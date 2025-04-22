import datetime
import os

import arakawa as ar  # type: ignore
import numpy as np
import pandas as pd
import plotly.express as px  # type: ignore

import ncbi_cluster_tracker.cluster as cluster

from ncbi_cluster_tracker.logger import logger

class ClusterReport:
    """
    Contains Arakawa Report for "Cluster details" tab with SNP cluster metrics,
    distance matrix heatmap, and isolate count over time.
    """
    PRIMARY_METADATA_COLS = [
        'isolate_id',
        'collection_date',
        'geo_loc_name',
    ]

    def __init__(
        self, 
        cluster: cluster.Cluster,
        metadata: pd.DataFrame,
        clusters_df: pd.DataFrame,
    ):

        self.cluster = cluster
        self.clusters_df = clusters_df
        self.metadata = self._truncate_metadata(metadata.copy()).fillna('')
        self.snp_matrix = self._create_snp_matrix()
        self.custom_labels = self._create_custom_labels()
        self.report = self._create_report()

    def _create_snp_matrix(self) -> ar.Group:
        """
        Create heatmap of SNP distance matrix annotated with metadata.
        """
        if self.cluster.filtered_matrix is None:
            return None

        matrix = self._add_metadata_to_matrix(self.cluster.filtered_matrix)
        matrix = self._rename_matrix_ids(matrix)
        # Some isolates may be duplicated due to multiple assemblies, keep first
        matrix = matrix[~matrix.index.duplicated(keep='first')]
        matrix = matrix.loc[:, ~matrix.columns.duplicated(keep='first')]

        style = (matrix
            .style
            .background_gradient()
            .set_table_styles([
                {
                    'selector': 'th',
                    'props': [
                        ('width', '40px'),
                        ('padding', '3px'),
                        ('font-size', '12px')
                    ]
                },
                {
                    'selector': 'th.col_heading',
                    'props': [
                        ('writing-mode', 'vertical-rl'),
                        ('transform', 'rotateZ(180deg)'), 
                        ('height', '120px'),
                        ('vertical-align', 'middle'),
                        ('horizontal-align', 'left'),
                    ]
                },
                # {
                #     'selector': 'thead th.blank',
                #     'props': [('display', 'none')]
                # },
                # {
                #     'selector': 'thead th.index_name',
                #     'props': [('display', 'none')]
                # },
                {
                    'selector': 'td',
                    'props': [
                        ('font-size', '12px'),
                        ('padding', '3px'),
                        ('white-space', 'nowrap'),
                        ('overflow', 'hidden'),
                        ('text-overflow', 'ellipsis'),
                    ]
                }
            ])
            .set_sticky(axis=0)
        )
        if self.cluster.filtered_matrix_message is not None:
            text = f'⚠️ WARNING: {self.cluster.filtered_matrix_message}'
            snp_matrix = ar.Group(
                ar.Text(text, name=self.cluster.name, label=self.cluster.name),
                ar.Table(style, label=self.cluster.name)
            )
        else:
            snp_matrix = ar.Group(ar.Table(style, label=self.cluster.name))
        return snp_matrix

    def _rename_matrix_ids(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Use 'id' column from sample sheet to name internal isolates in the
        matrix and also mark them with a star. Use BioSample ID for external
        isolates and for internal isolates if 'id' column not provided.
        """
        self.metadata = self.metadata.set_index('target_acc')
        if 'id' in self.metadata.columns:
            # indicate internal isolates with a star
            internal_renamer = self.metadata.loc[
                self.metadata['source'] == 'internal',
                'id'
            ] 
        else:
            internal_renamer = self.metadata.loc[
                self.metadata['source'] == 'internal',
                'biosample'
            ]
        external_renamer = self.metadata.loc[
            self.metadata['source'] == 'external',
            'biosample'
        ]
        internal_renamer = '⭐' + internal_renamer  
        renamer = pd.concat([internal_renamer, external_renamer])
        matrix = matrix.rename(columns=dict(renamer), index=dict(renamer))
        self.metadata = self.metadata.reset_index()
        return matrix
    
    def _add_metadata_to_matrix(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Add columns with basic metadata to the distance matrix.
        """
        matrix = (self.metadata[['target_acc', *self.PRIMARY_METADATA_COLS]]
            .set_index('target_acc')
            .merge(
                matrix,
                right_index=True,
                left_on='target_acc',
                how='right',
            )
        )
        matrix.index.name = ''
        return matrix

    def _truncate_metadata(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """
        Limit string length displayed in report to prevent text from wrapping
        and wide column widths.
        """
        def truncate(text: str) -> str:
            if isinstance(text, str) and len(text) > 15:
                return f'{text[:7]}...{text[-7:]}'
            return text

        metadata[self.PRIMARY_METADATA_COLS] = (
            metadata[self.PRIMARY_METADATA_COLS].map(lambda x: truncate(x))
        )
        return metadata 

    def _create_isolate_counts(self) -> ar.Group:
        """
        Summarize internal and external isolate counts.
        """
        internal_count = len(self.cluster.internal_isolates)
        if self.cluster.external_isolates is not None:
            external_count = len(self.cluster.external_isolates)
        else:
            external_count = 0
        total_count = internal_count + external_count
        count_blocks = ar.Group(
            ar.BigNumber(
                heading='Internal isolates',
                value=internal_count
            ),
            ar.BigNumber(
                heading='External isolates',
                value=external_count
            ),
            ar.BigNumber(
                heading='Total isolates',
                value=total_count
            ),
            columns=3,
        )
        return count_blocks

    def _create_isolate_count_graph(self) -> ar.Plot:
        """
        Histogram of isolate counts over creation date (similar to 'epi curve').
        """
        # TODO: show sample ID on hover
        # TODO: collection dates, date added, or both?
        # If isolate only has year or no collection date, use date_added,
        # otherwise, use the collection date?
        count_graph = ar.Plot(
            px.histogram(
                self.metadata,
                x='creation_date',
                color='source',
            )
        )
        return count_graph

    def _create_custom_labels(self) -> ar.Attachment:
        """
        Create custom labels file that can be uploaded to the NCBI Pathogen
        Detection SNP tree viewer
        """
        cols = self.PRIMARY_METADATA_COLS.copy()
        if 'id' in self.metadata.columns:
            cols.insert(0, 'id')

        # prefix label with start to avoid collision with Pathogen Detection
        star_cols = {k: f'*{k}' for k in cols}

        custom_labels = (
            self.metadata
            .rename(columns=star_cols)
            .melt(
                id_vars=['target_acc'],
                value_vars=list(star_cols.values()),
            )
        )
        subdir = os.path.join(os.environ['NCT_OUT_SUBDIR'], 'labels')
        os.makedirs(subdir, exist_ok=True)
        path = os.path.join(subdir, f'{self.cluster.name}_labels.txt')
        custom_labels.to_csv(path, sep='\t', index=False, header=False) 
        attachment = ar.Attachment(file=path)
        return attachment 

    def _create_report(self) -> ar.Group:
        """
        Combine tables and visualizations for the cluster into an Arawaka
        Group to be displayed together in the report.
        """
        taxgroup_name = self.clusters_df[
            self.clusters_df['cluster'] == self.cluster.name
        ]['taxgroup_name'].item()

        title = ar.HTML(f'<h2><i>{taxgroup_name}</i> cluster {self.cluster.name}</h2>')
        count_blocks = self._create_isolate_counts()
        tree_header = ar.HTML('<h3>NCBI Pathogen Detection</h3>')
        tree_url = self.clusters_df[
            self.clusters_df['cluster'] == self.cluster.name
        ]['tree_url'].item()

        # Modify URL so that internal isolates are highlighted red
        # FIXME: File server may update before web interface updates, leading
        # to a blank page when visiting the link. Back up link?
        tree_url_highlight = f'{tree_url}?accessions={','.join(self.cluster.internal_isolates)}'
        tree_link = ar.Text(f'[Link to tree]({tree_url_highlight})')

        snp_header = ar.HTML('<h3>SNP distance matrix</h3>')
        count_graph_header = ar.HTML('<h3>Isolates by date added</h3>')
        count_graph = self._create_isolate_count_graph()
        report = ar.Group(
            title,
            count_blocks,
            tree_header,
            tree_link,
            self.custom_labels,
            snp_header,
            self.snp_matrix,
            count_graph_header,
            count_graph,
            label=f'{self.cluster.name} - {taxgroup_name}'
        )
        return report



def combine_metadata(
    sample_sheet_df: pd.DataFrame,
    isolates_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Combine user-provided data in sample sheet with sample metadata
    from NCBI.
    """
    metadata = sample_sheet_df.merge(
        isolates_df,
        how='outer',
        on='biosample',
        suffixes=('__internal', '__external'),
        indicator='source'
    )
    # Use sample sheet values for matching columns in sample sheet and BigQuery
    internal_merged_cols = [c for c in metadata.columns if '__internal' in c]
    for internal_col in internal_merged_cols:
        base_col = internal_col[:-10]  # without __suffix
        external_col = f'{base_col}__external'
        metadata[base_col] = np.where(
            ~metadata[internal_col].isnull(),
            metadata[internal_col],
            metadata[external_col],
        )
        metadata = metadata.drop([internal_col, external_col], axis=1)

    metadata['source'] = metadata['source'].map({
        'left_only': 'internal',
        'right_only': 'external',
        'both': 'internal',
    })
    
    return metadata


def create_cluster_reports(
    clusters: list[cluster.Cluster],
    clusters_df: pd.DataFrame,
    metadata: pd.DataFrame
) -> list[ClusterReport]:
    """
    Create a ClusterReport for all clusters.
    """
    cluster_reports: list[ClusterReport] = []
    for cluster in clusters:
        cluster_report = ClusterReport(
        cluster,
        metadata[metadata['cluster'] == cluster.name],
        clusters_df,
    )
        cluster_reports.append(cluster_report)
    return cluster_reports


def add_counts(clusters_df: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Add external and total isolate counts to `clusters_df` DataFrame.
    """
    if 'total_count' not in clusters_df.columns:
        clusters_df['total_count'] = (
            clusters_df['internal_count'] + clusters_df['external_count']
        )
    else:
        internal_counts = (metadata
            .query('source == "internal"')
            .loc[:, 'cluster']
            .value_counts()
            .rename('internal_count')
        )
        clusters_df = clusters_df.merge(internal_counts, on='cluster')
        clusters_df['external_count'] =  (
            clusters_df['total_count'] - clusters_df['internal_count']
        )
    return clusters_df


def compare_counts(
    clusters_df: pd.DataFrame,
    old_clusters_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """
    Compare cluster isolate counts between previous report, if provided, and add
    differences to the `clusters_df` DataFrame.
    """
    clusters_df['cluster_base'] = clusters_df['cluster'].str.split('.').str[0]

    if old_clusters_df is None:
        clusters_df['change'] = 'new cluster'
        return clusters_df

    old_clusters_df['cluster_base'] = old_clusters_df['cluster'].str.split('.').str[0]
    clusters_df = clusters_df.astype({'internal_count': 'Int64', 'external_count': 'Int64'})
    old_clusters_df = old_clusters_df.astype({'internal_count': 'Int64', 'external_count': 'Int64'})

    compare_df = pd.merge(
        clusters_df[['cluster_base', 'internal_count', 'external_count']],
        old_clusters_df[['cluster_base', 'internal_count', 'external_count']],
        on='cluster_base',
        how='left',
        suffixes=('_new', '_old'),
    )
    compare_df['internal_change'] = compare_df['internal_count_new'] - compare_df['internal_count_old']
    compare_df['external_change'] = compare_df['external_count_new'] - compare_df['external_count_old']

    def create_change_column(row):
        if pd.isna(row['internal_change']) or pd.isna(row['external_change']):
            return 'new cluster'
        internal_prefix = '+' if row['internal_change'] >= 0 else ''
        external_prefix = '+' if row['external_change'] >= 0 else ''
        return f'{internal_prefix}{row["internal_change"]} / {external_prefix}{row["external_change"]}'

    compare_df['change'] = compare_df.apply(
        create_change_column,
        axis=1
    )
    if 'change' in clusters_df.columns:
        clusters_df = clusters_df.drop(columns='change')
    clusters_df = pd.merge(
        compare_df[['cluster_base', 'change']],
        clusters_df,
        how='right',
        on='cluster_base',
    )
    return clusters_df

def create_clusters_timeline_plot(metadata: pd.DataFrame) -> ar.Plot:
    """
    Create plot showing when each isolate was added to each cluster.
    """
    MAX_DISPLAY = 15

    # px.strip() wasn't jittering points with datetime, so manually jittering
    metadata_jittered = metadata.copy()
    metadata_jittered = metadata_jittered[metadata_jittered['cluster'].notna()]
    metadata_jittered = metadata_jittered.sort_values(by='creation_date', ascending=False)
    metadata_jittered['cluster_ticktext'] = (
        metadata_jittered['cluster'].str.cat(metadata_jittered['taxgroup_name'], sep='<br>')
    )
    message = None
    if metadata_jittered['cluster'].nunique() > MAX_DISPLAY:
        top_clusters = metadata_jittered['cluster'].unique()[:MAX_DISPLAY]
        metadata_jittered = metadata_jittered[metadata_jittered['cluster'].isin(top_clusters)]
        message = f'Top {MAX_DISPLAY} most recently updated clusters displayed'
    metadata_jittered['cluster_jittered'] = (
        pd.Categorical(metadata_jittered['cluster'], categories=metadata_jittered['cluster'].unique(), ordered=True).codes
        + np.random.uniform(-0.10, 0.10, size=len(metadata_jittered))
    )
    plot = ar.Plot(
        px.scatter(
            metadata_jittered,
            x='creation_date',
            y='cluster_jittered',
            color='source',
            height=90 * len(metadata_jittered['cluster'].unique()),
            hover_data=['cluster', 'isolate_id'],
        ).update_yaxes(
            tickmode='array',
            tickvals=list(range(len(metadata_jittered['cluster'].unique()))),
            ticktext=metadata_jittered['cluster_ticktext'].unique(),
            autorange='reversed',
        ).update_traces(
            hovertemplate='<b>Isolate ID:</b> %{customdata[1]}<br><b>Date:</b> %{x}'
        ).update_layout(
            yaxis_title='cluster',
            xaxis={'side': 'top'},
        )
    )
    return plot, message


def write_final_report(
     clusters_df: pd.DataFrame,
     old_clusters_df: pd.DataFrame | None,
     clusters: list[cluster.Cluster], 
     metadata: pd.DataFrame,
) -> None:
    """
    Output final, standalone HTML report with all tables and plots. This
    function also outputs the clusters CSV.
    """
    cluster_reports = create_cluster_reports(
        clusters,
        clusters_df,
        metadata
    )
    clusters_df = add_counts(clusters_df, metadata)
    clusters_df = compare_counts(clusters_df, old_clusters_df)
    keep_cols = [
        'cluster',
        'taxgroup_name',
        'internal_count',
        'external_count',
        'change',
        'latest_added',
        'earliest_added',
        'earliest_year_collected',
        'latest_year_collected',
        'tree_url',
    ]
    clusters_df = clusters_df[keep_cols]
    bq_clusters_csv = os.path.join(
        os.environ['NCT_OUT_SUBDIR'],
        f'bq_clusters_{os.environ["NCT_NOW"]}.csv'
    )
    clusters_df.to_csv(bq_clusters_csv, index=False)

    isolates_table = ar.DataTable(
        metadata
        .sort_values(by=['source', 'creation_date'], ascending=[False, False])
        .reset_index()
        .drop(columns='index')
    )

    clusters_table = ar.DataTable(
        clusters_df.sort_values('latest_added', ascending=False).reset_index().drop(columns='index')
    )

    clusters_timeline_plot, clusters_timeline_message = create_clusters_timeline_plot(metadata)
    cluster_report_blocks = [r.report for r in cluster_reports]
    cluster_page_blocks = [
        clusters_table,
        ar.HTML('<h3>Cluster timelines</h3>'),
    ]
    if clusters_timeline_message is not None: 
        cluster_page_blocks.append(clusters_timeline_message)
    cluster_page_blocks.append(clusters_timeline_plot)

    report = ar.Blocks(
        ar.Page(
            blocks=cluster_page_blocks,
            title='Clusters',
        ),
        ar.Page(isolates_table, title='Isolates'),
        ar.Page(
            ar.Select(
                blocks=cluster_report_blocks,
                type=ar.SelectType.DROPDOWN,
            ),
            title='Cluster details'
        )
    )

    ar.save_report(
        report,
        path=os.path.join(
            os.environ['NCT_OUT_SUBDIR'],
            f'clusters_{os.environ['NCT_NOW']}.html'
        ),
        standalone=True,
    )

