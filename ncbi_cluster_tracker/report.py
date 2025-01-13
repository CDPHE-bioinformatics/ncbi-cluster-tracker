import os

import arakawa as ar  # type: ignore
import numpy as np
import pandas as pd
import plotly.express as px  # type: ignore

import ncbi_cluster_tracker.cluster as cluster


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
        path = os.path.join(
            os.environ['NCT_OUT_DIR'],
            f'{self.cluster.name}_labels.txt'
        )
        custom_labels.to_csv(path, sep='\t', index=False, header=False) 
        attachment = ar.Attachment(file=path)
        return attachment 

    def _create_report(self) -> ar.Group:
        """
        Combine tables and visualizations for the cluster into an Arawaka
        Group to be displayed together in the report.
        """
        species = self.clusters_df[
            self.clusters_df['cluster'] == self.cluster.name
        ]['species'].item()

        title = ar.HTML(f'<h2><i>{species}</i> Cluster {self.cluster.name}</h2>')
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
            label=self.cluster.name
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


def create_final_report(
     clusters_df: pd.DataFrame,
     clusters: list[cluster.Cluster], 
     metadata: pd.DataFrame,
) -> None:
    cluster_reports = create_cluster_reports(
        clusters,
        clusters_df,
        metadata
    )

    internal_counts = (metadata
        .query('source == "internal"')
        .loc[:, 'cluster']
        .value_counts()
        .rename('internal_count')
    )
    clusters_df = clusters_df.merge(internal_counts, on='cluster')
    keep_cols = [
        'cluster',
        'internal_count',
        'total_count',
        'species',
        'earliest_added',
        'latest_added',
        'earliest_year_collected',
        'latest_year_collected',
        'tree_url',
    ]
    clusters_df = clusters_df[keep_cols]

    isolates_table = ar.DataTable(
        metadata[metadata['source'] == 'internal']
        .reset_index()
        .drop(columns='index')
    )

    clusters_table = ar.DataTable(clusters_df)
    clusters_gantt = ar.Plot(
        px.timeline(
            clusters_df.sort_values('latest_added', ascending=True),
            x_start='earliest_added',
            x_end='latest_added',
            y='cluster'
        )
    )

    cluster_report_blocks = [r.report for r in cluster_reports]

    report = ar.Blocks(
        ar.Page(clusters_table, clusters_gantt, title='Clusters'),
        ar.Page(isolates_table, title='Isolates'),
        ar.Page(
            ar.Select(
                blocks=cluster_report_blocks,
                type=ar.SelectType.DROPDOWN,
            ),
            title='Cluster details'
        )
    )

    #TODO Like multiqc, include JS directly in HTML (no href) so report
    # can be viewed offline (and security) 

    ar.save_report(report, path=os.path.join(os.environ['NCT_OUT_DIR'], 'clusters.html'))

