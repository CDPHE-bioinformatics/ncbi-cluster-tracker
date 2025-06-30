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
        self.metadata = metadata
        self.metadata_truncated = self._truncate_metadata(metadata.copy()).fillna('')
        self.snp_matrix = self._create_snp_matrix()
        self.custom_labels = self._create_custom_labels()
        self.report = self._create_report()

    def _create_snp_matrix(self) -> ar.Group:
        """
        Create heatmap of SNP distance matrix annotated with metadata.
        """
        if self.cluster.filtered_matrix is None:
            text = f'⚠️ WARNING: {self.cluster.filtered_matrix_message}'
            snp_matrix = ar.Group(
                ar.Text(text, name=self.cluster.name, label=self.cluster.name),
            )
            return snp_matrix

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
            internal = self.metadata.query('source == "internal"').copy()
            internal['label'] = internal['id']
        else:
            internal = self.metadata.query('source == "internal"').copy()
            internal['label'] = internal['biosample']

        external = self.metadata.query('source == "external"').copy()
        external['label'] = external['biosample']

        def apply_prefix(row):
            if row.get('is_new') == 'yes':
                if row.get('source') == 'internal':
                    return f'⭐🆕{row.label}'
                return f'🆕{row.label}'
            elif row.get('source') == 'internal':
                return f'⭐{row.label}'
            return row.label

        combined = pd.concat([
            internal,
            external if not external.empty else None,
        ])
        combined['renamed'] = combined.apply(apply_prefix, axis=1)

        renamer = combined['renamed']
        renamer.index = combined.index  # ensure target_acc is the index
        matrix = matrix.rename(columns=renamer.to_dict(), index=renamer.to_dict())

        self.metadata = self.metadata.reset_index()
        return matrix
    
    def _add_metadata_to_matrix(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Add columns with basic metadata to the distance matrix.
        """
        matrix = (self.metadata_truncated[['target_acc', *self.PRIMARY_METADATA_COLS]]
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

    def _count_isolates(self) -> tuple[ar.Group, ar.Text | None]:
        """
        Summarize internal and external isolate counts.
        """
        internal_count = np.int64(self.clusters_df['internal_count'].iloc[0].item())
        external_count = np.int64(self.clusters_df['external_count'].iloc[0].item())
        total_count = internal_count + external_count
        internal_change = np.int64(self.clusters_df['internal_change'].iloc[0].item())
        external_change = np.int64(self.clusters_df['external_change'].iloc[0].item())
        total_change = internal_change + external_change

        if internal_change > 0:
            internal_upward = True
        elif internal_change < 0:
            internal_upward = False
        else:
            internal_change = None
            internal_upward = None

        if external_change > 0:
            external_upward = True
        elif external_change < 0:
            external_upward = False
        else:
            external_change = None
            external_upward = None

        if total_change > 0:
            total_upward = True
        elif total_change < 0:
            total_upward = False
        else:
            total_change = None
            total_upward = None

        # Isolate counts may differ between FTP and BigQuery datasets
        # if the FTP site updated before BigQuery or the FTP downloaded data is
        # outdated (second case less likely). In either case the user should be
        # warned of the mismatch.
        warning_message = []
        if self.cluster.external_isolates is not None:
            external_count_ftp = len(self.cluster.external_isolates)
        else:
            external_count_ftp = 0

        if external_count != external_count_ftp:
            warning_message = [ar.Text(
                '⚠️ WARNING: A more up-to-date version of this cluster may be ' \
                'available on the Pathogen Detection site with more ' \
                'internal and/or external isolates (visit the backup link).' \
            )]

        # List specific new isolates
        MAX_DISPLAY = 5
        if internal_change is not None:
            if 'id' in self.metadata.columns:
                internal_new = self.metadata.query(
                    'source == "internal" and is_new == "yes"'
                )[['id', 'isolate_id']]
            else:
                internal_new = self.metadata.query(
                    'source == "internal" and is_new == "yes"'
                )[['biosample', 'isolate_id']]
            internal_new = internal_new.iloc[:, 0].str.cat(internal_new.iloc[:, 1], sep=' / ').tolist()
            internal_list = 'New internal isolates added:\n'
            internal_list = f"{internal_list} - {'\n - '.join(internal_new[:MAX_DISPLAY])}"
            if len(internal_new) > MAX_DISPLAY:
                internal_list = f'{internal_list}\nand {internal_change - MAX_DISPLAY} mor.'
        else:
            internal_list = 'No change in internal isolate count.'

        if external_change is not None:
            external_new = self.metadata.query(
                'source == "external" and is_new == "yes"'
            )[['biosample', 'isolate_id']]
            external_new = external_new.iloc[:, 0].str.cat(external_new.iloc[:, 1], sep=' / ').tolist()
            external_list = 'New external isolates added:\n'
            external_list = f"{external_list} - {'\n - '.join(external_new[:MAX_DISPLAY])}\n"
            if len(external_new) > MAX_DISPLAY:
                external_list = f'{external_list}- and {external_change - MAX_DISPLAY} more'
        else:
            external_list = 'No change in external isolate count.'

        count_blocks = ar.Group(
            blocks=[
                *warning_message,
                ar.Group(
                    ar.BigNumber(
                        heading='Internal isolates',
                        value=internal_count,
                        change=internal_change,
                        is_upward_change=internal_upward,
                    ),
                    ar.BigNumber(
                        heading='External isolates',
                        value=external_count,
                        change=external_change,
                        is_upward_change=external_upward,
                    ),
                    ar.BigNumber(
                        heading='Total isolates',
                        value=total_count,
                        change=total_change,
                        is_upward_change=total_upward,
                    ),
                    ar.Text(internal_list),
                    ar.Text(external_list),
                    columns=3,
                ),
            ]
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

        # prefix label with star to avoid collision with Pathogen Detection
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
        count_blocks = self._count_isolates()
        tree_header = ar.HTML('<h3>NCBI Pathogen Detection</h3>')
        tree_url = self.clusters_df[
            self.clusters_df['cluster'] == self.cluster.name
        ]['tree_url'].item()
        browser_link_base = 'https://www.ncbi.nlm.nih.gov/pathogens/isolates/#'

        # Modify URL so that internal isolates are highlighted red
        tree_url_all_internal = f'{tree_url}?accessions={','.join(self.cluster.internal_isolates)}'
        backup_link_all_internal = f'{browser_link_base}{'%20'.join(self.cluster.internal_isolates)}'
        tree_links = (
            f'Links to tree:\n\n' \
            f'- [Highlight all internal isolates]({tree_url_all_internal}) ([Backup link]({backup_link_all_internal}))\n'
        )
        new_internals = self.metadata.query('source == "internal" and is_new == "yes"')['target_acc'].tolist()
        if new_internals:
            tree_url_new_internal = f'{tree_url}?accessions={','.join(new_internals)}'
            backup_link_new_internal = f'{browser_link_base}{'%20'.join(new_internals)}'
            tree_links = f'{tree_links}- [Highlight new internal isolates]({tree_url_new_internal}) ([Backup link]({backup_link_new_internal}))\n'
        new_all = self.metadata.query('is_new == "yes"')['target_acc'].tolist()
        if new_all:
            tree_url_new_all = f'{tree_url}?accessions={','.join(new_all)}'
            backup_link_new_internal = f'{browser_link_base}{"%20".join(new_all)}'
            tree_links = f'{tree_links}- [Highlight all new isolates]({tree_url_new_all}) ([Backup link]({backup_link_new_internal}))\n'
        cluster_base = self.cluster.name.split('.')[0]
        tree_links_block = ar.Text(tree_links)

        snp_header = ar.HTML('<h3>SNP distance matrix</h3>')
        count_graph_header = ar.HTML('<h3>Isolates by date added</h3>')
        count_graph = self._create_isolate_count_graph()
        table_header = ar.HTML('<h3>Isolate table</h3>')
        table = ar.DataTable(
            self.metadata
            .sort_values(by=['source', 'creation_date'], ascending=[False, False])
            .reset_index()
            .drop(columns='index')
        )
        if (
            self.clusters_df['change'].iloc[0] == 'new cluster'
            or not self.clusters_df['change'].iloc[0].startswith('+0')
        ):
            new = '🆕'
        else:
            new = ''
        blocks = [
                title,
                count_blocks,
                tree_header,
                tree_links_block,
                self.custom_labels,
                snp_header,
                self.snp_matrix,
                count_graph_header,
                count_graph,
                table_header,
                table,
        ]
        report = ar.Group(
            blocks = [b for b in blocks if b is not None],
            label=f'{new} {self.cluster.name} - {taxgroup_name}'
        )
        return report

def mark_new_isolates(
    isolates_df: pd.DataFrame,
    old_isolates_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """
    Add "is_new" column to `isolates_df` to indicate any new isolates added
    compared to `old_cluster_df`.
    """
    if old_isolates_df is None:
        isolates_df['is_new'] = 'yes'
        return isolates_df
    old_biosamples = set(old_isolates_df['biosample'])
    isolates_df['is_new'] = isolates_df['biosample'].apply(
        lambda x: 'yes' if x not in old_biosamples else 'no'
    )
    return isolates_df

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
    optional_id = ['id'] if 'id' in metadata.columns else []
    column_order = [
        'biosample',
        *optional_id,
        'isolate_id',
        'cluster',
        'source',
        'is_new',
        'scientific_name',
        'collection_date',
        'creation_date',
        'geo_loc_name',
        'isolation_source',
        'taxgroup_name',
        'bioproject_acc',
        'target_acc',
        'sra_id',
    ]
    extra_internal_cols = [c for c in metadata.columns if c not in column_order]
    extra_internal_cols.sort()
    column_order.extend(extra_internal_cols)
    metadata = metadata[column_order]
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
            clusters_df[clusters_df['cluster'] == cluster.name],
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
    clusters_df = clusters_df.astype({'internal_count': 'Int64', 'external_count': 'Int64'})

    if old_clusters_df is None:
        clusters_df['internal_change'] = 0
        clusters_df['external_change'] = 0
        clusters_df['change'] = 'new cluster'
        return clusters_df

    old_clusters_df['cluster_base'] = old_clusters_df['cluster'].str.split('.').str[0]
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
        compare_df[['cluster_base', 'change', 'internal_change', 'external_change']],
        clusters_df,
        how='right',
        on='cluster_base',
    )
    clusters_df['internal_change'] = clusters_df['internal_change'].fillna(0)
    clusters_df['external_change'] = clusters_df['external_change'].fillna(0)
    return clusters_df

def create_clusters_timeline_plot(
    metadata: pd.DataFrame,
    previous_max_date: datetime.date | None
) -> ar.Plot:
    """
    Create plot showing when each isolate was added to each cluster. Add
    vertical red line to plot showing when the previous report's last
    isolate was
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
    fig = px.scatter(
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
    if previous_max_date is not None:
        fig.add_vline(
            x=previous_max_date.timestamp() * 1000,  # https://github.com/plotly/plotly.py/issues/3065
            line_dash='dash',
            line_color='indianred',
            annotation_text='previous report'
        )
    plot = ar.Plot(fig)
    return plot, message


def write_final_report(
     clusters_df: pd.DataFrame,
     old_clusters_df: pd.DataFrame | None,
     clusters: list[cluster.Cluster], 
     metadata: pd.DataFrame,
     sample_sheet: str,
     compare_dir: str | None,
     command: str
) -> None:
    """
    Output final, standalone HTML report with all tables and plots. This
    function also outputs the clusters CSV.
    """
    logger.info('Generating HTML report...')
    clusters_df = add_counts(clusters_df, metadata)
    clusters_df = compare_counts(clusters_df, old_clusters_df)
    cluster_reports = create_cluster_reports(
        clusters,
        clusters_df,
        metadata
    )
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
    clusters_csv = os.path.join(
        os.environ['NCT_OUT_SUBDIR'],
        f'clusters_{os.environ["NCT_NOW"]}.csv'
    )
    clusters_df.to_csv(clusters_csv, index=False)

    cluster_page_blocks = [ar.HTML(f'<h2>Cluster report {os.environ["NCT_NOW"]}</h2>')]
    command_header = ar.Text('Command: ')
    command_block = ar.Code(code=command, language='javascript')
    cluster_page_blocks.extend([command_header, command_block])
    isolate_page_blocks = []
    if compare_dir is not None:
        header_2 = ar.Text(f'↔️ Comparing to {compare_dir}')
        cluster_page_blocks.append(header_2)
        isolate_page_blocks.append(header_2)

    clusters_table = ar.DataTable(
        clusters_df.sort_values('latest_added', ascending=False).reset_index().drop(columns='index')
    )
    previous_max_date = (
        pd.to_datetime(old_clusters_df['latest_added'].max()) if old_clusters_df is not None else None
    )
    clusters_timeline_plot, clusters_timeline_message = create_clusters_timeline_plot(
        metadata,
        previous_max_date,
    )
    cluster_page_blocks.extend([
        clusters_table,
        ar.HTML('<h3>Cluster timelines</h3>'),
    ])
    if clusters_timeline_message is not None: 
        cluster_page_blocks.append(clusters_timeline_message)
    cluster_page_blocks.append(clusters_timeline_plot)

    isolates_table = ar.DataTable(
        metadata
        .sort_values(by=['source', 'creation_date'], ascending=[False, False])
        .reset_index()
        .drop(columns='index')
    )
    isolate_page_blocks.append(isolates_table)
    missing_isolates = metadata.query('source == "internal" and is_new.isna()')['biosample'].tolist()
    missing_isolates.sort(reverse=True)
    metadata['is_new'] = metadata['is_new'].fillna('(no data)')
    # TODO: List specific reason for missing data, at least when possible
    if missing_isolates:
        missing_message = ar.Text(
            f'⚠️ WARNING: No cluster data displayed for the internal '\
            f'isolates listed below, either because they have not finished '\
            f'processing, failed quality control, or are otherwise not '\
            f'included in the system. More information about these isolates '\
            f'may be available in the Isolates Browser search results '\
            f'[here](https://www.ncbi.nlm.nih.gov/pathogens/isolates/#{"%20".join(missing_isolates)})'
            f'\n\n{", ".join(missing_isolates)}'
        )
        isolate_page_blocks.append(missing_message)

    cluster_report_blocks = [r.report for r in cluster_reports]
    cluster_report_blocks.sort(key=lambda r: r.label, reverse=True)
    report = ar.Blocks(
        ar.Page(blocks=cluster_page_blocks, title='Clusters'),
        ar.Page(blocks=isolate_page_blocks, title='Isolates'),
        ar.Page(
            ar.Select(blocks=cluster_report_blocks, type=ar.SelectType.DROPDOWN),
            title='Cluster details',
        )
    )

    ar.save_report(
        report,
        path=os.path.join(
            os.environ['NCT_OUT_SUBDIR'],
            f'clusters_{os.environ['NCT_NOW']}.html'
        ),
        standalone=False,
    )

