import os

import arakawa as ar # type: ignore
import pandas as pd
import plotly.express as px # type: ignore

from ncbi_cluster_tracker import cli
from ncbi_cluster_tracker import cluster
from ncbi_cluster_tracker import download
from ncbi_cluster_tracker import query
from ncbi_cluster_tracker import report


def main() -> None:
    args = cli.parse_args()
    sample_sheet_df = (pd
        .read_csv(args.sample_sheet)
        .set_index('biosample', verify_integrity=True)
    )
    biosamples = sample_sheet_df.index.to_list()

    # TODO: make CLI argument
    os.environ['NCT_OUT_DIR'] = 'outputs'
    os.makedirs(os.environ['NCT_OUT_DIR'], exist_ok=True)

    isolates_df, clusters_df = get_clusters(biosamples, args.use_local)
    download.download_cluster_files(clusters_df)
    clusters = cluster.create_clusters(sample_sheet_df, isolates_df, clusters_df)
    metadata = report.combine_metadata(sample_sheet_df, isolates_df)
    report.create_final_report(clusters_df, clusters, metadata)


def get_clusters(
    biosamples: list[str],
    is_local: bool
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fetch cluster data from NCBI's BigQuery `pdbrowser` dataset for the given
    `biosamples`, or use existing data if `is_local` is True.

    Return `isolates_df` DataFrame with isolate-level metadata, and
    `clusters_df` DataFrame with cluster-level metadata. Additionally, the
    DataFrames' data is written to a CSV in the output directory.
    """
    bq_isolates_csv = os.path.join(os.environ['NCT_OUT_DIR'], 'bq_isolates.csv')
    bq_clusters_csv = os.path.join(os.environ['NCT_OUT_DIR'], 'bq_clusters.csv')

    if not is_local:
        clusters = query.query_set_of_clusters(biosamples)
        isolates_df = query.query_isolates(clusters, biosamples)
        clusters_df = query.query_clusters(biosamples)
        isolates_df.to_csv(bq_isolates_csv, index=False)
        clusters_df.to_csv(bq_clusters_csv, index=False)
    else:
        isolates_df = pd.read_csv(bq_isolates_csv)
        clusters_df = pd.read_csv(bq_clusters_csv)

    return (isolates_df, clusters_df)


if __name__ == '__main__':
    main()