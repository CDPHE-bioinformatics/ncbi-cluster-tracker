import datetime
import glob
import logging
import os
import re

import arakawa as ar # type: ignore
import pandas as pd
import plotly.express as px # type: ignore

from ncbi_cluster_tracker import cli
from ncbi_cluster_tracker import cluster
from ncbi_cluster_tracker import download
from ncbi_cluster_tracker import query
from ncbi_cluster_tracker import report

from ncbi_cluster_tracker.logger import logger

def main() -> None:
    logging.basicConfig(level=logging.INFO)
    args = cli.parse_args()
    sample_sheet_df = (pd
        .read_csv(args.sample_sheet)
        .set_index('biosample', verify_integrity=True)
    )
    biosamples = sample_sheet_df.index.to_list()

    compare_dir = args.compare_dir
    if compare_dir is None:
        try:
            # select directory with most recent timestamp
            dirs = glob.glob('outputs/*')
            valid_dirs = [d for d in dirs if re.search(r'\d{8}_\d{6}', d)]
            compare_dir = max(valid_dirs)
            logger.info(f'Comparing to {compare_dir}')
        except FileNotFoundError:
            logger.info('No comparison directory found.')
    else:
        if not os.path.isdir(compare_dir):
            raise ValueError(f'Directory "{compare_dir}" does not exist.')
    
    if compare_dir is None:
        old_clusters_df = None
    else:
        old_clusters_glob = glob.glob(os.path.join(compare_dir, 'bq_clusters*.csv'))
        if not old_clusters_glob:
            raise ValueError(f'Could not find bq_clusters file in {compare_dir}')
        if len(old_clusters_glob) > 1:
            raise ValueError(f'Multiple bq_clusters files found in {compare_dir}')
        old_clusters_df = pd.read_csv(old_clusters_glob[0])

    # TODO: make CLI argument
    os.environ['NCT_NOW'] = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    os.environ['NCT_OUT_DIR'] = os.path.join('outputs', os.environ['NCT_NOW'])
    os.makedirs(os.environ['NCT_OUT_DIR'], exist_ok=True)

    isolates_df, clusters_df = get_clusters(biosamples, args.use_local)
    download.download_cluster_files(clusters_df)
    clusters = cluster.create_clusters(sample_sheet_df, isolates_df, clusters_df)
    metadata = report.combine_metadata(sample_sheet_df, isolates_df)
    report.create_final_report(clusters_df, old_clusters_df, clusters, metadata)


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
    bq_isolates_csv = os.path.join(
        os.environ['NCT_OUT_DIR'],
        f'bq_isolates_{os.environ["NCT_NOW"]}.csv'
    )
    bq_clusters_csv = os.path.join(
        os.environ['NCT_OUT_DIR'],
        f'bq_clusters_{os.environ["NCT_NOW"]}.csv'
    )

    if not is_local:
        clusters = query.query_set_of_clusters(biosamples)
        isolates_df = query.query_isolates(clusters, biosamples)
        clusters_df = query.query_clusters(biosamples)
        isolates_df.to_csv(bq_isolates_csv, index=False)
    else:
        isolates_df = pd.read_csv(bq_isolates_csv)
        clusters_df = pd.read_csv(bq_clusters_csv)

    # compare to previous counts
    return (isolates_df, clusters_df)


if __name__ == '__main__':
    main()