
import argparse
import datetime
import glob
import os
import re
import sys

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
    command = f'{os.path.basename(sys.argv[0])} {" ".join(sys.argv[1:])}'
    args = cli.parse_args()
    sample_sheet_df = (pd
        .read_csv(args.sample_sheet, dtype={'id': 'string'})
        .set_index('biosample', verify_integrity=True)
    )
    biosamples = sample_sheet_df.index.to_list()

    out_dir = 'outputs' if args.out_dir is None else args.out_dir
    compare_dir, latest_dir = find_existing_dirs(args, out_dir)
    
    if compare_dir is None or args.no_compare:
        old_clusters_df = None
        old_isolates_df = None
    else:
        old_clusters_glob = glob.glob(os.path.join(compare_dir, '*clusters*.csv'))
        old_isolates_glob = glob.glob(os.path.join(compare_dir, '*isolates*.csv'))
        if not old_clusters_glob:
            raise FileNotFoundError(f'Could not find clusters CSV file in {compare_dir}')
        if len(old_clusters_glob) > 1:
            raise ValueError(f'Multiple clusters CSV files found in {compare_dir}')
        if not old_isolates_glob:
            raise FileNotFoundError(f'Could not find isolates CSV file in {compare_dir}')
        if len(old_isolates_glob) > 1:
            raise ValueError(f'Multiple isolates CSV files found in {compare_dir}')
        old_clusters_df = pd.read_csv(old_clusters_glob[0])
        old_isolates_df = pd.read_csv(old_isolates_glob[0])

    if not args.retry:
        os.environ['NCT_NOW'] = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        os.environ['NCT_OUT_SUBDIR'] = os.path.join(out_dir, os.environ['NCT_NOW'])
        os.makedirs(os.environ['NCT_OUT_SUBDIR'], exist_ok=True)
        if args.browser_file is not None:
            isolates_df, clusters_df = get_clusters(biosamples, args.browser_file)
        else:
            isolates_df, clusters_df = get_clusters(biosamples, 'bigquery')
        download.download_cluster_files(clusters_df)
    else:
        if latest_dir is None:
            raise FileNotFoundError(f'Could not find existing data at output directory {out_dir} for --retry')
        os.environ['NCT_OUT_SUBDIR'] = latest_dir
        os.environ['NCT_NOW'] = os.path.basename(os.environ['NCT_OUT_SUBDIR'])
        logger.info(f'Retrying with {os.environ["NCT_OUT_SUBDIR"]}')
        isolates_df, clusters_df = get_clusters(biosamples, 'local')
    
    clusters_df['tree_url'] = clusters_df.apply(download.build_tree_viewer_url, axis=1)
    clusters = cluster.create_clusters(sample_sheet_df, isolates_df, clusters_df)
    isolates_df = report.mark_new_isolates(isolates_df, old_isolates_df)
    metadata = report.combine_metadata(sample_sheet_df, isolates_df)
    report.write_final_report(
        clusters_df,
        old_clusters_df,
        clusters,
        metadata,
        args.sample_sheet,
        compare_dir,
        command,
    )

def find_existing_dirs(args: argparse.Namespace, out_dir: str) -> tuple[str, list[str]]:
    """
    Return sub-directory `compare_dir` to compare results (if any) and the 
    directory with most recent timestamp inside `out_dir`.
    """
    dirs = glob.glob(os.path.join(out_dir, '*'))
    valid_dirs = [d for d in dirs if re.search(r'\d{8}_\d{6}', d)]

    compare_dir = args.compare_dir
    if compare_dir is None and not args.no_compare:
        try:
            latest_dirs = sorted(valid_dirs)
            if args.retry:
                compare_dir = latest_dirs[-2]
            else:
                compare_dir = latest_dirs[-1]
            logger.info(f'Comparing to {compare_dir}')
        except (IndexError, ValueError):
            logger.info('No comparison directory found.')
            compare_dir = None
    elif compare_dir is not None and not os.path.isdir(compare_dir):
        raise ValueError(f'Directory "{compare_dir}" does not exist.')
    if not valid_dirs:
        return compare_dir, None
    return compare_dir, max(valid_dirs)


def get_clusters(
    biosamples: list[str],
    data_location: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fetch cluster data from NCBI's BigQuery `pdbrowser` dataset for the given
    `biosamples` if `data_location` is 'bigquery', or read from existing output
    CSV files if `data_location` is 'local`, otherwise read from --browser-tsv
    file.

    Return `isolates_df` DataFrame with isolate-level metadata, and
    `clusters_df` DataFrame with cluster-level metadata. Additionally, the
    DataFrames' data is written to a CSV in the output directory.
    """
    isolates_csv = os.path.join(
        os.environ['NCT_OUT_SUBDIR'],
        f'isolates_{os.environ["NCT_NOW"]}.csv'
    )
    clusters_csv = os.path.join(
        os.environ['NCT_OUT_SUBDIR'],
        f'clusters_{os.environ["NCT_NOW"]}.csv'
    )

    if data_location == 'bigquery':
        clusters = query.query_set_of_clusters(biosamples)
        isolates_df = query.query_isolates(clusters, biosamples)
        # TODO: query_clusters() should be replaceable with cluster_df_from_isolates_df()
        clusters_df = query.query_clusters(biosamples)
        isolates_df.to_csv(isolates_csv, index=False)
    elif data_location == 'local':
        isolates_df = pd.read_csv(isolates_csv)
        clusters_df = pd.read_csv(clusters_csv)
    else:
        if data_location.endswith('.tsv'):
            browser_df = pd.read_csv(data_location, sep='\t', low_memory=False, on_bad_lines='warn')
        elif data_location.endswith('.csv'):
            browser_df = pd.read_csv(data_location, low_memory=False, on_bad_lines='warn')
        else:
            raise ValueError(f'Invalid file type (must be .tsv or .csv): {data_location}')
        isolates_df = query.isolates_df_from_browser_df(browser_df)
        clusters_df = query.cluster_df_from_isolates_df(isolates_df)
        isolates_df.to_csv(isolates_csv, index=False)

    return (isolates_df, clusters_df)


if __name__ == '__main__':
    main()