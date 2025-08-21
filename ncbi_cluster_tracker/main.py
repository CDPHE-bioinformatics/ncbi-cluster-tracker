
import argparse
import datetime
import glob
import os
import re
import shutil
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
    args = cli.parse_args(sys.argv[1:])
    sample_sheet_df = (pd
        .read_csv(args.sample_sheet, dtype={'id': 'string'})
        .set_index('biosample', verify_integrity=True)
    )
    biosamples = sample_sheet_df.index.to_list()

    out_dir = 'outputs' if args.out_dir is None else args.out_dir
    
    if args.compare_dir is None:
        old_clusters_df = None
        old_isolates_df = None
    else:
        old_clusters_glob = glob.glob(os.path.join(args.compare_dir, '*clusters*.csv'))
        old_isolates_glob = glob.glob(os.path.join(args.compare_dir, '*isolates*.csv'))
        if not old_clusters_glob:
            raise FileNotFoundError(f'Could not find clusters CSV file in {args.compare_dir}')
        if len(old_clusters_glob) > 1:
            raise ValueError(f'Multiple clusters CSV files found in {args.compare_dir}')
        if not old_isolates_glob:
            raise FileNotFoundError(f'Could not find isolates CSV file in {args.compare_dir}')
        if len(old_isolates_glob) > 1:
            raise ValueError(f'Multiple isolates CSV files found in {args.compare_dir}')
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
        download.download_cluster_files(clusters_df, args.keep_snp_files)
    else:
        if not os.path.isdir(out_dir):
            raise FileNotFoundError(f'Could not find existing output directory {out_dir} for --retry')
        os.environ['NCT_OUT_SUBDIR'] = out_dir
        os.environ['NCT_NOW'] = os.path.basename(out_dir.rstrip(os.sep))
        logger.info(f'Retrying with {os.environ["NCT_OUT_SUBDIR"]}')
        isolates_df, clusters_df = get_clusters(biosamples, 'local')
    
    if args.amr:
        amr_ref_df = download.download_amr_reference_file()  
        amr_df = query.create_amr_df(isolates_df, amr_ref_df)
        if args.filter_amr:
            amr_df = query.filter_amr_df(amr_df, args.filter_amr)
    else:
        amr_df = None

    clusters_df['tree_url'] = clusters_df.apply(download.build_tree_viewer_url, axis=1)
    clusters = cluster.create_clusters(sample_sheet_df, isolates_df, clusters_df)
    if not args.keep_snp_files:
        shutil.rmtree(os.path.join(os.environ['NCT_OUT_SUBDIR'], 'snps'))
    isolates_df = report.mark_new_isolates(isolates_df, old_isolates_df)
    metadata = report.combine_metadata(sample_sheet_df, isolates_df, amr_df, args)
    report.write_final_report(
        clusters_df,
        old_clusters_df,
        clusters,
        metadata,
        amr_df,
        args.compare_dir,
        command,
    )


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
        try:
            isolates_df = pd.read_csv(isolates_csv)
            clusters_df = pd.read_csv(clusters_csv)
        except FileNotFoundError as e:
            message = f'Existing isolates/clusters CSV not found in {os.environ["NCT_OUT_SUBDIR"]}'
            raise Exception(message) from e
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

    return isolates_df, clusters_df


if __name__ == '__main__':
    main()