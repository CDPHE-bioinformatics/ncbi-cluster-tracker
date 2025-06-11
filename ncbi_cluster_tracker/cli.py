import argparse

from importlib.metadata import version

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments from the user.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'sample_sheet',
        help='Path to sample sheet CSV with required "biosample" column and any additional metadata columns. Use "id" column for alternate isolate IDs.',
    )
    parser.add_argument(
        '--out-dir', '-o',
        help='Path to directory to store outputs. Defaults to "./outputs/" if not specified.'
    )
    parser.add_argument(
        '--retry',
        help='Do not query BigQuery or NCBI, assumes data has already been downloaded to --out-dir or directory with most recent timestamp.',
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        '--browser-file',
        # TODO link to instructions
        help='Path to isolates TSV or CSV downloaded from the Pathogen Detection Isolates Browser with information for all internal and external isolates. When specified, data in file will be used instead of querying the BigQuery dataset.'
    )
    parser.add_argument(
        '--version', '-v',
        help='Print the version of ncbi_cluster_tracker and exit.',
        action='version',
        version=version('ncbi-cluster-tracker'),
    )
    mutex_group_compare = parser.add_mutually_exclusive_group()
    mutex_group_compare.add_argument(
        '--compare-dir',
        help='Path to previous output directory to detect and report new isolates. Defaults to directory inside --out-dir with most recent timestamp if not specified.',
    )
    mutex_group_compare.add_argument(
        '--no-compare',
        help='Do not compare to most recent output directory, all clusters and isolates will be considered "new".',
        action='store_true',
    )
    args = parser.parse_args()
    return args

