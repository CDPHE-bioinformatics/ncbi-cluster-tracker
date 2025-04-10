import argparse


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments from the user.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sample_sheet',
        help='CSV sample sheet with required "biosample" column and any additional metadata columns. Use "id" column for alternate isolate IDs',
    )
    parser.add_argument(
        '--compare-dir', '-c',
        help='Path to report directory to compare isolate counts, defaults to _current symlink'
    )
    parser.add_argument(
        '--use-local',
        help='Do not query BigQuery or NCBI, assumes data has already been downloaded',
        action=argparse.BooleanOptionalAction,
    )
    args = parser.parse_args()
    return args

