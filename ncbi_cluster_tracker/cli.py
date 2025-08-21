import argparse

from importlib.metadata import version

from typing import Sequence

def parse_args(command: Sequence[str]) -> argparse.Namespace:
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
        help='Do not query BigQuery or NCBI, assumes data has already been downloaded to --out-dir.',
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        '--browser-file',
        # TODO link to instructions
        help='Path to isolates TSV or CSV downloaded from the Pathogen Detection Isolates Browser with information for all internal and external isolates. When specified, data in file will be used instead of querying the BigQuery dataset.'
    )
    parser.add_argument(
        '--keep-snp-files',
        help='Keep downloaded SNP and tree files in the output directory. By default, files are deleted after processing.',
        action='store_true',
    )
    parser.add_argument(
        '--amr',
        help='Include AMR tab in report with antimicrobial resistance genes detected by AMRFinderPlus.',
        action='store_true',
    )
    parser.add_argument(
        '--filter-amr',
        help='Only include AMR genes in provided comma-separated list of CLASS:SUBCLASS pairs in the AMR tab. Also adds filtered_amr column to Isolates and Cluster details tab and matching genes to tree labels',
        type=lambda s: [i for i in s.split(',')],
    )
    parser.add_argument(
        '--version', '-v',
        help='Print the version of ncbi-cluster-tracker and exit.',
        action='version',
        version=version('ncbi-cluster-tracker'),
    )
    mutex_group_compare = parser.add_mutually_exclusive_group()
    mutex_group_compare.add_argument(
        '--compare-dir',
        help='Path to previous output directory to detect and report new isolates.',
    )
    args = parser.parse_args(command)

    if args.retry and not args.out_dir:
        parser.error('--retry flag requires --out_dir argument')

    if args.filter_amr:
        if not args.amr:
            parser.error('--filter-amr argument requires --amr flag')
        for item in args.filter_amr:
            if ':' not in item[1:-1]:
                parser.error('Each element in --filter-amr list must be in the form CLASS:SUBCLASS')
    
    return args

