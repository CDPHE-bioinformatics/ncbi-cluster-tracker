import unittest

import pandas as pd

from ncbi_cluster_tracker import report
from ncbi_cluster_tracker import query

class testCompareCounts(unittest.TestCase):
    COMPARE_COLS = ['cluster', 'cluster_base', 'internal_count', 'external_count', 'change']
    def test_no_compare(self):
        clusters_df = pd.DataFrame(
            [
                ['cluster1.1', 8, 0],
                ['cluster2.1', 3, 4],
                ['cluster3.1', 2, 9],
                ['cluster4.1', 1, 1],
            ],
            columns=['cluster', 'internal_count', 'external_count'],
        ).astype({'internal_count': 'Int64', 'external_count': 'Int64'})
        old_clusters_df = None

        actual_df = report.compare_counts(clusters_df, old_clusters_df)[self.COMPARE_COLS]
        expected_df = pd.DataFrame(
            [
                ['cluster1.1', 'cluster1', 8, 0, 'new cluster'],
                ['cluster2.1', 'cluster2', 3, 4, 'new cluster'],
                ['cluster3.1', 'cluster3', 2, 9, 'new cluster'],
                ['cluster4.1', 'cluster4', 1, 1, 'new cluster'],
            ],
            columns=['cluster', 'cluster_base', 'internal_count', 'external_count', 'change'],
        ).astype({'internal_count': 'Int64', 'external_count': 'Int64'})
        pd.testing.assert_frame_equal(actual_df, expected_df, check_like=True)

    def test_changes(self):
        clusters_df = pd.DataFrame(
            [
                ['cluster1.1', 8, 0],
                ['cluster2.1', 3, 4],
                ['cluster3.1', 2, 9],
                ['cluster4.1', 1, 1],
            ],
            columns=['cluster', 'internal_count', 'external_count'],
        ).astype({'internal_count': 'Int64', 'external_count': 'Int64'})
        old_clusters_df = pd.DataFrame(
            [
                ['cluster1.0', 7, 0],
                ['cluster2.0', 2, 5],
                ['cluster3.0', 2, 9],
                ['cluster5.0', 6, 2],
            ],
            columns=['cluster', 'internal_count', 'external_count'],
        ).astype({'internal_count': 'Int64', 'external_count': 'Int64'})
        
        actual_df = report.compare_counts(clusters_df, old_clusters_df)[self.COMPARE_COLS]
        expected_df = pd.DataFrame(
            [
                ['cluster1.1', 'cluster1', 8, 0, '+1 / +0'],
                ['cluster2.1', 'cluster2', 3, 4, '+1 / -1'],
                ['cluster3.1', 'cluster3', 2, 9, '+0 / +0'],
                ['cluster4.1', 'cluster4', 1, 1, 'new cluster'],
            ],
            columns=['cluster', 'cluster_base', 'internal_count', 'external_count', 'change'],
        ).astype({'internal_count': 'Int64', 'external_count': 'Int64'})
        pd.testing.assert_frame_equal(actual_df, expected_df, check_like=True)


class TestIsolatesDfFromBrowserDf(unittest.TestCase):
    COMPARE_COLS = ['isolate_id', 'biosample', 'target_acc', 'cluster',
                    'sra_id', 'isolation_source', 'geo_loc_name',
                    'collection_date', 'creation_date', 'taxgroup_name',
                    'scientific_name', 'bioproject_acc', 'amr_genotypes']
    
    def test_df_match_bigquery(self):
        browser_df = pd.read_csv('tests/data/pdbrowser_20250618.tsv', sep='\t')
        actual_df = query.isolates_df_from_browser_df(browser_df)
        actual_df = actual_df.sort_values(by='target_acc').reset_index(drop=True)
        expected_df = pd.read_csv(
            'tests/data/20250618_210234/isolates_20250618_210234.csv',
            dtype={'collection_date': 'string'},
        )
        expected_df = expected_df[self.COMPARE_COLS]
        expected_df = expected_df.sort_values(by='target_acc').reset_index(drop=True)
        pd.testing.assert_frame_equal(actual_df, expected_df, check_like=True)


class TestClusterDfFromIsolatesDf(unittest.TestCase):
    COMPARE_COLS = ['cluster', 'total_count', 'taxgroup_name',
                    'earliest_added', 'latest_added', 'earliest_year_collected',
                    'latest_year_collected']
    def test_df_match_bigquery(self):
        browser_df = pd.read_csv('tests/data/pdbrowser_20250618.tsv', sep='\t')
        isolates_df = query.isolates_df_from_browser_df(browser_df)
        actual_df = query.cluster_df_from_isolates_df(isolates_df)
        actual_df = actual_df.sort_values(by='cluster').reset_index(drop=True)
        expected_df = pd.read_csv(
            'tests/data/20250618_210234/clusters_20250618_210234.csv',
            dtype={
                'earliest_year_collected': 'string',
                'latest_year_collected': 'string',
            }
        )
        expected_df['total_count'] = expected_df['internal_count'] + expected_df['external_count']
        expected_df = expected_df[self.COMPARE_COLS]
        expected_df = expected_df.sort_values(by='cluster').reset_index(drop=True)
        pd.testing.assert_frame_equal(actual_df, expected_df, check_like=True)