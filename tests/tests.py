import unittest

import pandas as pd

from ncbi_cluster_tracker import report

class testCompareCounts(unittest.TestCase):
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

        actual_df = report.compare_counts(clusters_df, old_clusters_df)
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
        
        actual_df = report.compare_counts(clusters_df, old_clusters_df)
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
