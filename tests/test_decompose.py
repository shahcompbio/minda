import pytest
import pandas as pd
from minda.decompose import _get_paired_info_dfs


def _make_info_df(records):
    """Helper to build an info_df matching the shape expected by _get_paired_info_dfs."""
    return pd.DataFrame(records)


class TestGetPairedInfoDfsInsertions:
    """Test that INS records get END = POS + 1, not POS + abs(SVLEN)."""

    def test_ins_end_is_pos_plus_one(self):
        """For an insertion, the second breakpoint should be POS + 1."""
        info_df = _make_info_df([{
            '#CHROM': 'chr1', 'POS': 1000, 'ID': 'ins1',
            'ALT': '<INS>', 'FILTER': 'PASS',
            'INFO': 'SVTYPE=INS;SVLEN=500;END=1001',
            'SVTYPE': 'INS', 'VAF': '*',
        }])
        df_1, df_2 = _get_paired_info_dfs(info_df)
        assert df_2.iloc[0]['POS'] == 1001

    def test_ins_end_not_pos_plus_svlen(self):
        """Ensure INS END is not computed as POS + abs(SVLEN)."""
        info_df = _make_info_df([{
            '#CHROM': 'chr1', 'POS': 1000, 'ID': 'ins2',
            'ALT': '<INS>', 'FILTER': 'PASS',
            'INFO': 'SVTYPE=INS;SVLEN=500;END=1001',
            'SVTYPE': 'INS', 'VAF': '*',
        }])
        df_1, df_2 = _get_paired_info_dfs(info_df)
        assert df_2.iloc[0]['POS'] != 1500

    def test_ins_large_svlen(self):
        """A large SVLEN should still yield END = POS + 1 for INS."""
        info_df = _make_info_df([{
            '#CHROM': 'chr1', 'POS': 5000, 'ID': 'ins3',
            'ALT': '<INS>', 'FILTER': 'PASS',
            'INFO': 'SVTYPE=INS;SVLEN=10000;END=5001',
            'SVTYPE': 'INS', 'VAF': '*',
        }])
        df_1, df_2 = _get_paired_info_dfs(info_df)
        assert df_2.iloc[0]['POS'] == 5001


class TestGetPairedInfoDfsNonInsertions:
    """Test that non-INS SV types still use POS + abs(SVLEN)."""

    def test_del_end_is_pos_plus_svlen(self):
        """DEL should keep END = POS + abs(SVLEN)."""
        info_df = _make_info_df([{
            '#CHROM': 'chr1', 'POS': 1000, 'ID': 'del1',
            'ALT': '<DEL>', 'FILTER': 'PASS',
            'INFO': 'SVTYPE=DEL;SVLEN=-500;END=1500',
            'SVTYPE': 'DEL', 'VAF': '*',
        }])
        df_1, df_2 = _get_paired_info_dfs(info_df)
        assert df_2.iloc[0]['POS'] == 1500

    def test_dup_end_is_pos_plus_svlen(self):
        """DUP should keep END = POS + abs(SVLEN)."""
        info_df = _make_info_df([{
            '#CHROM': 'chr1', 'POS': 2000, 'ID': 'dup1',
            'ALT': '<DUP>', 'FILTER': 'PASS',
            'INFO': 'SVTYPE=DUP;SVLEN=300;END=2300',
            'SVTYPE': 'DUP', 'VAF': '*',
        }])
        df_1, df_2 = _get_paired_info_dfs(info_df)
        assert df_2.iloc[0]['POS'] == 2300

    def test_inv_end_is_pos_plus_svlen(self):
        """INV should keep END = POS + abs(SVLEN)."""
        info_df = _make_info_df([{
            '#CHROM': 'chr1', 'POS': 3000, 'ID': 'inv1',
            'ALT': '<INV>', 'FILTER': 'PASS',
            'INFO': 'SVTYPE=INV;SVLEN=200;END=3200',
            'SVTYPE': 'INV', 'VAF': '*',
        }])
        df_1, df_2 = _get_paired_info_dfs(info_df)
        assert df_2.iloc[0]['POS'] == 3200


class TestGetPairedInfoDfsMixed:
    """Test a mixed dataframe with INS and non-INS records together."""

    def test_mixed_sv_types(self):
        """INS should get POS+1 while DEL gets POS+abs(SVLEN) in the same df."""
        info_df = _make_info_df([
            {
                '#CHROM': 'chr1', 'POS': 1000, 'ID': 'del1',
                'ALT': '<DEL>', 'FILTER': 'PASS',
                'INFO': 'SVTYPE=DEL;SVLEN=-500;END=1500',
                'SVTYPE': 'DEL', 'VAF': '*',
            },
            {
                '#CHROM': 'chr1', 'POS': 2000, 'ID': 'ins1',
                'ALT': '<INS>', 'FILTER': 'PASS',
                'INFO': 'SVTYPE=INS;SVLEN=800;END=2001',
                'SVTYPE': 'INS', 'VAF': '*',
            },
        ])
        df_1, df_2 = _get_paired_info_dfs(info_df)
        del_row = df_2[df_2['ID'] == 'del1'].iloc[0]
        ins_row = df_2[df_2['ID'] == 'ins1'].iloc[0]
        assert del_row['POS'] == 1500
        assert ins_row['POS'] == 2001
