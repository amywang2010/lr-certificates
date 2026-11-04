"""Release checks that run on the shipped artifacts.

No vendor data is required. The checks compare the released row tables with the
counts reported in the paper, recompute the claim layer from the stored
permutation values, and exercise the multiple-testing routines on fixed inputs.

Run with:  python -m unittest discover -s tests -v
"""
import json
import os
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

import b03_claim_layer as claim  # noqa: E402

LEGS = {
    "breast": (os.path.join(ROOT, "results"), ("p_pos_raw", "p_neg_raw")),
    "lung": (os.path.join(ROOT, "lung"), ("p_pos", "p_neg")),
    "breast_s6": (os.path.join(ROOT, "breast_s6"), ("p_pos_raw", "p_neg_raw")),
}


def table(leg):
    return pd.read_csv(os.path.join(LEGS[leg][0], "scout_final_matched.csv"))


class TestRowTables(unittest.TestCase):
    def test_row_counts(self):
        for leg, n in (("breast", 75), ("lung", 54), ("breast_s6", 30)):
            self.assertEqual(len(table(leg)), n)

    def test_interval_signs(self):
        expected = {
            "breast": (39, 12, 24),
            "lung": (1, 32, 21),
            "breast_s6": (0, 17, 13),
        }
        for leg, (pos, neg, zero) in expected.items():
            df = table(leg)
            self.assertEqual(int((df.T_lo > 0).sum()), pos, leg)
            self.assertEqual(int((df.T_hi < 0).sum()), neg, leg)
            self.assertEqual(
                int(((df.T_lo <= 0) & (df.T_hi >= 0)).sum()), zero, leg)

    def test_claim_counts(self):
        expected = {
            "breast": (36, 9),
            "lung": (1, 4),
            "breast_s6": (0, 0),
        }
        for leg, (pos, neg) in expected.items():
            df = table(leg)
            self.assertEqual(int(df.claim_pos.astype(bool).sum()), pos, leg)
            self.assertEqual(int(df.claim_neg.astype(bool).sum()), neg, leg)

    def test_claims_require_interval_sign(self):
        for leg in LEGS:
            df = table(leg)
            self.assertFalse((df.claim_pos.astype(bool) & ~(df.T_lo > 0)).any(), leg)
            self.assertFalse((df.claim_neg.astype(bool) & ~(df.T_hi < 0)).any(), leg)


class TestClaimLayer(unittest.TestCase):
    def test_reproduces_released_claims(self):
        for leg, (_, (pos_col, neg_col)) in LEGS.items():
            df = table(leg)
            cl, _ = claim.claim_layer(df, pos_col, neg_col)
            self.assertTrue(
                (cl.claim_pos.to_numpy() == df.claim_pos.astype(bool).to_numpy()).all(),
                leg)
            self.assertTrue(
                (cl.claim_neg.to_numpy() == df.claim_neg.astype(bool).to_numpy()).all(),
                leg)

    def test_doubling_is_charged_before_by(self):
        p = np.array([0.002, 0.02, 0.5, 0.9])
        self.assertTrue(np.all(claim.by(p) >= claim.bh(p)))
        self.assertTrue(np.all(claim.by(np.minimum(1.0, 2 * p)) >= claim.by(p)))

    def test_bh_matches_hand_computation(self):
        p = np.array([0.01, 0.02, 0.03, 0.04])
        q = claim.bh(p)
        self.assertAlmostEqual(float(q[0]), 0.04)
        self.assertAlmostEqual(float(q[1]), 0.04)
        self.assertAlmostEqual(float(q[3]), 0.04)

    def test_by_harmonic_factor(self):
        p = np.array([0.01, 0.02, 0.03, 0.04])
        h4 = 1 + 1 / 2 + 1 / 3 + 1 / 4
        self.assertAlmostEqual(float(claim.by(p)[0]), 0.04 * h4)

    def test_tiny_family_claim_flags(self):
        df = pd.DataFrame({
            "pair": ["A->B", "A->B", "C->D", "C->D"],
            "region": [0, 0, 0, 0],
            "T_lo": [1.0, -2.0, 1.0, -2.0],
            "T_hi": [2.0, -1.0, 2.0, -1.0],
            "p_pos": [0.001, 0.9, 0.3, 0.9],
            "p_neg": [0.9, 0.001, 0.9, 0.3],
        })
        cl, _ = claim.claim_layer(df, "p_pos", "p_neg")
        self.assertEqual(int(cl.claim_pos.sum()), 1)
        self.assertEqual(int(cl.claim_neg.sum()), 1)


class TestExactEnds(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(ROOT, "results", "exact_end_summary.json")) as fh:
            self.summary = json.load(fh)

    def test_gene_ends(self):
        self.assertEqual(self.summary["n_gene_combos"], 150)
        self.assertEqual(self.summary["pos_exact_equal_count"], 90)

    def test_largest_shortfall(self):
        self.assertEqual(self.summary["pos_max_gap_count"], 361)

    def test_network_size(self):
        self.assertEqual(int(self.summary["nodes_median"]), 4576)
        self.assertEqual(int(self.summary["arcs_median"]), 7366)

    def test_capacity_moves_no_sign(self):
        self.assertEqual(self.summary["cap_sign_changes"], 0)
        self.assertEqual(self.summary["cap_rows_moved"], 46)


class TestAttainmentCounts(unittest.TestCase):
    def test_three_section_totals(self):
        parts = [
            ("results", 75),
            ("lung", 52),
            ("breast_s6", 30),
        ]
        total = 0
        for directory, expected in parts:
            df = pd.read_csv(os.path.join(ROOT, directory, "coupled_end.csv"))
            self.assertEqual(len(df), expected, directory)
            total += len(df)
        self.assertEqual(total, 157)


if __name__ == "__main__":
    unittest.main(verbosity=2)
