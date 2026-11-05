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
    "breast": (os.path.join(ROOT, "results", "breast"),
               "scout_final_matched.csv", ("p_pos_raw", "p_neg_raw")),
    "lung": (os.path.join(ROOT, "results", "lung"),
             "scout_final_lung.csv", ("p_pos", "p_neg")),
    "breast_s6": (os.path.join(ROOT, "results", "breast_s6"),
                  "scout_final_s6.csv", ("p_pos_raw", "p_neg_raw")),
}


def table(leg):
    directory, name, _ = LEGS[leg]
    return pd.read_csv(os.path.join(directory, name))


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
        for leg, (_, _, (pos_col, neg_col)) in LEGS.items():
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
        with open(os.path.join(ROOT, "results", "breast", "exact_end_summary.json")) as fh:
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


class TestPostClaimSummaries(unittest.TestCase):
    """Every downstream summary reports the selection-adjusted claim set."""

    def load(self, *parts):
        with open(os.path.join(ROOT, *parts)) as fh:
            return json.load(fh)

    def test_breast_coverage(self):
        s = self.load("coverage", "phase2_proseg_summary.json")
        self.assertEqual(s["certified_rows"], 45)
        self.assertEqual(s["certified_inside"], 42)
        self.assertEqual(s["undoubled_bh_rows"], 46)
        self.assertEqual(s["sign_agreement_on_certified"], "45/45")

    def test_breast_proseg_null(self):
        s = self.load("coverage", "phase3_summary.json")
        self.assertEqual(s["certified_rows"], 45)
        self.assertEqual(s["significant_in_certified_direction"], 39)
        self.assertAlmostEqual(s["fraction"], 39 / 45)

    def test_margins_use_claims(self):
        b = self.load("margins", "phase4_summary.json")
        self.assertEqual(b["n_certified"], 45)
        self.assertEqual(b["undoubled_bh_rows"], 46)
        self.assertAlmostEqual(b["frac_ge"]["1.25"], 0.6889, places=4)
        self.assertAlmostEqual(b["frac_ge"]["2.0"], 0.1556, places=4)
        self.assertEqual(int(b["c_star_min"] * 1000), 1048)
        lu = self.load("margins", "lung_phase4_summary.json")
        self.assertEqual(lu["n_certified"], 5)
        self.assertEqual(lu["undoubled_bh_rows"], 7)
        self.assertAlmostEqual(lu["c_star_median"], 20.4, places=1)
        self.assertAlmostEqual(lu["frac_ge"]["3.0"], 0.8)

    def test_lung_coverage(self):
        s = self.load("coverage", "phase5_e2_summary.json")
        self.assertEqual(s["certified_rows"], 5)
        self.assertEqual(s["primary"]["n_inside"], 4)
        self.assertEqual(s["undoubled_bh_inside"], 6)

    def test_joint_sign_changes_are_not_coupling(self):
        s = self.load("results", "coupled_end_summary.json")
        self.assertEqual(s["n_joint_sign_changes_from_coupling"], 0)
        self.assertEqual(s["n_joint_sign_changes_vs_sealed"], 9)
        self.assertEqual(s["n_lo_attained"], 157)
        self.assertEqual(s["n_hi_attained"], 149)

    def test_capped_builders_report_sound_bounds(self):
        for path in (("results", "lung", "capacity_check.json"),
                     ("results", "breast_s6", "capacity_check.json")):
            s = self.load(*path)
            self.assertTrue(s["capped_formula_applied"], path)
            self.assertTrue(s["sound_upper_bound"], path)


class TestAttainmentCounts(unittest.TestCase):
    def test_three_section_totals(self):
        parts = [
            (("results", "breast"), 75),
            (("results", "lung"), 52),
            (("results", "breast_s6"), 30),
        ]
        total = 0
        for directory, expected in parts:
            df = pd.read_csv(os.path.join(ROOT, *directory, "coupled_end.csv"))
            self.assertEqual(len(df), expected, directory)
            total += len(df)
        self.assertEqual(total, 157)


if __name__ == "__main__":
    unittest.main(verbosity=2)
