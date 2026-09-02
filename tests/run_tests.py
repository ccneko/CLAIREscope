import os
import sys
import unittest

# Ensure current directory is in sys.path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from clairescope.config import load_projects_config, load_settings_config, load_signatures_config
from clairescope.stats.hypothesis import run_mann_whitney, get_sig_label, format_sig_value
from clairescope.stats.correlation import compute_bivariate_correlation
from clairescope.stats.enrichment import run_hypergeometric_enrichment
from clairescope.core.schema import get_gene_display_mappings, resolve_gene_var_name
import numpy as np
import pandas as pd

class TestConfig(unittest.TestCase):
    def test_load_projects_config(self):
        projects = load_projects_config()
        self.assertIsInstance(projects, dict)
        self.assertIn("D001_Natsuga_JEB_snRNAseq", projects)
        d001 = projects["D001_Natsuga_JEB_snRNAseq"]
        self.assertEqual(d001["id"], "D001")
        self.assertIn("canonical_samples", d001)

    def test_load_settings_config(self):
        settings = load_settings_config()
        self.assertIsInstance(settings, dict)
        self.assertEqual(settings["ui"]["top_padding"], "3.0rem")

    def test_load_signatures_config(self):
        signatures = load_signatures_config()
        self.assertIn("human", signatures)
        self.assertIn("mouse", signatures)

class TestStats(unittest.TestCase):
    def test_sig_label(self):
        self.assertEqual(get_sig_label(0.00001), "****")
        self.assertEqual(get_sig_label(0.0005), "***")
        self.assertEqual(get_sig_label(0.005), "**")
        self.assertEqual(get_sig_label(0.03), "*")
        self.assertEqual(get_sig_label(0.12), "ns")

    def test_format_sig_value(self):
        self.assertEqual(format_sig_value(0.0000123), "1.23e-05")
        self.assertEqual(format_sig_value(1.234), "1.234")
        self.assertEqual(format_sig_value(1500), "1.50e+03")

    def test_mann_whitney(self):
        g1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        g2 = np.array([10.0, 11.0, 12.0, 13.0, 14.0])
        stat, pval, label = run_mann_whitney(g1, g2)
        self.assertFalse(np.isnan(pval))
        self.assertTrue(pval < 0.05)
        self.assertIn(label, ["*", "**"])

    def test_bivariate_correlation(self):
        x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=float)
        y = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float)
        res = compute_bivariate_correlation(x, y)
        self.assertEqual(res["n_cells"], 10)
        self.assertAlmostEqual(res["pearson_r"], 1.0, places=3)
        self.assertAlmostEqual(res["spearman_rho"], 1.0, places=3)
        self.assertAlmostEqual(res["slope"], 2.0, places=3)

    def test_hypergeometric_enrichment(self):
        bg = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10"]
        query = ["GENE1", "GENE2", "GENE3"]
        pathways = {
            "Test Pathway": ["GENE1", "GENE2", "GENE3", "GENE4"],
            "Other Pathway": ["GENE8", "GENE9", "GENE10"]
        }
        df = run_hypergeometric_enrichment(query, bg, pathways, min_overlap=2)
        self.assertFalse(df.empty)
        self.assertEqual(df.iloc[0]["Pathway"], "Test Pathway")
        self.assertEqual(df.iloc[0]["Overlap"], 3)

class TestSchema(unittest.TestCase):
    def test_gene_display_mappings(self):
        var_df = pd.DataFrame({
            "gene_name": ["CDH1", "COL17A1", "KRT14"],
            "gene_id": ["ENSG00000039068", "ENSG00000065618", "ENSG00000186847"]
        }, index=["CDH1_idx", "COL17A1_idx", "KRT14_idx"])
        options, disp_to_var, sym_to_disp, var_to_disp = get_gene_display_mappings(var_df, list(var_df.index))
        self.assertEqual(len(options), 3)
        self.assertIn("CDH1 (ENSG00000039068)", options)
        self.assertEqual(disp_to_var["CDH1 (ENSG00000039068)"], "CDH1_idx")
        self.assertEqual(sym_to_disp["CDH1"], "CDH1 (ENSG00000039068)")

if __name__ == "__main__":
    unittest.main(verbosity=2)
