import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from clairescope.config import (
    load_projects_config,
    load_settings_config,
    load_signatures_config,
    load_pathways_config,
    load_markers_config,
    load_css_styles,
    get_config_file_path,
)
from clairescope.stats.hypothesis import run_mann_whitney, get_sig_label, format_sig_value
from clairescope.stats.correlation import compute_bivariate_correlation
from clairescope.stats.enrichment import run_hypergeometric_enrichment
from clairescope.core.schema import get_gene_display_mappings, resolve_gene_var_name
import numpy as np
import pandas as pd

class TestConfig(unittest.TestCase):
    def test_config_resolution(self):
        p = get_config_file_path("settings.yaml")
        self.assertTrue(os.path.exists(p))

    def test_load_projects_config(self):
        projects = load_projects_config()
        self.assertIsInstance(projects, dict)
        self.assertTrue(len(projects) > 0)

    def test_load_settings_config(self):
        settings = load_settings_config()
        self.assertEqual(settings["ui"]["top_padding"], "3.0rem")

    def test_load_signatures_config(self):
        signatures = load_signatures_config()
        self.assertIn("human", signatures)

    def test_load_pathways_config(self):
        pathways = load_pathways_config()
        self.assertIn("Hallmark: Epithelial Mesenchymal Transition", pathways)

    def test_load_markers_config(self):
        markers = load_markers_config()
        self.assertIn("Human", markers)

    def test_load_css_styles(self):
        css = load_css_styles()
        self.assertIn("block-container", css)

class TestStats(unittest.TestCase):
    def test_sig_label(self):
        self.assertEqual(get_sig_label(0.00001), "****")

    def test_format_sig_value(self):
        self.assertEqual(format_sig_value(0.0000123), "1.23e-05")

    def test_mann_whitney(self):
        g1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        g2 = np.array([10.0, 11.0, 12.0, 13.0, 14.0])
        stat, pval, label = run_mann_whitney(g1, g2)
        self.assertTrue(pval < 0.05)

    def test_bivariate_correlation(self):
        x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=float)
        y = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float)
        res = compute_bivariate_correlation(x, y)
        self.assertAlmostEqual(res["pearson_r"], 1.0, places=3)

    def test_hypergeometric_enrichment(self):
        bg = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10"]
        query = ["GENE1", "GENE2", "GENE3"]
        pathways = {"Test Pathway": ["GENE1", "GENE2", "GENE3", "GENE4"]}
        df = run_hypergeometric_enrichment(query, bg, pathways, min_overlap=2)
        self.assertFalse(df.empty)

class TestSchema(unittest.TestCase):
    def test_gene_display_mappings(self):
        var_df = pd.DataFrame({"gene_name": ["CDH1"]}, index=["CDH1_idx"])
        options, disp_to_var, sym_to_disp, var_to_disp = get_gene_display_mappings(var_df, list(var_df.index))
        self.assertIn("CDH1", sym_to_disp)

if __name__ == "__main__":
    unittest.main(verbosity=2)
