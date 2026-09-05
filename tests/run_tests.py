import os
import sys
import unittest
import tempfile
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from clairescope.config import (
    load_projects_config,
    load_settings_config,
    load_signatures_config,
    load_pathways_config,
    load_markers_config,
    load_css_styles,
    get_config_file_path,
    save_user_project_config,
    get_next_new_project_name,
    scan_project_datasets,
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

    def test_get_next_new_project_name(self):
        dummy_projects = {
            "PROJ_1": {"name": "New Project 1"},
            "PROJ_2": {"name": "New Project 2"},
        }
        next_name, next_key = get_next_new_project_name(dummy_projects)
        self.assertEqual(next_name, "New Project 3")
        self.assertEqual(next_key, "PROJ_NEW_003")

    def test_scan_project_datasets(self):
        temp_dir = tempfile.mkdtemp()
        try:
            root_h5ad = os.path.join(temp_dir, "dataset_root.h5ad")
            with open(root_h5ad, "w") as f:
                f.write("mock")
                
            sub_dir = os.path.join(temp_dir, "out", "2026-09-05_test")
            os.makedirs(sub_dir, exist_ok=True)
            sub_h5ad = os.path.join(sub_dir, "dataset_sub.h5ad")
            with open(sub_h5ad, "w") as f:
                f.write("mock")

            active, all_ds = scan_project_datasets(temp_dir, scan_subdirs=["out/2026-09-05_test"])
            self.assertEqual(len(all_ds), 2)
            self.assertTrue(any("dataset_root" in k for k in all_ds.keys()))
            self.assertTrue(any("dataset_sub" in k for k in all_ds.keys()))
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

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

class TestGUI(unittest.TestCase):
    def test_gui_module_import(self):
        from clairescope.gui import ServerControllerGUI, is_port_in_use
        self.assertTrue(callable(is_port_in_use))

if __name__ == "__main__":
    unittest.main(verbosity=2)
