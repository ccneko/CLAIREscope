"""Unit tests for differential expression and subset filtering."""
import unittest
import numpy as np
import pandas as pd
import anndata as ad
from clairescope.stats.de import (
    get_de_candidate_columns,
    filter_anndata_for_de,
    run_differential_expression,
)


class TestDifferentialExpression(unittest.TestCase):
    def setUp(self):
        # Create a synthetic AnnData for testing
        np.random.seed(42)
        n_cells = 60
        n_genes = 20
        
        # 30 Basal, 30 Spinous; 30 Control, 30 Mutant
        lineages = ["Basal"] * 30 + ["Spinous"] * 30
        samples = (["Control"] * 15 + ["Mutant"] * 15) * 2
        
        obs_df = pd.DataFrame({
            "lineage": pd.Categorical(lineages),
            "sample": pd.Categorical(samples),
            "cluster": [f"c{i % 3}" for i in range(n_cells)],
        })
        
        # Gene expression: Gene 0 is upregulated in Mutant Basal
        X = np.random.poisson(lam=2.0, size=(n_cells, n_genes)).astype(float)
        # Boost gene 0 in Basal Mutant
        for idx in range(15, 30):
            X[idx, 0] += 10.0
            
        var_df = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
        self.adata = ad.AnnData(X=X, obs=obs_df, var=var_df)

    def test_get_de_candidate_columns(self):
        cols = get_de_candidate_columns(self.adata, sample_col="sample", selected_col="lineage")
        self.assertIn("sample", cols)
        self.assertIn("lineage", cols)
        self.assertIn("cluster", cols)

    def test_filter_anndata_for_de_no_filter(self):
        sub, n_filtered, n_tot = filter_anndata_for_de(self.adata, groupby_col="sample")
        self.assertEqual(n_filtered, 60)
        self.assertEqual(n_tot, 60)

    def test_filter_anndata_for_de_with_subset(self):
        sub, n_filtered, n_tot = filter_anndata_for_de(
            self.adata, groupby_col="sample", filter_col="lineage", filter_values=["Basal"]
        )
        self.assertEqual(n_filtered, 30)
        self.assertEqual(n_tot, 60)
        self.assertTrue(all(sub.obs["lineage"] == "Basal"))
        # Unused categories should be handled
        self.assertEqual(list(sub.obs["sample"].unique()), ["Control", "Mutant"])

    def test_run_differential_expression_within_lineage(self):
        df_deg = run_differential_expression(
            self.adata,
            groupby_col="sample",
            target_group="Mutant",
            reference_group="Control",
            filter_col="lineage",
            filter_values=["Basal"],
            method="wilcoxon"
        )
        self.assertFalse(df_deg.empty)
        self.assertIn("names", df_deg.columns)
        self.assertIn("logfoldchanges", df_deg.columns)
        self.assertIn("pvals_adj", df_deg.columns)
        # Top upregulated gene in Mutant Basal should be Gene_0
        top_gene = df_deg.iloc[0]["names"]
        self.assertEqual(top_gene, "Gene_0")
        self.assertGreater(df_deg.iloc[0]["logfoldchanges"], 0)


if __name__ == "__main__":
    unittest.main()
