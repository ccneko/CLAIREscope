"""Unit tests for statistical hypothesis testing, correlation, and ORA."""
import pytest
import numpy as np
import pandas as pd
from clairescope.stats.hypothesis import run_mann_whitney, get_sig_label, format_sig_value
from clairescope.stats.correlation import compute_bivariate_correlation
from clairescope.stats.enrichment import run_hypergeometric_enrichment

def test_sig_label():
    assert get_sig_label(0.00001) == "****"
    assert get_sig_label(0.0005) == "***"
    assert get_sig_label(0.005) == "**"
    assert get_sig_label(0.03) == "*"
    assert get_sig_label(0.12) == "ns"

def test_format_sig_value():
    assert format_sig_value(0.0000123) == "1.23e-05"
    assert format_sig_value(1.234) == "1.234"
    assert format_sig_value(1500) == "1.50e+03"

def test_mann_whitney():
    g1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    g2 = np.array([10.0, 11.0, 12.0, 13.0, 14.0])
    stat, pval, label = run_mann_whitney(g1, g2)
    assert not np.isnan(pval)
    assert pval < 0.05
    assert label in ["*", "**"]

def test_bivariate_correlation():
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=float)
    y = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float)
    res = compute_bivariate_correlation(x, y)
    assert res["n_cells"] == 10
    assert pytest.approx(res["pearson_r"], 0.001) == 1.0
    assert pytest.approx(res["spearman_rho"], 0.001) == 1.0
    assert pytest.approx(res["slope"], 0.001) == 2.0

def test_hypergeometric_enrichment():
    bg = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10"]
    query = ["GENE1", "GENE2", "GENE3"]
    pathways = {
        "Test Pathway": ["GENE1", "GENE2", "GENE3", "GENE4"],
        "Other Pathway": ["GENE8", "GENE9", "GENE10"]
    }
    df = run_hypergeometric_enrichment(query, bg, pathways, min_overlap=2)
    assert not df.empty
    assert df.iloc[0]["Pathway"] == "Test Pathway"
    assert df.iloc[0]["Overlap"] == 3
