"""Unit tests for schema mapping and column resolution."""
import pytest
import pandas as pd
from clairescope.core.schema import get_gene_display_mappings, resolve_gene_var_name

def test_gene_display_mappings():
    var_df = pd.DataFrame({
        "gene_name": ["CDH1", "COL17A1", "KRT14"],
        "gene_id": ["ENSG00000039068", "ENSG00000065618", "ENSG00000186847"]
    }, index=["CDH1_idx", "COL17A1_idx", "KRT14_idx"])
    
    options, disp_to_var, sym_to_disp, var_to_disp = get_gene_display_mappings(var_df, list(var_df.index))
    assert len(options) == 3
    assert "CDH1 (ENSG00000039068)" in options
    assert disp_to_var["CDH1 (ENSG00000039068)"] == "CDH1_idx"
    assert sym_to_disp["CDH1"] == "CDH1 (ENSG00000039068)"
