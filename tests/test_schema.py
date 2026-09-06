"""Unit tests for schema mapping and column resolution."""
import pytest
import pandas as pd
import anndata as ad
import numpy as np
from clairescope.core.schema import get_gene_display_mappings, resolve_gene_var_name, get_annotation_columns, get_sample_column

def test_gene_display_mappings_human_ensembl_index():
    var_df = pd.DataFrame({
        "gene_symbols": ["CDH1", "COL17A1", "KRT14"],
        "gene_ids": ["ENSG00000039068", "ENSG00000065618", "ENSG00000186847"]
    }, index=["ENSG00000039068", "ENSG00000065618", "ENSG00000186847"])
    
    options, disp_to_var, sym_to_disp, var_to_disp = get_gene_display_mappings(var_df, list(var_df.index))
    assert "None" in options
    assert "CDH1 (ENSG00000039068)" in options
    assert "COL17A1 (ENSG00000065618)" in options
    assert disp_to_var["CDH1 (ENSG00000039068)"] == "ENSG00000039068"
    assert sym_to_disp["CDH1"] == "CDH1 (ENSG00000039068)"
    assert sym_to_disp["ENSG00000039068"] == "CDH1 (ENSG00000039068)"

    # Test resolve_gene_var_name with AnnData
    adata = ad.AnnData(X=np.zeros((5, 3)), var=var_df)
    assert resolve_gene_var_name(adata, "CDH1", sym_to_disp, disp_to_var) == "ENSG00000039068"
    assert resolve_gene_var_name(adata, "cdh1", sym_to_disp, disp_to_var) == "ENSG00000039068"
    assert resolve_gene_var_name(adata, "ENSG00000039068", sym_to_disp, disp_to_var) == "ENSG00000039068"
    assert resolve_gene_var_name(adata, "CDH1 (ENSG00000039068)", sym_to_disp, disp_to_var) == "ENSG00000039068"
    assert resolve_gene_var_name(adata, "NONEXISTENT", sym_to_disp, disp_to_var) is None

def test_gene_display_mappings_mouse_symbols_index():
    # Project_B Lab B structure: var_names are mouse symbols, gene_ids are ENSMUSG
    var_df = pd.DataFrame({
        "gene_ids": ["ENSMUSG00000025064", "ENSMUSG00000025900", "ENSMUSG00000051951"]
    }, index=["Col17a1", "Rp1", "Xkr4"])
    
    options, disp_to_var, sym_to_disp, var_to_disp = get_gene_display_mappings(var_df, list(var_df.index))
    assert "Col17a1 (ENSMUSG00000025064)" in options
    assert disp_to_var["Col17a1 (ENSMUSG00000025064)"] == "Col17a1"
    assert sym_to_disp["COL17A1"] == "Col17a1 (ENSMUSG00000025064)"
    assert sym_to_disp["Col17a1"] == "Col17a1 (ENSMUSG00000025064)"
    assert sym_to_disp["ENSMUSG00000025064"] == "Col17a1 (ENSMUSG00000025064)"

    adata = ad.AnnData(X=np.zeros((5, 3)), var=var_df)
    assert resolve_gene_var_name(adata, "Col17a1", sym_to_disp, disp_to_var) == "Col17a1"
    assert resolve_gene_var_name(adata, "COL17A1", sym_to_disp, disp_to_var) == "Col17a1"
    assert resolve_gene_var_name(adata, "ENSMUSG00000025064", sym_to_disp, disp_to_var) == "Col17a1"
    assert resolve_gene_var_name(adata, "Col17a1 (ENSMUSG00000025064)", sym_to_disp, disp_to_var) == "Col17a1"
    # Unmatched Human Ensembl ID safely returns None without throwing KeyError
    assert resolve_gene_var_name(adata, "ENSG00000144749", sym_to_disp, disp_to_var) is None
