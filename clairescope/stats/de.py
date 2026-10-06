"""Differential expression analysis and subset filtering engine."""
from typing import Tuple, List, Optional, Any, Union, Dict
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad


def get_de_candidate_columns(
    adata: ad.AnnData,
    sample_col: Optional[str] = None,
    selected_col: Optional[str] = None,
    anno_cols: Optional[List[str]] = None
) -> List[str]:
    """Identify candidate columns in AnnData obs suitable for group comparison and subset filtering."""
    candidates = []
    
    # Priority columns
    preferred = [sample_col, selected_col, 'lineage', 'hp_lineage_relaxed', 'state_group', 
                 'cell_type', 'cell_states', 'predicted_labels', 'leiden_r02', 'kmeans_k4']
    for col in preferred:
        if col and col in adata.obs.columns and col not in candidates:
            candidates.append(col)
            
    # Include detected annotation columns
    if anno_cols:
        for col in anno_cols:
            if col and col in adata.obs.columns and col not in candidates:
                candidates.append(col)
                
    # Include other low-to-medium cardinality categorical/object columns
    for col in adata.obs.columns:
        if col not in candidates:
            s = adata.obs[col]
            if isinstance(s.dtype, pd.CategoricalDtype) or (s.dtype == object and 1 < s.nunique() <= 100):
                candidates.append(col)
                
    return candidates if candidates else list(adata.obs.columns)


def filter_anndata_for_de(
    adata: ad.AnnData,
    groupby_col: Optional[str] = None,
    filter_col: Optional[str] = None,
    filter_values: Optional[Union[List[Any], Tuple[Any, ...], Any]] = None
) -> Tuple[ad.AnnData, int, int]:
    """Filter AnnData by column values and clean up unused categorical levels.
    
    Returns:
        Tuple of (adata_sub, n_filtered_cells, n_total_cells)
    """
    n_total = len(adata)
    if not filter_col or filter_col in ("None", "None (All Cells)") or not filter_values:
        adata_sub = adata.copy()
        if groupby_col and groupby_col in adata_sub.obs.columns and hasattr(adata_sub.obs[groupby_col], "cat"):
            adata_sub.obs[groupby_col] = adata_sub.obs[groupby_col].cat.remove_unused_categories()
        return adata_sub, n_total, n_total

    if isinstance(filter_values, (str, int, float)):
        filter_values = [filter_values]
    filter_values = list(filter_values)

    mask = adata.obs[filter_col].isin(filter_values)
    adata_sub = adata[mask].copy()
    n_filtered = len(adata_sub)

    # Clean up categories in filtered subset
    for c in [groupby_col, filter_col]:
        if c and c in adata_sub.obs.columns and hasattr(adata_sub.obs[c], "cat"):
            adata_sub.obs[c] = adata_sub.obs[c].cat.remove_unused_categories()

    return adata_sub, n_filtered, n_total


def run_differential_expression(
    adata: ad.AnnData,
    groupby_col: str,
    target_group: str,
    reference_group: str = "Rest of Cells",
    filter_col: Optional[str] = None,
    filter_values: Optional[Union[List[Any], Tuple[Any, ...], Any]] = None,
    method: str = "wilcoxon"
) -> pd.DataFrame:
    """Compute differential expression (Wilcoxon rank-sum) between target and reference within optional subset."""
    adata_sub, _, _ = filter_anndata_for_de(
        adata,
        groupby_col=groupby_col,
        filter_col=filter_col,
        filter_values=filter_values
    )

    if groupby_col not in adata_sub.obs.columns:
        raise ValueError(f"Column '{groupby_col}' not found in AnnData obs.")

    available_groups = adata_sub.obs[groupby_col].dropna().unique().tolist()
    if target_group not in available_groups:
        raise ValueError(f"Target group '{target_group}' not present in subset for column '{groupby_col}'.")

    if reference_group != "Rest of Cells" and reference_group not in available_groups:
        raise ValueError(f"Reference group '{reference_group}' not present in subset for column '{groupby_col}'.")

    if reference_group == "Rest of Cells":
        if len(available_groups) < 2:
            raise ValueError(f"At least two groups are required to compare '{target_group}' against 'Rest of Cells'.")
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby=groupby_col,
            groups=[target_group],
            reference="rest",
            method=method,
            n_genes=None
        )
    else:
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby=groupby_col,
            groups=[target_group],
            reference=reference_group,
            method=method,
            n_genes=None
        )

    df_res = sc.get.rank_genes_groups_df(adata_sub, group=target_group)
    return df_res
