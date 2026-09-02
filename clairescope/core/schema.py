"""CLAIREscope schema standardizer, column resolvers, and gene mappers."""
import re
from typing import List, Tuple, Dict, Optional, Any
import pandas as pd

def get_annotation_columns(adata) -> List[str]:
    """Identify categorical annotation columns suitable for cell clustering."""
    candidates = []
    preferred = ["cell_type", "cell_states", "cell_state", "seurat_clusters", "leiden", "louvain", "cluster", "annotation", "Major_Cell_Type"]
    for col in preferred:
        if col in adata.obs.columns:
            candidates.append(col)
    for col in adata.obs.columns:
        if col not in candidates:
            if isinstance(adata.obs[col].dtype, pd.CategoricalDtype) or (adata.obs[col].dtype == object and adata.obs[col].nunique() < 100):
                candidates.append(col)
    return candidates if candidates else list(adata.obs.columns)

def get_sample_column(adata) -> str:
    """Identify the sample identifier column."""
    for col in ["sample", "Sample", "orig.ident", "condition", "Condition", "batch", "donor"]:
        if col in adata.obs.columns:
            return col
    return adata.obs.columns[0]

def get_gene_display_mappings(var_df: pd.DataFrame, var_names: List[str]) -> Tuple[List[str], Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Build bidirectional mappings between gene symbols, Ensembl IDs, and display labels."""
    symbol_cols = ["gene_name", "gene_symbols", "symbol", "feature_name", "symbols", "Gene"]
    id_cols = ["gene_id", "gene_ids", "ensembl_id", "id", "feature_id"]
    
    sym_col = next((c for c in symbol_cols if c in var_df.columns), None)
    id_col = next((c for c in id_cols if c in var_df.columns), None)
    
    display_options = []
    display_to_var = {}
    sym_to_display = {}
    var_to_display = {}
    
    for v in var_names:
        sym = str(var_df.loc[v, sym_col]) if sym_col else str(v)
        gid = str(var_df.loc[v, id_col]) if id_col else ""
        
        if gid and gid != sym and gid != "nan":
            disp = f"{sym} ({gid})"
        else:
            disp = sym
            
        display_options.append(disp)
        display_to_var[disp] = v
        var_to_display[v] = disp
        sym_to_display[sym.upper()] = disp
        if gid:
            sym_to_display[gid.upper()] = disp
            
    return display_options, display_to_var, sym_to_display, var_to_display

def resolve_gene_var_name(adata, gene_name: str, sym_to_display: Dict[str, str], display_to_var: Dict[str, str]) -> Optional[str]:
    """Resolve a case-insensitive user input string to the exact AnnData var_name index."""
    if not gene_name or gene_name == "None":
        return None
    q = gene_name.strip()
    if q in display_to_var:
        return display_to_var[q]
    if q in adata.var_names:
        return q
    q_upper = q.upper()
    if q_upper in sym_to_display:
        disp = sym_to_display[q_upper]
        return display_to_var.get(disp, None)
    return None
