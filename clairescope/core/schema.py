"""CLAIREscope schema standardizer, column resolvers, and gene mappers."""
import re
from typing import List, Tuple, Dict, Optional, Any
import pandas as pd

def get_annotation_columns(adata) -> List[str]:
    """Identify categorical annotation columns suitable for cell clustering."""
    candidates = []
    preferred = [
        "cell_type", "cell_states", "cell_state", "cell_type_annotated", "cell_state_annotated",
        "annot_res_0.5", "annot_res_0.7", "annot_res_0.8", "leiden_kc_res_0.5", "leiden_fb_res_0.5",
        "cluster.name", "cluster", "seurat_clusters", "leiden", "louvain", "annotation", "Major_Cell_Type"
    ]
    for col in preferred:
        if col in adata.obs.columns and col not in candidates:
            candidates.append(col)
    for col in adata.obs.columns:
        if col not in candidates:
            if isinstance(adata.obs[col].dtype, pd.CategoricalDtype) or (adata.obs[col].dtype == object and adata.obs[col].nunique() < 100):
                candidates.append(col)
    return candidates if candidates else list(adata.obs.columns)

def get_sample_column(adata) -> str:
    """Identify the sample identifier column."""
    for col in ["sample", "Sample", "wound_type", "orig.ident", "condition", "Condition", "batch", "donor", "tissue"]:
        if col in adata.obs.columns:
            return col
    return adata.obs.columns[0]

def get_gene_display_mappings(var_df: pd.DataFrame, var_names: Any) -> Tuple[List[str], Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Build bidirectional mappings between gene symbols, Ensembl IDs, and display labels."""
    symbol_cols = ["gene_name", "gene_symbols", "symbol", "feature_name", "symbols", "Gene", "gene"]
    id_cols = ["gene_id", "gene_ids", "ensembl_id", "id", "feature_id"]
    
    sym_col = next((c for c in symbol_cols if c in var_df.columns), None)
    id_col = next((c for c in id_cols if c in var_df.columns), None)
    
    display_options = ["None"]
    display_to_var = {}
    sym_to_display = {}
    var_to_display = {}
    
    records = []
    for v in var_names:
        v_str = str(v)
        if v in var_df.index:
            row = var_df.loc[v]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
        else:
            row = None
        
        # Determine symbol
        if sym_col and row is not None and pd.notna(row[sym_col]):
            sym = str(row[sym_col]).strip()
        else:
            sym = v_str.strip()
            
        # Determine ID
        if id_col and row is not None and pd.notna(row[id_col]):
            gid = str(row[id_col]).strip()
        elif v_str.startswith(('ENSG', 'ENSMUSG', 'ENS')):
            gid = v_str.strip()
        else:
            gid = ""
            
        if gid and gid != sym and gid.lower() != "nan" and gid.lower() != "none":
            disp = f"{sym} ({gid})"
        else:
            disp = sym
            
        records.append((disp, sym, gid, v_str))
        
    records.sort(key=lambda x: x[1].upper())
    
    for disp, sym, gid, v_str in records:
        display_options.append(disp)
        display_to_var[disp] = v_str
        var_to_display[v_str] = disp
        
        # Symbol mappings (case-insensitive and exact)
        sym_to_display[sym.upper()] = disp
        sym_to_display[sym] = disp
        
        # ID mappings
        if gid:
            sym_to_display[gid.upper()] = disp
            sym_to_display[gid] = disp
            
    return display_options, display_to_var, sym_to_display, var_to_display

def resolve_gene_var_name(adata, gene_name: str, sym_to_display: Dict[str, str], display_to_var: Dict[str, str]) -> Optional[str]:
    """Resolve a case-insensitive user input string to the exact AnnData var_name index."""
    if not gene_name or gene_name == "None" or not hasattr(adata, "var_names"):
        return None
    q = str(gene_name).strip()
    if not q:
        return None
    
    # 1. Direct match in display_to_var
    if q in display_to_var:
        res = display_to_var[q]
        if res in adata.var_names:
            return res
            
    # 2. Direct match in adata.var_names
    if q in adata.var_names:
        return q
        
    # 3. If string is formatted as "SYMBOL (ID)", parse parts
    if "(" in q and q.endswith(")"):
        parts = q.rsplit(" (", 1)
        clean_sym = parts[0].strip()
        clean_id = parts[1][:-1].strip()
        for cand in [clean_sym, clean_id, clean_sym.upper(), clean_id.upper()]:
            if cand in adata.var_names:
                return cand
            if cand in sym_to_display:
                disp = sym_to_display[cand]
                res = display_to_var.get(disp, None)
                if res is not None and res in adata.var_names:
                    return res
                    
    # 4. Symbol lookup (exact and uppercase)
    if q in sym_to_display:
        disp = sym_to_display[q]
        res = display_to_var.get(disp, None)
        if res is not None and res in adata.var_names:
            return res
            
    q_upper = q.upper()
    if q_upper in sym_to_display:
        disp = sym_to_display[q_upper]
        res = display_to_var.get(disp, None)
        if res is not None and res in adata.var_names:
            return res
            
    # 5. Check if case-insensitive match exists directly in adata.var_names
    # (Cached or linear scan for small/fallback cases)
    if hasattr(adata, "var") and "gene_symbols" in adata.var.columns:
        matches = adata.var[adata.var["gene_symbols"].astype(str).str.upper() == q_upper]
        if not matches.empty:
            return matches.index[0]
            
    return None
