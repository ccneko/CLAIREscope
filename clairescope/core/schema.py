"""CLAIREscope schema standardizer, column resolvers, gene mappers, and color palette engine."""
import re
from typing import List, Tuple, Dict, Optional, Any
import pandas as pd
import numpy as np

SEMANTIC_RULES = [
    (r'ife basal|basal 1|basal stem', '#2B4C7E'),         # Matte Slate Navy
    (r'basal 2|secretory basal', '#8A6D3B'),                 # Matte Ochre / Warm Umber
    (r'basal-spinous|transition', '#D97736'),                # Matte Terracotta / Clay
    (r'basal', '#2B4C7E'),                                       # Matte Slate Navy
    (r'spinous|suprabasal', '#5A8F76'),                      # Matte Sage Green
    (r'granular|terminally differentiated', '#437A75'),       # Matte Sea Pine
    (r'cycling|mitotic|cell cycle', '#935D8C'),          # Matte Dusty Plum / Mauve
    (r'isthmus|junctional zone', '#4A7C8C'),                 # Matte Nordic Blue
    (r'infundibulum|sg-opening', '#5E8B95'),                 # Matte Slate Teal
    (r'hfsc|bulge', '#8D6E63'),                              # Matte Warm Taupe
    (r'wound-activated.*migrating', '#B85D3B'),                  # Matte Burnt Sienna
    (r'wound-activated|hyperproliferative', '#C75D4D'),        # Matte Brick / Rust
    (r'channel', '#73628A'),                                     # Matte Slate Violet
    (r'follicular', '#A67B80'),                                  # Matte Dusty Rose
    (r'schwann', '#6E5D53'),                                     # Matte Bark Brown
    (r't cell|immune', '#7B904B'),                           # Matte Moss Green
    (r'lymphatic|endothelial', '#486581'),                   # Matte Steel Blue
    (r'pericyte|smc', '#5C4A72'),                            # Matte Deep Indigo
    (r'mel1|melanocyte', '#627D98'),                         # Matte Slate
    (r'mel2', '#486581'),                                        # Matte Dark Slate
    (r'wnt1', '#C08552'),                                        # Matte Caramel
    (r'low quality|mitochondrial', '#9E9E9E'),               # Matte Silver
]

PALETTE_POOL = [
    '#2B4C7E', '#D97736', '#5A8F76', '#935D8C', '#437A75', '#C75D4D', '#8D6E63', '#4A7C8C',
    '#8A6D3B', '#73628A', '#A67B80', '#7B904B', '#486581', '#5C4A72', '#627D98', '#C08552',
    '#5E8B95', '#6E5D53', '#A2678A', '#5B8266', '#B85D3B', '#607D8B', '#795548', '#856084',
    '#6B8E23', '#4682B4', '#D2691E', '#9370DB', '#3CB371', '#BC8F8F', '#708090', '#CD853F'
]

def rank_cell_state(cat_str: str) -> Tuple[int, str]:
    """Sort cell states in differentiation order."""
    c_low = str(cat_str).lower()
    if 'ife basal' in c_low or 'basal 1' in c_low or 'basal stem' in c_low:
        return (0, str(cat_str))
    elif 'basal 2' in c_low or 'secretory basal' in c_low:
        return (1, str(cat_str))
    elif 'basal' in c_low:
        return (2, str(cat_str))
    elif 'spinous' in c_low or 'suprabasal' in c_low:
        return (3, str(cat_str))
    elif 'granular' in c_low or 'terminally differentiated' in c_low:
        return (4, str(cat_str))
    elif 'cycling' in c_low or 'mitotic' in c_low:
        return (5, str(cat_str))
    elif 'wound' in c_low or 'migrating' in c_low:
        return (6, str(cat_str))
    elif 'isthmus' in c_low or 'infundibulum' in c_low or 'bulge' in c_low or 'hfsc' in c_low:
        return (7, str(cat_str))
    elif 'low quality' in c_low or 'contamination' in c_low:
        return (9, str(cat_str))
    else:
        return (8, str(cat_str))

def get_cluster_color_map(adata_obj, col_name: str) -> Tuple[Dict[str, str], List[str]]:
    """Generate high-contrast, unique colors for all categories in an annotation column."""
    if not col_name or col_name not in adata_obj.obs.columns:
        return {}, []
    if hasattr(adata_obj.obs[col_name], 'cat'):
        categories = adata_obj.obs[col_name].cat.categories.tolist()
    else:
        categories = sorted(adata_obj.obs[col_name].dropna().unique().tolist())
        
    if any(k in str(c).lower() for c in categories for k in ['basal', 'spinous', 'granular', 'differentiated']):
        categories = sorted(categories, key=rank_cell_state)
        
    color_key = f"{col_name}_colors"
    if col_name == 'predicted_labels':
        pred_map = {}
        for cat in categories:
            c_str = str(cat).lower()
            if 'undiff' in c_str:
                pred_map[cat] = '#b0bec5'
            elif 'diff' in c_str:
                pred_map[cat] = '#f5deb3'
            else:
                pred_map[cat] = '#b0bec5'
        colors = [pred_map[c] for c in categories]
        try:
            adata_obj.uns[color_key] = colors
        except Exception:
            pass
        return pred_map, categories

    # Build semantic and distinct palette
    result = {}
    used_colors = set()
    
    # 1. Semantic rules matching
    for cat in categories:
        c_low = str(cat).lower()
        matched_color = None
        for pattern, color in SEMANTIC_RULES:
            if re.search(pattern, c_low):
                if color not in used_colors:
                    matched_color = color
                    break
        if matched_color:
            result[cat] = matched_color
            used_colors.add(matched_color)
            
    # 2. Pool assignment for unmatched categories
    pool_idx = 0
    for cat in categories:
        if cat not in result:
            while pool_idx < len(PALETTE_POOL) and PALETTE_POOL[pool_idx] in used_colors:
                pool_idx += 1
            if pool_idx < len(PALETTE_POOL):
                color = PALETTE_POOL[pool_idx]
                pool_idx += 1
            else:
                import matplotlib.pyplot as plt
                import matplotlib.colors as mcolors
                cmap = plt.get_cmap('tab20b')
                color = mcolors.to_hex(cmap(len(result) % 20))
            result[cat] = color
            used_colors.add(color)

    if hasattr(adata_obj.obs[col_name], 'cat'):
        raw_cats = list(adata_obj.obs[col_name].cat.categories)
        colors_for_uns = [result.get(c, '#7f8c8d') for c in raw_cats]
    else:
        colors_for_uns = [result.get(c, '#7f8c8d') for c in categories]
    try:
        adata_obj.uns[color_key] = colors_for_uns
    except Exception:
        pass

    return result, categories

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
    if hasattr(adata, "var") and "gene_symbols" in adata.var.columns:
        matches = adata.var[adata.var["gene_symbols"].astype(str).str.upper() == q_upper]
        if not matches.empty:
            return matches.index[0]
            
    return None
