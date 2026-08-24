import os
import io
import zipfile
import datetime
import streamlit as st
import scanpy as sc
import pandas as pd
import yaml
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.sparse
from scipy.stats import mannwhitneyu, spearmanr, pearsonr
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# Global Plotly Typography Configuration
if "plotly_white" in pio.templates:
    pio.templates["plotly_white"].layout.font.family = "Segoe UI, Arial, sans-serif"
    pio.templates["plotly_white"].layout.font.size = 15
    pio.templates["plotly_white"].layout.title.font.size = 18
    pio.templates["plotly_white"].layout.legend.font.size = 14
    pio.templates["plotly_white"].layout.legend.title.font.size = 15
    pio.templates["plotly_white"].layout.xaxis.tickfont.size = 14
    pio.templates["plotly_white"].layout.yaxis.tickfont.size = 14
    pio.templates["plotly_white"].layout.xaxis.title.font.size = 16
    pio.templates["plotly_white"].layout.yaxis.title.font.size = 16
pio.templates.default = "plotly_white"

# Global Matplotlib Typography Configuration
matplotlib.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'font.sans-serif': ['Segoe UI', 'DejaVu Sans', 'Arial', 'Helvetica']
})

import streamlit.components.v1 as components

# Setup paths
APP_DIR = os.path.dirname(os.path.abspath(__file__))

# Custom Draggable Multiselect Component
COMPONENT_DIR = os.path.join(APP_DIR, "components", "draggable_multiselect")
if os.path.exists(COMPONENT_DIR):
    _draggable_multiselect_comp = components.declare_component("draggable_multiselect", path=COMPONENT_DIR)
else:
    _draggable_multiselect_comp = None

def draggable_multiselect(label, options, default=None, key=None):
    if default is None:
        default = []
    if _draggable_multiselect_comp is not None:
        val = _draggable_multiselect_comp(label=label, options=options, default=default, key=key)
        if val is None:
            return list(default)
        return list(val)
    else:
        return st.multiselect(label, options=options, default=default, key=key)

st.set_page_config(
    page_title="CLAIREscope | Single Cell Analysis Viewer",
    page_icon="🔬",
    layout="wide"
)

# -------------------------------------------------------------
# MULTI-PROJECT REGISTRY & DYNAMIC ENVIRONMENT RESOLUTION
# -------------------------------------------------------------
def get_platform_path(win_path, wsl_path):
    if os.path.exists(win_path):
        return win_path
    if os.path.exists(wsl_path):
        return wsl_path
    return win_path if os.name == 'nt' else wsl_path

PROJECT_REGISTRY = {
    "PROJ_001_HUMAN_EPIDERMAL": {
        "id": "PROJ_001",
        "name": "PROJ_001: Human Epidermal Single Cell",
        "desc": "Epidermal Differentiation Dynamics snRNA-seq & Single-Cell Epidermal Trajectory Analysis (Control, Phenotype_A, Phenotype_B, Rescued)",
        "win_base": r"G:\Data\SingleCell\PROJ_001_HUMAN_EPIDERMAL",
        "wsl_base": "/mnt/data/SingleCell/PROJ_001_HUMAN_EPIDERMAL",
        "scan_subdirs": [
            os.path.join("out", "2026-08-23_human_epidermal_subclustering"),
            os.path.join("out", "2026-07-06_adherens_junction_and_col17a1_correlation_analysis"),
            os.path.join("out", "2026-02-09_celltypist"),
            os.path.join("out", "2026-02-26_harmony"),
            os.path.join("data", "2025-03-18_Public_SingleCell")
        ],
        "default_preload": "adata_kc_norm_cell_typed.h5ad",
        "canonical_samples": ["Control", "Rescued_1", "Rescued_2", "Mutant", "Control_P4", "Sample_Rescued_1", "Sample_Rescued_2", "Sample_Mutant", "Normal", "JEB", "Revertant"],
        "sample_colors": {
            "Control": "#e74c3c", "Rescued_1": "#8e44ad", "Rescued_2": "#f1c40f", "Mutant": "#00a8ff",
            "Control_P4": "#e74c3c", "Sample_Rescued_1": "#8e44ad", "Sample_Rescued_2": "#f1c40f", "Sample_Mutant": "#00a8ff",
            "Normal": "#2ecc71", "JEB": "#e74c3c", "Revertant": "#3498db"
        },
        "default_pairs": [("Control", "Mutant"), ("Rescued_2", "Mutant"), ("Rescued_1", "Mutant")],
        "default_signatures": {
            "Adherens Junction Complex": ["CDH1", "CTNNB1", "CTNNA1", "CTNND1", "JUP"],
            "Desmosomes": ["DSP", "PKP1", "PKP3", "DSG1", "DSG3", "DSC1", "DSC3", "PPL", "EVPL"],
            "Hemidesmosome & Basement Membrane": ["COL17A1", "ITGB4", "ITGA6", "LAMA3", "LAMB3", "LAMC2", "DST"],
            "Cell Cycle / Proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA"]
        }
    },
    "PROJ_002_MURINE_WOUND_HEALING": {
        "id": "PROJ_002",
        "name": "PROJ_002: Murine Wound Healing",
        "desc": "Murine Cutaneous Wound Healing Single-Cell Multiome Dynamics (Full, Full_adj, SB, SB_adj)",
        "win_base": r"G:\Data\SingleCell\PROJ_002_MURINE_WOUND_HEALING",
        "wsl_base": "/mnt/data/SingleCell/PROJ_002_MURINE_WOUND_HEALING",
        "scan_subdirs": [
            os.path.join("out", "2026-08-19_subclustering"),
            os.path.join("out", "2026-08-20_expression_viewer"),
            "data"
        ],
        "default_preload": "adata_harmony.h5ad",
        "canonical_samples": ["Full_adj", "Full", "SB_adj", "SB"],
        "sample_colors": {"Full": "#c0392b", "Full_adj": "#e8a598", "SB": "#2465a8", "SB_adj": "#9cc3e0"},
        "default_pairs": [("Full_adj", "Full"), ("SB_adj", "SB")],
        "default_signatures": {
            "Adherens Junction Complex": ["Cdh1", "Ctnnb1", "Ctnna1", "Ctnnd1", "Dsp", "Jup", "Ppl", "Evpl", "Pkp1", "Pkp3", "Cdh2", "Ctnna2"],
            "Desmosomes": ["Dsp", "Pkp1", "Pkp3", "Dsg1a", "Dsg3", "Dsc1", "Dsc3", "Ppl", "Evpl"],
            "Hemidesmosome & Basement Membrane": ["Col17a1", "Itgb4", "Itga6", "Lama3", "Lamb3", "Lamc2", "Dst"],
            "Cell Cycle / Proliferation": ["Mki67", "Top2a", "Ccnb1", "Cdk1", "Pcna"],
            "Wound Activation": ["Krt6a", "Krt16", "Krt17", "Itgb6", "Sprr1b"]
        }
    }
}

# -------------------------------------------------------------
# SIDEBAR: CLAIREscope BRANDING & PROJECT SELECTOR
# -------------------------------------------------------------
with st.sidebar:
    st.markdown("""
    <div style="padding: 4px 0 12px 0;">
        <h1 style="font-size: 26px; font-weight: 800; color: #E11D48; margin: 0; padding: 0; letter-spacing: -0.5px;">
            🔬 CLAIREscope
        </h1>
        <div style="font-size: 11.5px; font-weight: 600; color: #475569; line-height: 1.25; margin-top: 2px;">
            Cellular Landscape Analysis, Interpretation & Results Explorer
        </div>
        <div style="display: inline-block; background-color: #FFE4E6; color: #BE123C; padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: 700; margin-top: 6px;">
            Single Cell Analysis Viewer
        </div>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("---")
    
    # Project Selector
    project_keys = list(PROJECT_REGISTRY.keys())
    current_proj_idx = 0
    if "selected_project_key" in st.session_state and st.session_state["selected_project_key"] in project_keys:
        current_proj_idx = project_keys.index(st.session_state["selected_project_key"])
        
    selected_project_key = st.selectbox(
        "📂 Select Active Project:",
        options=project_keys,
        index=current_proj_idx,
        format_func=lambda k: PROJECT_REGISTRY[k]["name"],
        key="project_picker"
    )
    
    if st.session_state.get("selected_project_key") != selected_project_key:
        st.session_state["selected_project_key"] = selected_project_key
        st.session_state.pop("selected_dataset_file", None)
        st.session_state.pop("selected_gene", None)
        st.rerun()

curr_proj = PROJECT_REGISTRY[selected_project_key]
PROJ_BASE = get_platform_path(curr_proj["win_base"], curr_proj["wsl_base"])
SCAN_DIRS = [os.path.join(PROJ_BASE, sub) for sub in curr_proj["scan_subdirs"]]
if not any(os.path.exists(d) for d in SCAN_DIRS):
    SCAN_DIRS = [PROJ_BASE]

canonical_samples = curr_proj["canonical_samples"]
sample_color_map = curr_proj["sample_colors"]
default_pairs = curr_proj["default_pairs"]
DEFAULT_SIGNATURES = curr_proj["default_signatures"]

YAML_PATH = os.path.join(APP_DIR, f"{curr_proj['id'].lower()}_cell_type_markers.yaml")
if not os.path.exists(YAML_PATH):
    YAML_PATH = os.path.join(APP_DIR, "cell_type_markers.yaml")

# -------------------------------------------------------------
# DYNAMIC DATASET SCANNER & LOADER
# -------------------------------------------------------------
@st.cache_data
def scan_for_datasets(scan_dirs):
    found = {}
    for s_dir in scan_dirs:
        if os.path.exists(s_dir):
            for root, _, files in os.walk(s_dir):
                for f in files:
                    if f.endswith('.h5ad') and not f.startswith('.'):
                        full_path = os.path.join(root, f)
                        rel_dir = os.path.basename(root)
                        label = f"{f} ({rel_dir})" if rel_dir else f
                        if label not in found:
                            found[label] = full_path
    return found

dataset_dict = scan_for_datasets(SCAN_DIRS)
if not dataset_dict:
    st.error(f"No `.h5ad` datasets found in {PROJ_BASE}. Please verify data directory paths.")
    st.stop()

# Helper for gene metadata & indexing
def prepare_gene_metadata(ad):
    symbols = []
    ids = []
    sym_to_var = {}
    id_to_var = {}
    sym_to_disp = {}
    var_to_disp = {}
    var_to_sym = {}
    
    col_sym = None
    col_id = None
    for c in ad.var.columns:
        c_low = c.lower()
        if c_low in ['gene_symbols', 'genesymbol', 'symbol', 'gene_name', 'genename', 'features']:
            col_sym = c
        if c_low in ['gene_ids', 'geneid', 'id', 'gene_id', 'ensembl_id', 'ensembl']:
            col_id = c
            
    for idx in ad.var_names:
        idx_str = str(idx)
        sym = None
        gid = None
        if col_sym and pd.notna(ad.var.loc[idx, col_sym]):
            sym = str(ad.var.loc[idx, col_sym]).strip()
        if col_id and pd.notna(ad.var.loc[idx, col_id]):
            gid = str(ad.var.loc[idx, col_id]).strip()
            
        if not sym and not idx_str.startswith('ENSG') and not idx_str.startswith('ENSMUSG'):
            sym = idx_str
        elif not sym and gid:
            sym = gid
        elif not sym:
            sym = idx_str
            
        if not gid and (idx_str.startswith('ENSG') or idx_str.startswith('ENSMUSG')):
            gid = idx_str
            
        if sym:
            symbols.append(sym)
            sym_to_var[sym] = idx_str
            sym_to_var[sym.upper()] = idx_str
            sym_to_var[sym.lower()] = idx_str
        if gid:
            ids.append(gid)
            id_to_var[gid] = idx_str
            id_to_var[gid.upper()] = idx_str
            
        if sym and gid and sym != gid:
            disp = f"{sym} ({gid})"
        elif sym:
            disp = sym
        else:
            disp = idx_str
            
        disp_clean = disp.split(" (")[0]
        sym_to_disp[sym] = disp
        sym_to_disp[disp_clean] = disp
        var_to_disp[idx_str] = disp
        var_to_sym[idx_str] = sym if sym else idx_str
        
    ad.uns["_gene_meta"] = {
        "symbols": sorted(list(set(symbols))),
        "ids": sorted(list(set(ids))),
        "sym_to_var": sym_to_var,
        "id_to_var": id_to_var,
        "sym_to_disp": sym_to_disp,
        "var_to_disp": var_to_disp,
        "var_to_sym": var_to_sym
    }

@st.cache_resource
def load_anndata(path):
    ad = sc.read_h5ad(path)
    prepare_gene_metadata(ad)
    return ad

# Pick default dataset
default_idx = 0
preferred_file = curr_proj.get("default_preload", "")
for i, k in enumerate(dataset_dict.keys()):
    if preferred_file in k:
        default_idx = i
        break

with st.sidebar:
    st.markdown("### 🗂️ Active Dataset")
    selected_dataset_label = st.selectbox("Choose Dataset (.h5ad):", list(dataset_dict.keys()), index=default_idx, key="dataset_selector")
    selected_dataset_path = dataset_dict[selected_dataset_label]
    selected_dataset_name = os.path.basename(selected_dataset_path).replace('.h5ad', '')

with st.spinner("Loading single-cell dataset..."):
    adata = load_anndata(selected_dataset_path)

# Auto-detect sample column
sample_col = None
for col in ["sample", "Sample", "orig.ident", "condition", "batch", "donor", "treatment"]:
    if col in adata.obs.columns:
        sample_col = col
        break

if sample_col:
    raw_samples = adata.obs[sample_col].unique().tolist()
    ordered_samples = [s for s in canonical_samples if s in raw_samples] + [s for s in raw_samples if s not in canonical_samples]
else:
    ordered_samples = []

# Auto-detect annotation column
categorical_cols = [c for c in adata.obs.columns if adata.obs[c].dtype.name in ['category', 'object'] and adata.obs[c].nunique() < 50]
default_annot_idx = 0
for i, col in enumerate(categorical_cols):
    if col in ['cell_type', 'celltype', 'cell_state', 'cheng2018_cell_state', 'cluster', 'leiden', 'seurat_clusters']:
        default_annot_idx = i
        break

with st.sidebar:
    selected_col = st.selectbox(
        "Annotation Column:",
        categorical_cols,
        index=default_annot_idx if categorical_cols else 0,
        help="Categorical column to use for cell state, cluster, or subtype coloring."
    ) if categorical_cols else None
    
    st.markdown("---")
    st.markdown("### 🔍 Gene & Marker Search")
    gene_meta = adata.uns.get("_gene_meta", {})
    var_to_disp = gene_meta.get("var_to_disp", {})
    sym_to_var = gene_meta.get("sym_to_var", {})
    sym_to_disp = gene_meta.get("sym_to_disp", {})
    all_display_options = sorted(list(set(var_to_disp.values()))) if var_to_disp else list(adata.var_names)
    
    default_gene_idx = 0
    for g_target in ["COL17A1", "Col17a1", "KRT14", "Krt14", "CDH1", "Cdh1"]:
        if g_target in sym_to_disp:
            target_disp = sym_to_disp[g_target]
            if target_disp in all_display_options:
                default_gene_idx = all_display_options.index(target_disp)
                break
                
    selected_display_name = st.selectbox(
        "Search & Select Gene:",
        all_display_options,
        index=default_gene_idx if all_display_options else 0,
        key="global_gene_picker"
    ) if all_display_options else None

# Helper to resolve gene
def resolve_gene_var(query_str):
    if not query_str:
        return None, None
    g_meta = adata.uns.get("_gene_meta", {})
    s_to_var = g_meta.get("sym_to_var", {})
    i_to_var = g_meta.get("id_to_var", {})
    v_to_disp = g_meta.get("var_to_disp", {})
    
    clean = query_str.split(" (")[0].strip()
    if clean in s_to_var:
        v = s_to_var[clean]
        return v, v_to_disp.get(v, clean)
    if clean.upper() in s_to_var:
        v = s_to_var[clean.upper()]
        return v, v_to_disp.get(v, clean)
    if clean in i_to_var:
        v = i_to_var[clean]
        return v, v_to_disp.get(v, clean)
    if query_str in adata.var_names:
        return query_str, v_to_disp.get(query_str, query_str)
    return None, None

resolved_var_name, resolved_display_name = resolve_gene_var(selected_display_name)

# Helper for cell state colors
def rank_cell_state(name):
    low = str(name).lower()
    if 'basal 1' in low or 'basal1' in low or 'quiescent' in low: return (0, str(name))
    if 'basal 2' in low or 'basal2' in low or 'secretory' in low or 'activated' in low: return (1, str(name))
    if 'spinous' in low or 'suprabasal' in low: return (2, str(name))
    if 'granular' in low or 'terminally' in low: return (3, str(name))
    return (4, str(name))

def get_cluster_color_map(ad, cat_col):
    if not cat_col or cat_col not in ad.obs.columns:
        return {}, []
    if hasattr(ad.obs[cat_col], "cat"):
        cats = list(ad.obs[cat_col].cat.categories)
    else:
        cats = sorted(ad.obs[cat_col].dropna().unique().tolist())
        
    has_basal_kw = any('basal' in str(c).lower() or 'spinous' in str(c).lower() for c in cats)
    if has_basal_kw:
        cats = sorted(cats, key=rank_cell_state)
        
    color_key = f"{cat_col}_colors"
    if color_key in ad.uns and len(ad.uns[color_key]) >= len(cats):
        colors = list(ad.uns[color_key])[:len(cats)]
    else:
        cmap = plt.get_cmap('tab20')
        colors = [matplotlib.colors.to_hex(cmap(i % 20)) for i in range(len(cats))]
    return dict(zip(cats, colors)), cats

def get_expression_vector(ad, var_id):
    if not var_id or var_id not in ad.var_names:
        return np.zeros(ad.n_obs)
    expr = ad[:, var_id].X
    if scipy.sparse.issparse(expr):
        return expr.toarray().flatten()
    return np.array(expr).flatten()

# CSV UMAP Embeddings Generator
@st.cache_data
def get_umap_embeddings_csv(ad_cache_key, s_col, c_col, v_name, v_disp):
    if 'X_umap' not in adata.obsm:
        return None
    umap_xy = adata.obsm['X_umap']
    df_umap = pd.DataFrame({
        "Barcode": adata.obs_names,
        "UMAP_1": np.round(umap_xy[:, 0], 4),
        "UMAP_2": np.round(umap_xy[:, 1], 4)
    })
    if s_col and s_col in adata.obs.columns:
        df_umap["Sample"] = adata.obs[s_col].values
    if c_col and c_col in adata.obs.columns:
        df_umap[c_col] = adata.obs[c_col].values
        
    if v_name and v_name in adata.var_names:
        raw_vals = get_expression_vector(adata, v_name)
        clean_name = v_disp.split(" (")[0] if v_disp else v_name
        df_umap[f"Expression_{clean_name}_Raw"] = np.round(raw_vals, 4)
        df_umap[f"Expression_{clean_name}_Log2"] = np.round(np.log2(raw_vals + 1), 4)
        
    return df_umap.to_csv(index=False).encode('utf-8')

# -------------------------------------------------------------
# MAIN HEADER & ACTIVE PROJECT METADATA
# -------------------------------------------------------------
st.title("🔬 CLAIREscope Single Cell Analysis Viewer")
st.markdown(f"""
<div style="background: linear-gradient(90deg, #FFF1F2 0%, #FFFFFF 100%); border-left: 5px solid #E11D48; padding: 10px 16px; border-radius: 6px; margin-bottom: 14px;">
    <div style="font-size: 16px; font-weight: 700; color: #9F1239; margin-bottom: 2px;">
        Active Project: {curr_proj['name']}
    </div>
    <div style="font-size: 13px; color: #475569;">
        {curr_proj['desc']}
    </div>
    <div style="font-size: 12.5px; color: #334155; margin-top: 6px;">
        Dataset: <code>{selected_dataset_label}</code> | Total Cells: <code>{adata.n_obs:,}</code> | Total Genes: <code>{adata.n_vars:,}</code> | Sample Column: <code>{sample_col if sample_col else 'None'}</code>
    </div>
</div>
""", unsafe_allow_html=True)

# -------------------------------------------------------------
# GLOBAL CONTROLS: COLORMAP, SCALE & CONTRAST
# -------------------------------------------------------------
expr_raw_all = get_expression_vector(adata, resolved_var_name) if resolved_var_name else np.zeros(1)
max_possible = float(np.percentile(expr_raw_all, 99.8)) if resolved_var_name and len(expr_raw_all) > 1 else 10.0
def_vmax_init = float(np.percentile(np.log2(expr_raw_all + 1), 99.0)) if resolved_var_name and len(expr_raw_all) > 1 else 3.5

use_log2 = True
chosen_vmax = def_vmax_init
chosen_scale_label = "Log2(Norm+1)"
cmap_choice = "viridis"
pt_size = 3.5
pt_alpha = 0.85

# -------------------------------------------------------------
# MAIN APP NAVIGATION TABS (INCLUDING BULK DOWNLOAD STUDIO)
# -------------------------------------------------------------
(
    tab_static,
    tab_interactive,
    tab_composition,
    tab_gene_violins,
    tab_score_violins,
    tab_scatter,
    tab_trajectory,
    tab_bulk_download
) = st.tabs([
    "🗺️ Static UMAP",
    "✨ Interactive UMAP",
    "📊 Sample Composition",
    "🎻 Gene Expression Violins",
    "📈 Signature & Pathway Scoring",
    "📉 Correlation & Scatter",
    "🌿 Trajectory Analysis",
    "📦 Bulk Download & Export"
])

# ---------------- TAB 1: STATIC UMAP ----------------
with tab_static:
    with st.expander("🎨 Colormap, Scale & Contrast Controls", expanded=False):
        c_s1, c_s2, c_s3, c_s4 = st.columns([1.2, 1.2, 1.2, 1.2])
        with c_s1:
            scale_mode = st.radio("Expression Scale:", ["Log2(Norm+1)", "Linear (Normalized)"], index=0, key="t1_scale")
            use_log2 = (scale_mode == "Log2(Norm+1)")
        with c_s2:
            cmap_choice = st.selectbox("Colormap Palette:", ["viridis", "inferno", "plasma", "magma", "turbo", "Reds", "YlOrRd"], index=0, key="t1_cmap")
        with c_s3:
            clip_perc = st.slider("Contrast Ceiling (%tile):", min_value=90.0, max_value=100.0, value=99.0, step=0.5, key="t1_perc")
        with c_s4:
            s_vals = np.log2(expr_raw_all + 1) if use_log2 else expr_raw_all
            s_vmax = float(np.percentile(s_vals, clip_perc)) if resolved_var_name and len(s_vals) > 1 else (3.5 if use_log2 else 10.0)
            custom_vmax = st.number_input("Colormap Max (vmax):", min_value=0.01, max_value=10000.0, value=round(s_vmax, 2), step=0.5 if use_log2 else 10.0, key="t1_vmax")
            chosen_vmax = float(custom_vmax)
            chosen_scale_label = "Log2(Norm+1)" if use_log2 else "Linear"

    with st.expander("⚙️ Static Grid Layout Controls", expanded=False):
        c_sg1, c_sg2 = st.columns(2)
        with c_sg1:
            stat_grid_cols = st.selectbox("Grid Columns:", [1, 2, 3, 4, 5, 6], index=2, key="stat_grid_cols")
        with c_sg2:
            stat_grid_rows = st.selectbox("Grid Rows:", ["Auto", 1, 2, 3, 4, 5, 6], index=0, key="stat_grid_rows")

    if resolved_var_name:
        with st.spinner("Generating static UMAP grid..."):
            n_samples = len(ordered_samples) if ordered_samples else 1
            n_cols = int(stat_grid_cols)
            n_rows = (n_samples + n_cols - 1) // n_cols if stat_grid_rows == "Auto" else int(stat_grid_rows)
            
            fig_grid, axes = plt.subplots(n_rows, n_cols, figsize=(4.8 * n_cols, 4.4 * n_rows), dpi=150)
            axes_flat = axes.flatten() if hasattr(axes, 'flatten') else [axes]
            umap_xy = adata.obsm['X_umap']
            
            clean_sym = resolved_display_name.split(" (")[0]
            expr_vals = np.log2(expr_raw_all + 1) if use_log2 else expr_raw_all
            
            for idx, s in enumerate(ordered_samples if ordered_samples else ["All"]):
                if idx < len(axes_flat):
                    ax = axes_flat[idx]
                    if s != "All" and sample_col:
                        mask = (adata.obs[sample_col] == s)
                    else:
                        mask = np.ones(adata.n_obs, dtype=bool)
                    
                    # Background unselected in grey
                    ax.scatter(umap_xy[~mask, 0], umap_xy[~mask, 1], c='#ECEFF1', s=1.2, alpha=0.5, rasterized=True)
                    sc_plot = ax.scatter(umap_xy[mask, 0], umap_xy[mask, 1], c=expr_vals[mask], cmap=cmap_choice, vmin=0, vmax=chosen_vmax, s=2.2, alpha=0.85, rasterized=True)
                    
                    s_color = sample_color_map.get(s, "#333333")
                    ax.set_title(f"{s} (N={np.sum(mask):,})", fontsize=11, fontweight='bold', color=s_color)
                    ax.set_aspect('equal', 'box')
                    ax.set_xticks([])
                    ax.set_yticks([])
                    
            for idx in range(len(ordered_samples if ordered_samples else ["All"]), len(axes_flat)):
                axes_flat[idx].axis('off')
                
            fig_grid.suptitle(f"{clean_sym} [{chosen_scale_label}] (vmax={chosen_vmax:.2f})", fontsize=14, fontweight='bold', y=0.99)
            fig_grid.tight_layout()
            st.pyplot(fig_grid)
            
            svg_grid_buf = io.BytesIO()
            fig_grid.savefig(svg_grid_buf, format="svg", bbox_inches="tight")
            
            c_dl1, c_dl2, _ = st.columns([0.28, 0.28, 0.44])
            with c_dl1:
                csv_data_t1 = get_umap_embeddings_csv(adata.shape, sample_col, selected_col, resolved_var_name, resolved_display_name)
                if csv_data_t1:
                    st.download_button("📥 Download UMAP Embeddings (CSV)", data=csv_data_t1, file_name=f"{selected_dataset_name}_umap_embeddings.csv", mime="text/csv", key="dl_tab1_umap_csv")
            with c_dl2:
                st.download_button("📥 Download Grid Plot as SVG", data=svg_grid_buf.getvalue(), file_name=f"{selected_dataset_name}_{clean_sym}_static_grid.svg", mime="image/svg+xml", key="dl_tab1_grid_svg")
            plt.close(fig_grid)
    else:
        st.info("💡 Select a gene in the sidebar to render the multi-condition static expression grid.")
        if 'X_umap' in adata.obsm:
            with st.spinner("Rendering reference UMAPs..."):
                c_ref1, c_ref2 = st.columns(2)
                umap_xy = adata.obsm['X_umap']
                with c_ref1:
                    fig_s, ax_s = plt.subplots(figsize=(6, 5), dpi=150)
                    if sample_col:
                        for s in ordered_samples:
                            m = adata.obs[sample_col] == s
                            ax_s.scatter(umap_xy[m, 0], umap_xy[m, 1], label=s, color=sample_color_map.get(s, "#7f8c8d"), s=1.8, alpha=0.85, rasterized=True)
                        ax_s.set_title(f"Samples ({sample_col})", fontsize=11, fontweight='bold')
                        ax_s.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=5, fontsize=8, frameon=False)
                    ax_s.set_aspect('equal', 'box')
                    ax_s.set_xticks([])
                    ax_s.set_yticks([])
                    fig_s.tight_layout()
                    st.pyplot(fig_s)
                    svg_s_buf = io.BytesIO()
                    fig_s.savefig(svg_s_buf, format="svg", bbox_inches="tight")
                    plt.close(fig_s)
                with c_ref2:
                    fig_c, ax_c = plt.subplots(figsize=(6, 5), dpi=150)
                    if selected_col:
                        cmap_dict, cats = get_cluster_color_map(adata, selected_col)
                        for cat in cats:
                            m = adata.obs[selected_col] == cat
                            ax_c.scatter(umap_xy[m, 0], umap_xy[m, 1], label=cat, color=cmap_dict.get(cat, "#7f8c8d"), s=1.8, alpha=0.85, rasterized=True)
                        ax_c.set_title(f"Cell States ({selected_col})", fontsize=11, fontweight='bold')
                        ax_c.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=5, fontsize=7.5, frameon=False)
                    ax_c.set_aspect('equal', 'box')
                    ax_c.set_xticks([])
                    ax_c.set_yticks([])
                    fig_c.tight_layout()
                    st.pyplot(fig_c)
                    svg_c_buf = io.BytesIO()
                    fig_c.savefig(svg_c_buf, format="svg", bbox_inches="tight")
                    plt.close(fig_c)
                    
                c_r_dl1, c_r_dl2, c_r_dl3, _ = st.columns([0.28, 0.25, 0.25, 0.22])
                with c_r_dl1:
                    csv_data_t1 = get_umap_embeddings_csv(adata.shape, sample_col, selected_col, resolved_var_name, resolved_display_name)
                    if csv_data_t1:
                        st.download_button("📥 Download UMAP Embeddings (CSV)", data=csv_data_t1, file_name=f"{selected_dataset_name}_umap_embeddings.csv", mime="text/csv", key="dl_tab1_umap_csv_ref")
                with c_r_dl2:
                    st.download_button("📥 Download Sample UMAP (SVG)", data=svg_s_buf.getvalue(), file_name=f"{selected_dataset_name}_sample_umap.svg", mime="image/svg+xml", key="dl_tab1_sample_svg")
                with c_r_dl3:
                    st.download_button("📥 Download Cell State UMAP (SVG)", data=svg_c_buf.getvalue(), file_name=f"{selected_dataset_name}_cell_state_umap.svg", mime="image/svg+xml", key="dl_tab1_cellstate_svg")

# ---------------- TAB 2: INTERACTIVE UMAP ----------------
with tab_interactive:
    if 'X_umap' not in adata.obsm:
        st.warning("This dataset does not contain UMAP coordinates ('X_umap').")
    else:
        with st.expander("Filter & Highlight Controls", expanded=True):
            col_f1, col_f2 = st.columns(2)
            all_samples = ordered_samples if ordered_samples else sorted(list(adata.obs[sample_col].unique())) if sample_col else []
            with col_f1:
                selected_samples = draggable_multiselect("Filter / Highlight Samples:", options=all_samples, default=all_samples, key="plotly_filter_samples")
            with col_f2:
                all_cats = list(get_cluster_color_map(adata, selected_col)[1]) if selected_col else []
                selected_cats = draggable_multiselect(f"Filter / Highlight {selected_col}:", options=all_cats, default=all_cats, key="plotly_filter_categories")
                
        with st.spinner("Rendering interactive Plotly UMAPs..."):
            max_cells = 10000
            if adata.n_obs > max_cells:
                np.random.seed(42)
                sub_idx = np.random.choice(adata.n_obs, max_cells, replace=False)
                adata_sub = adata[sub_idx].copy()
            else:
                adata_sub = adata
                
            umap_sub = adata_sub.obsm['X_umap']
            df_plot = pd.DataFrame({
                "UMAP 1": umap_sub[:, 0],
                "UMAP 2": umap_sub[:, 1],
                "Sample": adata_sub.obs[sample_col].values if sample_col else "All",
                "Cell State": adata_sub.obs[selected_col].values if selected_col else "All"
            })
            if resolved_var_name:
                raw_s = get_expression_vector(adata_sub, resolved_var_name)
                df_plot["Expression"] = np.log2(raw_s + 1) if use_log2 else raw_s
                
            c_p1, c_p2 = st.columns(2)
            with c_p1:
                fig_cats = px.scatter(
                    df_plot, x="UMAP 1", y="UMAP 2", color="Cell State",
                    title=f"Cell States ({selected_col})",
                    color_discrete_map=get_cluster_color_map(adata, selected_col)[0],
                    template="plotly_white"
                )
                fig_cats.update_traces(marker=dict(size=pt_size, opacity=0.85))
                fig_cats.update_xaxes(scaleanchor="x", scaleratio=1)
                st.plotly_chart(fig_cats, use_container_width=True)
                
            with c_p2:
                if resolved_var_name:
                    clean_sym = resolved_display_name.split(" (")[0]
                    fig_expr = px.scatter(
                        df_plot, x="UMAP 1", y="UMAP 2", color="Expression",
                        title=f"{clean_sym} ({chosen_scale_label})",
                        color_continuous_scale=cmap_choice,
                        range_color=[0, chosen_vmax],
                        template="plotly_white"
                    )
                    fig_expr.update_traces(marker=dict(size=pt_size, opacity=0.85))
                    fig_expr.update_xaxes(scaleanchor="x", scaleratio=1)
                    st.plotly_chart(fig_expr, use_container_width=True)
                    
            csv_data_t2 = get_umap_embeddings_csv(adata.shape, sample_col, selected_col, resolved_var_name, resolved_display_name)
            if csv_data_t2:
                st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
                c_t2_dl, _ = st.columns([0.28, 0.72])
                with c_t2_dl:
                    st.download_button("📥 Download UMAP Embeddings (CSV)", data=csv_data_t2, file_name=f"{selected_dataset_name}_umap_embeddings.csv", mime="text/csv", key="dl_tab2_umap_csv")

# ---------------- TAB 3: SAMPLE COMPOSITION ----------------
with tab_composition:
    if not selected_col or not sample_col:
        st.warning("Please configure Sample and Annotation columns in sidebar.")
    else:
        all_dataset_samples = ordered_samples if ordered_samples else list(adata.obs[sample_col].unique())
        selected_comp_samples = draggable_multiselect("Select & Reorder Samples to Display:", options=all_dataset_samples, default=all_dataset_samples, key="comp_samples_filter")
        
        df_comp = adata.obs[[sample_col, selected_col]].dropna().copy()
        df_comp = df_comp[df_comp[sample_col].isin(selected_comp_samples)]
        
        c_cp1, c_cp2 = st.columns(2)
        with c_cp1:
            ct = pd.crosstab(df_comp[sample_col], df_comp[selected_col], normalize='index') * 100
            ct = ct.reindex(selected_comp_samples).dropna(how='all')
            fig_bar = px.bar(
                ct.reset_index().melt(id_vars=sample_col, var_name="Cell State", value_name="Percentage"),
                x=sample_col, y="Percentage", color="Cell State",
                title="Cell State Distribution by Sample (%)",
                color_discrete_map=get_cluster_color_map(adata, selected_col)[0],
                template="plotly_white"
            )
            fig_bar.update_layout(barmode='stack', height=480)
            st.plotly_chart(fig_bar, use_container_width=True)
            
        with c_cp2:
            counts_total = df_comp["Cell State"].value_counts().reset_index()
            counts_total.columns = ["Cell State", "Cell Count"]
            fig_donut = px.pie(
                counts_total, names="Cell State", values="Cell Count",
                hole=0.45, title="Global Population Composition",
                color="Cell State",
                color_discrete_map=get_cluster_color_map(adata, selected_col)[0],
                template="plotly_white"
            )
            fig_donut.update_layout(height=480)
            st.plotly_chart(fig_donut, use_container_width=True)

# ---------------- TAB 4: GENE EXPRESSION VIOLINS ----------------
with tab_gene_violins:
    if not resolved_var_name:
        st.info("💡 Select a gene in the sidebar to view cross-condition expression violin distributions.")
    else:
        clean_sym = resolved_display_name.split(" (")[0]
        raw_vals_v = get_expression_vector(adata, resolved_var_name)
        df_v = pd.DataFrame({
            "Sample": adata.obs[sample_col].values if sample_col else "All",
            "Cell State": adata.obs[selected_col].values if selected_col else "All",
            "Expression": np.log2(raw_vals_v + 1) if use_log2 else raw_vals_v
        })
        all_states = list(get_cluster_color_map(adata, selected_col)[1]) if selected_col else ["All"]
        selected_states_v = draggable_multiselect("Select & Reorder Cell States to Include:", options=all_states, default=all_states, key="v_gene_states")
        
        if selected_states_v:
            plots_to_gen = ["All Cells (Global)"] + selected_states_v
            n_plots = len(plots_to_gen)
            fig_v, axes_v = plt.subplots(1, n_plots, figsize=(4.2 * n_plots, 4.5), dpi=150, sharey=True)
            axes_v_flat = axes_v if isinstance(axes_v, np.ndarray) else [axes_v]
            
            for idx, state in enumerate(plots_to_gen):
                ax = axes_v_flat[idx]
                if state == "All Cells (Global)":
                    sub_df = df_v[df_v["Sample"].isin(ordered_samples)]
                else:
                    sub_df = df_v[(df_v["Cell State"] == state) & (df_v["Sample"].isin(ordered_samples))]
                    
                sns.violinplot(
                    data=sub_df, x="Sample", y="Expression", order=ordered_samples,
                    palette=sample_color_map, ax=ax, inner="quartile", cut=0, linewidth=1.2
                )
                ax.set_title(f"{state}", fontsize=11, fontweight='bold')
                ax.set_xlabel("")
                ax.set_ylabel(f"Expression [{chosen_scale_label}]" if idx == 0 else "")
                ax.tick_params(axis='x', rotation=30)
                
            fig_v.suptitle(f"{clean_sym} Single-Cell Expression Distribution across Conditions", fontsize=13, fontweight='bold', y=1.03)
            fig_v.tight_layout()
            st.pyplot(fig_v)
            
            svg_v_buf = io.BytesIO()
            fig_v.savefig(svg_v_buf, format="svg", bbox_inches="tight")
            st.download_button("📥 Download Violins as SVG", data=svg_v_buf.getvalue(), file_name=f"{selected_dataset_name}_{clean_sym}_violins.svg", mime="image/svg+xml", key="dl_tab4_violins_svg")
            plt.close(fig_v)

# ---------------- TAB 5: SIGNATURE VIOLINS ----------------
with tab_score_violins:
    sig_choice = st.selectbox("Select Pathway / Signature:", list(DEFAULT_SIGNATURES.keys()), index=0, key="sig_picker")
    sig_genes = DEFAULT_SIGNATURES[sig_choice]
    
    found_genes = [resolve_gene_var(g)[0] for g in sig_genes if resolve_gene_var(g)[0] is not None]
    if not found_genes:
        st.warning(f"None of the genes for {sig_choice} were found in this dataset.")
    else:
        # Calculate mean expression score
        score_mat = np.column_stack([get_expression_vector(adata, g) for g in found_genes])
        sig_score = np.mean(np.log2(score_mat + 1), axis=1) if use_log2 else np.mean(score_mat, axis=1)
        
        df_sv = pd.DataFrame({
            "Sample": adata.obs[sample_col].values if sample_col else "All",
            "Cell State": adata.obs[selected_col].values if selected_col else "All",
            "Score": sig_score
        })
        all_states_s = list(get_cluster_color_map(adata, selected_col)[1]) if selected_col else ["All"]
        selected_states_sv = draggable_multiselect("Select & Reorder Cell States for Scoring Violins:", options=all_states_s, default=all_states_s, key="sv_states")
        
        if selected_states_sv:
            plots_to_gen_s = ["All Cells (Global)"] + selected_states_sv
            n_plots_s = len(plots_to_gen_s)
            fig_sv, axes_sv = plt.subplots(1, n_plots_s, figsize=(4.2 * n_plots_s, 4.5), dpi=150, sharey=True)
            axes_sv_flat = axes_sv if isinstance(axes_sv, np.ndarray) else [axes_sv]
            
            for idx, state in enumerate(plots_to_gen_s):
                ax = axes_sv_flat[idx]
                if state == "All Cells (Global)":
                    sub_df = df_sv[df_sv["Sample"].isin(ordered_samples)]
                else:
                    sub_df = df_sv[(df_sv["Cell State"] == state) & (df_sv["Sample"].isin(ordered_samples))]
                    
                sns.violinplot(
                    data=sub_df, x="Sample", y="Score", order=ordered_samples,
                    palette=sample_color_map, ax=ax, inner="quartile", cut=0, linewidth=1.2
                )
                ax.set_title(f"{state}", fontsize=11, fontweight='bold')
                ax.set_xlabel("")
                ax.set_ylabel(f"Mean Score [{chosen_scale_label}]" if idx == 0 else "")
                ax.tick_params(axis='x', rotation=30)
                
            fig_sv.suptitle(f"Pathway Activity: {sig_choice} ({len(found_genes)} Genes)", fontsize=13, fontweight='bold', y=1.03)
            fig_sv.tight_layout()
            st.pyplot(fig_sv)
            
            svg_sv_buf = io.BytesIO()
            fig_sv.savefig(svg_sv_buf, format="svg", bbox_inches="tight")
            st.download_button("📥 Download Scoring Violins as SVG", data=svg_sv_buf.getvalue(), file_name=f"{selected_dataset_name}_{sig_choice.replace(' ', '_')}_violins.svg", mime="image/svg+xml", key="dl_tab5_violins_svg")
            plt.close(fig_sv)

# ---------------- TAB 6: CORRELATION & SCATTER ----------------
with tab_scatter:
    c_sc1, c_sc2 = st.columns(2)
    with c_sc1:
        gene_x_disp = st.selectbox("Select Gene X:", all_display_options, index=0, key="scat_gene_x")
    with c_sc2:
        def_y = 1 if len(all_display_options) > 1 else 0
        gene_y_disp = st.selectbox("Select Gene Y:", all_display_options, index=def_y, key="scat_gene_y")
        
    var_x, disp_x = resolve_gene_var(gene_x_disp)
    var_y, disp_y = resolve_gene_var(gene_y_disp)
    
    if var_x and var_y:
        val_x = get_expression_vector(adata, var_x)
        val_y = get_expression_vector(adata, var_y)
        
        df_sc = pd.DataFrame({
            "Gene X": np.log2(val_x + 1) if use_log2 else val_x,
            "Gene Y": np.log2(val_y + 1) if use_log2 else val_y,
            "Sample": adata.obs[sample_col].values if sample_col else "All",
            "Cell State": adata.obs[selected_col].values if selected_col else "All"
        })
        
        # Subsample for rendering speed
        if len(df_sc) > 8000:
            df_sc = df_sc.sample(8000, random_state=42)
            
        fig_sc = px.scatter(
            df_sc, x="Gene X", y="Gene Y", color="Sample",
            title=f"Co-expression: {disp_x.split(' (')[0]} vs {disp_y.split(' (')[0]}",
            color_discrete_map=sample_color_map,
            labels={"Gene X": f"{disp_x.split(' (')[0]} [{chosen_scale_label}]", "Gene Y": f"{disp_y.split(' (')[0]} [{chosen_scale_label}]"},
            template="plotly_white"
        )
        fig_sc.update_traces(marker=dict(size=4, opacity=0.7))
        st.plotly_chart(fig_sc, use_container_width=True)

# ---------------- TAB 7: TRAJECTORY ANALYSIS ----------------
with tab_trajectory:
    pt_cols = [c for c in adata.obs.columns if 'pseudotime' in c.lower() or 'dpt' in c.lower() or 'latent_time' in c.lower()]
    if not pt_cols:
        st.warning("This dataset does not contain precomputed pseudotime coordinates (e.g. `dpt_pseudotime`).")
    else:
        chosen_pt = st.selectbox("Trajectory Pseudotime Variable:", pt_cols, index=0, key="traj_pt_col")
        traj_genes = ["COL17A1", "Col17a1", "KRT14", "Krt14", "KRT10", "Krt10", "AREG", "Areg"]
        avail_traj_genes = [sym_to_disp[g] for g in traj_genes if g in sym_to_disp]
        
        selected_multi_features = draggable_multiselect(
            "Select & Reorder Features along Pseudotime:",
            options=all_display_options,
            default=avail_traj_genes[:4] if avail_traj_genes else all_display_options[:4],
            key="traj_multi_features"
        )
        
        if selected_multi_features:
            with st.spinner("Calculating trajectory kinetics..."):
                n_feats = len(selected_multi_features)
                fig_t, axes_t = plt.subplots(1, n_feats, figsize=(5.2 * n_feats, 4.4), dpi=150, sharey=False)
                axes_t_flat = axes_t if isinstance(axes_t, np.ndarray) else [axes_t]
                
                pt_vals = adata.obs[chosen_pt].values
                valid_pt = ~np.isnan(pt_vals)
                
                for idx, feat_disp in enumerate(selected_multi_features):
                    ax = axes_t_flat[idx]
                    f_var, f_disp = resolve_gene_var(feat_disp)
                    clean_f = f_disp.split(" (")[0]
                    f_expr = get_expression_vector(adata, f_var)
                    f_scaled = np.log2(f_expr + 1) if use_log2 else f_expr
                    
                    for s in ordered_samples:
                        m = (adata.obs[sample_col] == s) & valid_pt
                        if np.sum(m) > 10:
                            s_pt = pt_vals[m]
                            s_exp = f_scaled[m]
                            sort_idx = np.argsort(s_pt)
                            s_pt_sorted = s_pt[sort_idx]
                            s_exp_sorted = s_exp[sort_idx]
                            
                            # Rolling window smooth
                            win = max(len(s_exp_sorted) // 25, 5)
                            smoothed = pd.Series(s_exp_sorted).rolling(win, center=True, min_periods=1).mean()
                            ax.plot(s_pt_sorted, smoothed, label=s, color=sample_color_map.get(s, "#333333"), linewidth=2.5)
                            
                    ax.set_title(f"{clean_f} along Pseudotime", fontsize=11, fontweight='bold')
                    ax.set_xlabel(chosen_pt, fontsize=9)
                    ax.set_ylabel(f"Expression [{chosen_scale_label}]" if idx == 0 else "", fontsize=9)
                    ax.legend(frameon=False, fontsize=8)
                    
                fig_t.tight_layout()
                st.pyplot(fig_t)
                svg_t_buf = io.BytesIO()
                fig_t.savefig(svg_t_buf, format="svg", bbox_inches="tight")
                st.download_button("📥 Download Trajectory Kinetics as SVG", data=svg_t_buf.getvalue(), file_name=f"{selected_dataset_name}_trajectory_kinetics.svg", mime="image/svg+xml", key="dl_tab7_traj_svg")
                plt.close(fig_t)

# -------------------------------------------------------------
# TAB 8: BULK DOWNLOAD & EXPORT STUDIO
# -------------------------------------------------------------
with tab_bulk_download:
    st.markdown("""
    <div style="background-color: #FEF3C7; border-left: 5px solid #F59E0B; padding: 12px 16px; border-radius: 6px; margin-bottom: 16px;">
        <div style="font-size: 15px; font-weight: 700; color: #92400E;">
            ⚠️ Bulk Generation & Package Export Studio
        </div>
        <div style="font-size: 13px; color: #78350F; margin-top: 2px;">
            Exporting high-resolution multi-panel figures and structured tabular datasets across multiple genes and pathways may take a few moments. All assets will be packaged into a single clean <code>.zip</code> file.
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # 1. Feature Selections
    st.markdown("### 1️⃣ Features to Include in Bulk Export")
    c_b1, c_b2 = st.columns(2)
    def_bulk_genes = ["COL17A1", "Col17a1", "KRT14", "Krt14", "KRT10", "Krt10", "CDH1", "Cdh1", "DSP", "Dsp"]
    valid_bulk_genes = [sym_to_disp[g] for g in def_bulk_genes if g in sym_to_disp]
    
    with c_b1:
        bulk_selected_genes = draggable_multiselect(
            "Genes to Export (One row per gene in table & individual plots):",
            options=all_display_options,
            default=valid_bulk_genes[:6] if valid_bulk_genes else all_display_options[:4],
            key="bulk_genes_input"
        )
    with c_b2:
        bulk_selected_pathways = draggable_multiselect(
            "Pathways / Signatures to Export:",
            options=list(DEFAULT_SIGNATURES.keys()),
            default=list(DEFAULT_SIGNATURES.keys())[:4],
            key="bulk_pathways_input"
        )
        
    st.markdown("---")
    st.markdown("### 2️⃣ Select Figures & Formats to Export")
    c_fig1, c_fig2, c_fig3 = st.columns(3)
    with c_fig1:
        inc_grid_umap = st.checkbox("Multi-Condition UMAP Expression Grids", value=True, key="b_inc_grid")
        inc_ref_umap = st.checkbox("Sample & Cell State Reference UMAPs", value=True, key="b_inc_ref")
    with c_fig2:
        inc_violins = st.checkbox("Gene Expression Statistical Violins", value=True, key="b_inc_violins")
        inc_sig_violins = st.checkbox("Signature & Pathway Score Violins", value=True, key="b_inc_sig_violins")
    with c_fig3:
        inc_comp_plots = st.checkbox("Population Composition (Bars & Donuts)", value=True, key="b_inc_comp")
        inc_traj_plots = st.checkbox("Trajectory Kinetics Curves", value=True if pt_cols else False, key="b_inc_traj")
        
    img_formats = st.multiselect("Image Formats to Generate:", ["SVG", "PNG", "PDF"], default=["SVG", "PNG"], key="bulk_img_formats")
    
    st.markdown("---")
    st.markdown("### 3️⃣ Select Tabular Summary Datasets")
    c_tab1, c_tab2 = st.columns(2)
    with c_tab1:
        inc_table_gene = st.checkbox("📊 Gene Expression Summary (1 row per gene, samples & cell states in columns)", value=True, key="b_tbl_gene")
        inc_table_sig = st.checkbox("📈 Pathway Score Summary (1 row per pathway, samples & cell states in columns)", value=True, key="b_tbl_sig")
    with c_tab2:
        inc_table_comp = st.checkbox("🥧 Cell Type Composition Table (Counts & Percentages)", value=True, key="b_tbl_comp")
        inc_table_umap = st.checkbox("🗺️ Full UMAP Embeddings & Coordinates Table (CSV)", value=True, key="b_tbl_umap")
        
    st.markdown("<div style='margin-top: 10px;'></div>", unsafe_allow_html=True)
    
    # ---------------------------------------------------------
    # BULK GENERATION ENGINE
    # ---------------------------------------------------------
    if st.button("⚡ Generate & Build Bulk Export Package (.ZIP)", type="primary", key="btn_run_bulk_export"):
        zip_buffer = io.BytesIO()
        progress_bar = st.progress(0, text="Initializing export package...")
        
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
            timestamp_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            proj_id = curr_proj["id"]
            
            # -----------------------------------------------------
            # A. TABULAR EXPORT 1: GENE EXPRESSION SUMMARY
            # (One file for all genes, 1 row per gene, sample & cell state attributes in columns)
            # -----------------------------------------------------
            if inc_table_gene and bulk_selected_genes:
                progress_bar.progress(10, text="Building Gene Expression Summary Table...")
                gene_rows = []
                for g_disp in bulk_selected_genes:
                    g_var, g_resolved_disp = resolve_gene_var(g_disp)
                    if g_var:
                        clean_sym = g_resolved_disp.split(" (")[0]
                        g_vals = get_expression_vector(adata, g_var)
                        g_log2 = np.log2(g_vals + 1)
                        
                        row = {
                            "Gene_Symbol": clean_sym,
                            "Gene_ID": g_var,
                            "Global_Mean_Raw": np.round(np.mean(g_vals), 4),
                            "Global_Mean_Log2": np.round(np.mean(g_log2), 4),
                            "Global_Pct_Expressing": np.round(np.mean(g_vals > 0) * 100, 2)
                        }
                        
                        # Per Sample columns
                        if sample_col:
                            for s in ordered_samples:
                                m_s = (adata.obs[sample_col] == s)
                                row[f"Sample_{s}_Mean_Log2"] = np.round(np.mean(g_log2[m_s]), 4) if np.sum(m_s) > 0 else 0.0
                                row[f"Sample_{s}_Pct_Expr"] = np.round(np.mean(g_vals[m_s] > 0) * 100, 2) if np.sum(m_s) > 0 else 0.0
                                
                        # Per Cell State columns
                        if selected_col:
                            _, state_cats = get_cluster_color_map(adata, selected_col)
                            for st_name in state_cats:
                                m_st = (adata.obs[selected_col] == st_name)
                                clean_st = str(st_name).replace(" ", "_")
                                row[f"State_{clean_st}_Mean_Log2"] = np.round(np.mean(g_log2[m_st]), 4) if np.sum(m_st) > 0 else 0.0
                                row[f"State_{clean_st}_Pct_Expr"] = np.round(np.mean(g_vals[m_st] > 0) * 100, 2) if np.sum(m_st) > 0 else 0.0
                                
                        gene_rows.append(row)
                        
                df_gene_summary = pd.DataFrame(gene_rows)
                zip_file.writestr(f"tables/gene_expression_summary.csv", df_gene_summary.to_csv(index=False))
                
            # -----------------------------------------------------
            # B. TABULAR EXPORT 2: PATHWAY / SIGNATURE SCORE SUMMARY
            # (One file for all pathways, 1 row per pathway, samples & states in columns)
            # -----------------------------------------------------
            if inc_table_sig and bulk_selected_pathways:
                progress_bar.progress(25, text="Building Pathway Score Summary Table...")
                pathway_rows = []
                for p_name in bulk_selected_pathways:
                    if p_name in DEFAULT_SIGNATURES:
                        p_genes = DEFAULT_SIGNATURES[p_name]
                        p_vars = [resolve_gene_var(g)[0] for g in p_genes if resolve_gene_var(g)[0] is not None]
                        if p_vars:
                            p_mat = np.column_stack([get_expression_vector(adata, g) for g in p_vars])
                            p_score = np.mean(np.log2(p_mat + 1), axis=1)
                            
                            p_row = {
                                "Pathway_Name": p_name,
                                "Genes_Included_Count": len(p_vars),
                                "Genes_List": "; ".join(p_genes),
                                "Global_Mean_Score": np.round(np.mean(p_score), 4)
                            }
                            
                            if sample_col:
                                for s in ordered_samples:
                                    m_s = (adata.obs[sample_col] == s)
                                    p_row[f"Sample_{s}_Mean_Score"] = np.round(np.mean(p_score[m_s]), 4) if np.sum(m_s) > 0 else 0.0
                                    
                            if selected_col:
                                _, state_cats = get_cluster_color_map(adata, selected_col)
                                for st_name in state_cats:
                                    m_st = (adata.obs[selected_col] == st_name)
                                    clean_st = str(st_name).replace(" ", "_")
                                    p_row[f"State_{clean_st}_Mean_Score"] = np.round(np.mean(p_score[m_st]), 4) if np.sum(m_st) > 0 else 0.0
                                    
                            pathway_rows.append(p_row)
                            
                df_pathway_summary = pd.DataFrame(pathway_rows)
                zip_file.writestr(f"tables/pathway_scores_summary.csv", df_pathway_summary.to_csv(index=False))

            # -----------------------------------------------------
            # C. TABULAR EXPORT 3: SAMPLE COMPOSITION TABLE
            # -----------------------------------------------------
            if inc_table_comp and sample_col and selected_col:
                progress_bar.progress(35, text="Building Sample Composition Table...")
                ct_counts = pd.crosstab(adata.obs[sample_col], adata.obs[selected_col])
                ct_pct = pd.crosstab(adata.obs[sample_col], adata.obs[selected_col], normalize='index') * 100
                
                df_comp_export = ct_counts.copy()
                df_comp_export.columns = [f"{c}_Count" for c in df_comp_export.columns]
                for c in ct_pct.columns:
                    df_comp_export[f"{c}_Percentage"] = np.round(ct_pct[c], 2)
                df_comp_export["Total_Cells"] = adata.obs[sample_col].value_counts()
                zip_file.writestr(f"tables/sample_composition_summary.csv", df_comp_export.to_csv())

            # -----------------------------------------------------
            # D. TABULAR EXPORT 4: UMAP EMBEDDINGS TABLE
            # -----------------------------------------------------
            if inc_table_umap and 'X_umap' in adata.obsm:
                progress_bar.progress(45, text="Exporting UMAP Coordinates & Cell Metadata...")
                umap_csv_bytes = get_umap_embeddings_csv(adata.shape, sample_col, selected_col, resolved_var_name, resolved_display_name)
                if umap_csv_bytes:
                    zip_file.writestr(f"tables/umap_embeddings_and_metadata.csv", umap_csv_bytes)

            # -----------------------------------------------------
            # E. FIGURE EXPORTS
            # -----------------------------------------------------
            def save_fig_to_zip(fig_obj, filename_base, fmts):
                for fmt in fmts:
                    buf = io.BytesIO()
                    dpi_val = 300 if fmt in ['png', 'pdf'] else None
                    fig_obj.savefig(buf, format=fmt.lower(), bbox_inches="tight", dpi=dpi_val)
                    zip_file.writestr(f"figures/{filename_base}.{fmt.lower()}", buf.getvalue())

            # Reference UMAPs
            if inc_ref_umap and 'X_umap' in adata.obsm:
                progress_bar.progress(55, text="Rendering Reference UMAP figures...")
                fig_ref_all, (ax_r1, ax_r2) = plt.subplots(1, 2, figsize=(12, 5), dpi=200)
                umap_xy = adata.obsm['X_umap']
                if sample_col:
                    for s in ordered_samples:
                        m = adata.obs[sample_col] == s
                        ax_r1.scatter(umap_xy[m, 0], umap_xy[m, 1], label=s, color=sample_color_map.get(s, "#7f8c8d"), s=2.0, alpha=0.85)
                    ax_r1.set_title(f"Samples ({sample_col})", fontsize=12, fontweight='bold')
                    ax_r1.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=4, fontsize=8.5, frameon=False)
                ax_r1.set_aspect('equal', 'box')
                ax_r1.set_xticks([])
                ax_r1.set_yticks([])
                
                if selected_col:
                    cmap_dict, cats = get_cluster_color_map(adata, selected_col)
                    for cat in cats:
                        m = adata.obs[selected_col] == cat
                        ax_r2.scatter(umap_xy[m, 0], umap_xy[m, 1], label=cat, color=cmap_dict.get(cat, "#7f8c8d"), s=2.0, alpha=0.85)
                    ax_r2.set_title(f"Cell States ({selected_col})", fontsize=12, fontweight='bold')
                    ax_r2.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=4, fontsize=8, frameon=False)
                ax_r2.set_aspect('equal', 'box')
                ax_r2.set_xticks([])
                ax_r2.set_yticks([])
                fig_ref_all.tight_layout()
                save_fig_to_zip(fig_ref_all, "reference_sample_and_cell_state_umaps", img_formats)
                plt.close(fig_ref_all)

            # Gene UMAP Grids
            if inc_grid_umap and bulk_selected_genes and 'X_umap' in adata.obsm:
                total_g = len(bulk_selected_genes)
                for idx_g, g_disp in enumerate(bulk_selected_genes):
                    progress_bar.progress(60 + int((idx_g / total_g) * 20), text=f"Rendering UMAP grid for {g_disp.split(' (')[0]}...")
                    g_var, g_resolved = resolve_gene_var(g_disp)
                    if g_var:
                        clean_sym = g_resolved.split(" (")[0]
                        g_raw = get_expression_vector(adata, g_var)
                        g_scaled = np.log2(g_raw + 1) if use_log2 else g_raw
                        g_vmax = float(np.percentile(g_scaled, 99.0))
                        
                        n_s = len(ordered_samples) if ordered_samples else 1
                        n_cols = 3 if n_s >= 3 else n_s
                        n_rows = (n_s + n_cols - 1) // n_cols
                        
                        fig_g, axes_g = plt.subplots(n_rows, n_cols, figsize=(4.8 * n_cols, 4.4 * n_rows), dpi=200)
                        axes_g_flat = axes_g.flatten() if hasattr(axes_g, 'flatten') else [axes_g]
                        
                        for idx_s, s in enumerate(ordered_samples if ordered_samples else ["All"]):
                            if idx_s < len(axes_g_flat):
                                ax = axes_g_flat[idx_s]
                                mask = (adata.obs[sample_col] == s) if (sample_col and s != "All") else np.ones(adata.n_obs, dtype=bool)
                                ax.scatter(umap_xy[~mask, 0], umap_xy[~mask, 1], c='#ECEFF1', s=1.2, alpha=0.5)
                                ax.scatter(umap_xy[mask, 0], umap_xy[mask, 1], c=g_scaled[mask], cmap="viridis", vmin=0, vmax=g_vmax, s=2.2, alpha=0.85)
                                ax.set_title(f"{s} (N={np.sum(mask):,})", fontsize=11, fontweight='bold', color=sample_color_map.get(s, "#333333"))
                                ax.set_aspect('equal', 'box')
                                ax.set_xticks([])
                                ax.set_yticks([])
                        for idx_rem in range(len(ordered_samples if ordered_samples else ["All"]), len(axes_g_flat)):
                            axes_g_flat[idx_rem].axis('off')
                        fig_g.suptitle(f"{clean_sym} [{chosen_scale_label}] (vmax={g_vmax:.2f})", fontsize=13, fontweight='bold', y=0.99)
                        fig_g.tight_layout()
                        save_fig_to_zip(fig_g, f"umap_grid_{clean_sym}", img_formats)
                        plt.close(fig_g)

            # Pathway Violins
            if inc_sig_violins and bulk_selected_pathways:
                for p_name in bulk_selected_pathways:
                    if p_name in DEFAULT_SIGNATURES:
                        p_genes = DEFAULT_SIGNATURES[p_name]
                        p_vars = [resolve_gene_var(g)[0] for g in p_genes if resolve_gene_var(g)[0] is not None]
                        if p_vars:
                            p_mat = np.column_stack([get_expression_vector(adata, g) for g in p_vars])
                            p_score = np.mean(np.log2(p_mat + 1), axis=1) if use_log2 else np.mean(p_mat, axis=1)
                            df_p_v = pd.DataFrame({
                                "Sample": adata.obs[sample_col].values if sample_col else "All",
                                "Cell State": adata.obs[selected_col].values if selected_col else "All",
                                "Score": p_score
                            })
                            _, cats = get_cluster_color_map(adata, selected_col)
                            plots_p = ["All Cells (Global)"] + cats
                            fig_pv, axes_pv = plt.subplots(1, len(plots_p), figsize=(4.0 * len(plots_p), 4.2), dpi=200, sharey=True)
                            axes_pv_flat = axes_pv if isinstance(axes_pv, np.ndarray) else [axes_pv]
                            for idx_pv, state in enumerate(plots_p):
                                ax = axes_pv_flat[idx_pv]
                                sub = df_p_v if state == "All Cells (Global)" else df_p_v[df_p_v["Cell State"] == state]
                                sns.violinplot(data=sub[sub["Sample"].isin(ordered_samples)], x="Sample", y="Score", order=ordered_samples, palette=sample_color_map, ax=ax, inner="quartile", cut=0)
                                ax.set_title(state, fontsize=10, fontweight='bold')
                                ax.set_xlabel("")
                                ax.set_ylabel("Score" if idx_pv == 0 else "")
                                ax.tick_params(axis='x', rotation=30)
                            fig_pv.suptitle(f"Pathway Activity: {p_name}", fontsize=12, fontweight='bold', y=1.02)
                            fig_pv.tight_layout()
                            save_fig_to_zip(fig_pv, f"violin_pathway_{p_name.replace(' ', '_')}", img_formats)
                            plt.close(fig_pv)

            progress_bar.progress(100, text="Export package successfully created!")
            
        zip_buffer.seek(0)
        st.success(f"🎉 Bulk export bundle generated successfully! (Archive contains figures in {', '.join(img_formats)} and structured summary tables)")
        st.download_button(
            label=f"📥 Download Complete {curr_proj['id']} Export Bundle (.ZIP)",
            data=zip_buffer.getvalue(),
            file_name=f"CLAIREscope_{curr_proj['id']}_BulkExport_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip",
            key="btn_download_zip_bundle"
        )
