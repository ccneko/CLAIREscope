import os
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

# Setup paths
APP_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_BASE = os.path.dirname(os.path.dirname(APP_DIR))
DATA_DIR = os.path.join(PROJ_BASE, "out", "2026-07-06_adherens_junction_and_col17a1_correlation_analysis")
SCAN_DIRS = [
    os.path.join(PROJ_BASE, "out", "2026-07-06_adherens_junction_and_col17a1_correlation_analysis"),
    os.path.join(PROJ_BASE, "out", "2026-02-09_celltypist"),
    os.path.join(PROJ_BASE, "out", "2026-02-26_harmony"),
    os.path.join(PROJ_BASE, "data", "2025-03-18_scRNAseq")
]
YAML_PATH = os.path.join(APP_DIR, "cell_type_markers.yaml")

st.set_page_config(page_title="D001 JEB Expression Viewer & Marker Editor", layout="wide")

# Global Matplotlib Typography Configuration
matplotlib.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 13
})

# Global Streamlit Typography & Readability CSS
st.markdown("""
<style>
    /* Restore Page Title & Heading Hierarchy */
    h1, div[data-testid="stMarkdownContainer"] h1 {
        font-size: 2.35rem !important;
        font-weight: 700 !important;
        line-height: 1.25 !important;
        margin-bottom: 0.4rem !important;
        color: #1e293b !important;
    }
    h2, div[data-testid="stMarkdownContainer"] h2 {
        font-size: 1.7rem !important;
        font-weight: 700 !important;
        margin-top: 1.2rem !important;
        margin-bottom: 0.5rem !important;
    }
    h3, div[data-testid="stMarkdownContainer"] h3 {
        font-size: 1.35rem !important;
        font-weight: 600 !important;
        margin-top: 1.0rem !important;
        margin-bottom: 0.4rem !important;
    }
    h4, div[data-testid="stMarkdownContainer"] h4 {
        font-size: 1.18rem !important;
        font-weight: 600 !important;
    }

    /* Enlarge Tab Names & Navigation */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        font-size: 17px !important;
        font-weight: 500 !important;
        padding: 10px 18px !important;
        border-radius: 4px 4px 0 0;
    }
    .stTabs [aria-selected="true"] {
        font-weight: 700 !important;
        color: inherit !important;
        border-bottom-color: #1e293b !important;
    }
    .stTabs [aria-selected="true"] p, .stTabs [aria-selected="true"] span {
        font-weight: 700 !important;
        color: inherit !important;
    }

    /* Enlarge Plot Controls, Widget Labels & Expander Headers */
    label[data-testid="stWidgetLabel"] p, label[data-testid="stWidgetLabel"] span {
        font-size: 16px !important;
        font-weight: 600 !important;
        color: #1e293b !important;
        margin-bottom: 4px !important;
    }
    .streamlit-expanderHeader p, details summary p {
        font-size: 16.5px !important;
        font-weight: 600 !important;
        color: #1e293b !important;
    }
    
    /* General Body, DataFrames, Inputs & Buttons */
    p, li {
        font-size: 15.5px !important;
    }
    .stCaption, caption, div[data-testid="stCaptionContainer"] p {
        font-size: 14px !important;
        line-height: 1.5 !important;
    }
    div[data-testid="stDataFrame"] {
        font-size: 15px !important;
    }
    .stSelectbox, .stMultiSelect, .stSlider, .stRadio, .stCheckbox {
        font-size: 15.5px !important;
    }
    .stButton button {
        font-size: 15px !important;
        font-weight: 600 !important;
        padding: 8px 18px !important;
    }
</style>
""", unsafe_allow_html=True)

DEFAULT_SIGNATURES = {
    "Adherens Junction Complex": ["CDH1", "CTNNB1", "CTNNA1", "CTNND1", "JUP"],
    "Desmosomes": ["DSP", "PKP1", "PKP3", "DSG1", "DSG3", "DSC1", "DSC3", "PPL", "EVPL"],
    "Hemidesmosome & Basement Membrane": ["COL17A1", "ITGB4", "ITGA6", "LAMA3", "LAMB3", "LAMC2", "DST"],
    "Cell Cycle / Proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "PCNA"]
}

# Default config if YAML is missing
DEFAULT_MARKERS = {
    "Fibroblasts": {
        "Activated / ECM-producing FB (Postn+)": ["Postn", "Fn1", "Col1a1", "Col3a1", "Cthrc1"],
        "Remodeling FB (Mmp3+)": ["Mmp3", "Mmp13", "Adamts4", "Col11a1"],
        "Inflammatory FB (C3+/Clu+)": ["C3", "Clu", "Apoe", "Apod", "Cxcl14"],
        "Inflammatory / Antigen-presenting FB (Cd74+)": ["Cd74", "H2-Ab1", "H2-Aa", "Cd86"],
        "Myofibroblast-like (Acta2+)": ["Acta2", "Tagln", "Myl9", "Tpm2"],
        "Reticular / Lower Dermal FB (Mest+)": ["Mest", "Dlk1", "Prg4", "Pcolce2"],
        "Pericyte / SMC (Contamination)": ["Rgs5", "Mcam", "Cspg4", "Pdgfrb", "Notch3"]
    },
    "Keratinocytes": {
        "IFE Basal / Basal stem-like (Krt14+/Krt5+)": ["Krt5", "Krt14", "Tp63", "Itga6", "Col17a1", "COL17A1", "KRT14", "KRT5", "TP63", "ITGA6", "ITGB4"],
        "IFE Suprabasal / Spinous (Krt10+/Krt1+)": ["Krt1", "Krt10", "Dsc1", "Dsg1a", "KRT1", "KRT10", "DSG1", "DSC1"],
        "Granular / Terminally Differentiated (Mt4+/Lor+)": ["Ivl", "Lor", "Flg", "Tgm1", "Krt2", "IVL", "LOR", "FLG", "TGM1"],
        "Adherens Junction Complex": ["Cdh1", "Ctnnb1", "Ctnna1", "Ctnnd1", "Dsp", "Jup", "CDH1", "CTNNB1", "CTNNA1", "CTNND1", "JUP"],
        "Desmosomes": ["Dsp", "Pkp1", "Pkp3", "Dsg3", "Dsc3", "Ppl", "Evpl", "DSP", "PKP1", "PKP3", "DSG3", "DSC3", "PPL", "EVPL"],
        "Hemidesmosome & Basement Membrane": ["Col17a1", "Itgb4", "Itga6", "Lama3", "Lamb3", "Lamc2", "Dst", "COL17A1", "ITGB4", "ITGA6", "LAMA3", "LAMB3", "LAMC2", "DST"],
        "Wound-activated (Migrating) KC": ["Krt6a", "Krt17", "Itgb6", "Sprr1b"],
        "Cycling KC (Top2a+/Mki67+)": ["Mki67", "Top2a", "Ccnb1", "Cdk1", "MKI67", "TOP2A", "CCNB1", "CDK1"]
    }
}

# Helper to load YAML
def load_yaml():
    if os.path.exists(YAML_PATH):
        try:
            with open(YAML_PATH, "r", encoding="utf-8") as f:
                content = yaml.safe_load(f)
                return content if content else DEFAULT_MARKERS
        except Exception:
            return DEFAULT_MARKERS
    return DEFAULT_MARKERS

# Helper to save YAML
def save_yaml(data):
    with open(YAML_PATH, "w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, default_flow_style=False, sort_keys=False)

# Scan folder dynamically for h5ad files
@st.cache_resource
def scan_datasets(data_dir):
    datasets = {}
    target_dirs = SCAN_DIRS if 'SCAN_DIRS' in globals() else [data_dir]
    for d in target_dirs:
        if not os.path.exists(d):
            continue
        for filename in sorted(os.listdir(d)):
            if filename.endswith(".h5ad") and not filename.startswith("."):
                filepath = os.path.join(d, filename)
                folder_tag = os.path.basename(d)
                try:
                    adata = sc.read_h5ad(filepath, backed='r')
                    dataset_name = adata.uns.get("dataset_name", f"{filename[:-5]} ({folder_tag})")
                    datasets[dataset_name] = filepath
                except Exception:
                    datasets[f"{filename[:-5]} ({folder_tag})"] = filepath
    return datasets

# Cache data loading
@st.cache_resource
def load_adata(h5ad_path):
    if not os.path.exists(h5ad_path):
        return None
    adata = sc.read_h5ad(h5ad_path)
    return adata

# Cache gene records: returns display_options with "gene_name (gene_id)" formatting and bidirectional lookups
@st.cache_data
def get_gene_display_mappings(_var_df, _var_names):
    records = []
    has_symbols = 'gene_symbols' in _var_df.columns
    has_ids = 'gene_ids' in _var_df.columns
    
    for idx_name in _var_names:
        row = _var_df.loc[idx_name] if idx_name in _var_df.index else None
        
        # Determine symbol
        if has_symbols and row is not None and pd.notna(row['gene_symbols']):
            sym = str(row['gene_symbols']).strip()
        else:
            sym = str(idx_name).strip()
            
        # Determine ID
        if has_ids and row is not None and pd.notna(row['gene_ids']):
            gid = str(row['gene_ids']).strip()
        elif idx_name.startswith(('ENSG', 'ENSMUSG', 'ENS')):
            gid = str(idx_name).strip()
        else:
            gid = None
            
        if gid and gid != sym:
            display_str = f"{sym} ({gid})"
        else:
            display_str = sym
            
        records.append((display_str, sym, gid, idx_name))
        
    records.sort(key=lambda x: x[1].upper())
    
    display_to_var = {r[0]: r[3] for r in records}
    sym_to_display = {r[1].upper(): r[0] for r in records}
    var_to_display = {r[3]: r[0] for r in records}
    display_options = ["None"] + [r[0] for r in records]
    
    return display_options, display_to_var, sym_to_display, var_to_display

# Helper to find categorical obs columns
def get_annotation_columns(adata):
    candidates = []
    for col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        if isinstance(dtype, pd.CategoricalDtype):
            candidates.append(col)
        elif dtype.name in ['object', 'string', 'int64', 'int32']:
            if adata.obs[col].nunique() < 100:
                candidates.append(col)
    preferred = ['predicted_labels', 'majority_voting', 'state_group', 'cell_state_annotated', 'cell_type', 'celltype', 'leiden', 'louvain', 'sample', 'condition']
    candidates = sorted(candidates, key=lambda x: (0 if x in preferred else 1, preferred.index(x) if x in preferred else x))
    return candidates

# Helper to find sample/condition column
def get_sample_column(adata):
    for col in ['sample', 'condition', 'donor_id', 'source', 'batch']:
        if col in adata.obs.columns:
            return col
    return adata.obs.columns[0] if len(adata.obs.columns) > 0 else 'sample'

# Significance label helper
def get_sig_label(p):
    if p < 0.0001: return "****"
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

# Format p-values and q-values with 4 sig figs and scientific notation for < 0.0001
def format_sig_value(val):
    if pd.isna(val):
        return "N/A"
    try:
        v = float(val)
        if v == 0.0:
            return "0.0"
        if abs(v) < 0.0001 or abs(v) >= 10000:
            return f"{v:.4e}"
        else:
            return f"{float(f'{v:.4g}'):g}"
    except Exception:
        return str(val)

# Resolve case-insensitive gene match to var_name index
def resolve_gene_var_name(adata, gene_name, sym_to_display, display_to_var):
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

# Sidebar menu
st.sidebar.title("Navigation")
app_mode = st.sidebar.selectbox("Choose the page:", [
    "Gene Expression UMAP", 
    "Cell-Type Marker Editor"
])

# Dynamic Dataset Picker
datasets_map = scan_datasets(DATA_DIR)
if not datasets_map:
    st.error("No `.h5ad` files found in specified directories.")
    st.stop()

selected_dataset_name = st.sidebar.selectbox("Select Dataset:", list(datasets_map.keys()))
h5ad_path = datasets_map[selected_dataset_name]

# Load selected dataset
adata = load_adata(h5ad_path)
if adata is None:
    st.error(f"Failed to load dataset at `{h5ad_path}`.")
    st.stop()

# Gene records with "gene_name (gene_id)" formatting
display_options, display_to_var, sym_to_display, var_to_display = get_gene_display_mappings(adata.var, adata.var_names)

# Dynamic Annotation & Sample Column Picker
anno_cols = get_annotation_columns(adata)
if not anno_cols:
    st.sidebar.warning("No categorical or annotation columns found in the dataset's `.obs`.")
    selected_col = None
else:
    selected_col = st.sidebar.selectbox("Annotation Column:", anno_cols)

sample_col = get_sample_column(adata)

# Determine dataset category
yaml_key = "Fibroblasts" if "fibroblast" in selected_dataset_name.lower() or "fb" in selected_dataset_name.lower() else "Keratinocytes"

# Canonical sample ordering
canonical_samples = ["NHEK", "RM1", "RM2", "nonRM", "NHEK_P4", "JEB_RM1", "JEB_RM2", "JEB_nonRM", "Normal", "JEB", "Revertant"]
if sample_col and sample_col in adata.obs.columns:
    unique_in_data = adata.obs[sample_col].dropna().unique().tolist()
    ordered_samples = [s for s in canonical_samples if s in unique_in_data] + [s for s in unique_in_data if s not in canonical_samples]
else:
    ordered_samples = []

sample_color_map = {"NHEK": "#d62728", "RM1": "#5c2d91", "RM2": "#b5a500", "nonRM": "#00a8e0", "NHEK_P4": "#d62728", "JEB_RM1": "#5c2d91", "JEB_RM2": "#b5a500", "JEB_nonRM": "#00a8e0", "Normal": "#2ecc71", "JEB": "#e74c3c", "Revertant": "#3498db"}
for idx, s in enumerate(ordered_samples):
    if s not in sample_color_map:
        cmap = plt.get_cmap('tab10')
        sample_color_map[s] = matplotlib.colors.to_hex(cmap(idx % 10))

# ----------------- PAGE 1: EXPRESSION VIEWER & ANALYSIS TABS -----------------
if app_mode == "Gene Expression UMAP":
    st.title("Single-Cell RNA-seq Expression, Scoring & Correlation Viewer")
    st.markdown(f'<div style="font-size: 20px; font-weight: 500; margin-top: 4px; margin-bottom: 18px; color: #1e293b; line-height: 1.5;">Active Dataset: <strong>{selected_dataset_name}</strong> | Total Cells: <strong>{adata.n_obs:,}</strong> | Total Genes: <strong>{adata.n_vars:,}</strong> | Sample Column: <strong>{sample_col}</strong></div>', unsafe_allow_html=True)
    
    # Initialize session state for selected gene display string
    if "selected_gene_display" not in st.session_state:
        st.session_state.selected_gene_display = "None"
    if "marker_select_key" not in st.session_state:
        st.session_state.marker_select_key = "None"
        
    markers_config = load_yaml()
    dataset_markers = markers_config.get(yaml_key, {})
    
    def on_marker_change():
        sel = st.session_state.marker_select_key
        if sel and sel != "None":
            disp = sym_to_display.get(sel.upper(), None)
            if disp:
                st.session_state.selected_gene_display = disp
                
    def on_gene_dropdown_change():
        sel = st.session_state.gene_dropdown_key
        st.session_state.selected_gene_display = sel
        st.session_state.marker_select_key = "None"
        
    col_gene, col_cell_type, col_marker = st.columns([1.6, 1.2, 1.2])
    
    with col_cell_type:
        cell_type_options = ["None"] + list(dataset_markers.keys())
        selected_ct = st.selectbox("Suggested Marker by Cell Type:", cell_type_options, index=0, key="ct_select_box")
        
    with col_marker:
        marker_options = ["None"]
        if selected_ct != "None":
            marker_options += dataset_markers.get(selected_ct, [])
        st.selectbox(
            "Select Marker:", 
            marker_options, 
            key="marker_select_key", 
            on_change=on_marker_change
        )
        
    with col_gene:
        curr_disp = st.session_state.selected_gene_display
        curr_idx = display_options.index(curr_disp) if curr_disp in display_options else 0
        
        selected_gene_box = st.selectbox(
            "Enter / Select Gene (gene_name (gene_id)):",
            options=display_options,
            index=curr_idx,
            key="gene_dropdown_key",
            on_change=on_gene_dropdown_change,
            help="Type to search by gene name or gene ID (Ensembl), then select to view plots."
        )
    
    resolved_var_name = display_to_var.get(selected_gene_box, None) if selected_gene_box != "None" else None
    resolved_display_name = selected_gene_box if selected_gene_box != "None" else None
    
    # Colormap & Scale Controls
    with st.expander("🎨 Colormap, Scale & Contrast Controls (Loupe-Style)", expanded=bool(resolved_var_name)):
        c_scale, c_cmap, c_pct, c_vmax = st.columns([1.2, 1.0, 1.8, 1.0])
        
        with c_scale:
            use_log2 = st.checkbox("Log2(Normalized + 1) Scale", value=True, help="Applies log2 transformation like Loupe Browser for balanced contrast.")
            
        with c_cmap:
            cmap_choice = st.selectbox("Colormap:", ["viridis", "YlOrRd", "Reds", "inferno", "plasma", "magma", "turbo"], index=0, key="global_cmap_choice")
            
        # Get raw expression values to calculate percentiles
        if resolved_var_name:
            if scipy.sparse.issparse(adata.X):
                raw_vals = adata[:, resolved_var_name].X.toarray().flatten()
            else:
                raw_vals = adata[:, resolved_var_name].X.flatten()
            
            expr_for_scale = np.log2(raw_vals + 1) if use_log2 else raw_vals
            max_possible = float(expr_for_scale.max())
            p80_val = float(np.percentile(expr_for_scale, 80))
            p90_val = float(np.percentile(expr_for_scale, 90))
            p95_val = float(np.percentile(expr_for_scale, 95))
            p99_val = float(np.percentile(expr_for_scale, 99))
        else:
            expr_for_scale = np.array([0.0])
            max_possible, p80_val, p90_val, p95_val, p99_val = 10.0, 4.0, 6.0, 8.0, 9.5
            
        with c_pct:
            pct_slider = st.select_slider(
                "Max Percentile Threshold (Anchors):",
                options=[50, 60, 70, 75, 80, 85, 90, 95, 98, 99, 99.5, 100],
                value=100,
                format_func=lambda x: f"{x}%" if x in [80, 90, 95, 99, 100] else f"{x}",
                help="Clip upper colormap limit to enhance visual contrast against outliers (anchors at 80%, 90%, 95%, 99%, 100%).",
                key="vmax_pct_slider"
            )
            
        suggested_vmax = float(np.percentile(expr_for_scale, pct_slider)) if resolved_var_name and len(expr_for_scale) > 1 else max_possible
        if suggested_vmax <= 0:
            suggested_vmax = max_possible if max_possible > 0 else 1.0
            
        with c_vmax:
            custom_vmax = st.number_input(
                "Colormap Max (vmax):",
                min_value=0.01,
                max_value=max(max_possible * 2.0, 10000.0),
                value=round(suggested_vmax, 2),
                step=0.5 if use_log2 else 10.0,
                help="Direct numeric limit for colormap maximum. Values above this will be saturated."
            )
            
        chosen_vmax = float(custom_vmax)
        chosen_scale_label = "Log2(Norm+1)" if use_log2 else "Linear"

    # Caching static multi-panel grid generator with Sample UMAP as 1st plot and custom vmax
    @st.cache_data
    def generate_static_grid(_adata, var_key, disp_title, col, s_col, dataset_name, log2_mode, v_max, cmap_name):
        if scipy.sparse.issparse(_adata.X):
            expr_raw = _adata[:, var_key].X.toarray().flatten()
        else:
            expr_raw = _adata[:, var_key].X.flatten()
            
        expr = np.log2(expr_raw + 1) if log2_mode else expr_raw
        
        umap_coords = _adata.obsm.get('X_umap', None)
        if umap_coords is None:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.text(0.5, 0.5, "No UMAP coordinates (X_umap) found in dataset", ha='center', va='center')
            ax.axis('off')
            return fig
            
        samples_list = [s for s in ordered_samples if s in _adata.obs[s_col].unique()][:4] if s_col else []
        num_splits = len(samples_list)
        total_plots = 3 + num_splits
        n_cols = 3 if total_plots >= 3 else total_plots
        n_rows = (total_plots + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 4.8 * n_rows))
        axes_flat = axes.flatten() if hasattr(axes, 'flatten') else [axes]
        
        # 1. Sample UMAP Plot (Leading 1st Plot)
        ax_samp = axes_flat[0]
        if s_col and s_col in _adata.obs.columns:
            for s in ordered_samples:
                if s in _adata.obs[s_col].values:
                    mask = _adata.obs[s_col] == s
                    coords = umap_coords[mask]
                    ax_samp.scatter(coords[:, 0], coords[:, 1], label=s, color=sample_color_map.get(s, "#7f8c8d"), s=1.5, alpha=0.85)
            ax_samp.set_title(f"Samples / Conditions ({s_col})", fontsize=12, fontweight='bold')
            ax_samp.set_xlabel("UMAP 1", fontsize=10)
            ax_samp.set_ylabel("UMAP 2", fontsize=10)
            ax_samp.tick_params(axis='both', which='major', labelsize=9.5)
            ax_samp.legend(title="Sample", bbox_to_anchor=(0.5, -0.2), loc="upper center", markerscale=6, fontsize=9.5, ncol=2, frameon=False)
        else:
            ax_samp.text(0.5, 0.5, "No Sample Column", ha='center', va='center', fontsize=11)
            ax_samp.axis('off')
        
        # 2. Cell States Reference
        ax_ct = axes_flat[1]
        if col:
            categories = _adata.obs[col].cat.categories.tolist() if hasattr(_adata.obs[col], "cat") else sorted(_adata.obs[col].dropna().unique().tolist())
            color_key = f"{col}_colors"
            if color_key in _adata.uns:
                colors = list(_adata.uns[color_key])
            else:
                cmap = plt.get_cmap('tab20')
                colors = [matplotlib.colors.to_hex(cmap(i % 20)) for i in range(len(categories))]
            color_map = dict(zip(categories, colors))
            
            for cat in categories:
                mask = _adata.obs[col] == cat
                coords = umap_coords[mask]
                ax_ct.scatter(coords[:, 0], coords[:, 1], label=cat, color=color_map.get(cat, "#7f8c8d"), s=1.5, alpha=0.8)
                
            ax_ct.set_title(f"Cell States ({col})", fontsize=12, fontweight='bold')
            ax_ct.set_xlabel("UMAP 1", fontsize=10)
            ax_ct.set_ylabel("UMAP 2", fontsize=10)
            ax_ct.tick_params(axis='both', which='major', labelsize=9.5)
            ax_ct.legend(title="Cell State", bbox_to_anchor=(0.5, -0.2), loc="upper center", markerscale=6, fontsize=9.0, ncol=2, frameon=False)
        else:
            ax_ct.text(0.5, 0.5, "No Annotation Column Selected", ha='center', va='center', fontsize=11)
            ax_ct.axis('off')

        # 3. Global Expression
        ax_all = axes_flat[2]
        sort_idx = np.argsort(expr)
        sc_all = ax_all.scatter(umap_coords[sort_idx, 0], umap_coords[sort_idx, 1], c=expr[sort_idx], cmap=cmap_name, s=1.5, alpha=0.85, vmin=0, vmax=v_max)
        clean_name = disp_title.split(" (")[0]
        scale_tag = "Log2" if log2_mode else "Linear"
        ax_all.set_title(f"{clean_name} - All Cells ({scale_tag}, max={v_max:.1f})", fontsize=12, fontweight='bold')
        ax_all.set_xlabel("UMAP 1", fontsize=10)
        ax_all.set_ylabel("UMAP 2", fontsize=10)
        ax_all.tick_params(axis='both', which='major', labelsize=9.5)
        cbar = fig.colorbar(sc_all, ax=ax_all, label=f"{scale_tag} Expr")
        cbar.ax.tick_params(labelsize=9.5)
        cbar.set_label(f"{scale_tag} Expr", fontsize=10)
        
        # 4..N. Splits per sample
        for idx, sample in enumerate(samples_list):
            ax_sub = axes_flat[3 + idx]
            mask = _adata.obs[s_col] == sample
            bg_x = umap_coords[~mask, 0]
            bg_y = umap_coords[~mask, 1]
            ax_sub.scatter(bg_x, bg_y, color='lightgrey', s=0.5, alpha=0.3)
            
            sub_expr = expr[mask]
            sub_coords = umap_coords[mask]
            sub_sort = np.argsort(sub_expr)
            sc_sub = ax_sub.scatter(sub_coords[sub_sort, 0], sub_coords[sub_sort, 1], c=sub_expr[sub_sort], cmap=cmap_name, s=1.5, alpha=0.85, vmin=0, vmax=v_max)
            ax_sub.set_title(f"{sample} only", fontsize=12, fontweight='bold')
            ax_sub.set_xlabel("UMAP 1", fontsize=10)
            ax_sub.set_ylabel("UMAP 2", fontsize=10)
            ax_sub.tick_params(axis='both', which='major', labelsize=9.5)
            cbar_s = fig.colorbar(sc_sub, ax=ax_sub)
            cbar_s.ax.tick_params(labelsize=9.5)
            
        for extra_idx in range(3 + num_splits, len(axes_flat)):
            axes_flat[extra_idx].axis('off')
            
        plt.tight_layout()
        return fig

    # 6 Main Tabs
    tab_static, tab_interactive, tab_composition, tab_gene_violin, tab_score_violin, tab_scatter = st.tabs([
        "Static Plots", 
        "Interactive Plots", 
        "Sample Composition",
        "Gene Expression Violins",
        "Signature & Pathway Scoring",
        "Correlation & Scatter Plots"
    ])
    
    # ---------------- TAB 1: STATIC PLOTS ----------------
    with tab_static:
        if resolved_var_name:
            with st.spinner("Generating static UMAP grid..."):
                fig_grid = generate_static_grid(adata, resolved_var_name, resolved_display_name, selected_col, sample_col, selected_dataset_name, use_log2, chosen_vmax, cmap_choice)
                st.pyplot(fig_grid)
                plt.close(fig_grid)
        else:
            st.info("💡 Select or search a gene above from the dropdown to view the expression comparison grid.")
            if 'X_umap' in adata.obsm:
                with st.spinner("Rendering reference UMAPs..."):
                    col_ref1, col_ref2 = st.columns(2)
                    umap_coords = adata.obsm['X_umap']
                    
                    # 1. Sample Reference
                    with col_ref1:
                        fig_s, ax_s = plt.subplots(figsize=(6, 5))
                        if sample_col and sample_col in adata.obs.columns:
                            for s in ordered_samples:
                                if s in adata.obs[sample_col].values:
                                    mask = adata.obs[sample_col] == s
                                    coords = umap_coords[mask]
                                    ax_s.scatter(coords[:, 0], coords[:, 1], label=s, color=sample_color_map.get(s, "#7f8c8d"), s=1.8, alpha=0.85)
                            ax_s.set_title(f"Samples / Conditions ({sample_col})", fontsize=11, fontweight='bold')
                            ax_s.set_xlabel("UMAP 1", fontsize=8)
                            ax_s.set_ylabel("UMAP 2", fontsize=8)
                            ax_s.legend(title="Sample", bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=5, fontsize=8, frameon=False)
                        plt.tight_layout()
                        st.pyplot(fig_s)
                        plt.close(fig_s)
                        
                    # 2. Cell States Reference
                    with col_ref2:
                        if selected_col:
                            fig_ref, ax_ref = plt.subplots(figsize=(6, 5))
                            categories = adata.obs[selected_col].cat.categories.tolist() if hasattr(adata.obs[selected_col], "cat") else sorted(adata.obs[selected_col].dropna().unique().tolist())
                            color_key = f"{selected_col}_colors"
                            if color_key in adata.uns:
                                colors = list(adata.uns[color_key])
                            else:
                                cmap = plt.get_cmap('tab20')
                                colors = [matplotlib.colors.to_hex(cmap(i % 20)) for i in range(len(categories))]
                            color_map = dict(zip(categories, colors))
                            
                            for cat in categories:
                                mask = adata.obs[selected_col] == cat
                                coords = umap_coords[mask]
                                ax_ref.scatter(coords[:, 0], coords[:, 1], label=cat, color=color_map.get(cat, "#7f8c8d"), s=1.8, alpha=0.8)
                                
                            ax_ref.set_title(f"Cell States ({selected_col})", fontsize=11, fontweight='bold')
                            ax_ref.set_xlabel("UMAP 1", fontsize=8)
                            ax_ref.set_ylabel("UMAP 2", fontsize=8)
                            ax_ref.legend(title="Cell State", bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=5, fontsize=7.5, frameon=False)
                            plt.tight_layout()
                            st.pyplot(fig_ref)
                            plt.close(fig_ref)
            
    # ---------------- TAB 2: INTERACTIVE PLOTS ----------------
    with tab_interactive:
        if 'X_umap' not in adata.obsm:
            st.warning("This dataset does not contain UMAP coordinates ('X_umap').")
        else:
            all_samples = ordered_samples if ordered_samples else (list(adata.obs[sample_col].unique()) if sample_col else [])
            if selected_col and hasattr(adata.obs[selected_col], "cat"):
                all_categories = list(adata.obs[selected_col].cat.categories)
            elif selected_col:
                all_categories = sorted(adata.obs[selected_col].dropna().unique().tolist())
            else:
                all_categories = []
            
            with st.expander("Filter & Highlight Controls", expanded=True):
                c_mode, c_size = st.columns([2, 1])
                with c_mode:
                    view_mode = st.radio(
                        "View / Highlight Mode:",
                        ["Color all cells", "Highlight selected (dim unselected in grey)", "Filter view (show selected only)"],
                        horizontal=True,
                        key="plotly_view_mode"
                    )
                with c_size:
                    pt_size = st.slider("Point Size:", min_value=1.5, max_value=8.0, value=3.5, step=0.5, key="plotly_pt_size")
                
                col_f1, col_f2 = st.columns(2)
                with col_f1:
                    selected_samples = st.multiselect(
                        "Filter / Highlight Samples:",
                        options=all_samples,
                        default=all_samples,
                        key="plotly_filter_samples"
                    )
                with col_f2:
                    if selected_col and all_categories:
                        selected_cats = st.multiselect(
                            f"Filter / Highlight {selected_col}:",
                            options=all_categories,
                            default=all_categories,
                            key="plotly_filter_categories"
                        )
                    else:
                        selected_cats = []
            
            with st.spinner("Rendering interactive Plotly UMAPs..."):
                max_cells = 10000
                if adata.n_obs > max_cells:
                    np.random.seed(42)
                    sub_idx = np.random.choice(adata.n_obs, max_cells, replace=False)
                    adata_sub = adata[sub_idx].copy()
                else:
                    adata_sub = adata
                    
                if resolved_var_name:
                    if scipy.sparse.issparse(adata_sub.X):
                        expr_sub_raw = adata_sub[:, resolved_var_name].X.toarray().flatten()
                    else:
                        expr_sub_raw = adata_sub[:, resolved_var_name].X.flatten()
                    expr_sub = np.log2(expr_sub_raw + 1) if use_log2 else expr_sub_raw
                else:
                    expr_sub = np.zeros(adata_sub.n_obs)
                    
                df_plotly = pd.DataFrame({
                    "UMAP 1": adata_sub.obsm['X_umap'][:, 0],
                    "UMAP 2": adata_sub.obsm['X_umap'][:, 1],
                    "Cell State": adata_sub.obs[selected_col].astype(str) if selected_col else "N/A",
                    "Sample": adata_sub.obs[sample_col].astype(str) if sample_col else "All",
                    "Expression": expr_sub
                })
                
                x_min, x_max = float(df_plotly["UMAP 1"].min()), float(df_plotly["UMAP 1"].max())
                y_min, y_max = float(df_plotly["UMAP 2"].min()), float(df_plotly["UMAP 2"].max())
                x_pad = (x_max - x_min) * 0.05
                y_pad = (y_max - y_min) * 0.05
                x_range = [x_min - x_pad, x_max + x_pad]
                y_range = [y_min - y_pad, y_max + y_pad]
                
                mask_samp = df_plotly["Sample"].isin(selected_samples) if selected_samples else pd.Series(True, index=df_plotly.index)
                if selected_col and selected_cats:
                    mask_c = df_plotly["Cell State"].isin(selected_cats)
                else:
                    mask_c = pd.Series(True, index=df_plotly.index)
                selected_mask = (mask_samp & mask_c)
                
                df_selected = df_plotly[selected_mask]
                df_unselected = df_plotly[~selected_mask]
                
                color_discrete_map = None
                if selected_col:
                    color_key = f"{selected_col}_colors"
                    if color_key in adata_sub.uns:
                        colors_list = list(adata_sub.uns[color_key])
                        unique_states = adata_sub.obs[selected_col].cat.categories.tolist() if hasattr(adata_sub.obs[selected_col], "cat") else sorted(adata_sub.obs[selected_col].dropna().unique().tolist())
                        color_discrete_map = dict(zip(unique_states, colors_list))
                        
                # 3 Side-by-Side Plots when gene selected, or 2 plots when no gene selected
                if resolved_var_name:
                    col_p0, col_p1, col_p2 = st.columns(3)
                else:
                    col_p0, col_p1 = st.columns(2)
                    col_p2 = None
                
                # Figure 0: Sample / Condition UMAP (Leading 1st Plot)
                with col_p0:
                    samp_title = f"Samples ({sample_col})" if sample_col else "Samples"
                    if view_mode == "Filter view (show selected only)":
                        fig_samp = px.scatter(
                            df_selected, x="UMAP 1", y="UMAP 2",
                            color="Sample",
                            hover_data=["Cell State", "Expression"],
                            color_discrete_map=sample_color_map,
                            title=f"{samp_title} (Filtered: {len(df_selected)} cells)",
                            template="plotly_white"
                        )
                        fig_samp.update_traces(marker=dict(size=pt_size, opacity=0.85))
                    elif view_mode == "Highlight selected (dim unselected in grey)":
                        fig_samp = go.Figure()
                        if not df_unselected.empty:
                            fig_samp.add_trace(go.Scattergl(
                                x=df_unselected["UMAP 1"], y=df_unselected["UMAP 2"],
                                mode='markers',
                                marker=dict(color='#F0F2F6', size=max(pt_size - 1.0, 1.5), opacity=0.9),
                                name='Other / Dimmed',
                                hoverinfo='skip'
                            ))
                        samps_to_plot = [s for s in ordered_samples if s in df_selected["Sample"].values]
                        for s in samps_to_plot:
                            sub_s_df = df_selected[df_selected["Sample"] == s]
                            if not sub_s_df.empty:
                                s_color = sample_color_map.get(s, None)
                                hover_txt = [
                                    f"Sample: {samp}<br>State: {st}<br>Expr: {e:.2f}"
                                    for samp, st, e in zip(sub_s_df["Sample"], sub_s_df["Cell State"], sub_s_df["Expression"])
                                ]
                                fig_samp.add_trace(go.Scattergl(
                                    x=sub_s_df["UMAP 1"], y=sub_s_df["UMAP 2"],
                                    mode='markers',
                                    marker=dict(color=s_color, size=pt_size, opacity=0.85),
                                    name=s,
                                    text=hover_txt,
                                    hoverinfo='text'
                                ))
                        fig_samp.update_layout(
                            title=f"{samp_title} (Highlighted: {len(df_selected)} / {len(df_plotly)} cells)",
                            template="plotly_white"
                        )
                    else:
                        fig_samp = px.scatter(
                            df_plotly, x="UMAP 1", y="UMAP 2",
                            color="Sample",
                            hover_data=["Cell State", "Expression"],
                            color_discrete_map=sample_color_map,
                            title=f"{samp_title} (All Cells)",
                            template="plotly_white"
                        )
                        fig_samp.update_traces(marker=dict(size=pt_size, opacity=0.85))
                        
                    fig_samp.update_xaxes(range=x_range, title_text="UMAP 1", zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_samp.update_yaxes(range=y_range, title_text="UMAP 2", scaleanchor="x", scaleratio=1, zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_samp.update_layout(
                        height=500,
                        margin=dict(l=10, r=10, t=40, b=10),
                        legend=dict(
                            itemsizing='constant',
                            font=dict(size=10, family="Segoe UI, sans-serif"),
                            title=dict(font=dict(size=11, family="Segoe UI, sans-serif")),
                            bgcolor="rgba(255,255,255,0.85)",
                            bordercolor="#E2E8F0",
                            borderwidth=1
                        )
                    )
                    st.plotly_chart(fig_samp, use_container_width=True)

                # Figure 1: Cell States
                with col_p1:
                    state_title = f"Cell States ({selected_col})" if selected_col else "Cell States"
                    if view_mode == "Filter view (show selected only)":
                        fig_states = px.scatter(
                            df_selected, x="UMAP 1", y="UMAP 2",
                            color="Cell State",
                            hover_data=["Sample", "Expression"],
                            color_discrete_map=color_discrete_map,
                            title=f"{state_title} (Filtered: {len(df_selected)} cells)",
                            template="plotly_white"
                        )
                        fig_states.update_traces(marker=dict(size=pt_size, opacity=0.85))
                    elif view_mode == "Highlight selected (dim unselected in grey)":
                        fig_states = go.Figure()
                        if not df_unselected.empty:
                            fig_states.add_trace(go.Scattergl(
                                x=df_unselected["UMAP 1"], y=df_unselected["UMAP 2"],
                                mode='markers',
                                marker=dict(color='#F0F2F6', size=max(pt_size - 1.0, 1.5), opacity=0.9),
                                name='Other / Dimmed',
                                hoverinfo='skip'
                            ))
                        cats_to_plot = all_categories if all_categories else sorted(df_selected["Cell State"].unique())
                        for cat in cats_to_plot:
                            sub_cat_df = df_selected[df_selected["Cell State"] == cat]
                            if not sub_cat_df.empty:
                                cat_color = color_discrete_map.get(cat, None) if color_discrete_map else None
                                hover_txt = [
                                    f"State: {s}<br>Sample: {samp}<br>Expr: {e:.2f}"
                                    for s, samp, e in zip(sub_cat_df["Cell State"], sub_cat_df["Sample"], sub_cat_df["Expression"])
                                ]
                                fig_states.add_trace(go.Scattergl(
                                    x=sub_cat_df["UMAP 1"], y=sub_cat_df["UMAP 2"],
                                    mode='markers',
                                    marker=dict(color=cat_color, size=pt_size, opacity=0.85),
                                    name=cat,
                                    text=hover_txt,
                                    hoverinfo='text'
                                ))
                        fig_states.update_layout(
                            title=f"{state_title} (Highlighted: {len(df_selected)} / {len(df_plotly)} cells)",
                            template="plotly_white"
                        )
                    else:
                        fig_states = px.scatter(
                            df_plotly, x="UMAP 1", y="UMAP 2",
                            color="Cell State",
                            hover_data=["Sample", "Expression"],
                            color_discrete_map=color_discrete_map,
                            title=f"{state_title} (All Cells)",
                            template="plotly_white"
                        )
                        fig_states.update_traces(marker=dict(size=pt_size, opacity=0.85))
                        
                    fig_states.update_xaxes(range=x_range, title_text="UMAP 1", zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_states.update_yaxes(range=y_range, title_text="UMAP 2", scaleanchor="x", scaleratio=1, zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_states.update_layout(
                        height=500,
                        margin=dict(l=10, r=10, t=40, b=10),
                        legend=dict(
                            itemsizing='constant',
                            font=dict(size=10, family="Segoe UI, sans-serif"),
                            title=dict(font=dict(size=11, family="Segoe UI, sans-serif")),
                            bgcolor="rgba(255,255,255,0.85)",
                            bordercolor="#E2E8F0",
                            borderwidth=1
                        )
                    )
                    st.plotly_chart(fig_states, use_container_width=True)

                # Figure 2: Gene Expression
                if resolved_var_name and col_p2:
                    clean_sym = resolved_display_name.split(" (")[0]
                    plotly_cs = cmap_choice if cmap_choice in ['viridis', 'inferno', 'plasma', 'magma', 'turbo', 'Reds', 'YlOrRd'] else 'viridis'
                    
                    with col_p2:
                        if view_mode == "Filter view (show selected only)":
                            fig_expr = px.scatter(
                                df_selected, x="UMAP 1", y="UMAP 2",
                                color="Expression",
                                hover_data=["Cell State", "Sample", "Expression"],
                                color_continuous_scale=plotly_cs,
                                range_color=[0, chosen_vmax],
                                title=f"{clean_sym} ({chosen_scale_label}, Filtered: {len(df_selected)})",
                                template="plotly_white"
                            )
                            fig_expr.update_traces(marker=dict(size=pt_size, opacity=0.85))
                        elif view_mode == "Highlight selected (dim unselected in grey)":
                            fig_expr = go.Figure()
                            if not df_unselected.empty:
                                fig_expr.add_trace(go.Scattergl(
                                    x=df_unselected["UMAP 1"], y=df_unselected["UMAP 2"],
                                    mode='markers',
                                    marker=dict(color='#F0F2F6', size=max(pt_size - 1.0, 1.5), opacity=0.9),
                                    name='Unselected',
                                    hoverinfo='skip'
                                ))
                            if not df_selected.empty:
                                hover_txt = [
                                    f"State: {s}<br>Sample: {samp}<br>Expr: {e:.2f}"
                                    for s, samp, e in zip(df_selected["Cell State"], df_selected["Sample"], df_selected["Expression"])
                                ]
                                fig_expr.add_trace(go.Scattergl(
                                    x=df_selected["UMAP 1"], y=df_selected["UMAP 2"],
                                    mode='markers',
                                    marker=dict(
                                        color=np.clip(df_selected["Expression"], 0, chosen_vmax),
                                        colorscale=plotly_cs,
                                        colorbar=dict(title=f"{chosen_scale_label}", len=0.8, thickness=14, tickfont=dict(size=10)),
                                        size=pt_size,
                                        opacity=0.85,
                                        cmin=0,
                                        cmax=chosen_vmax,
                                        showscale=True
                                    ),
                                    text=hover_txt,
                                    hoverinfo='text',
                                    name='Selected'
                                ))
                            fig_expr.update_layout(
                                title=f"{clean_sym} ({chosen_scale_label}, Highlighted: {len(df_selected)} / {len(df_plotly)})",
                                template="plotly_white"
                            )
                        else:
                            fig_expr = px.scatter(
                                df_plotly, x="UMAP 1", y="UMAP 2",
                                color="Expression",
                                hover_data=["Cell State", "Sample", "Expression"],
                                color_continuous_scale=plotly_cs,
                                range_color=[0, chosen_vmax],
                                title=f"{clean_sym} ({chosen_scale_label}, All Cells)",
                                template="plotly_white"
                            )
                            fig_expr.update_traces(marker=dict(size=pt_size, opacity=0.85))
                            
                        fig_expr.update_xaxes(range=x_range, title_text="UMAP 1", zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                        fig_expr.update_yaxes(range=y_range, title_text="UMAP 2", scaleanchor="x", scaleratio=1, zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                        fig_expr.update_layout(
                            height=500,
                            margin=dict(l=10, r=10, t=40, b=10),
                            legend=dict(itemsizing='constant', font=dict(size=10))
                        )
                        st.plotly_chart(fig_expr, use_container_width=True)

    # ---------------- TAB 3: SAMPLE COMPOSITION ----------------
    with tab_composition:
        if not selected_col or not sample_col:
            st.warning("Please select an Annotation Column in the sidebar to view cell population compositions.")
        else:
            all_dataset_samples = ordered_samples if ordered_samples else list(adata.obs[sample_col].unique())
            if hasattr(adata.obs[selected_col], "cat"):
                categories = adata.obs[selected_col].cat.categories.tolist()
            else:
                categories = sorted(adata.obs[selected_col].dropna().unique().tolist())
                
            color_key = f"{selected_col}_colors"
            if color_key in adata.uns:
                colors_list = list(adata.uns[color_key])
            else:
                cmap = plt.get_cmap('tab20')
                colors_list = [matplotlib.colors.to_hex(cmap(i % 20)) for i in range(len(categories))]
            color_map = dict(zip(categories, colors_list))
            
            with st.expander("Composition Display Options", expanded=True):
                c_bar_metric, c_samp_filter, c_donut_opt = st.columns([1.2, 2, 1])
                with c_bar_metric:
                    bar_metric = st.radio("Stacked Bar Metric:", ["Percentage (%)", "Absolute Counts"], horizontal=True, key="comp_bar_metric")
                with c_samp_filter:
                    selected_comp_samples = st.multiselect("Select Samples to Display:", options=all_dataset_samples, default=all_dataset_samples, key="comp_samples_filter")
                with c_donut_opt:
                    show_donut_pct = st.checkbox("Show % in Donut Slices", value=True, key="comp_show_pct")
                    
            if not selected_comp_samples:
                st.warning("Please select at least one sample to display composition charts.")
            else:
                ct_counts = pd.crosstab(adata.obs[selected_col], adata.obs[sample_col]).reindex(index=categories, columns=selected_comp_samples, fill_value=0)
                sample_totals = ct_counts.sum(axis=0)
                ct_pct = (ct_counts / sample_totals.replace(0, 1)) * 100
                
                col_chart1, col_chart2 = st.columns([1.1, 1.4])
                
                with col_chart1:
                    st.markdown("#### Stacked Population Composition")
                    fig_bar = go.Figure()
                    for cat in categories:
                        cat_col = color_map.get(cat, "#7f8c8d")
                        if bar_metric == "Percentage (%)":
                            y_vals = ct_pct.loc[cat]
                            hover_template = '<b>' + cat + '</b><br>Sample: %{x}<br>Percentage: %{y:.1f}%<br>Count: %{customdata:,} cells<extra></extra>'
                            custom_data = ct_counts.loc[cat]
                        else:
                            y_vals = ct_counts.loc[cat]
                            hover_template = '<b>' + cat + '</b><br>Sample: %{x}<br>Count: %{y:,} cells<br>Percentage: %{customdata:.1f}%<extra></extra>'
                            custom_data = ct_pct.loc[cat]
                            
                        fig_bar.add_trace(go.Bar(
                            x=selected_comp_samples,
                            y=y_vals,
                            name=cat,
                            marker=dict(color=cat_col, line=dict(color='#FFFFFF', width=0.5)),
                            customdata=custom_data,
                            hovertemplate=hover_template
                        ))
                    
                    yaxis_title = "Percentage (%)" if bar_metric == "Percentage (%)" else "Cell Count"
                    yaxis_range = [0, 100] if bar_metric == "Percentage (%)" else None
                    fig_bar.update_layout(
                        barmode='stack',
                        template='plotly_white',
                        height=480,
                        margin=dict(l=10, r=10, t=30, b=10),
                        xaxis=dict(title="Sample / Condition", tickfont=dict(size=11, family="Segoe UI, sans-serif")),
                        yaxis=dict(title=yaxis_title, range=yaxis_range, showgrid=True, gridcolor="#F1F5F9"),
                        legend=dict(
                            itemsizing='constant',
                            font=dict(size=10, family="Segoe UI, sans-serif"),
                            title=dict(font=dict(size=11, family="Segoe UI, sans-serif")),
                            bgcolor="rgba(255,255,255,0.85)",
                            bordercolor="#E2E8F0",
                            borderwidth=1
                        )
                    )
                    st.plotly_chart(fig_bar, use_container_width=True)
                    
                with col_chart2:
                    st.markdown("#### Sample Percentage Ring (Donut) Charts")
                    num_s = len(selected_comp_samples)
                    if num_s <= 4:
                        n_rows, n_cols = 1, num_s
                        donut_height = 480
                    else:
                        n_cols = min(num_s, 3)
                        n_rows = (num_s + n_cols - 1) // n_cols
                        donut_height = 250 * n_rows
                        
                    donut_specs = [[{"type": "domain"} for _ in range(n_cols)] for _ in range(n_rows)]
                    subplot_titles = [f"<b>{s}</b><br><sub>Total: {sample_totals[s]:,} cells</sub>" for s in selected_comp_samples]
                    
                    fig_donuts = make_subplots(
                        rows=n_rows, cols=n_cols,
                        subplot_titles=subplot_titles,
                        specs=donut_specs,
                        horizontal_spacing=0.04,
                        vertical_spacing=0.15
                    )
                    
                    for idx, s in enumerate(selected_comp_samples):
                        r_idx = (idx // n_cols) + 1
                        c_idx = (idx % n_cols) + 1
                        s_counts = ct_counts[s]
                        
                        fig_donuts.add_trace(
                            go.Pie(
                                labels=categories,
                                values=s_counts,
                                hole=0.55,
                                marker=dict(
                                    colors=[color_map.get(c, "#7f8c8d") for c in categories],
                                    line=dict(color='#FFFFFF', width=1.5)
                                ),
                                textinfo='percent' if show_donut_pct else 'none',
                                textposition='inside',
                                hovertemplate='<b>%{label}</b><br>Sample: ' + s + '<br>Count: %{value:,} cells<br>Percentage: %{percent}<extra></extra>',
                                showlegend=False,
                                sort=False
                            ),
                            row=r_idx, col=c_idx
                        )
                        
                    fig_donuts.update_layout(
                        template='plotly_white',
                        height=donut_height,
                        margin=dict(l=10, r=10, t=50, b=10)
                    )
                    st.plotly_chart(fig_donuts, use_container_width=True)
                    
                with st.expander("📊 View Composition Data Tables & Export"):
                    tab_cnt, tab_prc = st.tabs(["Cell Counts Table", "Percentages (%) Table"])
                    
                    with tab_cnt:
                        df_counts_display = ct_counts.copy()
                        df_counts_display.loc["Total (All Populations)"] = sample_totals
                        st.dataframe(df_counts_display, use_container_width=True)
                        
                    with tab_prc:
                        df_pct_display = ct_pct.round(2).copy()
                        st.dataframe(df_pct_display, use_container_width=True)
                        
                    combined_export_df = pd.concat({
                        "Cell Counts": ct_counts,
                        "Percentages (%)": ct_pct.round(3)
                    }, axis=1)
                    csv_data = combined_export_df.to_csv().encode('utf-8')
                    st.download_button(
                        label="📥 Download Composition Tables as CSV",
                        data=csv_data,
                        file_name=f"{selected_dataset_name}_{selected_col}_composition.csv",
                        mime="text/csv"
                    )

    # ---------------- TAB 4: GENE EXPRESSION VIOLINS (WITH STATS) ----------------
    with tab_gene_violin:
        st.markdown("### Gene Expression Significance Violins across Conditions")
        if not resolved_var_name:
            st.info("💡 Please select or search a gene above from the dropdown to generate expression violin plots with statistical significance testing.")
        elif not selected_col or not sample_col:
            st.warning("Please select an Annotation Column and verify sample columns.")
        else:
            if scipy.sparse.issparse(adata.X):
                expr_vals_raw = adata[:, resolved_var_name].X.toarray().flatten()
            else:
                expr_vals_raw = adata[:, resolved_var_name].X.flatten()
                
            expr_vals = np.log2(expr_vals_raw + 1) if use_log2 else expr_vals_raw
            clean_sym = resolved_display_name.split(" (")[0]
            df_gene_violin = pd.DataFrame({
                "Cell State": adata.obs[selected_col].astype(str),
                "Sample": adata.obs[sample_col].astype(str),
                "Expression": expr_vals
            })
            
            all_states = sorted(df_gene_violin["Cell State"].unique())
            available_samples = [s for s in ordered_samples if s in df_gene_violin["Sample"].unique()]
            if not available_samples:
                available_samples = sorted(df_gene_violin["Sample"].unique())
                
            with st.expander("Violin & Statistical Comparison Options", expanded=True):
                c_v0, c_v1, c_v2, c_v3 = st.columns([1.2, 1.8, 0.9, 1.1])
                with c_v0:
                    include_all_cells = st.checkbox("Include 'All Cells (Global)' as 1st Plot", value=True, key="v_include_all")
                    filter_zeros_v = st.checkbox("Remove Expression = 0 Cells", value=False, help="Restricts violin analysis and statistical testing to expressing cells only (>0).", key="v_filter_zeros")
                with c_v1:
                    selected_states_v = st.multiselect("Select Cell States to Include:", options=all_states, default=all_states, key="v_gene_states")
                with c_v2:
                    plot_ncols = st.selectbox("Subplot Grid Columns:", [2, 3, 4], index=1, key="v_gene_ncols")
                with c_v3:
                    gene_max_input = st.number_input(
                        "Expression Y-Axis Max (0 for Auto):",
                        min_value=0.0,
                        max_value=10000.0,
                        value=0.0,
                        step=0.5 if use_log2 else 10.0,
                        help="Set fixed upper limit for expression Y-axis across all subplots, or keep 0 for automatic."
                    )
                    
                if filter_zeros_v:
                    df_gene_violin = df_gene_violin[df_gene_violin["Expression"] > 0]
                    
                candidate_pairs = []
                for i in range(len(available_samples)):
                    for j in range(i+1, len(available_samples)):
                        candidate_pairs.append((available_samples[i], available_samples[j]))
                        
                default_active_pairs = [("NHEK", "nonRM"), ("RM2", "nonRM"), ("RM1", "nonRM")]
                valid_default_pairs = [p for p in default_active_pairs if p[0] in available_samples and p[1] in available_samples]
                if not valid_default_pairs and len(candidate_pairs) > 0:
                    valid_default_pairs = candidate_pairs[:2]
                    
                pair_options = [f"{p[0]} vs {p[1]}" for p in candidate_pairs]
                default_pair_strs = [f"{p[0]} vs {p[1]}" for p in valid_default_pairs]
                
                selected_pair_strs = st.multiselect("Select Comparison Pairs for Significance Brackets (S1 vs S2):", options=pair_options, default=default_pair_strs, key="v_gene_pairs")
                active_pairs = [tuple(s.split(" vs ")) for s in selected_pair_strs]

            plots_to_generate = (["All Cells (Global)"] if include_all_cells else []) + selected_states_v
            
            if not plots_to_generate:
                st.warning("Please select at least one plot to display.")
            else:
                with st.spinner("Generating statistical violin plots..."):
                    n_plots = len(plots_to_generate)
                    n_cols = min(plot_ncols, n_plots)
                    n_rows = (n_plots + n_cols - 1) // n_cols
                    
                    fig_v, axes_v = plt.subplots(n_rows, n_cols, figsize=(4.4 * n_cols, 3.8 * n_rows + 0.6))
                    axes_v_flat = axes_v.flatten() if hasattr(axes_v, 'flatten') else [axes_v]
                    
                    # Descriptive Super Title
                    fig_v.suptitle(f"{resolved_display_name} - Single-Cell Expression across Conditions ({chosen_scale_label})", fontsize=13, fontweight='bold', y=0.995)
                    
                    stats_records = []
                    
                    for idx, ct in enumerate(plots_to_generate):
                        ax = axes_v_flat[idx]
                        if ct == "All Cells (Global)":
                            df_sub = df_gene_violin
                            title_display = f"★ All Cells (Global, N={len(df_sub):,})"
                        else:
                            df_sub = df_gene_violin[df_gene_violin["Cell State"] == ct]
                            title_clean = ct.split("(")[0].strip()
                            title_display = title_clean[:26] + "..." if len(title_clean) > 28 else title_clean
                        
                        sns.violinplot(
                            data=df_sub, x="Sample", y="Expression", ax=ax,
                            palette=sample_color_map, order=available_samples, hue="Sample", legend=False,
                            inner=None, cut=0, alpha=0.85
                        )
                        sns.boxplot(
                            data=df_sub, x="Sample", y="Expression", ax=ax,
                            width=0.12, boxprops={'facecolor':'white', 'edgecolor':'black', 'linewidth':0.8},
                            whiskerprops={'color':'black', 'linewidth':0.8}, capprops={'color':'black', 'linewidth':0.8},
                            medianprops={'color':'red', 'linewidth':1.3},
                            showmeans=True,
                            meanprops={'marker':'D', 'markerfacecolor':'#0284c7', 'markeredgecolor':'black', 'markersize':4.5},
                            showfliers=False, order=available_samples
                        )
                        
                        ymin, ymax = ax.get_ylim()
                        y_range = max(ymax - ymin, 0.1)
                        y_curr = ymax + y_range * 0.05
                        h_bracket = y_range * 0.035
                        
                        for p_idx, (s1, s2) in enumerate(active_pairs):
                            if s1 in available_samples and s2 in available_samples:
                                idx1 = available_samples.index(s1)
                                idx2 = available_samples.index(s2)
                                
                                val1 = df_sub[df_sub["Sample"] == s1]["Expression"].values
                                val2 = df_sub[df_sub["Sample"] == s2]["Expression"].values
                                
                                p_val = 1.0
                                u_stat = np.nan
                                if len(val1) >= 3 and len(val2) >= 3:
                                    try:
                                        res = mannwhitneyu(val1, val2, alternative="two-sided")
                                        p_val = res.pvalue
                                        u_stat = res.statistic
                                    except Exception:
                                        p_val = 1.0
                                        
                                sig_str = get_sig_label(p_val)
                                m1, m2 = np.mean(val1) if len(val1)>0 else 0, np.mean(val2) if len(val2)>0 else 0
                                med1, med2 = np.median(val1) if len(val1)>0 else 0, np.median(val2) if len(val2)>0 else 0
                                stats_records.append({
                                    "Gene": resolved_display_name,
                                    "Population": ct,
                                    "Comparison (S1 vs S2)": f"{s1} vs {s2}",
                                    "N (S1)": len(val1),
                                    "Mean (S1)": round(float(m1), 3),
                                    "Median (S1)": round(float(med1), 3),
                                    "N (S2)": len(val2),
                                    "Mean (S2)": round(float(m2), 3),
                                    "Median (S2)": round(float(med2), 3),
                                    "Mean Diff (S1 - S2)": round(float(m1 - m2), 3),
                                    "Median Diff (S1 - S2)": round(float(med1 - med2), 3),
                                    "Mann-Whitney U": u_stat,
                                    "p-value_raw": p_val,
                                    "Significance": sig_str
                                })
                                
                                ax.plot([idx1, idx1, idx2, idx2], [y_curr, y_curr + h_bracket, y_curr + h_bracket, y_curr], lw=0.8, color="black")
                                ax.text((idx1 + idx2)/2.0, y_curr + h_bracket + y_range * 0.01, sig_str, ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="black")
                                y_curr += y_range * 0.16
                                
                        if gene_max_input > 0:
                            ax.set_ylim(ymin, gene_max_input)
                        else:
                            ax.set_ylim(ymin, y_curr + y_range * 0.05)
                        ax.set_title(title_display, fontsize=11.5, fontweight="bold")
                        ax.set_xlabel("")
                        ax.set_ylabel(f"{clean_sym} ({chosen_scale_label})" if (idx % n_cols) == 0 else "", fontsize=11)
                        ax.tick_params(axis='both', which='major', labelsize=10)
                        
                    for extra_idx in range(n_plots, len(axes_v_flat)):
                        axes_v_flat[extra_idx].axis('off')
                        
                    plt.tight_layout(rect=[0, 0, 1, 0.96])
                    st.pyplot(fig_v)
                    plt.close(fig_v)
                    st.caption("ℹ️ **Boxplot Legend**: 🔴 **Red line** = Median (50th percentile) &nbsp;|&nbsp; 🔷 **Blue diamond** = Mean (arithmetic average). For skewed single-cell distributions, the Mean is pulled by high-expressing cells.")
                    
                    # 1. Sample Mean & Median Summary Tables
                    st.markdown("#### Sample Expression Summary Tables")
                    mean_rows = []
                    median_rows = []
                    for pop in plots_to_generate:
                        df_pop = df_gene_violin if pop == "All Cells (Global)" else df_gene_violin[df_gene_violin["Cell State"] == pop]
                        r_mean = {"Population / Cluster": pop}
                        r_median = {"Population / Cluster": pop}
                        for s in available_samples:
                            s_vals = df_pop[df_pop["Sample"] == s]["Expression"].values
                            r_mean[s] = round(float(np.mean(s_vals)), 3) if len(s_vals) > 0 else np.nan
                            r_median[s] = round(float(np.median(s_vals)), 3) if len(s_vals) > 0 else np.nan
                        mean_rows.append(r_mean)
                        median_rows.append(r_median)
                        
                    df_sample_means = pd.DataFrame(mean_rows).set_index("Population / Cluster")
                    df_sample_medians = pd.DataFrame(median_rows).set_index("Population / Cluster")
                    
                    tab_m_tbl, tab_med_tbl = st.tabs(["Mean Expression Table", "Median Expression Table"])
                    with tab_m_tbl:
                        st.dataframe(df_sample_means, use_container_width=True)
                        csv_means = df_sample_means.to_csv().encode('utf-8')
                        st.download_button(
                            label=f"📥 Download {clean_sym} Sample Means CSV",
                            data=csv_means,
                            file_name=f"{selected_dataset_name}_{clean_sym}_sample_means.csv",
                            mime="text/csv",
                            key="dl_gene_means_csv"
                        )
                    with tab_med_tbl:
                        st.dataframe(df_sample_medians, use_container_width=True)
                        csv_medians = df_sample_medians.to_csv().encode('utf-8')
                        st.download_button(
                            label=f"📥 Download {clean_sym} Sample Medians CSV",
                            data=csv_medians,
                            file_name=f"{selected_dataset_name}_{clean_sym}_sample_medians.csv",
                            mime="text/csv",
                            key="dl_gene_medians_csv"
                        )
                    
                    # 2. Statistical Significance Table
                    if stats_records:
                        st.markdown("#### Statistical Significance Summary Table (Mann-Whitney U Test)")
                        st.caption("ℹ️ **Statistical Note**: Mann-Whitney U is a non-parametric rank-sum test evaluating whether S1 values tend to be higher than S2 (evaluating rank/median shifts). Both parametric Means and non-parametric Medians are reported for clarity.")
                        df_stats = pd.DataFrame(stats_records)
                        if len(df_stats) > 1:
                            from scipy.stats import false_discovery_control
                            try:
                                df_stats["FDR_raw"] = false_discovery_control(df_stats["p-value_raw"])
                            except Exception:
                                df_stats["FDR_raw"] = df_stats["p-value_raw"] * len(df_stats)
                        else:
                            df_stats["FDR_raw"] = df_stats["p-value_raw"]
                            
                        # Format 4 sig figs / scientific notation
                        df_stats["p-value"] = df_stats["p-value_raw"].apply(format_sig_value)
                        df_stats["FDR (q-value)"] = df_stats["FDR_raw"].apply(format_sig_value)
                        
                        cols_to_show = [
                            "Gene", "Population", "Comparison (S1 vs S2)", 
                            "N (S1)", "Mean (S1)", "Median (S1)", "N (S2)", "Mean (S2)", "Median (S2)",
                            "Mean Diff (S1 - S2)", "Median Diff (S1 - S2)", "Mann-Whitney U", "p-value", "FDR (q-value)", "Significance"
                        ]
                        df_stats_display = df_stats[cols_to_show]
                        st.dataframe(df_stats_display, use_container_width=True)
                        csv_stats = df_stats_display.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label=f"📥 Download {clean_sym} Statistical Summary CSV",
                            data=csv_stats,
                            file_name=f"{selected_dataset_name}_{clean_sym}_significance_stats.csv",
                            mime="text/csv",
                            key="dl_gene_stats_csv"
                        )

    # ---------------- TAB 5: SIGNATURE & PATHWAY SCORING (WITH STATS) ----------------
    with tab_score_violin:
        st.markdown("### Gene Signature & Pathway Scoring Violins (with Stats)")
        st.write("Calculate dynamic gene set scores (`sc.tl.score_genes`) and perform statistical comparisons across condition groups per cell state.")
        
        sig_options = list(DEFAULT_SIGNATURES.keys()) + ["Custom Gene List"]
        c_sig1, c_sig2 = st.columns([1.2, 2])
        with c_sig1:
            selected_sig_name = st.selectbox("Select Signature / Pathway:", sig_options, index=0, key="score_sig_choice")
        with c_sig2:
            if selected_sig_name == "Custom Gene List":
                custom_sig_str = st.text_input("Enter comma-separated gene list:", value="Cdh1, Ctnnb1, Dsp, Jup, Col17a1", key="score_custom_genes")
                target_genes = [g.strip() for g in custom_sig_str.split(",") if g.strip()]
            else:
                target_genes = DEFAULT_SIGNATURES[selected_sig_name]
                st.markdown(f'<div style="margin-top: 30px; font-size: 15.5px; color: #334155;"><strong>Genes in panel:</strong> {", ".join(target_genes)}</div>', unsafe_allow_html=True)
                
        resolved_sig_var_keys = []
        resolved_sig_names = []
        for g in target_genes:
            v_key = resolve_gene_var_name(adata, g, sym_to_display, display_to_var)
            if v_key and v_key not in resolved_sig_var_keys:
                resolved_sig_var_keys.append(v_key)
                disp_str = var_to_display.get(v_key, v_key).split(" (")[0]
                resolved_sig_names.append(disp_str)
                
        if not resolved_sig_var_keys:
            st.error(f"None of the genes in panel {target_genes} were found in this dataset.")
        else:
            st.success(f"Detected **{len(resolved_sig_var_keys)} / {len(target_genes)}** genes in dataset: `{', '.join(resolved_sig_names)}`")
            
            @st.cache_data
            def compute_score(_adata, gene_list, score_label, dataset_tag):
                adata_temp = _adata.copy()
                sc.tl.score_genes(adata_temp, gene_list=gene_list, score_name=score_label)
                return adata_temp.obs[score_label].values
                
            score_col_name = f"sig_score_{selected_sig_name.lower().replace(' ', '_')[:20]}"
            score_values = compute_score(adata, resolved_sig_var_keys, score_col_name, selected_dataset_name)
            
            df_score_violin = pd.DataFrame({
                "Cell State": adata.obs[selected_col].astype(str) if selected_col else "All",
                "Sample": adata.obs[sample_col].astype(str) if sample_col else "All",
                "Score": score_values
            })
            
            all_states_s = sorted(df_score_violin["Cell State"].unique())
            available_samples_s = [s for s in ordered_samples if s in df_score_violin["Sample"].unique()]
            if not available_samples_s:
                available_samples_s = sorted(df_score_violin["Sample"].unique())
                
            with st.expander("Scoring Violin Display Options", expanded=True):
                c_sv0, c_sv1, c_sv2, c_sv3 = st.columns([1.1, 1.8, 1.0, 1.0])
                with c_sv0:
                    include_all_cells_s = st.checkbox("Include 'All Cells (Global)' as 1st Plot", value=True, key="sv_include_all")
                with c_sv1:
                    selected_states_sv = st.multiselect("Select Cell States for Scoring Violins:", options=all_states_s, default=all_states_s, key="sv_states")
                with c_sv2:
                    plot_ncols_s = st.selectbox("Subplot Grid Columns:", [2, 3, 4], index=1, key="sv_ncols")
                with c_sv3:
                    score_max_input = st.number_input(
                        "Score Y-Axis Max (0 for Auto):",
                        min_value=0.0,
                        max_value=100.0,
                        value=0.0,
                        step=0.5,
                        help="Set fixed upper limit for score Y-axis across all subplots, or keep 0 for automatic."
                    )
                    
                candidate_pairs_s = []
                for i in range(len(available_samples_s)):
                    for j in range(i+1, len(available_samples_s)):
                        candidate_pairs_s.append((available_samples_s[i], available_samples_s[j]))
                        
                default_active_pairs_s = [("NHEK", "nonRM"), ("RM2", "nonRM"), ("RM1", "nonRM")]
                valid_default_pairs_s = [p for p in default_active_pairs_s if p[0] in available_samples_s and p[1] in available_samples_s]
                if not valid_default_pairs_s and len(candidate_pairs_s) > 0:
                    valid_default_pairs_s = candidate_pairs_s[:2]
                    
                pair_options_s = [f"{p[0]} vs {p[1]}" for p in candidate_pairs_s]
                default_pair_strs_s = [f"{p[0]} vs {p[1]}" for p in valid_default_pairs_s]
                
                selected_pair_strs_s = st.multiselect("Select Comparison Pairs for Significance Brackets (S1 vs S2):", options=pair_options_s, default=default_pair_strs_s, key="sv_pairs")
                active_pairs_s = [tuple(s.split(" vs ")) for s in selected_pair_strs_s]
                
            score_plots_to_generate = (["All Cells (Global)"] if include_all_cells_s else []) + selected_states_sv
            
            if not score_plots_to_generate:
                st.warning("Please select at least one cell state.")
            else:
                with st.spinner("Rendering scoring violin plots..."):
                    n_plots_s = len(score_plots_to_generate)
                    n_cols_s = min(plot_ncols_s, n_plots_s)
                    n_rows_s = (n_plots_s + n_cols_s - 1) // n_cols_s
                    
                    fig_sv, axes_sv = plt.subplots(n_rows_s, n_cols_s, figsize=(4.4 * n_cols_s, 3.8 * n_rows_s + 0.6))
                    axes_sv_flat = axes_sv.flatten() if hasattr(axes_sv, 'flatten') else [axes_sv]
                    
                    # Descriptive Super Title
                    fig_sv.suptitle(f"{selected_sig_name} Score", fontsize=13.5, fontweight='bold', y=0.99)
                    
                    score_stats_records = []
                    
                    for idx, ct in enumerate(score_plots_to_generate):
                        ax = axes_sv_flat[idx]
                        if ct == "All Cells (Global)":
                            df_sub = df_score_violin
                            title_display_s = f"★ All Cells (Global, N={len(df_sub):,})"
                        else:
                            df_sub = df_score_violin[df_score_violin["Cell State"] == ct]
                            title_clean = ct.split("(")[0].strip()
                            title_display_s = title_clean[:26] + "..." if len(title_clean) > 28 else title_clean
                        
                        sns.violinplot(
                            data=df_sub, x="Sample", y="Score", ax=ax,
                            palette=sample_color_map, order=available_samples_s, hue="Sample", legend=False,
                            inner=None, cut=0, alpha=0.85
                        )
                        sns.boxplot(
                            data=df_sub, x="Sample", y="Score", ax=ax,
                            width=0.12, boxprops={'facecolor':'white', 'edgecolor':'black', 'linewidth':0.8},
                            whiskerprops={'color':'black', 'linewidth':0.8}, capprops={'color':'black', 'linewidth':0.8},
                            medianprops={'color':'red', 'linewidth':1.3},
                            showmeans=True,
                            meanprops={'marker':'D', 'markerfacecolor':'#0284c7', 'markeredgecolor':'black', 'markersize':4.5},
                            showfliers=False, order=available_samples_s
                        )
                        
                        ymin, ymax = ax.get_ylim()
                        y_range = max(ymax - ymin, 0.1)
                        y_curr = ymax + y_range * 0.05
                        h_bracket = y_range * 0.035
                        
                        for p_idx, (s1, s2) in enumerate(active_pairs_s):
                            if s1 in available_samples_s and s2 in available_samples_s:
                                idx1 = available_samples_s.index(s1)
                                idx2 = available_samples_s.index(s2)
                                
                                val1 = df_sub[df_sub["Sample"] == s1]["Score"].values
                                val2 = df_sub[df_sub["Sample"] == s2]["Score"].values
                                
                                p_val = 1.0
                                u_stat = np.nan
                                if len(val1) >= 3 and len(val2) >= 3:
                                    try:
                                        res = mannwhitneyu(val1, val2, alternative="two-sided")
                                        p_val = res.pvalue
                                        u_stat = res.statistic
                                    except Exception:
                                        p_val = 1.0
                                        
                                sig_str = get_sig_label(p_val)
                                m1, m2 = np.mean(val1) if len(val1)>0 else 0, np.mean(val2) if len(val2)>0 else 0
                                med1, med2 = np.median(val1) if len(val1)>0 else 0, np.median(val2) if len(val2)>0 else 0
                                score_stats_records.append({
                                    "Signature": selected_sig_name,
                                    "Population": ct,
                                    "Comparison (S1 vs S2)": f"{s1} vs {s2}",
                                    "N (S1)": len(val1),
                                    "Mean Score (S1)": round(float(m1), 3),
                                    "Median Score (S1)": round(float(med1), 3),
                                    "N (S2)": len(val2),
                                    "Mean Score (S2)": round(float(m2), 3),
                                    "Median Score (S2)": round(float(med2), 3),
                                    "Mean Diff (S1 - S2)": round(float(m1 - m2), 3),
                                    "Median Diff (S1 - S2)": round(float(med1 - med2), 3),
                                    "Mann-Whitney U": u_stat,
                                    "p-value_raw": p_val,
                                    "Significance": sig_str
                                })
                                
                                ax.plot([idx1, idx1, idx2, idx2], [y_curr, y_curr + h_bracket, y_curr + h_bracket, y_curr], lw=0.8, color="black")
                                ax.text((idx1 + idx2)/2.0, y_curr + h_bracket + y_range * 0.01, sig_str, ha="center", va="bottom", fontsize=9.5, fontweight="bold", color="black")
                                y_curr += y_range * 0.16
                                
                        if score_max_input > 0:
                            ax.set_ylim(ymin, score_max_input)
                        else:
                            ax.set_ylim(ymin, y_curr + y_range * 0.05)
                            
                        ax.set_title(title_display_s, fontsize=11.5, fontweight="bold")
                        ax.set_xlabel("")
                        ax.set_ylabel(f"{selected_sig_name[:12]} Score" if (idx % n_cols_s) == 0 else "", fontsize=11)
                        ax.tick_params(axis='both', which='major', labelsize=10)
                        
                    for extra_idx in range(n_plots_s, len(axes_sv_flat)):
                        axes_sv_flat[extra_idx].axis('off')
                        
                    plt.tight_layout(rect=[0, 0, 1, 0.96])
                    st.pyplot(fig_sv)
                    plt.close(fig_sv)
                    st.caption("ℹ️ **Boxplot Legend**: 🔴 **Red line** = Median (50th percentile) &nbsp;|&nbsp; 🔷 **Blue diamond** = Mean (arithmetic average). For skewed single-cell distributions, the Mean is pulled by high-scoring cells.")
                    
                    # 1. Sample Mean & Median Score Summary Tables
                    st.markdown(f"#### {selected_sig_name} Sample Score Summary Tables")
                    score_mean_rows = []
                    score_median_rows = []
                    for pop in score_plots_to_generate:
                        df_pop = df_score_violin if pop == "All Cells (Global)" else df_score_violin[df_score_violin["Cell State"] == pop]
                        r_mean = {"Population / Cluster": pop}
                        r_med = {"Population / Cluster": pop}
                        for s in available_samples_s:
                            s_vals = df_pop[df_pop["Sample"] == s]["Score"].values
                            r_mean[s] = round(float(np.mean(s_vals)), 3) if len(s_vals) > 0 else np.nan
                            r_med[s] = round(float(np.median(s_vals)), 3) if len(s_vals) > 0 else np.nan
                        score_mean_rows.append(r_mean)
                        score_median_rows.append(r_med)
                        
                    df_score_sample_means = pd.DataFrame(score_mean_rows).set_index("Population / Cluster")
                    df_score_sample_medians = pd.DataFrame(score_median_rows).set_index("Population / Cluster")
                    
                    tab_sm_tbl, tab_smed_tbl = st.tabs(["Mean Score Table", "Median Score Table"])
                    with tab_sm_tbl:
                        st.dataframe(df_score_sample_means, use_container_width=True)
                        csv_score_means = df_score_sample_means.to_csv().encode('utf-8')
                        st.download_button(
                            label=f"📥 Download {selected_sig_name} Sample Means CSV",
                            data=csv_score_means,
                            file_name=f"{selected_dataset_name}_{selected_sig_name.replace(' ', '_')}_sample_means.csv",
                            mime="text/csv",
                            key="dl_score_means_csv"
                        )
                    with tab_smed_tbl:
                        st.dataframe(df_score_sample_medians, use_container_width=True)
                        csv_score_medians = df_score_sample_medians.to_csv().encode('utf-8')
                        st.download_button(
                            label=f"📥 Download {selected_sig_name} Sample Medians CSV",
                            data=csv_score_medians,
                            file_name=f"{selected_dataset_name}_{selected_sig_name.replace(' ', '_')}_sample_medians.csv",
                            mime="text/csv",
                            key="dl_score_medians_csv"
                        )
                    
                    # 2. Statistical Significance Table
                    if score_stats_records:
                        st.markdown("#### Statistical Significance Summary Table (Mann-Whitney U Test)")
                        st.caption("ℹ️ **Statistical Note**: Mann-Whitney U is a non-parametric rank-sum test evaluating whether S1 scores tend to be systematically higher than S2 (evaluating rank/median shifts). Both parametric Means and non-parametric Medians are reported for clarity.")
                        df_score_stats = pd.DataFrame(score_stats_records)
                        if len(df_score_stats) > 1:
                            from scipy.stats import false_discovery_control
                            try:
                                df_score_stats["FDR_raw"] = false_discovery_control(df_score_stats["p-value_raw"])
                            except Exception:
                                df_score_stats["FDR_raw"] = df_score_stats["p-value_raw"] * len(df_score_stats)
                        else:
                            df_score_stats["FDR_raw"] = df_score_stats["p-value_raw"]
                            
                        # Format 4 sig figs / scientific notation
                        df_score_stats["p-value"] = df_score_stats["p-value_raw"].apply(format_sig_value)
                        df_score_stats["FDR (q-value)"] = df_score_stats["FDR_raw"].apply(format_sig_value)
                        
                        cols_to_show_s = [
                            "Signature", "Population", "Comparison (S1 vs S2)", 
                            "N (S1)", "Mean Score (S1)", "Median Score (S1)", "N (S2)", "Mean Score (S2)", "Median Score (S2)", 
                            "Mean Diff (S1 - S2)", "Median Diff (S1 - S2)", "Mann-Whitney U", "p-value", "FDR (q-value)", "Significance"
                        ]
                        df_score_stats_display = df_score_stats[cols_to_show_s]
                        st.dataframe(df_score_stats_display, use_container_width=True)
                        csv_score_stats = df_score_stats_display.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label=f"📥 Download {selected_sig_name} Scoring Stats CSV",
                            data=csv_score_stats,
                            file_name=f"{selected_dataset_name}_{selected_sig_name.replace(' ', '_')}_score_stats.csv",
                            mime="text/csv",
                            key="dl_score_stats_csv"
                        )

    # ---------------- TAB 6: CORRELATION & SCATTER PLOTS (WITH STATS) ----------------
    with tab_scatter:
        st.markdown("### Co-expression & Correlation Scatter Plots (with Statistical Testing)")
        st.write("Analyze co-expression relationships between any two genes, gene signatures/pathways, or a gene vs. a pathway score across conditions and cell populations.")
        
        # Dual Variable Selectors (Gene vs Score)
        col_x_panel, col_y_panel = st.columns(2)
        
        # X-Axis Selector
        with col_x_panel:
            st.markdown("#### 🔹 X-Axis Variable")
            x_type = st.radio("X-Axis Type:", ["Gene Expression", "Signature / Pathway Score"], horizontal=True, key="scat_x_type")
            if x_type == "Gene Expression":
                x_marker_cat = st.selectbox("Filter Suggested Genes by Cell Type (X):", ["All Genes"] + list(dataset_markers.keys()), index=0, key="scat_x_ct_filter")
                if x_marker_cat != "All Genes":
                    x_gene_opts = ["None"] + [sym_to_display.get(g.upper(), g) for g in dataset_markers.get(x_marker_cat, []) if g.upper() in sym_to_display]
                else:
                    x_gene_opts = display_options
                selected_x_var_disp = st.selectbox("Select X-Axis Gene:", options=x_gene_opts, index=0, key="scat_x_gene_select")
                x_resolved_var = display_to_var.get(selected_x_var_disp, None) if selected_x_var_disp != "None" else None
                x_label_name = selected_x_var_disp.split(" (")[0] if selected_x_var_disp != "None" else None
                x_is_score = False
            else:
                selected_x_sig = st.selectbox("Select X-Axis Signature:", options=sig_options, index=0, key="scat_x_sig_select")
                if selected_x_sig == "Custom Gene List":
                    x_cust_genes = st.text_input("Custom X genes (comma-separated):", value="Cdh1, Ctnnb1, Col17a1", key="scat_x_cust")
                    x_sig_genes = [g.strip() for g in x_cust_genes.split(",") if g.strip()]
                else:
                    x_sig_genes = DEFAULT_SIGNATURES[selected_x_sig]
                x_resolved_var = [resolve_gene_var_name(adata, g, sym_to_display, display_to_var) for g in x_sig_genes]
                x_resolved_var = [v for v in x_resolved_var if v is not None]
                x_label_name = f"{selected_x_sig} Score"
                x_is_score = True

        # Y-Axis Selector
        with col_y_panel:
            st.markdown("#### 🔸 Y-Axis Variable")
            y_type = st.radio("Y-Axis Type:", ["Gene Expression", "Signature / Pathway Score"], horizontal=True, key="scat_y_type")
            if y_type == "Gene Expression":
                y_marker_cat = st.selectbox("Filter Suggested Genes by Cell Type (Y):", ["All Genes"] + list(dataset_markers.keys()), index=0, key="scat_y_ct_filter")
                if y_marker_cat != "All Genes":
                    y_gene_opts = ["None"] + [sym_to_display.get(g.upper(), g) for g in dataset_markers.get(y_marker_cat, []) if g.upper() in sym_to_display]
                else:
                    y_gene_opts = display_options
                selected_y_var_disp = st.selectbox("Select Y-Axis Gene:", options=y_gene_opts, index=0, key="scat_y_gene_select")
                y_resolved_var = display_to_var.get(selected_y_var_disp, None) if selected_y_var_disp != "None" else None
                y_label_name = selected_y_var_disp.split(" (")[0] if selected_y_var_disp != "None" else None
                y_is_score = False
            else:
                selected_y_sig = st.selectbox("Select Y-Axis Signature:", options=sig_options, index=1 if len(sig_options)>1 else 0, key="scat_y_sig_select")
                if selected_y_sig == "Custom Gene List":
                    y_cust_genes = st.text_input("Custom Y genes (comma-separated):", value="Mki67, Top2a, Pcna", key="scat_y_cust")
                    y_sig_genes = [g.strip() for g in y_cust_genes.split(",") if g.strip()]
                else:
                    y_sig_genes = DEFAULT_SIGNATURES[selected_y_sig]
                y_resolved_var = [resolve_gene_var_name(adata, g, sym_to_display, display_to_var) for g in y_sig_genes]
                y_resolved_var = [v for v in y_resolved_var if v is not None]
                y_label_name = f"{selected_y_sig} Score"
                y_is_score = True

        # Check if both variables are selected
        x_ready = (x_is_score and len(x_resolved_var) > 0) or (not x_is_score and x_resolved_var is not None)
        y_ready = (y_is_score and len(y_resolved_var) > 0) or (not y_is_score and y_resolved_var is not None)
        
        if not x_ready or not y_ready:
            st.info("💡 Please select both an X-axis and a Y-axis variable (gene or score) above to generate correlation scatter plots and statistics.")
        else:
            # Extract data values
            with st.spinner("Extracting expression/score data for correlation..."):
                # X values
                if x_is_score:
                    x_data = compute_score(adata, x_resolved_var, f"scat_x_{selected_x_sig[:15]}", selected_dataset_name)
                else:
                    if scipy.sparse.issparse(adata.X):
                        x_raw = adata[:, x_resolved_var].X.toarray().flatten()
                    else:
                        x_raw = adata[:, x_resolved_var].X.flatten()
                    x_data = np.log2(x_raw + 1) if use_log2 else x_raw

                # Y values
                if y_is_score:
                    y_data = compute_score(adata, y_resolved_var, f"scat_y_{selected_y_sig[:15]}", selected_dataset_name)
                else:
                    if scipy.sparse.issparse(adata.X):
                        y_raw = adata[:, y_resolved_var].X.toarray().flatten()
                    else:
                        y_raw = adata[:, y_resolved_var].X.flatten()
                    y_data = np.log2(y_raw + 1) if use_log2 else y_raw

                df_scatter_full = pd.DataFrame({
                    "X": x_data,
                    "Y": y_data,
                    "Sample": adata.obs[sample_col].astype(str) if sample_col else "All",
                    "Cell State": adata.obs[selected_col].astype(str) if selected_col else "All"
                })

            all_scat_samples = ordered_samples if ordered_samples else sorted(df_scatter_full["Sample"].unique())
            all_scat_states = sorted(df_scatter_full["Cell State"].unique())

            with st.expander("Scatter Plot Display & Subsetting Options", expanded=True):
                c_sc1, c_sc2, c_sc3, c_sc4 = st.columns([1.5, 1.5, 1.2, 1.2])
                with c_sc1:
                    filter_scat_samples = st.multiselect("Filter Samples:", options=all_scat_samples, default=all_scat_samples, key="scat_filter_s")
                with c_sc2:
                    filter_scat_states = st.multiselect("Filter Cell States / Populations:", options=all_scat_states, default=all_scat_states, key="scat_filter_st")
                with c_sc3:
                    split_mode = st.selectbox("Subplot Layout:", ["Split by Sample", "Split by Cell State", "Single Combined Overlay"], index=0, key="scat_split_mode")
                with c_sc4:
                    scat_ncols = st.selectbox("Grid Columns:", [2, 3, 4], index=1, key="scat_ncols_choice")
                    
                c_sc5, c_sc6, c_sc7, c_sc8 = st.columns([1.8, 1.0, 1.0, 1.0])
                with c_sc5:
                    zero_filter_mode = st.selectbox(
                        "Zero-Expression Filtering (Remove zeros):",
                        [
                            "Include all cells (Keep zeros)",
                            "Co-detected only: Remove cells with X = 0 OR Y = 0 (X > 0 and Y > 0)",
                            "Remove double-zeros (X > 0 or Y > 0)",
                            "Remove X = 0 cells only (X > 0)",
                            "Remove Y = 0 cells only (Y > 0)"
                        ],
                        index=0,
                        help="Filter out unexpressed / dropout cells across the selected axes.",
                        key="scat_zero_filter"
                    )
                with c_sc6:
                    show_regline = st.checkbox("Show Trendline", value=True, key="scat_reg_line")
                with c_sc7:
                    scat_pt_size = st.slider("Point Size:", min_value=1.0, max_value=8.0, value=2.5, step=0.5, key="scat_pt_size")
                with c_sc8:
                    scat_alpha = st.slider("Point Opacity:", min_value=0.1, max_value=1.0, value=0.6, step=0.05, key="scat_pt_alpha")

            # Filter data
            df_filtered = df_scatter_full[df_scatter_full["Sample"].isin(filter_scat_samples) & df_scatter_full["Cell State"].isin(filter_scat_states)]
            if zero_filter_mode == "Co-detected only: Remove cells with X = 0 OR Y = 0 (X > 0 and Y > 0)":
                df_filtered = df_filtered[(df_filtered["X"] > 0) & (df_filtered["Y"] > 0)]
            elif zero_filter_mode == "Remove double-zeros (X > 0 or Y > 0)":
                df_filtered = df_filtered[(df_filtered["X"] > 0) | (df_filtered["Y"] > 0)]
            elif zero_filter_mode == "Remove X = 0 cells only (X > 0)":
                df_filtered = df_filtered[df_filtered["X"] > 0]
            elif zero_filter_mode == "Remove Y = 0 cells only (Y > 0)":
                df_filtered = df_filtered[df_filtered["Y"] > 0]

            if df_filtered.empty:
                st.warning("No cells match the chosen filter criteria.")
            else:
                with st.spinner("Calculating correlations and rendering scatter plots..."):
                    # Prepare subplots
                    if split_mode == "Split by Sample":
                        groups_to_plot = ["★ All Selected Cells (Global)"] + [s for s in filter_scat_samples if s in df_filtered["Sample"].values]
                    elif split_mode == "Split by Cell State":
                        groups_to_plot = ["★ All Selected Cells (Global)"] + [st for st in filter_scat_states if st in df_filtered["Cell State"].values]
                    else:
                        groups_to_plot = ["★ All Selected Cells (Global)"]

                    n_scat_plots = len(groups_to_plot)
                    n_sc_cols = min(scat_ncols, n_scat_plots)
                    n_sc_rows = (n_scat_plots + n_sc_cols - 1) // n_sc_cols

                    fig_scat, axes_scat = plt.subplots(n_sc_rows, n_sc_cols, figsize=(4.8 * n_sc_cols, 4.3 * n_sc_rows + 0.6))
                    axes_scat_flat = axes_scat.flatten() if hasattr(axes_scat, 'flatten') else [axes_scat]

                    x_scale_tag = f" ({chosen_scale_label})" if not x_is_score else ""
                    y_scale_tag = f" ({chosen_scale_label})" if not y_is_score else ""
                    
                    # Cell state and zero-filtering status tags for super title
                    if split_mode == "Split by Cell State":
                        cell_state_title_tag = "Split by Cell State"
                    elif len(filter_scat_states) == len(all_scat_states):
                        cell_state_title_tag = "All Cell States"
                    elif len(filter_scat_states) == 1:
                        cell_state_title_tag = f"Cell State: {filter_scat_states[0]}"
                    elif len(filter_scat_states) <= 3:
                        cell_state_title_tag = f"Cell States: {', '.join(filter_scat_states)}"
                    else:
                        cell_state_title_tag = f"{len(filter_scat_states)} Cell States"

                    zero_tag_map = {
                        "Include all cells (Keep zeros)": "All Cells (Zeros Included)",
                        "Co-detected only: Remove cells with X = 0 OR Y = 0 (X > 0 and Y > 0)": "Co-detected Only (X > 0 & Y > 0)",
                        "Remove double-zeros (X > 0 or Y > 0)": "Double-zeros Removed (X > 0 | Y > 0)",
                        "Remove X = 0 cells only (X > 0)": "X > 0 Only",
                        "Remove Y = 0 cells only (Y > 0)": "Y > 0 Only"
                    }
                    zero_filter_title_tag = zero_tag_map.get(zero_filter_mode, zero_filter_mode)

                    fig_scat.suptitle(
                        f"Correlation: {x_label_name} vs. {y_label_name}\n({cell_state_title_tag}  |  {zero_filter_title_tag})", 
                        fontsize=13.0, fontweight='bold', y=0.99
                    )

                    scat_stats_records = []

                    for idx, grp in enumerate(groups_to_plot):
                        ax = axes_scat_flat[idx]
                        
                        if grp.startswith("★"):
                            df_grp = df_filtered
                            grp_title = f"★ All Cells (N={len(df_grp):,})"
                            grp_color = "#2c3e50"
                        elif split_mode == "Split by Sample":
                            df_grp = df_filtered[df_filtered["Sample"] == grp]
                            grp_title = f"{grp} (N={len(df_grp):,})"
                            grp_color = sample_color_map.get(grp, "#3498db")
                        else:  # Split by Cell State
                            df_grp = df_filtered[df_filtered["Cell State"] == grp]
                            grp_title = f"{grp[:24]} (N={len(df_grp):,})"
                            grp_color = "#16a085"

                        x_vals = df_grp["X"].values
                        y_vals = df_grp["Y"].values

                        # Compute Statistics
                        n_cells = len(x_vals)
                        if n_cells >= 3 and np.std(x_vals) > 0 and np.std(y_vals) > 0:
                            try:
                                sp_res = spearmanr(x_vals, y_vals)
                                sp_rho, sp_p = float(sp_res.statistic), float(sp_res.pvalue)
                            except Exception:
                                sp_rho, sp_p = np.nan, 1.0
                            try:
                                pe_res = pearsonr(x_vals, y_vals)
                                pe_r, pe_p = float(pe_res.statistic), float(pe_res.pvalue)
                            except Exception:
                                pe_r, pe_p = np.nan, 1.0
                        else:
                            sp_rho, sp_p, pe_r, pe_p = np.nan, 1.0, np.nan, 1.0

                        sig_label_sp = get_sig_label(sp_p)

                        scat_stats_records.append({
                            "Group / Subplot": grp,
                            "X Variable": f"{x_label_name}{x_scale_tag}",
                            "Y Variable": f"{y_label_name}{y_scale_tag}",
                            "N Cells": n_cells,
                            "Spearman rho": round(sp_rho, 4) if pd.notna(sp_rho) else np.nan,
                            "Spearman p-value_raw": sp_p,
                            "Pearson r": round(pe_r, 4) if pd.notna(pe_r) else np.nan,
                            "Pearson p-value_raw": pe_p,
                            "Significance": sig_label_sp
                        })

                        # Scatter Plot
                        if split_mode == "Single Combined Overlay":
                            for s in filter_scat_samples:
                                sub_s_mask = df_grp["Sample"] == s
                                if sub_s_mask.any():
                                    ax.scatter(
                                        df_grp.loc[sub_s_mask, "X"], 
                                        df_grp.loc[sub_s_mask, "Y"], 
                                        color=sample_color_map.get(s, "#7f8c8d"), 
                                        s=scat_pt_size, 
                                        alpha=scat_alpha, 
                                        label=s
                                    )
                            ax.legend(title="Sample", markerscale=4, fontsize=9.5, frameon=False, loc="best")
                        else:
                            ax.scatter(x_vals, y_vals, color=grp_color, s=scat_pt_size, alpha=scat_alpha)

                        # Regression Trendline
                        if show_regline and n_cells >= 3 and np.std(x_vals) > 0 and np.std(y_vals) > 0:
                            try:
                                poly_fit = np.polyfit(x_vals, y_vals, 1)
                                poly_fn = np.poly1d(poly_fit)
                                x_line = np.linspace(x_vals.min(), x_vals.max(), 100)
                                ax.plot(x_line, poly_fn(x_line), color="red" if split_mode != "Single Combined Overlay" else "black", lw=1.5, ls="--")
                            except Exception:
                                pass

                        # Correlation text on plot
                        rho_txt = f"Spearman ρ = {sp_rho:+.3f}" if pd.notna(sp_rho) else "Spearman ρ = N/A"
                        p_txt = f"p = {sp_p:.2e}" if sp_p < 0.001 else f"p = {sp_p:.3f}"
                        pe_txt = f"Pearson r = {pe_r:+.3f}" if pd.notna(pe_r) else "Pearson r = N/A"
                        stat_annotation = f"{rho_txt} ({p_txt}, {sig_label_sp})\n{pe_txt}"
                        ax.text(0.04, 0.95, stat_annotation, transform=ax.transAxes, fontsize=10, verticalalignment='top',
                                bbox=dict(boxstyle='round,pad=0.35', facecolor='white', alpha=0.85, edgecolor='#bdc3c7'))

                        ax.set_title(grp_title, fontsize=11.5, fontweight='bold')
                        ax.set_xlabel(f"{x_label_name}{x_scale_tag}", fontsize=11)
                        ax.set_ylabel(f"{y_label_name}{y_scale_tag}", fontsize=11)
                        ax.tick_params(axis='both', which='major', labelsize=10)
                        ax.grid(True, linestyle=':', alpha=0.5)

                    for extra_idx in range(n_scat_plots, len(axes_scat_flat)):
                        axes_scat_flat[extra_idx].axis('off')

                    plt.tight_layout(rect=[0, 0, 1, 0.94])
                    st.pyplot(fig_scat)
                    plt.close(fig_scat)

                    # Correlation Statistics Table
                    if scat_stats_records:
                        st.markdown("#### Correlation Statistical Summary Table")
                        st.caption("ℹ️ Computes both non-parametric **Spearman rank correlation** ($\rho$) and parametric **Pearson linear correlation** ($r$) with exact p-values and Benjamini-Hochberg FDR q-values.")
                        df_scat_stats = pd.DataFrame(scat_stats_records)
                        
                        if len(df_scat_stats) > 1:
                            from scipy.stats import false_discovery_control
                            try:
                                df_scat_stats["FDR_raw"] = false_discovery_control(df_scat_stats["Spearman p-value_raw"].fillna(1.0))
                            except Exception:
                                df_scat_stats["FDR_raw"] = df_scat_stats["Spearman p-value_raw"] * len(df_scat_stats)
                        else:
                            df_scat_stats["FDR_raw"] = df_scat_stats["Spearman p-value_raw"]

                        df_scat_stats["Spearman p-value"] = df_scat_stats["Spearman p-value_raw"].apply(format_sig_value)
                        df_scat_stats["Pearson p-value"] = df_scat_stats["Pearson p-value_raw"].apply(format_sig_value)
                        df_scat_stats["FDR (q-value)"] = df_scat_stats["FDR_raw"].apply(format_sig_value)

                        cols_to_display_scat = [
                            "Group / Subplot", "X Variable", "Y Variable", "N Cells", 
                            "Spearman rho", "Spearman p-value", "Pearson r", "Pearson p-value", 
                            "FDR (q-value)", "Significance"
                        ]
                        df_scat_display = df_scat_stats[cols_to_display_scat]
                        st.dataframe(df_scat_display, use_container_width=True)
                        
                        csv_scat = df_scat_display.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label="📥 Download Correlation Summary CSV",
                            data=csv_scat,
                            file_name=f"{selected_dataset_name}_{x_label_name}_vs_{y_label_name}_correlation_stats.csv",
                            mime="text/csv"
                        )

# ----------------- PAGE 2: MARKER EDITOR -----------------
elif app_mode == "Cell-Type Marker Editor":
    st.title("Cell-Type Population & Marker Gene Editor")
    st.write("View, add, modify, or delete cell populations and their key marker genes. Changes will be saved directly to the config YAML.")
    
    if "markers_config" not in st.session_state:
        st.session_state.markers_config = load_yaml()
        
    config = st.session_state.markers_config
    if yaml_key not in config:
        config[yaml_key] = {}
        
    class_config = config[yaml_key]
    
    if selected_col:
        file_categories = sorted(adata.obs[selected_col].cat.categories.tolist()) if hasattr(adata.obs[selected_col], 'cat') else sorted(adata.obs[selected_col].dropna().unique().tolist())
    else:
        file_categories = []
        
    st.subheader(f"Configure Markers for Dataset: {selected_dataset_name}")
    st.write(f"Annotation Column: `{selected_col}` | Found `{len(file_categories)}` categories in this file.")
    
    if len(file_categories) == 0:
        st.warning("No categories found. Please select an annotation column in the sidebar.")
    else:
        st.write("### Edit Markers per Cell State:")
        with st.form("marker_editor_form"):
            updated_class_config = {}
            cols = st.columns(2)
            for idx, cat in enumerate(file_categories):
                with cols[idx % 2]:
                    existing_markers = class_config.get(cat, [])
                    if not existing_markers and yaml_key in DEFAULT_MARKERS:
                        existing_markers = DEFAULT_MARKERS[yaml_key].get(cat, [])
                    marker_str = st.text_input(f"Markers for '{cat}':", value=", ".join(existing_markers))
                    parsed_markers = [g.strip() for g in marker_str.split(",") if g.strip()]
                    updated_class_config[cat] = parsed_markers
                    
            save_btn = st.form_submit_button("Save Configuration Changes")
            if save_btn:
                config[yaml_key] = updated_class_config
                st.session_state.markers_config = config
                save_yaml(config)
                st.success("Successfully updated and saved config YAML!")
                st.rerun()
                
        st.write("---")
        if st.button("Reset Markers to Default List"):
            if yaml_key in DEFAULT_MARKERS:
                config[yaml_key] = DEFAULT_MARKERS[yaml_key].copy()
                st.session_state.markers_config = config
                save_yaml(config)
                st.success("Reset markers to defaults.")
                st.rerun()
