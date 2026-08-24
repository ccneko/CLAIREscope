import os
import io
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
try:
    from streamlit_sortables import sort_items
except Exception:
    def sort_items(items, **kwargs):
        return items

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
PROJ_BASE = os.path.dirname(os.path.dirname(APP_DIR))
DATA_DIR = os.path.join(PROJ_BASE, "out", "2026-08-23_human_epidermal_subclustering")
SCAN_DIRS = [
    os.path.join(PROJ_BASE, "out", "2026-08-23_human_epidermal_subclustering"),
    os.path.join(PROJ_BASE, "out", "2026-07-06_adherens_junction_and_col17a1_correlation_analysis"),
    os.path.join(PROJ_BASE, "out", "2026-02-09_celltypist"),
    os.path.join(PROJ_BASE, "out", "2026-02-26_harmony"),
    os.path.join(PROJ_BASE, "data", "2025-03-18_Public_SingleCell")
]
YAML_PATH = os.path.join(APP_DIR, "cell_type_markers.yaml")

st.set_page_config(page_title="PROJ_001 JEB Expression Viewer & Marker Editor", layout="wide")

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
    iframe[title="streamlit.components.v1.html"], div[data-testid="stCustomComponentV1"] {
        display: none !important;
        height: 0 !important;
        margin: 0 !important;
        padding: 0 !important;
    }
</style>
""", unsafe_allow_html=True)

# Disable Streamlit 'C' / 'c' Clear Cache keyboard shortcut popup in browser
components.html("""
<script>
(function() {
    try {
        const targetWin = window.parent || window;
        const targetDoc = targetWin.document;
        targetWin.addEventListener('keydown', function(e) {
            const activeEl = targetDoc.activeElement;
            const isEditable = activeEl && (
                activeEl.tagName === 'INPUT' || 
                activeEl.tagName === 'TEXTAREA' || 
                activeEl.tagName === 'SELECT' || 
                activeEl.isContentEditable
            );
            if (!isEditable) {
                // Block standalone 'c' / 'C' (Streamlit's clear cache shortcut)
                if ((e.key === 'c' || e.key === 'C' || e.keyCode === 67) && !e.ctrlKey && !e.metaKey && !e.altKey) {
                    e.stopImmediatePropagation();
                    e.preventDefault();
                }
            }
        }, true);
    } catch (err) {
        console.warn("Unable to attach shortcut blocker:", err);
    }
})();
</script>
""", height=0, width=0)

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

DATASET_CONFIG_PATH = os.path.join(APP_DIR, "dataset_config.yaml")

# Helper to load dataset config
def load_dataset_config():
    if os.path.exists(DATASET_CONFIG_PATH):
        try:
            with open(DATASET_CONFIG_PATH, "r", encoding="utf-8") as f:
                cfg = yaml.safe_load(f)
                return cfg if isinstance(cfg, dict) else {}
        except Exception:
            return {}
    return {}

# Helper to save dataset config
def save_dataset_config(cfg):
    with open(DATASET_CONFIG_PATH, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, default_flow_style=False, sort_keys=False)

# Scan folder dynamically for h5ad files
def scan_datasets(data_dir):
    active_datasets = {}
    all_datasets = {}
    target_dirs = SCAN_DIRS if 'SCAN_DIRS' in globals() else [data_dir]
    cfg = load_dataset_config()
    hidden_list = cfg.get("hidden_datasets", [])
    
    for d in target_dirs:
        if not os.path.exists(d):
            continue
        for filename in sorted(os.listdir(d)):
            if filename.endswith(".h5ad") and not filename.startswith("."):
                filepath = os.path.join(d, filename)
                folder_tag = os.path.basename(d)
                ds_name = f"{filename[:-5]} ({folder_tag})"
                all_datasets[ds_name] = filepath
                if ds_name not in hidden_list:
                    active_datasets[ds_name] = filepath
    return active_datasets, all_datasets

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
    "Cell-Type Marker Editor",
    "Dataset Management & Launch Settings"
])

# Dynamic Dataset Picker
active_datasets, all_detected_datasets = scan_datasets(DATA_DIR)
if not active_datasets:
    if all_detected_datasets:
        st.sidebar.warning("All datasets are currently marked hidden in Dataset Settings.")
        active_datasets = all_detected_datasets
    else:
        st.error("No `.h5ad` files found in specified directories.")
        st.stop()

cfg = load_dataset_config()
default_ds_name = cfg.get("default_dataset", None)
active_keys = list(active_datasets.keys())
default_ds_idx = active_keys.index(default_ds_name) if default_ds_name in active_keys else 0

selected_dataset_name = st.sidebar.selectbox("Select Dataset:", active_keys, index=default_ds_idx)
h5ad_path = active_datasets[selected_dataset_name]

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
canonical_samples = ["Control", "Rescued_1", "Rescued_2", "Mutant", "Control_P4", "Sample_Rescued_1", "Sample_Rescued_2", "Sample_Mutant", "Normal", "JEB", "Revertant"]
if sample_col and sample_col in adata.obs.columns:
    unique_in_data = adata.obs[sample_col].dropna().unique().tolist()
    ordered_samples = [s for s in canonical_samples if s in unique_in_data] + [s for s in unique_in_data if s not in canonical_samples]
else:
    ordered_samples = []

sample_color_map = {"Control": "#e74c3c", "Rescued_1": "#8e44ad", "Rescued_2": "#f1c40f", "Mutant": "#00a8ff", "Control_P4": "#e74c3c", "Sample_Rescued_1": "#8e44ad", "Sample_Rescued_2": "#f1c40f", "Sample_Mutant": "#00a8ff", "Normal": "#2ecc71", "JEB": "#e74c3c", "Revertant": "#3498db"}
for idx, s in enumerate(ordered_samples):
    if s not in sample_color_map:
        cmap = plt.get_cmap('tab10')
        sample_color_map[s] = matplotlib.colors.to_hex(cmap(idx % 10))

# ----------------- PAGE 1: EXPRESSION VIEWER & ANALYSIS TABS -----------------
if app_mode == "Gene Expression UMAP":
    st.title("Single-Cell RNA-seq Expression, Scoring & Correlation Viewer")
    st.markdown(f'<div style="font-size: 20px; font-weight: 500; margin-top: 4px; margin-bottom: 18px; color: #1e293b; line-height: 1.5;">Active Dataset: <code>{selected_dataset_name}</code> | Total Cells: <code>{adata.n_obs:,}</code> | Total Genes: <code>{adata.n_vars:,}</code> | Sample Column: <code>{sample_col}</code></div>', unsafe_allow_html=True)
    
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
    if resolved_var_name and resolved_var_name not in adata.var_names:
        # Fallback search if var_names format differs
        if selected_gene_box in adata.var_names:
            resolved_var_name = selected_gene_box
        else:
            resolved_var_name = None
    resolved_display_name = selected_gene_box if resolved_var_name else None
    
    # Calculate raw expression values & percentiles in background
    if resolved_var_name:
        if scipy.sparse.issparse(adata.X):
            raw_vals = adata[:, resolved_var_name].X.toarray().flatten()
        else:
            raw_vals = adata[:, resolved_var_name].X.flatten()
        
        raw_log2_vals = np.log2(raw_vals + 1)
        max_possible_log2 = float(raw_log2_vals.max()) if len(raw_log2_vals) > 0 else 10.0
        max_possible_lin = float(raw_vals.max()) if len(raw_vals) > 0 else 100.0
    else:
        raw_vals = np.array([0.0])
        raw_log2_vals = np.array([0.0])
        max_possible_log2, max_possible_lin = 10.0, 100.0
        
    use_log2 = True
    cmap_choice = "viridis"
    pct_slider = 100
    chosen_vmax = max_possible_log2 if max_possible_log2 > 0 else 1.0
    chosen_scale_label = "Log2(Norm+1)"
    
    use_log2_t2 = True
    cmap_choice_t2 = "viridis"
    pct_slider_t2 = 100
    chosen_vmax_t2 = max_possible_log2 if max_possible_log2 > 0 else 1.0
    chosen_scale_label_t2 = "Log2(Norm+1)" 

    def rank_cell_state(c):
        c_low = str(c).lower()
        if "basal 1" in c_low or "basal_1" in c_low:
            return (0, c)
        elif "basal 2" in c_low or "basal_2" in c_low:
            return (1, c)
        elif "spinous" in c_low:
            return (2, c)
        elif "granular" in c_low:
            return (3, c)
        elif "mitotic" in c_low:
            return (4, c)
        elif "channel" in c_low:
            return (5, c)
        elif "wnti" in c_low or "bulge" in c_low:
            return (6, c)
        elif "follicular" in c_low or "sebaceous" in c_low:
            return (7, c)
        else:
            return (8, str(c))

    def get_cluster_color_map(adata_obj, col_name):
        if not col_name or col_name not in adata_obj.obs.columns:
            return {}, []
        if hasattr(adata_obj.obs[col_name], 'cat'):
            categories = adata_obj.obs[col_name].cat.categories.tolist()
        else:
            categories = sorted(adata_obj.obs[col_name].dropna().unique().tolist())
            
        if any("basal" in str(c).lower() or "spinous" in str(c).lower() or "granular" in str(c).lower() for c in categories):
            categories = sorted(categories, key=rank_cell_state)
        
        color_key = f"{col_name}_colors"
        if col_name == 'predicted_labels':
            pred_map = {}
            for cat in categories:
                c_str = str(cat).lower()
                if 'undiff' in c_str:
                    pred_map[cat] = '#c0c0c0'  # Silver for Undifferentiated
                elif 'diff' in c_str:
                    pred_map[cat] = '#f5deb3'  # Wheat for Differentiated
                else:
                    pred_map[cat] = '#c0c0c0'
            colors = [pred_map[c] for c in categories]
            try:
                adata_obj.uns[color_key] = colors
            except Exception:
                pass
            return pred_map, categories

        if any("basal" in str(c).lower() or "spinous" in str(c).lower() or "granular" in str(c).lower() for c in categories):
            fixed_palette = {
                'Basal 1': '#1f77b4',
                'Basal 2': '#ff7f0e',
                'Spinous': '#f1c40f',
                'Granular': '#2ca02c',
                'Mitotic': '#d62728',
                'Channel': '#9467bd',
                'Bulge': '#8c564b',
                'Follicular': '#e377c2'
            }
            colors = []
            for cat in categories:
                c_low = str(cat).lower()
                matched = False
                for k, v in fixed_palette.items():
                    if k.lower() in c_low:
                        colors.append(v)
                        matched = True
                        break
                if not matched:
                    colors.append('#7f8c8d')
            try:
                adata_obj.uns[color_key] = colors
            except Exception:
                pass
            return dict(zip(categories, colors)), categories

        if color_key in adata_obj.uns and len(adata_obj.uns[color_key]) >= len(categories):
            colors = list(adata_obj.uns[color_key])[:len(categories)]
        else:
            cmap = plt.get_cmap('tab20')
            colors = [matplotlib.colors.to_hex(cmap(i % 20)) for i in range(len(categories))]
            try:
                adata_obj.uns[color_key] = colors
            except Exception:
                pass
                
        return dict(zip(categories, colors)), categories

    # Caching static multi-panel grid generator with Sample UMAP as 1st plot and custom vmax
    @st.cache_data
    def generate_static_grid(_adata, var_key, disp_title, col, s_col, dataset_name, log2_mode, v_max, cmap_name, grid_cols=3, grid_rows="Auto", col_color_tag=""):
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
            
        x_min, x_max = float(np.min(umap_coords[:, 0])), float(np.max(umap_coords[:, 0]))
        y_min, y_max = float(np.min(umap_coords[:, 1])), float(np.max(umap_coords[:, 1]))
        x_pad = (x_max - x_min) * 0.05
        y_pad = (y_max - y_min) * 0.05
        u_xlim = (x_min - x_pad, x_max + x_pad)
        u_ylim = (y_min - y_pad, y_max + y_pad)

        samples_list = [s for s in ordered_samples if s in _adata.obs[s_col].unique()][:4] if s_col else []
        num_splits = len(samples_list)
        total_plots = 3 + num_splits
        n_cols = int(grid_cols)
        if grid_rows == "Auto" or grid_rows is None:
            n_rows = (total_plots + n_cols - 1) // n_cols
        else:
            n_rows = int(grid_rows)
        
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
            ax_samp.set_xlim(u_xlim)
            ax_samp.set_ylim(u_ylim)
            ax_samp.set_aspect("equal", adjustable="box")
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
            ax_ct.set_xlim(u_xlim)
            ax_ct.set_ylim(u_ylim)
            ax_ct.set_aspect("equal", adjustable="box")
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
        ax_all.set_xlim(u_xlim)
        ax_all.set_ylim(u_ylim)
        ax_all.set_aspect("equal", adjustable="box")
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
            ax_sub.set_xlim(u_xlim)
            ax_sub.set_ylim(u_ylim)
            ax_sub.set_aspect("equal", adjustable="box")
            ax_sub.tick_params(axis='both', which='major', labelsize=9.5)
            cbar_s = fig.colorbar(sc_sub, ax=ax_sub)
            cbar_s.ax.tick_params(labelsize=9.5)
            
        for extra_idx in range(3 + num_splits, len(axes_flat)):
            axes_flat[extra_idx].axis('off')
            
        plt.tight_layout()
        return fig

    @st.cache_data
    def get_umap_embeddings_csv(_adata, s_col, a_col, var_name, disp_name):
        if 'X_umap' not in _adata.obsm:
            return None
        df_out = pd.DataFrame({
            "Barcode": _adata.obs_names,
            "UMAP_1": _adata.obsm['X_umap'][:, 0],
            "UMAP_2": _adata.obsm['X_umap'][:, 1]
        })
        if s_col and s_col in _adata.obs:
            df_out["Sample"] = _adata.obs[s_col].values
        if a_col and a_col in _adata.obs:
            df_out[a_col] = _adata.obs[a_col].values
        if var_name and var_name in _adata.var_names:
            if scipy.sparse.issparse(_adata.X):
                g_raw = _adata[:, var_name].X.toarray().flatten()
            else:
                g_raw = _adata[:, var_name].X.flatten()
            clean_sym = disp_name.split(" (")[0] if disp_name else "Gene"
            df_out[f"Expression_{clean_sym}_Raw"] = g_raw
            df_out[f"Expression_{clean_sym}_Log2"] = np.log2(g_raw + 1)
        return df_out.to_csv(index=False).encode('utf-8')

    # 7 Main Tabs
    tab_static, tab_interactive, tab_composition, tab_gene_violin, tab_score_violin, tab_scatter, tab_trajectory = st.tabs([
        "Static UMAP", 
        "Interactive UMAP", 
        "Sample Composition",
        "Gene Expression Violins",
        "Signature & Pathway Scoring",
        "Correlation & Scatter Plots",
        "Trajectory Analysis"
    ])
    
    # ---------------- TAB 1: STATIC UMAP ----------------
    with tab_static:
        with st.expander("🎨 Colormap, Scale & Contrast Controls", expanded=bool(resolved_var_name)):
            c_scale, c_cmap, c_pct, c_vmax = st.columns([1.2, 1.0, 1.8, 1.0])
            with c_scale:
                use_log2 = st.checkbox("Log2(Normalized + 1) Scale", value=True, help="Applies log2 transformation like Loupe Browser for balanced contrast.", key="tab1_use_log2")
            with c_cmap:
                cmap_choice = st.selectbox("Colormap:", ["viridis", "YlOrRd", "Reds", "inferno", "plasma", "magma", "turbo"], index=0, key="tab1_cmap_choice")
                
            expr_for_scale = raw_log2_vals if use_log2 else raw_vals
            max_possible = max_possible_log2 if use_log2 else max_possible_lin
            
            with c_pct:
                pct_slider = st.select_slider(
                    "Max Percentile Threshold (Anchors):",
                    options=[50, 60, 70, 75, 80, 85, 90, 95, 98, 99, 99.5, 100],
                    value=100,
                    format_func=lambda x: f"{x}%" if x in [80, 90, 95, 99, 100] else f"{x}",
                    help="Clip upper colormap limit to enhance visual contrast against outliers.",
                    key="tab1_vmax_pct_slider"
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
                    help="Direct numeric limit for colormap maximum.",
                    key="tab1_custom_vmax"
                )
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
                color_tag_str = str(get_cluster_color_map(adata, selected_col)[0]) if selected_col else ''
                fig_grid = generate_static_grid(adata, resolved_var_name, resolved_display_name, selected_col, sample_col, selected_dataset_name, use_log2, chosen_vmax, cmap_choice, grid_cols=stat_grid_cols, grid_rows=stat_grid_rows, col_color_tag=color_tag_str)
                st.pyplot(fig_grid)
                
                svg_grid_buf = io.BytesIO()
                fig_grid.savefig(svg_grid_buf, format="svg", bbox_inches="tight")
                clean_sym_name = resolved_display_name.split(" (")[0] if resolved_display_name else "gene"
                
                c_dl1, c_dl2, _ = st.columns([0.28, 0.28, 0.44])
                with c_dl1:
                    csv_data_t1 = get_umap_embeddings_csv(adata, sample_col, selected_col, resolved_var_name, resolved_display_name)
                    if csv_data_t1:
                        st.download_button(
                            label="📥 Download UMAP Embeddings (CSV)",
                            data=csv_data_t1,
                            file_name=f"{selected_dataset_name}_umap_embeddings.csv",
                            mime="text/csv",
                            key="dl_tab1_umap_csv"
                        )
                with c_dl2:
                    st.download_button(
                        label="📥 Download Grid Plot as SVG",
                        data=svg_grid_buf.getvalue(),
                        file_name=f"{selected_dataset_name}_{clean_sym_name}_static_grid.svg",
                        mime="image/svg+xml",
                        key="dl_tab1_grid_svg"
                    )
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
                        svg_s_buf = io.BytesIO()
                        fig_s.savefig(svg_s_buf, format="svg", bbox_inches="tight")
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
                            svg_c_buf = io.BytesIO()
                            fig_ref.savefig(svg_c_buf, format="svg", bbox_inches="tight")
                            plt.close(fig_ref)

                    c_r_dl1, c_r_dl2, c_r_dl3, _ = st.columns([0.28, 0.25, 0.25, 0.22])
                    with c_r_dl1:
                        csv_data_t1 = get_umap_embeddings_csv(adata, sample_col, selected_col, resolved_var_name, resolved_display_name)
                        if csv_data_t1:
                            st.download_button(
                                label="📥 Download UMAP Embeddings (CSV)",
                                data=csv_data_t1,
                                file_name=f"{selected_dataset_name}_umap_embeddings.csv",
                                mime="text/csv",
                                key="dl_tab1_umap_csv_ref"
                            )
                    with c_r_dl2:
                        st.download_button(
                            label="📥 Download Sample UMAP (SVG)",
                            data=svg_s_buf.getvalue(),
                            file_name=f"{selected_dataset_name}_sample_umap.svg",
                            mime="image/svg+xml",
                            key="dl_tab1_sample_svg"
                        )
                    with c_r_dl3:
                        if selected_col:
                            st.download_button(
                                label="📥 Download Cell State UMAP (SVG)",
                                data=svg_c_buf.getvalue(),
                                file_name=f"{selected_dataset_name}_cell_state_umap.svg",
                                mime="image/svg+xml",
                                key="dl_tab1_cellstate_svg"
                            )
            
    # ---------------- TAB 2: INTERACTIVE UMAP ----------------
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
                    selected_samples = draggable_multiselect(
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
                    expr_sub = np.log2(expr_sub_raw + 1) if use_log2_t2 else expr_sub_raw
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
                        
                    fig_samp.update_xaxes(range=x_range, title_text="UMAP 1", title_font=dict(size=16, family="Segoe UI, sans-serif"), tickfont=dict(size=14, family="Segoe UI, sans-serif"), zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_samp.update_yaxes(range=y_range, title_text="UMAP 2", title_font=dict(size=16, family="Segoe UI, sans-serif"), tickfont=dict(size=14, family="Segoe UI, sans-serif"), scaleanchor="x", scaleratio=1, zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_samp.update_layout(
                        height=540,
                        title_font=dict(size=18, family="Segoe UI, sans-serif"),
                        margin=dict(l=10, r=10, t=50, b=10),
                        legend=dict(
                            itemsizing='constant',
                            font=dict(size=14, family="Segoe UI, sans-serif"),
                            title=dict(font=dict(size=15, family="Segoe UI, sans-serif")),
                            bgcolor="rgba(255,255,255,0.9)",
                            bordercolor="#CBD5E1",
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
                        
                    fig_states.update_xaxes(range=x_range, title_text="UMAP 1", title_font=dict(size=16, family="Segoe UI, sans-serif"), tickfont=dict(size=14, family="Segoe UI, sans-serif"), zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_states.update_yaxes(range=y_range, title_text="UMAP 2", title_font=dict(size=16, family="Segoe UI, sans-serif"), tickfont=dict(size=14, family="Segoe UI, sans-serif"), scaleanchor="x", scaleratio=1, zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                    fig_states.update_layout(
                        height=540,
                        title_font=dict(size=18, family="Segoe UI, sans-serif"),
                        margin=dict(l=10, r=10, t=50, b=10),
                        legend=dict(
                            itemsizing='constant',
                            font=dict(size=14, family="Segoe UI, sans-serif"),
                            title=dict(font=dict(size=15, family="Segoe UI, sans-serif")),
                            bgcolor="rgba(255,255,255,0.9)",
                            bordercolor="#CBD5E1",
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
                                range_color=[0, chosen_vmax_t2],
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
                                        color=np.clip(df_selected["Expression"], 0, chosen_vmax_t2),
                                        colorscale=plotly_cs,
                                        colorbar=dict(title=f"{chosen_scale_label}", len=0.8, thickness=14, tickfont=dict(size=10)),
                                        size=pt_size,
                                        opacity=0.85,
                                        cmin=0,
                                        cmax=chosen_vmax_t2,
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
                                range_color=[0, chosen_vmax_t2],
                                title=f"{clean_sym} ({chosen_scale_label}, All Cells)",
                                template="plotly_white"
                            )
                            fig_expr.update_traces(marker=dict(size=pt_size, opacity=0.85))
                            
                        fig_expr.update_xaxes(range=x_range, title_text="UMAP 1", title_font=dict(size=16, family="Segoe UI, sans-serif"), tickfont=dict(size=14, family="Segoe UI, sans-serif"), zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                        fig_expr.update_yaxes(range=y_range, title_text="UMAP 2", title_font=dict(size=16, family="Segoe UI, sans-serif"), tickfont=dict(size=14, family="Segoe UI, sans-serif"), scaleanchor="x", scaleratio=1, zeroline=False, showgrid=True, gridcolor="#F8FAFC")
                        fig_expr.update_layout(
                            height=540,
                            title_font=dict(size=18, family="Segoe UI, sans-serif"),
                            margin=dict(l=10, r=10, t=50, b=10),
                            legend=dict(itemsizing='constant', font=dict(size=14, family="Segoe UI, sans-serif")),
                            coloraxis_colorbar=dict(title_font=dict(size=15, family="Segoe UI, sans-serif"), tickfont=dict(size=13, family="Segoe UI, sans-serif"))
                        )
                        st.plotly_chart(fig_expr, use_container_width=True)

            # Download UMAP Embeddings (CSV) at the bottom of Interactive UMAP
            csv_data_t2 = get_umap_embeddings_csv(adata, sample_col, selected_col, resolved_var_name, resolved_display_name)
            if csv_data_t2:
                st.markdown("<div style='margin-top: 15px;'></div>", unsafe_allow_html=True)
                c_t2_dl, _ = st.columns([0.28, 0.72])
                with c_t2_dl:
                    st.download_button(
                        label="📥 Download UMAP Embeddings (CSV)",
                        data=csv_data_t2,
                        file_name=f"{selected_dataset_name}_umap_embeddings.csv",
                        mime="text/csv",
                        key="dl_tab2_umap_csv"
                    )

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
                    selected_comp_samples = draggable_multiselect("Select & Reorder Samples to Display:", options=all_dataset_samples, default=all_dataset_samples, key="comp_samples_filter")
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
                        height=540,
                        margin=dict(l=10, r=10, t=35, b=10),
                        xaxis=dict(title="Sample / Condition", title_font=dict(size=18, family="Segoe UI, sans-serif"), tickfont=dict(size=16, family="Segoe UI, sans-serif")),
                        yaxis=dict(title=yaxis_title, title_font=dict(size=18, family="Segoe UI, sans-serif"), tickfont=dict(size=15, family="Segoe UI, sans-serif"), range=yaxis_range, showgrid=True, gridcolor="#F1F5F9"),
                        legend=dict(
                            itemsizing='constant',
                            font=dict(size=14, family="Segoe UI, sans-serif"),
                            title=dict(font=dict(size=15, family="Segoe UI, sans-serif")),
                            bgcolor="rgba(255,255,255,0.9)",
                            bordercolor="#CBD5E1",
                            borderwidth=1
                        )
                    )
                    st.plotly_chart(fig_bar, use_container_width=True)
                    
                with col_chart2:
                    st.markdown("#### Sample Percentage Ring (Donut) Charts")
                    num_s = len(selected_comp_samples)
                    if num_s <= 4:
                        n_rows, n_cols = 1, num_s
                        donut_height = 520
                    else:
                        n_cols = min(num_s, 3)
                        n_rows = (num_s + n_cols - 1) // n_cols
                        donut_height = 270 * n_rows
                        
                    donut_specs = [[{"type": "domain"} for _ in range(n_cols)] for _ in range(n_rows)]
                    subplot_titles = [f"<b style='font-size: 20px; color: #1e293b;'>{s}</b><br><span style='font-size: 16px; color: #475569;'>Total: {sample_totals[s]:,} cells</span>" for s in selected_comp_samples]
                    
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
                                hole=0.48,
                                marker=dict(
                                    colors=[color_map.get(c, "#7f8c8d") for c in categories],
                                    line=dict(color='#FFFFFF', width=1.5)
                                ),
                                textinfo='percent' if show_donut_pct else 'none',
                                textposition='inside',
                                insidetextfont=dict(size=18, family="Segoe UI, Arial, sans-serif", color="#FFFFFF"),
                                textfont=dict(size=18, family="Segoe UI, Arial, sans-serif"),
                                insidetextorientation='horizontal',
                                hovertemplate='<b>%{label}</b><br>Sample: ' + s + '<br>Count: %{value:,} cells<br>Percentage: %{percent}<extra></extra>',
                                showlegend=False,
                                sort=False
                            ),
                            row=r_idx, col=c_idx
                        )
                        
                    fig_donuts.update_layout(
                        template='plotly_white',
                        height=donut_height,
                        margin=dict(l=10, r=10, t=65, b=10),
                        uniformtext_minsize=14,
                        uniformtext_mode='hide',
                        hoverlabel=dict(font_size=14)
                    )
                    for annotation in fig_donuts['layout']['annotations']:
                        annotation['font'] = dict(size=16, family="Segoe UI, Arial, sans-serif", color="#1e293b")
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
                    selected_states_v = draggable_multiselect("Select & Reorder Cell States to Include:", options=all_states, default=all_states, key="v_gene_states")
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
                        
                default_active_pairs = [("Control", "Mutant"), ("Rescued_2", "Mutant"), ("Rescued_1", "Mutant")]
                valid_default_pairs = [p for p in default_active_pairs if p[0] in available_samples and p[1] in available_samples]
                if not valid_default_pairs and len(candidate_pairs) > 0:
                    valid_default_pairs = candidate_pairs[:2]
                    
                pair_options = [f"{p[0]} vs {p[1]}" for p in candidate_pairs]
                default_pair_strs = [f"{p[0]} vs {p[1]}" for p in valid_default_pairs]
                
                selected_pair_strs = draggable_multiselect("Select Comparison Pairs for Significance Brackets (S1 vs S2):", options=pair_options, default=default_pair_strs, key="v_gene_pairs")
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
                                
                                if len(val1) < 3 or len(val2) < 3:
                                    # Do not show comparison lines when cells are missing / insufficient
                                    continue
                                
                                p_val = 1.0
                                u_stat = np.nan
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
                    selected_states_sv = draggable_multiselect("Select & Reorder Cell States for Scoring Violins:", options=all_states_s, default=all_states_s, key="sv_states")
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
                        
                default_active_pairs_s = [("Control", "Mutant"), ("Rescued_2", "Mutant"), ("Rescued_1", "Mutant")]
                valid_default_pairs_s = [p for p in default_active_pairs_s if p[0] in available_samples_s and p[1] in available_samples_s]
                if not valid_default_pairs_s and len(candidate_pairs_s) > 0:
                    valid_default_pairs_s = candidate_pairs_s[:2]
                    
                pair_options_s = [f"{p[0]} vs {p[1]}" for p in candidate_pairs_s]
                default_pair_strs_s = [f"{p[0]} vs {p[1]}" for p in valid_default_pairs_s]
                
                selected_pair_strs_s = draggable_multiselect("Select Comparison Pairs for Significance Brackets (S1 vs S2):", options=pair_options_s, default=default_pair_strs_s, key="sv_pairs")
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
                                
                                if len(val1) < 3 or len(val2) < 3:
                                    # Do not show comparison lines when cells are missing / insufficient
                                    continue
                                
                                p_val = 1.0
                                u_stat = np.nan
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
                    filter_scat_samples = draggable_multiselect("Filter & Reorder Samples:", options=all_scat_samples, default=all_scat_samples, key="scat_filter_s")
                with c_sc2:
                    filter_scat_states = draggable_multiselect("Filter & Reorder Cell States / Populations:", options=all_scat_states, default=all_scat_states, key="scat_filter_st")
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


    # ---------------- TAB 7: TRAJECTORY ANALYSIS ----------------
    with tab_trajectory:
        st.markdown("### Lineage Trajectory & Pseudotime Dynamics")
        st.write("Explore continuous differentiation trajectories, PAGA connectivity, and dynamic gene & pathway expression kinetics along developmental pseudotime.")
        
        # Detect pseudotime columns in adata.obs
        pt_cols = [c for c in adata.obs.columns if any(k in c.lower() for k in ['pseudotime', 'dpt', 'slingshot', 'latent_time', 'trajectory', 'time'])]
        
        if not pt_cols:
            st.warning("⚠️ No pre-computed pseudotime or trajectory variables found in this dataset (e.g. `dpt_pseudotime`, `dpt_landmark_aligned`, `latent_time`).")
            st.info("💡 You can compute Diffusion Pseudotime (DPT) on-the-fly using the button below.")
            
            if st.button("🚀 Compute Diffusion Pseudotime (DPT) on root Basal Cells", key="btn_compute_dpt"):
                with st.spinner("Computing neighborhood graph, diffusion map, and DPT (root=Basal)..."):
                    try:
                        if 'neighbors' not in adata.uns:
                            sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
                        if 'diffmap' not in adata.obsm:
                            sc.tl.diffmap(adata)
                        
                        # Select root cell in Basal cluster if available
                        root_idx = 0
                        if selected_col and selected_col in adata.obs.columns:
                            basal_mask = adata.obs[selected_col].astype(str).str.contains('Basal', case=False, na=False)
                            if basal_mask.sum() > 0:
                                root_idx = np.where(basal_mask)[0][0]
                        adata.uns['iroot'] = root_idx
                        sc.tl.dpt(adata)
                        st.success("✅ Diffusion Pseudotime computed successfully!")
                        st.rerun()
                    except Exception as err:
                        st.error(f"Error computing DPT: {err}")
        else:
            def format_pseudotime_label(pt_col):
                if pt_col == "dpt_landmark_aligned":
                    return "Pseudotime (Landmark Aligned)"
                elif pt_col == "dpt_sample_minmax":
                    return "Pseudotime (Within-Sample Min-Max [0, 1])"
                elif pt_col == "dpt_pseudotime":
                    return "Diffusion Pseudotime (Original Global DPT)"
                elif "aligned" in pt_col.lower():
                    return f"Pseudotime (Aligned: {pt_col})"
                elif "minmax" in pt_col.lower() or "norm" in pt_col.lower():
                    return f"Pseudotime (Normalized: {pt_col})"
                else:
                    return f"Pseudotime ({pt_col})"

            # Default to dpt_pseudotime (Original Global DPT) as requested
            pref_pt_idx = 0
            if "dpt_pseudotime" in pt_cols:
                pref_pt_idx = pt_cols.index("dpt_pseudotime")
            elif "dpt_landmark_aligned" in pt_cols:
                pref_pt_idx = pt_cols.index("dpt_landmark_aligned")

            selected_pt_col = st.selectbox(
                "Select Pseudotime Variable:",
                pt_cols,
                index=pref_pt_idx,
                format_func=format_pseudotime_label,
                key="traj_pt_col"
            )
            pt_label_str = format_pseudotime_label(selected_pt_col)
            
            with st.expander("Trajectory Overview Controls", expanded=True):
                c_tr1, c_tr2, c_tr3, c_tr4, c_tr5 = st.columns([1.2, 1.0, 0.9, 0.9, 1.1])
                with c_tr1:
                    tr_cmap = st.selectbox("Pseudotime Colormap:", ["plasma", "viridis", "inferno", "magma", "cividis", "turbo"], index=0, key="tr_cmap")
                with c_tr2:
                    tr_pt_size = st.slider("Trajectory Point Size:", min_value=0.5, max_value=6.0, value=2.0, step=0.5, key="tr_pt_size")
                with c_tr3:
                    tr_grid_cols = st.selectbox("Grid Columns:", [1, 2, 3, 4, 5, 6], index=2, key="tr_grid_cols")
                with c_tr4:
                    tr_grid_rows = st.selectbox("Grid Rows:", ["Auto", 1, 2, 3, 4, 5, 6], index=0, key="tr_grid_rows")
                with c_tr5:
                    show_densities = st.checkbox("Show Density by Sample", value=True, key="tr_show_dens")
            
            # 1. Trajectory Overview Multi-panel Grid
            if 'X_umap' in adata.obsm:
                with st.spinner("Rendering trajectory overview grid..."):
                    umap_xy = adata.obsm['X_umap']
                    pt_vals = adata.obs[selected_pt_col].values
                    tx_min, tx_max = float(np.min(umap_xy[:, 0])), float(np.max(umap_xy[:, 0]))
                    ty_min, ty_max = float(np.min(umap_xy[:, 1])), float(np.max(umap_xy[:, 1]))
                    tx_pad = (tx_max - tx_min) * 0.05
                    ty_pad = (ty_max - ty_min) * 0.05
                    tu_xlim = (tx_min - tx_pad, tx_max + tx_pad)
                    tu_ylim = (ty_min - ty_pad, ty_max + ty_pad)
                    
                    samples_list = [s for s in ordered_samples if s in adata.obs[sample_col].unique()][:4] if (sample_col and sample_col in adata.obs.columns) else []
                    n_sub_samples = len(samples_list)
                    
                    total_tr_plots = 2 + (1 if show_densities else 0) + n_sub_samples
                    n_tr_cols = int(tr_grid_cols)
                    if tr_grid_rows == "Auto" or tr_grid_rows is None:
                        n_tr_rows = (total_tr_plots + n_tr_cols - 1) // n_tr_cols
                    else:
                        n_tr_rows = int(tr_grid_rows)
                    
                    fig_tr, axes_tr = plt.subplots(n_tr_rows, n_tr_cols, figsize=(5.2 * n_tr_cols, 4.6 * n_tr_rows))
                    axes_tr_flat = axes_tr.flatten() if hasattr(axes_tr, 'flatten') else [axes_tr]
                    
                    # Plot 1: Continuous Pseudotime
                    ax1 = axes_tr_flat[0]
                    sort_pt = np.argsort(pt_vals)
                    sc_pt = ax1.scatter(umap_xy[sort_pt, 0], umap_xy[sort_pt, 1], c=pt_vals[sort_pt], cmap=tr_cmap, s=tr_pt_size, alpha=0.85)
                    ax1.set_title(f"{pt_label_str}", fontsize=11.5, fontweight='bold')
                    ax1.set_xlabel("UMAP 1", fontsize=10)
                    ax1.set_ylabel("UMAP 2", fontsize=10)
                    ax1.set_xlim(tu_xlim)
                    ax1.set_ylim(tu_ylim)
                    ax1.set_aspect("equal", adjustable="box")
                    ax1.tick_params(axis='both', which='major', labelsize=9.5)
                    cb_pt = fig_tr.colorbar(sc_pt, ax=ax1, label=pt_label_str)
                    cb_pt.ax.tick_params(labelsize=9)
                    
                    # Plot 2: Cell States
                    ax2 = axes_tr_flat[1]
                    if selected_col:
                        color_map, categories = get_cluster_color_map(adata, selected_col)
                        for cat in categories:
                            mask = adata.obs[selected_col] == cat
                            ax2.scatter(umap_xy[mask, 0], umap_xy[mask, 1], label=cat, color=color_map.get(cat, "#7f8c8d"), s=tr_pt_size, alpha=0.8)
                        ax2.set_title(f"Cell States ({selected_col})", fontsize=11.5, fontweight='bold')
                        ax2.set_xlabel("UMAP 1", fontsize=10)
                        ax2.set_ylabel("UMAP 2", fontsize=10)
                        ax2.set_xlim(tu_xlim)
                        ax2.set_ylim(tu_ylim)
                        ax2.set_aspect("equal", adjustable="box")
                        ax2.tick_params(axis='both', which='major', labelsize=9.5)
                        ax2.legend(title="Cell State", bbox_to_anchor=(0.5, -0.2), loc="upper center", markerscale=5, fontsize=8.5, ncol=2, frameon=False)
                    else:
                        ax2.text(0.5, 0.5, "No Annotation Column", ha='center', va='center')
                        ax2.axis('off')
                        
                    curr_plot_idx = 2
                    
                    # Plot 3: Pseudotime Distribution by Sample
                    if show_densities and sample_col and sample_col in adata.obs.columns:
                        ax3 = axes_tr_flat[curr_plot_idx]
                        curr_plot_idx += 1
                        for s in ordered_samples:
                            if s in adata.obs[sample_col].values:
                                s_pt = adata.obs[adata.obs[sample_col] == s][selected_pt_col].dropna().values
                                if len(s_pt) > 5:
                                    sns.kdeplot(s_pt, ax=ax3, label=s, color=sample_color_map.get(s, "#3498db"), lw=2.0)
                        ax3.set_title(f"Distribution by Sample", fontsize=11.5, fontweight='bold')
                        ax3.set_xlabel(pt_label_str, fontsize=10, fontweight='bold')
                        ax3.set_ylabel("Density", fontsize=10)
                        ax3.tick_params(axis='both', which='major', labelsize=9.5)
                        ax3.legend(title="Sample", fontsize=9, frameon=False)
                        ax3.grid(True, linestyle=':', alpha=0.5)
                        
                    # Split Subplots per Sample
                    for idx, sample in enumerate(samples_list):
                        if curr_plot_idx < len(axes_tr_flat):
                            ax_s = axes_tr_flat[curr_plot_idx]
                            curr_plot_idx += 1
                            mask = adata.obs[sample_col] == sample
                            bg_x = umap_xy[~mask, 0]
                            bg_y = umap_xy[~mask, 1]
                            ax_s.scatter(bg_x, bg_y, color='lightgrey', s=0.5, alpha=0.25)
                            
                            sub_pt = pt_vals[mask]
                            sub_coords = umap_xy[mask]
                            sub_sort = np.argsort(sub_pt)
                            sc_sub = ax_s.scatter(sub_coords[sub_sort, 0], sub_coords[sub_sort, 1], c=sub_pt[sub_sort], cmap=tr_cmap, s=tr_pt_size, alpha=0.85, vmin=0, vmax=1)
                            ax_s.set_title(f"{sample} ({selected_pt_col})", fontsize=11.5, fontweight='bold')
                            ax_s.set_xlabel("UMAP 1", fontsize=10)
                            ax_s.set_ylabel("UMAP 2", fontsize=10)
                            ax_s.set_xlim(tu_xlim)
                            ax_s.set_ylim(tu_ylim)
                            ax_s.set_aspect("equal", adjustable="box")
                            ax_s.tick_params(axis='both', which='major', labelsize=9.5)
                            cb_s = fig_tr.colorbar(sc_sub, ax=ax_s)
                            cb_s.ax.tick_params(labelsize=9)
                        
                    for extra in range(curr_plot_idx, len(axes_tr_flat)):
                        axes_tr_flat[extra].axis('off')
                        
                    plt.tight_layout()
                    st.pyplot(fig_tr)
                    plt.close(fig_tr)
            
            # 2. Dynamic Feature Kinetics along Pseudotime (Unified Genes & Pathway Scores Studio)
            st.markdown("---")
            st.markdown("### 📈 Dynamic Gene & Pathway Kinetics along Pseudotime")
            st.caption(f"ℹ️ Plot continuous gene expression and pathway/signature score kinetics along **{pt_label_str}** across cohorts.")
            
            # Build unified feature options: Scores + Signatures + Genes
            obs_scores = [c for c in adata.obs.columns if (c.startswith('score_') or 'score' in c.lower()) and pd.api.types.is_numeric_dtype(adata.obs[c])]
            score_labels = [f"Score: {c.replace('score_', '').replace('_', ' ')}" for c in obs_scores]
            score_col_map = dict(zip(score_labels, obs_scores))
            
            sig_labels = [f"Sig: {s}" for s in DEFAULT_SIGNATURES.keys()] if 'DEFAULT_SIGNATURES' in globals() or 'DEFAULT_SIGNATURES' in locals() else []
            gene_labels = display_options[1:]
            
            all_feature_options = score_labels + sig_labels + gene_labels
            
            # Helper function to extract array and labels for any selected feature (gene or score or signature)
            def get_feature_kinetics_data(feat_name):
                if feat_name in score_col_map:
                    col_name = score_col_map[feat_name]
                    vals = adata.obs[col_name].values
                    clean_title = feat_name.replace("Score: ", "") + " Score"
                    y_lbl = f"{clean_title} (Module Score)"
                    return vals, clean_title, y_lbl
                elif feat_name.startswith("Sig: "):
                    sig_name = feat_name.replace("Sig: ", "")
                    sig_genes = DEFAULT_SIGNATURES.get(sig_name, [])
                    sig_vars = [resolve_gene_var_name(adata, g, sym_to_display, display_to_var) for g in sig_genes]
                    sig_vars = [v for v in sig_vars if v and v in adata.var_names]
                    if sig_vars:
                        if scipy.sparse.issparse(adata.X):
                            sub_m = adata[:, sig_vars].X.toarray()
                        else:
                            sub_m = adata[:, sig_vars].X
                        vals = np.mean(sub_m, axis=1)
                    else:
                        vals = np.zeros(adata.n_obs)
                    clean_title = sig_name + " Signature"
                    y_lbl = f"{clean_title} (Mean Expr)"
                    return vals, clean_title, y_lbl
                else:
                    # Gene
                    var_g = resolve_gene_var_name(adata, feat_name, sym_to_display, display_to_var)
                    clean_g_sym = feat_name.split(" (")[0]
                    if var_g and var_g in adata.var_names:
                        if scipy.sparse.issparse(adata.X):
                            raw_g = adata[:, var_g].X.toarray().flatten()
                        else:
                            raw_g = adata[:, var_g].X.flatten()
                        vals = np.log2(raw_g + 1) if use_log2 else raw_g
                    else:
                        vals = np.zeros(adata.n_obs)
                    clean_title = clean_g_sym
                    y_lbl = f"Expression ({chosen_scale_label})"
                    return vals, clean_title, y_lbl

            with st.expander("⚙️ Dynamic Drawing Controls & Layout", expanded=True):
                c_mode, c_win, c_gcols, c_grows, c_scat = st.columns([1.3, 1.0, 0.8, 0.8, 1.1])
                with c_mode:
                    traj_plot_mode = st.radio("Display Mode:", ["Multi-Feature Kinetics Grid", "Single Feature (Samples & Cell States)"], horizontal=True, key="traj_plot_mode")
                with c_win:
                    traj_smooth_window = st.slider("Smoothing Window (Cells):", min_value=10, max_value=400, value=80, step=10, key="traj_smooth_win")
                with c_gcols:
                    mg_grid_cols = st.selectbox("Grid Columns:", [1, 2, 3, 4, 5, 6], index=1, key="mg_grid_cols")
                with c_grows:
                    mg_grid_rows = st.selectbox("Grid Rows:", ["Auto", 1, 2, 3, 4, 5, 6], index=0, key="mg_grid_rows")
                with c_scat:
                    traj_show_scatter = st.checkbox("Show Single-Cell Scatter Points", value=False, key="traj_show_scat")
                    
            if traj_plot_mode == "Multi-Feature Kinetics Grid":
                # Default selection (mix of key genes + stratum scores)
                def_genes = [g for g in ["COL17A1", "KRT14", "AREG", "KRT10"] if g in sym_to_display]
                default_disp = [sym_to_display[g] for g in def_genes if g in sym_to_display]
                if score_labels:
                    default_disp.extend(score_labels[:2])
                if not default_disp:
                    default_disp = all_feature_options[:4]
                    
                selected_multi_features = draggable_multiselect(
                    "Select & Reorder Features (Genes and/or Pathway Scores) along Pseudotime:",
                    options=all_feature_options,
                    default=default_disp[:4],
                    key="traj_multi_features"
                )
                
                if selected_multi_features:
                    with st.spinner("Generating multi-feature trajectory kinetics..."):
                        n_mg = len(selected_multi_features)
                        n_cols_mg = int(mg_grid_cols)
                        if mg_grid_rows == "Auto" or mg_grid_rows is None:
                            n_rows_mg = (n_mg + n_cols_mg - 1) // n_cols_mg
                        else:
                            n_rows_mg = int(mg_grid_rows)
                        
                        fig_mg, axes_mg = plt.subplots(n_rows_mg, n_cols_mg, figsize=(6.8 * n_cols_mg, 4.6 * n_rows_mg), dpi=200)
                        axes_mg_flat = axes_mg.flatten() if hasattr(axes_mg, 'flatten') else [axes_mg]
                        
                        for idx_g, feat_name in enumerate(selected_multi_features):
                            if idx_g < len(axes_mg_flat):
                                ax_g = axes_mg_flat[idx_g]
                                feat_vals, clean_feat_title, feat_y_lbl = get_feature_kinetics_data(feat_name)
                                
                                df_g = pd.DataFrame({
                                    "Pseudotime": adata.obs[selected_pt_col].values,
                                    "Value": feat_vals,
                                    "Sample": adata.obs[sample_col].astype(str) if (sample_col and sample_col in adata.obs.columns) else "All"
                                }).dropna(subset=["Pseudotime"]).sort_values("Pseudotime")
                                
                                if sample_col and sample_col in adata.obs.columns:
                                    for s in ordered_samples:
                                        if s in df_g["Sample"].values:
                                            df_gs = df_g[df_g["Sample"] == s]
                                            if len(df_gs) > 10:
                                                s_col = sample_color_map.get(s, "#3498db")
                                                if traj_show_scatter:
                                                    ax_g.scatter(df_gs["Pseudotime"], df_gs["Value"], color=s_col, s=1.2, alpha=0.2)
                                                num_gs = df_gs[['Pseudotime', 'Value']].copy()
                                                roll_win = max(min(traj_smooth_window, len(df_gs)//2), 5)
                                                df_gs_roll = num_gs.rolling(window=roll_win, min_periods=5)['Value'].mean()
                                                ax_g.plot(df_gs["Pseudotime"], df_gs_roll, label=s, color=s_col, lw=2.8)
                                                
                                ax_g.set_title(f"{clean_feat_title} along Pseudotime", fontsize=13, fontweight='bold')
                                ax_g.set_xlabel(pt_label_str, fontsize=11, fontweight='bold')
                                ax_g.set_ylabel(feat_y_lbl, fontsize=11)
                                ax_g.tick_params(axis='both', which='major', labelsize=10)
                                ax_g.grid(True, linestyle=':', alpha=0.5)
                                ax_g.legend(title="Cohort", fontsize=9.5, frameon=True, facecolor='white', framealpha=0.9)
                                
                        for extra in range(n_mg, len(axes_mg_flat)):
                            axes_mg_flat[extra].axis('off')
                            
                        plt.tight_layout()
                        st.pyplot(fig_mg)
                        plt.close(fig_mg)
                else:
                    st.info("Please select at least one gene or pathway score to display.")
                    
            else:
                # Single-feature detailed dual view (works for any gene or score)
                c_sg1, c_sg2 = st.columns([2, 1])
                with c_sg1:
                    traj_feat_choice = st.selectbox(
                        "Select Gene or Pathway Score for Trajectory Plot:", 
                        all_feature_options, 
                        index=all_feature_options.index(st.session_state.selected_gene_display) if st.session_state.selected_gene_display in all_feature_options else 0,
                        key="traj_single_feat"
                    )
                    
                with st.spinner("Plotting feature dynamics along pseudotime..."):
                    single_vals, single_title, single_y_lbl = get_feature_kinetics_data(traj_feat_choice)
                    
                    df_traj = pd.DataFrame({
                        "Pseudotime": adata.obs[selected_pt_col].values,
                        "Value": single_vals,
                        "Sample": adata.obs[sample_col].astype(str) if (sample_col and sample_col in adata.obs.columns) else "All",
                        "Cell State": adata.obs[selected_col].astype(str) if (selected_col and selected_col in adata.obs.columns) else "All"
                    }).dropna(subset=["Pseudotime"]).sort_values("Pseudotime")
                    
                    fig_dyn, (ax_dyn1, ax_dyn2) = plt.subplots(1, 2, figsize=(14, 5.0), dpi=200)
                    
                    # Subplot 1: Split by Condition / Sample
                    if sample_col and sample_col in adata.obs.columns:
                        for s in ordered_samples:
                            if s in df_traj["Sample"].values:
                                df_s = df_traj[df_traj["Sample"] == s]
                                if len(df_s) > 10:
                                    s_col = sample_color_map.get(s, "#3498db")
                                    if traj_show_scatter:
                                        ax_dyn1.scatter(df_s["Pseudotime"], df_s["Value"], color=s_col, s=1.5, alpha=0.25)
                                    num_s = df_s[['Pseudotime', 'Value']].copy()
                                    roll_win = max(min(traj_smooth_window, len(df_s)//2), 5)
                                    df_s_roll = num_s.rolling(window=roll_win, min_periods=5)['Value'].mean()
                                    ax_dyn1.plot(df_s["Pseudotime"], df_s_roll, label=s, color=s_col, lw=2.8)
                    ax_dyn1.set_title(f"{single_title} Kinetics by Sample", fontsize=12.5, fontweight='bold')
                    ax_dyn1.set_xlabel(pt_label_str, fontsize=11, fontweight='bold')
                    ax_dyn1.set_ylabel(single_y_lbl, fontsize=11)
                    ax_dyn1.tick_params(axis='both', which='major', labelsize=10)
                    ax_dyn1.legend(title="Cohort", fontsize=10, frameon=True, facecolor='white', framealpha=0.9)
                    ax_dyn1.grid(True, linestyle=':', alpha=0.5)
                    
                    # Subplot 2: Split by Cell State / Global Trend
                    if selected_col and selected_col in adata.obs.columns:
                        categories = adata.obs[selected_col].cat.categories.tolist() if hasattr(adata.obs[selected_col], "cat") else sorted(adata.obs[selected_col].dropna().unique().tolist())
                        if traj_show_scatter:
                            for cat in categories:
                                df_c = df_traj[df_traj["Cell State"] == cat]
                                if len(df_c) > 0:
                                    ax_dyn2.scatter(df_c["Pseudotime"], df_c["Value"], label=cat, color=color_map.get(cat, "#7f8c8d"), s=2.0, alpha=0.45)
                    num_glob = df_traj[['Pseudotime', 'Value']].copy()
                    df_glob_roll = num_glob.rolling(window=max(traj_smooth_window, 10), min_periods=5)['Value'].mean()
                    ax_dyn2.plot(df_traj["Pseudotime"], df_glob_roll, color="black", lw=3.2, label="Global Trend", ls="--")
                    ax_dyn2.set_title(f"{single_title} Kinetics by Cell State", fontsize=12.5, fontweight='bold')
                    ax_dyn2.set_xlabel(pt_label_str, fontsize=11, fontweight='bold')
                    ax_dyn2.set_ylabel(single_y_lbl, fontsize=11)
                    ax_dyn2.tick_params(axis='both', which='major', labelsize=10)
                    ax_dyn2.legend(title="Cell State / Trend", bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=4, fontsize=9.0, frameon=False)
                    ax_dyn2.grid(True, linestyle=':', alpha=0.5)
                    
                    plt.tight_layout()
                    st.pyplot(fig_dyn)
                    plt.close(fig_dyn)
            
            # 3. Precomputed Pipeline Trajectory Figures Gallery
            data_folder = os.path.dirname(h5ad_path)
            traj_pngs = [f for f in sorted(os.listdir(data_folder)) if f.endswith('.png') and any(k in f.lower() for k in ['trajectory', 'pseudotime', 'paga'])] if os.path.exists(data_folder) else []
            
            # Also check scan dirs for precomputed trajectory pngs
            for d in SCAN_DIRS:
                if os.path.exists(d) and d != data_folder:
                    extra_pngs = [os.path.join(d, f) for f in os.listdir(d) if f.endswith('.png') and any(k in f.lower() for k in ['trajectory', 'pseudotime', 'paga'])]
                    traj_pngs.extend(extra_pngs)
                    
            if traj_pngs:
                st.markdown("---")
                with st.expander("🖼️ View Pre-computed Pipeline Trajectory & PAGA Figures", expanded=False):
                    sel_fig = st.selectbox("Select Precomputed Figure:", traj_pngs)
                    full_fig_path = sel_fig if os.path.isabs(sel_fig) else os.path.join(data_folder, sel_fig)
                    st.image(full_fig_path, caption=os.path.basename(full_fig_path), use_container_width=True)

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

# ----------------- PAGE 3: DATASET MANAGEMENT & SETTINGS -----------------
elif app_mode == "Dataset Management & Launch Settings":
    st.title("📁 Dataset Management & Launch Preferences")
    st.write("Configure which single-cell datasets are active in dropdowns, choose the default dataset to load automatically upon app launch, or hide unneeded `.h5ad` files.")
    
    cfg = load_dataset_config()
    hidden_list = cfg.get("hidden_datasets", [])
    current_default = cfg.get("default_dataset", None)
    all_avail_keys = list(all_detected_datasets.keys())
    
    st.markdown("### 1. Default Dataset at Launch")
    default_idx = all_avail_keys.index(current_default) if current_default in all_avail_keys else 0
    chosen_default = st.selectbox("Select which dataset to load automatically at startup:", all_avail_keys, index=default_idx)
    
    st.markdown("### 2. Dataset Visibility in Dropdown")
    st.caption("ℹ️ Unchecking a dataset hides it from the main selector to keep the dropdown clean and accelerate startup.")
    
    new_hidden = []
    for ds_name in all_avail_keys:
        is_shown = ds_name not in hidden_list
        c_chk, c_info = st.columns([1.5, 3])
        with c_chk:
            show_box = st.checkbox(f"Show `{ds_name}`", value=is_shown, key=f"ds_vis_{ds_name}")
            if not show_box:
                new_hidden.append(ds_name)
        with c_info:
            st.caption(f"📁 Path: `{all_detected_datasets[ds_name]}`")
            
    st.write("")
    if st.button("💾 Save Dataset Preferences & Reload"):
        cfg["default_dataset"] = chosen_default
        cfg["hidden_datasets"] = new_hidden
        save_dataset_config(cfg)
        st.success("Dataset preferences saved successfully! Reloading...")
        st.rerun()
