# 🔬 Single-Cell RNA-seq Expression, Composition, Scoring & Correlation Viewer

An interactive, high-performance [Streamlit](https://streamlit.io/) application designed for comprehensive single-cell RNA sequencing (scRNA-seq / snRNA-seq) exploration, multi-dataset benchmarking, cell composition analytics, pathway scoring, and statistical correlation testing.

---

## 🌟 Key Features

### 1. Multi-Dataset & Dynamic Column Detection
- **Auto-discovery**: Dynamically discovers `.h5ad` AnnData files from local `data/` directories or paths specified via the `SC_DATA_DIR` environment variable.
- **Smart Annotation Detection**: Automatically identifies cell-type annotation columns (e.g. `predicted_labels`, `majority_voting`, `cell_type`, `leiden`, `louvain`) and experimental condition / sample columns.
- **Gene Search & ID Normalization**: Supports search by gene symbol, Ensembl ID, and curated cell-type marker suggestions with bidirectional symbol-to-ID mappings.

---

### 2. Tab 1: Static & Split Condition UMAPs
- **Loupe-Style Contrast Controls**: Linear vs. $\log_2(x+1)$ scaling, custom colormaps (`viridis`, `magma`, `plasma`, `inferno`, `cividis`, `YlOrRd`, `rocket`), text-input $v_{\max}$, and percentile clipping anchors ($80\%$, $90\%$, $95\%$, $99\%$, $100\%$).
- **Multi-Panel Grid**: Leading Sample UMAP $\rightarrow$ Cell State UMAP $\rightarrow$ Global Expression $\rightarrow$ Subplots split by each sample/condition.

---

### 3. Tab 2 & 3: Interactive Cell Exploration & Composition Analytics
- **Tab 2 (Interactive UMAP)**: Plotly-powered hover inspection showing Cell ID, Sample, Cell State, and raw/scaled gene expression.
- **Tab 3 (Cell Type Composition)**: Normalized proportions, absolute cell counts, stacked bar charts, and condition-specific donut distribution diagrams.

---

### 4. Tab 4: Gene Expression Violins & Statistical Tests
- Global ("All Cells") summary violin followed by split subplots per cell state.
- Embedded boxplots (medians, IQR) and pairwise condition significance brackets.
- Preceding **Sample Mean Expression Table** (clusters as rows, samples as columns) with CSV export.
- **Mann-Whitney U Test Summary Table** with two-sided $p$-values and Benjamini-Hochberg FDR $q$-values formatted to 4 significant figures.
- **Zero-Expression Filtering Toggle**: Option to exclude unexpressed cells (count $= 0$) from violin distributions and statistical tests.

---

### 5. Tab 5: Gene Signature & Pathway Scoring
- Calculates continuous gene set scores using `sc.tl.score_genes`.
- Built-in signatures (*Adherens Junctions*, *Desmosomes*, *Hemidesmosomes*, *Cell Cycle / Proliferation*) + **Custom Gene List Input**.
- Uniform Y-axis upper limit controls, sample mean score tables, and Mann-Whitney U statistical significance testing with CSV downloads.

---

### 6. Tab 6: Co-expression & Correlation Scatter Plots (with Stats)
- **Dual-Axis Flexible Selector**: Any combination of **Gene Expression** vs. **Gene Expression**, **Gene** vs. **Score**, or **Score** vs. **Score** on X and Y axes.
- **Subplot Split Modes**: Split by Sample, Split by Cell State, or Single Combined Overlay.
- **Zero-Expression Filtering Modes**:
  - `Include all cells (Keep zeros)`
  - `Co-detected only (X > 0 and Y > 0)`
  - `Remove double-zeros (X > 0 or Y > 0)`
  - `Remove X = 0 cells only (X > 0)`
  - `Remove Y = 0 cells only (Y > 0)`
- Automated non-parametric **Spearman rank correlation** ($\rho, p$) and parametric **Pearson linear correlation** ($r, p$) with linear regression trendlines and exportable summary tables.

---

### 7. Page 2: Cell-Type Marker YAML Editor
- Built-in GUI to view, add, edit, or delete cell populations and their key marker genes saved directly to `cell_type_markers.yaml`.

---

## 🚀 Quick Start

### 1. Clone & Install Dependencies

Using [`uv`](https://github.com/astral-sh/uv) (recommended):
```bash
git clone https://github.com/ccneko/sc-expression-viewer.git
cd sc-expression-viewer
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .
```

Or using `pip`:
```bash
git clone https://github.com/ccneko/sc-expression-viewer.git
cd sc-expression-viewer
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

### 2. Place Your Data Files

Place your `.h5ad` single-cell datasets into a `data/` directory (or set `SC_DATA_DIR`):
```bash
mkdir -p data
cp /path/to/your/dataset.h5ad data/
```

---

### 3. Launch the Viewer

```bash
streamlit run app.py
```
Open your browser at `http://localhost:8501`.

---

## 📁 Repository Structure

```
sc-expression-viewer/
├── app.py                   # Main Streamlit application
├── cell_type_markers.yaml   # Curated cell-type marker dictionary
├── pyproject.toml           # Project metadata and dependencies
├── requirements.txt         # Pip dependency manifest
├── .gitignore               # Strict ignore rules (excludes *.h5ad and datasets)
└── README.md                # Project documentation and guide
```

---

## 🔒 Data Privacy & License

- **Data Privacy**: All single-cell data objects (`*.h5ad`, `*.rds`, `*.loom`, `*.csv`) are strictly ignored via `.gitignore` and are never committed to version control.
- **License**: Private & Proprietary.
