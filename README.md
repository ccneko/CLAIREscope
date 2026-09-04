# CLAIREscope (Cellular Landscape Analysis, Interpretation & Results Explorer)

[![PyPI Version](https://img.shields.io/badge/pypi-v1.0.0-blue.svg)](https://pypi.org/project/clairescope/)
[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/Docs-ReadTheDocs-brightgreen.svg)](https://clairescope.readthedocs.io)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

An open-source, lightweight, and high-performance single-cell analysis platform built on **Python**, **Scanpy**, and **Streamlit**. CLAIREscope unifies end-to-end data preprocessing, live multi-project exploration, dynamic continuous trajectory kinetics modeling, on-the-fly differential expression (Volcano), hypergeometric pathway enrichment (ORA), and automated publication-grade packaging (300 DPI vector SVGs/PDFs and structured summary CSV matrices).

---

## 🚀 Quick Start

### 1. Installation
```bash
# Clone the repository
git clone https://github.com/ccneko/CLAIREscope.git
cd CLAIREscope

# Create virtual environment using uv (recommended)
uv venv --python 3.12
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
uv pip install -r requirements.txt
```

### 2. Launch the Application

#### Option A: Desktop Server Manager (Recommended)
Launch the graphical server controller with one-click Start, Stop, and browser management:

- **Windows (1-Click)**: Double-click `Server Manager.bat` (or `CLAIREscope Server Manager` desktop shortcut)
- **Cross-Platform / CLI**:
  ```bash
  python launcher.py
  # or after pip installation:
  clairescope-gui
  ```
> **Features**: Instant status badges, customizable ports, WSL instance management, and auto-browser launch.

#### Option B: Direct Streamlit CLI
```bash
streamlit run app.py
```
Open your browser and navigate to `http://localhost:8501`.

---

## 🌟 11 Specialized Interactive Analysis Studios & Real-Time Analytics

CLAIREscope provides an integrated suite of 11 high-performance analytical studios powered by real-time in-memory statistical and visualization engines:

1. **🗺️ Static UMAP Explorer (Tab 1)**:
   * **1:1 Isometric Square Geometry**: Mathematically enforced aspect ratio (`ax.set_box_aspect(1)`) eliminating coordinate distortion.
   * **Loupe-Calibrated Dynamic Scaling**: $\log_2(\text{Normalized}+1)$ transformation with upper-percentile clipping (95th/99th percentile) and continuous colormaps (`viridis`, `Reds`, `YlOrRd`, `turbo`, `inferno`).
   * **Multi-Sample Split Grids**: Side-by-side cohort expression comparison with 300 DPI vector SVG and CSV coordinate export.

2. **✨ Interactive Plotly UMAP Studio (Tab 2)**:
   * **Dynamic Cohort Filtering & Dimming**: Isolate active cell states/samples while dimming unselected populations in translucent light gray (`#F0F2F6`).
   * **Cellular Hover Tooltips**: Live cell-level inspection delivering exact expression values, cell state classifications, and donor metadata.

3. **📊 Sample Composition & Stratification (Tab 3)**:
   * **Automated Frequency Shifts**: Instant calculation of cell-state percentages (%) and absolute cell counts across conditions.
   * **Donut & Stacked Bar Charts**: Immediate visualization of lineage shifts and downloadable cross-tabulation CSVs.

4. **🎻 Gene Expression Violins with Statistical Testing (Tab 4)**:
   * **Real-Time Non-Parametric Hypothesis Testing**: Automated computation of **Mann-Whitney $U$ / Wilcoxon Rank-Sum tests** across custom condition pairs.
   * **Automated Significance Brackets**: Dynamic bracket placement (`ns`, `*`, `**`, `***`, `****`).
   * **Persistent Draggable Multiselect**: Reorderable comparison pairs with `⚡ Select all`, `✕ Clear`, and dropout ($E=0$) filtering.

5. **📈 Signature & Pathway Scoring Studio (Tab 5)**:
   * **Dynamic In-Memory Scoring**: `Scanpy` (`sc.tl.score_genes`) evaluation of composite gene modules (Adherens Junctions, Desmosomes, Hemidesmosomes, Cell Cycle G1/S vs G2/M, Custom Gene Lists).
   * **Cross-Condition Statistical Testing**: Automated Mann-Whitney tests comparing signature scores across patient groups per cell population.

6. **📉 Co-Expression & Correlation Scatter Studio (Tab 6)**:
   * **Multi-Variable Bivariate Analysis**: Real-time scatter plots for Gene vs. Gene, Gene vs. Pathway Score, or Score vs. Score.
   * **Dual Metric Computation**: Pearson linear correlation ($r, p$) and Spearman rank correlation ($\rho, p$) with regression trendlines and zero-dropout filtering.

7. **🌿 Trajectory Kinetics & Spline Modeling (Tab 7)**:
   * **Continuous Diffusion Pseudotime ($DPT$) Splines**: Third-order polynomial and B-spline regression curves along continuous developmental trajectories.
   * **State Transition Tracking**: Quantifies stem cell pool maintenance (Basal 1: $DPT < 0.15$), premature state transition dynamics, and phenotypic rescue dynamics.

8. **🌋 Differential Expression & Volcano Studio (Tab 8)**:
   * **Bidirectional Wilcoxon DEG**: Fast computation of $\log_2(\text{Fold Change})$ and Benjamini-Hochberg FDR-adjusted $p$-values between any two cohorts.
   * **Interactive Volcano Plot**: Dynamic significance thresholds ($\text{FDR} < 0.05, |\log_2\text{FC}| > 1.0$) with live top-gene labeling and scientific notation formatting.

9. **🔥 Clustered Heatmap Studio (Tab 9)**:
   * **SciPy $\ge$ 1.18.1 Hierarchical Clustering**: Pairwise Euclidean distance matrices (`pdist`) and average linkage dendrograms.
   * **Draggable Axis Reordering**: Sortable chip interface for customized sample and cell-state ordering with Z-score standardization.

10. **🧬 Pathway Over-Representation Analysis (ORA) Studio (Tab 10)**:
    * **Hypergeometric Enrichment Testing**: Tests for significant pathway over-representation against GO Biological Process, KEGG, and Reactome databases.
    * **Directional Enrichment**: Separates up-regulated and down-regulated gene sets with rich factor and $-\log_{10}(p\text{-adj})$ visualization.

11. **📦 Bulk Packaging & Provenance Studio (Tab 11)**:
    * **Automated Publication Package Generation**: One-click generation of a compressed ZIP bundle containing 300 DPI vector SVGs, high-resolution PNGs, and PDFs for all active figures.
    * **Standardized 1-Row-per-Feature CSV Matrices**: Comprehensive tabular statistical summaries ready for supplementary submission.
    * **Drag-and-Drop Config Importer**: Reusable YAML / CSV / Excel panel import ensuring 100% reproducible analysis provenance.

---

## 💻 Hardware & System Requirements

CLAIREscope is designed to be highly resource-efficient through dynamic on-demand memory management (`@st.cache_resource` and sparse CSR matrix support).

| Component | Minimum Specification (Exploratory / Small Atlases $< 20,000$ cells) | Recommended Specification (Standard Cohorts $20,000 - 100,000+$ cells) | Atlas-Scale / Production Server ($> 200,000$ cells) |
| :--- | :--- | :--- | :--- |
| **Operating System** | Linux (Ubuntu 20.04+), macOS (12+), Windows 10/11 (WSL2 or Native) | Linux (Ubuntu 22.04 / 24.04), macOS (Apple Silicon M1/M2/M3), Windows 11 WSL2 | Linux Server (Ubuntu / Debian / RHEL) / Docker |
| **Processor (CPU)** | Dual-Core (x86_64 or ARM64) | 8+ Cores (e.g. AMD Ryzen 7/9, Intel Core i7/i9, Apple M-Series) | 16–32+ Cores (e.g. Intel Xeon / AMD EPYC) |
| **Memory (RAM)** | **8 GB RAM** | **16 – 32 GB RAM** | **64 – 128+ GB RAM** |
| **Storage (Disk)** | 2 GB free SSD space for dependencies | 10 – 50 GB NVMe SSD (for cached `.h5ad` datasets & vector SVGs) | 100+ GB NVMe SSD |
| **Graphics (GPU)** | Not required (CPU accelerated via Scipy / Numpy) | Optional NVIDIA GPU (CUDA) for rapid UMAP/Harmony acceleration | NVIDIA RTX / A100 / V100 GPU |
| **Python Environment**| Python $\ge 3.10$ ($\le 3.13$) | Python 3.11 or 3.12 managed via `uv` or `conda` | Python 3.11/3.12 with `uv` virtual environment |
| **Web Browser** | Chrome, Firefox, Safari, Edge, Brave (HTML5 + WebSocket support) | Google Chrome, Mozilla Firefox, or Safari | Modern Chromium / WebKit engine |

---

## 📖 Documentation & Tutorials
For complete guides, configuration file specifications, and remote deployment protocols (NordVPN Meshnet, Cloudflare Tunnels), visit our [ReadTheDocs Manual](https://clairescope.readthedocs.io).

## 📖 Citation

If you use CLAIREscope in your research, please cite:

> Chung, C. (2026). *CLAIREscope: Cellular Landscape Analysis, Interpretation & Results Explorer for Single-Cell & Spatial Transcriptomics*. Department of Dermatology, Hokkaido University. https://github.com/ccneko/CLAIREscope

```bibtex
@software{chung2026clairescope,
  author = {Chung, Claire},
  title = {CLAIREscope: Cellular Landscape Analysis, Interpretation \& Results Explorer for Single-Cell \& Spatial Transcriptomics},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/ccneko/CLAIREscope}},
  institution = {Department of Dermatology, Hokkaido University}
}
```

## 📜 License
CLAIREscope is open-source software released under the [MIT License](LICENSE).
