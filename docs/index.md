# 🔬 CLAIREscope: Cellular Landscape Analysis, Interpretation & Results Explorer

[![Release](https://img.shields.io/github/v/release/ccneko/sc-expression-viewer?color=B32141&label=Release)](https://github.com/ccneko/sc-expression-viewer/releases)
[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**CLAIREscope** is a unified, publication-grade interactive single-cell analysis platform designed for exploratory transcriptomics, multi-project switching, developmental trajectory kinetics modeling, on-the-fly functional genomics, and automated multi-panel export.

---

## 🌟 Key Capabilities

1. **🚀 Dynamic Multi-Project Switching**: Instantly switch between heterogeneous single-cell cohorts (e.g. JEB snRNA-seq, Wound Healing scMultiome) with on-demand memory management.
2. **🗺️ Publication-Grade Isometric Visualizations**: 1:1 aspect ratio UMAPs with Loupe-like dynamic range contrast anchors, customizable colormaps, and one-click 300 DPI SVG/PNG/PDF downloads.
3. **🌿 Continuous Trajectory Expression Drawing Studio**: Multi-gene overlay kinetics along Diffusion Pseudotime (DPT), within-sample scaling, and landmark alignment.
4. **🌋 Differential Expression & Volcano Studio**: Bidirectional Wilcoxon statistical testing with tunable FDR/$\log_2	ext{FC}$ thresholds and live gene labeling.
5. **🔥 Hierarchical Heatmap Studio**: Group-averaged or single-cell clustered heatmaps with draggable sortable X-axis ordering and Z-score standardization.
6. **🧬 Pathway Over-Representation Analysis (ORA)**: Tunable top $X$ up/down gene enrichment across Hallmark and custom signature databases.
7. **📦 Bulk Export & Config Import Studio**: Generate comprehensive multi-figure packages and structured 1-row-per-feature CSV matrices with drag-and-drop Excel/CSV configuration importing.

---

## 🛠️ Quick Installation & Launch

```bash
# Clone the repository
git clone https://github.com/ccneko/sc-expression-viewer.git
cd sc-expression-viewer

# Run using uv / Python virtual environment
uv run streamlit run app.py --server.port 8501
```
