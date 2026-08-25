# 🔬 CLAIREscope: Cellular Landscape Analysis, Interpretation & Results Explorer

[![Release](https://img.shields.io/badge/Release-v1.0.0-crimson.svg)](https://github.com/ccneko/CLAIREscope/releases)
[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/Docs-ReadTheDocs-brightgreen.svg)](https://clairescope.readthedocs.io)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

An open-source, lightweight, and high-performance single-cell analysis platform built on **Python**, **Scanpy**, and **Streamlit**. CLAIREscope unifies end-to-end data preprocessing, live multi-project exploration, dynamic continuous trajectory kinetics modeling, on-the-fly differential expression (Volcano), hypergeometric pathway enrichment (ORA), and automated publication-grade packaging (300 DPI vector SVGs/PDFs and structured summary CSV matrices).

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

## 🌟 Key Functional Modules

1. **🧪 Preprocessing & Scanpy Pipeline Studio**: Diagnostic data checker, green readiness checklist badges, fuzzy column schema standardizer (`sample`, `cell_type`), and 11-stage automated Scanpy engine with zero-restart in-memory viewer launch.
2. **🔬 Static & Interactive UMAP Studios**: 1:1 isometric aspect ratio preservation, Loupe-like dynamic percentile contrast clipping, and dual-scatter synchronized hover inspection.
3. **📊 Sample Composition & Stratification Studios**: Reorderable cell-state donuts, stacked percentage bars, and multi-condition gene expression violins with automated Wilcoxon statistical hypothesis testing.
4. **📈 Dynamic Trajectory Kinetics Studio**: Continuous spline kinetics overlay across Diffusion Pseudotime ($DPT$), within-sample normalization, and landmark peak alignment.
5. **🌋 Differential Expression Volcano Studio**: Full-genome bidirectional testing with live gene label overlays and publication scientific notation formatting.
6. **🎨 Hierarchical Expression Heatmap Studio**: Clustered dendrograms with draggable sortable X-axis ordering and Z-score standardization.
7. **🧬 Pathway Over-Representation Analysis (ORA)**: Real-time Hypergeometric gene set enrichment for top $X$ up- and down-regulated genes.
8. **📦 Bulk Package Export & Config Importer**: Automated generation of 300 DPI vector SVGs, PNGs, PDFs, and structured 1-row-per-feature CSV summary matrices with drag-and-drop CSV/Excel configuration import.

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
```bash
streamlit run app.py
```
Open your browser and navigate to `http://localhost:8501`.

---

## 📖 Documentation & Tutorials
For complete guides, configuration file specifications, and remote deployment protocols (NordVPN Meshnet, Cloudflare Tunnels), visit our [ReadTheDocs Manual](https://clairescope.readthedocs.io).

## 📄 Citation
If you use CLAIREscope in your research, please cite our manuscript:
> **Chung C., et al.** (2026). *CLAIREscope: An Interactive, Lightweight Single-Cell Analysis Platform for Multi-Project Exploration, Dynamic Trajectory Kinetics, and Automated Publication Reporting.* (Under Review).

## 📜 License
CLAIREscope is open-source software released under the [MIT License](LICENSE).
