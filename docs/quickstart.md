# 🚀 Quickstart Guide

This guide will walk you through setting up and launching **CLAIREscope** on your local machine or high-performance compute server.

---

## 💻 Hardware & System Requirements

CLAIREscope is designed to be lightweight and resource-efficient. Memory footprint scales with the number of cells in your active AnnData (`.h5ad`) object:

| Cohort Scale | Number of Cells | Minimum RAM | Recommended RAM | Recommended CPU |
| :--- | :--- | :--- | :--- | :--- |
| **Small / Pilot** | $< 10,000$ cells | 4 – 8 GB | 8 GB | 4 Cores |
| **Standard Atlas** | $10,000 - 50,000$ cells | 8 – 16 GB | 16 GB | 8 Cores |
| **Large-Scale Atlas** | $50,000 - 150,000+$ cells | 16 – 32 GB | 32 – 64 GB | 8–16 Cores |
| **Million-Cell Atlas** | $> 200,000$ cells | 64 GB | 128 GB+ | 16–32 Cores |

* **Storage**: 2 GB free disk space for dependencies; NVMe SSD recommended for rapid `.h5ad` file I/O.
* **Operating Systems**: Linux (Ubuntu 20.04/22.04/24.04), macOS (Intel / Apple Silicon M1/M2/M3), Windows 10/11 (WSL2 or Native).
* **Python**: Version $\ge 3.10$ ($\le 3.13$).

---

## 📦 Installation via `uv` (Recommended)

[`uv`](https://github.com/astral-sh/uv) is an extremely fast Python package and environment manager.

```bash
# 1. Clone repository
git clone https://github.com/ccneko/CLAIREscope.git
cd CLAIREscope

# 2. Create isolated virtual environment
uv venv --python 3.12
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# 3. Install core single-cell dependencies
uv pip install scanpy anndata streamlit pandas numpy scipy matplotlib plotly pyyaml openpyxl
```

---

## 🔬 Launching the Platform

Run the unified launcher:
```bash
streamlit run app.py
```

CLAIREscope will initialize and open your default browser at **`http://localhost:8501`**.

---

## 🧪 Exploring with Demo Data
1. Navigate to the sidebar **🧭 Navigation** $ightarrow$ **`Single-Cell Preprocessing & Scanpy Pipeline`**.
2. Click **`💾 Download Demo AnnData (.h5ad)`** or click **`🚀 Load Synthetic Template into Pipeline`**.
3. Explore all 11 analysis studios immediately!
