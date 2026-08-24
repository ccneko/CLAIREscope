# 🚀 Quickstart Guide

## 1. System Requirements
- **Operating System**: Linux (Ubuntu/Debian), macOS, or Windows (WSL2 recommended).
- **Python**: Python 3.10 to 3.13.
- **Memory**: Minimum 8 GB RAM (16 GB+ recommended for datasets > 50,000 cells).

## 2. Launching the Application
Run the launcher command:
```bash
streamlit run app.py --server.port 8501
```
Navigate to `http://localhost:8501` in your browser.

## 3. Dataset Configuration (`dataset_config.yaml`)
CLAIREscope automatically reads project configurations from `dataset_config.yaml`:
```yaml
projects:
  D001:
    name: "D001: JEB snRNA-seq (Natsuga Lab)"
    desc: "Junctional Epidermolysis Bullosa (JEB) snRNA-seq & Revertant Mosaicism Analysis"
    data_dir: "data/D001"
    default_dataset: "adata_kc_norm_cell_typed.h5ad"
    sample_col: "sample"
    cell_state_col: "cell_type"
```
