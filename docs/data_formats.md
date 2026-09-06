# 📂 Supported Data Formats & Ingestion Guide

<span class="brand">CLAIREscope</span> natively integrates with standard single-cell genomics data structures, tabular matrices, and configuration formats.

---

## 1. AnnData (`.h5ad`) Data Structure

<span class="brand">CLAIREscope</span> reads standard **AnnData** (Annotated Data) `.h5ad` objects. Memory scaling is optimized through sparse CSR matrix representations and on-demand caching.

### Expected AnnData Slots

| Slot | Required | Description | Example Format |
| :--- | :---: | :--- | :--- |
| **`adata.X`** | **Yes** | Log-normalized or raw count expression matrix (cells × genes). Sparse `csr_matrix` or dense `numpy.ndarray`. | `log1p` normalized (CP10k / 10,000 counts per cell) |
| **`adata.obs`** | **Yes** | Cell-level observational metadata table. | `sample`, `cell_type`, `condition`, `donor`, `n_counts` |
| **`adata.var`** | **Yes** | Gene-level metadata table containing gene symbols and identifiers. | `gene_symbols`, `gene_ids` (Ensembl format: `ENSG...`) |
| **`adata.obsm['X_umap']`** | **Yes** | 2D or 3D coordinate embeddings for spatial/UMAP projections. | `numpy.ndarray` with shape `(n_cells, 2)` or `(n_cells, 3)` |
| **`adata.uns`** | Optional | Unstructured annotations, cluster palette hex codes, and dataset metadata. | `cell_type_colors: ['#1f77b4', '#ff7f0e', ...]` |
| **`adata.obs['dpt_pseudotime']`** | Optional | Precomputed Diffusion Pseudotime values (range: `0.0` to `1.0`). | Continuous float values for trajectory ordering |

---

## 2. Converting Seurat (R) Objects to AnnData (`.h5ad`)

If your analysis was performed in R with Seurat, you can export to `.h5ad` using **`zellkonverter`** (Bioconductor) or **`SeuratDisk`**:

### Using R (`zellkonverter` — Recommended)
```r
# Install zellkonverter if needed
BiocManager::install("zellkonverter")

library(Seurat)
library(zellkonverter)

# Convert Seurat object to SingleCellExperiment, then write .h5ad
sce <- as.SingleCellExperiment(seurat_obj)
writeH5AD(sce, file = "my_dataset.h5ad", X_name = "logcounts")
```

### Using Python (`scanpy` from 10x Genomics Outputs)
```python
import scanpy as sc

# Load standard 10x Genomics cellranger output folder
adata = sc.read_10x_mtx("filtered_feature_bc_matrix/", var_names="gene_symbols", cache=True)

# Standard preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="cell_type")

# Save processed AnnData
adata.write("my_dataset.h5ad")
```

---

## 3. Custom Gene Lists & Pathway Signatures (CSV / Excel)

<span class="brand">CLAIREscope</span> supports on-the-fly importing of gene sets and pathway signatures across multiple file formats in the **Signature & Pathway Scoring Studio**:

### A. Single-Column Gene List (`.csv` / `.txt`)
A simple one-column text or CSV file listing target gene symbols:

```text
Gene
COL17A1
KRT14
KRT5
ITGA6
ITGB4
```

### B. Multi-Column Pathway Signatures (`.csv`)
A CSV table where each column header represents a **Pathway / Signature Name**, and rows contain member genes:

| Basal_Stem | Spinous_Differentiation | Keratinization |
| :--- | :--- | :--- |
| COL17A1 | KRT1 | LOR |
| KRT14 | KRT10 | FLG |
| KRT5 | SBSN | IVL |
| ITGA6 | KRT16 | DSG1 |
| ITGB4 | CDH1 | TGM1 |

### C. Multi-Sheet Excel Workbook (`.xlsx`)
A `.xlsx` file containing designated sheets:
- **`Genes` Sheet**: Custom gene lists for expression violins and heatmaps.
- **`Pathways` Sheet**: Multi-column signature definitions for ssGSEA and Scanpy gene scoring.

---

## 4. Configuration YAML Files

<span class="brand">CLAIREscope</span> uses structured YAML configuration files under `config/` to manage active projects, cell-type marker databases, and signature catalogs:

### Project Configuration (`config/user/projects.yaml`)
```yaml
projects:
  Project_A_Lab A_Condition B:
    name: "Project_A: Condition B Somatic Condition A snRNA-seq (Lab A)"
    root: "data/Project_A"
    desc: "Single-nucleus RNA-seq of Junctional Epidermolysis Bullosa (Control, A1, A2, B)"
    sample_col: "sample"
    annotation_col: "cell_state_annotated"
```

### Marker Genes Configuration (`config/defaults/markers.yaml`)
```yaml
markers:
  Keratinocytes:
    Basal_1: ["COL17A1", "KRT14", "KRT5", "ITGA6"]
    Basal_2: ["POSTN", "MMP3", "PTGS2", "IL6"]
    Spinous: ["KRT1", "KRT10", "SBSN"]
    Granular: ["LOR", "FLG", "LCE1A", "IVL"]
```
