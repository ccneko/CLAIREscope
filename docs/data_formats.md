# 📂 Supported Data Formats

## 1. AnnData (`.h5ad`)
CLAIREscope natively reads standard AnnData `.h5ad` files:
- **`adata.X`**: Log-normalized expression matrix ($\log 1p$ normalized or CP10k).
- **`adata.obs`**: Cell-level metadata (must include sample columns and cell-type annotations).
- **`adata.var`**: Gene metadata with symbols and Ensembl IDs.
- **`adata.obsm['X_umap']`**: 2D UMAP coordinates for static and interactive scatter plots.
- **`adata.obs['dpt_pseudotime']`**: Optional trajectory diffusion pseudotime values.

## 2. Configuration Importer Formats
CLAIREscope supports importing gene lists and custom pathway signatures via `.csv` or `.xlsx`:
- **Gene List CSV**: A single column file containing gene names (e.g. `Gene`, `COL17A1`, `KRT14`).
- **Pathway CSV**: Multiple columns where the header row is the **Pathway Name** and rows below list member genes.
- **Multi-Sheet Excel**: A `.xlsx` workbook containing a `Genes` sheet and a `Pathways` sheet.
