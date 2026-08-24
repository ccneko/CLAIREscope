# 🗺️ Future Roadmap & Development Plans

## 1. SVG Diagram Coloring Engine
- **Interactive Anatomical & Pathway SVGs**: Upload custom SVG diagrams (e.g., epidermal basement membrane, desmosome junction complex, hair follicle niche).
- **Dynamic CSS/SVG Fill Coloring**: Dynamically color individual SVG path elements according to gene expression levels, condition fold-changes, or trajectory kinetic states.

## 2. Cloud Data Storage & Live Connectors
- **Amazon S3 Buckets**: Direct loading and streaming of remote AnnData objects via S3 URIs (`s3://...`) using `fsspec` / `s3fs`.
- **Google Drive Integration**: Direct cloud file picking and authentication for Google Drive-stored `.h5ad` datasets.
- **Zarr / H5Coro Streaming**: Sparse chunk streaming to allow visualization of multimillion-cell atlases with minimal local memory footprint.

## 3. Advanced Gene Regulatory Network (GRN) Modeling
- In silico perturbation simulations (TF knockdown / overexpression) overlaid on trajectory kinetics curves.
