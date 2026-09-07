# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.3] - 2026-09-08

### Added
- **Interactive First Experience**: Swapped tab hierarchy to place **✨ Interactive UMAP** as Tab 1 for immediate exploration, followed by **🗺️ Static UMAP** as Tab 2.
- **On-Data Centroid Labels Default**: Configured publication-grade on-data cluster centroid annotations as the default mode for Static UMAP.
- **Paul Tol Color System**: Integrated the Paul Tol Qualitative Palette into semantic schema rules for accessible, high-contrast, and colorblind-safe categorical palettes.
- **Live Demo Epidermal Atlas**: Bundled `demo_wang2020_epidermis.h5ad` with baked schema colors for zero-configuration testing and public showcase.

### Changed
- **Enlarged Interactive Legend Typography**: Scaled interactive Plotly UMAP legend font sizing (13–15pt) to optimize legibility on high-resolution displays.
- **Robust Streamlit Cloud Compatibility**: Standardized dependency bounds in `requirements.txt` to eliminate container build failures on Streamlit Community Cloud.

---

## [1.0.2] - 2026-09-07

### Added
- **Bidirectional Gene & Species Resolution**: Robust case-insensitive gene identifier mapping supporting Ensembl IDs, Gene Symbols, and mouse/human orthologs across multi-species datasets (e.g. D002 Watanabe Lab, D004 Wang 2020).
- **Synchronized Multiselect Legend Ordering**: Both Static UMAP (multi-panel grids and reference plots) and Interactive Plotly UMAPs strictly reflect user custom drag-and-drop sort order in legends.
- **Dynamic View Filtering & Highlighting on Static UMAP**: Extended "Color all cells", "Highlight selected", and "Filter view" modes with interactive background dimming to static grids and reference plots.
- **Configurable Plotly Legend Layouts**: Added customizable Interactive UMAP legend orientations (`Bottom (Horizontal)`, `Right Side (Compact)`, `On-Plot Centroids`, `Hide Legend`) with label truncation to prevent plot compression.
- **Comprehensive AST & Code Integrity Tests**: Added syntax and symbol presence verification in test suite (`tests/run_tests.py`), expanding suite to 21 unit & integrity tests.

### Changed
- **Overhauled Categorical Color Mapping Engine**: Implemented semantic skin lineage color mappings with a 50+ high-contrast palette pool in `clairescope/core/schema.py`, eliminating unintended grey collapses across cell states.
- **Point Size Slider Integration**: Unified `Point Size:` controls across all static grid subplots and reference plots with adaptive background point sizing.
- **Transient State Isolation**: Added automatic widget key cleanup upon dataset/column switching, preventing stale filter exceptions.

### Fixed
- **Dataset Switch KeyError**: Resolved transient filter state collisions when toggling between datasets with differing annotation columns and categories.

---

## [1.0.1] - 2026-09-06

### Added
- **Full-Width Bottom Legend Layout**: Default legend position placed below plots with dynamic multi-column wrapping to prevent coordinate axes shrinking on long cluster names.
- **On-Data Centroid Cluster Labels**: Added Seurat/Scanpy-style cluster label overlays directly at cluster centroids with white halo text outlines (`matplotlib.patheffects`).
- **Dynamic Dataset Scanner**: Automatic recursive directory scanner (`scan_project_datasets`) detecting `.h5ad` files in configured project root paths.
- **Interactive Project Ingestion**: Sidebar `"➕ New Project..."` option with drag-and-drop file ingestion, automatic default numbering, and root folder path detection.
- **Sortable Multiselect**: Sortable draggable chips with *Select All* and *Clear* controls for cell state and annotation filtering.
- **Intelligent Auto-Highlighting**: Subsetting clusters in Interactive UMAP automatically switches view mode to *"Highlight selected (dim unselected in grey)"*.

### Changed
- **Dynamic Colormap Reactivity**: Synchronized `Colormap Max (vmax)` directly with `Max Percentile Threshold (Anchors)` slider and selected gene expression across Static UMAP, Interactive UMAP, and Heatmap views.
- **Streamlined Layout Controls**: Unified Static Grid controls into a 3-column toolbar (`Grid Columns`, `Grid Rows`, `Legend Position`).

### Fixed
- **Scanpy Scoring Fallback**: Enforced `use_raw=False` and pre-filtered variable validation in `compute_score` to prevent `ValueError: No valid genes were passed for scoring`.
- **Trajectory Gallery Paths**: Fixed figure directory path resolution to dynamically scan project figure subdirectories.

---

## [1.0.0] - 2026-09-04

### Added
- **Static Grid UMAP Studio**: Multi-panel side-by-side comparison grids with Loupe Browser-style Log2 scaling and anchor contrast thresholds.
- **Interactive 2D & 3D UMAP Studio**: Hardware-accelerated Plotly WebGL scatter viewer with 3D camera controls and dynamic point sizing.
- **Sample Composition Studio**: Proportional donut charts and stacked distributions with Fisher's exact and Chi-square statistics.
- **Gene Expression Violins**: Multi-gene distribution violins with automated Mann-Whitney U test p-values.
- **Signature & Pathway Scoring**: Built-in score generator for custom gene sets and cutaneous biology signatures.
- **Bivariate Correlation Studio**: Co-expression scatter plots with linear regression and Pearson / Spearman correlation statistics.
- **Developmental Trajectory Analysis**: Diffusion Pseudotime (DPT) curves, vector trajectories, and precomputed figure galleries.
- **Differential Expression Studio**: Interactive volcano plots with customizable fold-change and p-value cutoffs.
- **Pathway Enrichment Studio**: Hypergeometric over-representation analysis against curated Gene Ontology and pathway databases.
- **Expression Heatmap Studio**: Publication-grade hierarchical clustered heatmaps with Z-score standardization.
- **Curated Cutaneous Biology Knowledgebase**: Integrated signatures spanning epidermal differentiation, hair follicle compartments, fibroblast subpopulations, immune infiltrates, and wound healing cascades.
- **Native Desktop Integration**: 1-click launcher (`CLAIREscope.bat`) and Python Tkinter Server Manager GUI.
- **Zenodo DOI Integration**: Official persistent identifier [`10.5281/zenodo.22308479`](https://doi.org/10.5281/zenodo.22308479).
