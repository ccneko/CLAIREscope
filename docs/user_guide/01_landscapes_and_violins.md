# Single-Cell Landscapes & Violins

## Static UMAP Grid Studio (Tab 1)
- **Colormap & Scale**: Toggle between $\log_2(\text{Normalized} + 1)$ and linear scales.
- **Percentile Threshold Anchors**: Adjust upper colormap percentiles ($50\%$ to $100\%$) to enhance contrast against outlier cells.
- **Layout Controls**: Choose 1–6 grid columns, auto-rows, and legend positioning (*Bottom (Full Width)*, *Right (Side)*, *On-Data Labels (Centroids)*, *Hidden*).

## Interactive 2D & 3D UMAP Studio (Tab 2)
- **Auto-Highlighting**: Subsetting clusters automatically switches view mode to *"Highlight selected (dim unselected in grey)"*.
- **Sortable Filters**: Reorder clusters and quickly select/clear categories using draggable chips.
- **3D Exploration**: Rotate, zoom, and inspect continuous manifold structures in 3D space.

## Gene Expression Violins (Tab 4)
- Compare single-cell expression distributions across conditions or cell types.
- Displays non-parametric Mann-Whitney U test statistics with significance asterisks ($p < 0.05^*, p < 0.01^{**}, p < 0.001^{***}$).
