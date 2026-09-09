# 🌋 Differential Expression & Volcano Studio

The **Differential Expression & Volcano Studio** (Tab 9) provides an interactive, publication-ready environment for pairwise and group-level statistical comparisons across cell types, patient cohorts, or experimental conditions.

---

## 🎯 Key Capabilities

- **Bidirectional Wilcoxon Rank-Sum Testing**: Rapid computation of Log2(Fold Change) and Benjamini-Hochberg FDR-adjusted $p$-values across thousands of genes.
- **Multi-Term Union (OR) Gene Search**: Real-time filtering and visualization of gene families or target marker lists (e.g. `ITGB, LAM, COL` or `ITGB.LAM.COL`).
- **Interactive SVG Volcano Studio**: High-precision interactive volcano plot with dynamic significance thresholds, pure SVG z-ordering, and smart anti-collision gene callouts.
- **Spatially Aligned DEG Tables**: Synchronized tables reflecting the volcano plot layout (Downregulated on the Left, Upregulated on the Right) with instant CSV export.

---

## 🔍 Multi-Term Union Gene Search & Filtering

Under the **Display Controls** section, use the search input to filter and highlight specific genes of interest in real time.

### Flexible Query Delimiters
The search parser supports multiple delimiter types, enabling quick copy-pasting from gene lists or literature:
- **Comma-separated**: `ITGB1, LAMB3, COL17A1`
- **Dot-separated**: `ITGB.LAM.COL`
- **Space or semicolon-separated**: `KRT14 KRT5; ACTB`

### Union Matching Logic
- The search evaluates each token independently across both **Gene Symbols** and **Ensembl / Gene IDs** (case-insensitive).
- Results represent the **union (OR)** of all matches: any gene matching at least one token is selected and displayed.
- The UI dynamically displays the active search query summary:
  ```text
  Highlighting & filtering union of matches for: 'ITGB', 'LAM', 'COL'
  ```

---

## 🌋 Interactive Volcano Plot Features

The Volcano Plot provides deep insight into global transcriptomic changes with full user interactivity:

### 1. Dynamic Cutoffs
- **Log2 Fold Change Threshold ($|\log_2\text{FC}|$ cut-off)**: Adjust via slider to filter for biologically meaningful effect sizes (default: `1.0`).
- **FDR Adjusted $p$-value Threshold ($-\log_{10}(p_{\text{adj}})$)**: Filter for statistical significance (default: `0.05`).
- Visualized as dashed reference threshold lines (red/blue for fold changes, grey for significance).

### 2. Pure SVG Z-Ordering & Highlight Rings
- Built with Plotly's explicit SVG rendering mode (`render_mode="svg"`), ensuring consistent DOM layering across all devices:
  1. **Background Layer**: Non-significant genes rendered in muted grey (`#CBD5E1`).
  2. **Significance Layer**: Upregulated genes in coral red (`#EF4444`) and downregulated genes in vibrant blue (`#3B82F6`).
  3. **Highlight Layer**: Matched genes are enclosed in crisp, hollow purple marker rings (`#8B5CF6`, 3px stroke width) strictly on the top layer, preserving each gene's true significance color inside the ring.

### 3. Anti-Collision Labels & Headroom Padding
- **Smart Dynamic Offsetting**: Gene label badges and arrows dynamically adjust their angles based on scatter density and axis bounds:
  - Points clustered at ceiling thresholds (e.g. capped at $y = 300$) are staggered upward and outward.
  - Points near the baseline ($y = 0$) are offset upward to prevent axis overlap.
- **Y-Axis Headroom**: Automatic padding on the top and bottom of the Y-axis ensures that extreme $p$-value dots and rings are never clipped by chart borders.

---

## 📋 Spatially Aligned DEG Result Tables

Beneath the Volcano Plot, the **Top Differentially Expressed Genes** section features dual mirrored data tables aligned with the volcano plot's spatial orientation:

| Column | Content | Significance Filter |
| :--- | :--- | :--- |
| **Left Table** | **Top Downregulated Genes** in target group | Negative $\log_2\text{FC}$, FDR $\le 0.05$ |
| **Right Table** | **Top Upregulated Genes** in target group | Positive $\log_2\text{FC}$, FDR $\le 0.05$ |

### Table Features
- **Real-Time Search Sync**: When a search filter is active, tables display only the matching subset (e.g. `Matches (18)`).
- **Formatted Statistics**: Displays Gene Symbol, Gene ID, Log2FC (3 decimals), FDR adjusted $p$-value in clean scientific notation (e.g. `< 1e-300`, `2.45e-12`), and Wilcoxon Z-score.
- **📥 Download Full Results**: One-click button to download the entire unthresholded differential expression result table as a CSV file (`DE_<Target>_vs_<Reference>_<GroupCol>.csv`).
