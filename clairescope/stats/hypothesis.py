"""Statistical hypothesis testing and significance annotation formatting."""
from typing import Tuple, Any
from scipy.stats import mannwhitneyu
import numpy as np

def get_sig_label(p: float) -> str:
    """Format p-value into publication significance asterisks."""
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"

def format_sig_value(val: Any) -> str:
    """Format statistical values cleanly for tooltips and tables."""
    try:
        v = float(val)
        if abs(v) >= 1000 or (abs(v) < 0.001 and v != 0):
            return f"{v:.2e}"
        else:
            return f"{float(f'{v:.4g}'):g}"
    except Exception:
        return str(val)

def run_mann_whitney(group1_vals: np.ndarray, group2_vals: np.ndarray) -> Tuple[float, float, str]:
    """Run two-sided Mann-Whitney U test between two numeric groups."""
    g1 = np.asarray(group1_vals).flatten()
    g2 = np.asarray(group2_vals).flatten()
    g1 = g1[~np.isnan(g1)]
    g2 = g2[~np.isnan(g2)]
    
    if len(g1) < 2 or len(g2) < 2:
        return np.nan, np.nan, "ns"
    try:
        stat, pval = mannwhitneyu(g1, g2, alternative='two-sided')
        return float(stat), float(pval), get_sig_label(pval)
    except Exception:
        return np.nan, np.nan, "ns"
