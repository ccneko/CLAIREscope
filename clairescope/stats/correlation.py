"""Bivariate co-expression and correlation analytical engines."""
from typing import Tuple, Dict, Any, Optional
import numpy as np
from scipy.stats import pearsonr, spearmanr

def compute_bivariate_correlation(x_vals: np.ndarray, y_vals: np.ndarray, filter_mode: str = "all") -> Dict[str, Any]:
    """Calculate Pearson and Spearman correlation coefficients with filtering."""
    x = np.asarray(x_vals).flatten()
    y = np.asarray(y_vals).flatten()
    
    valid_mask = ~(np.isnan(x) | np.isnan(y))
    if filter_mode == "x_pos":
        valid_mask &= (x > 0)
    elif filter_mode == "y_pos":
        valid_mask &= (y > 0)
    elif filter_mode == "co_expressed":
        valid_mask &= (x > 0) & (y > 0)
        
    x_filt = x[valid_mask]
    y_filt = y[valid_mask]
    
    n_cells = len(x_filt)
    if n_cells < 3 or np.all(x_filt == x_filt[0]) or np.all(y_filt == y_filt[0]):
        return {
            "n_cells": n_cells,
            "pearson_r": np.nan, "pearson_p": np.nan,
            "spearman_rho": np.nan, "spearman_p": np.nan,
            "slope": np.nan, "intercept": np.nan,
            "x": x_filt, "y": y_filt
        }
        
    pr_r, pr_p = pearsonr(x_filt, y_filt)
    sp_rho, sp_p = spearmanr(x_filt, y_filt)
    
    # Linear fit
    poly = np.polyfit(x_filt, y_filt, 1)
    
    return {
        "n_cells": n_cells,
        "pearson_r": float(pr_r), "pearson_p": float(pr_p),
        "spearman_rho": float(sp_rho), "spearman_p": float(sp_p),
        "slope": float(poly[0]), "intercept": float(poly[1]),
        "x": x_filt, "y": y_filt
    }
