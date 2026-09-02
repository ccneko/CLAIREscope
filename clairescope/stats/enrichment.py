"""Hypergeometric over-representation analysis (ORA) engine."""
from typing import List, Dict, Any
from scipy.stats import hypergeom
import pandas as pd

def run_hypergeometric_enrichment(query_genes: List[str], background_genes: List[str], pathway_dict: Dict[str, List[str]], min_overlap: int = 2) -> pd.DataFrame:
    """Run hypergeometric ORA test for a gene set against biological pathways."""
    bg_set = set(str(g).upper() for g in background_genes)
    q_set = set(str(g).upper() for g in query_genes).intersection(bg_set)
    
    N = len(bg_set)
    n = len(q_set)
    
    if N == 0 or n == 0:
        return pd.DataFrame(columns=["Pathway", "Overlap", "Pathway_Size", "Rich_Factor", "p_value", "Overlap_Genes"])
        
    records = []
    for p_name, p_genes in pathway_dict.items():
        p_set = set(str(g).upper() for g in p_genes).intersection(bg_set)
        K = len(p_set)
        if K == 0:
            continue
            
        overlap = q_set.intersection(p_set)
        k = len(overlap)
        
        if k >= min_overlap:
            pval = hypergeom.sf(k - 1, N, K, n)
            rich_factor = (k / n) / (K / N) if (K > 0 and n > 0) else 0.0
            records.append({
                "Pathway": p_name,
                "Overlap": k,
                "Pathway_Size": K,
                "Rich_Factor": round(rich_factor, 2),
                "p_value": float(pval),
                "Overlap_Genes": ", ".join(sorted(list(overlap)))
            })
            
    df = pd.DataFrame(records)
    if not df.empty:
        df = df.sort_values("p_value", ascending=True).reset_index(drop=True)
    return df
