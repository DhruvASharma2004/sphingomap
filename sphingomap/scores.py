"""Per-sample pathway module scores and the sphingolipid rheostat index."""
from __future__ import annotations

import pandas as pd

from .genesets import MODULES, RHEOSTAT


def zscore_rows(expr: pd.DataFrame, reference: list[str] | None = None) -> pd.DataFrame:
    """Gene-wise z-scores. If ``reference`` samples are given (e.g. normals),
    centre and scale on them so scores read as 'deviation from normal'."""
    ref = expr if reference is None else expr[reference]
    mu, sd = ref.mean(axis=1), ref.std(axis=1, ddof=1).replace(0, float("nan"))
    return expr.sub(mu, axis=0).div(sd, axis=0)


def module_scores(expr: pd.DataFrame, reference: list[str] | None = None,
                  min_genes: int = 2) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Mean z-score of each module's detected genes, per sample.

    Returns (samples x modules score table, genes actually used per module).
    """
    z = zscore_rows(expr, reference)
    scores, used = {}, {}
    for name, mod in MODULES.items():
        genes = [g for g in mod["genes"] if g in z.index and z.loc[g].notna().any()]
        used[name] = genes
        if len(genes) >= min_genes:
            scores[name] = z.loc[genes].mean(axis=0)
    return pd.DataFrame(scores), used


def rheostat_index(expr: pd.DataFrame, reference: list[str] | None = None) -> pd.Series:
    """Pro-survival (S1P/GSL) minus pro-ceramide mean z-score, per sample.

    Higher values = transcriptional profile tilted away from ceramide
    accumulation toward S1P/glycosphingolipid production. Exploratory only.
    """
    z = zscore_rows(expr, reference)
    up = [g for g in RHEOSTAT["pro_survival"] if g in z.index]
    down = [g for g in RHEOSTAT["pro_ceramide"] if g in z.index]
    if not up or not down:
        raise ValueError("Not enough rheostat genes present in the expression matrix")
    return (z.loc[up].mean(axis=0) - z.loc[down].mean(axis=0)).rename("rheostat_index")
