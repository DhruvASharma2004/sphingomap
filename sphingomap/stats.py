"""Differential expression and association statistics."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def bh_fdr(pvalues) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaNs preserved)."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    mask = ~np.isnan(p)
    pv = p[mask]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[mask] = adj
    return out


def differential_expression(expr: pd.DataFrame, groups: pd.Series,
                            case: str = "Tumor", control: str = "Normal") -> pd.DataFrame:
    """Unpaired tumor-vs-normal comparison on log2-scale expression.

    Uses Welch's t-test (robust to unequal variances / group sizes) plus a
    Mann-Whitney U test as a non-parametric check. log2FC is the difference
    of group means on the log2 scale.
    """
    a = expr.loc[:, groups[groups == case].index]
    b = expr.loc[:, groups[groups == control].index]
    rows = []
    for gene in expr.index:
        x, y = a.loc[gene].dropna().values, b.loc[gene].dropna().values
        if len(x) < 3 or len(y) < 3 or (np.var(x) == 0 and np.var(y) == 0):
            rows.append((gene, np.nan, np.nan, np.nan, np.nan, np.nan))
            continue
        _t, p = stats.ttest_ind(x, y, equal_var=False)
        _, pmw = stats.mannwhitneyu(x, y, alternative="two-sided")
        rows.append((gene, x.mean() - y.mean(), x.mean(), y.mean(), p, pmw))
    res = pd.DataFrame(rows, columns=["gene", "log2FC", f"mean_{case}", f"mean_{control}",
                                      "p_welch", "p_mannwhitney"]).set_index("gene")
    res["fdr"] = bh_fdr(res["p_welch"])
    res["fdr_mannwhitney"] = bh_fdr(res["p_mannwhitney"])
    return res.sort_values("p_welch")


def paired_differential_expression(expr: pd.DataFrame, groups: pd.Series,
                                   patients: pd.Series) -> pd.DataFrame:
    """Matched tumor/normal pairs from the same patient (paired t + Wilcoxon)."""
    df = pd.DataFrame({"group": groups, "patient": patients})
    tum = df[df.group == "Tumor"].reset_index(names="sample").drop_duplicates("patient").set_index("patient")["sample"]
    nor = df[df.group == "Normal"].reset_index(names="sample").drop_duplicates("patient").set_index("patient")["sample"]
    shared = tum.index.intersection(nor.index)
    if len(shared) < 3:
        return pd.DataFrame()
    t_samples, n_samples = tum.loc[shared].tolist(), nor.loc[shared].tolist()
    rows = []
    for gene in expr.index:
        x = expr.loc[gene, t_samples].values.astype(float)
        y = expr.loc[gene, n_samples].values.astype(float)
        d = x - y
        if np.allclose(d, 0) or np.isnan(d).all():
            rows.append((gene, np.nan, np.nan, np.nan))
            continue
        _, p = stats.ttest_rel(x, y, nan_policy="omit")
        _, pw = stats.wilcoxon(d[~np.isnan(d)])
        rows.append((gene, np.nanmean(d), p, pw))
    res = pd.DataFrame(rows, columns=["gene", "log2FC_paired", "p_paired_t", "p_wilcoxon"]).set_index("gene")
    res["fdr_paired"] = bh_fdr(res["p_paired_t"])
    res.attrs["n_pairs"] = len(shared)
    return res.sort_values("p_paired_t")


def stage_association(scores: pd.DataFrame, stage: pd.Series) -> pd.DataFrame:
    """Kruskal-Wallis across stages I-IV and Spearman trend for each score."""
    order = {"I": 1, "II": 2, "III": 3, "IV": 4}
    stage = stage.dropna()
    stage = stage[stage.isin(order)]
    common = scores.index.intersection(stage.index)
    rows = []
    for col in scores.columns:
        s = scores.loc[common, col]
        st = stage.loc[common]
        groups = [s[st == k].dropna().values for k in order if (st == k).sum() >= 3]
        p_kw = stats.kruskal(*groups).pvalue if len(groups) >= 2 else np.nan
        rho, p_rho = stats.spearmanr(st.map(order), s, nan_policy="omit")
        rows.append((col, len(common), p_kw, rho, p_rho))
    res = pd.DataFrame(rows, columns=["score", "n", "p_kruskal", "spearman_rho", "p_spearman"]).set_index("score")
    res["fdr_kruskal"] = bh_fdr(res["p_kruskal"])
    return res
