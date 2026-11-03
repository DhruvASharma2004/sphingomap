"""Publication-style figures."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from .genesets import MODULES, PATHWAY_LAYOUT, PATHWAY_NODES, gene_to_module

MODULE_COLORS = {
    "de_novo_synthesis": "#4C72B0", "sphingomyelin_cycle": "#55A868",
    "ceramide_catabolism": "#C44E52", "s1p_axis": "#8172B2",
    "glycosphingolipid_branch": "#CCB974", "transport": "#64B5CD",
}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 110, "savefig.bbox": "tight"})


def volcano(de: pd.DataFrame, path, fdr: float = 0.05, lfc: float = 0.5, title: str = ""):
    mod = gene_to_module()
    d = de.dropna(subset=["log2FC", "fdr"]).copy()
    d["nlp"] = -np.log10(d["fdr"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for name, color in MODULE_COLORS.items():
        sub = d[d.index.map(mod.get) == name]
        sig = (sub.fdr < fdr) & (sub.log2FC.abs() >= lfc)
        ax.scatter(sub.log2FC[~sig], sub.nlp[~sig], s=30, color=color, alpha=0.35, edgecolor="none")
        ax.scatter(sub.log2FC[sig], sub.nlp[sig], s=45, color=color, edgecolor="k", lw=0.5,
                   label=MODULES[name]["label"])
    for g, r in d[(d.fdr < fdr) & (d.log2FC.abs() >= lfc)].nlargest(14, "nlp").iterrows():
        ax.annotate(g, (r.log2FC, r.nlp), xytext=(4, 2), textcoords="offset points", fontsize=8)
    ax.axhline(-np.log10(fdr), ls="--", c="grey", lw=0.8)
    for x in (-lfc, lfc):
        ax.axvline(x, ls="--", c="grey", lw=0.8)
    ax.set(xlabel="log2 fold change (tumor vs normal)", ylabel="-log10 FDR", title=title)
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.savefig(path)
    plt.close(fig)


def module_heatmap(expr: pd.DataFrame, groups: pd.Series, path, title: str = ""):
    """Row-z-scored heatmap of pathway genes, samples ordered Normal | Tumor."""
    mod = gene_to_module()
    genes = [g for name in MODULES for g in MODULES[name]["genes"] if g in expr.index]
    order = list(groups[groups == "Normal"].index) + list(groups[groups == "Tumor"].index)
    z = expr.loc[genes, order]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0).clip(-3, 3)
    fig, (axc, ax) = plt.subplots(1, 2, figsize=(10, 0.18 * len(genes) + 1.5),
                                  gridspec_kw={"width_ratios": [0.15, 10], "wspace": 0.02})
    axc.imshow(np.array([[list(MODULE_COLORS).index(mod[g])] for g in genes]),
               aspect="auto", cmap=matplotlib.colors.ListedColormap(list(MODULE_COLORS.values())),
               vmin=0, vmax=len(MODULE_COLORS) - 1)
    axc.set_xticks([])
    axc.set_yticks(range(len(genes)), genes, fontsize=7)
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3, interpolation="nearest")
    n_norm = (groups == "Normal").sum()
    ax.axvline(n_norm - 0.5, c="k", lw=1)
    ax.set_yticks([])
    ax.set_xticks([n_norm / 2, n_norm + (len(order) - n_norm) / 2],
                  [f"Normal (n={n_norm})", f"Tumor (n={len(order) - n_norm})"])
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="row z-score")
    fig.savefig(path)
    plt.close(fig)


def pathway_map(de: pd.DataFrame, path, title: str = "Sphingolipid pathway: tumor vs normal"):
    """Schematic pathway with enzymes coloured by log2FC (FDR<0.05 outlined)."""
    fig, ax = plt.subplots(figsize=(10, 6.5))
    vmax = max(1.0, float(de["log2FC"].abs().reindex(PATHWAY_LAYOUT).max(skipna=True) or 1))
    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r")
    for label, (x, y) in PATHWAY_NODES.items():
        ax.text(x, y, label, ha="center", va="center", fontsize=10, weight="bold",
                bbox={"boxstyle": "round,pad=0.4", "fc": "#f3f3f3", "ec": "#666"})
    arrows = [((0.5, 6.9), (4.3, 4.8)), ((4.2, 4.3), (1.6, 3.2)), ((1.6, 2.8), (4.2, 4.2)),
              ((4.9, 4.4), (5.6, 4.05)), ((6.8, 4.0), (7.9, 4.0)), ((8.2, 3.7), (8.2, 1.85)),
              ((4.7, 4.9), (6.2, 8.1))]
    for (x0, y0), (x1, y1) in arrows:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops={"arrowstyle": "->", "color": "#999", "lw": 1.4, "connectionstyle": "arc3,rad=0.15"})
    for gene, (x, y) in PATHWAY_LAYOUT.items():
        if gene in de.index and pd.notna(de.loc[gene, "log2FC"]):
            lfc, q = de.loc[gene, "log2FC"], de.loc[gene, "fdr"]
            fc, ec, lw = cmap(norm(lfc)), ("k" if q < 0.05 else "#bbb"), (1.5 if q < 0.05 else 0.6)
            tc = "white" if abs(norm(lfc) - 0.5) > 0.33 else "black"
        else:
            fc, ec, lw, tc = "white", "#ddd", 0.6, "#999"
        ax.text(x, y - 0.35 if y > 7 else y, gene, ha="center", va="center", fontsize=8, color=tc,
                bbox={"boxstyle": "round,pad=0.25", "fc": fc, "ec": ec, "lw": lw})
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, label="log2FC (black outline: FDR < 0.05)")
    ax.set(xlim=(-0.8, 9.8), ylim=(0.6, 9), title=title)
    ax.axis("off")
    fig.savefig(path)
    plt.close(fig)


def score_boxplot(scores: pd.DataFrame, groups: pd.Series, path, title: str = ""):
    cols = list(scores.columns)
    fig, ax = plt.subplots(figsize=(1.3 * len(cols) + 1.5, 4.2))
    pos = np.arange(len(cols))
    for i, (grp, color) in enumerate((("Normal", "#7fa7d1"), ("Tumor", "#d6705f"))):
        data = [scores.loc[groups[groups == grp].index.intersection(scores.index), c].dropna() for c in cols]
        bp = ax.boxplot(data, positions=pos + (i - 0.5) * 0.35, widths=0.3, patch_artist=True,
                        showfliers=False, medianprops={"color": "k"})
        for b in bp["boxes"]:
            b.set(facecolor=color, alpha=0.8)
        ax.plot([], [], color=color, lw=8, label=grp)
    labels = [MODULES[c]["label"] if c in MODULES else c.replace("_", " ") for c in cols]
    ax.set_xticks(pos, labels, rotation=30, ha="right")
    ax.axhline(0, c="grey", lw=0.6, ls=":")
    ax.set(ylabel="score (z vs. normal)", title=title)
    ax.legend(frameon=False)
    fig.savefig(path)
    plt.close(fig)


def km_plot(fits: dict, stats: dict, path, title: str = ""):
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for kmf, color in zip(fits.values(), ("#d6705f", "#4C72B0"), strict=False):
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color)
    ax.set(xlabel="Days", ylabel="Overall survival probability", ylim=(0, 1.02), title=title)
    ax.text(0.03, 0.05, f"log-rank p = {stats['logrank_p']:.3g}\n"
            f"Cox HR = {stats['cox_hr_per_unit']:.2f} "
            f"({stats['cox_hr_ci_low']:.2f}-{stats['cox_hr_ci_high']:.2f})",
            transform=ax.transAxes, fontsize=9)
    fig.savefig(path)
    plt.close(fig)
