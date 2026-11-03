"""End-to-end analysis: expression -> DE -> module scores -> clinical links -> report."""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pandas as pd

from . import __version__, plots
from .genesets import MODULES, all_genes
from .io import load_clinical, load_expression, patient_id, sample_type
from .scores import module_scores, rheostat_index
from .stats import differential_expression, paired_differential_expression, stage_association
from .survival import km_curves, median_split_logrank, prepare


def run(expression: list[str | Path], outdir: str | Path, clinical: list[str | Path] | None = None,
        survival: list[str | Path] | None = None, label: str = "cohort") -> dict:
    out = Path(outdir)
    figs = out / "figures"
    figs.mkdir(parents=True, exist_ok=True)

    # --- load & merge cohorts ---------------------------------------------
    mats = [load_expression(p) for p in expression]
    expr = pd.concat(mats, axis=1, join="inner")
    expr = expr.loc[:, ~expr.columns.duplicated()]
    clinical, survival = clinical or [], survival or []
    clin_frames = [load_clinical(c, s) for c, s in zip(
        clinical + [None] * (len(expression) - len(clinical)),
        survival + [None] * (len(expression) - len(survival)), strict=True)]
    clin = pd.concat([c for c in clin_frames if not c.empty]) if any(not c.empty for c in clin_frames) \
        else pd.DataFrame(index=expr.columns)
    clin = clin[~clin.index.duplicated()]

    groups = pd.Series({s: sample_type(s) for s in expr.columns})
    expr = expr.loc[:, groups.isin(["Tumor", "Normal"])]
    groups = groups[expr.columns]
    patients = pd.Series({s: patient_id(s) for s in expr.columns})

    genes = [g for g in all_genes() if g in expr.index]
    missing = [g for g in all_genes() if g not in expr.index]
    sp = expr.loc[genes]
    # drop near-silent genes (mean log2 expression < 1 in both groups)
    expressed = sp[(sp.T.groupby(groups).mean().T > 1).any(axis=1)]

    summary = {"label": label, "version": __version__, "date": str(date.today()),
               "n_tumor": int((groups == "Tumor").sum()), "n_normal": int((groups == "Normal").sum()),
               "n_pathway_genes_found": len(genes), "missing_genes": missing,
               "n_expressed": len(expressed)}

    # --- differential expression -------------------------------------------
    de = differential_expression(expressed, groups)
    de.insert(0, "module", de.index.map({g: m for m, d in MODULES.items() for g in d["genes"]}))
    paired = paired_differential_expression(expressed, groups, patients)
    if not paired.empty:
        de = de.join(paired)
        summary["n_matched_pairs"] = paired.attrs["n_pairs"]
    de.to_csv(out / "differential_expression.csv")
    sig = de[(de.fdr < 0.05) & (de.log2FC.abs() >= 0.5)]
    summary["n_significant"] = len(sig)
    summary["top_up"] = sig[sig.log2FC > 0].sort_values("log2FC", ascending=False).index[:8].tolist()
    summary["top_down"] = sig[sig.log2FC < 0].sort_values("log2FC").index[:8].tolist()

    # --- scores -----------------------------------------------------------
    normals = list(groups[groups == "Normal"].index)
    scores, _used = module_scores(expressed, reference=normals)
    scores["rheostat_index"] = rheostat_index(expressed, reference=normals)
    scores.insert(0, "group", groups)
    scores.to_csv(out / "sample_scores.csv")
    num = scores.drop(columns="group")
    score_de = differential_expression(num.T, groups)[["log2FC", "p_welch", "fdr"]]
    score_de = score_de.rename(columns={"log2FC": "mean_difference"})
    score_de.to_csv(out / "score_tumor_vs_normal.csv")
    summary["score_tests"] = score_de.round(6).reset_index().to_dict(orient="records")

    # --- figures ------------------------------------------------------------
    plots.volcano(de, figs / "volcano.png", title=f"{label}: sphingolipid genes, tumor vs normal")
    plots.pathway_map(de, figs / "pathway_map.png", title=f"{label}: sphingolipid pathway, tumor vs normal")
    plots.module_heatmap(expressed, groups, figs / "heatmap.png", title=f"{label}: pathway gene expression")
    plots.score_boxplot(num, groups, figs / "module_scores.png", title=f"{label}: module scores")

    # --- clinical associations (tumors only) -----------------------------------
    tumors = groups[groups == "Tumor"].index
    if "stage_group" in clin.columns:
        st = stage_association(num.loc[tumors], clin["stage_group"])
        st.to_csv(out / "stage_association.csv")
        summary["stage_tests"] = st.round(6).reset_index().to_dict(orient="records")
    if {"OS", "OS.time"} <= set(clin.columns):
        df = prepare(num.loc[tumors, "rheostat_index"], clin)
        if len(df) >= 20 and df.event.sum() >= 5:
            km = median_split_logrank(df)
            plots.km_plot(km_curves(df), km, figs / "km_rheostat.png",
                          title=f"{label}: OS by rheostat index (median split)")
            summary["survival_rheostat"] = km

    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    write_report(summary, de, out)
    return summary


def write_report(summary: dict, de: pd.DataFrame, out: Path) -> None:
    lines = [f"# Sphingolipid pathway report — {summary['label']}", "",
             f"_sphingomap v{summary['version']}, generated {summary['date']}_", "",
             "## Cohort", "",
             f"- Tumor samples: **{summary['n_tumor']}**  ",
             f"- Normal samples: **{summary['n_normal']}**"
             + (f" ({summary['n_matched_pairs']} matched pairs)" if "n_matched_pairs" in summary else ""),
             f"- Pathway genes found: {summary['n_pathway_genes_found']} "
             f"({summary['n_expressed']} expressed)"
             + (f"; not in matrix: {', '.join(summary['missing_genes'])}" if summary['missing_genes'] else ""),
             "", "## Differential expression (Welch t-test, BH FDR)", "",
             f"{summary['n_significant']} genes with FDR < 0.05 and |log2FC| ≥ 0.5.", "",
             f"- Up in tumor: {', '.join(summary['top_up']) or '—'}",
             f"- Down in tumor: {', '.join(summary['top_down']) or '—'}", "",
             "![pathway](figures/pathway_map.png)", "", "![volcano](figures/volcano.png)", "",
             "| Gene | Module | log2FC | FDR |", "|---|---|---:|---:|"]
    for g, r in de.dropna(subset=["fdr"]).head(20).iterrows():
        lines.append(f"| {g} | {MODULES.get(r.module, {}).get('label', '')} | {r.log2FC:+.2f} | {r.fdr:.2e} |")
    lines += ["", "## Module scores and rheostat index", "",
              "Scores are mean z-scores relative to normal tissue.", "",
              "| Score | Tumor − Normal | FDR |", "|---|---:|---:|"]
    for t in summary["score_tests"]:
        lines.append(f"| {t['gene']} | {t['mean_difference']:+.2f} | {t['fdr']:.2e} |")
    lines += ["", "![scores](figures/module_scores.png)", "", "![heatmap](figures/heatmap.png)", ""]
    if "stage_tests" in summary:
        lines += ["## Association with pathologic stage (tumors)", "",
                  "| Score | n | Kruskal p | Spearman ρ | Spearman p |", "|---|---:|---:|---:|---:|"]
        for t in summary["stage_tests"]:
            lines.append(f"| {t['score']} | {t['n']} | {t['p_kruskal']:.3g} | "
                         f"{t['spearman_rho']:+.2f} | {t['p_spearman']:.3g} |")
        lines.append("")
    if "survival_rheostat" in summary:
        s = summary["survival_rheostat"]
        lines += ["## Overall survival", "",
                  (f"Median split on rheostat index: n = {s['n']}, events = {s['events']}, "
                  f"log-rank p = {s['logrank_p']:.3g}; Cox HR per unit = {s['cox_hr_per_unit']:.2f} "
                  f"(95% CI {s['cox_hr_ci_low']:.2f}–{s['cox_hr_ci_high']:.2f}, p = {s['cox_p']:.3g})."),
                  "", "![km](figures/km_rheostat.png)", ""]
    lines += ["## Caveats", "",
              "- Transcript abundance is not enzyme activity or lipid level; confirm with lipidomics.",
              "- The rheostat index is an exploratory heuristic, not a validated biomarker.",
              "- Unpaired tests treat tumor and normal as independent; see paired columns where available.", ""]
    (out / "REPORT.md").write_text("\n".join(lines))
