"""Synthetic TCGA-format data for tests, CI and offline demos.

The simulated cohort mimics the *shape* of UCSC Xena TCGA files (barcodes,
log2(norm+1) values, clinical and survival columns) so the full pipeline can
be exercised without network access. The planted effects are arbitrary and
must not be interpreted biologically.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .genesets import ALIASES, all_genes

PLANTED = {"SPHK1": 1.2, "UGCG": 0.9, "B4GALT5": 0.8, "ACER2": 0.7, "CERS6": -0.6,
           "SGPL1": -0.9, "S1PR1": -0.8, "SMPD3": -0.7, "DEGS1": 0.5}


def simulate_cohort(n_tumor: int = 120, n_normal: int = 30, n_pairs: int = 20,
                    seed: int = 7, legacy_symbols: bool = True):
    rng = np.random.default_rng(seed)
    genes = all_genes() + [f"CTRL{i:03d}" for i in range(40)]
    reverse = {v: k for k, v in ALIASES.items()}
    tumors = [f"TCGA-SM-{i:04d}-01" for i in range(n_tumor)]
    normals = [f"TCGA-SM-{i:04d}-11" for i in range(n_pairs)] + \
              [f"TCGA-SN-{i:04d}-11" for i in range(n_normal - n_pairs)]
    base = rng.uniform(4, 12, len(genes))
    patient_effect = {p: rng.normal(0, 0.25, len(genes)) for p in range(n_tumor)}
    cols = {}
    rheo = rng.normal(0, 1, n_tumor)  # latent per-tumor rheostat tilt
    for i, s in enumerate(tumors):
        shift = np.array([PLANTED.get(g, 0.0) for g in genes])
        tilt = np.array([0.4 if g in ("SPHK1", "UGCG", "B4GALT5") else
                         -0.4 if g in ("SGPL1", "CERS6") else 0 for g in genes]) * rheo[i]
        cols[s] = base + shift + tilt + patient_effect[i] + rng.normal(0, 0.6, len(genes))
    for j, s in enumerate(normals):
        pe = patient_effect[j] if s.startswith("TCGA-SM") else rng.normal(0, 0.25, len(genes))
        cols[s] = base + pe + rng.normal(0, 0.45, len(genes))
    names = [reverse.get(g, g) if legacy_symbols else g for g in genes]
    expr = pd.DataFrame(cols, index=pd.Index(names, name="sample")).clip(lower=0).round(4)

    stages = rng.choice(["Stage I", "Stage IIA", "Stage IIIB", "Stage IV"], n_tumor, p=[.2, .35, .3, .15])
    clin = pd.DataFrame({"pathologic_stage": stages,
                         "gender": rng.choice(["MALE", "FEMALE"], n_tumor),
                         "age_at_initial_pathologic_diagnosis": rng.integers(35, 90, n_tumor)},
                        index=pd.Index(tumors, name="sampleID"))
    hazard = 0.0004 * np.exp(0.45 * rheo)
    t_event = rng.exponential(1 / hazard)
    t_cens = rng.uniform(200, 3000, n_tumor)
    surv = pd.DataFrame({"OS": (t_event <= t_cens).astype(int),
                         "OS.time": np.minimum(t_event, t_cens).round(0),
                         "_PATIENT": [s[:12] for s in tumors]},
                        index=pd.Index(tumors, name="sample"))
    return expr, clin, surv


def write_demo(outdir: str | Path, **kw) -> dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    expr, clin, surv = simulate_cohort(**kw)
    paths = {"expression": outdir / "DEMO_expression.tsv", "clinical": outdir / "DEMO_clinical.tsv",
             "survival": outdir / "DEMO_survival.tsv"}
    expr.to_csv(paths["expression"], sep="\t")
    clin.to_csv(paths["clinical"], sep="\t")
    surv.to_csv(paths["survival"], sep="\t")
    return paths
