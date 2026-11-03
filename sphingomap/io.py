"""Data loading for UCSC Xena TCGA matrices and generic expression tables."""
from __future__ import annotations

import gzip
import shutil
import urllib.request
import warnings
from pathlib import Path

import pandas as pd

from .genesets import ALIASES

XENA_HUB = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download"

# TCGA sample-type codes (barcode characters 14-15).
SAMPLE_TYPES = {"01": "Tumor", "02": "Tumor", "06": "Metastatic", "11": "Normal"}


def xena_urls(cohort: str) -> dict[str, str]:
    """URLs for a TCGA cohort on the UCSC Xena TCGA hub (e.g. 'COAD', 'READ')."""
    c = cohort.upper()
    return {
        "expression": f"{XENA_HUB}/TCGA.{c}.sampleMap%2FHiSeqV2.gz",
        "clinical": f"{XENA_HUB}/TCGA.{c}.sampleMap%2F{c}_clinicalMatrix",
        "survival": f"{XENA_HUB}/survival%2F{c}_survival.txt",
    }


def fetch_cohort(cohort: str, outdir: str | Path, overwrite: bool = False) -> dict[str, Path]:
    """Download expression, clinical and survival tables for a TCGA cohort."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for kind, url in xena_urls(cohort).items():
        dest = outdir / f"{cohort.upper()}_{kind}.tsv"
        if dest.exists() and not overwrite:
            paths[kind] = dest
            continue
        tmp = dest.with_suffix(".download")
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
                shutil.copyfileobj(r, fh)
        except Exception as exc:  # clinical/survival are optional; expression is not
            tmp.unlink(missing_ok=True)
            if kind == "expression":
                raise RuntimeError(f"Could not download {url}: {exc}") from exc
            warnings.warn(f"Skipping {cohort} {kind}: {exc}", stacklevel=2)
            continue
        if url.endswith(".gz"):
            with gzip.open(tmp, "rb") as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst)
            tmp.unlink()
        else:
            tmp.rename(dest)
        paths[kind] = dest
    return paths


def harmonize_symbols(expr: pd.DataFrame) -> pd.DataFrame:
    """Rename retired gene symbols to current HGNC names; average duplicates."""
    expr = expr.rename(index=lambda g: ALIASES.get(str(g).upper(), str(g).upper()))
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    return expr


def load_expression(path: str | Path) -> pd.DataFrame:
    """Load a genes x samples matrix (TSV/CSV, optionally gzipped).

    The first column must hold gene symbols. Values are expected on a log
    scale (Xena HiSeqV2 is log2(RSEM normalized count + 1)).
    """
    path = Path(path)
    sep = "," if path.suffix in {".csv"} or path.name.endswith(".csv.gz") else "\t"
    expr = pd.read_csv(path, sep=sep, index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return harmonize_symbols(expr)


def sample_type(barcode: str) -> str:
    """Map a TCGA barcode (e.g. TCGA-AA-3511-01) to Tumor/Normal/Other."""
    parts = str(barcode).split("-")
    if len(parts) >= 4 and parts[0] == "TCGA":
        return SAMPLE_TYPES.get(parts[3][:2], "Other")
    return "Other"


def patient_id(barcode: str) -> str:
    return "-".join(str(barcode).split("-")[:3])


def load_clinical(path: str | Path | None, survival_path: str | Path | None = None) -> pd.DataFrame:
    """Load Xena clinical + survival tables indexed by sample barcode."""
    frames = []
    if path:
        clin = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
        keep = [c for c in ["pathologic_stage", "gender", "age_at_initial_pathologic_diagnosis",
                            "anatomic_neoplasm_subdivision", "MSI_updated_Oct62011",
                            "CDE_ID_3226963"] if c in clin.columns]
        frames.append(clin[keep])
    if survival_path:
        surv = pd.read_csv(survival_path, sep="\t", index_col=0)
        frames.append(surv[[c for c in ["OS", "OS.time", "PFI", "PFI.time"] if c in surv.columns]])
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, axis=1)
    if "pathologic_stage" in out.columns:
        out["stage_group"] = out["pathologic_stage"].map(simplify_stage)
    return out


def simplify_stage(stage) -> str | None:
    """'Stage IIIB' -> 'III'."""
    if not isinstance(stage, str):
        return None
    s = stage.upper().replace("STAGE", "").strip()
    for roman in ("IV", "III", "II", "I"):
        if s.startswith(roman):
            return roman
    return None
