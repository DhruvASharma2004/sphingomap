import json

import numpy as np
import pandas as pd
import pytest

from sphingomap.genesets import ALIASES, MODULES, PATHWAY_LAYOUT, RHEOSTAT, all_genes
from sphingomap.io import harmonize_symbols, patient_id, sample_type, simplify_stage
from sphingomap.scores import module_scores, rheostat_index
from sphingomap.simulate import PLANTED, simulate_cohort, write_demo
from sphingomap.stats import bh_fdr, differential_expression, paired_differential_expression


def test_genesets_consistent():
    genes = all_genes()
    assert len(genes) == len(set(genes))
    assert all(g in genes for g in RHEOSTAT["pro_survival"] + RHEOSTAT["pro_ceramide"])
    assert all(g in genes for g in PATHWAY_LAYOUT)
    assert all(v in genes for v in ALIASES.values())
    assert not set(RHEOSTAT["pro_survival"]) & set(RHEOSTAT["pro_ceramide"])


@pytest.mark.parametrize("barcode,expected", [
    ("TCGA-AA-3511-01", "Tumor"), ("TCGA-AA-3511-11", "Normal"),
    ("TCGA-AA-3511-01A-01R", "Tumor"), ("TCGA-AA-3511-06", "Metastatic"), ("sample_1", "Other")])
def test_sample_type(barcode, expected):
    assert sample_type(barcode) == expected


def test_patient_and_stage():
    assert patient_id("TCGA-AA-3511-11") == "TCGA-AA-3511"
    assert simplify_stage("Stage IIIB") == "III"
    assert simplify_stage("Stage IV") == "IV"
    assert simplify_stage("Stage I") == "I"
    assert simplify_stage(np.nan) is None


def test_alias_harmonization_merges_duplicates():
    df = pd.DataFrame({"s1": [1.0, 3.0, 5.0]}, index=["LASS2", "CERS2", "gba"])
    out = harmonize_symbols(df)
    assert out.loc["CERS2", "s1"] == 2.0
    assert "GBA1" in out.index


def test_bh_fdr_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.005, np.nan])
    q = bh_fdr(p)
    # reference values from R: p.adjust(c(.01,.04,.03,.005), "BH")
    np.testing.assert_allclose(q[:4], [0.02, 0.04, 0.04, 0.02])
    assert np.isnan(q[4])
    assert np.all(np.diff(np.sort(q[:4])) >= 0)


@pytest.fixture(scope="module")
def cohort():
    expr, clin, surv = simulate_cohort(seed=1)
    expr = harmonize_symbols(expr)
    groups = pd.Series({s: sample_type(s) for s in expr.columns})
    return expr, clin, surv, groups


def test_de_recovers_planted_effects(cohort):
    expr, _, _, groups = cohort
    de = differential_expression(expr.loc[[g for g in all_genes() if g in expr.index]], groups)
    for gene, eff in PLANTED.items():
        assert np.sign(de.loc[gene, "log2FC"]) == np.sign(eff)
        assert abs(de.loc[gene, "log2FC"] - eff) < 0.35
    assert de.loc["SPHK1", "fdr"] < 1e-5
    null = [g for g in all_genes() if g not in PLANTED and g in de.index]
    assert (de.loc[null, "fdr"] < 0.05).mean() < 0.15


def test_paired_de(cohort):
    expr, _, _, groups = cohort
    patients = pd.Series({s: patient_id(s) for s in expr.columns})
    res = paired_differential_expression(expr.loc[["SPHK1", "SGPL1", "CERS2"]], groups, patients)
    assert res.attrs["n_pairs"] == 20
    assert res.loc["SPHK1", "log2FC_paired"] > 0.6
    assert res.loc["SGPL1", "log2FC_paired"] < -0.4


def test_scores_and_rheostat(cohort):
    expr, _, _, groups = cohort
    normals = list(groups[groups == "Normal"].index)
    scores, _used = module_scores(expr, reference=normals)
    assert set(scores.columns) <= set(MODULES)
    assert abs(scores.loc[normals].mean().mean()) < 0.1  # centred on normals
    r = rheostat_index(expr, reference=normals)
    assert r[groups == "Tumor"].mean() > r[groups == "Normal"].mean() + 0.5


def test_full_pipeline(tmp_path):
    from sphingomap.pipeline import run
    paths = write_demo(tmp_path / "data", n_tumor=60, n_normal=15, n_pairs=10)
    s = run([paths["expression"]], tmp_path / "out", [paths["clinical"]], [paths["survival"]], label="t")
    assert s["n_tumor"] == 60 and s["n_normal"] == 15
    assert s["missing_genes"] == []  # legacy symbols in the demo are harmonised
    for f in ["REPORT.md", "differential_expression.csv", "sample_scores.csv", "summary.json",
              "figures/volcano.png", "figures/pathway_map.png", "figures/heatmap.png"]:
        assert (tmp_path / "out" / f).exists(), f
    assert json.loads((tmp_path / "out" / "summary.json").read_text())["label"] == "t"


def test_pipeline_without_clinical(tmp_path):
    from sphingomap.pipeline import run
    paths = write_demo(tmp_path / "data", n_tumor=30, n_normal=10, n_pairs=5)
    s = run([paths["expression"]], tmp_path / "out")
    assert "survival_rheostat" not in s and s["n_significant"] >= 1


def test_cli_genes(capsys):
    from sphingomap.cli import main
    assert main(["genes"]) == 0
    assert "SPHK1" in capsys.readouterr().out
