"""Command-line interface: ``sphingomap fetch | demo | run | genes``."""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="sphingomap",
                                description="Sphingolipid pathway expression analysis for TCGA-style cohorts")
    p.add_argument("--version", action="version", version=f"sphingomap {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="download TCGA cohort(s) from the UCSC Xena hub")
    f.add_argument("cohorts", nargs="+", help="e.g. COAD READ")
    f.add_argument("-o", "--outdir", default="data/raw")

    d = sub.add_parser("demo", help="write a synthetic TCGA-format cohort and analyse it")
    d.add_argument("-o", "--outdir", default="results/demo")

    r = sub.add_parser("run", help="analyse expression (+ optional clinical/survival) tables")
    r.add_argument("-e", "--expression", nargs="+", required=True)
    r.add_argument("-c", "--clinical", nargs="*", default=[])
    r.add_argument("-s", "--survival", nargs="*", default=[])
    r.add_argument("-o", "--outdir", default="results")
    r.add_argument("-l", "--label", default="cohort")

    t = sub.add_parser("tcga", help="fetch + run in one step (e.g. sphingomap tcga COAD READ)")
    t.add_argument("cohorts", nargs="+")
    t.add_argument("--datadir", default="data/raw")
    t.add_argument("-o", "--outdir", default=None)

    sub.add_parser("genes", help="print the curated gene sets as JSON")
    a = p.parse_args(argv)

    if a.cmd == "genes":
        from .genesets import MODULES, RHEOSTAT
        print(json.dumps({"modules": MODULES, "rheostat": RHEOSTAT}, indent=2))
        return 0
    if a.cmd == "fetch":
        from .io import fetch_cohort
        for c in a.cohorts:
            for k, v in fetch_cohort(c, a.outdir).items():
                print(f"{c} {k}: {v}")
        return 0

    from .pipeline import run
    if a.cmd == "demo":
        from .simulate import write_demo
        paths = write_demo(f"{a.outdir}/data")
        s = run([paths["expression"]], a.outdir, [paths["clinical"]], [paths["survival"]],
                label="Synthetic demo")
    elif a.cmd == "tcga":
        from .io import fetch_cohort
        got = [fetch_cohort(c, a.datadir) for c in a.cohorts]
        label = "TCGA-" + "+".join(c.upper() for c in a.cohorts)
        s = run([g["expression"] for g in got], a.outdir or f"results/{label}",
                [g.get("clinical") for g in got], [g.get("survival") for g in got], label=label)
    else:
        s = run(a.expression, a.outdir, a.clinical, a.survival, label=a.label)
    print(f"Done: {s['n_tumor']} tumor / {s['n_normal']} normal, "
          f"{s['n_significant']} significant pathway genes. Report in {a.outdir or 'results/'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
