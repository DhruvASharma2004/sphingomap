"""Curated sphingolipid metabolism gene sets.

Genes are grouped into functional modules that follow the canonical
sphingolipid pathway: de novo ceramide synthesis in the ER, the
sphingomyelin cycle, ceramide catabolism, the sphingosine-1-phosphate (S1P)
signalling axis and the glycosphingolipid (GSL) branch that starts with
glucosylceramide (GlcCer) and lactosylceramide (LacCer).

Symbols are current HGNC symbols. Older expression matrices (e.g. TCGA
releases built on hg19 annotations) sometimes use retired symbols, so
``ALIASES`` maps those back to the current name.
"""
from __future__ import annotations

MODULES: dict[str, dict] = {
    "de_novo_synthesis": {
        "label": "De novo synthesis",
        "description": "Serine palmitoyltransferase -> 3-ketosphinganine -> dihydroceramide -> ceramide",
        "genes": ["SPTLC1", "SPTLC2", "SPTLC3", "KDSR", "CERS1", "CERS2",
                  "CERS3", "CERS4", "CERS5", "CERS6", "DEGS1", "DEGS2"],
    },
    "sphingomyelin_cycle": {
        "label": "Sphingomyelin cycle",
        "description": "Sphingomyelin synthases and sphingomyelinases interconverting SM and ceramide",
        "genes": ["SGMS1", "SGMS2", "SMPD1", "SMPD2", "SMPD3", "SMPD4"],
    },
    "ceramide_catabolism": {
        "label": "Ceramide catabolism",
        "description": "Ceramidases releasing sphingosine",
        "genes": ["ASAH1", "ASAH2", "ACER1", "ACER2", "ACER3"],
    },
    "s1p_axis": {
        "label": "S1P axis",
        "description": "Sphingosine kinases, S1P phosphatases/lyase, S1P transporter and receptors",
        "genes": ["SPHK1", "SPHK2", "SGPP1", "SGPP2", "SGPL1", "SPNS2",
                  "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5"],
    },
    "glycosphingolipid_branch": {
        "label": "Glycosphingolipid branch",
        "description": "GlcCer/LacCer synthesis and lysosomal GSL degradation",
        "genes": ["UGCG", "B4GALT5", "B4GALT6", "ST3GAL5", "A4GALT",
                  "GBA1", "GBA2", "GLA", "GLB1", "HEXA", "HEXB"],
    },
    "transport": {
        "label": "Ceramide transport",
        "description": "ER-to-Golgi ceramide transfer",
        "genes": ["CERT1"],
    },
}

# Retired / alternative symbols seen in public expression matrices.
ALIASES: dict[str, str] = {
    "LASS1": "CERS1", "LASS2": "CERS2", "LASS3": "CERS3",
    "LASS4": "CERS4", "LASS5": "CERS5", "LASS6": "CERS6",
    "GBA": "GBA1", "COL4A3BP": "CERT1", "SIAT9": "ST3GAL5",
    "EDG1": "S1PR1", "EDG5": "S1PR2", "EDG3": "S1PR3",
    "EDG6": "S1PR4", "EDG8": "S1PR5", "FVT1": "KDSR",
    "ASAHL": "ASAH2", "ASAH3": "ACER1", "ASAH3L": "ACER2", "PHCA": "ACER3",
}

# "Sphingolipid rheostat" (Cuvillier et al., Nature 1996): the balance between
# pro-apoptotic ceramide and pro-survival S1P / glycosphingolipids. The index is
# an exploratory, transcript-level proxy — NOT a measurement of lipid levels.
RHEOSTAT = {
    "pro_survival": ["SPHK1", "SPHK2", "S1PR1", "SPNS2", "UGCG",
                     "B4GALT5", "B4GALT6", "ACER2"],
    "pro_ceramide": ["CERS2", "CERS5", "CERS6", "SMPD1", "SMPD3",
                     "SGPP1", "SGPL1", "GBA1"],
}

# Approximate layout of the pathway diagram (x, y) used by plots.pathway_map.
PATHWAY_LAYOUT: dict[str, tuple[float, float]] = {
    "SPTLC1": (0, 6), "SPTLC2": (1, 6), "KDSR": (2, 6),
    "CERS2": (3, 6.5), "CERS5": (3, 5.5), "CERS6": (4, 6.5), "DEGS1": (4, 5.5),
    "SGMS1": (2, 3.5), "SGMS2": (2, 2.5), "SMPD1": (3, 3.5), "SMPD3": (3, 2.5),
    "ASAH1": (5.3, 5.0), "ACER2": (5.3, 3.2),
    "SPHK1": (7.1, 5.0), "SPHK2": (7.1, 3.2),
    "SGPP1": (9.0, 5.2), "SGPL1": (9.0, 2.8), "S1PR1": (8.2, 1.5), "SPNS2": (7.1, 1.5),
    "UGCG": (5, 7.5), "GBA1": (6, 7.5), "B4GALT5": (7, 7.5), "B4GALT6": (8, 7.5),
    "ST3GAL5": (9, 7.5),
}

PATHWAY_NODES: dict[str, tuple[float, float]] = {
    "Serine +\nPalmitoyl-CoA": (0.5, 7.2),
    "Ceramide": (4.5, 4.5),
    "Sphingomyelin": (1.2, 3.0),
    "Sphingosine": (6.2, 4.0),
    "S1P": (8.2, 4.0),
    "GlcCer -> LacCer\n-> GM3": (7.0, 8.4),
}


def all_genes() -> list[str]:
    """Unique gene list across modules, in module order."""
    seen: list[str] = []
    for mod in MODULES.values():
        for g in mod["genes"]:
            if g not in seen:
                seen.append(g)
    return seen


def gene_to_module() -> dict[str, str]:
    return {g: name for name, mod in MODULES.items() for g in mod["genes"]}
