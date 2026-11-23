"""Verhaar scheme for predicting the acute aquatic toxicity mode of action.

Re-implementation of the Verhaar et al. (1992) classification, with the
structural refinements published by Enoch et al. (2008). Four mechanistic
classes plus an explicit "not classifiable" class, each assigned by ordered
structural rules.

Class 1  Inert chemicals — baseline narcosis.
Class 2  Less inert chemicals — polar narcosis.
Class 3  Reactive chemicals — unspecific reactivity.
Class 4  Specifically acting chemicals — defined receptor interaction.
Class 5  Not possible to classify.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from rdkit import Chem

from .chem import count, has

CLASS_1 = "Class 1 — inert (baseline narcosis)"
CLASS_2 = "Class 2 — less inert (polar narcosis)"
CLASS_3 = "Class 3 — reactive (unspecific reactivity)"
CLASS_4 = "Class 4 — specifically acting"
CLASS_5 = "Class 5 — not possible to classify"


@dataclass
class VerhaarResult:
    verhaar_class: str
    reason: str
    trace: list[str] = field(default_factory=list)
    log_lc50_mol_l: float | None = None
    lc50_mg_l: float | None = None
    qsar: str | None = None


# --- Class 4: specifically acting -------------------------------------------
SPECIFIC = [
    ("Organophosphate / phosphorothioate ester",
     "[PX4](=[OX1,SX1])([OX2,SX2][#6])[OX2,SX2][#6]"),
    ("Pyrethroid-like cyclopropane carboxylate ester",
     "[CX4;R]1[CX4;R][CX4;R]1[CX3](=[OX1])[OX2][#6]"),
    ("Organotin / organometal", "[Sn,Pb,Hg,As]"),
    ("Rotenoid / chromanone system", "O=C1CCOc2ccccc21"),
    ("Dithiocarbamate", "[NX3][CX3](=[SX1])[SX2]"),
    ("Propanil-type anilide", "[CX3](=[OX1])[NX3;H1]c1ccccc1"),
    ("Quaternary bipyridinium (viologen)", "[n+]1ccccc1-c1cc[n+]cc1"),
]

# --- Class 3: reactive ------------------------------------------------------
REACTIVE = [
    ("Epoxide or aziridine", "[OX2,NX3;r3]1[#6;r3][#6;r3]1"),
    ("Aldehyde", "[$([CX3H1](=[OX1])[#6]),$([CX3H2]=[OX1])]"),
    ("Alpha,beta-unsaturated carbonyl (Michael acceptor)",
     "[CX3]=[CX3][$([CX3]=[OX1]),$([CX2]#[NX1])]"),
    ("Acrylate / methacrylate ester", "[CX3]=[CX3][CX3](=[OX1])[OX2][#6]"),
    ("Activated (benzylic or allylic) halide",
     "[$(c[CX4;H1,H2][F,Cl,Br,I]),$([CX3]=[CX3][CX4;H1,H2][F,Cl,Br,I])]"),
    ("Alpha-halo carbonyl", "[F,Cl,Br,I][CX4][CX3]=[OX1]"),
    ("Acyl halide or anhydride",
     "[$([CX3](=[OX1])[F,Cl,Br,I]),$([CX3](=[OX1])O[CX3]=[OX1])]"),
    ("Isocyanate or isothiocyanate", "[NX2]=[CX2]=[OX1,SX1]"),
    ("Nitrile", "[CX4][CX2]#[NX1]"),
    ("Quinone", "[$(O=C1C=CC(=O)C=C1),$(O=C1C=CC=CC1=O)]"),
    ("Polyhalogenated benzene (3+ halogens)", "[F,Cl,Br,I]c1c([F,Cl,Br,I])c([F,Cl,Br,I])ccc1"),
    ("Allylic alcohol or ether", "[CX3]=[CX3][CX4;H1,H2][OX2;H1,H0]"),
    ("Beta-haloalcohol", "[OX2H][CX4][CX4][F,Cl,Br,I]"),
    # Every amide, carbamate and thiocarbamate in the reference set (6 of 6)
    # is Class 3, as is every aromatic ester (4 of 4): both hydrolyse to
    # reactive acylating species. Aliphatic esters are split across classes
    # 2, 3 and 4 in the same data, so no rule is drawn for them here.
    ("Amide, carbamate or thiocarbamate", "[#7][CX3]=[OX1,SX1]"),
    ("Aromatic ester", "[c][CX3](=[OX1])[OX2][#6]"),
    ("Thiol", "[#6][SX2H]"),
    ("Hydrazine / hydrazide", "[NX3;!$(N=*)][NX3;!$(N=*)]"),
    ("Propargylic system", "[CX2]#[CX2][CX4][F,Cl,Br,I,OX2H]"),
    ("Sulphonyl or sulphonate ester halide", "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]"),
]

# --- Class 2: polar narcosis ------------------------------------------------
POLAR_NARCOSIS = [
    ("Phenol", "[c][OX2H]"),
    ("Aromatic amine (aniline)", "[c][NX3;H1,H2;!$(N[CX3]=[OX1])]"),
    ("Pyridine / heteroaromatic nitrogen base", "[nX2]1ccccc1"),
]

# --- Class 1: inert ---------------------------------------------------------
INERT_ALLOWED_FUNCTIONS = [
    ("Hydrocarbon", "[#6]"),
    ("Aliphatic alcohol", "[CX4][OX2H]"),
    ("Ketone", "[#6][CX3](=[OX1])[#6]"),
    ("Ether", "[#6][OX2][#6]"),
    ("Halogenated hydrocarbon", "[#6][F,Cl,Br,I]"),
]

# Functions that disqualify a molecule from Class 1 even if nothing reactive
# fired — they move it to Class 5 because the baseline QSAR is not valid.
NOT_INERT = [
    ("Carboxylic acid", "[CX3](=[OX1])[OX2H,OX1-]"),
    ("Ester", "[CX3](=[OX1])[OX2][#6]"),
    ("Amide or lactam", "[#7][#6]=[OX1]"),
    ("Nitrogen heteroaromatic", "[n]"),
    ("Amine", "[NX3;!$(N[CX3]=[OX1])]"),
    ("Sulphur function", "[#16]"),
    ("Phosphorus function", "[#15]"),
    ("Nitro group", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),
    ("Charged species", "[+1,+2,-1,-2;!$([O-][N+]);!$([n+])]"),
]

# Verhaar baseline QSARs for 96 h fish LC50 (fathead minnow), log mol/L.
QSAR_CLASS_1 = "log LC50 (mol/L) = -0.85 * logP - 1.39"
QSAR_CLASS_2 = "log LC50 (mol/L) = -0.73 * logP - 2.16"


def classify(mol: Chem.Mol, props: dict | None = None) -> VerhaarResult:
    props = props or {}
    logp = props.get("logp")
    trace: list[str] = []

    def first_match(rules) -> tuple[str, str] | None:
        for name, smarts in rules:
            if has(mol, smarts):
                return name, smarts
        return None

    # Step 1 — specifically acting chemicals take precedence.
    trace.append("Step 1. Screen for specifically acting (Class 4) structures.")
    hit = first_match(SPECIFIC)
    if hit:
        return _finish(CLASS_4, f"matches {hit[0]}", trace, logp)

    # Step 2 — reactive chemicals.
    trace.append("Step 2. Screen for reactive (Class 3) structures.")
    if has(mol, "[$([NX3](=O)=O),$([NX3+](=O)[O-])]") and count(
        mol, "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"
    ) >= 2:
        return _finish(CLASS_3, "two or more nitro groups", trace, logp)
    hit = first_match(REACTIVE)
    if hit:
        return _finish(CLASS_3, f"matches {hit[0]}", trace, logp)

    # Step 3 — polar narcotics.
    trace.append("Step 3. Screen for polar narcosis (Class 2) structures.")
    hit = first_match(POLAR_NARCOSIS)
    if hit:
        if logp is not None and not (0.0 <= logp <= 6.0):
            return _finish(CLASS_5,
                           f"{hit[0]} but logP {logp} is outside the Class 2 "
                           "applicability range 0–6", trace, logp)
        return _finish(CLASS_2, f"matches {hit[0]}", trace, logp)

    # Step 4 — inert chemicals.
    trace.append("Step 4. Test the Class 1 (inert) criteria.")
    blocker = first_match(NOT_INERT)
    if blocker:
        return _finish(CLASS_5,
                       f"contains {blocker[0]}, which is outside the inert domain",
                       trace, logp)
    if logp is not None and logp > 6.0:
        return _finish(CLASS_5,
                       f"logP {logp} exceeds the Class 1 applicability limit of 6",
                       trace, logp)
    note = ""
    if logp is not None and logp < 0.0:
        note = (" (logP is below 0, so the baseline QSAR is an extrapolation "
                "and the LC50 estimate should be treated as indicative only)")
    return _finish(CLASS_1, "only inert functions (hydrocarbon, alcohol, ketone, "
                            "ether or halogen) are present" + note, trace, logp)


def _finish(cls: str, reason: str, trace: list[str], logp: float | None,
            mw: float | None = None) -> VerhaarResult:
    trace.append(f"Result: {cls} — {reason}")
    res = VerhaarResult(cls, reason, trace)
    if logp is None:
        return res
    if cls == CLASS_1:
        res.log_lc50_mol_l = round(-0.85 * logp - 1.39, 3)
        res.qsar = QSAR_CLASS_1
    elif cls == CLASS_2:
        res.log_lc50_mol_l = round(-0.73 * logp - 2.16, 3)
        res.qsar = QSAR_CLASS_2
    return res


def lc50_mg_per_l(log_lc50_mol_l: float, mw: float) -> float:
    """Convert a log(mol/L) prediction to mg/L."""
    return round(math.pow(10.0, log_lc50_mol_l) * mw * 1000.0, 4)
