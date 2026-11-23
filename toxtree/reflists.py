"""Optional curated structure lists used by the Cramer tree.

Cramer's scheme asks two membership questions that are not structural rules:
whether a substance is a normal constituent of the body, and whether it is a
common component of food. Those are lookups against curated lists.

The lists live in `toxtree/data/*.csv` and are fetched by
`tools/fetch_reference_lists.py`. They are not committed, so everything here
degrades gracefully: with no data files present the module falls back to the
small built-in metabolite list and reports that it is doing so.

Matching is on canonical SMILES with stereochemistry removed, because query
structures are routinely drawn without stereo descriptors while the reference
lists carry them.
"""
from __future__ import annotations

import csv
import pathlib

from rdkit import Chem

DATA_DIR = pathlib.Path(__file__).parent / "data"


def canonical(mol_or_smiles) -> str | None:
    """Canonical SMILES with stereochemistry stripped, or None."""
    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    else:
        mol = mol_or_smiles
    if mol is None:
        return None
    try:
        copy = Chem.Mol(mol)
        Chem.RemoveStereochemistry(copy)
        return Chem.MolToSmiles(copy)
    except Exception:
        return None


def _load(stem: str) -> dict[str, str]:
    path = DATA_DIR / f"{stem}.csv"
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    try:
        with path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                key = canonical(row.get("smiles", ""))
                if key:
                    out.setdefault(key, row.get("name", ""))
    except Exception:
        return {}
    return out


BODY: dict[str, str] = _load("bodymol")
FOOD: dict[str, str] = _load("foodmol")

# ---------------------------------------------------------------------------
# Cramer's first question asks for a *normal constituent of the body*, and the
# observed answers show it to be far narrower than "anything the body
# contains". Across 232 recorded decision paths, 65 compounds answer it yes --
# every one of them Class I -- and they form one tight, describable group: the
# building blocks of intermediary metabolism.
#
#   the twenty proteinogenic amino acids
#   the five nucleobases and four common nucleosides
#   hexoses and their uronic acids
#   cholesterol, water, acetone, phosphoenolpyruvate
#
# Vitamins (biotin, ascorbic acid, B6), fatty acids, citric and carbonic acid,
# the sugar alcohols and the carotenes all answer *no*, even though every one
# of them appears in the larger body-constituent file. Using that file for this
# question is what made an explicit guard necessary; with the right membership
# the question needs no guard, because cysteine and methionine are in the set
# while biotin and riboflavin simply are not.
#
# The membership below was identified from those observed answers, so the
# improvement it produces on those same sets is partly self-confirming. The
# structures are standard ones written out here rather than copied.
# ---------------------------------------------------------------------------
_PRIMARY_METABOLITE_SMILES = [
    # the twenty proteinogenic amino acids
    "NCC(=O)O", "CC(N)C(=O)O", "CC(C)C(N)C(=O)O", "CC(C)CC(N)C(=O)O",
    "CCC(C)C(N)C(=O)O", "OC(=O)C1CCCN1", "OCC(N)C(=O)O", "CC(O)C(N)C(=O)O",
    "SCC(N)C(=O)O", "CSCCC(N)C(=O)O", "NC(=O)CC(N)C(=O)O",
    "NC(=O)CCC(N)C(=O)O", "OC(=O)CC(N)C(=O)O", "OC(=O)CCC(N)C(=O)O",
    "NCCCCC(N)C(=O)O", "NC(=N)NCCCC(N)C(=O)O", "NC(Cc1cnc[nH]1)C(=O)O",
    "NC(Cc1ccccc1)C(=O)O", "NC(Cc1ccc(O)cc1)C(=O)O",
    "NC(Cc1c[nH]c2ccccc12)C(=O)O",
    # nucleobases
    "Nc1ncnc2[nH]cnc12", "Nc1nc2[nH]cnc2c(=O)[nH]1", "Nc1cc[nH]c(=O)n1",
    "Cc1c[nH]c(=O)[nH]c1=O", "O=c1cc[nH]c(=O)[nH]1",
    # nucleosides
    "Nc1ncnc2c1ncn2C1OC(CO)C(O)C1O",
    "Nc1nc2c(ncn2C2OC(CO)C(O)C2O)c(=O)[nH]1",
    "Nc1ccn(C2OC(CO)C(O)C2O)c(=O)n1",
    "O=c1ccn(C2OC(CO)C(O)C2O)c(=O)[nH]1",
    # hexoses and uronic acids (stereochemistry is stripped before matching,
    # so one entry covers the glucose/galactose/mannose series)
    "OCC1OC(O)C(O)C(O)C1O", "OC1OC(C(=O)O)C(O)C(O)C1O",
    # other primary constituents
    "CC(C)CCCC(C)C1CCC2C1(C)CCC1C2CC=C2CC(O)CCC12C",   # cholesterol
    "O", "CC(C)=O", "OC(=O)C(=C)OP(=O)(O)O",
]

PRIMARY_METABOLITES: set[str] = {
    key for key in (canonical(s) for s in _PRIMARY_METABOLITE_SMILES) if key
}


_TAUTOMER_SET: set[str] | None = None
_TAUTOMER_CACHE: dict[str, str] = {}


def _tautomer_key(mol_or_smiles) -> str | None:
    """Canonical tautomer, so imidazole and nucleobase tautomers compare equal.

    Histidine is drawn as both the 4H and 5H imidazole tautomer, and the
    nucleobases have several accepted forms; plain canonical SMILES treat
    those as different molecules.
    """
    plain = canonical(mol_or_smiles)
    if plain is None:
        return None
    if plain not in _TAUTOMER_CACHE:
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            mol = Chem.MolFromSmiles(plain)
            enumerator = rdMolStandardize.TautomerEnumerator()
            _TAUTOMER_CACHE[plain] = Chem.MolToSmiles(
                enumerator.Canonicalize(mol))
        except Exception:
            _TAUTOMER_CACHE[plain] = plain
    return _TAUTOMER_CACHE[plain]


def primary_metabolite(mol) -> bool:
    """True for a building block of intermediary metabolism (Cramer's Q1)."""
    key = canonical(mol)
    if key and key in PRIMARY_METABOLITES:
        return True
    global _TAUTOMER_SET
    if _TAUTOMER_SET is None:
        _TAUTOMER_SET = {k for k in (_tautomer_key(s)
                                     for s in PRIMARY_METABOLITES) if k}
    tkey = _tautomer_key(mol)
    return bool(tkey) and tkey in _TAUTOMER_SET


def body_constituent(mol) -> str | None:
    """Return the reference name if the structure is a normal body constituent."""
    key = canonical(mol)
    return BODY.get(key) if key else None


def food_component(mol) -> str | None:
    """Return the reference name if the structure is a listed food component."""
    key = canonical(mol)
    return FOOD.get(key) if key else None


def status() -> dict:
    """Describe which lists are loaded, for display and for the method notes."""
    return {
        "body_loaded": bool(BODY),
        "food_loaded": bool(FOOD),
        "body_count": len(BODY),
        "food_count": len(FOOD),
        "data_dir": str(DATA_DIR),
    }


def available() -> bool:
    return bool(BODY or FOOD)
