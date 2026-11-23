"""Molecule parsing, descriptor calculation and depiction helpers.

Pure-Python layer over RDKit used by every decision-tree module.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")

ORGANIC_ELEMENTS = {"C", "H", "N", "O", "S", "F", "Cl", "Br", "I", "P"}
CRAMER_ALLOWED = {"C", "H", "N", "O", "S", "F", "Cl", "Br", "I"}

# Alkali and alkaline-earth metals. These are benign as counterions and are
# stripped as such, so one still present in the parent after stripping is
# covalently bound -- calcium cyanamide rather than sodium benzoate -- which
# puts the compound outside the domain of every scheme here.
LIGHT_METALS = {"Li", "Na", "K", "Rb", "Cs", "Mg", "Ca"}

HEAVY_METALS = {
    "Be", "Al", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "Ba", "La", "Hf", "Ta", "W", "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Ra", "Th", "U",
}


@dataclass
class MolRecord:
    """A parsed structure plus the descriptors the rulebases need."""

    identifier: str
    raw_input: str
    mol: Chem.Mol | None
    error: str | None = None
    props: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.mol is not None

    @property
    def smiles(self) -> str:
        return Chem.MolToSmiles(self.mol) if self.mol is not None else ""


# Counterions that Cramer's scheme treats as not altering the class of the
# organic parent: alkali and alkaline-earth cations, ammonium, and the common
# mineral anions. Anything else (organotin, heavy metals) is kept, because the
# metal is then the toxicologically relevant part.
BENIGN_COUNTERIONS = {
    "[Na+]", "[K+]", "[Li+]", "[Ca+2]", "[Mg+2]", "[Zn+2]", "[NH4+]",
    "[Cl-]", "[Br-]", "[I-]", "[F-]", "[OH-]", "Cl", "Br", "I",
    "O=S(=O)(O)O", "OS(=O)(=O)O", "O=C(O)O", "OC(=O)O", "O",
    "[O-]S(=O)(=O)[O-]", "CC(=O)[O-]", "CC(=O)O",
}


def strip_counterions(mol: Chem.Mol) -> tuple[Chem.Mol, list[str]]:
    """Remove benign counterions, returning the parent and what was removed.

    Cramer's scheme classifies the sodium or ammonium salt of an organic acid
    on the organic parent — the counterion is not what carries the toxicity.
    Fragments are only dropped when every one of them is a recognised benign
    counterion and at least one carbon-containing fragment remains.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) < 2:
        return mol, []

    organic, removed = [], []
    for frag in frags:
        smi = Chem.MolToSmiles(frag)
        has_carbon = any(a.GetSymbol() == "C" for a in frag.GetAtoms())
        if smi in BENIGN_COUNTERIONS and not has_carbon:
            removed.append(smi)
        elif smi in BENIGN_COUNTERIONS:
            removed.append(smi)
        else:
            organic.append(frag)

    if not organic or not removed:
        return mol, []

    # Keep the largest remaining fragment; if several survive, this is a true
    # multi-component substance and the scheme does not apply cleanly.
    organic.sort(key=lambda m: m.GetNumHeavyAtoms(), reverse=True)
    parent = organic[0]
    if len(organic) > 1:
        removed.extend(Chem.MolToSmiles(f) for f in organic[1:])
    return parent, removed


def parse_structure(text: str, identifier: str | None = None) -> MolRecord:
    """Parse SMILES, InChI or a molblock into a sanitised RDKit molecule."""
    raw = (text or "").strip()
    ident = identifier or raw
    if not raw:
        return MolRecord(ident, raw, None, "empty input")

    mol = None
    if raw.lower().startswith("inchi="):
        mol = Chem.MolFromInchi(raw, sanitize=False)
    elif "\n" in raw and ("V2000" in raw or "V3000" in raw):
        mol = Chem.MolFromMolBlock(raw, sanitize=False)
    else:
        mol = Chem.MolFromSmiles(raw, sanitize=False)

    if mol is None:
        return MolRecord(ident, raw, None, "could not be parsed as SMILES/InChI/molblock")

    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # pragma: no cover - depends on input
        return MolRecord(ident, raw, None, f"sanitisation failed: {exc}")

    # Structures written with explicit hydrogens ([H]O[H] for water) would
    # otherwise carry them as heavy atoms through every downstream count.
    try:
        mol = Chem.RemoveHs(mol)
    except Exception:
        pass

    parent, removed = strip_counterions(mol)
    if removed:
        try:
            Chem.SanitizeMol(parent)
        except Exception:
            parent, removed = mol, []

    rec = MolRecord(ident, raw, parent)
    rec.props = compute_descriptors(parent)
    rec.props["counterions"] = removed
    rec.props["salt_stripped"] = bool(removed)
    return rec


def parse_many(lines: Iterable[str]) -> list[MolRecord]:
    out: list[MolRecord] = []
    for line in lines:
        line = (line or "").strip()
        if not line or line.startswith("#"):
            continue
        # allow "SMILES<tab or 2+ spaces>name"
        name = None
        for sep in ("\t", "  ", " "):
            if sep in line:
                head, _, tail = line.partition(sep)
                if Chem.MolFromSmiles(head.strip(), sanitize=False) is not None:
                    line, name = head.strip(), tail.strip() or None
                    break
        out.append(parse_structure(line, name))
    return out


def domain_exclusion(rec: "MolRecord", for_tree: bool = False) -> str | None:
    """Why this structure is outside a scheme's applicability domain.

    A metal compound is outside every rulebase here: no SMARTS pattern in the
    library addresses a metal centre, and neither tree asks a question that
    applies to one. Scoring those counts a deliberate refusal as a miss.

    Carbon-free compounds are different, and the distinction is `for_tree`.
    The Cramer and Verhaar trees ask about organic functional groups and
    cannot answer for hydrazine or sodium nitrite at all. The *alert* sets
    can: SA13 is a hydrazine pattern and SA34 a peroxide one. Excluding
    inorganics from an alert run would therefore be ducking -- potassium
    bromate and sodium nitrite are genuine mutagens that this library has no
    pattern for, and they should be counted as the misses they are.

    Benign counterions have already been stripped by the time this is asked,
    so a sodium carboxylate is in domain and sodium nitrite is not.

    Multi-fragment records are deliberately *not* excluded, although the
    README lists mixtures as out of domain. Counting fragments does not
    identify a mixture: the reference files write racemates as two
    enantiomers in one record (DL-menthol, DL-ethionine) and commercial
    isomer blends as their components (toluene diisocyanate), and the rules
    answer the same for each fragment. Excluding those would quietly drop
    genuine actives, so the only exclusions here are the two that are real.
    """
    p = rec.props
    if not rec.ok:
        return "unparsable"
    metals = sorted(set(p.get("elements", [])) & (HEAVY_METALS | LIGHT_METALS))
    if metals:
        return f"metal compound ({', '.join(metals)})"
    if for_tree and not p.get("is_organic"):
        return "inorganic (no carbon)"
    return None


def compute_descriptors(mol: Chem.Mol) -> dict:
    elements = sorted({a.GetSymbol() for a in mol.GetAtoms()})
    # Fragments handed back by GetMolFrags(sanitizeFrags=False) carry no ring
    # perception, and asking such a molecule for NumRings() raises rather than
    # returning zero. Reached by multi-fragment records in an SDF.
    ri = mol.GetRingInfo()
    try:
        ri.NumRings()
    except RuntimeError:
        Chem.FastFindRings(mol)
        ri = mol.GetRingInfo()
    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 3),
        "logp": round(Crippen.MolLogP(mol), 3),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings": ri.NumRings(),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "elements": elements,
        "has_non_cramer_element": bool(set(elements) - CRAMER_ALLOWED),
        "has_heavy_metal": bool(set(elements) & HEAVY_METALS),
        "is_organic": "C" in elements,
        "charge": Chem.GetFormalCharge(mol),
        "fragments": len(Chem.GetMolFrags(mol)),
    }


def has(mol: Chem.Mol, smarts: str) -> bool:
    patt = _pattern(smarts)
    return patt is not None and mol.HasSubstructMatch(patt)


def count(mol: Chem.Mol, smarts: str) -> int:
    patt = _pattern(smarts)
    if patt is None:
        return 0
    return len(mol.GetSubstructMatches(patt, uniquify=True))


_CACHE: dict[str, Chem.Mol | None] = {}


def _pattern(smarts: str) -> Chem.Mol | None:
    if smarts not in _CACHE:
        _CACHE[smarts] = Chem.MolFromSmarts(smarts)
    return _CACHE[smarts]


def to_svg(mol: Chem.Mol, width: int = 420, height: int = 300,
           highlight: list[int] | None = None, dark: bool = False) -> str:
    """Render a molecule as inline SVG, optionally highlighting matched atoms."""
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    if dark:
        opts.setAtomPalette({-1: (0.85, 0.87, 0.9)})
    mc = Draw.PrepareMolForDrawing(mol)
    drawer.DrawMolecule(mc, highlightAtoms=highlight or [])
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
