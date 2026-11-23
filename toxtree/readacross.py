"""Read-across neighbours and an evidential-support index.

Every verdict this app gives is a rule firing, with no sense of whether the
molecule is one the rules have ever been checked against. A Cramer Class III
for cinnamaldehyde and a Class III for something unlike anything in any
reference set read identically, and they should not.

This module closes that gap by looking up the query's nearest neighbours in
the reference sets used for validation and asking whether their *experimental*
outcomes agree with what the rules just said.

A note on what the index means
------------------------------
The idea is borrowed from VEGA's applicability domain index, but the meaning
is deliberately different. VEGA scores statistical models, where the index
asks whether a prediction is an interpolation or an extrapolation. The
decision trees here are deterministic rules: they will answer identically for
a molecule in the reference set and one unlike anything ever tested, and the
rules are not fitted to these compounds at all.

So this is **not** an applicability domain in the QSAR sense. It measures
*evidential support*: how much checked evidence stands behind this particular
answer. A low score does not mean the rule misfired -- it means nothing in the
reference data would have caught it if it had.

The components are reported separately, because the single number is the least
informative part of it.
"""
from __future__ import annotations

import csv
import pathlib
from dataclasses import dataclass, field

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

DATA_DIR = pathlib.Path(__file__).parent / "data"

# Each entry: file, the column holding the outcome, and how to render it.
SOURCES: dict[str, dict] = {
    "cramer": {
        "files": ["cramer_reference.csv", "foodmol_validation.csv"],
        "column": "cramer",
        "label": "Cramer class",
        "render": lambda v: {"1": "Class I", "2": "Class II",
                             "3": "Class III"}.get(str(v), str(v)),
        "origin": "Cramer et al. (1978) Appendix 1 and the reference implementation",
    },
    "verhaar": {
        "files": ["verhaar_reference.csv"],
        "column": "verhaar",
        "label": "Verhaar class",
        "render": lambda v: f"Class {v}",
        "origin": "Verhaar et al.",
    },
    "genotox": {
        "files": ["ames_reference.csv"],
        "column": "active",
        "label": "Ames result",
        "render": lambda v: "mutagenic" if str(v) == "1" else "not mutagenic",
        "origin": "ISSCAN Salmonella (Ames) data",
    },
    "sensitisation": {
        "files": ["llna_reference.csv"],
        "column": "active",
        "label": "LLNA result",
        "render": lambda v: "sensitiser" if str(v) == "1" else "non-sensitiser",
        "origin": "local lymph node assay data",
    },
}

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


@dataclass
class Neighbour:
    similarity: float
    name: str
    smiles: str
    outcome: str
    raw: str


@dataclass
class Support:
    """Evidential support for one verdict."""

    endpoint: str
    available: bool
    n_reference: int = 0
    neighbours: list[Neighbour] = field(default_factory=list)
    concordance: float | None = None
    similarity: float = 0.0
    in_range: bool = True
    novelty: float = 0.0
    index: float = 0.0
    notes: list[str] = field(default_factory=list)

    @property
    def band(self) -> str:
        if not self.available:
            return "no reference data"
        if self.index >= 0.70:
            return "well supported"
        if self.index >= 0.45:
            return "partly supported"
        return "little support"


class _Reference:
    def __init__(self, endpoint: str, spec: dict):
        self.endpoint = endpoint
        self.spec = spec
        self.rows: list[dict] = []
        self.fps: list = []
        self.mw: list[float] = []
        self.logp: list[float] = []
        self.bits: set[int] = set()
        self._load()

    def _load(self) -> None:
        from .chem import compute_descriptors

        seen = set()
        for filename in self.spec["files"]:
            path = DATA_DIR / filename
            if not path.exists():
                continue
            try:
                with path.open(newline="", encoding="utf-8") as fh:
                    for row in csv.DictReader(fh):
                        smiles = (row.get("smiles") or "").strip()
                        value = row.get(self.spec["column"])
                        if not smiles or value in (None, ""):
                            continue
                        mol = Chem.MolFromSmiles(smiles)
                        if mol is None:
                            continue
                        key = Chem.MolToSmiles(mol)
                        if key in seen:
                            continue
                        seen.add(key)
                        fp = _GEN.GetFingerprint(mol)
                        props = compute_descriptors(mol)
                        self.rows.append({"smiles": key, "raw": str(value),
                                          "name": (row.get("name") or "").strip()})
                        self.fps.append(fp)
                        self.mw.append(props["mw"])
                        self.logp.append(props["logp"])
                        self.bits.update(fp.GetOnBits())
            except Exception:
                continue


_CACHE: dict[str, _Reference] = {}


def _reference(endpoint: str) -> _Reference | None:
    if endpoint not in SOURCES:
        return None
    if endpoint not in _CACHE:
        _CACHE[endpoint] = _Reference(endpoint, SOURCES[endpoint])
    ref = _CACHE[endpoint]
    return ref if ref.rows else None


def available_endpoints() -> list[str]:
    return [e for e in SOURCES if _reference(e) is not None]


def support(mol: Chem.Mol, endpoint: str, predicted: str | None = None,
            k: int = 5) -> Support:
    """Nearest reference neighbours and an evidential-support index."""
    ref = _reference(endpoint)
    if ref is None:
        return Support(endpoint=endpoint, available=False,
                       notes=["No reference data installed for this endpoint."])

    fp = _GEN.GetFingerprint(mol)
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref.fps)
    order = sorted(range(len(sims)), key=lambda i: -sims[i])[:k]

    render = ref.spec["render"]
    neighbours = [
        Neighbour(similarity=round(sims[i], 3),
                  name=ref.rows[i]["name"] or ref.rows[i]["smiles"][:28],
                  smiles=ref.rows[i]["smiles"],
                  outcome=render(ref.rows[i]["raw"]),
                  raw=ref.rows[i]["raw"])
        for i in order
    ]

    s = Support(endpoint=endpoint, available=True, n_reference=len(ref.rows),
                neighbours=neighbours,
                similarity=round(max(sims) if sims else 0.0, 3))

    # Concordance: do the neighbours' measured outcomes match the verdict,
    # weighted by how similar each one actually is? VEGA calls the equivalent
    # check the most critical of its components, and the same holds here.
    if predicted is not None and neighbours:
        weight = sum(n.similarity for n in neighbours)
        if weight > 0:
            agree = sum(n.similarity for n in neighbours
                        if n.outcome.lower() == predicted.lower())
            s.concordance = round(agree / weight, 3)

    # Descriptor range, and how much of the molecule is structurally new.
    from .chem import compute_descriptors
    props = compute_descriptors(mol)
    if ref.mw:
        s.in_range = (min(ref.mw) <= props["mw"] <= max(ref.mw)
                      and min(ref.logp) <= props["logp"] <= max(ref.logp))
    on = set(fp.GetOnBits())
    s.novelty = round(len(on - ref.bits) / len(on), 3) if on else 0.0

    # The weighting is this module's own, not VEGA's: concordance counts most,
    # then similarity, with range and novelty as smaller corrections.
    concordance = s.concordance if s.concordance is not None else s.similarity
    s.index = round(0.40 * concordance + 0.30 * s.similarity
                    + 0.15 * (1.0 if s.in_range else 0.0)
                    + 0.15 * (1.0 - s.novelty), 3)

    if s.similarity < 0.75:
        s.notes.append(
            f"Nearest reference compound is only {s.similarity:.2f} similar, "
            "so the comparison is weak.")
    if not s.in_range:
        s.notes.append("Molecular weight or logP falls outside the range of "
                       "the reference set.")
    if s.novelty > 0.3:
        s.notes.append(f"{100 * s.novelty:.0f}% of the molecule's structural "
                       "fragments do not appear anywhere in the reference set.")
    if s.concordance is not None and s.concordance < 0.5:
        s.notes.append("The most similar reference compounds were measured to "
                       "behave differently from this verdict.")
    return s


def status() -> dict:
    out = {}
    for endpoint in SOURCES:
        ref = _reference(endpoint)
        out[endpoint] = len(ref.rows) if ref else 0
    return out
