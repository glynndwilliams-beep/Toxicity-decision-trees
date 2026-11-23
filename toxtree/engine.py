"""Runs every rulebase over a structure and collects the results.

Also implements the Kroes et al. (2004) tiered threshold-of-toxicological-
concern (TTC) decision tree, which combines the genotoxicity alerts with the
Cramer class.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem

from . import cramer, verhaar
from .alerts import (
    COMPUTED,
    GENOTOX,
    IRRITATION,
    NONGENOTOX,
    SENSITISATION,
    Alert,
    screen,
)
from .chem import MolRecord, has

SCHEME_LABELS = {
    "cramer": "Cramer rules (oral toxicity class)",
    "verhaar": "Verhaar scheme (aquatic mode of action)",
    "genotox": "Structural alerts for genotoxic carcinogenicity",
    "nongenotox": "Structural alerts for non-genotoxic carcinogenicity",
    "sensitisation": "Skin sensitisation (protein-binding mechanisms)",
    "irritation": "Skin and eye irritation / corrosion",
    "ttc": "Kroes TTC decision tree",
}


@dataclass
class AlertHit:
    code: str
    name: str
    scheme: str
    rationale: str
    atoms: list[int] = field(default_factory=list)


@dataclass
class Assessment:
    record: MolRecord
    cramer: cramer.CramerResult | None = None
    verhaar: verhaar.VerhaarResult | None = None
    genotox: list[AlertHit] = field(default_factory=list)
    nongenotox: list[AlertHit] = field(default_factory=list)
    sensitisation: list[AlertHit] = field(default_factory=list)
    irritation: list[AlertHit] = field(default_factory=list)
    ttc_class: str = ""
    ttc_value_ug_day: float | None = None
    ttc_trace: list[str] = field(default_factory=list)

    # -- derived verdicts ---------------------------------------------------
    @property
    def genotox_verdict(self) -> str:
        return "Alerting" if self.genotox else "No alert"

    @property
    def nongenotox_verdict(self) -> str:
        return "Alerting" if self.nongenotox else "No alert"

    @property
    def sensitisation_verdict(self) -> str:
        if not self.sensitisation:
            return "No alert"
        n = len(self.sensitisation)
        return f"Alerting ({n} mechanism{'s' if n > 1 else ''})"

    @property
    def irritation_verdict(self) -> str:
        n = len(self.irritation)
        if n == 0:
            return "No alert"
        if n >= 3:
            return "Alerting — multiple irritation triggers"
        return "Alerting"

    def all_hits(self) -> list[AlertHit]:
        return self.genotox + self.nongenotox + self.sensitisation + self.irritation

    def to_row(self) -> dict:
        p = self.record.props
        row = {
            "Input": self.record.raw_input,
            "Name": self.record.identifier if self.record.identifier != self.record.raw_input else "",
            "SMILES": self.record.smiles,
            "Formula": p.get("formula", ""),
            "MW": p.get("mw"),
            "logP": p.get("logp"),
            "Cramer class": self.cramer.cramer_class if self.cramer else "",
            "Verhaar class": self.verhaar.verhaar_class if self.verhaar else "",
            "Genotox alerts": "; ".join(h.name for h in self.genotox),
            "Non-genotox alerts": "; ".join(h.name for h in self.nongenotox),
            "Sensitisation alerts": "; ".join(h.name for h in self.sensitisation),
            "Irritation alerts": "; ".join(h.name for h in self.irritation),
            "TTC tier": self.ttc_class,
            "TTC (µg/person/day)": self.ttc_value_ug_day,
        }
        if self.verhaar and self.verhaar.log_lc50_mol_l is not None:
            row["Predicted log LC50 (mol/L)"] = self.verhaar.log_lc50_mol_l
            if p.get("mw"):
                row["Predicted LC50 (mg/L)"] = verhaar.lc50_mg_per_l(
                    self.verhaar.log_lc50_mol_l, p["mw"])
        return row


def _hits(mol: Chem.Mol, alerts: list[Alert], scheme: str) -> list[AlertHit]:
    out: list[AlertHit] = []
    for alert, matches in screen(mol, alerts):
        atoms = sorted({i for m in matches for i in m})
        out.append(AlertHit(alert.code, alert.name, scheme, alert.rationale, atoms))
    for comp in COMPUTED:
        if comp.scheme == scheme and comp.fires(mol):
            out.append(AlertHit(comp.code, comp.name, scheme, comp.rationale, []))
    return out


# --- Kroes TTC tree ---------------------------------------------------------
P_ACHE = "[$([PX4](=[OX1,SX1])([OX2,SX2][#6])[OX2,SX2][#6]),$([NX3;H1,H2][CX3](=[OX1])[OX2][#6])]"
P_AFLATOXIN_LIKE = "[$(O=C1C=COC2=C1C=CC=C2),$(c1cc2ccoc2cc1)]"

TTC_GENOTOX = 0.15      # µg/person/day
TTC_ACHE = 18.0
TTC_CLASS_III = 90.0
TTC_CLASS_II = 540.0
TTC_CLASS_I = 1800.0


def kroes_ttc(mol: Chem.Mol, cram: cramer.CramerResult,
              genotox_hits: list[AlertHit], props: dict) -> tuple[str, float | None, list[str]]:
    trace: list[str] = []

    trace.append("Q1. Is the substance a metal, organometallic, dioxin-like or "
                 "a known high-potency carcinogen class?")
    if props.get("has_heavy_metal"):
        trace.append("  → Yes. The TTC approach is not applicable; a substance-"
                     "specific risk assessment is required.")
        return "Not applicable (metal/organometallic)", None, trace
    if has(mol, P_AFLATOXIN_LIKE):
        trace.append("  → Yes. Aflatoxin-like furocoumarin cohort of concern; "
                     "TTC is not applicable.")
        return "Not applicable (cohort of concern)", None, trace
    trace.append("  → No.")

    trace.append("Q2. Does the structure raise a genotoxicity/DNA-reactivity alert?")
    if genotox_hits:
        names = ", ".join(h.name for h in genotox_hits[:4])
        trace.append(f"  → Yes ({names}). Apply the TTC for potential "
                     "genotoxic carcinogens.")
        return "Potential genotoxic carcinogen", TTC_GENOTOX, trace
    trace.append("  → No.")

    trace.append("Q3. Is it an organophosphate or carbamate (cholinesterase "
                 "inhibitor structural class)?")
    if has(mol, P_ACHE):
        trace.append("  → Yes. Apply the organophosphate/carbamate TTC.")
        return "Organophosphate / carbamate", TTC_ACHE, trace
    trace.append("  → No.")

    trace.append(f"Q4. Apply the Cramer class TTC. Cramer result: {cram.cramer_class}.")
    value = cramer.TTC_UG_PER_PERSON_PER_DAY[cram.cramer_class]
    trace.append(f"  → TTC = {value:g} µg/person/day.")
    return f"Cramer {cram.cramer_class}", value, trace


def assess(record: MolRecord, schemes: set[str] | None = None) -> Assessment:
    """Run the selected rulebases over one structure."""
    schemes = schemes or set(SCHEME_LABELS)
    a = Assessment(record=record)
    if not record.ok:
        return a
    mol, props = record.mol, record.props

    if "genotox" in schemes or "ttc" in schemes:
        a.genotox = _hits(mol, GENOTOX, "genotox")
    if "nongenotox" in schemes:
        a.nongenotox = _hits(mol, NONGENOTOX, "nongenotox")
    if "sensitisation" in schemes:
        a.sensitisation = _hits(mol, SENSITISATION, "sensitisation")
    if "irritation" in schemes:
        a.irritation = _hits(mol, IRRITATION, "irritation")
    if "cramer" in schemes or "ttc" in schemes:
        a.cramer = cramer.classify(mol, props)
    if "verhaar" in schemes:
        a.verhaar = verhaar.classify(mol, props)
    if "ttc" in schemes and a.cramer is not None:
        a.ttc_class, a.ttc_value_ug_day, a.ttc_trace = kroes_ttc(
            mol, a.cramer, a.genotox, props)
    if "genotox" not in schemes:
        a.genotox = []
    return a


def assess_many(records: list[MolRecord], schemes: set[str] | None = None) -> list[Assessment]:
    return [assess(r, schemes) for r in records]
