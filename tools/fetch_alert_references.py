"""Fetch the experimental datasets for the structural-alert rulebases.

Two sets, both distributed as test data with the Toxtree project:

  LLNA   209 compounds with a local lymph node assay result. Sensitisers carry
         an EC3 potency value; non-sensitisers are marked NC. Each also carries
         the mechanistic reaction domain it was assigned, which can be compared
         against the mechanism this implementation's alerts name.

  ISSCAN 925 records with Salmonella (Ames) mutagenicity calls and a
         carcinogenicity call in the same encoding: **3 is active, 1 is
         inactive**, 2 marginal, ND not determined. That direction is worth
         stating plainly -- assuming the intuitive 1=positive inverts every
         label and produces a validation that looks catastrophic but means
         nothing. It is confirmed by the compounds themselves: benzo[a]pyrene
         and benzidine are 3, thiourea and caffeine are 1.

    python tools/fetch_alert_references.py out/
    python validate_alerts.py out/llna_reference.csv --scheme sensitisation
    python validate_alerts.py out/ames_reference.csv --scheme genotox

Neither file is committed; both are GPL-licensed where they are distributed.
"""
from __future__ import annotations

import csv
import io
import pathlib
import sys
import urllib.request
from collections import Counter

BASE = ("https://raw.githubusercontent.com/ideaconsult/toxtree-mirror/HEAD/"
        "toxtree/toxtree-plugins")
LLNA_URL = (f"{BASE}/toxtree-michaelacceptors/src/test/resources/"
            "toxtree/michaelacceptors/LLNA_3D.sdf")
AMES_URL = f"{BASE}/toxtree-ames/src/test/resources/data/iss2_rid.sdf"


def _download(url: str) -> bytes | None:
    """Fetch and normalise to UTF-8.

    These SDFs carry cp1252 bytes in some chemical names, which makes RDKit's
    property reader raise on an otherwise perfectly good record.
    """
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            raw = resp.read()
    except Exception as exc:
        print(f"  FAILED {url.rsplit('/', 1)[-1]}: {exc}")
        return None
    try:
        raw.decode("utf-8")
    except UnicodeDecodeError:
        raw = raw.decode("cp1252", errors="replace").encode("utf-8")
    return raw


def main(outdir: str) -> int:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    written = 0

    # --- skin sensitisation -------------------------------------------------
    raw = _download(LLNA_URL)
    if raw:
        rows = []
        for mol in Chem.ForwardSDMolSupplier(io.BytesIO(raw), sanitize=True,
                                             removeHs=True):
            if mol is None:
                continue
            props = mol.GetPropsAsDict()
            ec3 = str(props.get("LLNA EC3 %", "")).strip()
            if not ec3:
                continue
            # a numeric EC3 is a sensitiser; NC means not classified as one
            active = 0 if ec3.upper() == "NC" else 1
            rows.append({
                "smiles": Chem.MolToSmiles(mol),
                "active": active,
                "name": str(props.get("Names", ""))[:60],
                "domain": str(props.get("ReactionDomain", "")).strip(),
                "potency": "" if active == 0 else ec3,
            })
        target = out / "llna_reference.csv"
        with target.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["smiles", "active", "name",
                                               "domain", "potency"])
            w.writeheader()
            w.writerows(rows)
        counts = Counter(r["active"] for r in rows)
        print(f"  LLNA: {len(rows)} compounds -> {target}")
        print(f"        {counts[1]} sensitisers, {counts[0]} non-sensitisers")
        print(f"        domains: {dict(Counter(r['domain'] for r in rows).most_common(8))}")
        written += len(rows)

    # --- Ames mutagenicity --------------------------------------------------
    raw = _download(AMES_URL)
    if raw:
        rows = []
        for mol in Chem.ForwardSDMolSupplier(io.BytesIO(raw), sanitize=True,
                                             removeHs=True):
            if mol is None:
                continue
            props = mol.GetPropsAsDict()
            call = str(props.get("SAL_iss2", "")).strip()
            if call not in {"1", "3"}:      # skip marginal (2) and ND
                continue
            rows.append({
                "smiles": Chem.MolToSmiles(mol),
                "active": 1 if call == "3" else 0,
                "name": str(props.get("ChemName_iss2", ""))[:60],
                "carcinogenicity": str(props.get("Canc_iss2", "")).strip(),
            })
        target = out / "ames_reference.csv"
        with target.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["smiles", "active", "name",
                                               "carcinogenicity"])
            w.writeheader()
            w.writerows(rows)
        counts = Counter(r["active"] for r in rows)
        print(f"  Ames: {len(rows)} compounds -> {target}")
        print(f"        {counts[1]} mutagens, {counts[0]} non-mutagens")
        written += len(rows)

    return 0 if written else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "out"))
