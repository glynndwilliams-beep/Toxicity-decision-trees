"""Fetch the Verhaar aquatic mode-of-action validation set.

121 compounds with the mode-of-action class assigned in the Verhaar et al.
evaluation, distributed as test data with the Toxtree project. Structures are
held as molblocks rather than SMILES, so RDKit is used to convert them.

    python tools/fetch_verhaar_reference.py out/
    python validate_verhaar.py out/verhaar_reference.csv

As with the Cramer sets, the data is not committed -- it is GPL-licensed where
it is distributed -- so fetch it locally when you want to re-validate.
"""
from __future__ import annotations

import csv
import pathlib
import sys
import urllib.request
from collections import Counter

URL = ("https://raw.githubusercontent.com/ideaconsult/toxtree-mirror/HEAD/"
       "toxtree/toxtree-plugins/toxtree-verhaar2/src/test/resources/"
       "toxtree/plugins/verhaar2/Verhaar2000.sdf")


def main(outdir: str) -> int:
    import io

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    try:
        with urllib.request.urlopen(URL, timeout=60) as resp:
            raw = resp.read()
    except Exception as exc:
        print(f"FAILED: {exc}")
        return 1

    # RDKit's SDF reader is used rather than splitting the text by hand: a
    # molblock's first three lines are header, so naive splitting on the
    # record separator misaligns the counts line and almost nothing parses.
    rows, skipped = [], 0
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(raw), sanitize=True,
                                         removeHs=True)
    for mol in supplier:
        if mol is None:
            skipped += 1
            continue
        props = mol.GetPropsAsDict()
        try:
            number = int(float(props.get("Class")))
        except (TypeError, ValueError):
            skipped += 1
            continue
        try:
            logkow = f"{float(props.get('log Kowa')):.3f}"
        except (TypeError, ValueError):
            logkow = ""
        rows.append({"smiles": Chem.MolToSmiles(mol), "verhaar": number,
                     "name": str(props.get("Name", "")),
                     "cas": str(props.get("CAS", "")), "log_kow": logkow})

    target = out / "verhaar_reference.csv"
    with target.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["smiles", "verhaar", "name",
                                                "cas", "log_kow"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"{len(rows)} compounds -> {target}"
          + (f"  ({skipped} skipped)" if skipped else ""))
    print("by class:", dict(sorted(Counter(r["verhaar"] for r in rows).items())))
    return 0 if rows else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "out"))
