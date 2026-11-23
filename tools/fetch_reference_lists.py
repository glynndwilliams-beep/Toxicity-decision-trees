"""Fetch the curated structure lists that the Cramer scheme relies on.

Cramer's tree asks two questions that are not structural rules but membership
tests against curated lists:

  Q2  -- is the substance a normal constituent of the body?
  Q13 -- is it a common component of food?

Those lists are data, not logic, and they are what decides the Class I / Class
II split. They are distributed with the Toxtree project under the GPL, so they
are NOT committed to this repository. Run this script to fetch them locally:

    python tools/fetch_reference_lists.py toxtree/data/

The application works without them -- it falls back to the small built-in
metabolite list and reports reduced confidence on the Class I / II split -- so
this is an optional enhancement rather than a dependency.
"""
from __future__ import annotations

import csv
import pathlib
import sys
import urllib.request
from collections import Counter

BASE = ("https://raw.githubusercontent.com/ideaconsult/toxtree-mirror/HEAD/"
        "toxtree/toxtree-core/src/main/resources")
SOURCES = {
    "bodymol": ("bodymol.sdf", "normal constituents of the body"),
    "foodmol": ("foodmol.sdf", "common components of food"),
}


def sdf_records(text: str):
    """Yield the data fields of each record in an SDF."""
    for block in text.split("$$$$"):
        if not block.strip():
            continue
        fields, key = {}, None
        for line in block.splitlines():
            if line.startswith("> "):
                key = line.split("<", 1)[-1].split(">", 1)[0] if "<" in line else None
            elif key is not None:
                value = line.strip()
                if value:
                    fields.setdefault(key, value)
                else:
                    key = None
        if fields:
            yield fields


def main(outdir: str) -> int:
    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    total = 0

    for stem, (filename, description) in SOURCES.items():
        try:
            with urllib.request.urlopen(f"{BASE}/{filename}", timeout=60) as resp:
                text = resp.read().decode("utf-8", errors="replace")
        except Exception as exc:
            print(f"  FAILED {filename}: {exc}")
            continue

        rows = []
        for fields in sdf_records(text):
            smiles = fields.get("SMILES", "").strip()
            if not smiles:
                continue
            rows.append({
                "smiles": smiles,
                "name": fields.get("NAME", "").strip(),
                "cas": fields.get("CAS", "").strip(),
            })

        target = out / f"{stem}.csv"
        with target.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["smiles", "name", "cas"])
            w.writeheader()
            w.writerows(rows)
        print(f"  {filename}: {len(rows)} structures ({description}) -> {target}")
        total += len(rows)

        # These records also carry the class each structure was assigned, which
        # makes them a second validation set independent of Cramer Appendix 1.
        graded = []
        for fields in sdf_records(text):
            smiles = fields.get("SMILES", "").strip()
            verdict = fields.get("toxTree.tree.cramer.CramerRules", "")
            num = {"I": 1, "II": 2, "III": 3}.get(
                verdict.split("Class ")[-1].rstrip(")").strip())
            if smiles and num:
                graded.append({"smiles": smiles, "cramer": num,
                               "name": fields.get("NAME", "").strip(),
                               "source": filename})
        if graded:
            vtarget = out / f"{stem}_validation.csv"
            with vtarget.open("w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=["smiles", "cramer", "name",
                                                   "source"])
                w.writeheader()
                w.writerows(graded)
            counts = Counter(g["cramer"] for g in graded)
            print(f"     + {len(graded)} graded for validation "
                  f"{dict(sorted(counts.items()))} -> {vtarget}")

    print(f"\n{total} structures written to {out}")
    return 0 if total else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "toxtree/data"))
