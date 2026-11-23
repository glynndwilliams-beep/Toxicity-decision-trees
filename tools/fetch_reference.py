"""Rebuild the Cramer validation set.

The reference compounds are those of Cramer, Ford & Hall (1978), Appendix 1,
together with the class each was assigned there. They are distributed as test
fixtures in the Toxtree project, which is GPL-licensed, so the extracted CSV is
deliberately *not* committed to this repository — run this script to produce it
locally when you want to re-validate.

    python tools/fetch_reference.py out/
    python validate_cramer.py out/cramer_reference.csv

Only the data columns are read (structure, assigned class, name). No rule logic
from that project is read or reproduced here.
"""
from __future__ import annotations

import csv
import pathlib
import re
import sys
import urllib.request
from collections import Counter

BASE = ("https://raw.githubusercontent.com/ideaconsult/toxtree-mirror/HEAD/"
        "toxtree/toxtree-plugins/toxtree-cramer/src/test/java/toxTree/cramer")
FILES = ["CramerDataTest.java", "MunroTest.java"]

ENTRY = re.compile(r'\{\s*("(?:[^"\\]|\\.)*"(?:\s*,\s*"(?:[^"\\]|\\.)*")*)\s*\}')
STR = re.compile(r'"((?:[^"\\]|\\.)*)"')


def unescape(s: str) -> str:
    return s.replace('\\\\', '\\').replace('\\"', '"')


def main(outdir: str) -> int:
    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []

    for fn in FILES:
        try:
            with urllib.request.urlopen(f"{BASE}/{fn}", timeout=30) as resp:
                txt = resp.read().decode("utf-8", errors="replace")
        except Exception as exc:
            print(f"  skipped {fn}: {exc}")
            continue
        found = 0
        for m in ENTRY.finditer(txt):
            parts = [unescape(x) for x in STR.findall(m.group(1))]
            if len(parts) < 2:
                continue
            smiles, cls = parts[0].strip(), parts[1].strip()
            if cls not in {"1", "2", "3"} or not smiles:
                continue
            rows.append({
                "smiles": smiles,
                "cramer": int(cls),
                "name": parts[2] if len(parts) > 2 else "",
                "source": fn,
            })
            found += 1
        print(f"  {fn}: {found} entries")

    seen, uniq = set(), []
    for r in rows:
        key = (r["smiles"], r["cramer"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(r)

    target = out / "cramer_reference.csv"
    with target.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["smiles", "cramer", "name", "source"])
        w.writeheader()
        w.writerows(uniq)

    print(f"\n{len(uniq)} unique compounds -> {target}")
    print("by class:", dict(sorted(Counter(r['cramer'] for r in uniq).items())))
    return 0 if uniq else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "refdata"))
