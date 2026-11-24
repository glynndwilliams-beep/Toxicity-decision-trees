"""Build the Toxtree regression oracle that test_rules.py checks against.

The test suite pins about thirty verdicts by hand, and twice those pins have
turned out to encode a guess rather than a fact -- acrolein was pinned to
Cramer Class III on an assumption the authors contradict, and methyl acrylate
was pinned as genotoxicity-alerting on an alert that measurement later showed
was in the wrong rulebase. A hand-written expectation cannot catch that,
because it shares the author's assumptions.

This writes a fixture holding, for every compound in the reference files,
what Toxtree says and what this implementation said at the time the fixture
was built. `test_against_the_toxtree_oracle` then recomputes our side and
fails on any compound whose verdict has moved, naming it and saying whether
it moved toward Toxtree or away. Changing a rule therefore means
regenerating this file deliberately, which is the step that was missing.

Usage:

    python tools/build_oracle.py <refdata-dir> [toxtree/data/toxtree_oracle.csv]

<refdata-dir> holds the SDFs fetched from the Toxtree project (see the README
section "Comparing against Toxtree"). Toxtree is located as compare_toxtree.py
locates it: --jar/--java, else $TOXTREE_JAR/$TOXTREE_JAVA, else $TOXTREE_HOME,
else a Toxtree-*.jar under ~/toxtree-java.

The output is GPL-derived, like everything else the tools/ scripts produce, so
it is written into toxtree/data/ and is not committed.
"""
from __future__ import annotations

import csv
import pathlib
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from rdkit import RDLogger  # noqa: E402

from compare_toxtree import (SCHEMES, locate_toxtree,  # noqa: E402
                             run_toxtree, sdf_blocks)
from toxtree import assess, parse_structure  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Which reference file each scheme is scored over. Several schemes share a
# Toxtree module, and the module is the slow part -- mutant.BB_CarcMutRules
# runs semi-empirical QM descriptors per molecule -- so runs are grouped by
# (module, file) and each scheme reads its own field out of the one output.
SOURCES = {
    "cramer": "foodmol.sdf",
    "verhaar": "Verhaar2000.sdf",
    "genotox": "iss2_rid.sdf",
    "nongenotox": "iss2_rid.sdf",
    "sensitisation": "LLNA_3D.sdf",
    "irritation": "LLNA_3D.sdf",
    "eye": "LLNA_3D.sdf",
}


def our_verdict(scheme: str, molblock: str, name: str):
    """What this implementation says, and the SMILES the test will re-read."""
    spec = SCHEMES[scheme]
    rec = parse_structure(molblock, name)
    if not rec.ok:
        return None, None, rec.error
    verdict = spec["ours"](assess(rec, spec["our_schemes"]))
    # The fixture stores SMILES rather than molblocks, so confirm the round
    # trip does not change the answer. Where it does, the compound is
    # dropped: a regression test that cannot reproduce its own baseline is
    # worse than no test.
    round_tripped = parse_structure(rec.smiles, name)
    if not round_tripped.ok:
        return None, None, "SMILES round trip failed to parse"
    if spec["ours"](assess(round_tripped, spec["our_schemes"])) != verdict:
        return None, None, "verdict differs between molblock and SMILES"
    return verdict, rec.smiles, None


def main(argv: list[str]) -> int:
    if not argv:
        raise SystemExit(__doc__)
    refdir = pathlib.Path(argv[0])
    out = pathlib.Path(argv[1] if len(argv) > 1
                       else "toxtree/data/toxtree_oracle.csv")
    jar, java = locate_toxtree(None, None)
    print(f"Toxtree : {jar}\nJava    : {java}\nOutput  : {out}\n")

    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for scheme, filename in SOURCES.items():
        groups[(SCHEMES[scheme]["module"], filename)].append(scheme)

    rows, skipped = [], defaultdict(int)
    for (module, filename), schemes in groups.items():
        src = refdir / filename
        if not src.exists():
            print(f"  skipping {filename}: not found in {refdir}")
            continue
        print(f"  {module}\n    over {filename} for {', '.join(schemes)}")
        # The extension matters: Toxtree picks its output writer from it, so
        # anything that is not a chemical format gives an unreadable file.
        dst = pathlib.Path(tempfile.mkdtemp()) / f"{src.stem}_out.sdf"
        run_toxtree(java, jar, module, src.resolve(), dst, verbose=False)

        produced = list(sdf_blocks(dst))
        records = list(sdf_blocks(src))
        if len(produced) != len(records):
            print(f"    note: Toxtree returned {len(produced)} of "
                  f"{len(records)}; pairing by position is unsafe, so this "
                  f"group is matched on its identifier field")
        by_name = {}
        for mb, f in produced:
            key = (f.get("ChemName_iss2") or f.get("NAME") or f.get("Name")
                   or f.get("Names") or f.get("cdk:Title"))
            if key:
                by_name.setdefault(key, f)

        for molblock, fields in records:
            name = (fields.get("ChemName_iss2") or fields.get("NAME")
                    or fields.get("Name") or fields.get("Names")
                    or fields.get("cdk:Title") or "")
            out_fields = by_name.get(name)
            if out_fields is None:
                skipped["no Toxtree record"] += 1
                continue
            for scheme in schemes:
                spec = SCHEMES[scheme]
                tt = spec["read"](out_fields.get(spec["verdict"], ""))
                ours, smiles, why = our_verdict(scheme, molblock, name)
                if why is not None:
                    skipped[why] += 1
                    continue
                rows.append({"scheme": scheme, "name": name, "smiles": smiles,
                             "toxtree": "" if tt is None else tt,
                             "ours": "" if ours is None else ours})

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["scheme", "name", "smiles",
                                           "toxtree", "ours"])
        w.writeheader()
        w.writerows(rows)

    per = defaultdict(int)
    agree = defaultdict(int)
    for r in rows:
        per[r["scheme"]] += 1
        if r["toxtree"] and r["ours"] and r["toxtree"] == r["ours"]:
            agree[r["scheme"]] += 1
    print(f"\n{len(rows)} rows written to {out}")
    for scheme in sorted(per):
        print(f"  {scheme:15s} {per[scheme]:5d} compounds, "
              f"{agree[scheme]:5d} agree with Toxtree")
    if skipped:
        print("\nnot recorded:")
        for reason, n in sorted(skipped.items(), key=lambda x: -x[1]):
            print(f"  {n:5d}  {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
