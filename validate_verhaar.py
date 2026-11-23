"""Validate the Verhaar scheme against a reference set.

Usage:
    .venv/Scripts/python.exe validate_verhaar.py out/verhaar_reference.csv

The CSV needs a `smiles` column and a `verhaar` column holding 1-4. Class 5
("not possible to classify") is a legitimate answer from this implementation
but never appears in the reference set, so predicting it always counts as a
miss -- which is the honest way to score it, since the reference did manage to
assign a mechanism.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter

from toxtree import assess, domain_exclusion, parse_structure


def main(path: str, show: int = 30) -> int:
    rows = [r for r in csv.DictReader(open(path, newline="", encoding="utf-8"))
            if r.get("smiles")]
    matrix: Counter[tuple[int, int]] = Counter()
    errors: list[dict] = []
    unparsed = 0
    excluded: list[tuple[str, str]] = []

    for r in rows:
        expected = int(r["verhaar"])
        rec = parse_structure(r["smiles"], r.get("name") or None)
        if not rec.ok:
            unparsed += 1
            continue
        # As in validate_cramer.py: the scheme cannot answer for a metal
        # compound or a carbon-free inorganic, so it is not scored.
        reason = domain_exclusion(rec, for_tree=True)
        if reason is not None:
            excluded.append((r.get("name") or r["smiles"], reason))
            continue
        result = assess(rec, {"verhaar"}).verhaar
        got = int(result.verhaar_class.split()[1])
        matrix[(expected, got)] += 1
        if got != expected:
            errors.append({"name": r.get("name", ""), "smiles": r["smiles"],
                           "expected": expected, "got": got,
                           "reason": result.reason})

    scored = sum(matrix.values())
    correct = sum(n for (e, g), n in matrix.items() if e == g)
    classes = [1, 2, 3, 4, 5]

    print(f"Reference set: {len(rows)} compounds "
          f"({scored} classified, {unparsed} unparsable)\n")
    print("Confusion matrix (rows = reference, columns = this implementation)")
    print(f"{'':>11}" + "".join(f"{'-> ' + str(c):>8}" for c in classes) + f"{'n':>7}")
    for e in classes:
        n = sum(matrix[(e, g)] for g in classes)
        if not n:
            continue
        cells = "".join(f"{matrix[(e, g)]:>8}" for g in classes)
        print(f"  Class {e}  {cells}{n:>7}")

    print(f"\nOverall concordance: {correct}/{scored} = {100 * correct / scored:.1f}%")
    for e in classes:
        n = sum(matrix[(e, g)] for g in classes)
        if n:
            print(f"  Class {e}: {matrix[(e, e)]}/{n} = {100 * matrix[(e, e)] / n:.0f}%")

    unclassified = sum(n for (e, g), n in matrix.items() if g == 5)
    print(f"\n  Answered 'not possible to classify': {unclassified}"
          f" ({100 * unclassified / scored:.0f}% of the set)")

    if errors:
        print(f"\n{'=' * 72}\nDisagreements ({len(errors)}), showing up to {show}"
              f"\n{'=' * 72}")
        pairs = Counter((e["expected"], e["got"]) for e in errors)
        print("\nmost common confusions:")
        for (exp, got), n in pairs.most_common(6):
            print(f"  reference {exp} -> got {got}:  {n}")
        print()
        for err in errors[:show]:
            print(f"- {err['name'][:34] or err['smiles'][:34]}  "
                  f"reference {err['expected']} -> got {err['got']}")
            print(f"    {err['smiles'][:66]}")
            print(f"    because: {err['reason'][:88]}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
