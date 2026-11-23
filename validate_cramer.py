"""Validate the Cramer implementation against a reference set.

Usage:
    .venv/Scripts/python.exe validate_cramer.py path/to/cramer_reference.csv

The CSV needs a `smiles` column and a `cramer` column holding 1, 2 or 3.
A `name` column is used in the error report if present.

Prints a confusion matrix, the overall and per-class concordance, and every
disagreement with the decision path that produced it, so each one can be judged
on its merits rather than counted.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict

from toxtree import assess, parse_structure
from toxtree.cramer import CLASS_I, CLASS_II, CLASS_III

CLASS_OF = {CLASS_I: 1, CLASS_II: 2, CLASS_III: 3}


def load(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh) if r.get("smiles")]


def main(path: str, show_errors: int = 25) -> int:
    rows = load(path)
    matrix: Counter[tuple[int, int]] = Counter()
    errors: list[dict] = []
    unparsed: list[dict] = []

    for r in rows:
        expected = int(r["cramer"])
        rec = parse_structure(r["smiles"], r.get("name") or None)
        if not rec.ok:
            unparsed.append({"smiles": r["smiles"], "name": r.get("name", ""),
                             "error": rec.error})
            continue
        got = CLASS_OF[assess(rec, {"cramer"}).cramer.cramer_class]
        matrix[(expected, got)] += 1
        if got != expected:
            errors.append({
                "name": r.get("name", ""), "smiles": r.get("smiles"),
                "expected": expected, "got": got,
                "trace": assess(rec, {"cramer"}).cramer.trace,
            })

    scored = sum(matrix.values())
    correct = sum(n for (e, g), n in matrix.items() if e == g)

    print(f"Reference set: {len(rows)} compounds "
          f"({scored} classified, {len(unparsed)} unparsable)\n")

    print("Confusion matrix (rows = reference, columns = this implementation)")
    print(f"{'':>12}" + "".join(f"{'-> ' + str(c):>10}" for c in (1, 2, 3)) + f"{'n':>8}")
    for e in (1, 2, 3):
        n = sum(matrix[(e, g)] for g in (1, 2, 3))
        cells = "".join(f"{matrix[(e, g)]:>10}" for g in (1, 2, 3))
        print(f"  Class {e}   {cells}{n:>8}")

    print(f"\nOverall concordance: {correct}/{scored} = {100 * correct / scored:.1f}%")
    for e in (1, 2, 3):
        n = sum(matrix[(e, g)] for g in (1, 2, 3))
        if n:
            print(f"  Class {e}: {matrix[(e, e)]}/{n} = {100 * matrix[(e, e)] / n:.0f}%")

    # Conservatism matters as much as accuracy: under-calling is the unsafe error.
    under = sum(n for (e, g), n in matrix.items() if g < e)
    over = sum(n for (e, g), n in matrix.items() if g > e)
    print(f"\n  Under-called (less strict than reference): {under}"
          f"   <- the unsafe direction")
    print(f"  Over-called  (more strict than reference): {over}")

    if unparsed:
        print(f"\nUnparsable ({len(unparsed)}):")
        for u in unparsed[:10]:
            print(f"  {u['name'][:18].ljust(18)} {u['error']}  {u['smiles'][:50]}")

    if errors:
        print(f"\n{'=' * 72}\nDisagreements ({len(errors)}), "
              f"showing up to {show_errors}\n{'=' * 72}")
        by_question: Counter[str] = Counter()
        def deciding_line(trace):
            # the decision is the last line that concludes with a class
            for line in reversed(trace):
                if "⇒" in line:
                    return line
            return trace[-1] if trace else "?"
        for err in errors:
            by_question[deciding_line(err["trace"]).split(".")[0]] += 1
        print("\nDeciding question for each disagreement:")
        for q, n in by_question.most_common():
            print(f"  {q}: {n}")
        print()
        for err in errors[:show_errors]:
            print(f"- {err['name'] or err['smiles'][:40]}  "
                  f"reference {err['expected']} -> got {err['got']}")
            print(f"    {err['smiles'][:72]}")
            print(f"    decided by: {deciding_line(err['trace'])}")

    return 0 if scored else 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
