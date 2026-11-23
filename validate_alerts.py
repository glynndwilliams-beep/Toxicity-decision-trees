"""Validate a structural-alert rulebase against experimental data.

    .venv/Scripts/python.exe validate_alerts.py out/llna_reference.csv --scheme sensitisation
    .venv/Scripts/python.exe validate_alerts.py out/ames_reference.csv --scheme genotox

Structural alerts are screening filters, not balanced classifiers: they are
built to miss as few actives as possible, and they buy that with false
positives. So the headline number here is **sensitivity** (what fraction of
actives are flagged), reported alongside specificity and balanced accuracy
rather than a single concordance figure that would hide the trade-off.

With `--per-alert`, each alert's own precision is listed, which is what shows
whether a particular pattern is pulling its weight or just adding noise.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict

from toxtree import assess, parse_structure


def main(path: str, scheme: str, per_alert: bool = False) -> int:
    rows = [r for r in csv.DictReader(open(path, newline="", encoding="utf-8"))
            if r.get("smiles")]

    tp = fp = tn = fn = 0
    unparsed = 0
    alert_counts: dict[str, Counter] = defaultdict(Counter)
    domain_hits: dict[str, Counter] = defaultdict(Counter)
    missed: list[dict] = []

    for r in rows:
        rec = parse_structure(r["smiles"], r.get("name") or None)
        if not rec.ok:
            unparsed += 1
            continue
        hits = getattr(assess(rec, {scheme}), scheme)
        flagged = bool(hits)
        active = r["active"] == "1"

        if active and flagged:
            tp += 1
        elif active and not flagged:
            fn += 1
            missed.append(r)
        elif not active and flagged:
            fp += 1
        else:
            tn += 1

        for h in hits:
            alert_counts[f"{h.code} {h.name}"][active] += 1
            if r.get("domain"):
                domain_hits[r["domain"]][h.code] += 1

    n = tp + fp + tn + fn
    if not n:
        print("Nothing scored.")
        return 1

    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    ppv = tp / (tp + fp) if tp + fp else 0.0

    print(f"{path.rsplit('/', 1)[-1]} -- {scheme} alerts")
    print(f"{n} compounds scored ({tp + fn} active, {tn + fp} inactive)"
          + (f", {unparsed} unparsable" if unparsed else "") + "\n")
    print("                 flagged   not flagged")
    print(f"  active        {tp:>8}{fn:>14}")
    print(f"  inactive      {fp:>8}{tn:>14}")
    print(f"\n  Sensitivity (actives flagged):   {100 * sens:5.1f}%   <- what an alert set is for")
    print(f"  Specificity (inactives cleared): {100 * spec:5.1f}%")
    print(f"  Balanced accuracy:               {100 * (sens + spec) / 2:5.1f}%")
    print(f"  Precision of a positive call:    {100 * ppv:5.1f}%")
    print(f"  Overall concordance:             {100 * (tp + tn) / n:5.1f}%")

    if per_alert and alert_counts:
        print(f"\n{'alert':44}{'active':>8}{'inactive':>10}{'precision':>11}")
        print("-" * 73)
        for label, c in sorted(alert_counts.items(),
                               key=lambda kv: -(kv[1][True] + kv[1][False])):
            a, i = c[True], c[False]
            prec = f"{100 * a / (a + i):.0f}%" if a + i else "-"
            print(f"  {label[:42]:42}{a:>8}{i:>10}{prec:>11}")

    if domain_hits:
        print("\nReference mechanistic domain vs the alert this implementation fires:")
        for domain, codes in sorted(domain_hits.items(),
                                    key=lambda kv: -sum(kv[1].values())):
            top = ", ".join(f"{c}x{n}" for c, n in codes.most_common(3))
            print(f"  {domain or '(blank)':10} {top}")

    if missed:
        print(f"\nActives with no alert ({len(missed)}), showing up to 15:")
        for m in missed[:15]:
            print(f"  {m.get('name', '')[:34]:36}{m['smiles'][:40]}")

    return 0


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__)
        raise SystemExit(2)
    scheme = "genotox"
    if "--scheme" in sys.argv:
        scheme = sys.argv[sys.argv.index("--scheme") + 1]
    raise SystemExit(main(args[0], scheme, per_alert="--per-alert" in sys.argv))
