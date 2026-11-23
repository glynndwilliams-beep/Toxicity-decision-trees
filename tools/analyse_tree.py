"""Recover the Cramer tree topology from recorded decision paths.

Both reference sources record, for each compound, the sequence of questions
asked and the Y/N answer given -- caffeine, for instance, is
``1N,2N,3N,5N,6N,7Y,8N,10N,11Y,33N``. Roughly 190 such paths are enough to
reconstruct the shape of the original 33-question tree without reading any
implementation of it: which question follows which on each answer, and which
questions terminate in which class.

Usage:

    python tools/fetch_reference.py out/
    python tools/fetch_reference_lists.py toxtree/data/
    python tools/analyse_tree.py out/ toxtree/data/          # topology
    python tools/analyse_tree.py out/ toxtree/data/ --q 22   # one question

Passing ``--q N`` lists the compounds that answered question N yes and no, so
the rule behind it can be inferred from the structures rather than guessed.
That is how Q22 was identified as the food-component question (answering yes
gives Class II, not Class I) and Q11 as the acylated-ring-heteroatom question
that gates access to it.
"""
from __future__ import annotations

import csv
import pathlib
import re
import sys
from collections import Counter, defaultdict

TOKEN = re.compile(r"(\d+)([YN])")
ENTRY = re.compile(r'\{\s*("(?:[^"\\]|\\.)*"(?:\s*,\s*"(?:[^"\\]|\\.)*")*)\s*\}')
STR = re.compile(r'"((?:[^"\\]|\\.)*)"')


def steps(path: str) -> list[str]:
    """Flatten a recorded path to ordered <number><Y|N> tokens."""
    return [n + yn for n, yn in TOKEN.findall(path)]


def load_paths(*dirs: str) -> list[dict]:
    rows: list[dict] = []
    for d in (pathlib.Path(x) for x in dirs):
        for sdf in d.glob("*.sdf"):
            rows.extend(_from_sdf(sdf.read_text(encoding="utf-8", errors="replace")))
        for java in d.glob("*.java"):
            rows.extend(_from_java(java.read_text(encoding="utf-8", errors="replace")))
        for csvf in d.glob("*_paths.csv"):
            rows.extend(csv.DictReader(csvf.open(encoding="utf-8")))
    return [r for r in rows if r.get("path") and r.get("cramer")]


def _from_sdf(text: str) -> list[dict]:
    out = []
    for block in text.split("$$$$"):
        if not block.strip():
            continue
        fields, key = {}, None
        for line in block.splitlines():
            if line.startswith("> "):
                key = line.split("<", 1)[-1].split(">", 1)[0] if "<" in line else None
            elif key is not None:
                v = line.strip()
                fields.setdefault(key, v) if v else None
                if not v:
                    key = None
        verdict = fields.get("toxTree.tree.cramer.CramerRules", "")
        num = {"I": 1, "II": 2, "III": 3}.get(
            verdict.split("Class ")[-1].rstrip(")").strip())
        path = fields.get("toxTree.tree.cramer.CramerTreeResult", "")
        if num and path:
            out.append({"smiles": fields.get("SMILES", ""), "cramer": num,
                        "path": path, "name": fields.get("NAME", "")})
    return out


def _from_java(text: str) -> list[dict]:
    out = []
    for m in ENTRY.finditer(text):
        parts = [p.replace('\\\\', '\\') for p in STR.findall(m.group(1))]
        if len(parts) >= 4 and parts[1] in {"1", "2", "3"} and TOKEN.match(parts[3]):
            out.append({"smiles": parts[0], "cramer": int(parts[1]),
                        "path": parts[3], "name": parts[2]})
    return out


def topology(rows: list[dict]) -> None:
    succ: dict[str, Counter] = defaultdict(Counter)
    terminal: dict[str, Counter] = defaultdict(Counter)
    seen: Counter = Counter()

    for r in rows:
        st = steps(r["path"])
        for i, s in enumerate(st):
            seen[s[:-1]] += 1
            if i + 1 < len(st):
                succ[s][st[i + 1][:-1]] += 1
            else:
                terminal[s][int(r["cramer"])] += 1

    print(f"{len(rows)} decision paths\n")
    print(f"{'question':10}{'n':<6}{'on Y ->':<24}on N ->")
    print("-" * 74)
    for q in sorted(seen, key=int):
        def fmt(branch: str) -> str:
            nx, tm = succ.get(q + branch, Counter()), terminal.get(q + branch, Counter())
            bits = []
            if nx:
                bits.append("Q" + ",Q".join(k for k, _ in nx.most_common(3)))
            if tm:
                bits.append(" ".join(f"[CLASS {c}]x{n}" for c, n in sorted(tm.items())))
            return " ".join(bits) or "-"
        print(f"  Q{q:<8}{seen[q]:<6}{fmt('Y'):<24}{fmt('N')}")


def inspect(rows: list[dict], target: str, limit: int = 12) -> None:
    yes, no = [], []
    for r in rows:
        for s in steps(r["path"]):
            if s[:-1] == target:
                (yes if s.endswith("Y") else no).append(r)
                break
    print(f"Q{target}: {len(yes)} answered Y, {len(no)} answered N\n")
    for label, group in (("YES", yes), ("NO", no)):
        classes = dict(sorted(Counter(int(r["cramer"]) for r in group).items()))
        print(f"--- {label} ({len(group)})  resulting classes: {classes}")
        for r in group[:limit]:
            print(f"    c{r['cramer']}  {r['name'][:26]:28s} {r['smiles'][:44]}")
        if len(group) > limit:
            print(f"    ... and {len(group) - limit} more")
        print()


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        print(__doc__)
        raise SystemExit(2)
    data = load_paths(*args)
    if not data:
        print("No decision paths found. Run the two fetch scripts first.")
        raise SystemExit(1)
    if "--q" in sys.argv:
        inspect(data, sys.argv[sys.argv.index("--q") + 1])
    else:
        topology(data)
