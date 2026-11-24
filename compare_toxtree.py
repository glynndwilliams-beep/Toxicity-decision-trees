"""Three-way comparison against the Java Toxtree reference implementation.

For each structure in an SDF this reports three labels side by side:

    reference   the experimental or published class recorded in the file
    ours        what the rules in this repository say
    toxtree     what Toxtree 3.1.0 says, run headless on the same file

and sorts the disagreements into the four buckets that matter:

    both agree with the reference        -- nothing to learn
    both differ, identically             -- the *scheme's* limit, not a bug here
    Toxtree right, ours wrong            -- a real gap in this implementation
    ours right, Toxtree wrong            -- a defensible improvement

Structures are passed to both sides as the molblock drawn in the file, never
as regenerated SMILES. That matters: Toxtree answers two Cramer questions by
InChI lookup against its curated lists, and the same compound drawn
differently can resolve differently (ascorbic acid drawn by RDKit answers Q1
yes and lands in Class I; as drawn in foodmol.sdf it answers no and lands in
Class III). Re-deriving structures would measure the drawing, not the rules.

Usage
-----
    python compare_toxtree.py cramer  path/to/foodmol.sdf
    python compare_toxtree.py verhaar path/to/Verhaar2000.sdf
    python compare_toxtree.py genotox path/to/iss2_rid.sdf --limit 100

Toxtree is located from --jar/--java, else $TOXTREE_JAR/$TOXTREE_JAVA, else
$TOXTREE_HOME, else a Toxtree-*.jar under ~/toxtree-java. It is GPL-licensed
and not distributed here; see the README for how to obtain it.
"""
from __future__ import annotations

import argparse
import csv
import os
import pathlib
import re
import subprocess
import sys
import tempfile
from collections import Counter

from rdkit import Chem, RDLogger

from toxtree import assess, parse_structure
from toxtree.cramer import CLASS_I, CLASS_II, CLASS_III
from toxtree import verhaar as verhaar_mod

RDLogger.DisableLog("rdApp.*")

# Reference files carry mojibake in some compound names; a cp1252 console
# would otherwise abort the whole report on one unprintable character.
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

CRAMER_OF = {CLASS_I: 1, CLASS_II: 2, CLASS_III: 3}
VERHAAR_OF = {
    verhaar_mod.CLASS_1: 1, verhaar_mod.CLASS_2: 2, verhaar_mod.CLASS_3: 3,
    verhaar_mod.CLASS_4: 4, verhaar_mod.CLASS_5: 5,
}
ROMAN = {"I": 1, "II": 2, "III": 3}

ALERT = "alerting"
CLEAR = "no alert"


# --------------------------------------------------------------------------
# label readers: raw SDF text -> a comparable label, or None for "no answer"
# --------------------------------------------------------------------------

def _cramer_label(text: str):
    m = re.search(r"Class\s+(I{1,3})\b", text or "")
    if m:
        return ROMAN[m.group(1)]
    m = re.fullmatch(r"\s*([123])(?:\.0)?\s*", text or "")
    return int(m.group(1)) if m else None


def _verhaar_label(text: str):
    m = re.search(r"Class\s*(\d)", text or "")
    if m:
        return int(m.group(1))
    m = re.fullmatch(r"\s*([1-5])(?:\.0)?\s*", text or "")
    return int(m.group(1)) if m else None


def _yesno_negative_label(text: str):
    """Toxtree phrases these as 'Negative for X' = YES/NO."""
    t = (text or "").strip().upper()
    if t == "YES":
        return CLEAR
    if t == "NO":
        return ALERT
    return None


def _counter_label(text: str):
    try:
        return ALERT if float(text) > 0 else CLEAR
    except (TypeError, ValueError):
        return None


def _sicret_label(text: str):
    t = (text or "").strip()
    if not t or t.lower().startswith("unknown"):
        return None
    if t.lower().startswith("not"):
        return CLEAR
    if "irritat" in t.lower() or "corrosiv" in t.lower():
        return ALERT
    return None


def _eye_label(text: str):
    t = (text or "").strip()
    if not t or t.lower().startswith("unknown"):
        return None
    return CLEAR if t.lower().startswith("not") else ALERT


def _ames_label(text: str):
    """ISSCAN mutagenicity: 3 is active, 1 is inactive. 2 and ND are dropped."""
    t = (text or "").strip()
    return {"3": ALERT, "1": CLEAR}.get(t)


def _ec3_label(text: str):
    """LLNA: an EC3 value means a sensitiser, 'NC' means not classified."""
    t = (text or "").strip()
    try:
        float(t)
        return ALERT
    except ValueError:
        return CLEAR if t.upper() in {"NC", "NS", "NON"} else None


# --------------------------------------------------------------------------
# our side
# --------------------------------------------------------------------------

def _ours_cramer(a):
    return CRAMER_OF.get(a.cramer.cramer_class) if a.cramer else None


def _ours_verhaar(a):
    return VERHAAR_OF.get(a.verhaar.verhaar_class) if a.verhaar else None


def _ours_alerts(attr):
    def get(a):
        return ALERT if getattr(a, attr) else CLEAR
    return get


SCHEMES = {
    "cramer": dict(
        module="toxTree.tree.cramer.CramerRules",
        verdict="Cramer rules",
        path="toxTree.tree.cramer.CramerTreeResult",
        read=_cramer_label,
        ref_field="toxTree.tree.cramer.CramerRules",
        ref_read=_cramer_label,
        ours=_ours_cramer,
        our_schemes={"cramer"},
        labels=[1, 2, 3],
    ),
    "verhaar": dict(
        module="verhaar.VerhaarScheme",
        verdict="Verhaar scheme",
        path="Verhaar scheme#explanation",
        read=_verhaar_label,
        ref_field="Class",
        ref_read=_verhaar_label,
        ours=_ours_verhaar,
        our_schemes={"verhaar"},
        labels=[1, 2, 3, 4, 5],
        caveat=(
            "Toxtree 3.1.0 ships only verhaar.VerhaarScheme, the original\n"
            "  1992 scheme; the verhaar2 plugin that carries the Enoch (2008)\n"
            "  refinements is not in the distribution, and Verhaar2000.sdf is\n"
            "  verhaar2's own test fixture. This implementation follows the\n"
            "  refined scheme, so the two sides are NOT like for like: expect\n"
            "  Toxtree to answer Class 5 (cannot classify) far more often,\n"
            "  which is the difference the refinements were made to fix."),
    ),
    "genotox": dict(
        module="mutant.BB_CarcMutRules",
        verdict="Negative for genotoxic carcinogenicity",
        path=None,
        read=_yesno_negative_label,
        ref_field="SAL_iss2",
        ref_read=_ames_label,
        ours=_ours_alerts("genotox"),
        our_schemes={"genotox"},
        labels=[ALERT, CLEAR],
        cross_check="ALERTSCounter",
    ),
    "nongenotox": dict(
        module="mutant.BB_CarcMutRules",
        verdict="Negative for nongenotoxic carcinogenicity",
        path=None,
        read=_yesno_negative_label,
        ref_field=None,
        ref_read=None,
        ours=_ours_alerts("nongenotox"),
        our_schemes={"nongenotox"},
        labels=[ALERT, CLEAR],
    ),
    "sensitisation": dict(
        module="toxtree.plugins.skinsensitisation.SkinSensitisationPlugin",
        verdict="ALERTSCounter",
        path=None,
        read=_counter_label,
        ref_field="LLNA EC3 %",
        ref_read=_ec3_label,
        ours=_ours_alerts("sensitisation"),
        our_schemes={"sensitisation"},
        labels=[ALERT, CLEAR],
    ),
    "irritation": dict(
        module="sicret.SicretRules",
        verdict="Skin irritation / skin corrosion",
        path="Skin irritation / skin corrosion#explanation",
        read=_sicret_label,
        ref_field=None,
        ref_read=None,
        ours=_ours_alerts("irritation"),
        our_schemes={"irritation"},
        labels=[ALERT, CLEAR],
    ),
    "eye": dict(
        module="eye.EyeIrritationRules",
        verdict="Eye irritation and corrosion",
        path="Eye irritation and corrosion#explanation",
        read=_eye_label,
        ref_field=None,
        ref_read=None,
        ours=_ours_alerts("irritation"),
        our_schemes={"irritation"},
        labels=[ALERT, CLEAR],
    ),
}


# --------------------------------------------------------------------------
# SDF handling
# --------------------------------------------------------------------------

def sdf_blocks(path: pathlib.Path):
    """Yield (molblock, fields) for each record, keeping the drawing intact."""
    text = path.read_text(encoding="utf-8", errors="replace")
    for idx, chunk in enumerate(text.split("$$$$")):
        if not chunk.strip():
            continue
        # Every chunk but the first opens with the newline that terminated the
        # previous "$$$$". Drop exactly that one: a molblock's first line is
        # the title and is often blank, so stripping blank lines generally
        # would shift the whole record and break parsing.
        if idx and chunk.startswith("\n"):
            chunk = chunk[1:]
        lines = chunk.splitlines()
        end = next((i for i, l in enumerate(lines) if l.rstrip() == "M  END"),
                   None)
        mol_lines = lines[:end + 1] if end is not None else []
        fields, i = {}, (end + 1 if end is not None else 0)
        while i < len(lines):
            m = re.match(r"^>\s+<([^>]+)>", lines[i])
            if m:
                key, vals = m.group(1), []
                i += 1
                while i < len(lines) and lines[i].strip():
                    vals.append(lines[i].strip())
                    i += 1
                fields[key] = " ".join(vals)
            else:
                i += 1
        if mol_lines or fields:
            yield "\n".join(mol_lines), fields


def locate_toxtree(jar_arg, java_arg):
    jar = jar_arg or os.environ.get("TOXTREE_JAR")
    java = java_arg or os.environ.get("TOXTREE_JAVA") or "java"
    if not jar:
        roots = []
        if os.environ.get("TOXTREE_HOME"):
            roots.append(pathlib.Path(os.environ["TOXTREE_HOME"]))
        roots.append(pathlib.Path.home() / "toxtree-java")
        for root in roots:
            found = sorted(root.rglob("Toxtree-*.jar"))
            found = [f for f in found
                     if "sources" not in f.name and "javadoc" not in f.name]
            if found:
                jar = str(found[0])
                break
    if not jar:
        raise SystemExit(
            "Toxtree jar not found. Pass --jar, or set TOXTREE_JAR.\n"
            "See the README section 'Comparing against Toxtree' for how to "
            "obtain it.")
    if not java_arg and not os.environ.get("TOXTREE_JAVA"):
        local = sorted(pathlib.Path(jar).parents[2].rglob("bin/java.exe")) \
            if len(pathlib.Path(jar).parents) > 2 else []
        if local:
            java = str(local[0])
    return pathlib.Path(jar), java


def run_toxtree(java, jar: pathlib.Path, module: str,
                src: pathlib.Path, dst: pathlib.Path, verbose: bool) -> None:
    """Run Toxtree headless. cwd must be the jar's folder so `ext` resolves."""
    cmd = [str(java), "-jar", jar.name, "--noui",
           "-i", str(src), "-o", str(dst), "-m", module]
    if verbose:
        print("  $ " + " ".join(cmd))
    proc = subprocess.run(cmd, cwd=str(jar.parent), capture_output=True,
                          text=True, errors="replace")
    log = (proc.stdout or "") + (proc.stderr or "")
    done = re.search(r"Records processed (\d+)\s", log)
    if not dst.exists():
        tail = "\n".join(l for l in log.splitlines()
                         if "SEVERE" in l or "Exception" in l)[-800:]
        raise SystemExit(f"Toxtree produced no output for {module}.\n{tail}")
    if done:
        print(f"  Toxtree processed {done.group(1)} records")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scheme", choices=sorted(SCHEMES))
    p.add_argument("sdf")
    p.add_argument("--jar")
    p.add_argument("--java")
    p.add_argument("--reference-field",
                   help="SDF field holding the reference label "
                        "(default: the scheme's usual one)")
    p.add_argument("--limit", type=int, help="only the first N records")
    p.add_argument("--show", type=int, default=25,
                   help="disagreements to print (default 25)")
    p.add_argument("--out", help="write the full three-way join to this CSV")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args(argv)

    spec = SCHEMES[args.scheme]
    src = pathlib.Path(args.sdf)
    if not src.exists():
        raise SystemExit(f"no such file: {src}")
    jar, java = locate_toxtree(args.jar, args.java)
    print(f"Toxtree : {jar}")
    print(f"Java    : {java}")
    print(f"Module  : {spec['module']}")
    if spec.get("caveat"):
        print(f"\n  CAVEAT: {spec['caveat']}\n")

    records = list(sdf_blocks(src))
    if args.limit:
        records = records[:args.limit]
        trimmed = pathlib.Path(tempfile.mkdtemp()) / src.name
        trimmed.write_text(
            "".join(f"{mb}\n" + "".join(f">  <{k}>\n{v}\n\n"
                                        for k, v in f.items()) + "$$$$\n"
                    for mb, f in records),
            encoding="utf-8")
        src_for_run = trimmed
    else:
        src_for_run = src
    print(f"Input   : {src}  ({len(records)} records)")

    out_sdf = pathlib.Path(tempfile.mkdtemp()) / f"{args.scheme}_out.sdf"
    run_toxtree(java, jar, spec["module"], src_for_run.resolve(),
                out_sdf, args.verbose)
    produced = list(sdf_blocks(out_sdf))

    # Join. Same file in, same order out, so position is the natural key;
    # Toxtree silently drops records it cannot process, so verify.
    if len(produced) == len(records):
        pairs = list(zip(records, produced))
        join = "position"
    else:
        # Toxtree silently drops records it cannot process, so position no
        # longer lines up. It copies the input's own data fields through to
        # the output, so join on whichever identifier field is unique.
        candidates = ["ID", "Substance ID", "CAS", "CasRN", "NAME", "Name",
                      "Names", "ChemName_iss2", "SMILES", "cdk:Title"]
        best, join_field, best_index = -1, None, None
        for cand in candidates:
            ins = Counter(f[cand] for _, f in records if f.get(cand))
            outs = Counter(f[cand] for _, f in produced if f.get(cand))
            if not ins:
                continue
            # A reference file may repeat a substance, so tolerate a few
            # ambiguous values and simply leave those records unjoined
            # rather than discard an otherwise good identifier.
            usable = {v for v in ins if ins[v] == 1 and outs.get(v) == 1}
            if len(usable) < 0.9 * len(ins):
                continue
            index = {f[cand]: f for _, f in produced
                     if f.get(cand) in usable}
            matched = len(index)
            if matched > best:
                best, join_field, best_index = matched, cand, index
        if join_field is None or best < 0.9 * len(produced):
            raise SystemExit(
                f"Toxtree returned {len(produced)} of {len(records)} records "
                f"and no shared identifier field matches enough of them "
                f"(best: {join_field} with {best}). Tried: "
                + ", ".join(candidates))
        index = best_index
        pairs = [((mb, f), (mb, index.get(f.get(join_field))))
                 for mb, f in records]
        join = f"field <{join_field}>"
        print(f"  note: Toxtree returned {len(produced)} of {len(records)} "
              f"records; joining on {join}")

    ref_field = args.reference_field or spec["ref_field"]
    ref_read = spec["ref_read"] if not args.reference_field else (
        spec["ref_read"] or spec["read"])
    rows, errors = [], Counter()

    for (molblock, in_fields), (_, out_fields) in pairs:
        name = (in_fields.get("NAME") or in_fields.get("Name")
                or in_fields.get("Names") or in_fields.get("ChemName_iss2")
                or in_fields.get("cdk:Title") or "?")
        row = {"name": name, "reference": None, "ours": None,
               "toxtree": None, "toxtree_path": "", "note": ""}

        if ref_field and ref_read:
            row["reference"] = ref_read(in_fields.get(ref_field, ""))

        rec = parse_structure(molblock, name)
        if not rec.ok:
            row["note"] = f"unparsed here: {rec.error}"
            errors["unparsed by RDKit"] += 1
        else:
            row["ours"] = spec["ours"](assess(rec, spec["our_schemes"]))

        if out_fields is None:
            row["note"] = (row["note"] + "; " if row["note"] else "") + \
                "no Toxtree record"
            errors["dropped by Toxtree"] += 1
        else:
            if out_fields.get("Error when applying the decision tree",
                              "").strip().upper() == "YES":
                errors["Toxtree internal error"] += 1
                row["note"] = (row["note"] + "; " if row["note"] else "") + \
                    "Toxtree reported a tree error"
            row["toxtree"] = spec["read"](out_fields.get(spec["verdict"], ""))
            if spec.get("path"):
                row["toxtree_path"] = out_fields.get(spec["path"], "")
            cc = spec.get("cross_check")
            if cc and row["toxtree"] is not None:
                alt = _counter_label(out_fields.get(cc, ""))
                if alt is not None and alt != row["toxtree"]:
                    errors[f"{spec['verdict']!r} disagrees with {cc}"] += 1
        rows.append(row)

    # ---- report ----------------------------------------------------------
    print(f"\nJoined on {join}.")
    if errors:
        print("Excluded or flagged:")
        for k, v in errors.most_common():
            print(f"  {v:5d}  {k}")

    both = [r for r in rows if r["ours"] is not None and r["toxtree"] is not None]
    print(f"\n=== ours vs Toxtree ({len(both)} comparable) ===")
    agree = sum(1 for r in both if r["ours"] == r["toxtree"])
    pct = 100 * agree / len(both) if both else 0
    print(f"agree {agree}/{len(both)} = {pct:.1f}%")
    labels = spec["labels"]
    seen = sorted({r["ours"] for r in both} | {r["toxtree"] for r in both},
                  key=lambda x: (labels.index(x) if x in labels else 99, str(x)))
    if seen:
        head = "".join(f"{str(s)[:12]:>14s}" for s in seen)
        print(f"\n{'ours \\ toxtree':>20s}{head}")
        for a in seen:
            cells = "".join(
                f"{sum(1 for r in both if r['ours'] == a and r['toxtree'] == b):14d}"
                for b in seen)
            print(f"{str(a)[:20]:>20s}{cells}")

    graded = [r for r in both if r["reference"] is not None]
    if graded:
        print(f"\n=== against the reference ({len(graded)} graded) ===")
        buckets = Counter()
        for r in graded:
            ok_ours = r["ours"] == r["reference"]
            ok_tt = r["toxtree"] == r["reference"]
            buckets["both right" if ok_ours and ok_tt else
                    "ours right, Toxtree wrong" if ok_ours else
                    "Toxtree right, ours wrong" if ok_tt else
                    "both wrong, same answer" if r["ours"] == r["toxtree"] else
                    "both wrong, differently"] += 1
        for k in ("both right", "both wrong, same answer",
                  "Toxtree right, ours wrong", "ours right, Toxtree wrong",
                  "both wrong, differently"):
            if buckets[k]:
                print(f"  {buckets[k]:5d}  ({100*buckets[k]/len(graded):4.1f}%)  {k}")
        o = sum(1 for r in graded if r["ours"] == r["reference"])
        t = sum(1 for r in graded if r["toxtree"] == r["reference"])
        print(f"\n  ours    {o}/{len(graded)} = {100*o/len(graded):.1f}%")
        print(f"  Toxtree {t}/{len(graded)} = {100*t/len(graded):.1f}%")

        # For the alert rulebases a single concordance number hides the
        # screening trade-off, so report the same figures validate_alerts.py
        # does: these sets are built to miss as few actives as possible.
        if set(spec["labels"]) == {ALERT, CLEAR}:
            print(f"\n  {'':10s}{'sens':>8s}{'spec':>8s}{'bal acc':>9s}"
                  f"{'prec':>8s}")
            for who in ("ours", "toxtree"):
                tp = sum(1 for r in graded
                         if r["reference"] == ALERT and r[who] == ALERT)
                fn = sum(1 for r in graded
                         if r["reference"] == ALERT and r[who] == CLEAR)
                fp = sum(1 for r in graded
                         if r["reference"] == CLEAR and r[who] == ALERT)
                tn = sum(1 for r in graded
                         if r["reference"] == CLEAR and r[who] == CLEAR)
                sens = tp / (tp + fn) if tp + fn else 0
                spc = tn / (tn + fp) if tn + fp else 0
                prec = tp / (tp + fp) if tp + fp else 0
                print(f"  {who:10s}{100*sens:7.0f}%{100*spc:7.0f}%"
                      f"{100*(sens+spc)/2:8.0f}%{100*prec:7.0f}%")

        interesting = [r for r in graded
                       if (r["ours"] == r["reference"]) != (r["toxtree"] == r["reference"])]
        if interesting:
            print(f"\n  where exactly one is right ({len(interesting)}; "
                  f"showing {min(args.show, len(interesting))}):")
            for r in interesting[:args.show]:
                who = "ours" if r["ours"] == r["reference"] else "Toxtree"
                print(f"    {r['name'][:34]:34s} ref {str(r['reference']):9s} "
                      f"ours {str(r['ours']):9s} tt {str(r['toxtree']):9s} "
                      f"-> {who}")
                if r["toxtree_path"]:
                    print(f"      toxtree path: {r['toxtree_path'][:110]}")
    else:
        disagree = [r for r in both if r["ours"] != r["toxtree"]]
        if disagree:
            print(f"\nno reference labels; showing "
                  f"{min(args.show, len(disagree))} of {len(disagree)} "
                  f"disagreements:")
            for r in disagree[:args.show]:
                print(f"  {r['name'][:34]:34s} ours {str(r['ours']):9s} "
                      f"tt {str(r['toxtree']):9s}")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"\nfull join written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
