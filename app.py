"""Streamlit front end for the pure-Python toxicity decision-tree engine."""
from __future__ import annotations

import io

import pandas as pd
import streamlit as st

from toxtree import (SCHEME_LABELS, assess, parse_many, parse_structure,
                     readacross, reflists)
from toxtree.readacross import SOURCES
from toxtree.alerts import (
    ALL_ALERTS,
    COMPUTED,
    GENOTOX,
    IRRITATION,
    NONGENOTOX,
    SENSITISATION,
)
from toxtree.chem import to_svg
from toxtree.cramer import CLASS_I, CLASS_II, CLASS_III
from toxtree.examples import EXAMPLES, example_text
from toxtree.verhaar import lc50_mg_per_l

st.set_page_config(page_title="Toxicity Decision Trees", page_icon="⚗️",
                   layout="wide", initial_sidebar_state="expanded")

CSS = """
<style>
  .block-container {padding-top: 2.2rem; max-width: 1400px;}
  .verdict {border-radius: 10px; padding: 0.85rem 1rem; margin-bottom: 0.6rem;
            border: 1px solid rgba(128,128,128,0.35);}
  .verdict h4 {margin: 0 0 0.3rem 0; font-size: 0.78rem; letter-spacing: 0.08em;
               text-transform: uppercase; opacity: 0.7; font-weight: 600;}
  .verdict .val {font-size: 1.05rem; font-weight: 600; line-height: 1.3;}
  .v-low  {background: rgba(34,160,90,0.12);  border-color: rgba(34,160,90,0.5);}
  .v-mid  {background: rgba(220,160,20,0.12); border-color: rgba(220,160,20,0.5);}
  .v-high {background: rgba(220,60,60,0.12);  border-color: rgba(220,60,60,0.5);}
  .v-none {background: rgba(128,128,128,0.08);}
  .trace {font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          font-size: 0.82rem; line-height: 1.75;}
  .pill {display: inline-block; padding: 0.12rem 0.55rem; border-radius: 999px;
         font-size: 0.72rem; background: rgba(128,128,128,0.18); margin-right: 0.3rem;}
  .molbox svg {max-width: 100%; height: auto;}
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)

TONE = {
    CLASS_I: "v-low", CLASS_II: "v-mid", CLASS_III: "v-high",
    "No alert": "v-low",
}


def verdict_card(title: str, value: str, tone: str = "v-none", sub: str = "") -> str:
    extra = f'<div style="opacity:0.7;font-size:0.8rem;margin-top:0.25rem">{sub}</div>' if sub else ""
    return (f'<div class="verdict {tone}"><h4>{title}</h4>'
            f'<div class="val">{value}</div>{extra}</div>')


def tone_for(value: str) -> str:
    if value in TONE:
        return TONE[value]
    # The reactive-group flags deliberately read amber, not red: they say a
    # reactive fragment is present, not that the substance is an irritant.
    if value.startswith("None flagged"):
        return "v-low"
    if value.startswith("Flagged"):
        return "v-mid"
    if value.startswith("No alert"):
        return "v-low"
    if value.startswith("Alerting"):
        return "v-high"
    if "Class 1" in value or "Class I (" in value:
        return "v-low"
    if "Class 2" in value or "Class II " in value:
        return "v-mid"
    if "Class 3" in value or "Class 4" in value or "Class III" in value:
        return "v-high"
    return "v-none"


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.title("⚗️ Decision trees")
    st.caption("Structure-based toxicity classification in pure Python "
               "(RDKit + Streamlit).")

    st.subheader("Rulebases")
    chosen: set[str] = set()
    defaults = {"cramer", "verhaar", "genotox", "sensitisation", "irritation",
                "nongenotox", "ttc"}
    for key, label in SCHEME_LABELS.items():
        if st.checkbox(label, value=key in defaults, key=f"scheme_{key}"):
            chosen.add(key)
    if not chosen:
        st.warning("Select at least one rulebase.")

    st.divider()
    st.subheader("Library")
    st.metric("SMARTS alerts", len(ALL_ALERTS) + len(COMPUTED))
    st.caption(
        f"{len(GENOTOX)} genotoxicity · {len(NONGENOTOX)} non-genotoxicity · "
        f"{len(SENSITISATION)} sensitisation · {len(IRRITATION)} reactive-group"
    )

    rl = reflists.status()
    st.caption(f"✓ {len(reflists.PRIMARY_METABOLITES)} primary metabolites "
               "built in (Cramer Q1).")
    ra_status = readacross.status()
    ra_total = sum(ra_status.values())
    if ra_total:
        st.caption(f"✓ Read-across data loaded — {ra_total} reference "
                   "compounds across "
                   f"{sum(1 for v in ra_status.values() if v)} endpoints.")
    else:
        st.caption("Read-across data not installed, so no evidential-support "
                   "scores. Fetch it with the scripts in `tools/`.")

    if rl["food_loaded"]:
        st.caption(f"✓ Food-component list loaded — {rl['food_count']} "
                   "structures, used at the Q22 question.")
    else:
        st.caption("Food-component list not installed, so the Q22 question "
                   "cannot be answered and Class II will be under-reported. "
                   "Fetch it with "
                   "`python tools/fetch_reference_lists.py toxtree/data/`.")

st.title("Toxicity decision trees")
st.caption("Cramer classes · Verhaar aquatic modes of action · genotoxicity "
           "and sensitisation alerts · reactive-group flags · Kroes TTC")

tab_single, tab_batch, tab_rules, tab_about = st.tabs(
    ["Single structure", "Batch screening", "Rulebase reference", "Method notes"]
)

# ---------------------------------------------------------------------------
# Single structure
# ---------------------------------------------------------------------------
with tab_single:
    c1, c2 = st.columns([3, 2], gap="large")
    with c1:
        names = [f"{n} — {why}" for n, _, why in EXAMPLES]
        pick = st.selectbox("Load an example", ["(type your own)"] + names, index=1)
        if pick != "(type your own)":
            default_smiles = EXAMPLES[names.index(pick)][1]
            default_name = EXAMPLES[names.index(pick)][0]
        else:
            default_smiles, default_name = "", ""
        query = st.text_input("SMILES or InChI", value=default_smiles,
                              placeholder="e.g. COc1cc(C=O)ccc1O",
                              key=f"q_{pick}")
    with c2:
        st.write("")
        st.write("")
        exposure = st.number_input(
            "Optional: human exposure for the TTC comparison (µg/person/day)",
            min_value=0.0, value=0.0, step=1.0,
            help="Leave at 0 to skip the margin-of-exposure comparison.")

    if query.strip():
        rec = parse_structure(query, default_name or None)
        if not rec.ok:
            st.error(f"Could not read that structure — {rec.error}")
        else:
            a = assess(rec, chosen)
            p = rec.props

            left, right = st.columns([2, 3], gap="large")
            with left:
                hits = a.all_hits()
                options = ["(no highlight)"] + [f"{h.code} · {h.name}" for h in hits
                                                if h.atoms]
                sel = st.selectbox("Highlight an alert", options,
                                   label_visibility="collapsed") if len(options) > 1 else None
                highlight: list[int] = []
                if sel and sel != "(no highlight)":
                    code = sel.split(" · ")[0]
                    highlight = next((h.atoms for h in hits if h.code == code), [])
                st.markdown(f'<div class="molbox">{to_svg(rec.mol, 440, 320, highlight)}</div>',
                            unsafe_allow_html=True)
                st.code(rec.smiles, language=None)

            with right:
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("MW", f"{p['mw']:.1f}")
                m2.metric("logP", f"{p['logp']:.2f}")
                m3.metric("TPSA", f"{p['tpsa']:.0f}")
                m4.metric("Rings", p["rings"])
                st.caption(f"Formula {p['formula']} · {p['heavy_atoms']} heavy atoms · "
                           f"elements {', '.join(p['elements'])}")
                if p.get("counterions"):
                    st.caption("⚗️ Counterions removed before classification: "
                               + ", ".join(p["counterions"])
                               + " — the organic parent is assessed.")

                cols = st.columns(2)
                i = 0
                if a.cramer:
                    cols[i % 2].markdown(
                        verdict_card("Cramer class", a.cramer.cramer_class,
                                     tone_for(a.cramer.cramer_class),
                                     f"TTC basis {a.cramer.ttc_ug_day:g} µg/person/day"
                                     + (f"<br>Listed food component: "
                                        f"{a.cramer.food_component.title()}"
                                        if a.cramer.food_component else "")),
                        unsafe_allow_html=True)
                    i += 1
                if a.verhaar:
                    sub = ""
                    if a.verhaar.log_lc50_mol_l is not None and p.get("mw"):
                        sub = (f"Predicted 96 h fish LC50 ≈ "
                               f"{lc50_mg_per_l(a.verhaar.log_lc50_mol_l, p['mw']):g} mg/L")
                    cols[i % 2].markdown(
                        verdict_card("Verhaar mode of action", a.verhaar.verhaar_class,
                                     tone_for(a.verhaar.verhaar_class), sub),
                        unsafe_allow_html=True)
                    i += 1
                if "genotox" in chosen:
                    cols[i % 2].markdown(
                        verdict_card("Genotoxic carcinogenicity", a.genotox_verdict,
                                     tone_for(a.genotox_verdict),
                                     ", ".join(h.code for h in a.genotox)),
                        unsafe_allow_html=True)
                    i += 1
                if "nongenotox" in chosen:
                    cols[i % 2].markdown(
                        verdict_card("Non-genotoxic carcinogenicity", a.nongenotox_verdict,
                                     tone_for(a.nongenotox_verdict),
                                     ", ".join(h.code for h in a.nongenotox)),
                        unsafe_allow_html=True)
                    i += 1
                if "sensitisation" in chosen:
                    cols[i % 2].markdown(
                        verdict_card("Skin sensitisation", a.sensitisation_verdict,
                                     tone_for(a.sensitisation_verdict),
                                     ", ".join(h.code for h in a.sensitisation)),
                        unsafe_allow_html=True)
                    i += 1
                if "irritation" in chosen:
                    cols[i % 2].markdown(
                        verdict_card("Reactive group flags", a.irritation_verdict,
                                     tone_for(a.irritation_verdict),
                                     ", ".join(h.code for h in a.irritation)),
                        unsafe_allow_html=True)
                    i += 1

            st.divider()

            if a.ttc_value_ug_day is not None or a.ttc_trace:
                st.subheader("Threshold of toxicological concern")
                tc1, tc2 = st.columns([2, 3], gap="large")
                with tc1:
                    if a.ttc_value_ug_day is None:
                        st.warning(f"**{a.ttc_class}** — the TTC approach does not "
                                   "apply to this structure.")
                    else:
                        st.metric(f"TTC — {a.ttc_class}",
                                  f"{a.ttc_value_ug_day:g} µg/person/day")
                        if exposure > 0:
                            ratio = a.ttc_value_ug_day / exposure
                            if ratio >= 1:
                                st.success(f"Exposure of {exposure:g} µg/day is "
                                           f"{ratio:.1f}× below the threshold.")
                            else:
                                st.error(f"Exposure of {exposure:g} µg/day exceeds "
                                         f"the threshold by {1 / ratio:.1f}×.")
                with tc2:
                    st.markdown('<div class="trace">' +
                                "<br>".join(line.replace("  ", "&nbsp;&nbsp;&nbsp;")
                                            for line in a.ttc_trace) +
                                "</div>", unsafe_allow_html=True)

            if a.cramer:
                with st.expander("Cramer decision path", expanded=False):
                    st.markdown('<div class="trace">' +
                                "<br>".join(a.cramer.trace) + "</div>",
                                unsafe_allow_html=True)
                    if a.cramer.out_of_domain:
                        st.warning("This structure is outside the applicability "
                                   "domain of the Cramer tree, which was derived "
                                   "for organic substances without metals.")

            if a.verhaar:
                with st.expander("Verhaar decision path", expanded=False):
                    st.markdown('<div class="trace">' +
                                "<br>".join(a.verhaar.trace) + "</div>",
                                unsafe_allow_html=True)
                    if a.verhaar.qsar:
                        st.caption(f"Baseline QSAR applied: {a.verhaar.qsar}")

            # --- evidential support -------------------------------------
            predicted_for = {}
            if a.cramer:
                predicted_for["cramer"] = {
                    CLASS_I: "Class I", CLASS_II: "Class II",
                    CLASS_III: "Class III"}[a.cramer.cramer_class]
            if a.verhaar:
                predicted_for["verhaar"] = (
                    f"Class {a.verhaar.verhaar_class.split()[1]}")
            if "genotox" in chosen:
                predicted_for["genotox"] = ("mutagenic" if a.genotox
                                            else "not mutagenic")
            if "sensitisation" in chosen:
                predicted_for["sensitisation"] = ("sensitiser"
                                                  if a.sensitisation
                                                  else "non-sensitiser")

            usable = [e for e in readacross.available_endpoints()
                      if e in predicted_for]
            if usable:
                st.subheader("Evidential support")
                st.caption(
                    "How much checked evidence stands behind each verdict: the "
                    "nearest compounds in the reference sets used for "
                    "validation, and whether their **measured** outcomes agree. "
                    "These rules are not fitted to these compounds, so this is "
                    "not an applicability domain in the QSAR sense — a low "
                    "score means nothing in the reference data would have "
                    "caught this answer if it were wrong."
                )
                tabs = st.tabs([SCHEME_LABELS[e].split(" (")[0] for e in usable])
                for tab, endpoint in zip(tabs, usable):
                    with tab:
                        sup = readacross.support(rec.mol, endpoint,
                                                 predicted_for[endpoint])
                        sc1, sc2 = st.columns([2, 3], gap="large")
                        with sc1:
                            tone = ("v-low" if sup.index >= 0.70
                                    else "v-mid" if sup.index >= 0.45
                                    else "v-high")
                            st.markdown(
                                verdict_card(
                                    "Support index", f"{sup.index:.2f}", tone,
                                    f"{sup.band} · {sup.n_reference} reference "
                                    "compounds"),
                                unsafe_allow_html=True)
                            st.caption(
                                f"nearest {sup.similarity:.2f} · concordance "
                                + (f"{sup.concordance:.2f}"
                                   if sup.concordance is not None else "n/a")
                                + f" · new fragments {100 * sup.novelty:.0f}%"
                                + (" · in descriptor range" if sup.in_range
                                   else " · **outside** descriptor range"))
                            for note in sup.notes:
                                st.warning(note)
                        with sc2:
                            st.caption(f"Nearest reference compounds "
                                       f"({SOURCES[endpoint]['origin']})")
                            st.dataframe(pd.DataFrame([{
                                "Similarity": n.similarity,
                                "Compound": n.name,
                                "Measured": n.outcome,
                                "Agrees": ("yes" if n.outcome.lower()
                                           == predicted_for[endpoint].lower()
                                           else "no"),
                            } for n in sup.neighbours]),
                                hide_index=True, width='stretch')
                st.divider()

            hits = a.all_hits()
            if hits:
                with st.expander(f"Structural alerts that fired ({len(hits)})",
                                 expanded=True):
                    st.dataframe(
                        pd.DataFrame([{
                            "Code": h.code,
                            "Alert": h.name,
                            "Rulebase": SCHEME_LABELS.get(h.scheme, h.scheme),
                            "Mechanistic rationale": h.rationale,
                            "Matched atoms": len(h.atoms),
                        } for h in hits]),
                        hide_index=True, width='stretch')
            elif chosen - {"cramer", "verhaar", "ttc"}:
                st.info("No structural alerts fired for the selected rulebases.")

# ---------------------------------------------------------------------------
# Batch screening
# ---------------------------------------------------------------------------
with tab_batch:
    st.write("Paste one structure per line (`SMILES` optionally followed by a tab "
             "and a name), or upload a `.smi`, `.txt` or `.csv` file.")
    bc1, bc2 = st.columns([3, 2], gap="large")
    # bc2 is written first so the example loader can seed the text area before
    # that widget is instantiated.
    with bc2:
        upload = st.file_uploader("…or upload a file", type=["smi", "txt", "csv"])
        smiles_col = None
        if upload is not None and upload.name.lower().endswith(".csv"):
            raw = upload.getvalue().decode("utf-8", errors="replace")
            df_in = pd.read_csv(io.StringIO(raw))
            smiles_col = st.selectbox("SMILES column", list(df_in.columns))
        if st.button(f"Load the {len(EXAMPLES)} example structures",
                     width='stretch'):
            st.session_state["batch_text"] = example_text()
    with bc1:
        text = st.text_area("Structures", height=220, key="batch_text",
                            placeholder="CCO\tEthanol\nOc1ccccc1\tPhenol")

    lines: list[str] = []
    if text.strip():
        lines = text.splitlines()
    elif upload is not None:
        raw = upload.getvalue().decode("utf-8", errors="replace")
        if upload.name.lower().endswith(".csv") and smiles_col:
            df_in = pd.read_csv(io.StringIO(raw))
            lines = [str(v) for v in df_in[smiles_col].tolist()]
        else:
            lines = raw.splitlines()

    if st.button("Run screening", type="primary", disabled=not lines):
        records = parse_many(lines)
        failed = [r for r in records if not r.ok]
        results = [assess(r, chosen) for r in records if r.ok]
        rows = [a.to_row() for a in results]

        if not rows:
            st.error("Nothing could be parsed.")
        else:
            df = pd.DataFrame(rows)
            st.success(f"Screened {len(rows)} structure(s)."
                       + (f" {len(failed)} could not be parsed." if failed else ""))

            k = st.columns(4)
            if "cramer" in chosen:
                high = int((df["Cramer class"] == CLASS_III).sum())
                k[0].metric("Cramer Class III", high)
            if "genotox" in chosen:
                k[1].metric("Genotoxicity alerts",
                            int((df["Genotox alerts"] != "").sum()))
            if "sensitisation" in chosen:
                k[2].metric("Sensitisation alerts",
                            int((df["Sensitisation alerts"] != "").sum()))
            if "irritation" in chosen:
                k[3].metric("Reactive group flags",
                            int((df["Reactive group flags"] != "").sum()))

            st.dataframe(df, hide_index=True, width='stretch', height=460)

            d1, d2 = st.columns(2, gap="large")
            if "cramer" in chosen:
                with d1:
                    st.caption("Cramer class distribution")
                    st.bar_chart(df["Cramer class"].value_counts())
            if "verhaar" in chosen:
                with d2:
                    st.caption("Verhaar class distribution")
                    st.bar_chart(df["Verhaar class"].value_counts())

            st.download_button("Download results as CSV",
                               df.to_csv(index=False).encode("utf-8"),
                               file_name="toxicity_screening.csv",
                               mime="text/csv")
            if failed:
                with st.expander(f"{len(failed)} structure(s) could not be parsed"):
                    st.dataframe(pd.DataFrame(
                        [{"Input": r.raw_input, "Problem": r.error} for r in failed]),
                        hide_index=True, width='stretch')

# ---------------------------------------------------------------------------
# Rulebase reference
# ---------------------------------------------------------------------------
with tab_rules:
    st.write("Every alert in the library, with the SMARTS pattern it uses and the "
             "mechanism it stands for. You can test a pattern against any structure "
             "in the **Single structure** tab by highlighting it.")
    scheme_names = sorted({a.scheme for a in ALL_ALERTS})
    pick_schemes = st.multiselect("Filter by rulebase", scheme_names,
                                  default=scheme_names,
                                  format_func=lambda s: SCHEME_LABELS.get(s, s))
    search = st.text_input("Search alert names and rationales", "")

    rows = []
    for al in ALL_ALERTS:
        if al.scheme not in pick_schemes:
            continue
        blob = f"{al.code} {al.name} {al.rationale}".lower()
        if search and search.lower() not in blob:
            continue
        rows.append({"Code": al.code, "Alert": al.name,
                     "Rulebase": SCHEME_LABELS.get(al.scheme, al.scheme),
                     "SMARTS": al.smarts, "Rationale": al.rationale})
    for comp in COMPUTED:
        if comp.scheme in pick_schemes:
            rows.append({"Code": comp.code, "Alert": comp.name,
                         "Rulebase": SCHEME_LABELS.get(comp.scheme, comp.scheme),
                         "SMARTS": "(computed from ring topology)",
                         "Rationale": comp.rationale})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width='stretch',
                 height=620)

# ---------------------------------------------------------------------------
# Method notes
# ---------------------------------------------------------------------------
with tab_about:
    st.markdown(
        """
### What this is

A pure-Python implementation of the classical structure–activity decision trees
used in regulatory toxicology screening. Everything runs locally: RDKit does the
cheminformatics, Streamlit draws the interface, and there is no Java, no server
call and no external service.

The rules are written from the published literature, not transcribed from any
existing program:

| Rulebase | Source |
| --- | --- |
| Cramer classes | Cramer, Ford & Hall (1978), *Fd Cosmet. Toxicol.* 16, 255 |
| Verhaar scheme | Verhaar *et al.* (1992), *Chemosphere* 25, 471; refinements from Enoch *et al.* (2008) |
| Genotoxicity alerts | Benigni & Bossa structural-alert set for mutagenicity and carcinogenicity |
| Sensitisation alerts | Mechanistic protein-binding domains (Michael addition, Schiff base, acyl transfer, SN2, SNAr) |
| Reactive group flags | Protein-reactive and corrosive functional groups. **Not** a prediction of irritancy — see below |
| TTC tiers | Kroes *et al.* (2004), *Food Chem. Toxicol.* 42, 65 |

### How to read a result

Each scheme returns both a verdict **and** the path that produced it. The
decision path is the useful part: it tells you which question decided the class,
so you can judge whether the rule really applies to your structure. Expand the
*Cramer decision path* or *Verhaar decision path* panels to see it.

Structural alerts are **hypothesis generators**. An alert that fires says the
molecule contains a fragment associated with a mechanism — not that the molecule
is toxic. Absence of an alert is weaker still: it means none of the fragments in
this library matched.

### Validation

The Cramer implementation has been checked against the 83 compounds of
Cramer, Ford & Hall (1978) Appendix 1 — the classifications of the authors who
devised the scheme. Two of those have malformed structures in the published
fixture and cannot be read, leaving 81 scored:

| | → Class I | → Class II | → Class III | n |
| --- | --- | --- | --- | --- |
| **Reference Class I** | **31** | 0 | 1 | 32 |
| **Reference Class II** | 4 | **3** | 0 | 7 |
| **Reference Class III** | 1 | 1 | **36** | 38 |

**Overall concordance 70/77 = 90.9%** — Class I 97%, Class II 43%, Class III 95%.
Six compounds are placed in a *lower* class than the reference (the unsafe
direction) and one in a higher one.

A second set of 108 food-related structures carries the class each was assigned
by the reference implementation. That is a weaker oracle — it measures agreement
with a program rather than with the scheme's authors — but it contains 30
Class II compounds against Appendix 1's seven:

| | → Class I | → Class II | → Class III | n |
| --- | --- | --- | --- | --- |
| **Reference Class I** | **66** | 1 | 0 | 67 |
| **Reference Class II** | 2 | **27** | 1 | 30 |
| **Reference Class III** | 0 | 0 | **10** | 10 |

**102/106 = 96.2%** — Class I 98%, Class II 90%, Class III 100%.

Run both yourself: `python tools/fetch_reference.py out/` and
`python tools/fetch_reference_lists.py toxtree/data/`, then
`python validate_cramer.py out/cramer_reference.csv`.

### The Verhaar scheme

The aquatic scheme has its own reference set of 121 compounds, and checking it
moved it from **73.5% to 87.2%** — Class 1 100%, Class 2 82%, Class 3 92%,
Class 4 53%.

Two of the fixes were bugs in rules written here: an alert labelled "two or
more nitro groups" whose pattern matched a single one, so every nitrophenol
read as reactive rather than a polar narcotic; and a trifluoromethyl group
being matched as an activated benzylic halide. Three more were asserted
without evidence and contradicted by the data — triazines are Class 2,
dinitrophenols and thiocarbamates Class 3, and carbamates Class 3 rather than
specifically acting.

Class 4 is the weak spot and cannot be improved from this data. Six of the
eight misses are organochlorine insecticides, and no halogen-count rule
separates them: pentachlorophenol is Class 3 while tetrachlorophenol is
Class 4, and heptachlor is Class 3 while endrin is Class 4. Such a rule would
break seven compounds to fix six. The organophosphates are all correct.

### Evidential support, and what it is not

The **Evidential support** panel finds the query's nearest neighbours in
the reference sets used for validation and reports whether their measured
outcomes agree with the verdict. The index weights concordance 0.40,
similarity 0.30, descriptor range 0.15 and fragment novelty 0.15.

The idea is borrowed from VEGA's applicability domain index, but the
meaning differs in a way worth stating. VEGA scores statistical models,
where the question is whether a prediction interpolates or extrapolates.
These are deterministic rules, not fitted to these compounds, and they
answer identically for a molecule in the reference set and one unlike
anything ever tested. So this is not an applicability domain in the QSAR
sense: it measures how much checked evidence stands behind an answer. A
low score does not mean the rule misfired — it means nothing in the
reference data would have caught it if it had.

### Why this is "reactive group flags" and not "irritation"

This panel used to be labelled skin and eye irritation. It is not that,
and it has been renamed rather than left to imply more than it does.

What it computes is a SMARTS list of protein-reactive and corrosive
functional groups. The published approach to local irritation is a
different kind of thing: the BfR decision support systems work largely
from **physicochemical cut-offs** — logP, molecular weight, melting point,
aqueous solubility and pKa — with structural exclusions layered on top.

Comparing the two settles it. Run against the same compounds, this list
and Toxtree's SICRET skin rulebase agree on 50% of the 129 cases SICRET
will answer for; for eye irritation the agreement is 60% of 37. That is
about what two unrelated classifiers would manage, and it says the two are
measuring different things.

**Implementing the BfR rules properly is not available here**, which is
the other half of the reason for renaming rather than rewriting. Toxtree
itself returned "Unknown" for 80 of 209 compounds on SICRET and 172 of 209
on the eye rulebase, because the measured properties its rules need were
not in the file. Melting point, aqueous solubility and pKa cannot be
computed reliably from a structure, and this app has no source for them.
A partial implementation would be the same overstatement in new clothes.

So: read a flag here as "this molecule contains a reactive fragment",
which is a useful thing to know and is why the panel is kept. Do not read
it as a prediction of irritancy or corrosivity.

What could be checked without labels was. Across the 941 unique compounds
from the other reference sets, at least one flag fires on 29% of them and
no single alert exceeds an 11% rate, so the set discriminates rather than
flagging everything. That found no pathology, which is a negative result
rather than a validation. It did expose one defect that needs no labels:
the amino acids were firing both the acid and the base alert, though an
alpha-amino acid is zwitterionic at skin pH. Both now exclude that case.

### The structural alerts

Alert sets are screening filters, so the number that matters is
sensitivity. Against a 209-compound local lymph node assay set the
sensitisation alerts went from **52% to 92% sensitivity** (balanced
accuracy 60% to 70%, precision 88%); against the Ames set the
genotoxicity alerts reach **90% sensitivity**, 75% balanced accuracy.
Both sets have been compared against the Java Toxtree implementation
compound by compound; see the project README.

The sensitisation gaps grouped by mechanism: azlactones and lactones were
missing from acyl transfer, the SN2 pattern accepted bromide and iodide
but not chloride, the dicarbonyls were absent from Schiff base formation,
and the eugenol family of quinone-methide precursors was unreachable from
the catechol pattern alone.

The genotoxicity run first scored far **worse than chance**, which is the
signature of inverted labels rather than a bad model: the ISSCAN field
encodes 3 as active and 1 as inactive. Chemistry gave it away — 42
N-nitroso compounds had come out non-mutagenic. Corrected, the per-alert
precisions read as they should: N-nitroso 93%, nitroaromatic 86%,
polycyclic aromatic 85%.

### Recovering questions from decision paths

Both reference sources record the path each compound took through the tree, not
just its class — caffeine is `1N,2N,3N,5N,6N,7Y,8N,10N,11Y,33N`. Around 190 such
paths are enough to reconstruct the tree's shape without reading any
implementation of it, using `tools/analyse_tree.py`. Three questions were
identified this way and implemented:

- **Q22, "is it a common component of food?"** — answering **yes gives Class
  II**, not Class I. All 28 observed yes-answers are Class II.
- **Q11, "does a ring heteroatom bear an acyl group?"** — a cyclic amide, urea,
  imide or sulfonamide, or an acyl hydrazide. Yes *bypasses* the food question.
  This is the gate that separates **maltol from caffeine**, both food-listed
  heterocycles but Class II and Class III. The pattern splits all 13 observed
  yes-answers from all 14 no-answers.
- **Q33, the sulphonate-salt question**, which the original tree asks *last*.

Four more followed from the same analysis: **Q1**, the normal-body-constituent
question, which turns out to mean the building blocks of intermediary
metabolism (amino acids, nucleobases and nucleosides, hexoses, cholesterol) and
*not* the 439-structure body-constituent file — vitamins, fatty acids, citric
acid and the sugar alcohols all answer no; a **carbohydrate** question, since
cyclic sugars are Class I while the open-chain sugar alcohols are Class II;
**oxygen density** in place of an arbitrary size cutoff, one or two oxygen
functions reading as readily metabolised where three or more do not; and
**tautomer-aware matching**, because histidine is drawn either way.

A further pass fixed the remaining clusters: open-chain sugars (fructose is
drawn acyclic and is Class I, while the fully reduced sugar alcohols are
Class II — the carbonyl is the distinction); explicit hydrogens, since water
written `[H]O[H]` was counted as two heavy atoms and rejected as inorganic;
question order, because the body-constituent question genuinely comes first
and water answers it yes; two of my own sulphonate conditions, which were
wrong (a carboxylate alongside a sulphonate is fine, and the sulphonation
density requirement belongs only to large aromatic dyes); and routing sulphur
compounds away from the food question, which is why thiophene is Class III
despite being listed.

Together these took Class II from 22% to 81% across both validation sets
(8/37 → 30/37), and removed the ad-hoc guard the old body-constituent list
had needed.

**On self-confirmation.** The Q1 membership was identified from the recorded
paths, which include both validation sets, so that part is partly circular.
The effect is measurable and small: on the Appendix 1 set *no compound at all*
is decided by that list, so its 90.9% is entirely structural; on the food set
23 of 106 are, and excluding them the remaining 83 score 95.2% against 96.2%
overall.

### What the curated lists can and cannot do

Cramer's tree asks two questions that are list lookups rather than structural
rules: whether a substance is a normal constituent of the body, and whether it
is a common component of food. Both lists are available (439 and 108
structures) and both were tested as class-deciding rules. The results are worth
stating plainly, because the intuitive reading of each is wrong:

- **Food-list membership carries no class information on its own.** Every
  food-listed compound in the Appendix 1 set is Class II, which makes a
  "food ⇒ Class II" rule look excellent — it scores 86.4% there. It does not
  generalise: the same list holds ethanol, vanillin, benzoic acid and limonene
  (Class I) and caffeine (Class III). The rule was tested, rejected, and
  membership is now reported as context only.
- **Body-constituent membership needs a guard.** Applied as a plain
  "on the list ⇒ Class I" shortcut it lifts overall agreement from 55.6% to
  68.5%, but it cuts Class III accuracy from 100% to 73% and doubles the
  under-calls, because biotin and ascorbic acid are body constituents that the
  scheme still treats as Class III and Class II. Blocking the shortcut whenever
  the structure carries an enhanced-toxicity group, sulphur or a lactone gives
  69.4% while keeping Class III at 91%. That guarded form is what ships.

Neither list is a lever to pull. Each answers one narrow question inside a
branch of the tree, and using it as a global override trades accuracy in one
class for accuracy in another.

### Known limitations

- **Class II is now 30/37 across both sets** (81%, up from 22%), after
  recovering seven questions from the recorded decision paths. The weakest
  remaining area is the Appendix 1 Class II set (3/7), which contains the four
  compounds that turn on question 18 — and **question 18 could not be pinned
  down**. It decides Class I against Class II for 56 compounds, but offers only
  four positive examples against fifty-six negatives with no clean structural
  separator: malonic acid is Class I and oxalic acid Class II, one CH₂ apart;
  acetoin is Class I and 2-heptanone Class II, both methyl ketones. Enlarging
  the path set to 232 produced the same four positives. Guessing it would be
  exactly the overfitting the food-list rule was rejected for, so it is left
  unimplemented.
- The Cramer tree is implemented as 14 ordered questions covering the main
  branch points. The original 33-question tree distinguishes some narrow
  sub-cases (particular terpene skeletons, specific food heterocycles) that are
  approximated here rather than enumerated.
- The sulphonate-salt rule uses sulphonation density (two or more groups, or a
  small molecule) as a proxy for whether the substance is absorbed. This is a
  heuristic standing in for expert judgement, and it is what places one
  reference Class III dye in Class I.
- Normal-body-constituent recognition in Q2 uses an explicit list of about 50
  metabolites, matched on canonical SMILES. Close derivatives of a metabolite
  fall through to the structural questions instead.
- The Verhaar baseline QSARs are valid for neutral organics within roughly
  0 < logP < 6. Outside that window the app reports Class 5 or flags the
  estimate as an extrapolation.
- Inorganics, metals, organometallics, salts and mixtures are outside the
  applicability domain of every scheme here. They are flagged rather than
  silently classified.
- Nothing here is a substitute for measured data or for a regulatory
  assessment. Treat the output as a screening-level prioritisation aid.
        """
    )
