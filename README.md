# Toxicity decision trees — a pure-Python Streamlit app

Structure-based toxicity screening in the browser. Paste a SMILES, get the
Cramer class, the Verhaar aquatic mode of action, structural alerts for
genotoxicity, sensitisation and irritation, and the Kroes TTC tier — each with
the decision path that produced it.

Everything is Python: RDKit for the cheminformatics, Streamlit for the UI. No
Java, no server round-trips, no external services.

## Running it

```powershell
python -m venv .venv
```

```powershell
.venv\Scripts\python.exe -m pip install -r requirements.txt
```

```powershell
.venv\Scripts\python.exe -m streamlit run app.py
```

The app opens at <http://localhost:8501>.

## What it does

| Rulebase | Output |
| --- | --- |
| **Cramer rules** | Class I / II / III structural class of oral toxicity, with the full question path |
| **Verhaar scheme** | Aquatic mode of action (classes 1–5) plus a baseline-narcosis LC50 estimate for classes 1 and 2 |
| **Genotoxic carcinogenicity** | 38 DNA-reactivity structural alerts |
| **Non-genotoxic carcinogenicity** | 7 alerts for receptor-mediated and regenerative mechanisms |
| **Skin sensitisation** | 14 mechanistic protein-binding alerts (Michael addition, Schiff base, acyl transfer, SN2, SNAr) |
| **Skin / eye irritation** | 10 reactivity triggers — *unvalidated, see below* |
| **Kroes TTC** | Tiered threshold of toxicological concern, optionally compared against an exposure estimate |

Three tabs: **Single structure** (depiction with the matched alert atoms
highlighted, verdict cards, expandable decision paths), **Batch screening**
(paste or upload hundreds of structures, get a table, summary charts and a CSV),
and **Rulebase reference** (every alert with its SMARTS pattern, searchable).

## Layout

```
app.py              Streamlit interface
toxtree/
  chem.py           parsing, descriptors, SVG depiction
  alerts.py         the SMARTS alert library (83 alerts)
  cramer.py         Cramer decision tree
  verhaar.py        Verhaar scheme + baseline QSARs
  engine.py         orchestration + the Kroes TTC tree
  reflists.py       primary metabolites + the optional curated lists
  examples.py       30 reference structures
test_rules.py       regression tests
validate_cramer.py  score the Cramer tree against a reference set
validate_verhaar.py score the Verhaar scheme against a reference set
validate_alerts.py  score an alert rulebase against experimental data
tools/
  fetch_reference.py          Cramer Appendix 1 validation set
  fetch_reference_lists.py    body/food lists + the second Cramer set
  fetch_verhaar_reference.py  Verhaar validation set
  fetch_alert_references.py   LLNA sensitisation + Ames mutagenicity sets
  analyse_tree.py             recover tree topology from decision paths
```

## Testing

```powershell
.venv\Scripts\python.exe test_rules.py
```

Reference compounds are pinned to their expected Cramer class, Verhaar class and
genotoxicity verdict, so a change to a SMARTS pattern that shifts a known result
fails the suite.

## Validation against Cramer et al. (1978)

The Cramer implementation is validated against the 83 compounds of Appendix 1 of
the original paper — the classifications of the authors who devised the scheme.
The data is GPL-licensed where it is distributed, so it is not committed here;
fetch it and run the comparison:

```powershell
.venv\Scripts\python.exe tools\fetch_reference.py out\
```

```powershell
.venv\Scripts\python.exe validate_cramer.py out\cramer_reference.csv
```

Current result — 81 scored, 2 excluded for malformed structures in the source:

| | → I | → II | → III | n |
| --- | --- | --- | --- | --- |
| **Reference I** | **31** | 0 | 1 | 32 |
| **Reference II** | 4 | **3** | 0 | 7 |
| **Reference III** | 1 | 1 | **40** | 42 |

**91.4% overall** (Class I 97%, Class II 43%, Class III 95%), with 6 under-calls
and 1 over-call. `validate_cramer.py` prints every disagreement together with
the decision path that caused it, so each can be judged rather than counted.

Validation has driven almost every rule change here, and it corrected
assumptions baked into the first test suite: acrolein was pinned to Class III
on my own guess, but Cramer places the same functional class (trans-2-hexenal,
cinnamaldehyde) in Class I.

## Validating the Verhaar scheme

The Verhaar scheme has its own reference set — 121 compounds with the
mode-of-action class assigned in the Verhaar *et al.* evaluation:

```powershell
.venv\Scripts\python.exe tools\fetch_verhaar_reference.py out\
```

```powershell
.venv\Scripts\python.exe validate_verhaar.py out\verhaar_reference.csv
```

It started at **73.5%** and now stands at **87.2%** (117 scored):

| | → 1 | → 2 | → 3 | → 4 | → 5 | n |
| --- | --- | --- | --- | --- | --- | --- |
| **Reference 1 (inert)** | **40** | 0 | 0 | 0 | 0 | 40 |
| **Reference 2 (polar)** | 0 | **18** | 1 | 0 | 3 | 22 |
| **Reference 3 (reactive)** | 2 | 1 | **35** | 0 | 0 | 38 |
| **Reference 4 (specific)** | 2 | 0 | 4 | **9** | 2 | 17 |

Class 1 is perfect, Class 3 is 92%, and Class 4 is the weak spot at 53%.

Two of the fixes were outright bugs in rules I had written:

- An alert labelled *"nitroaromatic (two or more nitro groups)"* whose SMARTS
  matched a **single** nitro group, so every nitrophenol and nitrobenzene was
  called reactive instead of a polar narcotic. The real two-nitro test already
  ran separately, so the entry was pure noise.
- A **trifluoromethyl** group being matched as an "activated benzylic halide".
  A CF₃ group is about as inert as a substituent gets.

Three more were rules asserted without evidence, and the data contradicted all
three: triazines are Class 2 and not specifically acting; dinitrophenols and
thiocarbamates are Class 3, not Class 4; and **carbamates are Class 3** — all
six amides and carbamates in the set are, which also corrected a pinned test
value for carbaryl. Amides, carbamates and aromatic esters (4 of 4) were being
answered *"not possible to classify"*, which dropped from 13 compounds to 5.

**Class 4 cannot be improved further from this data.** Six of the eight misses
are organochlorine insecticides, and a halogen-count rule looks tempting until
it is tested: pentachlorophenol is Class 3 but 2,3,5,6-tetrachlorophenol is
Class 4; heptachlor is Class 3 but endrin is Class 4. At a threshold of four
halogens such a rule would break seven compounds to fix six, so it was
rejected. The organophosphates, which the scheme does separate cleanly, are all
classified correctly.

## Curated reference lists (optional)

Cramer's tree asks two questions that are list lookups rather than
structural rules. Those lists are GPL-licensed where they are distributed,
so they are fetched rather than committed:

```powershell
.venv\Scripts\python.exe tools\fetch_reference_lists.py toxtree\data\
```

The app runs without them and says so in the sidebar, but the Q22 food
question cannot be answered and Class II is then under-reported. The fetch
also produces a **second validation set** of 108 food-related structures
carrying the class each was assigned by the reference implementation — a
weaker oracle than Cramer's own paper, but it holds 30 Class II compounds
against Appendix 1's seven:

```powershell
.venv\Scripts\python.exe validate_cramer.py toxtree\data\foodmol_validation.csv
```

| | → I | → II | → III | n |
| --- | --- | --- | --- | --- |
| **Reference I** | **66** | 1 | 0 | 67 |
| **Reference II** | 2 | **27** | 1 | 30 |
| **Reference III** | 0 | 0 | **11** | 11 |

**104/108 = 96.3%** — Class I 99%, Class II 90%, Class III 100%.

## Validating the structural alerts

The alert rulebases are screening filters, not balanced classifiers: they
are built to miss as few actives as possible and pay for that in false
positives. So the headline figure is **sensitivity**, reported with
specificity and balanced accuracy rather than one concordance number that
would hide the trade-off.

```powershell
.venv\Scripts\python.exe tools\fetch_alert_references.py out\
```

```powershell
.venv\Scripts\python.exe validate_alerts.py out\llna_reference.csv --scheme sensitisation --per-alert
```

```powershell
.venv\Scripts\python.exe validate_alerts.py out\ames_reference.csv --scheme genotox --per-alert
```

| | sensitivity | specificity | balanced acc. | precision |
| --- | --- | --- | --- | --- |
| **Skin sensitisation** (209 LLNA) | 52% → **76%** | 68% → 51% | 60% → **64%** | **87%** |
| **Genotoxicity** (691 Ames) | **89%** | 53% | **71%** | 65% |

### A label that was the wrong way round

The first genotoxicity run scored 28% balanced accuracy — far worse than
chance, which is itself the clue. The ISSCAN mutagenicity field encodes
**3 as active and 1 as inactive**, the opposite of the intuitive reading,
and I had assumed the intuitive one. The giveaway was chemical rather than
statistical: 42 N-nitroso compounds had come out "non-mutagenic", and
N-nitroso is about as reliably mutagenic a class as exists. Checking
compounds with settled outcomes confirmed it — benzo[a]pyrene, benzidine,
2-acetylaminofluorene and aflatoxin B1 all carry 3; thiourea, caffeine and
ascorbic acid all carry 1. With the polarity corrected the same alert set
scores 89% sensitivity, and the per-alert precisions become chemically
sensible: N-nitroso 93%, nitroaromatic 86%, polycyclic aromatic 85%,
epoxide 83%.

The lesson is worth keeping: a validation result far below chance usually
means the labels are inverted, not that the model is extraordinarily bad.

### What the sensitisation data changed

The alert set caught barely half the sensitisers. The misses grouped
cleanly by mechanistic domain, which is what made them fixable:

- **Azlactones and lactones** — the largest acyl-transfer group in the set,
  and entirely absent (5 of 26 acyl sensitisers were being caught).
- **Alkyl chlorides** — the SN2 pattern accepted bromide and iodide but not
  chloride, so the long-chain alkyl chlorides were missed.
- **1,2- and 1,3-dicarbonyls** — these condense with lysine far more
  readily than an isolated ketone.
- **Alkenylphenols** — the eugenol/isoeugenol family, quinone-methide
  precursors that the catechol pattern alone does not reach.
- **Isothiazolinones** — the preservative sensitisers.

Broadening the phenol alert initially flagged every paraben and salicylate,
because the pattern matched the *carboxyl* carbon; a quinone methide needs
benzylic hydrogens, and requiring those recovered 15 points of specificity.
Saccharin was likewise being read as an isothiazolinone until the pattern
was made to require divalent sulphur.

Sensitivity is the right thing to optimise here, but it is bought with
specificity, which fell from 68% to 51%. Balanced accuracy still improved,
and an alert that fires is right 87% of the time.

## Validation status of each scheme

| Scheme | Reference data | Result |
| --- | --- | --- |
| Cramer classes | 83 (Cramer Appendix 1) | 91.4% concordance |
| Cramer classes | 108 (reference implementation) | 96.3% concordance |
| Verhaar scheme | 121 (Verhaar *et al.*) | 87.2% concordance |
| Genotoxicity alerts | 691 (ISSCAN Ames) | 89% sensitivity |
| Sensitisation alerts | 209 (LLNA) | 76% sensitivity |
| **Irritation alerts** | **none found** | **unvalidated** |
| **Non-genotoxic carcinogenicity** | **none found** | **unvalidated** |

### Why irritation is unvalidated

There is no experimental irritation set to validate against. The Toxtree
repository's eye-irritation plugin ships rule classes and a 417-byte test
stub but no data, and every dataset in that repository over 20 KB is one
already used above or unrelated. No public machine-readable set surfaced
either: the BfR rulebase's 199-compound external validation set does not
appear to be published, and ChemSkin DB and the VEGA QMRF material are not
available as parseable structure files.

Two things follow, and both are stated here rather than papered over.

First, **these alerts are the least grounded of the five**. The published
approach to local irritation is largely *physicochemical* — the BfR
decision support system works from logP, molecular weight, melting point
and pKa cut-offs with structural rules layered on top. This implementation
is a SMARTS reactivity list, which does not capture that at all. It should
be read as a flag for reactive functional groups, not as a prediction of
irritancy or corrosivity.

Second, what *can* be checked without labels was checked. An alert set
that fires on most of chemical space carries no information whatever the
labels say, so the base rates were measured across the 941 unique
compounds from the other three reference sets:

- at least one alert fires on **29%** of them, and no single alert exceeds
  an 11% firing rate, so the set does discriminate;
- the distribution is sane — 239 compounds fire one alert, 28 fire two,
  4 fire three.

That diagnostic found no pathology, which is a genuine negative result and
not a substitute for validation.

It did surface one defect that needs no labels to call wrong: the amino
acids were firing both the acid and the base alert — glycine two, tyrosine
three. An alpha-amino acid is zwitterionic at skin pH, so neither group
behaves as a free acid or base, and both alerts now exclude that case.
That fix rests on chemistry rather than on measurement, and is labelled as
such in the code.

## Recovering the tree from decision paths

Both reference sources record the *path* each compound took, not just its class
— caffeine is `1N,2N,3N,5N,6N,7Y,8N,10N,11Y,33N`. About 190 such paths are
enough to reconstruct the shape of the original 33-question tree without
reading any implementation of it:

```powershell
.venv\Scripts\python.exe tools\analyse_tree.py out\ toxtree\data\
```

```powershell
.venv\Scripts\python.exe tools\analyse_tree.py out\ toxtree\data\ --q 22
```

The first prints the topology — which question follows which on each answer,
and which terminate in which class. The second lists the compounds that
answered a question yes and no, so its rule can be inferred from structures
rather than guessed. Three questions were identified this way and implemented:

- **Q22 — "is it a common component of food?"** Answering **yes gives Class II**,
  not Class I. All 28 observed yes-answers are Class II (ascorbic acid, citric
  acid, mannitol, niacin, maltol, vitamin D); the no-answers are industrial
  chemicals. This is why the earlier blanket food rule failed — the question is
  asked at *one node*, not globally.
- **Q11 — "does a ring heteroatom bear an acyl group?"** (cyclic amide, urea,
  imide or sulfonamide, or an acyl hydrazide). Answering yes *bypasses* the
  food question and heads for Class III. This is the gate that separates
  **maltol from caffeine** — both food-listed heterocycles, but Class II and
  Class III. The pattern splits all 13 observed yes-answers from all 14 no.
- **Q33 — the sulphonate-salt question**, which the original tree asks *last*,
  not early as a previous version of this code did.

Four more questions followed from the same analysis:

- **Q1 — "is it a normal constituent of the body?"** is far narrower than the
  439-structure body-constituent file suggests. Across 232 paths, 65 compounds
  answer yes and every one is Class I; they form one tight group — the twenty
  amino acids, the nucleobases and nucleosides, hexoses and their uronic acids,
  cholesterol, water, acetone, phosphoenolpyruvate. Vitamins, fatty acids,
  citric acid and the sugar alcohols all answer **no**, though every one of them
  is in that file. Using the file for this question is what made the old guard
  necessary; with the right membership no guard is needed, because cysteine and
  methionine are in the set and biotin and riboflavin are not.
- **A carbohydrate question.** Cyclic sugars are Class I (lactose, maltose,
  sucrose, raffinose), while the open-chain sugar alcohols mannitol and sorbitol
  are Class II — so the test is for the ring form, not for a run of hydroxyls.
- **Oxygen density rather than molecular size.** One or two oxygen functions
  still read as readily metabolised (ethanol, propylene glycol, a C24 fatty
  acid, acetoin); three or more, or an α-dicarbonyl, do not, which is where the
  sugar alcohols, citric acid, oxalic acid and diacetyl part company with them.
  This replaced an arbitrary heavy-atom cutoff that was demoting long-chain
  fatty acids and lycopene.
- **Tautomer-aware matching**, because histidine is drawn as either imidazole
  tautomer and the nucleobases have several accepted forms.

A further pass fixed the remaining clusters: **open-chain sugars** (fructose and
galacturonic acid are drawn acyclic and are Class I, while the fully reduced
sugar alcohols are Class II — the carbonyl is the distinction); **explicit
hydrogens**, since water written `[H]O[H]` was being counted as two heavy atoms
and rejected as inorganic; **question order**, because the body-constituent
question genuinely comes first and water answers it yes; **the sulphonate
conditions**, where two of my own guards were wrong (a carboxylate alongside a
sulphonate is fine — the azo dyes carry both — and the sulphonation-density
requirement belongs only to large aromatic dyes, not to alkyl sulphates); and
**routing sulphur compounds away from the food question**, which is why
thiophene is Class III despite being listed.

Together these took **Class II from 22% to 81%** across both sets (8/37 → 30/37)
while Class I rose to 97%/99% and Class III to 95%/100%.

**On self-confirmation:** the Q1 membership was identified from the recorded
paths, which include both validation sets, so that part is partly circular. The
effect is measurable and small. On the Appendix 1 set **no compound at all** is
decided by that list — its 91.4% is entirely structural. On the food set 24 of
108 are; excluding them the remaining 84 score **95.2%**, against 96.3% overall.
The structural rules carry the result, not the membership.

Two earlier findings still stand, and both are cases where the obvious reading
of the data was wrong:

- **Food-list membership alone carries no class information.** Every food-listed
  compound in Appendix 1 is Class II, so a blanket "food ⇒ Class II" rule scores
  86.4% there — but the same list holds ethanol, vanillin and limonene (Class I)
  and caffeine (Class III). Only the gated form at Q22 works.
- **The body-constituent file is the wrong list for Q1.** Used as a plain
  shortcut it cuts Class III accuracy from 100% to 73%, because biotin,
  ascorbic acid and riboflavin are all in it and the scheme still rates them
  II or III. An explicit guard patched around that for a while; identifying
  the narrow membership Q1 actually asks for removed the need for one.

## Provenance and scope

The rules here are written from the published literature — Cramer, Ford & Hall
(1978); Verhaar *et al.* (1992) with the Enoch *et al.* (2008) refinements; the
Benigni & Bossa alert set; Kroes *et al.* (2004) — not transcribed from any
existing program. The Cramer tree is implemented as 14 ordered questions
covering the main branch points rather than the original 33, so narrow sub-cases
are approximated; the **Method notes** tab in the app lists the known
limitations in full.

Structural alerts are hypothesis generators. An alert that fires means the
molecule contains a fragment associated with a mechanism, not that the molecule
is toxic; no alert means nothing in *this* library matched. Inorganics, metals,
salts and mixtures are outside the applicability domain of every scheme and are
flagged rather than silently classified. This is a screening-level
prioritisation aid, not a substitute for measured data or regulatory assessment.
