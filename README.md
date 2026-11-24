# Toxicity decision trees — a pure-Python Streamlit app

Structure-based toxicity screening in the browser. Paste a SMILES, get the
Cramer class, the Verhaar aquatic mode of action, structural alerts for
genotoxicity and sensitisation, flags for reactive groups, and the Kroes TTC
tier — each with
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
| **Genotoxic carcinogenicity** | 35 DNA-reactivity structural alerts |
| **Non-genotoxic carcinogenicity** | 10 alerts for receptor-mediated, hormonal and regenerative mechanisms |
| **Skin sensitisation** | 23 mechanistic protein-binding alerts (Michael addition, Schiff base, acyl transfer, SN2, SNAr) plus the pro-haptens |
| **Reactive group flags** | 10 protein-reactive and corrosive triggers. **Not** an irritancy prediction — see below |
| **Kroes TTC** | Tiered threshold of toxicological concern, optionally compared against an exposure estimate |

Three tabs: **Single structure** (depiction with the matched alert atoms
highlighted, verdict cards, expandable decision paths), **Batch screening**
(paste or upload hundreds of structures, get a table, summary charts and a CSV),
and **Rulebase reference** (every alert with its SMARTS pattern, searchable).

## Layout

```
LICENSE             MIT
app.py              Streamlit interface
toxtree/
  chem.py           parsing, descriptors, SVG depiction
  alerts.py         the SMARTS alert library (92 alerts)
  cramer.py         Cramer decision tree
  verhaar.py        Verhaar scheme + baseline QSARs
  engine.py         orchestration + the Kroes TTC tree
  reflists.py       primary metabolites + the optional curated lists
  readacross.py     nearest reference neighbours + evidential-support index
  examples.py       30 reference structures
test_rules.py       regression tests
validate_cramer.py  score the Cramer tree against a reference set
validate_verhaar.py score the Verhaar scheme against a reference set
validate_alerts.py  score an alert rulebase against experimental data
compare_toxtree.py  three-way comparison against the Java Toxtree
tools/
  fetch_reference.py          Cramer Appendix 1 validation set
  fetch_reference_lists.py    body/food lists + the second Cramer set
  fetch_verhaar_reference.py  Verhaar validation set
  fetch_alert_references.py   LLNA sensitisation + Ames mutagenicity sets
  build_oracle.py             build the Toxtree regression fixture
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

Current result — 77 scored; 2 excluded for malformed structures in the source
and 4 for being outside the domain of the scheme (two organotins, a beryllium
salt and a mercurial):

| | → I | → II | → III | n |
| --- | --- | --- | --- | --- |
| **Reference I** | **31** | 0 | 1 | 32 |
| **Reference II** | 4 | **3** | 0 | 7 |
| **Reference III** | 1 | 1 | **36** | 38 |

**90.9% overall** (Class I 97%, Class II 43%, Class III 95%), with 6 under-calls
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
| **Reference I** | **65** | 1 | 0 | 66 |
| **Reference II** | 2 | **27** | 1 | 30 |
| **Reference III** | 0 | 0 | **10** | 10 |

**102/106 = 96.2%** — Class I 98%, Class II 90%, Class III 100%. Sodium
chloride and water are not scored; they are not organic molecules and the
scheme does not claim to classify them.

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
| **Skin sensitisation** (209 LLNA) | 52% → 76% → **92%** | 68% → 51% → 49% | 60% → 64% → **70%** | **88%** |
| **Genotoxicity** (672 Ames) | 89% → **90%** | 53% → **59%** | 71% → **75%** | 65% → **68%** |

### Two alerts that were in the wrong rulebase

Per-alert precision is the useful diagnostic here, because an alert set is a
union and a bad member costs specificity without buying sensitivity. Two
members were doing exactly that on the 691 Ames compounds:

| | true pos | false pos | precision |
| --- | --- | --- | --- |
| SA35 Thiourea / thioamide | **0** | **10** | **0%** |
| SA33 Halogenated dibenzodioxin / biphenyl | 1 | 5 | 17% |

A pattern that fires ten times and is wrong ten times carries no information,
and the diagnosis was mechanistic rather than statistical: both rationales, as
written in the code, described mechanisms that are *not* DNA-reactive. A
thiocarbamide inhibits thyroid peroxidase; a halogenated dioxin is a persistent
AhR ligand and tumour promoter. Both are real carcinogenic mechanisms and
neither is genotoxicity, so they were in the wrong one of the two rulebases.
The thiourea case was already settled elsewhere in this file — thiourea is one
of the compounds used above to confirm that the ISSCAN labels were inverted,
*because* it is a known non-mutagen.

They now sit in the non-genotoxic set, SA33 as NG8 and SA35 folded into NG6,
whose thiocarbamide arm was widened from `[NX3]` to `[NX3,#6]` so the
thioamides keep their coverage. That widening rests on chemistry rather than
measurement — thioacetamide and its relatives are non-genotoxic
hepatocarcinogens — and is labelled as such in the code, since the
non-genotoxic set still has no reference data.

The result is 12 false positives removed at **no cost in sensitivity at all**:
the misses stay at 36 of 340, exactly the same compounds. Specificity goes
53.0% → 56.4%, precision 64.8% → 66.5%, balanced accuracy 71.2% → 72.9%. The
non-genotoxic set, which cannot be scored directly, was checked for shape the
same way the irritation set is: across the 831 unique compounds of the other
three reference sets it fires on 21.7%, no single alert above 10.0%, with the
later arrivals at 2.6% (NG6, widened), 0.8% (NG8) and 0.5% each for NG9 and
NG10.

### Counting only the alerts that actually cost something

Per-alert precision is still the wrong measure for deciding what to remove,
because an alert set is a union: a pattern with many false positives that are
*also* flagged by something else costs nothing, and tightening it changes no
verdict. What matters is the **solo** count — the compounds an alert alone is
responsible for. The two measures disagree sharply. SA38 (propargyl / allyl
halide) sits at 27% precision with 8 false positives and **none** of them solo;
it is free. Meanwhile:

| alert | solo actives | solo inactives | if removed |
| --- | --- | --- | --- |
| SA10 α,β-unsaturated carbonyl / nitrile / sulphone | 3 | 14 | +11 spec, −3 sens |
| SA29 aromatic mono- or dialkylamine | 5 | 8 | +3 spec, −5 sens |
| SA4 monohaloalkene | 5 | 7 | +2 spec, −5 sens |
| SA2b dialkyl sulfate / trialkyl phosphate | 5 | 7 | +2 spec, −5 sens |
| SA31 aromatic diazo | 7 | 8 | +1 spec, −7 sens |

Only the first is worth acting on, and for the same reason as the two above:
**a Michael acceptor is protein-reactive, not DNA-reactive.** Thiols and
glutathione consume it long before DNA does, which is why the identical pattern
scores 96% precision as the sensitisation alert SS1 and 36% here. The 14
compounds it alone condemned are acrylamide, acrylic acid, the acrylate and
methacrylate esters, ascorbic acid, curcumin, caffeic acid, carvone, isophorone
and the tetracyclines — none mutagenic.

Two of its three actives were pyrrolizidine alkaloids, and it caught them for
the wrong reason: what makes lasiocarpine and riddelline genotoxic is hepatic
dehydrogenation to a pyrrolic ester, not the enoate. **SA39** now describes
that, and keeps them. Its one Ames-negative hit, monocrotaline, is itself an
established genotoxic hepatocarcinogen the bacterial assay misses. The third,
acrylonitrile, is given up: it is genuinely Ames-positive through its epoxide
and no pattern separates it from methacrylonitrile, which is not.

The rows below SA10 were all rejected. Trading five or seven actives for two or
three inactives is the wrong direction for a screening filter.

**SA16 was kept against the measurement**, which is the one case here where the
data should lose. Removing it buys 2 points of specificity, because the only
two compounds it alone flags — urethane and methyl carbamate — are both
Ames-negative. But urethane is the textbook Ames-negative genotoxic
*carcinogen*: it is activated in the liver to vinyl carbamate epoxide, which a
bacterial assay cannot model. Dropping it would fit this rulebase to the Ames
endpoint at the expense of the endpoint it is for.

Net: sensitivity 89.4% → 90.1%, specificity 56.4% → 58.9%, balanced accuracy
72.9% → 74.5%.

### Compounds the schemes never claimed to answer for

The README has always said that inorganics, metals, salts and mixtures are
outside the applicability domain of every scheme here, and that they are
flagged rather than silently classified. The validators did not know that, and
scored them anyway — so a deliberate refusal was being counted as a wrong
answer. `domain_exclusion()` now reports them and they are left unscored, with
the count and the reason printed, as the unparsable structures already were.

The rule is narrower than the sentence it implements, in two ways that matter.

**Metal compounds are excluded everywhere.** No SMARTS pattern in the library
addresses a metal centre and neither tree asks a question that applies to one.
This is what removes the two organotins, the beryllium salt and the mercurial
from Appendix 1, and potassium dichromate from the LLNA set.

**Carbon-free compounds are excluded only from the trees.** Cramer and Verhaar
ask about organic functional groups and can say nothing about hydrazine or
sodium nitrite. The alert sets can: SA13 *is* a hydrazine pattern and SA34 a
peroxide one. Excluding inorganics from an alert run would be ducking, since
potassium bromate and sodium nitrite are real mutagens this library has no
pattern for, and they stay counted as the misses they are.

**Mixtures are not excluded at all**, despite the sentence. Counting fragments
does not identify a mixture: the reference files write racemates as two
enantiomers in one record (DL-menthol, DL-ethionine) and commercial blends as
their components (toluene diisocyanate), and the rules answer identically for
each fragment. Excluding those would have quietly dropped genuine actives.

Worth noting what this did *not* do. Both Cramer figures went **down** —
91.4% → 90.9% on Appendix 1 and 96.3% → 96.2% on the food set — because the
compounds removed were ones the tree had been getting right. The alert
sensitivities rose, 89.4% → 90.1% for genotoxicity and 91.7% → 92.2% for
sensitisation. A change that only ever improved the numbers would be a reason
for suspicion.

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

### The pro-haptens, found by comparison

At 76% sensitivity this set was still missing a quarter of the sensitisers,
and nothing in the LLNA data said what they had in common. Running Toxtree's
own rulebase over the same file did: it caught 34 that this missed, and they
grouped by *mechanism*. Every alert above is an electrophile as drawn. None of
the misses was. They become one in the skin — a **pro-hapten** oxidised by
cutaneous enzymes, or a **pre-hapten** oxidised by air — and that whole
category was absent:

- **Alkenylphenols and their ethers** (SS13) — anethole, the isoeugenol
  ethers, 4-allylanisole. Quinone-methide precursors.
- **Allylic and benzylic primary alcohols** (SS14) — cinnamic alcohol and
  geraniol act through their aldehydes.
- **Polycyclic aromatics** (SS15) — benzo[a]pyrene and
  dimethylbenz[a]anthracene, by the same activation that makes them genotoxic.
- **Beta-lactams** (SS16), **aliphatic polyamines** (SS17),
  **N-nitrosoureas** (SS18), **N-methylol formaldehyde releasers** (SS19),
  **heteroaromatic thiols** (SS20) and **activated esters** (SS21).

Sensitivity went **76% → 92%** for **one** new false positive in 41 inactives:
specificity 51% → 49%, precision 87% → 88%, balanced accuracy 64% → 70%.

Two details are worth keeping, because both are about restraint.

**The phenol alert had a SMARTS precedence bug, and the comment hid it.** It
claimed to cover "the eugenol/isoeugenol family"; it covered eugenol only. The
benzylic branch read `[CX4;!$(C=O),CX3;!$([CX3]=[OX1])]`, and in SMARTS `,`
binds *tighter* than `;`, so that parses as `CX4 AND (not-C=O OR CX3) AND
not-CX3=O` — the leading `CX4` makes the `CX3` alternative unreachable.
Eugenol's benzylic carbon is saturated and matched; isoeugenol's propenyl
carbon is not and never could. A test now pins isoeugenol.

**The ester alert is deliberately narrower than the misses would justify.**
Four of the remaining misses are esters, and the obvious pattern — any aroyl
ester — gains two of them and costs four inactives: methyl paraben, propyl
paraben, methyl salicylate and diethyl phthalate, the same false positives the
phenol alert threw up when it was first broadened. Restricting the alert to
the ester *of a phenol* gains three and costs none, which is also the
chemistry: a phenoxide is a far better leaving group than an alkoxide, so an
aryl ester acylates lysine and an alkyl benzoate does not. Benzyl benzoate is
an alkyl ester and is given up as a miss rather than paid for with the
parabens. A cyclic-aryl-ester pattern was rejected the same way, for costing
coumarin and 6-methylcoumarin.

What this does cost is breadth of firing. Across the 831 unique compounds in
the *other* three reference sets the sensitisation rulebase now fires on 42.5%
rather than 35.0%, with no single alert above 9.9%. That is the price of
reading pro-haptens as alerting, and it is the same trade the specificity
figure reports.

## Comparing against Toxtree

Everything above scores this implementation against reference data. It does not
say how the Java [Toxtree](https://sourceforge.net/projects/toxtree/) answers
the same compounds, which is a different question: a miss both programs make is
the *scheme's* limit, while a miss only this one makes is a gap here.

`compare_toxtree.py` answers it. Toxtree has an undocumented headless mode, so
it can be driven over a file and its verdicts read back:

```powershell
java -jar Toxtree-3.1.0.1851.jar --noui -i in.sdf -o out.sdf -m toxTree.tree.cramer.CramerRules
```

The output SDF carries the class *and* the decision path, in the same format
`analyse_tree.py` parses — so paths can now be generated for any compound, not
only the ~190 that happen to be recorded in the test fixtures.

Toxtree is GPL-licensed and is not distributed here. Fetch the plain zip (no
installer, no Maven build needed) from SourceForge, unpack it, and point the
script at it with `TOXTREE_JAR`; `TOXTREE_JAVA` sets the JRE if `java` is not on
the path. Toxtree needs Java 7 or newer — use 11 rather than 17+, since the
stack is 2014-era CDK and strong module encapsulation can break its reflective
plugin loading.

```powershell
$env:TOXTREE_JAR = "C:\path\to\Toxtree\Toxtree-3.1.0.1851.jar"
```

```powershell
.venv\Scripts\python.exe compare_toxtree.py cramer path\to\foodmol.sdf
```

Structures are handed to both sides as **the molblock drawn in the file**, never
as regenerated SMILES. That is not fussiness: Toxtree answers two Cramer
questions by InChI lookup against curated lists, and the same compound drawn
differently resolves differently. Ascorbic acid drawn by RDKit answers Q1 *yes*
and lands in Class I; as drawn in `foodmol.sdf` it answers *no* and lands in
Class III. Re-deriving structures would measure the drawing, not the rules.

### Results

| scheme | set | ours | Toxtree | module |
| --- | --- | --- | --- | --- |
| Cramer | 108 foodmol | **96.3%** | 77.8% | `toxTree.tree.cramer.CramerRules` |
| Cramer | 78 Appendix 1 † | **91.0%** | 84.6% | `toxTree.tree.cramer.CramerRules` |
| Verhaar | 117 Verhaar2000 | **87.1%** | 53.4% | `verhaar.VerhaarScheme` |
| Sensitisation | 208 LLNA | **92% sens** | 91% sens | `…SkinSensitisationPlugin` |
| Genotoxicity | 665 ISSCAN Ames | **89% sens** | 87% sens | `mutant.BB_CarcMutRules` |
| Non-genotox carcinogenicity | 886 comparable | — | — | `mutant.BB_CarcMutRules` |
| Skin irritation | 129 comparable | — | — | `sicret.SicretRules` |
| Eye irritation | 37 comparable | — | — | `eye.EyeIrritationRules` |

† Appendix 1 is distributed as Java test fixtures rather than an SDF, so its
structures have to be rebuilt from SMILES. That is the representation caveat
below in action, and it makes this row weaker evidence than the others; it is
used mainly to read *which questions* decided a compound, which are
structural.

These figures are deliberately scored on every record both programs answer
for, including the metal compounds and inorganics that the validation section
above leaves unscored — the point here is to compare two programs, and Toxtree
answers for sodium chloride as readily as this does. That is why the Cramer
figure reads 96.3% of 108 here and 96.2% of 106 there.

Five of those seven need reading carefully, and the caveats matter more than
the numbers.

**Cramer is the one clean win.** Like for like, same file, same module. The
reference labels are the classes stored in `foodmol.sdf`, and Toxtree 3.1.0
reproduces only 84 of its own 108 stored labels — so that oracle is a
*historical* run of the program, not the program as it ships. 22 compounds are
right here and wrong in current Toxtree; 2 the other way (cystine among them).

**Verhaar is not like for like, so the 34-point gap is mostly not merit.**
Toxtree 3.1.0 ships only `verhaar.VerhaarScheme`, the original 1992 scheme; the
`verhaar2` plugin carrying the Enoch (2008) refinements is absent from the
distribution even though it is listed in `toxtree-plugins.properties` — and
`Verhaar2000.sdf` is verhaar2's own test fixture. This implementation follows
the refined scheme. Toxtree answers Class 5, *cannot classify*, for 30 of the
117, the phenols among them, which is precisely what the refinements were made
to fix. The script prints this caveat on every Verhaar run.

**Sensitisation was a loss, and that is what the comparison was for.** Toxtree
caught 91% of the LLNA sensitisers against 76% here, and the 34 compounds it
flagged and this did not turned out to group by mechanism rather than scatter,
which made them fixable — see *The pro-haptens, found by comparison* above.
With those mechanisms added the two are level and this is marginally ahead:

| | sensitivity | specificity | balanced acc. | precision |
| --- | --- | --- | --- | --- |
| ours (23 alerts) | **92%** | **49%** | **70%** | **88%** |
| Toxtree (OASIS domains) | 91% | 46% | 69% | 87% |

The compounds Toxtree still catches and this does not fell from 34 to 8.

**Genotoxicity is a near-draw on sensitivity and a loss on specificity.** Over
665 graded ISSCAN compounds:

| | sensitivity | specificity | balanced acc. | precision |
| --- | --- | --- | --- | --- |
| ours (35 alerts) | **89%** | 62% | 75% | 69% |
| Toxtree (Benigni/Bossa) | 87% | **72%** | **79%** | **75%** |

Sensitivity is what this rulebase is tuned for, and 35 alerts match the full
rulebase on it — 89% against 87%. The cost is the compounds flagged here that
Toxtree leaves alone, which is where the remaining 14 points of specificity go;
two of the worst offenders were removed from this set entirely, as described
under *Two alerts that were in the wrong rulebase*. Toxtree's
module also runs semi-empirical QM descriptors per molecule, making it far
slower, and reports an internal tree error on a minority of records.

The most interesting number in that run is neither side's: **106 of 665 (15.9%)
are called wrong by both programs, identically.** That is the Benigni/Bossa
approach's own floor on this data, not a defect in either implementation — and
it is the figure no amount of work on these SMARTS patterns will move.

**Non-genotoxic carcinogenicity gets its first external check**, and it comes
with a warning about how to read it. The two agree on 87.8% of 886 compounds —
but 670 of those are cases where *neither* fires, so the headline is mostly
mutual silence. Among the 216 compounds where either one fires they agree on
108, which is **50.0%**. That is the number worth quoting. It began at 37%,
and the climb is the whole content of this section: the comparison was the
first evidence of any kind about this rulebase, which until then had none.

The misses came in families, which is what made them fixable. Of the 44 that
Toxtree flagged and this did not:

| family | alert | recovered | newly disagreeing |
| --- | --- | --- | --- |
| 2-ethylhexyl esters | NG9 | 3 | 0 |
| thiouracils and thioguanines | NG6 (widened) | 5 | 0 |
| safroles and piperonyl synergists | NG10 | 5 | 0 |
| chlorinated ethylenes, hexachlorocyclohexanes | NG2 (widened) | 8 | 0 |
| ortho-phenylphenol | NG11 | 2 | 0 |
| dithiocarbamates and thiurams | NG6 (widened) | 3 | **3** |

**NG6 was simply too narrow for its own mechanism.** Its thiocarbamide arm
required an aliphatic `[NX3][CX3](=[SX1])`, but the thiouracils and
thioguanines are written as aromatic thiones — `c(=S)[nH]` — which that can
never reach. Half the thyroid story was unreachable by construction.

**NG2 was the wrong shape.** Two halogens on one saturated carbon misses both
families: the chlorinated ethylenes are sp², so no `CX4` matches, and the
hexachlorocyclohexanes carry one halogen per carbon rather than two on any
one of them. Three halogens across a double bond is the threshold that
separates tri- and tetrachloroethylene, hexachlorobutadiene and
tetrafluoroethylene from vinylidene chloride, dichlorvos and deltamethrin —
which carry two, and which Toxtree does not flag either.

**The dithiocarbamates are the one case where agreeing with Toxtree got
worse,** and they were added anyway. The arm gains three compounds Toxtree
flags and newly flags three it does not: sulfallate, sodium diethyldithio-
carbamate and dimethyldithiocarbamic acid. Those three are themselves
dithiocarbamates. The disagreement therefore looks like inconsistency on the
other side rather than over-reach here, and treating a mechanistic class
uniformly was preferred to scoring a point. Agreement is not ground truth,
and this is the case that makes the difference visible.

What remains is the other direction. This implementation flags 93 compounds
Toxtree does not, most of them on NG1 (halogenated benzene or phenol) and NG2,
which are broad: 1,2-dichlorobenzene, chlorodibromomethane,
1,1,1-trichloroethane. Whether that is over-reach cannot be settled, because
nothing can score precision on this endpoint — there is still no experimental
reference set, only a second opinion.

**Irritation is the answer to the open question above, and it is not
reassuring.** There are still no experimental labels, so this is agreement
between two programs with neither known to be right — but they agree on 50% of
the 129 compounds SICRET will answer for, and 60% of the 37 the eye rulebase
will. Toxtree declines to answer ("Unknown") for 80 and 172 respectively,
because its rules want physicochemical properties the file does not carry.
Chance-level agreement is what the README predicted from first principles: a
SMARTS reactivity list and a logP/MW/pKa decision system are not two
implementations of one thing. Treat that tab accordingly.

### The comparison as a regression test

The expectations in `test_rules.py` were written by the same person who wrote
the rules, so they encode the same assumptions — and twice they have been
wrong in the same direction. Acrolein was pinned to Cramer Class III against
the classification its own authors give; methyl acrylate was pinned as
genotoxicity-alerting on an alert that measurement later moved out of the
rulebase entirely. Thirty hand-written pins cannot catch that.

`tools/build_oracle.py` records, for every compound in the reference files,
what Toxtree says and what this implementation said when the fixture was
built:

```powershell
.venv\Scripts\python.exe tools\build_oracle.py out\ toxtree\data\toxtree_oracle.csv
```

`test_against_the_toxtree_oracle` then recomputes this side and fails on any
compound whose verdict has moved, naming it and saying whether it moved
*toward* Toxtree or away:

```
sensitisation/Isoeugenol: no alert -> alerting (Toxtree says alerting, toward Toxtree)
```

So changing a rule now means regenerating the fixture deliberately. That
deliberate step is exactly what was missing when the two bad pins went in.

The fixture stores canonical SMILES rather than molblocks, which is a
compromise: Toxtree answers two Cramer questions by InChI lookup and can give
a different answer for the same compound drawn differently. The build guards
against it by computing each verdict from the molblock *and* from the round
-tripped SMILES and dropping any compound where the two disagree — a
regression test that cannot reproduce its own baseline is worse than none.

A regression test that has never failed is a test nobody has checked, so this
one was verified by breaking a rule on purpose: narrowing NG10 so the
benzodioxoles stop matching. It named all five affected compounds and the
direction —

```
nongenotox/Safrole: alerting -> no alert (Toxtree says alerting, away from Toxtree)
```

— and passed again once the rule was restored. The fixture currently holds
**2624 verdicts** across seven schemes; 3 compounds were dropped because their
SMILES round trip would not reproduce the molblock answer, and 40 because
Toxtree returned no record for them.

Like everything else `tools/` produces, the fixture is GPL-derived and is not
committed. Without it the test reports itself skipped rather than failing.

### What the comparison found in this code

Feeding real SDF files through the rulebases surfaced two defects that the
SMILES-only test suite could not reach:

- `compute_descriptors` raised `RingInfo not initialized` on any multi-fragment
  record, because `GetMolFrags(sanitizeFrags=False)` returns fragments without
  ring perception and `NumRings()` throws rather than returning zero. Any SDF
  containing a salt crashed the batch tab. Fixed in `chem.py`.
- `parse_structure` stripped its input before testing for `V2000`, so a molblock
  whose title line is blank had that line removed, every following line shifted
  by one, and the block failed to parse. This was not an edge case: a molblock's
  title is blank unless a writer fills it in, and RDKit's own `MolToMolBlock`
  leaves it empty — so *every* RDKit-generated molblock pasted into the app was
  rejected. Fixed in `chem.py`: the input is classified before it is touched,
  and a molblock has only its trailing end stripped. The title line, when a
  writer did fill it in, now names the record.

## Read-across and evidential support

Every verdict above is a rule firing. On its own a Class III for
cinnamaldehyde and a Class III for something unlike anything ever tested
read identically, and they should not. The **Evidential support** panel
closes that gap: for each scheme it finds the query's nearest neighbours
in the reference sets used for validation and asks whether their
*measured* outcomes agree with what the rules just said.

It needs the reference data, which the `tools/` scripts fetch into
`toxtree/data/`. Without it the app runs and says so in the sidebar.

The index combines four components, weighted as shown:

| component | weight | what it asks |
| --- | --- | --- |
| Concordance | 0.40 | do the nearest compounds' measured outcomes match this verdict? |
| Similarity | 0.30 | how close is the nearest reference compound? |
| Descriptor range | 0.15 | are MW and logP inside the reference range? |
| Fragment novelty | 0.15 | how much of the structure never appears in the data? |

Above 0.70 is reported as well supported, 0.45 to 0.70 as partly
supported, below that as little support, with the components always shown
separately — the single number is the least informative part.

### What it is not

The idea comes from VEGA's applicability domain index, but the meaning is
deliberately different and the difference matters. VEGA scores statistical
models, where the index asks whether a prediction interpolates or
extrapolates from training data. These decision trees are deterministic
rules that are **not fitted to these compounds at all**, and they answer
identically for a molecule in the reference set and one unlike anything
ever tested.

So this is not an applicability domain in the QSAR sense. It measures how
much checked evidence stands behind a particular answer. A low score does
not mean the rule misfired; it means nothing in the reference data would
have caught it if it had.

The feature earned its place on first run by finding a mistake in this
repository's own example set: the structure labelled benzo[a]pyrene was
matching *pyrene* at similarity 1.00 and the real benzo[a]pyrene at 0.50.
It was pyrene — C16H10 where benzo[a]pyrene is C20H12 — and had been
mislabelled since the examples were written.

## Validation status of each scheme

| Scheme | Reference data | Result |
| --- | --- | --- |
| Cramer classes | 77 of 83 (Cramer Appendix 1) | 90.9% concordance |
| Cramer classes | 108 (reference implementation) | 96.2% concordance — but see *Comparing against Toxtree*: that oracle is a historical Toxtree run, which current Toxtree reproduces 84/108 |
| Verhaar scheme | 121 (Verhaar *et al.*) | 87.2% concordance |
| Genotoxicity alerts | 672 of 691 (ISSCAN Ames) | 90% sensitivity, 59% specificity |
| Sensitisation alerts | 208 of 209 (LLNA) | 92% sensitivity, 49% specificity |
| **Reactive group flags** | **none applicable** | **not a prediction** — renamed; agrees with Toxtree's SICRET/eye rulebases at about chance |
| **Non-genotoxic carcinogenicity** | **none found** | **unvalidated** — agrees with Toxtree on 50% of the compounds where either fires |

### Why this is called reactive group flags

There is no experimental irritation set to validate against. The Toxtree
repository's eye-irritation plugin ships rule classes and a 417-byte test
stub but no data, and every dataset in that repository over 20 KB is one
already used above or unrelated. No public machine-readable set surfaced
either: the BfR rulebase's 199-compound external validation set does not
appear to be published, and ChemSkin DB and the VEGA QMRF material are not
available as parseable structure files.

Two things follow, and both are stated here rather than papered over.

First, **these alerts are not an irritation rulebase and are no longer
labelled as one.** The published approach to local irritation is largely
*physicochemical* — the BfR decision support systems work from logP,
molecular weight, melting point, aqueous solubility and pKa cut-offs with
structural rules layered on top. This implementation is a SMARTS
reactivity list, which does not capture that at all.

The comparison settles what used to be an argument. Run over the same
compounds, this list agrees with Toxtree's SICRET skin rulebase on 50% of
the 129 cases SICRET will answer for, and with its eye rulebase on 60% of
37 — about what two unrelated classifiers would manage.

**And the BfR route is not open here**, which is why the panel was renamed
rather than rewritten. Toxtree itself answered "Unknown" for 80 of 209
compounds on SICRET and 172 of 209 on the eye rules, because the measured
properties those rules need were absent from the file. Melting point,
aqueous solubility and pKa cannot be computed reliably from a structure
and this app has no source for them, so a partial implementation would be
the same overstatement in new clothes. The panel is kept because "this
molecule contains a reactive fragment" is worth knowing; it is labelled
for what it is.

The verdict wording and colour changed with the name: a flag now reads
amber rather than red, since it reports a fragment rather than a hazard.

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
decided by that list — its 90.9% is entirely structural. On the food set 23 of
106 are; excluding them the remaining 83 score **95.2%**, against 96.2% overall.
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

### Pricing the questions that are not implemented

The Cramer tree here is 14 questions standing in for the original 33, and
Class II on Appendix 1 is the worst per-class figure in the project at 43%.
Whether finishing the tree would fix that used to be an argument. It is now a
measurement, because Toxtree emits a decision path for any compound, so the
question can be asked of structures the recorded paths never covered.

Run over Appendix 1, Toxtree gets six compounds right that this
implementation gets wrong. **Three of them turn on the same question, Q18**,
and its behaviour across both reference sets is unambiguous:

| | n | class |
| --- | --- | --- |
| Q18 answered **yes** | 5 | all Class II |
| Q18 answered **no** | 64 | all Class I |

It is the terminal Class I/II gate on the branch this implementation exits
early, which is precisely where the Class II accuracy is being lost. So the
answer is yes: the missing questions are worth implementing, and Q18 is the
one to do first.

What Q18 *asks* took one wrong hypothesis to find. The obvious guess, given
that answering yes yields Class II, was the food-list lookup that Q22 already
performs. The data refuses it: none of the five yes-answers is in the food
list, while 18 of the 64 no-answers are. That is the wrong way round for a
food test.

Reading the two groups together gives it. The yes-answers are 1,3-butanediol,
acetoin, BHT, 2-heptanone and allyl heptanoate — flavour and additive
chemicals. The no-answers on the same branch are, almost without exception,
fatty acids and amino acids: alanine, valine, leucine, isoleucine, aspartic
and glutamic acid, and essentially every straight-chain fatty acid from
hexanoic to tetracosanoic, saturated and unsaturated alike, plus ethanol,
propylene glycol and the simple diacids. Q18 separates **normal dietary and
metabolic constituents from everything else**, and it is a membership
question rather than a structural one — which is why no SMARTS pattern was
going to find it, and why implementing it will need the curated lists this
project already fetches rather than a new rule.

That is as far as this goes for now: the question is identified and priced,
not implemented. Adding a terminal gate to a branch is a change with real
regression risk, and it belongs in its own pass with the oracle in place.

## Licence

The code in this repository is **MIT licensed** (see `LICENSE`). The rules
are written from the published literature, so the implementation is
original work and carries no obligations from any other project.

**The reference data is a separate matter.** Everything the `tools/`
scripts fetch into `toxtree/data/` — the body and food constituent lists,
and the Cramer, Verhaar, LLNA and Ames validation sets — is GPL-licensed
where it is distributed. None of it is committed here, which is why those
paths are in `.gitignore` and why the app is built to run without them.
Fetching it for your own use is fine; redistributing it under this MIT
licence would not be, so do not add it to a fork.

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
