"""Regression tests for the rulebases.

Run with:  .venv/Scripts/python.exe -m pytest -q
or simply: .venv/Scripts/python.exe test_rules.py
"""
from __future__ import annotations

from toxtree import assess, parse_structure
from toxtree.alerts import validate_library
from toxtree.cramer import CLASS_I, CLASS_II, CLASS_III
from toxtree.examples import EXAMPLES

# (SMILES, expected Cramer class, expected Verhaar class number, genotox alert?)
CASES = [
    ("CCO", CLASS_I, 1, False),                       # ethanol
    ("OCCO", CLASS_I, 1, False),                      # ethylene glycol
    ("c1ccccc1", CLASS_I, 1, False),                  # benzene
    ("Oc1ccccc1", CLASS_I, 2, False),                 # phenol
    ("Nc1ccccc1", CLASS_III, 2, True),                # aniline
    ("C1OC1c1ccccc1", CLASS_III, 3, True),            # styrene oxide
    # Acrolein was pinned to Class III here on my own assumption. Validation
    # against Cramer et al. (1978) Appendix 1 corrected it: that paper places
    # trans-2-hexenal and cinnamaldehyde, the same functional class, in Class I,
    # so a simple alpha,beta-unsaturated aldehyde is Class I in this scheme.
    # The genotoxicity and Verhaar reactivity alerts still fire, which is where
    # the concern for acrolein actually shows up.
    ("C=CC=O", CLASS_I, 3, True),                     # acrolein
    ("O=CC=CCCC", CLASS_I, 3, True),                  # trans-2-hexenal (ref: 1)
    ("O=CC=Cc1ccccc1", CLASS_I, 3, True),             # cinnamaldehyde (ref: 1)
    ("c1cc(ccc1OC)C=CC", CLASS_III, 1, False),        # anethole (ref: 3)
    # Methyl acrylate was pinned here as genotoxicity-alerting, on SA10. It no
    # longer alerts, and the Ames data says that is right: the acrylates and
    # methacrylates in the reference set are all inactive. Its reactivity is
    # with protein, so it alerts for sensitisation instead, and Verhaar still
    # calls it Class 3.
    ("COC(=O)C=C", CLASS_III, 3, False),              # methyl acrylate
    ("OCCCl", CLASS_III, 3, True),                    # 2-chloroethanol (ref V3)
    # Sulphonate and sulphamate salts are Class I: not absorbed, excreted intact
    ("O=S(=O)([O-].[Na+])NC1CCCCC1", CLASS_I, 5, False),   # cyclamate (ref: 1)
    ("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]", CLASS_I, 5, False),  # lauryl sulfate
    ("CN(C)N=O", CLASS_III, 5, True),                 # NDMA
    ("CCOP(=S)(OCC)Oc1ccc(cc1)[N+](=O)[O-]", CLASS_III, 4, True),   # parathion
    ("CNC(=O)Oc1cccc2ccccc12", CLASS_III, 3, False),  # carbaryl (ref V3)
    ("CCNc1nc(Cl)nc(NC(C)C)n1", CLASS_III, 2, False), # atrazine (ref V2)
    ("OCC1OC(O)C(O)C(O)C1O", CLASS_I, 1, False),      # glucose
    ("NCC(=O)O", CLASS_I, 5, False),                  # glycine
    ("CC(=C)C1CCC(C)=CC1", CLASS_I, 1, False),        # limonene
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", CLASS_III, 5, False),  # caffeine
]


def test_smarts_library_compiles():
    assert validate_library() == []


def test_cramer_and_verhaar_classes():
    for smi, cramer_class, verhaar_no, genotox in CASES:
        rec = parse_structure(smi)
        assert rec.ok, f"{smi} failed to parse: {rec.error}"
        a = assess(rec)
        assert a.cramer.cramer_class == cramer_class, (
            f"{smi}: Cramer {a.cramer.cramer_class}, expected {cramer_class}")
        got = int(a.verhaar.verhaar_class.split()[1])
        assert got == verhaar_no, (
            f"{smi}: Verhaar class {got}, expected {verhaar_no}")
        assert bool(a.genotox) == genotox, (
            f"{smi}: genotox alerts {[h.code for h in a.genotox]}, "
            f"expected {'some' if genotox else 'none'}")


def test_every_example_parses_and_assesses():
    for name, smi, _ in EXAMPLES:
        rec = parse_structure(smi, name)
        assert rec.ok, f"{name} ({smi}) failed to parse: {rec.error}"
        a = assess(rec)
        assert a.cramer is not None and a.verhaar is not None
        assert a.to_row()["SMILES"]


def test_ttc_tiers():
    # a DNA-reactive alert collapses the TTC to the genotoxic-carcinogen tier
    a = assess(parse_structure("CN(C)N=O"))
    assert a.ttc_value_ug_day == 0.15
    # a carbamate without a genotox alert uses the cholinesterase tier
    a = assess(parse_structure("CNC(=O)Oc1cccc2ccccc12"))
    assert a.ttc_value_ug_day == 18.0
    # an inert Class I substance uses the Cramer tier
    a = assess(parse_structure("CCO"))
    assert a.ttc_value_ug_day == 1800.0
    # metals fall outside the TTC domain entirely
    a = assess(parse_structure("[Hg]"))
    assert a.ttc_value_ug_day is None


def test_primary_metabolite_question():
    """Cramer's Q1 is narrow: metabolic building blocks, not "in the body".

    An earlier version used the full 439-structure body-constituent file and
    needed an explicit guard to stop biotin and riboflavin reaching Class I.
    With the right membership no guard is needed -- cysteine and methionine
    are in the set, biotin and ascorbic acid simply are not -- and the guard
    was removed.
    """
    from toxtree import reflists

    # The app must work whether or not the optional lists are installed.
    assert isinstance(reflists.status(), dict)

    def cls(smiles):
        return assess(parse_structure(smiles), {"cramer"}).cramer.cramer_class

    for smiles in ("NCC(=O)O",                      # glycine
                   "N[C@@H](CS)C(O)=O",             # cysteine (sulphur)
                   "CSCC[C@H](N)C(O)=O",            # methionine (sulphur)
                   "[H][C@@]1(CCCN1)C(O)=O",        # proline (2nd amine)
                   "N[C@@H](CC1=CNC=N1)C(O)=O"):    # histidine (tautomer)
        assert cls(smiles) == CLASS_I, smiles

    # In the body, but not metabolic building blocks: these must fall through
    # to the structural questions rather than short-circuit to Class I.
    for smiles in ("OC(=O)CCCCC1SCC2NC(=O)NC12",            # biotin
                   "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",        # ascorbic acid
                   "OC(=O)CC(O)(CC(O)=O)C(O)=O",            # citric acid
                   "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO"):  # sorbitol
        assert cls(smiles) != CLASS_I, smiles


def test_carbohydrates_and_oxygen_density():
    """Cyclic sugars are Class I; open-chain polyols and polyacids are not."""
    def cls(smiles):
        return assess(parse_structure(smiles), {"cramer"}).cramer.cramer_class

    # ring form -> Class I
    assert cls("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O[C@@H]2O[C@H](CO)"
               "[C@H](O)[C@H](O)[C@H]2O") == CLASS_I          # lactose
    assert cls("OC[C@H]1O[C@H](O[C@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)"
               "[C@H](O)[C@@H](O)[C@@H]1O") == CLASS_I        # sucrose

    # one or two oxygen functions still read as readily metabolised
    assert cls("CCO") == CLASS_I                              # ethanol
    assert cls("OCC(O)C") == CLASS_I                          # propylene glycol
    assert cls("CCCCCCCCCCCCCCCCCCCCCCCC(O)=O") == CLASS_I    # lignoceric acid
    # three or more, or an alpha-dicarbonyl, do not
    assert cls("OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO") != CLASS_I  # mannitol
    assert cls("CC(=O)C(C)=O") != CLASS_I                     # diacetyl
    assert cls("OC(=O)C(O)=O") != CLASS_I                     # oxalic acid


def test_food_question_is_gated():
    from toxtree import reflists
    if reflists.FOOD:
        # Food-list membership is context only: it must not change the class.
        # Caffeine is on the list and is Class III; ethanol is on it and is
        # Class I. Any rule keying off membership alone breaks one of them.
        caffeine = assess(parse_structure("Cn1cnc2c1c(=O)n(C)c(=O)n2C"))
        assert caffeine.cramer.cramer_class == CLASS_III
        assert caffeine.cramer.food_component
        ethanol = assess(parse_structure("CCO"))
        assert ethanol.cramer.cramer_class == CLASS_I
        assert ethanol.cramer.food_component


def test_class_two_discrimination():
    """The food question (Q22) and the acylated-heteroatom question (Q11).

    These were recovered from ~190 recorded decision paths. Q22 asks whether a
    substance is a common component of food and answers Class II, but it is
    only reached when the ring heteroatoms are unacylated. That single gate is
    what separates maltol from caffeine -- both food-listed heterocycles, but
    Class II and Class III respectively -- and it is the whole reason Class II
    is reproducible at all.
    """
    from toxtree import reflists
    if not reflists.FOOD:
        return  # reference lists not installed

    def cls(smiles):
        return assess(parse_structure(smiles), {"cramer"}).cramer.cramer_class

    # Acylated ring heteroatom -> bypasses the food question -> Class III
    assert cls("Cn1cnc2c1c(=O)n(C)c(=O)n2C") == CLASS_III          # caffeine
    assert cls("CN1C=NC2=C1C(=O)NC(=O)N2C") == CLASS_III           # theobromine
    assert cls("CC1=CC2=C(C=C1C)N(C[C@@H](O)[C@@H](O)[C@@H](O)CO)"
               "C3=NC(=O)NC(=O)C3=N2") == CLASS_III                # riboflavin

    # No acylated heteroatom -> reaches the food question -> Class II
    assert cls("CC1=C(O)C(=O)C=CO1") == CLASS_II                   # maltol
    assert cls("CCC1=C(O)C(=O)C=CO1") == CLASS_II                  # ethyl maltol
    assert cls("OC(=O)C1=CC=CN=C1") == CLASS_II                    # niacin
    assert cls("CC1=CN=CC=N1") == CLASS_II                         # methylpyrazine

    # Compounds that exit as Class I earlier must never reach the food
    # question, even though they are on the food list.
    assert cls("COc1cc(C=O)ccc1O") == CLASS_I                      # vanillin
    assert cls("CCO") == CLASS_I                                   # ethanol


def test_molblock_with_a_blank_title_line_parses():
    """A molblock's title line is usually blank, and that must still parse.

    parse_structure() used to strip its input before deciding what it was
    looking at. A molblock's first line is the title and is very often empty
    -- RDKit's own MolToMolBlock writes it that way -- so stripping deleted
    it, shifted every following line, and the block failed to parse. Pasting
    any RDKit-generated molblock into the app hit this.
    """
    from rdkit import Chem

    block = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
    assert block.split("\n")[0].strip() == "", "test assumes a blank title"

    for label, text, ident in (
            ("blank title", block, None),
            ("named title", "ethanol\n" + block.split("\n", 1)[1], None),
            ("explicit identifier", block, "my ethanol"),
            ("trailing whitespace", block + "\n\n  ", None)):
        rec = parse_structure(text, ident)
        assert rec.ok, f"{label}: {rec.error}"
        assert rec.smiles == "CCO", label

    # the title names the record when there is one, and never the whole block
    assert parse_structure(block).identifier == "molblock"
    assert parse_structure("ethanol\n" + block.split("\n", 1)[1]
                           ).identifier == "ethanol"

    # the other input forms are unaffected
    assert parse_structure("CCO").ok
    assert parse_structure("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3").ok
    assert not parse_structure("   ").ok
    assert not parse_structure("not a molecule").ok


def test_question_order_and_parsing_edge_cases():
    """Ordering and input quirks that silently produced wrong classes."""
    def cls(smiles):
        return assess(parse_structure(smiles), {"cramer"}).cramer.cramer_class

    # Water is a normal body constituent and Class I. It only reaches that
    # answer because the body-constituent question precedes the
    # organic/metal applicability check -- reversed, water is "inorganic".
    assert cls("O") == CLASS_I
    # ... and it must survive being written with explicit hydrogens, which
    # would otherwise be counted as heavy atoms throughout.
    assert cls("[H]O[H]") == CLASS_I

    # Open-chain sugars are carbohydrates; the fully reduced polyols are not.
    assert cls("OCC(=O)[C@@H](O)[C@H](O)[C@H](O)CO") == CLASS_I   # fructose
    assert cls("OC(=O)[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)C=O") == CLASS_I

    # A single ionised head group on an alkyl chain is enough; the density
    # requirement applies only to large aromatic dyes.
    assert cls("O=S(=O)([O-])OCCOCCOCCOCCCCCCCCCCCC.[Na+]") == CLASS_I

    # Sulphur compounds are routed away from the food question, so a listed
    # food component that reaches Class III through sulphur stays there.
    from toxtree import reflists
    if reflists.FOOD:
        assert cls("S1C=CC=C1") == CLASS_III                      # thiophene


def test_verhaar_rules():
    """Verhaar behaviours fixed by validation against the 121-compound set."""
    def vh(smiles):
        return int(assess(parse_structure(smiles),
                          {"verhaar"}).verhaar.verhaar_class.split()[1])

    # Two bugs: an alert labelled "two or more nitro groups" whose pattern
    # matched a single one, and a trifluoromethyl group being read as an
    # activated halide. Both made polar narcotics look reactive.
    assert vh("O=[N+]([O-])c1ccc(O)cc1") == 2              # 4-nitrophenol
    assert vh("O=[N+]([O-])c1ccc(O)cc1C(F)(F)F") == 2      # and its CF3 analogue
    # genuine dinitro compounds are still reactive
    assert vh("O=[N+]([O-])c1ccc(O)c([N+](=O)[O-])c1") == 3

    # Amides, carbamates and aromatic esters hydrolyse to acylating species:
    # unanimous in the reference data, and previously answered "class 5".
    assert vh("CN(C)C=O") == 3                              # dimethylformamide
    assert vh("CCCCOC(=O)c1ccccc1C(=O)OCCCC") == 3          # dibutyl phthalate
    assert vh("CNC(=O)Oc1cccc2ccccc12") == 3                # carbaryl

    # Reactive aliphatics that were being called inert
    assert vh("C=CCO") == 3                                 # allyl alcohol
    assert vh("OCCCl") == 3                                 # 2-chloroethanol

    # Organophosphates remain specifically acting; triazines are not
    assert vh("CCOP(=S)(OCC)Oc1ccc(cc1)[N+](=O)[O-]") == 4  # parathion
    assert vh("CCNc1nc(Cl)nc(NC(C)C)n1") == 2               # atrazine

    # Baseline narcosis is untouched
    assert vh("CCO") == 1 and vh("c1ccccc1") == 1


def test_alert_rulebases():
    """Alert behaviours fixed by validating against LLNA and Ames data."""
    def codes(smiles, scheme):
        return {h.code for h in getattr(
            assess(parse_structure(smiles), {scheme}), scheme)}

    # --- skin sensitisation (LLNA) -------------------------------------
    # Gaps that cost half the sensitivity of the original alert set.
    assert "SS3b" in codes("CCCCC1=NC(C)(C)C(=O)O1", "sensitisation")  # azlactone
    assert "SS3b" in codes("O=C1CCc2ccccc2O1", "sensitisation")        # lactone
    assert "SS4" in codes("CCCCCCCCCCCCCCCCCl", "sensitisation")       # alkyl chloride
    assert "SS2b" in codes("CC(=O)C(C)=O", "sensitisation")            # 1,2-diketone
    assert "SS2b" in codes("CC(=O)CC(=O)c1ccccc1", "sensitisation")    # 1,3-diketone
    assert "SS6" in codes("C=CCc1ccc(O)c(OC)c1", "sensitisation")      # eugenol
    assert "SS11" in codes("Cn1sccc1=O", "sensitisation")              # methylisothiazolinone

    # A quinone methide needs benzylic hydrogens, so a carboxyl or ester
    # carbon on the ring must not fire SS6 -- that was flagging every paraben
    # and salicylate as a sensitiser.
    for non_sensitiser in ("CCCOC(=O)c1ccc(O)cc1",     # propyl paraben
                           "COC(=O)c1ccccc1O",          # methyl salicylate
                           "O=C(O)c1ccc(O)cc1"):        # 4-hydroxybenzoic acid
        assert "SS6" not in codes(non_sensitiser, "sensitisation"), non_sensitiser
    # saccharin's sulphonyl sulphur is not a reactive isothiazolinone
    assert "SS11" not in codes("O=C1NS(=O)(=O)c2ccccc21", "sensitisation")

    # --- genotoxicity (Ames) -------------------------------------------
    # Fully esterified sulfates and phosphates alkylate; the salts do not.
    assert "SA2b" in codes("CCOS(=O)(=O)OCC", "genotox")        # diethyl sulfate
    assert "SA2b" in codes("COP(=O)(OC)OC", "genotox")          # trimethyl phosphate
    assert not codes("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]", "genotox")

    # This used to assert that acrylonitrile fires SA10, the alpha,beta-
    # unsaturated acceptor. SA10 has since been removed: measured over the
    # Ames set it alone flagged 14 inactives against 3 actives, because a
    # Michael acceptor is protein-reactive rather than DNA-reactive. The same
    # pattern earns 96% precision as a sensitisation alert, so acrylonitrile
    # is now expected to alert there and not here. Giving it up as a
    # genotoxicity miss was the deliberate price.
    assert not codes("C=CC#N", "genotox")                       # acrylonitrile
    assert "SS1" in codes("C=CC#N", "sensitisation")
    # ...but the pyrrolizidine alkaloids SA10 was catching incidentally are
    # kept, by the alert that describes what actually makes them genotoxic.
    assert "SA39" in codes(
        r"C/C=C(\C)C(=O)OC1CCN2CC=C(COC(=O)C(C)(O)C(C)(O)C)C12",
        "genotox")                                              # lasiocarpine


def test_pro_hapten_sensitisers_and_the_esters_they_must_not_catch():
    """The pro-hapten alerts, and the boundary that keeps them affordable.

    Isoeugenol is the regression that matters: the phenol alert claimed to
    cover the eugenol/isoeugenol family but a SMARTS precedence slip ("," binds
    tighter than ";") made its unsaturated branch unreachable, so it caught
    eugenol and missed isoeugenol for as long as the comment said otherwise.

    The esters are the other boundary. An ester of a phenol acylates lysine and
    sensitises; an alkyl benzoate does not, and a pattern loose enough to take
    in both flags every paraben. These four inactives are what that costs.
    """
    def codes(smiles):
        rec = parse_structure(smiles)
        assert rec.ok, smiles
        return {h.code for h in assess(rec, {"sensitisation"}).sensitisation}

    assert "SS6" in codes(r"C/C=C\c1ccc(O)c(OC)c1")      # isoeugenol
    assert "SS6" in codes("C=CCc1ccc(O)c(OC)c1")          # eugenol
    assert "SS13" in codes("C/C=C/c1ccc(OC)cc1")          # trans-anethole
    assert "SS14" in codes(r"OC/C=C\c1ccccc1")            # cinnamic alcohol
    assert "SS16" in codes("CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)"
                           "C(=O)N2[C@H]1C(=O)O")         # penicillin G
    assert "SS17" in codes("NCCN")                        # ethylenediamine
    assert "SS18" in codes("CN(N=O)C(N)=O")               # N-methyl-N-nitrosourea
    assert "SS21" in codes("O=C(Oc1ccccc1)c1ccccc1")      # phenyl benzoate

    # ...and the inactives that a broader ester pattern would have flagged
    for ester in ("O=C(OC)c1ccc(O)cc1",                   # methyl paraben
                  "O=C(OCCC)c1ccc(O)cc1",                 # propyl paraben
                  "O=C(OC)c1ccccc1O",                     # methyl salicylate
                  "O=C(OCC)c1ccccc1C(=O)OCC"):            # diethyl phthalate
        assert "SS21" not in codes(ester), ester

    # an alpha-amino acid is not a polyamine and must not fire SS17
    assert "SS17" not in codes("NCCCCC(N)C(=O)O")         # lysine


def test_peroxisome_proliferator_ester_needs_a_hydrolysable_link():
    """NG9, and the boundary that makes it a mechanism rather than a shape.

    The alert exists because the 2-ethylhexyl esters are hydrolysed to
    2-ethylhexanol and oxidised to 2-ethylhexanoic acid, a PPAR-alpha agonist.
    That requires a linkage that hydrolyses, so the alcohol on its own and the
    corresponding ether must not fire: neither releases the metabolite. DEHP
    fires NG7 and NG9 both, which is correct -- they describe one mechanism
    from two directions.
    """
    def codes(smiles):
        rec = parse_structure(smiles)
        assert rec.ok, smiles
        return {h.code for h in assess(rec, {"nongenotox"}).nongenotox}

    assert "NG9" in codes("CCCCC(CC)COC(=O)CCCCC(=O)OCC(CC)CCCC")   # DEHA
    assert "NG9" in codes("CCCCC(CC)COP(=O)(OCC(CC)CCCC)OCC(CC)CCCC")
    assert "NG9" in codes("CCCCC(CC)COS(=O)(=O)[O-].[Na+]")
    assert {"NG7", "NG9"} <= codes("CCCCC(CC)COC(=O)c1ccccc1C(=O)OCC(CC)CCCC")

    assert "NG9" not in codes("CCCCC(CC)CO")                        # the alcohol
    assert "NG9" not in codes("CCCCC(CC)COCC(CC)CCCC")              # the ether
    assert "NG9" not in codes("CCOC(C)=O")                          # ethyl acetate
    # a phthalate of an unbranched alcohol is NG7's business, not NG9's
    assert codes("CCCCOC(=O)c1ccccc1C(=O)OCCCC") == {"NG7"}


def test_non_dna_reactive_mechanisms_are_not_genotoxicity_alerts():
    """Thyroid and AhR mechanisms belong in the non-genotoxic rulebase.

    Both of these were genotoxicity alerts until they were measured against
    the ISSCAN Ames set, where the thiourea pattern scored 0 true positives
    against 10 false and the dioxin pattern 1 against 5. Neither mechanism is
    DNA-reactive -- a thiocarbamide inhibits thyroid peroxidase and an AhR
    ligand promotes tumours -- so a firing here would be the old mistake
    coming back rather than a new finding.
    """
    def codes(smiles, scheme):
        a = assess(parse_structure(smiles), {scheme})
        return {h.code for h in getattr(a, scheme)}

    for smiles in ("NC(N)=S",                    # thiourea
                   "CC(N)=S",                    # thioacetamide
                   "Clc1ccc(-c2ccccc2)cc1"):     # 4-chlorobiphenyl
        assert not codes(smiles, "genotox"), smiles
        assert codes(smiles, "nongenotox"), smiles

    # the genotoxicity set must still catch what is genuinely DNA-reactive
    assert codes("CN(C)N=O", "genotox")          # N-nitrosodimethylamine
    assert codes("C1CO1", "genotox")             # ethylene oxide


def test_irritation_alerts_are_not_indiscriminate():
    """The irritation set has no reference data, so this guards its shape.

    No experimental irritation set is available (the Toxtree repository has
    none, and no public machine-readable one surfaced), so these alerts are
    unvalidated. What can be checked without labels is that they discriminate
    at all, and that they do not fire on compounds whose chemistry rules it
    out. An alpha-amino acid is zwitterionic at skin pH, so neither its
    carboxyl nor its amine behaves as a free acid or base.
    """
    def codes(smiles):
        return {h.code for h in
                assess(parse_structure(smiles), {"irritation"}).irritation}

    for amino_acid in ("NCC(=O)O",              # glycine
                       "CC(N)C(=O)O",           # alanine
                       "CC(C)CC(N)C(=O)O"):     # leucine
        assert not codes(amino_acid) & {"IR1", "IR4"}, amino_acid

    # but a free acid, a free base and a phenol must still fire
    assert "IR1" in codes("CC(=O)O")            # acetic acid
    assert "IR4" in codes("CCN")                # ethylamine
    assert "IR7" in codes("Oc1ccccc1")          # phenol
    assert "IR3" in codes("CC(=O)Cl")           # acetyl chloride

    # and the set must not simply flag everything
    for benign in ("CCCC", "OCC1OC(O)C(O)C(O)C1O", "c1ccccc1"):
        assert not codes(benign), benign


def test_read_across_support():
    """Evidential support: nearest reference compounds and the index."""
    from toxtree import readacross

    assert isinstance(readacross.status(), dict)
    if not readacross.available_endpoints():
        return  # reference data not installed in this checkout

    # A compound in the reference set, predicted as the data measured it,
    # should come out well supported.
    rec = parse_structure("Cn1cnc2c1c(=O)n(C)c(=O)n2C")        # caffeine
    strong = readacross.support(rec.mol, "cramer", "Class III")
    assert strong.available and strong.neighbours
    assert strong.similarity > 0.95          # it is in the set
    assert strong.index > 0.70 and strong.band == "well supported"

    # Something unlike anything in the reference data must say so, whatever
    # the rules decided -- this is the whole point of the panel.
    rec = parse_structure("C1C2C3C4C1C5C2C3C45C(=O)N6CCOCC6")
    weak = readacross.support(rec.mol, "cramer", "Class III")
    assert weak.similarity < 0.5
    assert weak.index < 0.45 and weak.band == "little support"
    assert any("similar" in n for n in weak.notes)

    # Concordance must respond to the verdict it is given, not ignore it.
    rec = parse_structure("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    agreeing = readacross.support(rec.mol, "cramer", "Class III")
    disagreeing = readacross.support(rec.mol, "cramer", "Class I")
    assert agreeing.concordance > disagreeing.concordance

    # An endpoint with no data installed must degrade, not raise.
    missing = readacross.support(rec.mol, "not_an_endpoint", "whatever")
    assert not missing.available and missing.band == "no reference data"


def test_bad_input_is_reported_not_raised():
    rec = parse_structure("not a molecule at all")
    assert not rec.ok and rec.error


if __name__ == "__main__":
    passed = 0
    for fn in [v for k, v in sorted(globals().items()) if k.startswith("test_")]:
        fn()
        print(f"  ok  {fn.__name__}")
        passed += 1
    print(f"\n{passed} test group(s) passed.")
