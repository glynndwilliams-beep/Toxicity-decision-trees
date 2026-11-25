"""Cramer classification of structural classes of oral toxicity.

A Python re-implementation of the Cramer, Ford & Hall (1978) decision tree,
written from the published rule descriptions. Each molecule walks an ordered
sequence of questions; the first question that resolves decides the class and
the full path is returned so the answer can be audited.

Classes
-------
I    Substances of simple chemical structure with efficient metabolism,
     suggesting a low order of oral toxicity.
II   Substances that are intermediate: less innocuous than Class I but
     without the positive indication of toxicity of Class III.
III  Substances with a structure that permits no strong initial presumption
     of safety, or that suggests significant toxicity.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem

from . import reflists
from .alerts import CRAMER_HIGH_TOX, screen
from .chem import CRAMER_ALLOWED, count, has

CLASS_I = "Class I (low)"
CLASS_II = "Class II (intermediate)"
CLASS_III = "Class III (high)"

TTC_UG_PER_PERSON_PER_DAY = {CLASS_I: 1800.0, CLASS_II: 540.0, CLASS_III: 90.0}


@dataclass
class CramerResult:
    cramer_class: str
    trace: list[str] = field(default_factory=list)
    alerts: list[str] = field(default_factory=list)
    out_of_domain: bool = False
    food_component: str | None = None

    @property
    def ttc_ug_day(self) -> float:
        return TTC_UG_PER_PERSON_PER_DAY[self.cramer_class]



# ---------------------------------------------------------------------------
# Structural question patterns
# ---------------------------------------------------------------------------
P_BENZENE = "c1ccccc1"
P_HETEROAROMATIC = "[n,o,s;a]"
P_HETEROCYCLE_SAT = "[NX3,OX2,SX2;R]"
P_LACTONE = "[#6;R][CX3;R](=[OX1])[OX2;R]"
P_CYCLIC_DIESTER = "[OX2;R][CX3;R](=[OX1])[#6;R][CX3;R](=[OX1])[OX2;R]"
P_SULFONATE = "[SX4](=[OX1])(=[OX1])[OX2,OX1-]"
P_SULFAMATE = "[NX3][SX4](=[OX1])(=[OX1])[OX2,OX1-]"
# Salt forms only: the acidic oxygen must be free (O- or OH), not esterified
# to carbon. A fully esterified sulphonate is an alkylating agent, not an
# unabsorbed salt, so it must not reach this question as a Class I.
P_SULFONATE_SALT = "[#6][SX4](=[OX1])(=[OX1])[OX2H,OX1-]"
P_SULFAMATE_SALT = "[NX3][SX4](=[OX1])(=[OX1])[OX2H,OX1-]"
P_SULFATE_SALT = "[#6][OX2][SX4](=[OX1])(=[OX1])[OX2H,OX1-]"
P_PRIMARY_AROMATIC_AMINE = "[c][NX3;H2]"
# Recovered Q11: a ring heteroatom bearing an acyl group -- a cyclic amide,
# urea, imide or sulfonamide, or an acyl hydrazide. This is what separates
# caffeine, theobromine, biotin, saccharin, riboflavin and folic acid (which
# bypass the food question and end at Class III) from maltol, ascorbic acid,
# the pyrazines, niacin and the tocopherols (which reach it and end at
# Class II). It splits all 13 observed Y answers from all 14 N answers.
# Element-based carbon (#6X3, not CX3) is essential: RDKit aromatises the
# xanthine and pteridine rings, so the lactam carbonyls of caffeine,
# theobromine, folic acid and riboflavin sit on aromatic carbons that an
# aliphatic-carbon pattern silently misses.
P_ACYL_HETEROATOM = (
    "[$([#7;R][#6X3]=[OX1]),$([#7;R][#6X3]=[SX1]),"
    "$([#7;R][SX4](=[OX1])=[OX1]),$([#7][#7X3][#6X3]=[OX1]),"
    "$([#7;R][#6X3](=[OX1])[#7])]"
)
P_PHOSPHORUS_ACID = "[PX4](=[OX1])[OX2H,OX1-]"
P_ALPHA_DICARBONYL = "[CX3](=[OX1])[CX3]=[OX1]"
P_ALIPHATIC_HALIDE = "[CX4;!$(C(F)(F)F)][F,Cl,Br,I]"
P_SULFUR_ANY = "[#16]"
P_ALCOHOL = "[CX4][OX2H]"
P_ALDEHYDE = "[CX3H1]=[OX1]"
P_KETONE = "[#6][CX3](=[OX1])[#6]"
P_ACID = "[CX3](=[OX1])[OX2H,OX1-]"
P_ESTER = "[CX3](=[OX1])[OX2][#6]"
P_ETHER = "[#6][OX2][#6]"
P_AMIDE = "[NX3][CX3](=[OX1])[#6]"
# Substituents that Cramer regards as readily metabolised. The alkenyl entry
# requires an oxidised terminus (-CH=CH-CHO, -CH=CH-CH2OH, -CH=CH-COOH): a
# cinnamyl side chain is oxidised to cinnamic acid and excreted, which is why
# cinnamaldehyde is Class I, whereas a bare propenyl group as in anethole can
# be epoxidised and is not innocuous.
P_HARMLESS_AR_SUBST = (
    "[$([c][OX2H]),$([c][OX2][CX4]),$([c][CX4]),$([c][CX3](=O)[OX2H]),"
    "$([c][CX3](=O)[OX2][#6]),$([c][CX3H]=O),$([c][OX2][CX3]=O),"
    "$([c][CX3]=[CX3][$([CX3H]=[OX1]),$([CX4][OX2H]),$([CX3](=[OX1])[OX2H])])]"
)
P_TERPENE_SKELETON = "[CX4]([CX4])([CX4])[CX4]"
P_POLYCYCLIC_NONAROMATIC = "[#6;R2;!a]"
P_FURAN = "c1ccoc1"
P_PYRIDINE = "c1ccncc1"
# Q6 calls a lactone "simple" below this size, and sends the larger ones to
# Class III. The number is a choice, not the paper's, and the sweep in the
# README's Q6 section shows it discriminates nothing: every value from 0 to 11
# scores identically, every value from 12 up scores identically and one worse,
# because the only two lactones in either reference set that the cut could
# separate -- dehydroacetic acid (Class III) and ascorbic acid (Class II) --
# have exactly 12 heavy atoms each. 11 is used because it gets dehydroacetic
# acid right and costs nothing, not because 11 is calibrated; ascorbic acid
# then reaches Class II through the food question rather than through this
# rule. Three lactones in 183 compounds is the whole evidence base.
Q6_SIMPLE_LACTONE_MAX_HEAVY = 11

# Q8 calls a heterocycle on the common-in-food list "simple" below this size.
# Swept in the README's Q8 section.
Q8_SIMPLE_HETEROCYCLE_MAX_HEAVY = 12

# Q12 reads a chain as readily metabolised up to this many oxygen functions,
# and treats it as simple while it has no more rings than this. Both are swept
# in the README's Q12 section, and both sit at an optimum the data constrains
# from each side.
Q12_MAX_OXYGEN_FUNCTIONS = 2
Q12_MAX_RINGS = 1

P_COMMON_FOOD_HETEROCYCLE = (
    "[$(c1ccoc1),$(c1ccsc1),$(c1ccncc1),$(C1CCOC1),$(C1CCCO1),$(O=C1CCCO1)]"
)




def _is_carbohydrate(mol: Chem.Mol) -> bool:
    """True for a cyclic sugar: a pyranose or furanose ring, or a glycoside.

    Cramer places the sugars in Class I. The ring matters: lactose, maltose,
    sucrose and raffinose are Class I, while the open-chain sugar alcohols
    mannitol and sorbitol are Class II, so this tests for the ring form and
    not merely for a run of hydroxyls. A ring carbonyl disqualifies, which
    keeps ascorbic acid (a lactone, Class II) out.
    """
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(a.GetIsAromatic() for a in atoms):
            continue
        if sum(1 for a in atoms if a.GetSymbol() == "O") != 1:
            continue
        if any(a.GetSymbol() not in ("C", "O") for a in atoms):
            continue
        if any(b.GetBondTypeAsDouble() == 2 for a in atoms for b in a.GetBonds()):
            continue
        oxygens = 0
        for a in atoms:
            if a.GetSymbol() != "C":
                continue
            for nb in a.GetNeighbors():
                if nb.GetIdx() not in ring and nb.GetSymbol() == "O":
                    oxygens += 1
        if oxygens >= 2:
            return True

    # Open-chain form: a sugar is a polyhydroxy aldehyde or ketone. Fructose
    # and galacturonic acid are drawn this way and are Class I, while the
    # fully reduced sugar alcohols mannitol and sorbitol -- same hydroxyl
    # count, no carbonyl -- are Class II. The carbonyl is the distinction.
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    if elements <= {"C", "H", "O"}:
        carbonyl = has(mol, P_ALDEHYDE) or has(mol, P_KETONE)
        if carbonyl and count(mol, P_ALCOHOL) >= 3:
            return True
    return False


# --- Q18 ---------------------------------------------------------------
# Cramer's Q18 is a nine-branch list, not a single test, which is why it
# resisted every attempt to recover it from decision paths alone. The wording
# is reproduced in the JRC's analysis of the scheme (EUR 24898 EN, 2011);
# answering yes gives Class II where the branch would otherwise end at
# Class I. Branches (a) to (h) are implemented here.
#
# Branch (i), "a substance in which the functional groups are all sterically
# hindered", is undefined in the paper, and Phillips et al. (1987) deleted
# the phrase as ambiguous when they revised the scheme. It is implemented
# here under a definition this file had to choose rather than read off the
# paper: see _q18_all_hindered. The definition was simulated against both
# reference sets before being written in, and kept only because it moved
# one compound and that compound moved the right way.
Q18_BRANCHES = (
    # A diketone, so both carbonyls must carry carbon on the outside. The
    # looser [CX3](=O)[CX3]=O also matches oxalic acid and the other vicinal
    # diacids: they reach Class II anyway through the alpha-dicarbonyl rule,
    # but naming branch (a) as the reason would be wrong, and the decision
    # path is shown to the user.
    ("a", "vicinal diketone, or a ketone attached to a terminal vinyl group",
     "[$([#6][CX3](=[OX1])[CX3](=[OX1])[#6]),"
     "$([CH2]=[CH][CX3](=[OX1])[#6]),"
     "$([CH2]=[CH][CX4]([OX2])([OX2])[#6])]"),
    ("b", "secondary alcohol, or its ester, attached to a terminal vinyl group",
     "[$([CH2]=[CH][CX4H1]([#6])[OX2H]),"
     "$([CH2]=[CH][CX4H1]([#6])[OX2][CX3]=[OX1])]"),
    ("c", "allyl alcohol or its acetal, ketal or ester derivative",
     "[$([CH2]=[CH][CH2][OX2H]),"
     "$([CH2]=[CH][CH2][OX2][CX3]=[OX1]),"
     "$([CH2]=[CH][CH2][OX2][CX4][OX2])]"),
    ("d", "allyl mercaptan, sulphide, thioester or amine",
     "[CH2]=[CH][CH2][SX2,NX3]"),
    ("e", "acrolein, methacrolein or their acetals",
     "[$([CH2]=[CX3][CX3H1]=[OX1]),"
     "$([CH2]=[CX3][CX4H1]([OX2])[OX2])]"),
    # The acid only. Cramer assigns acrylic acid to Class II but leaves the
    # acrylate esters in Class I; the JRC report flags that as an
    # inconsistency in the scheme, but it is the scheme's, not this file's.
    ("f", "acrylic or methacrylic acid",
     "[CH2]=[CX3][CX3](=[OX1])[OX2H,OX1-]"),
    ("g", "acetylenic compound", "[CX2]#[CX2]"),
)


def _q18_large_ketone(mol: Chem.Mol) -> bool:
    """Branch (h): an acyclic aliphatic ketone with a long chain on one side.

    "Four or more carbons on either side of the keto group" is read as
    written: four on *one* side suffices. This reverses the strict reading
    (four on *both* sides) that this file shipped first, and the reversal was
    forced by measurement, not by rereading the sentence.

    The strict reading was justified here on the grounds that the loose one
    is where Bhatia et al. report Toxtree's Q18 over-calling, and that it
    would otherwise catch "every methyl ketone". The second half of that is
    simply false, and measuring it said so: the four-carbon requirement still
    excludes acetone, 2-butanone, 2-pentanone and 3-methyl-2-butanone, and
    the loose reading adds only the methyl ketones from 2-hexanone upward --
    which is precisely the class of large, poorly metabolised ketones the
    question is about. Across both reference sets it fires on exactly one
    compound that the strict reading missed, 2-heptanone, whose reference
    class is II; it breaks nothing in either set. The strict reading fires on
    no compound in either set at all, so it was carrying its justification
    and no evidence.

    Bhatia et al.'s warning is still live, and this rule is bounded by the
    rest of the question rather than by the carbon count: "no other
    functional groups" is taken as carbon, hydrogen and oxygen only, acyclic,
    a single ketone, and no oxygen function beyond the keto group itself and
    at most one hydroxyl (the "ketoalcohol" the question allows). Dropping
    that restriction as well was measured and rejected: it breaks raspberry
    ketone, 4-(4-hydroxyphenyl)butan-2-one, a reference Class I that the
    acyclic condition is what keeps out, for the one compound it gains.

    What the loose reading does widen, and what no reference compound tests,
    is the branched ketoalcohols: diacetone alcohol reaches four carbons on
    one side through two branch methyls rather than through chain length,
    and it answers yes here where the strict reading said no. Neither
    reference set contains it, so that verdict is unmeasured rather than
    confirmed, and it is recorded in the README as the loose reading's known
    exposure.
    """
    if any(a.GetSymbol() not in ("C", "H", "O") for a in mol.GetAtoms()):
        return False
    if mol.GetRingInfo().NumRings():
        return False
    hydroxyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]")))
    if hydroxyls > 1:
        return False
    ketones = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[#6][CX3](=[OX1])[#6]"))
    if len(ketones) != 1:
        return False
    total_o = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX1,OX2]")))
    if total_o != 1 + hydroxyls:
        return False

    carbonyl = ketones[0][1]
    for side in (ketones[0][0], ketones[0][3]):
        seen, stack = {carbonyl}, [side]
        carbons = 0
        while stack:
            idx = stack.pop()
            if idx in seen:
                continue
            seen.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "C":
                carbons += 1
            stack.extend(n.GetIdx() for n in atom.GetNeighbors()
                         if n.GetIdx() not in seen)
        # One side reaching four carbons is enough; see the docstring for
        # why the "both sides" reading was measured and dropped.
        if carbons >= 4:
            return True
    return False


def _is_quaternary_carbon(atom: Chem.Mol) -> bool:
    """A fully alkylated sp3 carbon: the central carbon of a tert-butyl group.

    Heteroatom substitution disqualifies, so the quaternary carbon of a
    hemiketal or an amino acid does not count as a blocking group.
    """
    if atom.GetSymbol() != "C" or atom.GetIsAromatic():
        return False
    if atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False
    neighbours = list(atom.GetNeighbors())
    if any(n.GetSymbol() not in ("C", "H") for n in neighbours):
        return False
    return sum(1 for n in neighbours if n.GetSymbol() == "C") == 4


# Groups that count as "functional" for branch (i), each with the index of
# the match atom the group hangs off. The list is deliberately generous: the
# question says *all* groups must be hindered, so a group this list misses
# can only make the branch fire too often.
_Q18_HINDERED_GROUPS = (
    ("[OX2H]", 0),                     # alcohol or phenol
    ("[CX3](=[OX1])[OX2H,OX1-]", 0),   # carboxylic acid
    ("[CX3H1]=[OX1]", 0),              # aldehyde
    ("[#6][CX3](=[OX1])[#6]", 1),      # ketone
    ("[CX3](=[OX1])[OX2][#6]", 0),     # ester
    ("[NX3;!$(N[CX3]=[OX1])]", 0),     # amine
    ("[NX3][CX3](=[OX1])[#6]", 0),     # amide
    ("[#6][OX2][#6]", 1),              # ether
    ("[SX2]", 0),                      # thiol or thioether
    ("[CX2]#[CX2]", 0),
    ("[NX1]#[CX2]", 0),                # nitrile
)


def _q18_all_hindered(mol: Chem.Mol) -> bool:
    """Branch (i): every functional group in the substance is hindered.

    Cramer does not define "sterically hindered", and Phillips et al. struck
    the phrase out for that reason, so the definition below is this file's
    choice and not the paper's. It is read as narrowly as the one worked
    example allows: a group is hindered when its attachment position is
    flanked on *two* sides by a quaternary carbon, which is the
    2,6-di-tert-butylphenol pattern of BHT. The flanking carbon may be the
    neighbouring atom itself or a substituent on it, so an ortho tert-butyl
    on a benzene ring blocks the phenol between two of them.

    The word "all" in the question is enforced: one free functional group
    anywhere in the molecule disqualifies it, so a hindered phenol carrying
    an unhindered carboxylic acid elsewhere stays in Class I.

    Measured, before it was written here, against both reference sets: it
    moves one compound, BHT, from Class I to its reference Class II, and
    breaks nothing in either set. Allowing tert-alkyl (isopropyl) flanks as
    well as fully quaternary ones was tried and changed no verdict in either
    set, so the narrower reading is kept. The looser readings that would also
    reach di-tert-butyl-4-hydroxymethylphenol all require treating a bare
    benzylic -CH2OH as hindered, which contradicts "all"; that compound is
    therefore a known miss, recorded in the README.
    """
    sites: list[tuple[int, set[int]]] = []
    for smarts, anchor in _Q18_HINDERED_GROUPS:
        patt = Chem.MolFromSmarts(smarts)
        for match in mol.GetSubstructMatches(patt):
            sites.append((match[anchor], set(match)))
    if not sites:
        return False

    for anchor, group in sites:
        atom = mol.GetAtomWithIdx(anchor)
        if atom.GetSymbol() == "C":
            attachment, blocked = anchor, group
        else:
            # A heteroatom group (-OH, -NH2, -SH) hangs off its own carbon.
            carbons = [n.GetIdx() for n in atom.GetNeighbors()
                       if n.GetSymbol() == "C"]
            if len(carbons) != 1:
                return False
            attachment, blocked = carbons[0], group | {anchor}

        seen = blocked | {attachment}
        flanks = 0
        for nb in mol.GetAtomWithIdx(attachment).GetNeighbors():
            if nb.GetIdx() in seen or nb.GetSymbol() == "H":
                continue
            if _is_quaternary_carbon(nb) or any(
                    _is_quaternary_carbon(n2) for n2 in nb.GetNeighbors()
                    if n2.GetIdx() not in seen):
                flanks += 1
        if flanks < 2:
            return False
    return True


def q18_feature(mol: Chem.Mol) -> str | None:
    """Which branch of Q18 the substance answers yes on, if any."""
    for letter, description, smarts in Q18_BRANCHES:
        if has(mol, smarts):
            return f"({letter}) {description}"
    if _q18_large_ketone(mol):
        return ("(h) acyclic aliphatic ketone with four or more carbons on "
                "one side of the keto group and no other functional groups")
    if _q18_all_hindered(mol):
        return ("(i) every functional group is sterically hindered, each "
                "flanked by two quaternary carbons")
    return None


def _oxygen_function_count(mol: Chem.Mol) -> int:
    """Number of readily metabolised oxygen functions."""
    return sum(count(mol, p) for p in
               (P_ALCOHOL, P_ACID, P_KETONE, P_ALDEHYDE))


def _sulphonate_salt(mol: Chem.Mol) -> bool:
    """True for a salt of a sulphonic, sulphamic or sulphuric acid.

    Cramer's final question: such salts are ionised at gut pH, are not
    absorbed to any appreciable extent, and are excreted unchanged. The
    acidic oxygen must be free rather than esterified to carbon -- a fully
    esterified sulphonate is an alkylating agent, not an unabsorbed salt.
    """
    groups = (count(mol, P_SULFONATE_SALT) + count(mol, P_SULFAMATE_SALT)
              + count(mol, P_SULFATE_SALT))
    if not groups:
        return False
    # A carboxylate alongside the sulphonate does not change the picture --
    # the sulphonated azo dyes carry both and are Class I -- but a phosphorus
    # acid does.
    if count(mol, P_PHOSPHORUS_ACID):
        return False
    # The three conditions below stand in for questions 16-32 of the original
    # tree, which are not implemented here. Without them a benzidine azo dye,
    # a sparsely sulphonated polycyclic dye and a chlorinated aryl sulphate
    # all escape to Class I, because nothing earlier in this shortened tree
    # catches them.
    if has(mol, "[#6][F,Cl,Br,I]"):
        return False
    if has(mol, P_PRIMARY_AROMATIC_AMINE):
        return False
    # Sulphonation has to be dense enough to prevent absorption, but that is
    # a question about large aromatic dyes: one sulphonate on a polycyclic dye
    # does not stop uptake, whereas a single ionised head group on an alkyl
    # chain does, which is why the laureth and lauryl sulphates are Class I at
    # any chain length.
    if not mol.GetRingInfo().NumRings() or not any(
            a.GetIsAromatic() for a in mol.GetAtoms()):
        return True
    return groups >= 2 or mol.GetNumHeavyAtoms() <= 25


def classify(mol: Chem.Mol, props: dict | None = None) -> CramerResult:
    """Walk the Cramer tree and return the class plus the decision path."""
    props = props or {}
    trace: list[str] = []
    elements = {a.GetSymbol() for a in mol.GetAtoms()}
    fired = [a.name for a, _ in screen(mol, CRAMER_HIGH_TOX)]

    def done_q18(q: str, answer: str) -> CramerResult:
        """A Class I verdict, unless Q18 intercepts it and makes it Class II.

        Q18 is the gate the original tree puts last on the open-chain and
        mononuclear branches, and it is where this implementation was losing
        Class II compounds: both branches ended at Class I without asking it.
        """
        feature = q18_feature(mol)
        if feature is None:
            return done(CLASS_I, q, answer)
        trace.append(f"{q} → {answer}")
        return done(CLASS_II,
                    "Q18. Is it one of the structures Cramer lists as "
                    "carrying enhanced toxicity?", f"Yes — {feature}")

    def done(cls: str, q: str, answer: str, ood: bool = False,
             food_question: bool = True) -> CramerResult:
        # Every Class III conclusion passes through the two questions the
        # original tree asks last. Their recovered positions matter: asking
        # the sulphonate question first (as an earlier version of this file
        # did) lets sulphonated dyes escape the structural questions that
        # should catch them, which is why it needs no ad-hoc guards here.
        if cls in (CLASS_II, CLASS_III) and not ood:
            food = reflists.food_component(mol)
            unusual_elements = bool(elements - CRAMER_ALLOWED)
            if (cls == CLASS_III and food and food_question
                    and not unusual_elements
                    and not has(mol, P_ACYL_HETEROATOM)):
                trace.append(f"{q} → {answer}")
                cls, q, answer = (
                    CLASS_II,
                    "Q22. Is it a common component of food?",
                    f"Yes — {food.title()} (reference list)")
            elif _sulphonate_salt(mol):
                trace.append(f"{q} → {answer}")
                cls, q, answer = (
                    CLASS_I,
                    "Q33. Is it a salt of a sulphonic, sulphamic or "
                    "sulphuric acid, readily excreted unchanged?",
                    "Yes")
        trace.append(f"{q} → {answer} ⇒ {cls}")
        # Food-list membership is recorded as context, never as a class
        # modifier. Tempting though it is -- every food-listed compound in the
        # Cramer Appendix 1 set is Class II -- it does not generalise: the same
        # list contains ethanol, vanillin and limonene (Class I) and caffeine
        # (Class III). Membership is necessary context for the questions
        # Cramer asks about food components, not sufficient to answer them.
        food = reflists.food_component(mol)
        if food:
            trace.append(f"Note: listed as a common component of food "
                         f"({food.title()}). Recorded as context; the Cramer "
                         f"class is unchanged by it.")
        return CramerResult(cls, trace, fired, ood, food_component=food)

    def no(q: str, answer: str = "No") -> None:
        trace.append(f"{q} → {answer}")

    if props.get("counterions"):
        trace.append("Q0. Benign counterions removed before classification: "
                     + ", ".join(props["counterions"])
                     + " — the organic parent is classified.")

    # Q2 — normal body constituent. Uses the curated reference list when it has
    # been fetched, and the small built-in metabolite list otherwise.
    # The shortcut is guarded: being a body constituent does not override the
    # tree's own positive-concern questions. Biotin and ascorbic acid are both
    # body constituents, but the first carries sulphur in a heterocycle and the
    # second is a lactone, and the reference implementation classifies them II
    # and III rather than I. Without this guard the curated list raises overall
    # agreement but halves Class III accuracy and doubles the under-calls.
    # Q1 — a normal constituent of the body. This is genuinely the tree's
    # first question, and it has to precede the applicability check below:
    # water answers it yes, and would otherwise be rejected as "inorganic".
    q = "Q1. Is it a normal constituent of the body?"
    if reflists.primary_metabolite(mol):
        return done(CLASS_I, q, "Yes — a building block of intermediary "
                                "metabolism")
    no(q)

    # Q1b — applicability domain.
    q = "Q1b. Is the substance organic (carbon-containing) and metal-free?"
    if "C" not in elements:
        return done(CLASS_III, q, "No — inorganic", ood=True)
    if props.get("has_heavy_metal"):
        return done(CLASS_III, q, "No — contains a metal/metalloid", ood=True)
    no(q, "Yes")

    # Q3 — elements outside the Cramer set
    q = "Q3. Does it contain elements other than C, H, N, O, S and halogen?"
    extra = sorted(elements - CRAMER_ALLOWED)
    if extra:
        return done(CLASS_III, q, f"Yes — {', '.join(extra)}")
    no(q)

    # Q4 — functional groups associated with enhanced toxicity
    q = "Q4. Does it contain a functional group associated with enhanced toxicity?"
    if fired:
        return done(CLASS_III, q, f"Yes — {', '.join(fired)}")
    no(q)

    # Q4b — carbohydrates
    q = "Q4b. Is it a carbohydrate (a cyclic sugar or a glycoside of one)?"
    if _is_carbohydrate(mol):
        return done(CLASS_I, q, "Yes")
    no(q)

    # Q5 — common terpene / terpenoid of food origin
    q = "Q5. Is it a common terpene or terpenoid hydrocarbon of food origin?"
    is_terpene = (
        elements <= {"C", "H", "O"}
        and has(mol, P_TERPENE_SKELETON)
        and 10 <= mol.GetNumHeavyAtoms() <= 20
        and not has(mol, P_BENZENE)
        and mol.GetRingInfo().NumRings() <= 2
    )
    if is_terpene:
        return done(CLASS_I, q, "Yes")
    no(q)

    # Q6 — lactone or cyclic diester
    q = "Q6. Is it a lactone or a cyclic diester?"
    if has(mol, P_CYCLIC_DIESTER):
        return done(CLASS_III, q, "Yes — cyclic diester")
    if has(mol, P_LACTONE):
        small = (mol.GetNumHeavyAtoms() <= Q6_SIMPLE_LACTONE_MAX_HEAVY
                 and not has(mol, P_BENZENE))
        return done(CLASS_II if small else CLASS_III, q,
                    "Yes — simple lactone" if small else "Yes — complex lactone")
    no(q)

    # Q7 — sulphur other than sulphonate/sulphamate
    q = "Q7. Does it contain sulphur other than as sulphonate or sulphamate?"
    if "S" in elements and not (has(mol, P_SULFONATE) or has(mol, P_SULFAMATE)):
        # The recovered topology routes sulphur compounds away from the food
        # question, which is why thiophene is Class III despite being listed.
        return done(CLASS_III, q, "Yes", food_question=False)
    no(q)

    # Q8 — heterocyclic systems
    q = "Q8. Is it a heterocyclic compound?"
    heterocyclic = has(mol, P_HETEROAROMATIC) or has(mol, P_HETEROCYCLE_SAT)
    if heterocyclic:
        if (has(mol, P_COMMON_FOOD_HETEROCYCLE)
                and mol.GetNumHeavyAtoms() <= Q8_SIMPLE_HETEROCYCLE_MAX_HEAVY):
            return done(CLASS_II, q, "Yes — simple heterocycle common in food")
        return done(CLASS_III, q, "Yes — other heterocycle")
    no(q)

    # Q9 — benzene derivatives
    q = "Q9. Is it a benzene derivative bearing only innocuous substituents?"
    if has(mol, P_BENZENE):
        # Cramer places benzene derivatives in Class I on the nature of the
        # substituents, not on molecular size: propyl paraben and propylene
        # glycol dibenzoate are both Class I. Only a substituent of concern,
        # or fusion to a second ring, moves the compound up.
        if _only_harmless_aryl_substituents(mol):
            return done_q18(q, "Yes — only readily metabolised substituents")
        return done(CLASS_III, q, "No — bears a substituent of concern")
    no(q)

    # Q10 — polycyclic non-aromatic (fused/bridged) carbocycle
    q = "Q10. Is it a fused or bridged polycyclic carbocycle?"
    if mol.GetRingInfo().NumRings() >= 2 and has(mol, P_POLYCYCLIC_NONAROMATIC):
        return done(CLASS_II, q, "Yes")
    no(q)

    # Q11 — alpha,beta-unsaturation and other activated aliphatic systems.
    # Cramer treats these separately from the readily metabolised oxygen
    # functions of Q12 because conjugation to the carbonyl makes them
    # electrophilic rather than simply oxidisable.
    q = ("Q11. Is it an activated unsaturated system (acrylate, unsaturated "
         "nitrile), an acetylenic compound, or an aliphatic halide?")
    activated = [
        ("acrylate or methacrylate ester", "[CX3]=[CX3][CX3](=[OX1])[OX2][#6]"),
        ("alpha,beta-unsaturated nitrile or sulphone",
         "[CX3]=[CX3][$([CX2]#[NX1]),$([SX4](=O)=O)]"),
        ("allylic or benzylic leaving group",
         "[$([CX3]=[CX3][CX4][F,Cl,Br,I]),$(c[CX4][F,Cl,Br,I])]"),
        ("vinyl ether or enol ether", "[CX3]=[CX3][OX2][#6]"),
        ("acetylenic system", "[CX2]#[CX2]"),
        ("alpha-halo carbonyl", "[F,Cl,Br,I][CX4][CX3]=[OX1]"),
        ("aliphatic halide", P_ALIPHATIC_HALIDE),
    ]
    for label, patt in activated:
        if has(mol, patt):
            return done(CLASS_III, q, f"Yes — {label}")
    no(q)

    # Q12 — simple aliphatic with readily metabolised oxygen function
    q = ("Q12. Is it an acyclic or monocyclic aliphatic compound whose only "
         "functions are alcohol, aldehyde, ketone, acid, ester or ether?")
    oxygen_funcs = any(has(mol, p) for p in
                       (P_ALCOHOL, P_ALDEHYDE, P_KETONE, P_ACID, P_ESTER, P_ETHER))
    only_cho = elements <= {"C", "H", "O"}
    if only_cho and (oxygen_funcs or elements <= {"C", "H"}):
        # No size gate: the reference data puts straight-chain fatty acids up
        # to C24, and lycopene, in Class I. Chain length is not what Cramer
        # asks about here -- it is how densely oxygenated the chain is. One
        # or two oxygen functions still read as readily metabolised (ethanol,
        # propylene glycol, a fatty acid, acetoin); three or more, or an
        # alpha-dicarbonyl, do not, which is where the sugar alcohols, citric
        # acid, oxalic acid and diacetyl part company with them.
        oxygens = _oxygen_function_count(mol)
        alpha_dicarbonyl = has(mol, P_ALPHA_DICARBONYL)
        if oxygens >= Q12_MAX_OXYGEN_FUNCTIONS + 1:
            no(q, f"Yes, but it carries {oxygens} oxygen functions")
        elif alpha_dicarbonyl:
            no(q, "Yes, but it carries an alpha-dicarbonyl")
        elif mol.GetRingInfo().NumRings() <= Q12_MAX_RINGS:
            return done_q18(q, "Yes")
        else:
            return done(CLASS_II, q, "Yes, but polycyclic")
    else:
        no(q)

    # Q13 — nitrogen-containing aliphatics and everything else
    q = "Q13. Does it contain nitrogen in an aliphatic amine or amide function?"
    if "N" in elements:
        if has(mol, P_AMIDE) and mol.GetNumHeavyAtoms() <= 12:
            return done(CLASS_II, q, "Yes — simple aliphatic amide")
        return done(CLASS_III, q, "Yes — other nitrogen function")
    no(q)

    q = "Q14. Default — no positive presumption of safety was established."
    return done(CLASS_II, q, "Default")


def _only_harmless_aryl_substituents(mol: Chem.Mol) -> bool:
    """True when every non-ring neighbour of the aromatic ring is innocuous."""
    benz = Chem.MolFromSmarts(P_BENZENE)
    harmless = Chem.MolFromSmarts(P_HARMLESS_AR_SUBST)
    harmless_atoms = {m[0] for m in mol.GetSubstructMatches(harmless)}
    for match in mol.GetSubstructMatches(benz):
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in match or nb.GetSymbol() == "H":
                    continue
                if nb.IsInRing():
                    return False
                if idx not in harmless_atoms:
                    return False
    return True
