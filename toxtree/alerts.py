"""SMARTS alert library.

Each :class:`Alert` is one structural alert: a name, a SMARTS pattern, the
scheme it belongs to and a short mechanistic rationale. The alert sets follow
the published rulebases (Benigni & Bossa for carcinogenicity/mutagenicity,
the mechanistic protein-binding domains for skin sensitisation, Cramer's
"functional groups of enhanced toxicity", and the Verhaar mode-of-action
classes). They are re-derived from the literature, not copied from any
existing implementation.
"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from .chem import count, has


@dataclass(frozen=True)
class Alert:
    code: str
    name: str
    smarts: str
    scheme: str
    rationale: str = ""

    def matches(self, mol: Chem.Mol) -> list[tuple[int, ...]]:
        patt = Chem.MolFromSmarts(self.smarts)
        if patt is None:
            return []
        return list(mol.GetSubstructMatches(patt, uniquify=True))


# --------------------------------------------------------------------------
# Genotoxic carcinogenicity / mutagenicity (Benigni–Bossa style)
# --------------------------------------------------------------------------
GENOTOX: list[Alert] = [
    Alert("SA1", "Acyl halide", "[CX3](=[OX1])[F,Cl,Br,I]", "genotox",
          "Direct acylating agent; reacts with nucleophilic DNA/protein sites."),
    Alert("SA2", "Alkyl ester of sulphonic or phosphonic acid",
          "[$([#6][SX4](=O)(=O)[OX2][CX4]),$([#6][PX4](=O)([OX2,OX1-])[OX2][CX4])]", "genotox",
          "Alkylating agent via nucleophilic substitution at carbon."),
    Alert("SA2b", "Dialkyl sulfate or trialkyl phosphate ester",
          "[$([#6][OX2][SX4](=[OX1])(=[OX1])[OX2][#6]),"
          "$([#6][OX2][PX4](=[OX1,SX1])([OX2][#6])[OX2,SX2][#6])]", "genotox",
          "Fully esterified: a direct alkylating agent at carbon. The "
          "mono-ester salts are not, and are excluded by requiring carbon on "
          "every acidic oxygen."),
    Alert("SA3", "N-methylol derivative", "[NX3;!$(N=*)][CX4H2][OX2H]", "genotox",
          "Releases formaldehyde, a DNA-protein cross-linker."),
    Alert("SA4", "Monohaloalkene", "[CX3]=[CX3][F,Cl,Br,I]", "genotox",
          "Metabolically activated to a reactive epoxide."),
    Alert("SA5", "S or N mustard", "[Cl,Br,I][CX4][CX4][S,N;!$(N=*)]", "genotox",
          "Forms a cyclic onium ion that bifunctionally alkylates DNA."),
    Alert("SA6", "Beta-lactone / sultone", "[$(O=C1OCC1),$(O=S1(=O)OCC1)]", "genotox",
          "Strained ring opens to acylate or sulphonate nucleophiles."),
    Alert("SA7", "Epoxide or aziridine", "[OX2,NX3;r3]1[#6;r3][#6;r3]1", "genotox",
          "Strained three-membered ring is a direct-acting alkylator."),
    Alert("SA8", "Aliphatic halide", "[CX4;!$(C(F)(F)F)][F,Cl,Br,I]", "genotox",
          "Potential SN2 alkylating agent."),
    Alert("SA9", "Alkyl nitrite", "[CX4]O[NX2]=[OX1]", "genotox",
          "Nitrosating agent; can generate alkylating nitrosamines."),
    Alert("SA10", "Alpha,beta-unsaturated carbonyl, nitrile or sulphone",
          "[CX3]=[CX3][$([CX3]=[OX1]),$([CX2]#[NX1]),$([SX4](=O)=O)]", "genotox",
          "Michael acceptor; alkylates DNA bases at nucleophilic nitrogen."),
    Alert("SA11", "Simple aldehyde", "[$([CX3H1](=[OX1])[#6;!$([CX3]=[OX1])]),$([CX3H2]=[OX1])]", "genotox",
          "Forms Schiff bases and DNA-protein cross-links."),
    Alert("SA12", "Quinone", "[$(O=C1C=CC(=O)C=C1),$(O=C1C=CC=CC1=O)]", "genotox",
          "Redox cycling plus Michael-type addition to DNA."),
    Alert("SA13", "Hydrazine", "[NX3;!$(N=*)][NX3;!$(N=*)]", "genotox",
          "Metabolised to carbonium/radical species that methylate DNA."),
    Alert("SA14", "Aliphatic azo or azoxy", "[CX4][NX2]=[NX2,NX3][CX4]", "genotox",
          "Activated to alkyl diazonium ions."),
    Alert("SA15", "Isocyanate or isothiocyanate", "[NX2]=[CX2]=[OX1,SX1]", "genotox",
          "Direct carbamoylating agent."),
    Alert("SA16", "Alkyl carbamate or thiocarbamate",
          "[NX3][CX3](=[OX1,SX1])[OX2][CX4]", "genotox",
          "Carbamoylates DNA and protein nucleophiles."),
    Alert("SA18", "Polycyclic aromatic hydrocarbon", "c1ccc2c(c1)ccc1ccccc12", "genotox",
          "Bay/fjord-region diol epoxides form bulky DNA adducts."),
    Alert("SA19", "Heterocyclic polycyclic aromatic", "[$(c1ccc2c(c1)ccc1ccncc12),$(c1ccc2c(c1)[nH]c1ccccc12)]",
          "genotox", "Heteroaromatic analogue of PAH activation."),
    Alert("SA21", "N-nitroso", "[NX3;!$(N=*)][NX2]=[OX1]", "genotox",
          "Alpha-hydroxylation yields a potent alkyldiazonium alkylator."),
    Alert("SA22", "Azide or triazene", "[$([NX2]=[NX2+]=[NX1-]),$(N=NN)]", "genotox",
          "Decomposes to diazonium/carbonium electrophiles."),
    Alert("SA23", "Aliphatic N-nitro", "[CX4][NX3][$([NX3](=O)=O),$([NX3+](=O)[O-])]", "genotox",
          "Nitramine; alkylating on metabolic activation."),
    Alert("SA24", "Alpha,beta-unsaturated alkoxy / vinyl ether",
          "[CX3]=[CX3][OX2][#6]", "genotox",
          "Oxidised to a reactive epoxide."),
    Alert("SA25", "Aromatic nitroso", "[c;!$(c~[n,o,s])][NX2]=[OX1]", "genotox",
          "Reactive nitrogen electrophile; forms nitrenium ions."),
    Alert("SA26", "Aromatic ring N-oxide", "[n+][O-]", "genotox",
          "Reduced to a nitrenium-type electrophile."),
    Alert("SA27", "Nitroaromatic", "[c][$([NX3](=O)=O),$([NX3+](=O)[O-])]", "genotox",
          "Nitroreduction gives a nitrenium ion that binds DNA."),
    Alert("SA28", "Primary aromatic amine / hydroxylamine",
          "[c;!$(c~[n,o,s])][NX3;H2,H1;!$(N[CX3]=[OX1])]", "genotox",
          "N-hydroxylation and conjugation give a DNA-binding nitrenium ion."),
    Alert("SA29", "Aromatic mono- or dialkylamine", "[c;!$(c~[n,o,s])][NX3]([CX4])[#6,#1]", "genotox",
          "N-dealkylation/N-oxidation route to the same nitrenium species."),
    Alert("SA30", "Aromatic N-acylamine", "[c;!$(c~[n,o,s])][NX3;H1][CX3](=[OX1])[#6]", "genotox",
          "Deacetylated to the primary aromatic amine."),
    Alert("SA31", "Aromatic diazo", "[c;!$(c~[n,o,s])][NX2]=[NX2][#6]", "genotox",
          "Azo reduction releases aromatic amines."),
    Alert("SA32", "Coumarin / furocoumarin", "[$(O=c1ccc2ccccc2o1),$(c1cc2ccoc2cc1)]", "genotox",
          "Photoactivated or epoxidised to a DNA-reactive species."),
    Alert("SA33", "Halogenated dibenzodioxin / biphenyl core",
          "[$(c1ccc2Oc3ccccc3Oc2c1),$(c1ccc(-c2ccccc2)cc1[F,Cl,Br])]", "genotox",
          "Persistent aryl-hydrocarbon receptor ligand; tumour promoter."),
    Alert("SA34", "Alkyl hydroperoxide or peroxide", "[OX2][OX2]", "genotox",
          "Generates oxygen radicals causing oxidative DNA damage."),
    Alert("SA35", "Thiourea / thioamide", "[NX3][CX3](=[SX1])[NX3,#6]", "genotox",
          "Oxidised to reactive sulphene/sulphine metabolites."),
    Alert("SA36", "Aromatic sulphonate ester", "[c][SX4](=[OX1])(=[OX1])[OX2][CX4]", "genotox",
          "Direct alkylating agent."),
    Alert("SA37", "Nitrogen or sulphur ylide / diazoketone", "[CX3](=[OX1])[CX2]=[NX2+]=[NX1-]",
          "genotox", "Carbene/diazo precursor; direct DNA alkylation."),
    Alert("SA38", "Propargyl / allyl halide", "[$([CX2]#[CX2][CX4][F,Cl,Br,I]),$([CX3]=[CX3][CX4][F,Cl,Br,I])]",
          "genotox", "Highly activated SN2/SN2' alkylating agent."),
]

# --------------------------------------------------------------------------
# Non-genotoxic carcinogenicity
# --------------------------------------------------------------------------
NONGENOTOX: list[Alert] = [
    Alert("NG1", "Halogenated benzene / phenol", "[F,Cl,Br,I]c1ccccc1", "nongenotox",
          "Peroxisome proliferation and sustained hepatocyte stress."),
    Alert("NG2", "Chlorinated alkane (C1-C4)", "[CX4]([F,Cl,Br,I])([F,Cl,Br,I])", "nongenotox",
          "Cytotoxicity with regenerative hyperplasia."),
    Alert("NG4", "Phenoxy or chlorophenoxy acid", "[c][OX2][CX4][CX3](=[OX1])[OX2H,OX1-]", "nongenotox",
          "Peroxisome proliferator-activated receptor agonism."),
    Alert("NG5", "Perfluoroalkyl acid/sulphonate",
          "[$(C(F)(F)C(F)(F)C(=O)[OX2H,OX1-]),$(C(F)(F)S(=O)(=O)[OX2H,OX1-])]", "nongenotox",
          "Persistent peroxisome proliferator."),
    Alert("NG6", "Thiocarbamide / imidazole antithyroid motif",
          "[$([NX3][CX3](=[SX1])[NX3]),$(c1cnc[nH]1)]", "nongenotox",
          "Thyroid peroxidase inhibition leading to follicular tumours."),
    Alert("NG7", "Long-chain phthalate ester",
          "[CX4][OX2][CX3](=O)c1ccccc1[CX3](=O)[OX2][CX4]", "nongenotox",
          "Hepatic peroxisome proliferation."),
]

# --------------------------------------------------------------------------
# Skin sensitisation — mechanistic protein-binding domains
# --------------------------------------------------------------------------
SENSITISATION: list[Alert] = [
    Alert("SS1", "Michael acceptor", "[CX3]=[CX3][$([CX3]=[OX1]),$([SX4](=O)=O),$([CX2]#[NX1])]",
          "sensitisation", "Conjugate addition of skin-protein thiol/amine nucleophiles."),
    Alert("SS2", "Schiff base former (aldehyde)",
          "[$([CX3H1](=[OX1])[#6]),$([CX3H2]=[OX1]),$([CX3](=[OX1])([#6])[CX4][F,Cl,Br,I])]",
          "sensitisation", "Condenses with lysine amino groups to form an imine."),
    Alert("SS2b", "1,2- or 1,3-dicarbonyl",
          "[$([#6][CX3](=[OX1])[CX3](=[OX1])[#6]),$([CX3](=[OX1])[CX4;H1,H2][CX3]=[OX1])]",
          "sensitisation",
          "Alpha- and beta-dicarbonyls condense with lysine far more readily "
          "than an isolated ketone; the diketone sensitisers of the reference "
          "set are all of this type."),
    Alert("SS3", "Acyl transfer agent",
          "[$([CX3](=[OX1])[F,Cl,Br,I]),$([CX3](=[OX1])O[CX3]=[OX1]),$([NX2]=[CX2]=[OX1])]",
          "sensitisation", "Acylates lysine side chains."),
    Alert("SS3b", "Lactone or azlactone",
          "[$([#6;R][CX3;R](=[OX1])[OX2;R]),$([#6;R]1=[#7;R][#6;R][CX3;R](=[OX1])[OX2;R]1)]",
          "sensitisation",
          "Ring-opening acyl transfer to protein nucleophiles. The azlactones "
          "are the single largest acyl-transfer group in the reference data."),
    Alert("SS4", "SN2 electrophile",
          "[$([CX4;H1,H2][Cl,Br,I]),$(c[CX4;H1,H2][Cl,Br,I]),$([CX4][OX2][SX4](=O)=O)]",
          "sensitisation",
          "Direct alkylation of cysteine thiol. Primary and benzylic chlorides "
          "count: the long-chain alkyl chlorides are sensitisers here."),
    Alert("SS5", "SNAr electrophile",
          "[$(c([F,Cl,Br,I])1c([$([NX3](=O)=O),$([NX3+](=O)[O-])])cccc1),"
          "$(c([F,Cl,Br,I])1cccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]),"
          "$([F,Cl,Br,I]c1nc([F,Cl,Br,I])nc([F,Cl,Br,I])n1)]",
          "sensitisation", "Nucleophilic aromatic substitution on an activated ring."),
    Alert("SS6", "Pro-electrophilic phenol (quinone or quinone-methide precursor)",
          "[$(c1cc([OX2H])ccc1[OX2H]),$(c1cccc([OX2H])c1[OX2H]),"
          "$([OX2H]c1ccc([CX4;!$(C=O),CX3;!$([CX3]=[OX1])])cc1),"
          "$([OX2H]c1ccccc1[CX4;!$(C=O),CX3;!$([CX3]=[OX1])]),"
          "$([OX2H]c1ccc([NX3])cc1),$([OX2H]c1ccccc1[NX3])]",
          "sensitisation",
          "Oxidised in skin to a quinone or quinone methide. This covers the "
          "eugenol/isoeugenol family, which the catechol pattern alone misses."),
    Alert("SS7", "Epoxide / aziridine", "[OX2,NX3;r3]1[#6;r3][#6;r3]1", "sensitisation",
          "Ring-opening alkylation of protein nucleophiles."),
    Alert("SS8", "Quinone / quinone imine", "[$(O=C1C=CC(=O)C=C1),$(O=C1C=CC(=N)C=C1)]",
          "sensitisation", "Strong Michael acceptor formed in situ."),
    Alert("SS9", "Primary aromatic amine", "[c;!$(c~[n,o,s])][NX3;H2]", "sensitisation",
          "Oxidised to a quinone imine or nitroso sensitiser."),
    Alert("SS10", "Isothiocyanate / thiol-reactive sulphur", "[NX2]=[CX2]=[SX1]",
          "sensitisation", "Thiocarbamoylates protein nucleophiles."),
    Alert("SS11", "Isothiazolinone",
          "[#16X2;r5][#7;r5][#6;r5]=[OX1]", "sensitisation",
          "The ring opens at the S-N bond on attack by a protein thiol; the "
          "isothiazolinone preservatives are potent sensitisers."),
    Alert("SS12", "Disulfide or thiuram",
          "[#16X2][#16X2]", "sensitisation",
          "Thiol-disulfide exchange with cysteine residues."),
]

# --------------------------------------------------------------------------
# Local irritation / corrosion
# --------------------------------------------------------------------------
IRRITATION: list[Alert] = [
    Alert("IR1", "Carboxylic acid",
          "[$([CX3](=[OX1])[OX2H]);!$([CX3](=[OX1])([OX2H])[CX4][NX3;H1,H2])]", "irritation",
          "Low pH denatures skin and ocular protein. An alpha-amino acid is "
          "excluded: it is zwitterionic at skin pH, so neither the carboxyl "
          "nor the amine behaves as a free acid or base."),
    Alert("IR2", "Sulphonic / sulphuric acid", "[SX4](=[OX1])(=[OX1])[OX2H]", "irritation",
          "Strongly acidic; corrosive."),
    Alert("IR3", "Acyl halide or anhydride",
          "[$([CX3](=[OX1])[F,Cl,Br,I]),$([CX3](=[OX1])O[CX3]=[OX1])]", "irritation",
          "Hydrolyses to acid in tissue; acylates protein."),
    Alert("IR4", "Aliphatic amine (base)", "[NX3;H2,H1,H0;!$(N[CX3]=[OX1]);!$(N=*);!$(Nc);!$(N[#7]);!$(N[SX4]);!$(N[PX4]);!$(N[CX4][CX3](=[OX1])[OX2H,OX1-])][CX4]",
          "irritation", "Alkaline; saponifies membrane lipids."),
    Alert("IR5", "Quaternary ammonium surfactant", "[NX4+]([CX4])([CX4])([CX4])[CX4]", "irritation",
          "Cationic surfactant disrupts cell membranes."),
    Alert("IR6", "Epoxide", "[OX2;r3]1[#6;r3][#6;r3]1", "irritation",
          "Direct alkylation of skin protein."),
    Alert("IR7", "Phenol", "[c][OX2H]", "irritation",
          "Protein denaturant and membrane disruptor."),
    Alert("IR8", "Alpha,beta-unsaturated carbonyl", "[CX3]=[CX3][CX3]=[OX1]", "irritation",
          "Reactive Michael acceptor; cytotoxic to epithelium."),
    Alert("IR9", "Peroxide", "[OX2][OX2]", "irritation",
          "Oxidative tissue damage."),
    Alert("IR10", "Isocyanate", "[NX2]=[CX2]=[OX1]", "irritation",
          "Acylating respiratory and dermal irritant."),
]

# --------------------------------------------------------------------------
# Cramer "functional groups associated with enhanced toxicity" (Q2 of the tree)
# --------------------------------------------------------------------------
CRAMER_HIGH_TOX: list[Alert] = [
    Alert("CQ2a", "Aliphatic secondary amine", "[NX3;H1]([CX4])[CX4]", "cramer", ""),
    Alert("CQ2b", "N-nitroso", "[NX3][NX2]=[OX1]", "cramer", ""),
    Alert("CQ2c", "Diazo or triazene", "[$([NX2]=[NX2+]=[NX1-]),$(N=NN)]", "cramer", ""),
    Alert("CQ2d", "Quaternary nitrogen (not trimethylammonium)",
          "[NX4+;!$([NX4+](C)(C)(C)[CX4])]", "cramer", ""),
    Alert("CQ2e", "N or S mustard", "[Cl,Br,I][CX4][CX4][S,N;!$(N=*)]", "cramer", ""),
    Alert("CQ2f", "Three-membered heterocycle", "[OX2,NX3,SX2;r3]1[#6;r3][#6;r3]1", "cramer", ""),
    Alert("CQ2g", "Aliphatic azo or azoxy", "[CX4][NX2]=[NX2,NX3][CX4]", "cramer", ""),
    Alert("CQ2h", "Hydrazine", "[NX3;!$(N=*)][NX3;!$(N=*)]", "cramer", ""),
    Alert("CQ2i", "Aliphatic N-nitro", "[CX4][NX3][$([NX3](=O)=O),$([NX3+](=O)[O-])]", "cramer", ""),
    Alert("CQ2j", "Nitrile / cyanide", "[CX2]#[NX1]", "cramer", ""),
    Alert("CQ2k", "Isocyanate or isothiocyanate", "[NX2]=[CX2]=[OX1,SX1]", "cramer", ""),
    Alert("CQ2l", "Carbamate or thiocarbamate ester",
          "[NX3][CX3](=[OX1,SX1])[OX2,SX2][#6]", "cramer", ""),
    Alert("CQ2m", "Organophosphate / phosphonate / phosphorothioate",
          "[PX4](=[OX1,SX1])([OX2,SX2,NX3])[OX2,SX2,NX3]", "cramer", ""),
    Alert("CQ2n", "Aromatic nitro", "[c][$([NX3](=O)=O),$([NX3+](=O)[O-])]", "cramer", ""),
]

ALL_ALERTS: list[Alert] = (
    GENOTOX + NONGENOTOX + SENSITISATION + IRRITATION + CRAMER_HIGH_TOX
)


def screen(mol: Chem.Mol, alerts: list[Alert]) -> list[tuple[Alert, list[tuple[int, ...]]]]:
    """Return every alert in *alerts* that fires, with its atom matches."""
    out = []
    for alert in alerts:
        m = alert.matches(mol)
        if m:
            out.append((alert, m))
    return out


def validate_library() -> list[str]:
    """Return the codes of any alerts whose SMARTS does not compile."""
    bad = []
    for alert in ALL_ALERTS:
        if Chem.MolFromSmarts(alert.smarts) is None:
            bad.append(f"{alert.code} ({alert.name}): {alert.smarts}")
    return bad


# --------------------------------------------------------------------------
# Computed (topology-based) alerts — expressed in code because ring-fusion
# patterns are clearer and safer than deeply nested SMARTS.
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class ComputedAlert:
    code: str
    name: str
    scheme: str
    rationale: str
    fn: "object"

    def fires(self, mol: Chem.Mol) -> bool:
        return bool(self.fn(mol))


def _fused_ring_sizes(mol: Chem.Mol) -> list[list[int]]:
    """Group rings into fused systems and return the ring sizes of each."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    systems: list[list[set[int]]] = []
    for ring in rings:
        joined = [s for s in systems if any(len(ring & r) >= 2 for r in s)]
        if not joined:
            systems.append([ring])
        else:
            merged = [ring]
            for s in joined:
                merged.extend(s)
                systems.remove(s)
            systems.append(merged)
    return [sorted(len(r) for r in s) for s in systems]


def is_steroid_like(mol: Chem.Mol) -> bool:
    """Four fused rings with the 6-6-6-5 cyclopenta[a]phenanthrene pattern."""
    return any(sizes == [5, 6, 6, 6] for sizes in _fused_ring_sizes(mol))


def fused_aromatic_count(mol: Chem.Mol) -> int:
    """Largest number of aromatic rings in a single fused ring system."""
    ri = mol.GetRingInfo()
    arom = [set(r) for r in ri.AtomRings()
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    best = 0
    for seed in arom:
        system = {id(seed): seed}
        changed = True
        while changed:
            changed = False
            for other in arom:
                if id(other) in system:
                    continue
                if any(len(other & r) >= 2 for r in system.values()):
                    system[id(other)] = other
                    changed = True
        best = max(best, len(system))
    return best


COMPUTED: list[ComputedAlert] = [
    ComputedAlert("NG3", "Steroid-like fused tetracycle", "nongenotox",
                  "Hormonal (receptor-mediated) promotion of tumour growth.",
                  is_steroid_like),
    ComputedAlert("SA17", "Three or more fused aromatic rings", "genotox",
                  "Polycyclic aromatic system; bay-region diol epoxide formation.",
                  lambda m: fused_aromatic_count(m) >= 3),
]
