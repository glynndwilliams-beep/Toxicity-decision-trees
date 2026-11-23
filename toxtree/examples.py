"""Reference structures used to populate the UI and the self-test."""
from __future__ import annotations

EXAMPLES: list[tuple[str, str, str]] = [
    # (name, SMILES, what it is meant to demonstrate)
    ("Ethanol", "CCO", "Simple aliphatic alcohol — Cramer I, Verhaar inert"),
    ("Vanillin", "COc1cc(C=O)ccc1O", "Flavouring with an aldehyde and phenol"),
    ("Benzene", "c1ccccc1", "Aromatic hydrocarbon — baseline narcosis"),
    ("Phenol", "Oc1ccccc1", "Polar narcotic, irritant"),
    ("Aniline", "Nc1ccccc1", "Aromatic amine — genotox alert"),
    ("4-Nitroaniline", "Nc1ccc(cc1)[N+](=O)[O-]", "Nitroaromatic plus aromatic amine"),
    ("Benzo[a]pyrene", "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",
     "Polycyclic aromatic hydrocarbon"),
    ("Styrene oxide", "C1OC1c1ccccc1", "Epoxide — direct alkylator"),
    ("Acrolein", "C=CC=O", "Michael acceptor and aldehyde"),
    ("Formaldehyde", "C=O", "Simple aldehyde"),
    ("N-Nitrosodimethylamine", "CN(C)N=O", "N-nitroso genotoxic carcinogen"),
    ("Parathion", "CCOP(=S)(OCC)Oc1ccc(cc1)[N+](=O)[O-]",
     "Organophosphate — Verhaar class 4, AChE TTC tier"),
    ("Carbaryl", "CNC(=O)Oc1cccc2ccccc12", "Carbamate insecticide"),
    ("Atrazine", "CCNc1nc(Cl)nc(NC(C)C)n1", "Triazine herbicide"),
    ("2,4-Dinitrophenol", "Oc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]", "Uncoupler"),
    ("Hexachlorobenzene", "Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl", "Polyhalogenated aromatic"),
    ("Glucose", "OCC1OC(O)C(O)C(O)C1O", "Normal body constituent — Cramer I"),
    ("Glycine", "NCC(=O)O", "Amino acid — normal body constituent"),
    ("Citric acid", "OC(=O)CC(O)(CC(=O)O)C(=O)O", "Intermediary metabolite"),
    ("Limonene", "CC(=C)C1CCC(C)=CC1", "Terpene of food origin"),
    ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Purine alkaloid heterocycle"),
    ("Cinnamaldehyde", "O=C/C=C/c1ccccc1", "Sensitiser — Michael acceptor plus aldehyde"),
    ("2,4-Dinitrochlorobenzene", "Clc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]",
     "Strong sensitiser — SNAr electrophile"),
    ("Hydroquinone", "Oc1ccc(O)cc1", "Pro-electrophilic sensitiser"),
    ("Glyceryl trioleate", "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(COC(=O)CCCCCCC/C=C\\CCCCCCCC)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
     "Large lipophilic ester — outside the aquatic QSAR domain"),
    ("Methyl acrylate", "COC(=O)C=C", "Acrylate ester — reactive"),
    ("Dimethyl phthalate", "COC(=O)c1ccccc1C(=O)OC", "Phthalate ester"),
    ("Sodium lauryl sulfate", "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]", "Anionic surfactant"),
    ("Thiourea", "NC(N)=S", "Thyroid-active non-genotoxic alert"),
    ("Ethylene glycol", "OCCO", "Aliphatic polyol"),
]


def example_text() -> str:
    """The examples as a SMILES + name list suitable for the batch text area."""
    return "\n".join(f"{smi}\t{name}" for name, smi, _ in EXAMPLES)
