"""A pure-Python toxicity decision-tree engine.

Modules
-------
chem      structure parsing, descriptors and depiction
alerts    the SMARTS structural-alert library
cramer    Cramer classification of oral toxicity classes
verhaar   Verhaar aquatic mode-of-action scheme
engine    orchestration plus the Kroes TTC decision tree
examples  a small set of reference structures
"""

__version__ = "1.0.0"

from .chem import MolRecord, parse_many, parse_structure  # noqa: F401
from .engine import Assessment, SCHEME_LABELS, assess, assess_many  # noqa: F401
