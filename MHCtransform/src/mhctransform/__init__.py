"""Human/mouse MHC full, AF3-chain and 34-aa pseudo conversion."""

from .errors import MHCTransformError
from .models import (
    AF3MoleculeCandidate,
    AF3SequenceCandidate,
    MoleculeCandidate,
    PseudoCandidate,
    QueryResult,
    SequenceCandidate,
)
from .registry import MHCRegistry

__all__ = [
    "MHCRegistry",
    "MHCTransformError",
    "AF3MoleculeCandidate",
    "AF3SequenceCandidate",
    "MoleculeCandidate",
    "PseudoCandidate",
    "QueryResult",
    "SequenceCandidate",
]

__version__ = "0.2.0"

