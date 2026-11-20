"""Auditable human/mouse MHC exact, partial and imputed conversion."""

from .errors import MHCTransformError
from .models import (
    AF3MoleculeCandidate,
    AF3SequenceCandidate,
    FallbackAF3SequenceCandidate,
    FallbackPBRSequenceCandidate,
    FallbackSequenceCandidate,
    MoleculeCandidate,
    PBRMoleculeCandidate,
    PBRSequenceCandidate,
    PseudoCandidate,
    QueryResult,
    SequenceCandidate,
)
from .registry import MHCRegistry
from .pbr import (
    PBR_DEFINITIONS,
    PBR_I_DEFINITION,
    PBR_II_ALPHA_DEFINITION,
    PBR_II_BETA_DEFINITION,
)

__all__ = [
    "MHCRegistry",
    "MHCTransformError",
    "AF3MoleculeCandidate",
    "AF3SequenceCandidate",
    "FallbackAF3SequenceCandidate",
    "FallbackPBRSequenceCandidate",
    "FallbackSequenceCandidate",
    "MoleculeCandidate",
    "PBRMoleculeCandidate",
    "PBRSequenceCandidate",
    "PseudoCandidate",
    "QueryResult",
    "SequenceCandidate",
    "PBR_DEFINITIONS",
    "PBR_I_DEFINITION",
    "PBR_II_ALPHA_DEFINITION",
    "PBR_II_BETA_DEFINITION",
]

__version__ = "0.5.0"

