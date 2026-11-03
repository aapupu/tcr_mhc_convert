"""Package-specific exceptions."""


class TCRTransformError(ValueError):
    """A TCR cannot be transformed without making an unsafe assumption."""


class TCRReferenceError(TCRTransformError):
    """A bundled or user-supplied TCR reference is invalid or incomplete."""


class TCRParseError(TCRTransformError):
    """A full TCR chain cannot be parsed unambiguously."""

