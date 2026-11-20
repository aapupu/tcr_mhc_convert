"""Peptide-binding-region representation identifiers and validation."""

from __future__ import annotations

from .errors import MHCTransformError


PBR_I_DEFINITION = "mhci-alpha1-alpha2-v1"
PBR_II_ALPHA_DEFINITION = "mhcii-alpha1-v1"
PBR_II_BETA_DEFINITION = "mhcii-beta1-v1"
PBR_DEFINITIONS = (
    PBR_I_DEFINITION,
    PBR_II_ALPHA_DEFINITION,
    PBR_II_BETA_DEFINITION,
)


def pbr_definition(mhc_class: str, chain_role: str) -> str:
    """Return the only valid PBR definition for a class/chain combination."""

    if mhc_class == "I" and chain_role == "heavy":
        return PBR_I_DEFINITION
    if mhc_class == "II" and chain_role == "alpha":
        return PBR_II_ALPHA_DEFINITION
    if mhc_class == "II" and chain_role == "beta":
        return PBR_II_BETA_DEFINITION
    raise MHCTransformError(
        f"no peptide-binding-region definition for MHC-{mhc_class} {chain_role}"
    )
