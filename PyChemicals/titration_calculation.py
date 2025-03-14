"""
Titration Calculation Module

This module provides classes and functions for titration calculations,
including analyte and titrant handling as well as pH calculations at various
points in the titration process.
"""

import numpy as np
from .chemicals import Acid, Base, Chemical
from .chemicals_db import get_valid_acids, get_valid_bases

valid_acids = get_valid_acids()
valid_bases = get_valid_bases()


def _validate_positive(value, name: str):
    if value is None:
        raise ValueError(f"{name} must be provided.")
    if value <= 0:
        raise ValueError(f"{name} must be a positive number.")


class Analyte(Chemical):
    """Class representing an analyte in a titration process."""

    def __new__(
        cls, name: str, volume: float = None, concentration: float = None, titrant=None
    ):
        # Choose the appropriate subclass based on the name.
        if name in valid_acids:
            return Acid(
                name,
                properties=valid_acids[name],
                concentration=concentration,
                volume=volume,
            )
        if name in valid_bases:
            return Base(
                name,
                properties=valid_bases[name],
                concentration=concentration,
                volume=volume,
            )
        return super().__new__(cls)

    def __init__(
        self, name: str, volume: float = None, concentration: float = None, titrant=None
    ):
        _validate_positive(volume, "Volume")
        _validate_positive(concentration, "Concentration")
        # Call Chemical's initializer
        super().__init__(name=name, concentration=concentration, volume=volume)
        if titrant is not None:
            self.validate_titrant(titrant)

    def validate_titrant(self, titrant: Chemical):
        """
        Validate that the titrant is chemically compatible with the analyte.
        Raises a ValueError if both the analyte and titrant are of the same type.
        """
        if isinstance(titrant, Acid) and isinstance(self, Acid):
            raise ValueError("Both Analyte and Titrant cannot be acids.")
        if isinstance(titrant, Base) and isinstance(self, Base):
            raise ValueError("Both Analyte and Titrant cannot be bases.")


class Titrant(Chemical):
    """Class representing a titrant in a titration process."""

    def __new__(
        cls, name: str, concentration: float = None, volume: float = None, analyte=None
    ):
        # Choose the appropriate subclass based on the name.
        if name in valid_acids:
            return Acid(
                name,
                properties=valid_acids[name],
                concentration=concentration,
                volume=volume,
            )
        if name in valid_bases:
            return Base(
                name,
                properties=valid_bases[name],
                concentration=concentration,
                volume=volume,
            )
        return super().__new__(cls)

    def __init__(
        self, name: str, concentration: float = None, volume: float = None, analyte=None
    ):
        _validate_positive(volume, "Volume")
        _validate_positive(concentration, "Concentration")
        super().__init__(name=name, concentration=concentration, volume=volume)
        if analyte is not None:
            self.validate_analyte(analyte)

    def validate_analyte(self, analyte: Chemical):
        """
        Validate that the analyte is chemically compatible with the titrant.
        Raises a ValueError if both the analyte and titrant are of the same type.
        """
        if isinstance(analyte, Acid) and isinstance(self, Acid):
            raise ValueError("Both Analyte and Titrant cannot be acids.")
        if isinstance(analyte, Base) and isinstance(self, Base):
            raise ValueError("Both Analyte and Titrant cannot be bases.")


def strength(analyte: Chemical) -> str:
    """
    Determine the strength of the analyte.
    Returns a string such as "strong acid", "weak acid", "strong base", "weak base",
    or None if the analyte is neither an acid nor a base.
    """
    if isinstance(analyte, Acid):
        if hasattr(analyte, "ka") and analyte.ka > 1:
            return "strong acid"
        return "weak acid"
    if isinstance(analyte, Base):
        if hasattr(analyte, "kb") and analyte.kb > 1:
            return "strong base"
        return "weak base"
    return None


def volume_of_titrant(
    analyte: Chemical, titrant: Chemical, stoichiometric_ratio: float = 1.0
) -> float:
    """
    Calculate the volume of titrant needed to reach equivalence.

    Args:
        analyte: The analyte being titrated.
        titrant: The titrant used in the titration.
        stoichiometric_ratio: The stoichiometric ratio of the analyte to titrant.

    Returns:
        The volume of titrant needed in liters.

    Raises:
        ValueError: If required concentrations, volumes, or
        stoichiometric ratio are missing or invalid.
    """
    _validate_positive(analyte.volume, "Analyte volume")
    _validate_positive(analyte.concentration, "Analyte concentration")
    _validate_positive(titrant.concentration, "Titrant concentration")
    if stoichiometric_ratio <= 0:
        raise ValueError("Stoichiometric ratio must be greater than zero.")

    moles_analyte = analyte.concentration * analyte.volume
    moles_titrant_needed = moles_analyte / stoichiometric_ratio

    # If the analyte is diprotic/dibasic, adjust based on proticity (if defined)
    proticity = getattr(analyte, "proticity", 1)
    if proticity == 2 and stoichiometric_ratio == 2.0:
        moles_titrant_needed *= 2

    volume_titrant_needed = moles_titrant_needed / titrant.concentration
    return volume_titrant_needed


def calculate_analyte_concentration(
    analyte: Chemical, titrant: Chemical, stoichiometric_ratio: float
) -> float:
    """
    Calculate the concentration of the analyte from titration data.

    Args:
        analyte: The analyte object.
        titrant: The titrant object.
        stoichiometric_ratio: Stoichiometric ratio between the analyte and titrant.

    Returns:
        The calculated concentration of the analyte (in mol/L).

    Raises:
        ValueError: If volumes, concentrations, or stoichiometric ratio are invalid.
    """
    _validate_positive(analyte.volume, "Analyte volume")
    _validate_positive(titrant.volume, "Titrant volume")
    _validate_positive(titrant.concentration, "Titrant concentration")
    if stoichiometric_ratio <= 0:
        raise ValueError("Stoichiometric ratio must be greater than zero.")

    moles_titrant = titrant.concentration * titrant.volume
    moles_analyte = moles_titrant * stoichiometric_ratio
    concentration_analyte = moles_analyte / analyte.volume
    return concentration_analyte


def ph_at_equivalence(analyte: Chemical, titrant: Chemical) -> float:
    """
    Calculate the pH at the equivalence point of the titration.

    Note: This is a simplified placeholder. For weak acid/strong base (or vice versa)
    titrations, the equivalence pH will differ from 7.

    Returns:
        pH at equivalence or None if not applicable.
    """
    if isinstance(analyte, Acid) and isinstance(titrant, Base):
        return 7.0  # For a strong acid + strong base titration.
    if isinstance(analyte, Base) and isinstance(titrant, Acid):
        return 7.0
    return None


def calculate_half_equivalence(analyte: Chemical, titrant: Chemical):
    """
    Calculate the pH at the half-equivalence point.

    Raises:
        NotImplementedError: This function is not implemented yet.
    """
    raise NotImplementedError("calculate_half_equivalence is not implemented yet.")


def _handle_acid_ph(
    acid: Acid, base: Base, excess_volume: float, total_volume: float
) -> float:
    """Helper function to calculate pH for acids."""
    moles_acid = acid.concentration * acid.volume
    if hasattr(acid, "ka") and acid.ka < 1:
        if excess_volume == acid.volume / 2:
            return acid.pka() if hasattr(acid, "pka") else -np.log10(acid.ka)

        moles_base_added = base.concentration * excess_volume
        remaining_acid = moles_acid - moles_base_added
        if remaining_acid <= 0:
            excess_oh = (moles_base_added - moles_acid) / total_volume
            if excess_oh <= 0:
                raise ValueError("Invalid OH concentration computed.")
            return 14 + np.log10(excess_oh)

        conjugate_base = moles_base_added
        conc_acid = remaining_acid / total_volume
        conc_base = conjugate_base / total_volume
        if conc_acid <= 0 or conc_base <= 0:
            raise ValueError("Invalid concentrations for Henderson-Hasselbalch.")
        return -np.log10(acid.ka) + np.log10(conc_base / conc_acid)

    moles_base = base.concentration * excess_volume
    if moles_base > moles_acid:
        excess_oh = (moles_base - moles_acid) / total_volume
        if excess_oh <= 0:
            raise ValueError("Calculated excess OH concentration is non-positive.")
        return 14 + np.log10(excess_oh)
    if moles_acid > moles_base:
        excess_h = (moles_acid - moles_base) / total_volume
        if excess_h <= 0:
            raise ValueError("Calculated excess H+ concentration is non-positive.")
        return -np.log10(excess_h)
    return ph_at_equivalence(acid, base)


def ph_after_equivalence(acid: Acid, base: Base, excess_volume: float) -> float:
    """Calculate pH after equivalence point."""
    _validate_positive(acid.volume, "Acid volume")
    _validate_positive(acid.concentration, "Acid concentration")
    _validate_positive(base.concentration, "Base concentration")
    if excess_volume < 0:
        raise ValueError("Excess volume cannot be negative.")

    total_volume = acid.volume + excess_volume
    return _handle_acid_ph(acid, base, excess_volume, total_volume)


# Example usage of the functions
if __name__ == "__main__":
    hcl = Analyte("Hydrochloric Acid", concentration=0.1, volume=0.025)
    naoh = Titrant("Sodium Hydroxide", concentration=0.1, volume=0.025)

    print("pH at Equivalence:", ph_at_equivalence(hcl, naoh))
    print("Volume of Titrant Needed:", volume_of_titrant(hcl, naoh))
    print("pH after Equivalence:", ph_after_equivalence(hcl, naoh, 0.01))
