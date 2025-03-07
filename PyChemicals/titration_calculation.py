"""
Titration Calculation Module

This module provides classes and functions for titration calculations,
including analyte and titrant handling as well as pH calculations at various
points in the titration process.
"""

import numpy as np
from .chemicals import Acid, Base, Chemical
from .predefined_chemicals import valid_acids, valid_bases


class Analyte(Chemical):
    """Class representing an analyte in a titration process."""

    def __init__(
        self, name: str, volume: float = None, concentration: float = None, titrant=None
    ):
        """
        Initialize an Analyte instance.

        Args:
            name (str): Name of the analyte.
            volume (float, optional): Volume of the analyte.
            concentration (float, optional): Concentration of the analyte.
            titrant (optional): Titrant used in the titration.
        """
        super().__init__(name=name, concentration=concentration, volume=volume)
        if name in valid_acids:
            self.__class__ = Acid
            self.__init__(name=name, concentration=concentration, volume=volume)
        elif name in valid_bases:
            self.__class__ = Base
            self.__init__(name=name, concentration=concentration, volume=volume)

        if titrant is not None:
            self.validate_titrant(titrant)

    def validate_titrant(self, titrant):
        """
        Validate that the titrant is chemically compatible with the analyte.

        Args:
            titrant: The titrant to be validated.

        Raises:
            ValueError: If both analyte and titrant are acids or both are bases.
        """
        if isinstance(titrant, Acid) and isinstance(self, Acid):
            raise ValueError("Both Analyte and Titrant cannot be acids.")
        if isinstance(titrant, Base) and isinstance(self, Base):
            raise ValueError("Both Analyte and Titrant cannot be bases.")


class Titrant(Chemical):
    """Class representing a titrant in a titration process."""

    def __init__(
        self, name: str, concentration: float = None, volume: float = None, analyte=None
    ):
        """
        Initialize a Titrant instance.

        Depending on the name provided, the instance's class type may change to Acid
        or Base.

        Args:
            name (str): Name of the titrant.
            concentration (float, optional): Concentration of the titrant.
            volume (float, optional): Volume of the titrant.
            analyte (optional): Analyte for cross-validation.
        """
        if name in valid_acids:
            super().__init__(name=name, concentration=concentration, volume=volume)
            self.__class__ = Acid  # Change the class type to Acid
        elif name in valid_bases:
            super().__init__(name=name, concentration=concentration, volume=volume)
            self.__class__ = Base  # Change the class type to Base
        else:
            super().__init__(name=name, concentration=concentration, volume=volume)

        if analyte is not None:
            self.validate_analyte(analyte)

    def validate_analyte(self, analyte):
        """
        Validate that the analyte is chemically compatible with the titrant.

        Args:
            analyte: The analyte to be validated.

        Raises:
            ValueError: If both analyte and titrant are acids or both are bases.
        """
        if isinstance(analyte, Acid) and isinstance(self, Acid):
            raise ValueError("Both Analyte and Titrant cannot be acids.")
        if isinstance(analyte, Base) and isinstance(self, Base):
            raise ValueError("Both Analyte and Titrant cannot be bases.")


def strength(analyte: Analyte):
    """
    Determine the strength of the analyte.

    Args:
        analyte (Analyte): The analyte whose strength is to be determined.

    Returns:
        str or None: "strong acid", "weak acid", "strong base", "weak base", or
        None if the analyte is neither an acid nor a base.
    """
    if isinstance(analyte, Acid):
        if analyte.ka > 1:
            return "strong acid"
        return "weak acid"
    if isinstance(analyte, Base):
        if analyte.kb > 1:
            return "strong base"
        return "weak base"
    return None


def volume_of_titrant(
    analyte: Analyte, titrant: Titrant, stoichiometric_ratio: float = 1.0
) -> float:
    """
    Calculate the volume of titrant needed to reach equivalence.

    Args:
        analyte (Analyte): The analyte being titrated.
        titrant (Titrant): The titrant used in the titration.
        stoichiometric_ratio (float): The stoichiometric ratio of the analyte to titrant.

    Returns:
        float: Volume of titrant needed in liters.

    Raises:
        ValueError: If the analyte's or titrant's concentration/volume is not provided.
    """
    if analyte.concentration is None or analyte.volume is None:
        raise ValueError("Analyte concentration and volume must be provided.")
    if titrant.concentration is None:
        raise ValueError("Titrant concentration must be provided.")

    moles_analyte = analyte.concentration * analyte.volume
    moles_titrant_needed = moles_analyte / stoichiometric_ratio

    # For diprotic acids/bases, adjust the volume based on the dissociation step
    if isinstance(analyte, (Acid, Base)) and analyte.proticity == 2:
        if stoichiometric_ratio == 2.0:
            moles_titrant_needed *= (
                2  # Double the moles needed for complete neutralization
            )

    volume_titrant_needed = moles_titrant_needed / titrant.concentration
    return volume_titrant_needed


def calculate_analyte_concentration(
    analyte: Analyte, titrant: Titrant, stoichiometric_ratio: float
):
    """
    Calculate the concentration of the analyte using the titrant's data.

    Args:
        analyte (Analyte): The analyte object.
        titrant (Titrant): The titrant object.
        stoichiometric_ratio (float): Stoichiometric ratio between the analyte and titrant.

    Returns:
        float: The calculated concentration of the analyte (in mol/L).

    Raises:
        ValueError: If required volumes or concentrations are missing or invalid.
    """
    if analyte.volume is None:
        raise ValueError("Volume of analyte must be provided.")
    if titrant.volume is None:
        raise ValueError("Volume of titrant must be provided.")
    if titrant.concentration is None:
        raise ValueError("Concentration of titrant must be provided.")
    if stoichiometric_ratio <= 0:
        raise ValueError("Stoichiometric ratio must be greater than zero.")

    moles_titrant = titrant.concentration * titrant.volume
    moles_analyte = moles_titrant * stoichiometric_ratio
    concentration_analyte = moles_analyte / analyte.volume

    return concentration_analyte


def ph_at_equivalence(analyte: Analyte, titrant: Titrant):
    """
    Calculate the pH at the equivalence point of the titration.

    This is a placeholder for actual pH calculation logic based on the types of
    analyte and titrant.

    Args:
        analyte (Analyte): The analyte being titrated.
        titrant (Titrant): The titrant used in the titration.

    Returns:
        float or None: pH at equivalence or None if not applicable.
    """
    if isinstance(analyte, Acid) and isinstance(titrant, Base):
        return 7.0  # Neutralization point for strong acid + strong base
    if isinstance(analyte, Base) and isinstance(titrant, Acid):
        return 7.0  # Neutralization point for strong base + strong acid
    return None


def calculate_half_equivalence(analyte, titrant):
    """
    Calculate the pH at the half-equivalence point.

    Args:
        analyte: The analyte object.
        titrant: The titrant object.

    Raises:
        NotImplementedError: This function is not implemented yet.
    """
    raise NotImplementedError("calculate_half_equivalence is not implemented yet.")


def ph_after_equivalence(acid: Acid, base: Base, excess_volume: float) -> float:
    """
    Calculate pH at any point during titration.

    Args:
        acid (Acid): The acid in the titration.
        base (Base): The base titrant.
        excess_volume (float): Volume of base added (in L).

    Returns:
        float: The calculated pH.
    """
    moles_acid = acid.concentration * acid.volume

    # For weak acids, use Ka to calculate the pH
    if acid.ka < 1:  # Weak acid
        if excess_volume == acid.volume / 2:  # At half-equivalence
            return acid.pka()

        # Calculate concentration of remaining acid and conjugate base
        total_volume = acid.volume + excess_volume
        moles_base_added = base.concentration * excess_volume
        remaining_acid = (moles_acid - moles_base_added) / total_volume
        conjugate_base = moles_base_added / total_volume
        if remaining_acid > 0:
            return -np.log10((acid.ka * remaining_acid) / conjugate_base)

    # For strong acids and bases, use excess concentration
    moles_base = base.concentration * excess_volume
    total_volume = acid.volume + excess_volume

    if moles_base > moles_acid:  # After equivalence
        excess_oh = (moles_base - moles_acid) / total_volume
        return 14 + np.log10(excess_oh)

    # Before equivalence
    excess_h = (moles_acid - moles_base) / total_volume
    return -np.log10(excess_h)


# Example usage of the functions
if __name__ == "__main__":
    hcl = Acid("Hydrochloric Acid", concentration=0.1, volume=0.025)  # 0.1 M HCl, 25 mL
    naoh = Base(
        "Sodium Hydroxide", concentration=0.1, volume=0.025
    )  # 0.1 M NaOH, 25 mL

    print("pH at Equivalence:", ph_at_equivalence(hcl, naoh))
    print("Volume of Titrant Needed:", volume_of_titrant(hcl, naoh))
    # Example call for ph_after_equivalence with an excess volume of 0.01 L
    print("pH after Equivalence:", ph_after_equivalence(hcl, naoh, 0.01))
