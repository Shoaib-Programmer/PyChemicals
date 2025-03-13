"""
This module plots titration curves for monoprotic acid-base titrations.
"""

import matplotlib  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from .chemicals import Acid, Base
from .titration_calculation import volume_of_titrant, Analyte, Titrant

# Set the matplotlib backend (adjust or remove as needed)
matplotlib.use("TkAgg")


def calculate_ph_strong(moles_primary, moles_titrant, total_volume) -> np.ndarray:
    """
    Vectorized pH calculation for strong acid or base titrations.

    Args:
        moles_primary (float or np.ndarray): Moles of the analyte (acid or base).
        moles_titrant (float or np.ndarray): Moles of titrant added.
        total_volume (float or np.ndarray): Total volume of the solution.

    Returns:
        np.ndarray: Calculated pH values.
    """
    delta = moles_primary - moles_titrant
    # For acid excess: pH = -log10([H+]); for base excess: pH = 14 + log10([OH-]);
    # At equivalence (delta == 0), pH is set to 7.
    return np.where(
        delta > 0,
        -np.log10(delta / total_volume),
        np.where(delta < 0, 14 + np.log10((-delta) / total_volume), 7.0),
    )


def calculate_ph_acid(
    acid_conc: float, base_conc: float, acid_volume: float, v_titrant
) -> np.ndarray:
    """
    Calculate pH for a strong acid titrated with a strong base using vectorized computation.

    Args:
        acid_conc (float): Concentration of the acid.
        base_conc (float): Concentration of the base.
        acid_volume (float): Volume of the acid.
        v_titrant (float or np.ndarray): Volume(s) of base added.

    Returns:
        np.ndarray: Array of calculated pH values.
    """
    moles_acid = acid_conc * acid_volume
    moles_base = base_conc * np.array(v_titrant)
    total_volume = acid_volume + np.array(v_titrant)
    return calculate_ph_strong(moles_acid, moles_base, total_volume)


def calculate_ph_base(
    base_conc: float, acid_conc: float, base_volume: float, v_titrant
) -> np.ndarray:
    """
    Calculate pH for a strong base titrated with a strong acid using vectorized computation.

    Args:
        base_conc (float): Concentration of the base.
        acid_conc (float): Concentration of the acid.
        base_volume (float): Volume of the base.
        v_titrant (float or np.ndarray): Volume(s) of acid added.

    Returns:
        np.ndarray: Array of calculated pH values.
    """
    moles_base = base_conc * base_volume
    moles_acid = acid_conc * np.array(v_titrant)
    total_volume = base_volume + np.array(v_titrant)
    return calculate_ph_strong(moles_base, moles_acid, total_volume)


def titrate_curve_monoprotic(analyte: Analyte, titrant: Titrant) -> None:
    """
    Plots the titration curve for a monoprotic acid-base titration.

    Args:
        analyte (Analyte): The acid or base being titrated.
        titrant (Titrant): The titrant added to the analyte.

    Raises:
        ValueError: If any required property (volume, concentration) is missing or if
                    the analyte is not an Acid or a Base.
    """
    # Validate necessary properties
    if analyte.volume is None:
        raise ValueError("Analyte must have a defined volume.")
    if analyte.concentration is None:
        raise ValueError("Analyte must have a concentration.")
    if titrant.concentration is None:
        raise ValueError("Titrant must have a concentration.")
    if titrant.volume is None:
        raise ValueError("Titrant must have a defined volume.")

    # Calculate the volume of titrant needed to reach equivalence
    v_eqv = volume_of_titrant(analyte, titrant)
    # Define a range of titrant volumes from 0 to twice the equivalence volume
    v_added = np.linspace(0, 2 * v_eqv, 500)

    if isinstance(analyte, Acid):
        acid_conc = analyte.concentration
        base_conc = titrant.concentration
        acid_volume = analyte.volume
        ph_values = calculate_ph_acid(acid_conc, base_conc, acid_volume, v_added)
        label = (
            f"{acid_volume * 1000:.0f} mL {analyte.name} ({acid_conc} M) titrated with "
            f"{titrant.name} ({base_conc} M)"
        )
    elif isinstance(analyte, Base):
        base_conc = analyte.concentration
        acid_conc = titrant.concentration
        base_volume = analyte.volume
        ph_values = calculate_ph_base(base_conc, acid_conc, base_volume, v_added)
        label = (
            f"{base_volume * 1000:.0f} mL {analyte.name} ({base_conc} M) titrated with "
            f"{titrant.name} ({acid_conc} M)"
        )
    else:
        raise ValueError("Analyte must be either an Acid or a Base.")

    # Plot the titration curve
    plt.figure(figsize=(8, 5))
    plt.plot(v_added * 1000, ph_values, label=label, color="blue")
    plt.axvline(x=v_eqv * 1000, color="red", linestyle="--", label="Equivalence Point")
    plt.axhline(y=7, color="green", linestyle="--", label="pH = 7")
    plt.xlabel(f"Volume of {titrant.name} added (mL)")
    plt.ylabel("pH")
    plt.title("Monoprotic Acid-Base Titration Curve")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
