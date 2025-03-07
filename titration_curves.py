"""
This module plots titration curves for monoprotic acid-base titrations.
"""

import matplotlib  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
from chemicals import Acid, Base
from titration_calculation import volume_of_titrant, Analyte, Titrant

matplotlib.use("tkAgg")  # or use 'Qt5Agg'


def calculate_ph_acid(
    acid_conc: float, base_conc: float, acid_volume: float, v_titrant: float
) -> float:
    """
    Calculate pH for a strong acid titrated with a strong base.

    Args:
        acid_conc: Concentration of the acid.
        base_conc: Concentration of the base.
        acid_volume: Volume of the acid.
        v_titrant: Volume of base added.

    Returns:
        Calculated pH value.
    """
    moles_acid = acid_conc * acid_volume
    moles_base = base_conc * v_titrant
    total_volume = acid_volume + v_titrant
    if moles_base < moles_acid:
        h_plus = (moles_acid - moles_base) / total_volume
        return -np.log10(h_plus)
    if moles_base == moles_acid:
        return 7.0
    oh_minus = (moles_base - moles_acid) / total_volume
    return 14 + np.log10(oh_minus)


def calculate_ph_base(
    base_conc: float, acid_conc: float, base_volume: float, v_titrant: float
) -> float:
    """
    Calculate pH for a strong base titrated with a strong acid.

    Args:
        base_conc: Concentration of the base.
        acid_conc: Concentration of the acid.
        base_volume: Volume of the base.
        v_titrant: Volume of acid added.

    Returns:
        Calculated pH value.
    """
    moles_base = base_conc * base_volume
    moles_acid = acid_conc * v_titrant
    total_volume = base_volume + v_titrant
    if moles_acid < moles_base:
        oh_minus = (moles_base - moles_acid) / total_volume
        return 14 + np.log10(oh_minus)
    if moles_acid == moles_base:
        return 7.0
    h_plus = (moles_acid - moles_base) / total_volume
    return -np.log10(h_plus)


def titrate_curve_monoprotic(analyte: Analyte, titrant: Titrant):
    """
    Plots the titration curve for a monoprotic acid-base titration.

    Parameters:
        analyte (Analyte): The acid or base being titrated.
        titrant (Titrant): The titrant added to the analyte.

    Raises:
        ValueError: If any required property (volume, concentration) is missing or if
                    the analyte is not an Acid or Base.
    """
    if analyte.volume is None:
        raise ValueError("Analyte must have a defined volume.")
    if analyte.concentration is None:
        raise ValueError("Analyte must have a concentration.")
    if titrant.concentration is None:
        raise ValueError("Titrant must have a concentration.")
    if titrant.volume is None:
        raise ValueError("Titrant must have a defined volume.")

    # Calculate volume of titrant needed to reach the equivalence point
    v_eqv = volume_of_titrant(analyte, titrant)
    # Define volumes to add titrant (from 0 to twice the equivalence volume)
    v_added = np.linspace(0, 2 * v_eqv, 500)

    if isinstance(analyte, Acid):
        acid_conc = analyte.concentration
        base_conc = titrant.concentration
        acid_volume = analyte.volume
        ph_values = [
            calculate_ph_acid(acid_conc, base_conc, acid_volume, v) for v in v_added
        ]
        label = (
            f"{acid_volume * 1000:.0f} mL {analyte.name} ({acid_conc} M) titrated with "
            f"{titrant.name} ({base_conc} M)"
        )
    elif isinstance(analyte, Base):
        base_conc = analyte.concentration
        acid_conc = titrant.concentration
        base_volume = analyte.volume
        ph_values = [
            calculate_ph_base(base_conc, acid_conc, base_volume, v) for v in v_added
        ]
        label = (
            f"{base_volume * 1000:.0f} mL {analyte.name} ({base_conc} M) titrated with "
            f"{titrant.name} ({acid_conc} M)"
        )
    else:
        raise ValueError("Analyte must be either an Acid or a Base.")

    plt.figure(figsize=(8, 5))
    plt.plot(v_added * 1000, ph_values, label=label)
    plt.axvline(x=v_eqv * 1000, color="r", linestyle="--", label="Equivalence Point")
    plt.axhline(y=7, color="g", linestyle="--", label="pH = 7 at Equivalence")
    plt.xlabel(f"Volume of {titrant.name} added (mL)")
    plt.ylabel("pH")
    plt.title("Titration Curve")
    plt.legend()
    plt.grid(True)
    plt.show()
