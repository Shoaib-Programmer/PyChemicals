"""Test module for titration calculations."""

from PyChemicals.chemicals import Acid, Base
from PyChemicals.titration_calculation import (
    volume_of_titrant,
    calculate_analyte_concentration,
    ph_at_equivalence,
    ph_after_equivalence,
)
from .conftest import hcl_properties, naoh_properties


def test_volume_of_titrant():
    """Test calculation of required titrant volume."""
    hcl = Acid(
        "Hydrogen chloride", properties=hcl_properties, concentration=0.1, volume=0.025
    )
    naoh = Base("Sodium hydroxide", properties=naoh_properties, concentration=0.1)
    volume = volume_of_titrant(hcl, naoh)
    assert (
        abs(volume - 0.025) < 0.001
    )  # Should need equal volumes for equal concentrations


def test_ph_at_equivalence():
    """Test pH calculation at equivalence point."""
    hcl = Acid(
        "Hydrogen chloride", properties=hcl_properties, concentration=0.1, volume=0.025
    )
    naoh = Base("Sodium hydroxide", properties=naoh_properties, concentration=0.1)
    ph = ph_at_equivalence(hcl, naoh)
    assert abs(ph - 7.0) < 0.1  # Should be neutral at equivalence


def test_ph_after_equivalence():
    """Test pH calculation after equivalence point."""
    hcl = Acid(
        "Hydrogen chloride", properties=hcl_properties, concentration=0.1, volume=0.025
    )
    naoh = Base("Sodium hydroxide", properties=naoh_properties, concentration=0.1)

    # Test after equivalence - use 120% of equivalence volume
    v_eq = volume_of_titrant(hcl, naoh)
    v_after = 1.2 * v_eq  # 20% excess
    ph = ph_after_equivalence(hcl, naoh, v_after)
    assert ph > 7.0  # Should be basic after equivalence


def test_analyte_concentration():
    """Test calculation of analyte concentration."""
    unknown_acid = Acid("Hydrogen chloride", properties=hcl_properties, volume=0.025)
    naoh = Base(
        "Sodium hydroxide", properties=naoh_properties, concentration=0.1, volume=0.025
    )
    conc = calculate_analyte_concentration(unknown_acid, naoh, 1.0)
    assert abs(conc - 0.1) < 0.001  # Should be 0.1M
