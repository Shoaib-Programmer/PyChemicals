"""Advanced tests for titration calculations."""

import pytest
from PyChemicals.titration_calculation import (
    Analyte,
    Titrant,
    volume_of_titrant,
    calculate_analyte_concentration,
    ph_at_equivalence,
    ph_after_equivalence,
)
from PyChemicals.chemicals import Acid


def test_titration_stoichiometry():
    """Test stoichiometric calculations in titrations."""
    analyte = Analyte("Hydrogen chloride", volume=0.025, concentration=0.1)
    titrant = Titrant("Sodium hydroxide", concentration=0.1)

    # Test 1:1 stoichiometry
    v_eq = volume_of_titrant(analyte, titrant)
    assert abs(v_eq - 0.025) < 0.001

    # Test 2:1 stoichiometry
    v_eq = volume_of_titrant(analyte, titrant, stoichiometric_ratio=2.0)
    assert abs(v_eq - 0.0125) < 0.001


def test_titration_curves():
    """Test various points on titration curves."""
    acid = Analyte("Hydrogen chloride", volume=0.025, concentration=0.1)
    base = Titrant("Sodium hydroxide", concentration=0.1)

    # Before equivalence (half-way)
    v_eq = volume_of_titrant(acid, base)
    v_before = 0.5 * v_eq  # Use 50% of equivalence volume
    ph = ph_after_equivalence(acid, base, v_before)  # Pass actual volume added
    assert 1.0 < ph < 7.0

    # At equivalence
    ph_eq = ph_at_equivalence(acid, base)
    assert abs(ph_eq - 7.0) < 0.1

    # After equivalence
    v_after = 1.5 * v_eq  # Use 150% of equivalence volume
    ph = ph_after_equivalence(acid, base, v_after)  # Pass actual volume added
    assert ph > 7.0


def test_concentration_determination():
    """Test determination of unknown concentrations."""
    unknown = Analyte("Hydrogen chloride", volume=0.025)
    titrant = Titrant("Sodium hydroxide", concentration=0.1, volume=0.025)

    conc = calculate_analyte_concentration(unknown, titrant, 1.0)
    assert abs(conc - 0.1) < 0.001


def test_invalid_titrations():
    """Test invalid titration scenarios."""
    acid1 = Analyte("Hydrogen chloride")
    acid2 = Titrant("Sulfuric acid")

    with pytest.raises(ValueError):
        # Cannot titrate acid with acid
        volume_of_titrant(acid1, acid2)


def test_diprotic_titrations():
    """Test titrations involving diprotic acids."""
    acid = Analyte("Sulfuric acid", concentration=0.1, volume=0.025)
    base = Titrant("Sodium hydroxide", concentration=0.1)

    # First equivalence point (1:1 ratio)
    v_eq1 = volume_of_titrant(acid, base)
    assert abs(v_eq1 - 0.025) < 0.001

    # Second equivalence point (1:2 ratio, requires double the volume)
    v_eq2 = volume_of_titrant(
        acid, base, stoichiometric_ratio=0.5
    )  # Changed from 2.0 to 0.5
    assert abs(v_eq2 - 0.05) < 0.001  # Should be double v_eq1


def test_buffer_regions():
    """Test pH calculations in buffer regions."""
    weak_acid = Analyte("Ethanoic acid", volume=0.025, concentration=0.1)
    strong_base = Titrant("Sodium hydroxide", concentration=0.1)

    # At half equivalence (buffer region)
    v_eq = volume_of_titrant(weak_acid, strong_base)
    v_half = v_eq / 2  # Half-way to equivalence
    ph = ph_after_equivalence(
        weak_acid, strong_base, v_half
    )  # Pass actual volume added

    # For weak acid buffer, pH ≈ pKa at half-neutralization
    if isinstance(weak_acid, Acid):
        assert abs(ph - weak_acid.pka()) < 0.2
