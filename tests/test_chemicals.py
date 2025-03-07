"""Test module for chemical classes including Chemical, Acid, and Base."""

import pytest
from PyChemicals.chemicals import Chemical, Acid, Base

# Chemical class tests
def test_chemical_basic_properties():
    """Test basic properties of Chemical class."""
    chem = Chemical("Test", concentration=0.1, volume=0.5, mass=2.0)
    assert chem.name == "Test"
    assert chem.concentration == 0.1
    assert chem.volume == 0.5
    assert chem.mass == 2.0

def test_chemical_validation():
    """Test validation of chemical properties."""
    with pytest.raises(ValueError, match="Concentration must be positive"):
        Chemical("Test", concentration=-1)
    with pytest.raises(ValueError, match="Volume must be positive"):
        Chemical("Test", volume=-1)
    with pytest.raises(ValueError, match="Mass must be positive"):
        Chemical("Test", mass=-1)

def test_chemical_moles_calculation():
    """Test moles calculation."""
    chem = Chemical("Test", concentration=0.1, volume=0.5)
    assert chem.calculate_moles() == 0.05

# Acid class tests
def test_acid_strength_classification(strong_acid, weak_acid):
    """Test acid strength classification."""
    assert strong_acid.ka > 1  # Strong acid
    assert weak_acid.ka < 1  # Weak acid

def test_acid_ph_calculations():
    """Test pH calculations for different acid concentrations."""
    acids = [
        Acid("Hydrochloric Acid", concentration=c) 
        for c in [0.1, 0.01, 0.001]
    ]
    expected_phs = [1.0, 2.0, 3.0]
    for acid, expected_ph in zip(acids, expected_phs):
        assert abs(acid.ph() - expected_ph) < 0.1

def test_acid_edge_cases():
    """Test edge cases for acids."""
    # Very dilute acid
    dilute = Acid("Hydrochloric Acid", concentration=1e-7)
    assert dilute.ph() > 6.9  # Should be close to neutral

    # Very concentrated acid
    conc = Acid("Hydrochloric Acid", concentration=12.0)
    assert conc.ph() < 0  # Negative pH is possible for concentrated acids

def test_diprotic_acid_properties(diprotic_acid):
    """Test properties specific to diprotic acids."""
    assert diprotic_acid.proticity == 2
    # For sulfuric acid: ka1 = 1e3, ka2 = 1.2e-2
    assert diprotic_acid.ka1 == 1e3
    assert diprotic_acid.ka2 == 1.2e-2
    assert diprotic_acid.ka1 > diprotic_acid.ka2

# Base class tests
def test_base_strength_classification(strong_base, weak_base):
    """Test base strength classification."""
    assert strong_base.kb > 1  # Strong base
    assert weak_base.kb < 1  # Weak base

def test_base_ph_calculations():
    """Test pH calculations for different base concentrations."""
    bases = [
        Base("Sodium Hydroxide", concentration=c)
        for c in [0.1, 0.01, 0.001]
    ]
    expected_phs = [13.0, 12.0, 11.0]
    for base, expected_ph in zip(bases, expected_phs):
        assert abs(base.ph() - expected_ph) < 0.1

def test_base_edge_cases():
    """Test edge cases for bases."""
    # Very dilute base
    dilute = Base("Sodium Hydroxide", concentration=1e-7)
    assert dilute.ph() < 7.1  # Should be close to neutral

    # Very concentrated base
    conc = Base("Sodium Hydroxide", concentration=10.0)
    assert conc.ph() > 14  # Can exceed pH 14 for concentrated bases

def test_diprotic_base_properties(diprotic_base):
    """Test properties specific to diprotic bases."""
    assert diprotic_base.proticity == 2
    # For calcium hydroxide: kb1 = 5.5e-2, kb2 = 1.3e-5
    assert diprotic_base.kb1 == 5.5e-2
    assert diprotic_base.kb2 == 1.3e-5
    assert diprotic_base.kb1 > diprotic_base.kb2

# Temperature dependence tests
def test_temperature_effects():
    """Test temperature dependence of equilibrium constants."""
    acid = Acid("Acetic Acid", concentration=0.1)
    ka_25c = acid.ka
    ka_50c = Acid.calculate_temp_dependence(ka_25c, 25, 50, 50000)  # Example ΔH
    assert ka_50c > ka_25c  # Endothermic reaction

# Mass and concentration tests
def test_mass_concentration_relationships():
    """Test relationships between mass, concentration, and volume."""
    acid = Acid("Hydrochloric Acid", concentration=0.1, volume=1.0)
    calc_mass = acid.calculate_mass()
    expected_mass = 0.1 * 1.0 * acid.molar_mass
    assert abs(calc_mass - expected_mass) < 0.001

def test_invalid_chemical_combinations():
    """Test invalid chemical property combinations."""
    with pytest.raises(ValueError):
        Acid("Invalid Acid")
    with pytest.raises(ValueError):
        Base("Invalid Base")
    with pytest.raises(ValueError):
        Acid("Hydrochloric Acid", concentration=0.1, volume=1.0, mass=100)  # Inconsistent mass
