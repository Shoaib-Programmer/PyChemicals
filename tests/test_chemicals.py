"""Test module for chemical classes including Chemical, Acid, and Base."""

import pytest
import numpy as np
from PyChemicals.chemicals import Chemical, Acid, Base

def test_chemical_initialization():
    """Test basic initialization of Chemical class."""
    chem = Chemical(name="Test Chemical", concentration=0.1, volume=0.5)
    assert chem.name == "Test Chemical"
    assert chem.concentration == 0.1
    assert chem.volume == 0.5

def test_chemical_validation():
    """Test validation of chemical properties."""
    with pytest.raises(ValueError):
        Chemical(name="Test", concentration=-1)
    with pytest.raises(ValueError):
        Chemical(name="Test", volume=-1)

def test_acid_initialization():
    """Test initialization of Acid class."""
    hcl = Acid("Hydrochloric Acid", concentration=0.1)
    assert hcl.name == "Hydrochloric Acid"
    assert hcl.concentration == 0.1
    assert hcl.proticity == 1

def test_acid_ph():
    """Test pH calculation for strong acid."""
    hcl = Acid("Hydrochloric Acid", concentration=0.1)
    assert abs(hcl.ph() - 1.0) < 0.1  # pH should be around 1 for 0.1M HCl

def test_base_initialization():
    """Test initialization of Base class."""
    naoh = Base("Sodium Hydroxide", concentration=0.1)
    assert naoh.name == "Sodium Hydroxide"
    assert naoh.concentration == 0.1
    assert naoh.proticity == 1

def test_base_ph():
    """Test pH calculation for strong base."""
    naoh = Base("Sodium Hydroxide", concentration=0.1)
    assert abs(naoh.ph() - 13.0) < 0.1  # pH should be around 13 for 0.1M NaOH

def test_invalid_acid():
    """Test error handling for invalid acid names."""
    with pytest.raises(ValueError):
        Acid("Invalid Acid")

def test_invalid_base():
    """Test error handling for invalid base names."""
    with pytest.raises(ValueError):
        Base("Invalid Base")

def test_acid_calculations():
    """Test various acid calculations including pKa."""
    acetic = Acid("Acetic Acid", concentration=0.1)
    assert abs(acetic.pka() - (-np.log10(1.8e-5))) < 0.01
    assert acetic.verify_acid_type() == "Monoprotic"

def test_base_calculations():
    """Test various base calculations including pKb."""
    ammonia = Base("Ammonia", concentration=0.1)
    assert abs(ammonia.pkb() - (-np.log10(1.76e-5))) < 0.01
    assert ammonia.verify_base_type() == "Monoprotic"
