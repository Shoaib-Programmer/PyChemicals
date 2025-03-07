"""Test module for predefined chemicals data."""

from PyChemicals.predefined_chemicals import valid_acids, valid_bases

def test_valid_acids_data():
    """Test validity and structure of predefined acids data."""
    assert "Hydrochloric Acid" in valid_acids
    assert "Acetic Acid" in valid_acids
    assert valid_acids["Hydrochloric Acid"]["proticity"] == 1
    assert valid_acids["Sulfuric Acid"]["proticity"] == 2

def test_valid_bases_data():
    """Test validity and structure of predefined bases data."""
    assert "Sodium Hydroxide" in valid_bases
    assert "Ammonia" in valid_bases
    assert valid_bases["Sodium Hydroxide"]["proticity"] == 1
    assert valid_bases["Calcium Hydroxide"]["proticity"] == 2

def test_acid_properties():
    """Test specific properties of predefined acids."""
    hcl = valid_acids["Hydrochloric Acid"]
    assert hcl["Ka"] > 1  # Strong acid
    assert hcl["molar_mass"] > 0
    assert "Ka_temp" in hcl

def test_base_properties():
    """Test specific properties of predefined bases."""
    naoh = valid_bases["Sodium Hydroxide"]
    assert naoh["Kb"] > 1  # Strong base
    assert naoh["molar_mass"] > 0
    assert "Kb_temp" in naoh
