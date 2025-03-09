"""Test module for predefined chemicals data."""

from PyChemicals.chemicals import valid_acids, valid_bases


def test_valid_acids_data():
    """Test validity and structure of predefined acids data."""
    assert "Hydrogen chloride" in valid_acids
    assert "Ethanoic acid" in valid_acids
    assert valid_acids["Hydrogen chloride"]["proticity"] == 1
    assert valid_acids["Sulfuric acid"]["proticity"] == 2


def test_valid_bases_data():
    """Test validity and structure of predefined bases data."""
    assert "Sodium hydroxide" in valid_bases
    assert "Ammonia" in valid_bases
    assert valid_bases["Sodium hydroxide"]["proticity"] == 1
    assert valid_bases["Calcium hydroxide"]["proticity"] == 2


def test_acid_properties():
    """Test specific properties of predefined acids."""
    hcl = valid_acids["Hydrogen chloride"]
    assert hcl["Ka"] > 1  # Strong acid
    assert hcl["molar_mass"] > 0
    assert "Ka_temp" in hcl


def test_base_properties():
    """Test specific properties of predefined bases."""
    naoh = valid_bases["Sodium hydroxide"]
    assert naoh["Kb"] > 1  # Strong base
    assert naoh["molar_mass"] > 0
    assert "Kb_temp" in naoh
