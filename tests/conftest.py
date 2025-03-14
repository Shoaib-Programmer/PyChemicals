"""Pytest fixtures for chemical testing."""

import pytest
from PyChemicals.chemicals import Acid, Base

from PyChemicals.chemicals_db import get_valid_acids, get_valid_bases

valid_acids = get_valid_acids()
valid_bases = get_valid_bases()

hcl_properties = valid_acids["Hydrogen chloride"]
naoh_properties = valid_bases["Sodium hydroxide"]
ch3cooh_properties = valid_acids["Ethanoic acid"]
nh3_properties = valid_bases["Ammonia"]
h2so4_properties = valid_acids["Sulfuric acid"]
caoh2_properties = valid_bases["Calcium hydroxide"]


@pytest.fixture
def strong_acid():
    """Fixture providing a strong acid (HCl) for testing."""
    return Acid(
        "Hydrogen chloride", properties=hcl_properties, concentration=0.1, volume=0.025
    )


@pytest.fixture
def strong_base():
    """Fixture providing a strong base (NaOH) for testing."""
    return Base(
        "Sodium hydroxide", properties=naoh_properties, concentration=0.1, volume=0.025
    )


@pytest.fixture
def weak_acid():
    """Fixture providing a weak acid (CH3COOH) for testing."""
    return Acid(
        "Ethanoic acid", properties=ch3cooh_properties, concentration=0.1, volume=0.025
    )


@pytest.fixture
def weak_base():
    """Fixture providing a weak base (NH3) for testing."""
    return Base("Ammonia", properties=nh3_properties, concentration=0.1, volume=0.025)


@pytest.fixture
def diprotic_acid():
    """Fixture providing a diprotic acid (H2SO4) for testing."""
    return Acid(
        "Sulfuric acid", properties=h2so4_properties, concentration=0.1, volume=0.025
    )


@pytest.fixture
def diprotic_base():
    """Fixture providing a diprotic base (Ca(OH)2) for testing."""
    return Base(
        "Calcium hydroxide",
        properties=caoh2_properties,
        concentration=0.1,
        volume=0.025,
    )
