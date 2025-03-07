"""Pytest fixtures for chemical testing."""

import pytest
from PyChemicals.chemicals import Acid, Base

@pytest.fixture
def strong_acid():
    """Fixture providing a strong acid (HCl) for testing."""
    return Acid("Hydrochloric Acid", concentration=0.1, volume=0.025)

@pytest.fixture
def strong_base():
    """Fixture providing a strong base (NaOH) for testing."""
    return Base("Sodium Hydroxide", concentration=0.1, volume=0.025)

@pytest.fixture
def weak_acid():
    """Fixture providing a weak acid (CH3COOH) for testing."""
    return Acid("Acetic Acid", concentration=0.1, volume=0.025)

@pytest.fixture
def weak_base():
    """Fixture providing a weak base (NH3) for testing."""
    return Base("Ammonia", concentration=0.1, volume=0.025)
