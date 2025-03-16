# pylint: disable=missing-module-docstring
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring

try:
    from .chemicals import Chemical, Acid, Base
except ImportError:
    print("ImportError: Unable to import chemicals module.")
    from chemicals import Chemical, Acid, Base

# Import all contents from titration_calculation and titration_curves
from .titration_calculation import *  # noqa: F403
from .titration_curves import *  # noqa: F403

# Define __all__ to explicitly state what gets exported
__all__ = [
    "Chemical",
    "Acid",
    "Base",
]

# Get all public attributes from titration_calculation and titration_curves
from . import titration_calculation, titration_curves

__all__ += [name for name in dir(titration_calculation) if not name.startswith("_")]
__all__ += [name for name in dir(titration_curves) if not name.startswith("_")]
