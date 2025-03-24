"""
PyChemicals

This library provides tools, classes, and methods for handling
common chemical entities and operations.

It is designed to aid in conducting calculations and analyses for substances
such as chemicals, acids, bases, and gases.

You can find more information about the project at:
    https://github.com/shoaib-quantumcalc/PyChemicals

You can find the documentation for this library at:
    https://shoaib-quantumcalc.github.io/PyChemicals/

Author: Shoaib Nigam Shaik <shoaibnigam422@gmail.com>
License: BSD-3-Clause
Version: 0.0.1
    Initial release.
"""
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring

try:
    from .chemicals import Chemical, Acid, Base
except ImportError:
    print("ImportError: Unable to import chemicals module.")
    from chemicals import Chemical, Acid, Base

# Import all contents from titration_calculation and titration_curves
from . import titration_calculation
from .titration_curves import titrate_curve_monoprotic, volume_of_titrant  # noqa: F403

# Define __all__ to explicitly state what gets exported
__all__ = [
    "Chemical",
    "Acid",
    "Base",
    "titration_calculation",
    "titrate_curve_monoprotic",
    "volume_of_titrant",
]
