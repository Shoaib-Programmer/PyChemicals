# pylint: disable=missing-module-docstring
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring

try:
    from .chemicals import Chemical, Acid, Base
except ImportError:
    print("ImportError: Unable to import chemicals module.")
    from chemicals import Chemical, Acid, Base

from . import titration_calculation
from . import titration_curves

__all__ = [
    "Chemical",
    "Acid",
    "Base",
    "titration_calculation",
    "titration_curves",
]


def main():
    acetic_acid = Acid(name="Acetic Acid", volume=0.25, concentration=0.1)
    print(acetic_acid)


if __name__ == "__main__":
    main()
