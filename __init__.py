# pylint: disable=missing-module-docstring
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring

try:
    from .predefined_chemicals import valid_acids, valid_bases
    from .chemicals import Chemical, Acid, Base
except ImportError:
    from predefined_chemicals import valid_acids, valid_bases
    from chemicals import Chemical, Acid, Base

__all__ = ["Chemical", "Acid", "Base", "valid_acids", "valid_bases"]


def main():
    acetic_acid = Acid(name="Acetic Acid", volume=0.25, concentration=0.1)
    print(acetic_acid)


if __name__ == "__main__":
    main()
