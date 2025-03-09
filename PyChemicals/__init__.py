# pylint: disable=missing-module-docstring
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring

try:
    from .chemicals import Chemical, Acid, Base, valid_acids, valid_bases, valid_gases
except ImportError:
    from chemicals import Chemical, Acid, Base, valid_acids, valid_bases, valid_gases

__all__ = ["Chemical", "Acid", "Base", "valid_acids", "valid_bases", "valid_gases"]


def main():
    acetic_acid = Acid(name="Acetic Acid", volume=0.25, concentration=0.1)
    print(acetic_acid)


if __name__ == "__main__":
    main()
