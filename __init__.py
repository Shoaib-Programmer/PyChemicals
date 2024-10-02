# __init__.py

from .predefined_chemicals import valid_acids, valid_bases
from .chemicals import Chemical, Acid, Base

# Example entry point
if __name__ == "__main__":
    acetic_acid = Acid(name="Acetic Acid", volume=0.25, concentration=0.1)
    print(acetic_acid)
