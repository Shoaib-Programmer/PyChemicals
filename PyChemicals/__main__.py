# pylint: disable=missing-module-docstring
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring


from .chemicals import Acid, Base
from .titration_curves import titrate_curve_monoprotic
from .titration_calculation import Analyte, Titrant
from .chemicals_db import get_valid_acids, get_valid_bases

hcl_properties = get_valid_acids()["Hydrogen chloride"]
naoh_properties = get_valid_bases()["Sodium hydroxide"]


def main():
    # Example usage
    print("Testing PyChemicals module...")

    # Create and test an acid
    hcl = Acid(name="Hydrogen chloride", concentration=0.1, properties=hcl_properties)
    print(f"HCl pH: {hcl.ph():.2f}")

    # Create and test a base
    naoh = Base(name="Sodium hydroxide", concentration=0.1, properties=naoh_properties)
    print(f"NaOH pH: {naoh.ph():.2f}")
    # Example: Titrate 25 mL of 0.1 M HCl (acid) with 0.1 M NaOH (base)
    # Note: The titrant's volume is provided (50 mL here) to satisfy validation,
    # though the plotting function uses a range of added volumes.
    # Create instances using the factory behavior in Analyte and Titrant.
    # Note: If "Hydrochloric Acid" and "Sodium Hydroxide" are in valid_acids and valid_bases,
    # they will be returned as Acid and Base instances respectively.
    hcl = Analyte("Hydrochloric Acid", concentration=0.1, volume=0.025)
    naoh = Titrant("Sodium Hydroxide", concentration=0.1, volume=0.025)

    # Plot the titration curve
    titrate_curve_monoprotic(hcl, naoh)


if __name__ == "__main__":
    main()
