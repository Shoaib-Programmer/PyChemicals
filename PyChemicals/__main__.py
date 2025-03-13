# pylint: disable=missing-module-docstring
# pylint: disable=invalid-name
# pylint: disable=missing-function-docstring


from .chemicals import Acid, Base
from .titration_curves import titrate_curve_monoprotic
from .titration_calculation import Analyte, Titrant


def main():
    # Example usage
    print("Testing PyChemicals module...")

    # Create and test an acid
    hcl = Acid(name="Hydrogen chloride", concentration=0.1)
    print(f"HCl pH: {hcl.ph():.2f}")

    # Create and test a base
    naoh = Base(name="Sodium hydroxide", concentration=0.1)
    print(f"NaOH pH: {naoh.ph():.2f}")
    # Example: Titrate 25 mL of 0.1 M HCl (acid) with 0.1 M NaOH (base)
    # Note: The titrant's volume is provided (50 mL here) to satisfy validation,
    # though the plotting function uses a range of added volumes.
    hcl = Analyte("Hydrogen chloride", volume=0.025, concentration=0.1)
    naoh = Titrant("Sodium hydroxide", volume=0.05, concentration=0.1)
    
    # Plot the titration curve
    titrate_curve_monoprotic(hcl, naoh)


if __name__ == "__main__":
    main()
