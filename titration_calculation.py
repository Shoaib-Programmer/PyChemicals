import numpy as np
from chemicals import Acid, Base, Chemical
from predefined_chemicals import valid_acids, valid_bases
from chempy import balance_stoichiometry

class Analyte(Chemical):
    def __init__(self, name: str, volume: float = None, titrant=None):
        if name in valid_acids:
            super().__init__(name=name, concentration=None, volume=volume)
            self.__class__ = Acid  # Change the class type to Acid
        elif name in valid_bases:
            super().__init__(name=name, concentration=None, volume=volume)
            self.__class__ = Base  # Change the class type to Base
        else:
            super().__init__(name=name, concentration=None, volume=volume)

        if titrant is not None:
            self.validate_titrant(titrant)

    def validate_titrant(self, titrant):
        if isinstance(titrant, Acid) and isinstance(self, Acid):
            raise ValueError("Both Analyte and Titrant cannot be acids.")
        elif isinstance(titrant, Base) and isinstance(self, Base):
            raise ValueError("Both Analyte and Titrant cannot be bases.")


class Titrant(Chemical):
    def __init__(self, name: str, concentration: float = None, volume: float = None, analyte=None):
        if name in valid_acids:
            super().__init__(name=name, concentration=concentration, volume=volume)
            self.__class__ = Acid  # Change the class type to Acid
        elif name in valid_bases:
            super().__init__(name=name, concentration=concentration, volume=volume)
            self.__class__ = Base  # Change the class type to Base
        else:
            super().__init__(name=name, concentration=concentration, volume=volume)

        if analyte is not None:
            self.validate_analyte(analyte)

    def validate_analyte(self, analyte):
        if isinstance(analyte, Acid) and isinstance(self, Acid):
            raise ValueError("Both Analyte and Titrant cannot be acids.")
        elif isinstance(analyte, Base) and isinstance(self, Base):
            raise ValueError("Both Analyte and Titrant cannot be bases.")



def volume_of_titrant(analyte: Analyte, titrant: Titrant, stoichiometric_ratio: float = 1.0):
    """Calculate the volume of titrant needed to reach equivalence.
    
    Args:
        analyte (Analyte): The analyte being titrated.
        titrant (Titrant): The titrant used in the titration.
        stoichiometric_ratio (float): The stoichiometric ratio of the analyte to the titrant.

    Returns:
        float: Volume of titrant needed in liters.
    
    Raises:
        ValueError: If concentration or volume of analyte or titrant is not provided.
    """
    # Check that analyte has a valid concentration and volume
    if analyte.concentration is None or analyte.volume is None:
        raise ValueError("Analyte concentration and volume must be provided.")
    
    # Check that titrant has a valid concentration
    if titrant.concentration is None:
        raise ValueError("Titrant concentration must be provided.")

    # Calculate moles of analyte
    moles_analyte = analyte.concentration * analyte.volume

    # Calculate moles of titrant needed based on the stoichiometric ratio
    moles_titrant_needed = moles_analyte / stoichiometric_ratio

    # Calculate the volume of titrant needed (V = n / C)
    volume_titrant_needed = moles_titrant_needed / titrant.concentration

    # Output results
    return volume_titrant_needed

def C_analyte(analyte: Analyte, titrant: Titrant, stoichiometric_ratio: float):
    """
    Calculate the concentration of the analyte using the volumes and concentration of the titrant.

    Parameters:
    - analyte: Analyte object
    - titrant: Titrant object
    - volume_titrant: Volume of the titrant used (in liters)
    - stoichiometric_ratio: Stoichiometric ratio between the analyte and titrant

    Returns:
    - Concentration of the analyte (in mol/L)
    """

    # Validate input parameters
    if analyte.volume is None:
        raise ValueError("Volume of analyte must be provided.")
    if titrant.volume is None:
        raise ValueError("Volume of titrant must be provided.")
    if titrant.concentration is None:
        raise ValueError("Concentration of titrant must be provided.")
    if stoichiometric_ratio <= 0:
        raise ValueError("Stoichiometric ratio must be greater than zero.")

    # Calculate moles of titrant used
    moles_titrant = titrant.concentration * titrant.volume

    # Calculate moles of analyte using stoichiometric ratio
    moles_analyte = moles_titrant * (stoichiometric_ratio)

    # Calculate concentration of analyte (C = n / V)
    concentration_analyte = moles_analyte / analyte.volume

    return concentration_analyte


def pH_at_equivalence(analyte: Analyte, titrant: Titrant):
    # This is a placeholder for the actual pH calculation logic.
    # The pH can depend on whether the titrant is a strong or weak acid/base.
    # You can implement more advanced logic based on the types of analyte and titrant.

    if isinstance(analyte, Acid) and isinstance(titrant, Base):
        # Example calculation for strong acid + strong base reaction
        return 7.0  # Neutralization point
    elif isinstance(analyte, Base) and isinstance(titrant, Acid):
        return 7.0  # Neutralization point
    # Implement other cases based on weak acids/bases
    return None  # Default case, return None if no specific pH calculation is done


def calculate_pH(acid, volume_added):
    """Calculate pH at different volumes of titrant added."""
    # Implement using the Henderson-Hasselbalch equation for weak acids/bases
    if isinstance(acid, Acid):
        # Assuming strong acid vs weak base or similar logic
        return acid.calculate_pH()
    return None

def calculate_half_equivalence(acid):
    """Calculate pH at the half-equivalence point."""
    ...


def pH_after_equivalence(acid : Acid, base: Base, excess_volume):
    """Calculate pH after reaching the equivalence point."""
    if isinstance(acid, Acid) and isinstance(base, Base):
        if base.concentration > acid.concentration:
            # Calculate excess OH⁻ from the base
            excess_oh_concentration = base.concentration - (acid.concentration * acid.volume / base.volume)
            return 14 + np.log10(excess_oh_concentration)  # pOH to pH conversion
    return None

# Example of using the functions
if __name__ == "__main__":
    hcl = Acid("Hydrochloric Acid", concentration=0.1, volume=0.025)  # 0.1 M HCl, 25 mL
    naoh = Base("Sodium Hydroxide", concentration=0.1, volume=0.025)  # 0.1 M NaOH, 25 mL

    print("pH at Equivalence:", pH_at_equivalence(hcl, naoh))
    print("Volume of Titrant Needed:", volume_of_titrant(hcl, naoh))
    print("pH after Equivalence:", pH_after_equivalence(hcl, naoh, 0.01))  # excess volume of 10 mL
