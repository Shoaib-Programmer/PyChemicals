import numpy as np
from .predefined_chemicals import valid_acids, valid_bases

class Chemical:
    def __init__(self, name: str, concentration: float = None, volume: float = None, mass : float = None):
        self.name = name
        self.concentration = concentration
        self.volume = volume
        self.mass = mass
        self.validate()

    def validate(self):
        if self.concentration is not None and self.concentration <= 0:
            raise ValueError(f"{self.name}: Concentration must be positive.")
        if self.volume is not None and self.volume <= 0:
            raise ValueError(f"{self.name}: Volume must be positive.")
        if self.mass is not None and self.mass <= 0:
            raise ValueError(f"{self.name}: Mass must be positive.")

    def __repr__(self):
        return f"{self.name}: {self.concentration} M, {self.volume} L"
    
    def calculate_moles(self):
        ... # Add the method to calculate the moles when you have the database for the Mr
        if self.concentration is not None and self.volume is not None:
            return self.concentration * self.mass

class Acid(Chemical):
    def __init__(self, name: str, concentration: float = None, volume: float = None, mass: float = None):
        if name not in valid_acids:
            raise ValueError(f"Unknown acid: {name}")
        super().__init__(name, concentration, volume, mass)
        self.proticity = valid_acids[name]["proticity"]
        self.ka = valid_acids[name]["Ka"]
        self.Mr = valid_acids[name]["molar_mass"]
        self.validate()
        self.validate_acid()

    def mass(self):
        if self.mass is None:
            return self.concentration * self.volume * self.Mr


    def H_plus(self):

        if self.ka > 1:
            return self.concentration

        else:
            if self.concentration is None:
                return None # You need the concentration to calculate H_plus
                
            C = self.concentration
            Ka = self.ka

            # Solving the quadratic equation for [H+]
            discriminant = Ka**2 + 4 * Ka * C
            if discriminant < 0:
                raise ValueError("Discriminant is negative; check the values of Ka and concentration.")
            
            # Positive root
            H_plus = (-Ka + np.sqrt(discriminant)) / 2
            return H_plus
        

    def calculate_temp_dependence(K_initial, T_initial, T_final, delta_H):
        """
        Calculate the new equilibrium constant (Ka or Kb) at a different temperature using van 't Hoff equation.
        
        Parameters:
        K_initial (float): The initial equilibrium constant (Ka or Kb) at T_initial.
        T_initial (float): The initial temperature in Celsius.
        T_final (float): The final temperature in Celsius.
        delta_H (float): The enthalpy change (ΔH) of dissociation in J/mol.

        Returns:
        float: The equilibrium constant at T_final.
        """
        # Convert Celsius to Kelvin
        T_initial_K = T_initial + 273.15
        T_final_K = T_final + 273.15
        
        # Universal gas constant R = 8.314 J/(mol·K)
        R = 8.314

        # van 't Hoff equation to calculate new Ka or Kb
        ln_K_ratio = (-delta_H / R) * ((1 / T_final_K) - (1 / T_initial_K))
        K_final = K_initial * np.exp(ln_K_ratio)
        
        return K_final



    def validate_acid(self):
        if self.concentration is not None and self.volume is not None and self.mass is not None:
            if not self.mass / self.Mr == self.concentration * self.volume:
                raise ValueError(f"{self.volume}L of {self.concentration}M {self.name} with {self.mass}g cannot physically exist")
            else:
                pass

    def verify_acid_type(self):
        if self.proticity == 1:
            return "Monoprotic"
        elif self.proticity == 2:
            return "Diprotic"
        else:
            return "Polyprotic"

    def calculate_pKa(self):
        """Calculates the pKa from the dissociation constant Ka."""
        return -np.log10(self.ka)

    def calculate_pH(self):
        """Calculates pH for strong and weak acids."""
        return -np.log10(self.H_plus())
            
    def calculate_moles(self):
        if self.concentration is not None and self.volume is not None:
            return super().calculate_moles()
        elif self.mass is not None:
            return self.mass / self.Mr
        else:
            return None

    def __repr__(self):
        acid_type = self.verify_acid_type()
        return f"{self.name} ({acid_type}): {self.concentration} M, pKa = {self.calculate_pKa():.2f}"

class Base(Chemical):
    def __init__(self, name: str, concentration : float = None, volume : float = None, mass: float = None):
        if name not in valid_bases:
            raise ValueError(f"Unknown base: {name}")
        self.proticity = valid_bases[name]["proticity"]
        self.kb = valid_bases[name]["Kb"]
        super().__init__(name, concentration=concentration, volume=volume, mass=mass)
        self.validate()
        self.Mr = valid_bases[name]["molar_mass"]
        self.validate_base()

    def validate_base(self):
        if self.concentration is not None and self.volume is not None and self.mass is not None:
            if not self.mass / self.Mr == self.concentration * self.volume:
                raise ValueError(f"{self.volume}L of {self.concentration}M {self.name} with {self.mass}g cannot physically exist")
            else:
                pass

    def verify_base_type(self):
        if self.proticity == 1:
            return "Monoprotic"
        elif self.proticity == 2:
            return "Diprotic"
        else:
            return "Polyprotic"
        
    def mass(self):
        if self.mass is None:
            return self.concentration * self.volume * self.Mr
        
    def OH_minus(self):
        if self.kb > 1:
            return self.concentration  # Strong base, complete dissociation

        else:
            if self.concentration is None:
                return None  # Need concentration to calculate OH_minus

            C = self.concentration
            Kb = self.kb

            # Solving the quadratic equation for [OH-]
            discriminant = Kb**2 + 4 * Kb * C
            if discriminant < 0:
                raise ValueError("Discriminant is negative; check the values of Kb and concentration.")
            
            # Positive root
            OH_minus = (-Kb + np.sqrt(discriminant)) / 2
            return OH_minus
        
    def calculate_temp_dependence(K_initial, T_initial, T_final, delta_H):
        """
        Calculate the new equilibrium constant (Ka or Kb) at a different temperature using van 't Hoff equation.
        
        Parameters:
        K_initial (float): The initial equilibrium constant (Ka or Kb) at T_initial.
        T_initial (float): The initial temperature in Celsius.
        T_final (float): The final temperature in Celsius.
        delta_H (float): The enthalpy change (ΔH) of dissociation in J/mol.

        Returns:
        float: The equilibrium constant at T_final.
        """
        # Convert Celsius to Kelvin
        T_initial_K = T_initial + 273.15
        T_final_K = T_final + 273.15
        
        # Universal gas constant R = 8.314 J/(mol·K)
        R = 8.314

        # van 't Hoff equation to calculate new Ka or Kb
        ln_K_ratio = (-delta_H / R) * ((1 / T_final_K) - (1 / T_initial_K))
        K_final = K_initial * np.exp(ln_K_ratio)
        
        return K_final


    def calculate_pKb(self):
        """Calculates the pKb from the dissociation constant Kb."""
        return -np.log10(self.kb)

    def calculate_pOH(self, concentration: float = None, volume: float = None):
        """Calculates pOH for strong and weak bases."""
        self.concentration = concentration
        if self.kb >= 1e-10:  # For strong bases (complete dissociation)
            oh_concentration = self.concentration * self.proticity  # 1 mol NaOH gives 1 mol OH-
            return -np.log10(oh_concentration)
        else:  # For weak bases, using the approximation [OH-] = sqrt(Kb * [B])
            if self.concentration:
                oh_concentration = np.sqrt(self.kb * self.concentration)
                return -np.log10(oh_concentration)
            else:
                raise ValueError("Concentration must be provided to calculate pOH.")

    def calculate_pH(self):
        """Calculates pH from pOH using the relationship pH + pOH = 14."""
        return 14 - self.calculate_pOH()
    
    def calculate_moles(self):
        if self.concentration is not None and self.volume is not None:
            return super().calculate_moles()
        elif self.mass is not None:
            return self.mass / self.Mr
        else:
            return None
        

    def __repr__(self):
        base_type = self.verify_base_type()
        return f"{self.name} ({base_type}): {self.concentration} M, pKb = {self.calculate_pKb():.2f}"

if __name__ == "__main__":
    ... # Implement some sample runs
    # Create an instance of an acid with known concentration and Ka
    hcl = Acid(name="Hydrochloric Acid", concentration=0.1)  # 0.1 M HCL    


