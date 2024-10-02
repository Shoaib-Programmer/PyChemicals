import numpy as np
from .predefined_chemicals import valid_acids, valid_bases

class Chemical:
    def __init__(self, name: str, concentration: float = None, volume: float = None, mass: float = None):
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
        if self.concentration is not None and self.volume is not None:
            return self.concentration * self.volume
        return None

class Acid(Chemical):
    def __init__(self, name: str, concentration: float = None, volume: float = None, mass: float = None):
        if name not in valid_acids:
            raise ValueError(f"Unknown acid: {name}")
        super().__init__(name, concentration, volume, mass)
        self.proticity = valid_acids[name]["proticity"]
        self.ka = valid_acids[name]["Ka"]
        self.Mr = valid_acids[name]["molar_mass"]
        self.validate_acid()

    def calculate_mass(self):
        if self.mass is None:
            return self.concentration * self.volume * self.Mr

    def H_plus(self):
        if self.ka > 1:
            return self.concentration
        else:
            if self.concentration is None:
                return None  # Need the concentration to calculate H_plus

            C = self.concentration
            Ka = self.ka

            # Solving the quadratic equation for [H+]
            discriminant = Ka**2 + 4 * Ka * C
            if discriminant < 0:
                raise ValueError("Discriminant is negative; check the values of Ka and concentration.")

            # Positive root
            H_plus = (-Ka + np.sqrt(discriminant)) / 2
            return H_plus

    @staticmethod
    def calculate_temp_dependence(K_initial, T_initial, T_final, delta_H):
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

    def verify_acid_type(self):
        if self.proticity == 1:
            return "Monoprotic"
        elif self.proticity == 2:
            return "Diprotic"
        else:
            return "Polyprotic"

    def pKa(self):
        """Calculates the pKa from the dissociation constant Ka."""
        return -np.log10(self.ka)

    def pH(self):
        """Calculates pH for strong and weak acids."""
        return -np.log10(self.H_plus())

    def moles(self):
        if self.concentration is not None and self.volume is not None:
            return super().calculate_moles()
        elif self.mass is not None:
            return self.mass / self.Mr
        return None

    def __repr__(self):
        acid_type = self.verify_acid_type()
        return f"{self.name} ({acid_type}): {self.concentration} M, pKa = {self.pKa():.2f}"

class Base(Chemical):
    def __init__(self, name: str, concentration: float = None, volume: float = None, mass: float = None):
        if name not in valid_bases:
            raise ValueError(f"Unknown base: {name}")
        self.proticity = valid_bases[name]["proticity"]
        self.kb = valid_bases[name]["Kb"]
        super().__init__(name, concentration=concentration, volume=volume, mass=mass)
        self.Mr = valid_bases[name]["molar_mass"]
        self.validate_base()

    def validate_base(self):
        if self.concentration is not None and self.volume is not None and self.mass is not None:
            if not self.mass / self.Mr == self.concentration * self.volume:
                raise ValueError(f"{self.volume}L of {self.concentration}M {self.name} with {self.mass}g cannot physically exist")

    def verify_base_type(self):
        if self.proticity == 1:
            return "Monoprotic"
        elif self.proticity == 2:
            return "Diprotic"
        else:
            return "Polyprotic"

    def calculate_mass(self):
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

    @staticmethod
    def calculate_temp_dependence(K_initial, T_initial, T_final, delta_H):
        # Convert Celsius to Kelvin
        T_initial_K = T_initial + 273.15
        T_final_K = T_final + 273.15

        # Universal gas constant R = 8.314 J/(mol·K)
        R = 8.314

        # van 't Hoff equation to calculate new Ka or Kb
        ln_K_ratio = (-delta_H / R) * ((1 / T_final_K) - (1 / T_initial_K))
        K_final = K_initial * np.exp(ln_K_ratio)

        return K_final

    def pKb(self):
        """Calculates the pKb from the dissociation constant Kb."""
        return -np.log10(self.kb)

    def pOH(self):
        """Calculates pOH for strong and weak bases."""
        if self.kb >= 1e-10:  # For strong bases (complete dissociation)
            oh_concentration = self.concentration * float(self.proticity)  # 1 mol NaOH gives 1 mol OH-
            return -np.log10(oh_concentration)
        else:  # For weak bases, using the approximation [OH-] = sqrt(Kb * [B])
            if self.concentration:
                oh_concentration = np.sqrt(self.kb * self.concentration)
                return -np.log10(oh_concentration)
            else:
                raise ValueError("Concentration must be provided to calculate pOH.")

    def pH(self):
        """Calculates pH from pOH using the relationship pH + pOH = 14."""
        return 14 - self.pOH()

    def moles(self):
        if self.concentration is not None and self.volume is not None:
            return super().calculate_moles()
        elif self.mass is not None:
            return self.mass / self.Mr
        return None

    def __repr__(self):
        base_type = self.verify_base_type()
        return f"{self.name} ({base_type}): {self.concentration} M, pKb = {self.pKb():.2f}"
    

class Gas(Chemical):
    R = 0.0821  # Ideal gas constant in L·atm/(K·mol)

    def __init__(self, name: str, pressure: float = None, temperature: float = None, volume: float = None, mass: float = None):
        super().__init__(name, volume=volume, mass=mass)  # Call to the parent constructor
        self.pressure = pressure  # Pressure in atm
        self.temperature = temperature  # Temperature in Kelvin

    def volume_from_moles(self, moles: float) -> float:
        """Calculate volume using the ideal gas law: V = nRT/P."""
        if self.pressure is not None and self.temperature is not None:
            return (moles * self.R * self.temperature) / self.pressure  # V = nRT/P
        raise ValueError("Pressure and temperature must be defined to calculate volume.")

    def moles_from_volume(self, volume: float) -> float:
        """Calculate moles using the ideal gas law: n = PV/RT."""
        if self.pressure is not None and self.temperature is not None:
            return (self.pressure * volume) / (self.R * self.temperature)  # n = PV/RT
        raise ValueError("Pressure and temperature must be defined to calculate moles.")
    
    def density(self) -> float:
        """Calculate the density of the gas."""
        if self.volume == 0:
            raise ValueError("Volume must be greater than zero to calculate density.")
        return self.mass / self.volume  # Density = mass / volume

    def molar_mass(self) -> float:
        """Calculate the molar mass of the gas."""
        # Assuming 'mass' is the mass of the gas in grams and 'moles' is the number of moles.
        if self.moles == 0:
            raise ValueError("Moles must be greater than zero to calculate molar mass.")
        return self.mass / self.moles  # Molar mass = mass / moles

    def partial_pressure(self, total_pressure: float, mole_fraction: float) -> float:
        """Calculate the partial pressure of the gas in a mixture."""
        return total_pressure * mole_fraction  # Partial pressure = total pressure * mole fraction
    
    def grahams_law_effusion(self, molar_mass_other: float) -> float:
        """Calculate the rate of effusion relative to another gas using Graham's Law."""
        return (self.molar_mass / molar_mass_other) ** 0.5  # Rate of effusion ∝ 1/√(molar mass)

    def compressibility_factor(self, pressure: float, volume: float, temperature: float) -> float:
        """Calculate the compressibility factor Z of the gas."""
        R = 0.0821  # Ideal gas constant in L·atm/(K·mol)
        # Z = PV/nRT, where n = moles
        return (pressure * volume) / (self.moles * R * temperature)  # Z = PV/nRT
    
    def work_done(self, n: float, T: float, V_i: float, V_f: float) -> float:
        """Calculate the work done by the gas during isothermal expansion."""
        R = 0.0821  # Ideal gas constant in L·atm/(K·mol)
        return -n * R * T * (np.log(V_f / V_i))  # Work = -nRT ln(V_f/V_i)

    def __repr__(self):
        return f"{self.name}: {self.pressure} atm, {self.temperature} K, {self.volume} L"
