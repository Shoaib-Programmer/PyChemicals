"""
This module contains all the classes for the chemicals in the PyChemicals package.
"""

import numpy as np
from cs50 import SQL  # type: ignore

# Connect to the chemicals database (this file should be in your current working directory)
db = SQL("sqlite:///chemicals.db")


def get_acids():
    """Retrieve acids data from the database and return a dictionary."""
    rows = db.execute("SELECT * FROM acids")
    acids = {}
    for row in rows:
        acid_id = row["id"]
        chem_name = row["name"]
        temp_rows = db.execute(
            "SELECT temperature, ka FROM acid_temperatures WHERE acid_id = ?", acid_id
        )
        ka_temp = {t["temperature"]: t["ka"] for t in temp_rows}
        acids[chem_name] = {
            "proticity": row["proticity"],
            "Ka": row["Ka"],
            "molar_mass": row["molar_mass"],
            "ka1": row["ka1"],
            "ka2": row["ka2"],
            "Ka_temp": ka_temp,
        }
    return acids


def get_bases():
    """Retrieve bases data from the database and return a dictionary."""
    rows = db.execute("SELECT * FROM bases")
    bases = {}
    for row in rows:
        base_id = row["id"]
        chem_name = row["name"]
        temp_rows = db.execute(
            "SELECT temperature, kb FROM base_temperatures WHERE base_id = ?", base_id
        )
        kb_temp = {t["temperature"]: t["kb"] for t in temp_rows}
        bases[chem_name] = {
            "proticity": row["proticity"],
            "Kb": row["Kb"],
            "molar_mass": row["molar_mass"],
            "kb1": row["kb1"],
            "kb2": row["kb2"],
            "Kb_temp": kb_temp,
        }
    return bases


def get_gases():
    """Retrieve gases data from the database and return a dictionary."""
    rows = db.execute("SELECT * FROM gases")
    gases = {}
    for row in rows:
        chem_name = row["name"]
        gases[chem_name] = {
            "molar_mass": row["molar_mass"],
            "density": row[
                "density"
            ],  # Note: This is a reference density (typically at standard conditions)
        }
    return gases


valid_acids = get_acids()
valid_bases = get_bases()
valid_gases = get_gases()


class Chemical:
    """Base class for all chemicals."""

    def __init__(
        self,
        name: str,
        concentration: float = None,
        volume: float = None,
        mass: float = None,
    ):
        """
        Initialize a Chemical instance.

        Args:
            name (str): Name of the chemical.
            concentration (float, optional): Concentration in molarity.
            volume (float, optional): Volume in liters.
            mass (float, optional): Mass in grams.
        """
        self.name = name
        self.concentration = concentration
        self.volume = volume
        self.mass = mass
        self.validate()

    def validate(self):
        """Validate chemical properties."""
        if self.concentration is not None and self.concentration <= 0:
            raise ValueError(f"{self.name}: Concentration must be positive.")
        if self.volume is not None and self.volume <= 0:
            raise ValueError(f"{self.name}: Volume must be positive.")
        if self.mass is not None and self.mass <= 0:
            raise ValueError(f"{self.name}: Mass must be positive.")

    def __repr__(self):
        return f"{self.name}: {self.concentration} M, {self.volume} L"

    def calculate_moles(self):
        """
        Calculate moles from concentration and volume.

        Returns:
            float or None: Calculated moles or None if data is missing.
        """
        if self.concentration is not None and self.volume is not None:
            return self.concentration * self.volume
        return None


class Acid(Chemical):
    """Class representing an acid chemical."""

    def __init__(
        self,
        name: str,
        concentration: float = None,
        volume: float = None,
        mass: float = None,
    ):
        """
        Initialize an Acid instance.

        Args:
            name (str): Name of the acid.
            concentration (float, optional): Concentration in molarity.
            volume (float, optional): Volume in liters.
            mass (float, optional): Mass in grams.
        """
        if name not in valid_acids:
            raise ValueError(f"Unknown acid: {name}")
        super().__init__(name, concentration, volume, mass)
        self.proticity = valid_acids[name]["proticity"]
        self.ka = valid_acids[name]["Ka"]
        self.ka1 = valid_acids[name]["ka1"]  # For diprotic acids
        self.ka2 = valid_acids[name]["ka2"]
        self.molar_mass = valid_acids[name]["molar_mass"]
        self.validate_acid()

    def calculate_mass(self):
        """
        Calculate the mass if not provided.

        Returns:
            float: Calculated mass.
        """
        if self.mass is None:
            return self.concentration * self.volume * self.molar_mass
        return self.mass

    def h_plus(self):
        """
        Calculate the concentration of H+ ions.

        Returns:
            float or None: H+ concentration.
        """
        if self.ka > 1:
            return self.concentration
        if self.concentration is None:
            return None  # Need the concentration to calculate H+
        conc = self.concentration
        ka_value = self.ka
        # Coefficients for the quadratic equation:
        # [H+]^2 + ka_value * [H+] - ka_value * conc = 0
        a = 1
        b = ka_value
        c = -ka_value * conc
        discriminant = b**2 - 4 * a * c
        if discriminant < 0:
            raise ValueError(
                "Discriminant is negative; check the values of Ka and concentration."
            )
        h_plus_conc = (-b + np.sqrt(discriminant)) / (2 * a)
        return h_plus_conc

    @staticmethod
    def calculate_temp_dependence(k_initial, t_initial, t_final, delta_h):
        """
        Calculate temperature dependence of the equilibrium constant using the van 't Hoff equation.

        Args:
            k_initial (float): Initial equilibrium constant.
            t_initial (float): Initial temperature in Celsius.
            t_final (float): Final temperature in Celsius.
            delta_h (float): Enthalpy change in J/mol.

        Returns:
            float: New equilibrium constant.
        """
        t_initial_k = t_initial + 273.15
        t_final_k = t_final + 273.15
        gas_constant = 8.314  # J/(mol·K)
        ln_k_ratio = (-delta_h / gas_constant) * ((1 / t_final_k) - (1 / t_initial_k))
        k_final = k_initial * np.exp(ln_k_ratio)
        return k_final

    def validate_acid(self):
        """Validate the acid properties."""
        if (
            self.concentration is not None
            and self.volume is not None
            and self.mass is not None
        ):
            if not self.mass / self.molar_mass == self.concentration * self.volume:
                raise ValueError(
                    f"{self.volume}L of {self.concentration}M {self.name} with {self.mass}g "
                    "cannot physically exist"
                )

    def verify_acid_type(self):
        """
        Verify the type of acid based on proticity.

        Returns:
            str: Acid type as 'Monoprotic', 'Diprotic', or 'Polyprotic'.
        """
        if self.proticity == 1:
            return "Monoprotic"
        if self.proticity == 2:
            return "Diprotic"
        return "Polyprotic"

    def pka(self):
        """
        Calculate the pKa value.

        Returns:
            float: pKa value.
        """
        return -np.log10(self.ka)

    def ph(self):
        """
        Calculate the pH of the acid.

        Returns:
            float: pH value.
        """
        return -np.log10(self.h_plus())

    def moles(self):
        """
        Calculate the number of moles of the acid.

        Returns:
            float or None: Number of moles.
        """
        if self.concentration is not None and self.volume is not None:
            return super().calculate_moles()
        if self.mass is not None:
            return self.mass / self.molar_mass
        return None

    def __repr__(self):
        acid_type = self.verify_acid_type()
        return (
            f"{self.name} ({acid_type}): {self.concentration} M, pKa = {self.pka():.2f}"
        )


class Base(Chemical):
    """Class representing a base chemical."""

    def __init__(
        self,
        name: str,
        concentration: float = None,
        volume: float = None,
        mass: float = None,
    ):
        """
        Initialize a Base instance.

        Args:
            name (str): Name of the base.
            concentration (float, optional): Concentration in molarity.
            volume (float, optional): Volume in liters.
            mass (float, optional): Mass in grams.
        """
        if name not in valid_bases:
            raise ValueError(f"Unknown base: {name}")
        self.proticity = valid_bases[name]["proticity"]
        self.kb = valid_bases[name]["Kb"]
        self.kb1 = valid_bases[name]["kb1"]  # For diprotic bases
        self.kb2 = valid_bases[name]["kb2"]
        super().__init__(name, concentration=concentration, volume=volume, mass=mass)
        self.molar_mass = valid_bases[name]["molar_mass"]
        self.validate_base()

    def validate_base(self):
        """Validate the base properties."""
        if (
            self.concentration is not None
            and self.volume is not None
            and self.mass is not None
        ):
            if not self.mass / self.molar_mass == self.concentration * self.volume:
                raise ValueError(
                    f"{self.volume}L of {self.concentration}M {self.name} with {self.mass}g "
                    "cannot physically exist"
                )

    def verify_base_type(self):
        """
        Verify the type of base based on proticity.

        Returns:
            str: Base type as 'Monoprotic', 'Diprotic', or 'Polyprotic'.
        """
        if self.proticity == 1:
            return "Monoprotic"
        if self.proticity == 2:
            return "Diprotic"
        return "Polyprotic"

    def calculate_mass(self):
        """
        Calculate the mass if not provided.

        Returns:
            float: Calculated mass.
        """
        if self.mass is None:
            return self.concentration * self.volume * self.molar_mass
        return self.mass

    def oh_minus(self):
        """
        Calculate the concentration of OH- ions.

        Returns:
            float or None: OH- concentration.
        """
        if self.kb > 1:
            return self.concentration  # Strong base, complete dissociation
        if self.concentration is None:
            return None  # Need concentration to calculate OH-
        conc = self.concentration
        kb_value = self.kb
        discriminant = kb_value**2 + 4 * kb_value * conc
        if discriminant < 0:
            raise ValueError(
                "Discriminant is negative; check the values of Kb and concentration."
            )
        oh_minus_conc = (-kb_value + np.sqrt(discriminant)) / 2
        return oh_minus_conc

    @staticmethod
    def calculate_temp_dependence(k_initial, t_initial, t_final, delta_h):
        """
        Calculate temperature dependence of the equilibrium constant using the van 't Hoff equation.

        Args:
            k_initial (float): Initial equilibrium constant.
            t_initial (float): Initial temperature in Celsius.
            t_final (float): Final temperature in Celsius.
            delta_h (float): Enthalpy change in J/mol.

        Returns:
            float: New equilibrium constant.
        """
        t_initial_k = t_initial + 273.15
        t_final_k = t_final + 273.15
        gas_constant = 8.314  # J/(mol·K)
        ln_k_ratio = (-delta_h / gas_constant) * ((1 / t_final_k) - (1 / t_initial_k))
        k_final = k_initial * np.exp(ln_k_ratio)
        return k_final

    def pkb(self):
        """
        Calculate the pKb value.

        Returns:
            float: pKb value.
        """
        return -np.log10(self.kb)

    def poh(self):
        """
        Calculate the pOH of the base.

        Returns:
            float: pOH value.
        """
        if self.kb >= 1e-10:  # For strong bases (complete dissociation)
            oh_conc = self.concentration * float(self.proticity)
            return -np.log10(oh_conc)
        if self.concentration:
            oh_conc = np.sqrt(self.kb * self.concentration)
            return -np.log10(oh_conc)
        raise ValueError("Concentration must be provided to calculate pOH.")

    def ph(self):
        """
        Calculate the pH of the base.

        Returns:
            float: pH value.
        """
        return 14 - self.poh()

    def moles(self):
        """
        Calculate the number of moles of the base.

        Returns:
            float or None: Number of moles.
        """
        if self.concentration is not None and self.volume is not None:
            return super().calculate_moles()
        if self.mass is not None:
            return self.mass / self.molar_mass
        return None

    def __repr__(self):
        base_type = self.verify_base_type()
        return (
            f"{self.name} ({base_type}): {self.concentration} M, pKb = {self.pkb():.2f}"
        )


class Gas(Chemical):
    """Class representing a gas."""

    R = 0.0821  # Ideal gas constant in L·atm/(K·mol)

    def __init__(  # pylint: disable=too-many-arguments, too-many-positional-arguments
        self,
        name: str,
        pressure: float = None,
        temperature: float = None,
        volume: float = None,
        mass: float = None,
    ):
        """
        Initialize a Gas instance.

        Args:
            name (str): Name of the gas.
            pressure (float, optional): Pressure in atm.
            temperature (float, optional): Temperature in Kelvin.
            volume (float, optional): Volume in liters.
            mass (float, optional): Mass in grams.
        """
        if name not in valid_gases:
            raise ValueError(f"Invalid gas: {name}")
        super().__init__(name, volume=volume, mass=mass)
        self.pressure = pressure
        self.temperature = temperature

    @property
    def moles(self):
        """
        Calculate the number of moles of the gas using the ideal gas law.

        Returns:
            float: Number of moles.
        """
        if (
            self.pressure is not None
            and self.temperature is not None
            and self.volume is not None
        ):
            return (self.pressure * self.volume) / (self.R * self.temperature)
        raise ValueError("Insufficient data to calculate moles for gas.")

    def volume_from_moles(self, moles: float) -> float:
        """
        Calculate volume using the ideal gas law: V = nRT/P.

        Args:
            moles (float): Number of moles.

        Returns:
            float: Volume in liters.
        """
        if self.pressure is not None and self.temperature is not None:
            return (moles * self.R * self.temperature) / self.pressure
        raise ValueError(
            "Pressure and temperature must be defined to calculate volume."
        )

    def moles_from_volume(self, volume: float) -> float:
        """
        Calculate moles using the ideal gas law: n = PV/RT.

        Args:
            volume (float): Volume in liters.

        Returns:
            float: Number of moles.
        """
        if self.pressure is not None and self.temperature is not None:
            return (self.pressure * volume) / (self.R * self.temperature)
        raise ValueError("Pressure and temperature must be defined to calculate moles.")

    def density(self) -> float:
        """
        Calculate the density of the gas.

        Returns:
            float: Density in g/L.
        """
        if self.volume == 0:
            raise ValueError("Volume must be greater than zero to calculate density.")
        return self.mass / self.volume

    def molar_mass_calc(self) -> float:
        """
        Calculate the molar mass of the gas.

        Returns:
            float: Molar mass in g/mol.
        """
        n = self.moles
        if n == 0:
            raise ValueError("Moles must be greater than zero to calculate molar mass.")
        return self.mass / n

    def partial_pressure(self, total_pressure: float, mole_fraction: float) -> float:
        """
        Calculate the partial pressure of the gas in a mixture.

        Args:
            total_pressure (float): Total pressure in atm.
            mole_fraction (float): Mole fraction of the gas.

        Returns:
            float: Partial pressure in atm.
        """
        return total_pressure * mole_fraction

    def grahams_law_effusion(self, molar_mass_other: float) -> float:
        """
        Calculate the rate of effusion relative to another gas using Graham's Law.

        Args:
            molar_mass_other (float): Molar mass of the other gas in g/mol.

        Returns:
            float: Relative rate of effusion.
        """
        return (self.molar_mass_calc() / molar_mass_other) ** 0.5

    def compressibility_factor(
        self, pressure: float, volume: float, temperature: float
    ) -> float:
        """
        Calculate the compressibility factor Z of the gas: Z = PV/nRT.

        Args:
            pressure (float): Pressure in atm.
            volume (float): Volume in liters.
            temperature (float): Temperature in Kelvin.

        Returns:
            float: Compressibility factor.
        """
        return (pressure * volume) / (self.moles * self.R * temperature)

    def work_done(self, v_initial: float, v_final: float) -> float:
        """
        Calculate the work done by the gas during isothermal expansion.

        Args:
            v_initial (float): Initial volume in liters.
            v_final (float): Final volume in liters.

        Returns:
            float: Work done in L·atm.
        """
        return -self.moles * self.R * self.temperature * np.log(v_final / v_initial)

    def __repr__(self):
        return (
            f"{self.name}: {self.pressure} atm, {self.temperature} K, {self.volume} L"
        )


if __name__ == "__main__":
    # Implement test runs
    hcl = Acid(name="Hydrochloric Acid", concentration=0.3)
    print(hcl.ph())

    print("Acids:")
    for acid_name, data in valid_acids.items():  # Changed name to acid_name
        print(f"  {acid_name}: {data}")

    print("\nBases:")
    for base_name, data in valid_bases.items():  # Changed name to base_name
        print(f"  {base_name}: {data}")

    print("\nGases:")
    for gas_name, data in valid_gases.items():  # Changed name to gas_name
        print(f"  {gas_name}: {data}")
