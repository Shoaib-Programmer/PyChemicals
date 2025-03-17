"""
This module contains all the classes for the chemicals in the PyChemicals package.
Database handling is decoupled; if the user does not supply chemical properties explicitly,
the classes will retrieve them based on the chemical name from the set data loader.
"""

from typing import Dict, Any
import numpy as np
from .chemicals_db import get_valid_acids, get_valid_bases, get_valid_gases


class Chemical:
    """Base class for all chemicals."""

    def __init__(self, name: str, **kwargs: Dict[str, Any]):
        """
        Initialize a Chemical instance.

        Args:
            name (str): Name of the chemical.
            **kwargs: Arbitrary keyword arguments including:
                     concentration, volume, mass
        """
        self.name = name
        self.concentration = kwargs.get("concentration")
        self.volume = kwargs.get("volume")
        self.mass = kwargs.get("mass")
        self.validate()

    def validate(self) -> None:
        """Validate chemical properties."""
        if self.concentration is not None and self.concentration <= 0:
            raise ValueError(f"{self.name}: Concentration must be positive.")
        if self.volume is not None and self.volume <= 0:
            raise ValueError(f"{self.name}: Volume must be positive.")
        if self.mass is not None and self.mass <= 0:
            raise ValueError(f"{self.name}: Mass must be positive.")

    def __repr__(self) -> str:
        return f"{self.name}: {self.concentration} M, {self.volume} L"

    def calculate_moles(self) -> float:
        """
        Calculate moles from concentration and volume.

        Returns:
            float: Calculated moles.
        """
        if self.concentration is not None and self.volume is not None:
            return self.concentration * self.volume
        raise ValueError(f"Insufficient data to calculate moles for {self.name}.")


class Acid(Chemical):
    """Class representing an acid chemical."""

    def __init__(self, name: str, **kwargs: Dict[str, Any]):
        """
        Initialize an Acid instance.

        Args:
            name (str): Name of the acid.
            **kwargs: Arbitrary keyword arguments including:
                     concentration, volume, mass, properties
        """
        # If the user did not supply a properties dict, automatically load it.
        properties = kwargs.pop("properties", None)
        if properties is None:
            acids_data = get_valid_acids()
            if name not in acids_data:
                raise ValueError(f"Unknown acid: {name} in the current data source.")
            properties = acids_data[name]
        self.properties = properties

        super().__init__(name, **kwargs)
        self.proticity = properties["proticity"]
        self.ka = properties["Ka"]
        self.ka1 = properties.get("ka1")
        self.ka2 = properties.get("ka2")
        self.molar_mass = properties["molar_mass"]
        self.validate_acid()

    def calculate_mass(self) -> float:
        """Calculate the mass if not provided."""
        if self.mass:
            return self.mass
        if self.concentration is None or self.volume is None:
            raise ValueError(
                "Concentration and volume must be provided to calculate mass."
            )
        return self.concentration * self.volume * self.molar_mass

    def h_plus(self) -> float:
        """Calculate the concentration of H+ ions."""
        if self.ka > 1:
            return self.concentration
        if self.concentration is None:
            raise ValueError(
                "Concentration must be provided to calculate H+ concentration."
            )
        conc = self.concentration
        # Solve quadratic: [H+]^2 + ka*[H+] - ka*conc = 0
        a = 1
        b = self.ka
        c = -self.ka * conc
        discriminant = b**2 - 4 * a * c
        if discriminant < 0:
            raise ValueError(
                "Negative discriminant; check Ka and concentration values."
            )
        return (-b + np.sqrt(discriminant)) / (2 * a)

    @staticmethod
    def calculate_temp_dependence(
        k_initial: float, t_initial: float, t_final: float, delta_h: float
    ) -> float:
        """
        Calculate the temperature dependence of the equilibrium constant
        using the van't Hoff equation.
        """
        t_initial_k = t_initial + 273.15
        t_final_k = t_final + 273.15
        gas_constant = 8.314  # J/(mol·K)
        ln_k_ratio = (-delta_h / gas_constant) * ((1 / t_final_k) - (1 / t_initial_k))
        return k_initial * np.exp(ln_k_ratio)

    def validate_acid(self) -> None:
        """Validate acid properties using a tolerance for float comparisons."""
        if (
            self.concentration is not None
            and self.volume is not None
            and self.mass is not None
        ):
            expected_mass = self.concentration * self.volume * self.molar_mass
            if not np.isclose(self.mass, expected_mass, rtol=1e-3):
                raise ValueError(
                    f"{self.volume} L of {self.concentration}M {self.name} "
                    "(expected mass {expected_mass:.2f}g) "
                    f"does not match provided mass {self.mass}g."
                )

    def verify_acid_type(self) -> str:
        """Verify the type of acid based on proticity."""
        if self.proticity == 1:
            return "Monoprotic"
        if self.proticity == 2:
            return "Diprotic"
        return "Polyprotic"

    def pka(self) -> float:
        """Calculate the pKa value."""
        return -np.log10(self.ka)

    def ph(self) -> float:
        """Calculate the pH of the acid."""
        h_conc = self.h_plus()
        ph_value = -np.log10(h_conc)
        return 0.0 if np.isclose(ph_value, 0.0) else ph_value

    def moles(self) -> float:
        """Calculate the number of moles of the acid."""
        if self.concentration is not None and self.volume is not None:
            return self.calculate_moles()
        if self.mass is not None:
            return self.mass / self.molar_mass
        raise ValueError("Insufficient data to calculate moles.")

    def __repr__(self) -> str:
        acid_type = self.verify_acid_type()
        return (
            f"{self.name} ({acid_type}): {self.concentration} M, pKa = {self.pka():.2f}"
        )


class Base(Chemical):
    """Class representing a base chemical."""

    def __init__(self, name: str, **kwargs: Dict[str, Any]):
        """
        Initialize a Base instance.

        Args:
            name (str): Name of the base.
            **kwargs: Arbitrary keyword arguments including:
                     concentration, volume, mass, properties
        """
        properties = kwargs.pop("properties", None)
        if properties is None:
            bases_data = get_valid_bases()
            if name not in bases_data:
                raise ValueError(f"Unknown base: {name} in the current data source.")
            properties = bases_data[name]
        self.properties = properties

        super().__init__(name, **kwargs)
        self.proticity = properties["proticity"]
        self.kb = properties["Kb"]
        self.kb1 = properties.get("kb1")
        self.kb2 = properties.get("kb2")
        self.molar_mass = properties["molar_mass"]
        self.validate_base()

    def validate_base(self) -> None:
        """Validate base properties using a tolerance for float comparisons."""
        if (
            self.concentration is not None
            and self.volume is not None
            and self.mass is not None
        ):
            expected_mass = self.concentration * self.volume * self.molar_mass
            if not np.isclose(self.mass, expected_mass, rtol=1e-3):
                raise ValueError(
                    f"{self.volume} L of {self.concentration}M {self.name} "
                    "(expected mass {expected_mass:.2f}g) "
                    f"does not match provided mass {self.mass}g."
                )

    def verify_base_type(self) -> str:
        """Verify the type of base based on proticity."""
        if self.proticity == 1:
            return "Monoprotic"
        if self.proticity == 2:
            return "Diprotic"
        return "Polyprotic"

    def calculate_mass(self) -> float:
        """Calculate the mass if not provided."""
        if self.mass:
            return self.mass
        if self.concentration is None or self.volume is None:
            raise ValueError(
                "Concentration and volume must be provided to calculate mass."
            )
        return self.concentration * self.volume * self.molar_mass

    def oh_minus(self) -> float:
        """Calculate the concentration of OH- ions."""
        if self.kb > 1:
            return self.concentration  # Strong base: complete dissociation.
        if self.concentration is None:
            raise ValueError(
                "Concentration must be provided to calculate OH- concentration."
            )
        conc = self.concentration
        discriminant = self.kb**2 + 4 * self.kb * conc
        if discriminant < 0:
            raise ValueError(
                "Negative discriminant; check Kb and concentration values."
            )
        return (-self.kb + np.sqrt(discriminant)) / 2

    @staticmethod
    def calculate_temp_dependence(
        k_initial: float, t_initial: float, t_final: float, delta_h: float
    ) -> float:
        """Calculate temperature dependence using the van 't Hoff equation."""
        t_initial_k = t_initial + 273.15
        t_final_k = t_final + 273.15
        gas_constant = 8.314
        ln_k_ratio = (-delta_h / gas_constant) * ((1 / t_final_k) - (1 / t_initial_k))
        return k_initial * np.exp(ln_k_ratio)

    def pkb(self) -> float:
        """Calculate the pKb value."""
        return -np.log10(self.kb)

    def poh(self) -> float:
        """Calculate the pOH of the base."""
        if self.kb >= 1e-10:
            oh_conc = self.concentration * float(self.proticity)
            return -np.log10(oh_conc)
        if self.concentration is not None:
            oh_conc = np.sqrt(self.kb * self.concentration)
            return -np.log10(oh_conc)
        raise ValueError("Concentration must be provided to calculate pOH.")

    def ph(self) -> float:
        """Calculate the pH of the base."""
        return 14 - self.poh()

    def moles(self) -> float:
        """Calculate the number of moles of the base."""
        if self.concentration is not None and self.volume is not None:
            return self.calculate_moles()
        if self.mass is not None:
            return self.mass / self.molar_mass
        raise ValueError("Insufficient data to calculate moles.")

    def __repr__(self) -> str:
        base_type = self.verify_base_type()
        return (
            f"{self.name} ({base_type}): {self.concentration} M, pKb = {self.pkb():.2f}"
        )


class Gas(Chemical):
    """Class representing a gas."""

    R = 0.0821  # Ideal gas constant in L·atm/(K·mol)

    def __init__(self, name: str, **kwargs: Dict[str, Any]):
        """
        Initialize a Gas instance.

        Args:
            name (str): Name of the gas.
            **kwargs: Arbitrary keyword arguments including:
                     pressure, temperature, volume, mass, properties
        """
        properties = kwargs.pop("properties", None)
        if properties is None:
            gases_data = get_valid_gases()
            if name not in gases_data:
                raise ValueError(f"Unknown gas: {name} in the current data source.")
            properties = gases_data[name]
        self.properties = properties

        super().__init__(name, **kwargs)
        self.pressure = kwargs.get("pressure")
        self.temperature = kwargs.get("temperature")

    @property
    def moles(self) -> float:
        """Calculate the number of moles of the gas using the ideal gas law."""
        if (
            self.pressure is not None
            and self.temperature is not None
            and self.volume is not None
        ):
            return (self.pressure * self.volume) / (self.R * self.temperature)
        raise ValueError("Insufficient data to calculate moles for gas.")

    def volume_from_moles(self, moles: float) -> float:
        """Calculate volume using the ideal gas law: V = nRT/P."""
        if self.pressure is not None and self.temperature is not None:
            return (moles * self.R * self.temperature) / self.pressure
        raise ValueError(
            "Pressure and temperature must be defined to calculate volume."
        )

    def moles_from_volume(self, volume: float) -> float:
        """Calculate moles using the ideal gas law: n = PV/RT."""
        if self.pressure is not None and self.temperature is not None:
            return (self.pressure * volume) / (self.R * self.temperature)
        raise ValueError("Pressure and temperature must be defined to calculate moles.")

    def density(self) -> float:
        """Calculate the density of the gas in g/L."""
        if self.volume is None or self.volume == 0:
            raise ValueError("Volume must be greater than zero to calculate density.")
        if self.mass is None:
            raise ValueError("Mass must be provided to calculate density.")
        return self.mass / self.volume

    def molar_mass_calc(self) -> float:
        """Calculate the molar mass of the gas in g/mol."""
        n = self.moles
        if n == 0:
            raise ValueError("Moles must be greater than zero to calculate molar mass.")
        if self.mass is None:
            raise ValueError("Mass must be provided to calculate molar mass.")
        return self.mass / n

    def partial_pressure(self, total_pressure: float, mole_fraction: float) -> float:
        """Calculate the partial pressure of the gas in a mixture."""
        return total_pressure * mole_fraction

    def grahams_law_effusion(self, molar_mass_other: float) -> float:
        """Calculate the rate of effusion relative to another gas using Graham's Law."""
        return (self.molar_mass_calc() / molar_mass_other) ** 0.5

    def compressibility_factor(
        self, pressure: float, volume: float, temperature: float
    ) -> float:
        """Calculate the compressibility factor Z of the gas: Z = PV/nRT."""
        n = self.moles
        return (pressure * volume) / (n * self.R * temperature)

    def work_done(self, v_initial: float, v_final: float) -> float:
        """Calculate the work done by the gas during isothermal expansion in L·atm."""
        n = self.moles
        return -n * self.R * self.temperature * np.log(v_final / v_initial)

    def __repr__(self) -> str:
        return (
            f"{self.name}: {self.pressure} atm, {self.temperature} K, {self.volume} L"
        )
