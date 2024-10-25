import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from chemicals import Chemical, Acid, Base
from titration_calculation import volume_of_titrant, Analyte, Titrant

matplotlib.use("tkAgg")  # or use 'Qt5Agg'

def titrate_curve_monoprotic(analyte: Analyte, titrant: Titrant):
    if analyte.volume is None:
        raise ValueError("Analyte must have a defined volume.")
    if analyte.concentration is None:
        raise ValueError("Analyte must have a concentration.")
    if titrant.concentration is None:
        raise ValueError("Titrant must have a concentration.")
    if titrant.volume is None:
        raise ValueError("Titrant must have a defined volume.")

    # Calculate volume of titrant needed to reach equivalence point
    V_eqv = volume_of_titrant(analyte, titrant)

    # Define volumes to add titrant, from 0 to a bit beyond equivalence point
    V_added = np.linspace(0, 2 * V_eqv, 500)

    pH = []

    if isinstance(analyte, Acid):
        # Strong Acid - Strong Base Titration
        C_acid = analyte.concentration
        C_base = titrant.concentration
        V_acid = analyte.volume

        for V_titrant in V_added:
            moles_acid = C_acid * V_acid
            moles_base = C_base * V_titrant

            total_volume = V_acid + V_titrant

            if moles_base < moles_acid:
                # Before equivalence point (excess acid)
                H_plus = (moles_acid - moles_base) / total_volume
                pH.append(-np.log10(H_plus))
            elif moles_base == moles_acid:
                # At equivalence point (neutral, pH = 7 for strong acid-strong base titration)
                pH.append(7)
            else:
                # After equivalence point (excess base)
                OH_minus = (moles_base - moles_acid) / total_volume
                pH.append(14 + np.log10(OH_minus))

    elif isinstance(analyte, Base):
        # Strong Base - Strong Acid Titration
        C_base = analyte.concentration
        C_acid = titrant.concentration
        V_base = analyte.volume

        for V_titrant in V_added:
            moles_base = C_base * V_base
            moles_acid = C_acid * V_titrant

            total_volume = V_base + V_titrant

            if moles_acid < moles_base:
                # Before equivalence point (excess base)
                OH_minus = (moles_base - moles_acid) / total_volume
                pH.append(14 + np.log10(OH_minus))
            elif moles_acid == moles_base:
                # At equivalence point (neutral, pH = 7)
                pH.append(7)
            else:
                # After equivalence point (excess acid)
                H_plus = (moles_acid - moles_base) / total_volume
                pH.append(-np.log10(H_plus))

    # Plot the titration curve
    plt.figure(figsize=(8, 5))
    plt.plot(V_added * 1000, pH, label=f"{analyte.volume * 1000}mL {analyte.name} ({C_acid} M) titrated with {titrant.name} ({C_base} M)")
    plt.axvline(x=V_eqv * 1000, color='r', linestyle='--', label='Equivalence Point')
    plt.axhline(y=7, color='g', linestyle='--', label='pH = 7 at Equivalence')
    plt.xlabel(f'Volume of {titrant.name} added (mL)')
    plt.ylabel('pH')
    plt.title('Titration Curve')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    hcl = Acid(name="Hydrochloric Acid", concentration=0.1, volume=0.2)
    naoh = Base(name="Sodium Hydroxide", concentration=0.2, volume=0.1)
    titrate_curve_monoprotic(hcl, naoh)
