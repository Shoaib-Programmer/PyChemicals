valid_acids = {
    "Hydrochloric Acid": {
        "proticity": 1, 
        "Ka": 1e6, 
        "molar_mass": 36.461, 
        "Ka_temp": {"25C": 1e6, "50C": 2e6},  # Temperature dependence
        "ka1": 1e6, "ka2": None
    },  # Strong acid
    "Nitric Acid": {
        "proticity": 1, 
        "Ka": 1e3, 
        "molar_mass": 63.012, 
        "Ka_temp": {"25C": 1e3, "50C": 1.5e3},
        "ka1": 1e3, "ka2": None
    },  # Strong acid
    "Acetic Acid": {
        "proticity": 1, 
        "Ka": 1.8e-5, 
        "molar_mass": 60.052, 
        "Ka_temp": {"25C": 1.8e-5, "50C": 2.1e-5},
        "ka1": 1.8e-5, "ka2": None
    },  # Weak acid
    "Hydrobromic Acid": {
        "proticity": 1, 
        "Ka": 1e6, 
        "molar_mass": 80.911, 
        "Ka_temp": {"25C": 1e6, "50C": 1.9e6},
        "ka1": 1e6, "ka2": None
    },  # Strong acid
    "Hydroiodic Acid": {
        "proticity": 1, 
        "Ka": 1e6, 
        "molar_mass": 127.913, 
        "Ka_temp": {"25C": 1e6, "50C": 2e6},
        "ka1": 1e6, "ka2": None
    },  # Strong acid
    "Formic Acid": {
        "proticity": 1, 
        "Ka": 1.77e-4, 
        "molar_mass": 46.025, 
        "Ka_temp": {"25C": 1.77e-4, "50C": 2.0e-4},
        "ka1": 1.77e-4, "ka2": None
    },  # Weak acid
    "Sulfuric Acid": {
        "proticity": 2, 
        "Ka": 1e3, 
        "molar_mass": 98.079, 
        "Ka_temp": {"25C": 1e3, "50C": 1.2e3},
        "ka1": 1e3, "ka2": 1.2e-2
    },  # Strong acid (1st dissociation)
    "Carbonic Acid": {
        "proticity": 2, 
        "Ka": 4.3e-7, 
        "molar_mass": 62.033, 
        "Ka_temp": {"25C": 4.3e-7, "50C": 5e-7},
        "ka1": 4.3e-7, "ka2": 5.6e-11
    },  # Weak acid
    "Oxalic Acid": {
        "proticity": 2, 
        "Ka": 5.37e-2, 
        "molar_mass": 90.031, 
        "Ka_temp": {"25C": 5.37e-2, "50C": 6.0e-2},
        "ka1": 5.37e-2, "ka2": 6.4e-5
    },  # Weak acid
    "Phosphoric Acid": {
        "proticity": 3, 
        "Ka": 7.5e-3, 
        "molar_mass": 97.994, 
        "Ka_temp": {"25C": 7.5e-3, "50C": 9e-3},
        "ka1": 7.5e-3, "ka2": 6.2e-8
    },  # Weak acid
    "Citric Acid": {
        "proticity": 3, 
        "Ka": 7.4e-4, 
        "molar_mass": 192.124, 
        "Ka_temp": {"25C": 7.4e-4, "50C": 8.5e-4},
        "ka1": 7.4e-4, "ka2": 1.7e-5
    },  # Weak acid
    "Arsenic Acid": {
        "proticity": 3, 
        "Ka": 5.5e-3, 
        "molar_mass": 177.833, 
        "Ka_temp": {"25C": 5.5e-3, "50C": 6.1e-3},
        "ka1": 5.5e-3, "ka2": 1.2e-7
    }  # Weak acid
}

valid_bases = {
    "Sodium Hydroxide": {
        "proticity": 1, 
        "Kb": 1e14, 
        "molar_mass": 40.00, 
        "Kb_temp": {"25C": 1e14, "50C": 1.1e14},
        "kb1": 1e14, "kb2": None
    },  # Strong base
    "Potassium Hydroxide": {
        "proticity": 1, 
        "Kb": 1e14, 
        "molar_mass": 56.105, 
        "Kb_temp": {"25C": 1e14, "50C": 1.1e14},
        "kb1": 1e14, "kb2": None
    },  # Strong base
    "Ammonia": {
        "proticity": 1, 
        "Kb": 1.76e-5, 
        "molar_mass": 17.031, 
        "Kb_temp": {"25C": 1.76e-5, "50C": 2.0e-5},
        "kb1": 1.76e-5, "kb2": None
    },  # Weak base
    "Calcium Hydroxide": {
        "proticity": 2, 
        "Kb": 5.5e-2, 
        "molar_mass": 74.092, 
        "Kb_temp": {"25C": 5.5e-2, "50C": 6.0e-2},
        "kb1": 5.5e-2, "kb2": 1.3e-5
    },  # Strong base
    "Barium Hydroxide": {
        "proticity": 2, 
        "Kb": 5.5e-2, 
        "molar_mass": 171.34, 
        "Kb_temp": {"25C": 5.5e-2, "50C": 6.5e-2},
        "kb1": 5.5e-2, "kb2": 1.5e-5
    },  # Strong base
    "Methylamine": {
        "proticity": 1, 
        "Kb": 4.38e-4, 
        "molar_mass": 31.057, 
        "Kb_temp": {"25C": 4.38e-4, "50C": 5.0e-4},
        "kb1": 4.38e-4, "kb2": None
    }  # Weak base
}

valid_gases = {}

if __name__ == "__main__":
    print(valid_acids, valid_bases)