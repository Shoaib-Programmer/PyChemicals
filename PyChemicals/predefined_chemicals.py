"""
Module containing predefined chemical data and database setup functionality.
This module handles initialization and population of the SQLite database
with predefined acids, bases, and gases.
"""

from cs50 import SQL  # type: ignore

# Predefined chemical data with IUPAC names
valid_acids = {
    "Hydrogen chloride": {
        "proticity": 1,
        "Ka": 1e6,
        "molar_mass": 36.461,
        "Ka_temp": {"25C": 1e6, "50C": 2e6},
        "ka1": 1e6,
        "ka2": None,
    },
    "Nitric acid": {
        "proticity": 1,
        "Ka": 1e3,
        "molar_mass": 63.012,
        "Ka_temp": {"25C": 1e3, "50C": 1.5e3},
        "ka1": 1e3,
        "ka2": None,
    },
    "Ethanoic acid": {
        "proticity": 1,
        "Ka": 1.8e-5,
        "molar_mass": 60.052,
        "Ka_temp": {"25C": 1.8e-5, "50C": 2.1e-5},
        "ka1": 1.8e-5,
        "ka2": None,
    },
    "Hydrogen bromide": {
        "proticity": 1,
        "Ka": 1e6,
        "molar_mass": 80.911,
        "Ka_temp": {"25C": 1e6, "50C": 1.9e6},
        "ka1": 1e6,
        "ka2": None,
    },
    "Hydrogen iodide": {
        "proticity": 1,
        "Ka": 1e6,
        "molar_mass": 127.913,
        "Ka_temp": {"25C": 1e6, "50C": 2e6},
        "ka1": 1e6,
        "ka2": None,
    },
    "Methanoic acid": {
        "proticity": 1,
        "Ka": 1.77e-4,
        "molar_mass": 46.025,
        "Ka_temp": {"25C": 1.77e-4, "50C": 2.0e-4},
        "ka1": 1.77e-4,
        "ka2": None,
    },
    "Sulfuric acid": {
        "proticity": 2,
        "Ka": 1e3,
        "molar_mass": 98.079,
        "Ka_temp": {"25C": 1e3, "50C": 1.2e3},
        "ka1": 1e3,
        "ka2": 1.2e-2,
    },
    "Carbonic acid": {
        "proticity": 2,
        "Ka": 4.3e-7,
        "molar_mass": 62.033,
        "Ka_temp": {"25C": 4.3e-7, "50C": 5e-7},
        "ka1": 4.3e-7,
        "ka2": 5.6e-11,
    },
    "Ethanedioic acid": {
        "proticity": 2,
        "Ka": 5.37e-2,
        "molar_mass": 90.031,
        "Ka_temp": {"25C": 5.37e-2, "50C": 6.0e-2},
        "ka1": 5.37e-2,
        "ka2": 6.4e-5,
    },
    "Phosphoric acid": {
        "proticity": 3,
        "Ka": 7.5e-3,
        "molar_mass": 97.994,
        "Ka_temp": {"25C": 7.5e-3, "50C": 9e-3},
        "ka1": 7.5e-3,
        "ka2": 6.2e-8,
    },
    "2-hydroxypropane-1,2,3-tricarboxylic acid": {
        "proticity": 3,
        "Ka": 7.4e-4,
        "molar_mass": 192.124,
        "Ka_temp": {"25C": 7.4e-4, "50C": 8.5e-4},
        "ka1": 7.4e-4,
        "ka2": 1.7e-5,
    },
    "Arsenic acid": {
        "proticity": 3,
        "Ka": 5.5e-3,
        "molar_mass": 177.833,
        "Ka_temp": {"25C": 5.5e-3, "50C": 6.1e-3},
        "ka1": 5.5e-3,
        "ka2": 1.2e-7,
    },
}

valid_bases = {
    "Sodium hydroxide": {
        "proticity": 1,
        "Kb": 1e14,
        "molar_mass": 40.00,
        "Kb_temp": {"25C": 1e14, "50C": 1.1e14},
        "kb1": 1e14,
        "kb2": None,
    },
    "Potassium hydroxide": {
        "proticity": 1,
        "Kb": 1e14,
        "molar_mass": 56.105,
        "Kb_temp": {"25C": 1e14, "50C": 1.1e14},
        "kb1": 1e14,
        "kb2": None,
    },
    "Ammonia": {
        "proticity": 1,
        "Kb": 1.76e-5,
        "molar_mass": 17.031,
        "Kb_temp": {"25C": 1.76e-5, "50C": 2.0e-5},
        "kb1": 1.76e-5,
        "kb2": None,
    },
    "Calcium hydroxide": {
        "proticity": 2,
        "Kb": 5.5e-2,
        "molar_mass": 74.092,
        "Kb_temp": {"25C": 5.5e-2, "50C": 6.0e-2},
        "kb1": 5.5e-2,
        "kb2": 1.3e-5,
    },
    "Barium hydroxide": {
        "proticity": 2,
        "Kb": 5.5e-2,
        "molar_mass": 171.34,
        "Kb_temp": {"25C": 5.5e-2, "50C": 6.5e-2},
        "kb1": 5.5e-2,
        "kb2": 1.5e-5,
    },
    "Methanamine": {
        "proticity": 1,
        "Kb": 4.38e-4,
        "molar_mass": 31.057,
        "Kb_temp": {"25C": 4.38e-4, "50C": 5.0e-4},
        "kb1": 4.38e-4,
        "kb2": None,
    },
}

# Add some sample valid gases (names remain unchanged)
valid_gases = {
    "Carbon dioxide": {
        "molar_mass": 44.01,
        "density": 1.98,  # example density (in appropriate units)
    },
    "Oxygen": {
        "molar_mass": 32.00,
        "density": 1.43,
    },
    "Nitrogen": {
        "molar_mass": 28.02,
        "density": 1.25,
    },
    "Methane": {
        "molar_mass": 16.04,
        "density": 0.656,
    },
}

# Connect to (or create) the database using cs50's SQL module.
# The connection string "sqlite:///chemicals.db" creates the file in your current working directory.
db = SQL("sqlite:///chemicals.db")


def create_tables():
    """
    Create tables in the database if they don't already exist.
    """
    # Create acids table
    db.execute("""
        CREATE TABLE IF NOT EXISTS acids (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            proticity INTEGER,
            Ka REAL,
            molar_mass REAL,
            ka1 REAL,
            ka2 REAL
        )
    """)
    # Create acid temperature table
    db.execute("""
        CREATE TABLE IF NOT EXISTS acid_temperatures (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            acid_id INTEGER,
            temperature TEXT,
            ka REAL,
            FOREIGN KEY (acid_id) REFERENCES acids(id)
        )
    """)
    # Create bases table
    db.execute("""
        CREATE TABLE IF NOT EXISTS bases (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            proticity INTEGER,
            Kb REAL,
            molar_mass REAL,
            kb1 REAL,
            kb2 REAL
        )
    """)
    # Create base temperature table
    db.execute("""
        CREATE TABLE IF NOT EXISTS base_temperatures (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            base_id INTEGER,
            temperature TEXT,
            kb REAL,
            FOREIGN KEY (base_id) REFERENCES bases(id)
        )
    """)
    # Create gases table
    db.execute("""
        CREATE TABLE IF NOT EXISTS gases (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE,
            molar_mass REAL,
            density REAL
        )
    """)


def insert_acids(acids):
    """
    Insert acid data into the acids table.

    Args:
        acids (dict): Dictionary containing acid data with acid names as keys.
    """
    for name, data in acids.items():
        db.execute(
            """
            INSERT OR IGNORE INTO acids (name, proticity, Ka, molar_mass, ka1, ka2)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            name,
            data["proticity"],
            data["Ka"],
            data["molar_mass"],
            data["ka1"],
            data["ka2"],
        )
        result = db.execute("SELECT id FROM acids WHERE name = ?", name)
        acid_id = result[0]["id"]
        for temp, ka in data.get("Ka_temp", {}).items():
            db.execute(
                """
                INSERT OR IGNORE INTO acid_temperatures (acid_id, temperature, ka)
                VALUES (?, ?, ?)
                """,
                acid_id,
                temp,
                ka,
            )


def insert_bases(bases):
    """
    Insert base data into the bases table.

    Args:
        bases (dict): Dictionary containing base data with base names as keys.
    """
    for name, data in bases.items():
        db.execute(
            """
            INSERT OR IGNORE INTO bases (name, proticity, Kb, molar_mass, kb1, kb2)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            name,
            data["proticity"],
            data["Kb"],
            data["molar_mass"],
            data["kb1"],
            data["kb2"],
        )
        result = db.execute("SELECT id FROM bases WHERE name = ?", name)
        base_id = result[0]["id"]
        for temp, kb in data.get("Kb_temp", {}).items():
            db.execute(
                """
                INSERT OR IGNORE INTO base_temperatures (base_id, temperature, kb)
                VALUES (?, ?, ?)
                """,
                base_id,
                temp,
                kb,
            )


def insert_gases(gases):
    """
    Insert gas data into the gases table.

    Args:
        gases (dict): Dictionary containing gas data with gas names as keys.
    """
    for name, data in gases.items():
        db.execute(
            """
            INSERT OR IGNORE INTO gases (name, molar_mass, density)
            VALUES (?, ?, ?)
            """,
            name,
            data["molar_mass"],
            data["density"],
        )


def main():
    """
    Initialize database tables and populate them with predefined chemical data.
    Creates tables if they don't exist and inserts predefined acids, bases, and gases.
    """
    create_tables()
    insert_acids(valid_acids)
    insert_bases(valid_bases)
    insert_gases(valid_gases)
    print(
        "Database setup complete using the cs50 library. Data inserted into chemicals.db."
    )


if __name__ == "__main__":
    main()
