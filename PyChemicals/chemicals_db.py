"""
This module handles database access for chemical data.
It decouples database operations from the chemical logic.
"""

import os
from typing import Dict, Any
from cs50 import SQL

# Get the absolute path to the current file's directory.
module_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.join(module_dir, "..")
db_path = os.path.join(project_root, "chemicals.db")
db = SQL(f"sqlite:///{db_path}")

def get_valid_acids() -> Dict[str, Dict[str, Any]]:
    """Retrieve acids data from the database and return a dictionary."""
    try:
        rows = db.execute("SELECT * FROM acids")
    except Exception as e:
        raise RuntimeError(f"Error retrieving acids from database: {e}") from e
    result = {}
    for row in rows:
        acid_id = row["id"]
        chem_name = row["name"]
        temp_rows = db.execute(
            "SELECT temperature, ka FROM acid_temperatures WHERE acid_id = ?", acid_id
        )
        ka_temp = {t["temperature"]: t["ka"] for t in temp_rows}
        result[chem_name] = {
            "proticity": row["proticity"],
            "Ka": row["Ka"],
            "molar_mass": row["molar_mass"],
            "ka1": row["ka1"],
            "ka2": row["ka2"],
            "Ka_temp": ka_temp,
        }
    return result

def get_valid_bases() -> Dict[str, Dict[str, Any]]:
    """Retrieve bases data from the database and return a dictionary."""
    try:
        rows = db.execute("SELECT * FROM bases")
    except Exception as e:
        raise RuntimeError(f"Error retrieving bases from database: {e}") from e
    result = {}
    for row in rows:
        base_id = row["id"]
        chem_name = row["name"]
        temp_rows = db.execute(
            "SELECT temperature, kb FROM base_temperatures WHERE base_id = ?", base_id
        )
        kb_temp = {t["temperature"]: t["kb"] for t in temp_rows}
        result[chem_name] = {
            "proticity": row["proticity"],
            "Kb": row["Kb"],
            "molar_mass": row["molar_mass"],
            "kb1": row["kb1"],
            "kb2": row["kb2"],
            "Kb_temp": kb_temp,
        }
    return result

def get_valid_gases() -> Dict[str, Dict[str, Any]]:
    """Retrieve gases data from the database and return a dictionary."""
    try:
        rows = db.execute("SELECT * FROM gases")
    except Exception as e:
        raise RuntimeError(f"Error retrieving gases from database: {e}") from e
    result = {}
    for row in rows:
        chem_name = row["name"]
        result[chem_name] = {
            "molar_mass": row["molar_mass"],
            "density": row["density"],  # Reference density (typically at standard conditions)
        }
    return result

if __name__ == "__main__":
    acids = get_valid_acids()
    bases = get_valid_bases()
    gases = get_valid_gases()
    print("Acids:", acids)
    print("Bases:", bases)
    print("Gases:", gases)
