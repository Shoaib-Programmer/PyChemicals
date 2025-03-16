"""
This module handles database access for chemical data.
It decouples database operations from the chemical logic.
"""

import os
import csv
from typing import Dict, Any, List, Union
from cs50 import SQL

# Global variables to hold the data loader and source type.
_data_loader = None  # type: Union[None, callable] # pylint: disable=invalid-name
_source_type = None  # "sql" or "csv" # pylint: disable=invalid-name


def use_source(source: Union[str, os.PathLike]) -> None:
    """
    Set the data source for chemical data based on the provided file or directory.

    Compatible extensions: .db, .sqlite3, .sql for SQL sources,
                           .csv for CSV files or a directory containing CSV files.

    Args:
        source: A string or os.PathLike object representing a file or directory path.
    """
    global _data_loader, _source_type  # pylint: disable=global-statement
    source = os.fspath(source)

    if os.path.isdir(source):
        # CSV directory mode: Expect files like acids.csv, acid_temperatures.csv, bases.csv,
        # base_temperatures.csv, gases.csv.
        _source_type = "csv"

        def loader(
            table_name: str,
            query_params: List[Any] = None,  # pylint: disable=unused-argument
        ) -> List[Dict[str, Any]]:
            filename = os.path.join(source, f"{table_name}.csv")
            if not os.path.exists(filename):
                raise FileNotFoundError(
                    f"Expected CSV file {filename} not found in directory {source}."
                )
            with open(filename, newline="", encoding="utf-8") as f:
                return list(csv.DictReader(f))

        _data_loader = loader

    else:
        # Single file mode. Check extension.
        ext = os.path.splitext(source)[1].lower()
        if ext in (".db", ".sqlite3", ".sql"):
            # SQL file mode.
            connection = SQL(f"sqlite:///{source}")
            _source_type = "sql"

            def loader(
                table_name: str,
                query_params: List[Any] = None,  # pylint: disable=unused-argument
            ) -> List[Dict[str, Any]]:
                if table_name == "acids":
                    return connection.execute("SELECT * FROM acids")
                if table_name == "acid_temperatures":
                    return connection.execute("SELECT * FROM acid_temperatures")
                if table_name == "bases":
                    return connection.execute("SELECT * FROM bases")
                if table_name == "base_temperatures":
                    return connection.execute("SELECT * FROM base_temperatures")
                if table_name == "gases":
                    return connection.execute("SELECT * FROM gases")
                raise ValueError(f"Unknown table: {table_name}")

            _data_loader = loader

        elif ext == ".csv":
            # Single CSV file mode. We assume that the file name (without path) indicates the table.
            _source_type = "csv"

            def loader(
                table_name: str,
                query_params: List[Any] = None,  # pylint: disable=unused-argument
            ) -> List[Dict[str, Any]]:
                basename = os.path.basename(source).lower()
                if not basename.startswith(table_name.lower()):
                    raise ValueError(
                        f"CSV file {source} does not appear to be for table '{table_name}'."
                    )
                with open(source, newline="", encoding="utf-8") as f:
                    return list(csv.DictReader(f))

            _data_loader = loader

        else:
            raise ValueError(f"Unsupported file extension: '{ext}'")


def _load_table(table_name: str) -> List[Dict[str, Any]]:
    """
    Helper function to load a table using the global data loader.

    Args:
        table_name: The name of the table to load.

    Returns:
        A list of dictionaries representing rows.
    """
    if _data_loader is None:
        raise RuntimeError("No data source has been set. Call use_source() first.")
    return _data_loader(table_name)


def get_valid_acids() -> Dict[str, Dict[str, Any]]:
    """
    Retrieve acids data from the data source and return a dictionary.

    Returns:
        A dictionary keyed by acid name.
    """
    try:
        rows = _load_table("acids")
    except Exception as e:
        raise RuntimeError(f"Error retrieving acids: {e}") from e
    result = {}
    for row in rows:
        acid_id = row["id"]
        chem_name = row["name"]
        temp_rows = _load_table("acid_temperatures")
        ka_temp = {
            t["temperature"]: t["ka"]
            for t in temp_rows
            if t.get("acid_id") == acid_id or t.get("acid_id") == str(acid_id)
        }
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
    """
    Retrieve bases data from the data source and return a dictionary.

    Returns:
        A dictionary keyed by base name.
    """
    try:
        rows = _load_table("bases")
    except Exception as e:
        raise RuntimeError(f"Error retrieving bases: {e}") from e
    result = {}
    for row in rows:
        base_id = row["id"]
        chem_name = row["name"]
        temp_rows = _load_table("base_temperatures")
        kb_temp = {
            t["temperature"]: t["kb"]
            for t in temp_rows
            if t.get("base_id") == base_id or t.get("base_id") == str(base_id)
        }
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
    """
    Retrieve gases data from the data source and return a dictionary.

    Returns:
        A dictionary keyed by gas name.
    """
    try:
        rows = _load_table("gases")
    except Exception as e:
        raise RuntimeError(f"Error retrieving gases: {e}") from e
    result = {}
    for row in rows:
        chem_name = row["name"]
        result[chem_name] = {
            "molar_mass": row["molar_mass"],
            "density": row[
                "density"
            ],  # Reference density (typically at standard conditions)
        }
    return result


# Default behavior: if no use_source() is called explicitly, use the built-in chemicals.db file.
if _data_loader is None:
    default_db_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "chemicals.db"
    )
    use_source(default_db_path)

if __name__ == "__main__":
    try:
        acids = get_valid_acids()
        bases = get_valid_bases()
        gases = get_valid_gases()
        print("Acids:", acids)
        print("Bases:", bases)
        print("Gases:", gases)
    except Exception as e:  # pylint: disable=broad-exception-caught
        print("Error:", e)
