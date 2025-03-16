"""
Test module for using different data sources, in this case `.csv`.
"""
import csv
import pytest
from PyChemicals.chemicals_db import use_source, get_valid_acids, get_valid_bases, get_valid_gases

@pytest.fixture
def csv_data_source(tmp_path): # pylint: disable=redefined-outer-name
    """
    Creates a temporary directory with CSV files for acids, acid_temperatures,
    bases, base_temperatures, and gases.
    """
    # Define file paths.
    acids_csv = tmp_path / "acids.csv"
    acid_temps_csv = tmp_path / "acid_temperatures.csv"
    bases_csv = tmp_path / "bases.csv"
    base_temps_csv = tmp_path / "base_temperatures.csv"
    gases_csv = tmp_path / "gases.csv"

    # Create sample acids.csv.
    with acids_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["id", "name", "proticity", "Ka", "molar_mass", "ka1", "ka2"]
        )
        writer.writeheader()
        writer.writerow({
            "id": "1",
            "name": "Test Acid",
            "proticity": "1",
            "Ka": "1e-3",
            "molar_mass": "36.46",
            "ka1": "1e-3",
            "ka2": "0"
        })

    # Create sample acid_temperatures.csv.
    with acid_temps_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["acid_id", "temperature", "ka"])
        writer.writeheader()
        writer.writerow({
            "acid_id": "1",
            "temperature": "25",
            "ka": "1e-3"
        })

    # Create sample bases.csv.
    with bases_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["id", "name", "proticity", "Kb", "molar_mass", "kb1", "kb2"]
        )
        writer.writeheader()
        writer.writerow({
            "id": "1",
            "name": "Test Base",
            "proticity": "1",
            "Kb": "1e-5",
            "molar_mass": "40.00",
            "kb1": "1e-5",
            "kb2": "0"
        })

    # Create sample base_temperatures.csv.
    with base_temps_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["base_id", "temperature", "kb"])
        writer.writeheader()
        writer.writerow({
            "base_id": "1",
            "temperature": "25",
            "kb": "1e-5"
        })

    # Create sample gases.csv.
    with gases_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "molar_mass", "density"])
        writer.writeheader()
        writer.writerow({
            "name": "Test Gas",
            "molar_mass": "28.01",
            "density": "1.25"
        })

    return tmp_path

def test_csv_source(csv_data_source):  # pylint: disable=redefined-outer-name
    """
    Test the CSV data source functionality by verifying acid, base, and gas data retrieval.
    """
    # Set the CSV directory as the data source.
    use_source(csv_data_source)

    # Retrieve data from the CSV source.
    acids = get_valid_acids()
    bases = get_valid_bases()
    gases = get_valid_gases()

    # Verify acids data.
    assert "Test Acid" in acids
    acid = acids["Test Acid"]
    assert acid["proticity"] == "1"
    assert acid["Ka"] == "1e-3"
    assert acid["molar_mass"] == "36.46"
    assert "25" in acid["Ka_temp"]
    assert acid["Ka_temp"]["25"] == "1e-3"

    # Verify bases data.
    assert "Test Base" in bases
    base = bases["Test Base"]
    assert base["proticity"] == "1"
    assert base["Kb"] == "1e-5"
    assert base["molar_mass"] == "40.00"
    assert "25" in base["Kb_temp"]
    assert base["Kb_temp"]["25"] == "1e-5"

    # Verify gases data.
    assert "Test Gas" in gases
    gas = gases["Test Gas"]
    assert gas["molar_mass"] == "28.01"
    assert gas["density"] == "1.25"
