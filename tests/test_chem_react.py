# tests/test_chem_react.py
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_react.py"
# COOH → COCl
ACID_TO_CHLORIDE = "[C:1](=O)[OH]>>[C:1](=O)Cl"
# Ester hydrolysis: ester → COOH + alcohol
ESTER_HYDROLYSIS = "[C:1](=O)[O:2][C:3]>>[C:1](=O)O.[O:2][C:3]"
# Non-matching reaction (nitrogen to nothing — won't match aspirin)
NO_MATCH_RXN = "[NH2:1]>>[N+:1]"


def test_acid_to_chloride():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--reaction", ACID_TO_CHLORIDE])
    assert result["n_products"] >= 1
    assert any("Cl" in p for p in result["products"])


def test_no_match_returns_empty():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--reaction", NO_MATCH_RXN])
    assert result["n_products"] == 0
    assert result["products"] == []


def test_products_are_valid_smiles():
    from rdkit import Chem
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--reaction", ACID_TO_CHLORIDE])
    for smi in result["products"]:
        assert Chem.MolFromSmiles(smi) is not None


def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--reaction", ACID_TO_CHLORIDE])
    for f in ["reactant", "reaction", "n_products", "products"]:
        assert f in result


def test_products_deduplicated():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--reaction", ACID_TO_CHLORIDE])
    assert len(result["products"]) == len(set(result["products"]))
