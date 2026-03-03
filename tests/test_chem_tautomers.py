# tests/test_chem_tautomers.py
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_tautomers.py"
PHENOL = "Oc1ccccc1"          # keto-enol tautomerism
ACETYLACETONE = "CC(=O)CC(=O)C"  # classic keto-enol


def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for f in ["input_smiles", "canonical_tautomer", "n_tautomers", "tautomers"]:
        assert f in result


def test_tautomers_is_list():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert isinstance(result["tautomers"], list)
    assert result["n_tautomers"] == len(result["tautomers"])


def test_canonical_in_tautomers():
    result = run_script(SCRIPT, ["--smiles", PHENOL])
    assert result["canonical_tautomer"] in result["tautomers"]


def test_acetylacetone_multiple_tautomers():
    result = run_script(SCRIPT, ["--smiles", ACETYLACETONE])
    assert result["n_tautomers"] >= 2


def test_aspirin_canonical_valid():
    from rdkit import Chem
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert Chem.MolFromSmiles(result["canonical_tautomer"]) is not None
