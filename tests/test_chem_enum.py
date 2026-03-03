# tests/test_chem_enum.py
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_enum.py"
IBUPROFEN  = "CC(Cc1ccc(CC(C)C(=O)O)cc1)C(=O)O"   # 2 stereocenters (unassigned)
ALANINE    = "N[C@@H](C)C(=O)O"                     # 1 assigned stereocenter


def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for f in ["input_smiles", "stereocenters_in_input", "n_isomers", "isomers"]:
        assert f in result


def test_aspirin_no_stereocenters():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["stereocenters_in_input"] == 0
    assert result["n_isomers"] == 1  # only itself


def test_ibuprofen_two_centers():
    result = run_script(SCRIPT, ["--smiles", IBUPROFEN])
    assert result["stereocenters_in_input"] >= 1
    assert result["n_isomers"] >= 2


def test_alanine_one_assigned():
    result = run_script(SCRIPT, ["--smiles", ALANINE])
    assert result["stereocenters_in_input"] == 1
    # already fully specified → only 1 isomer
    assert result["n_isomers"] == 1


def test_isomers_are_valid_smiles():
    from rdkit import Chem
    result = run_script(SCRIPT, ["--smiles", IBUPROFEN])
    for smi in result["isomers"]:
        assert Chem.MolFromSmiles(smi) is not None


def test_isomers_deduplicated():
    result = run_script(SCRIPT, ["--smiles", IBUPROFEN])
    assert len(result["isomers"]) == len(set(result["isomers"]))
