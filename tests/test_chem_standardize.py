# tests/test_chem_standardize.py
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_standardize.py"

# Sodium benzoate: salt + charged
SODIUM_BENZOATE = "[Na+].[O-]C(=O)c1ccccc1"
# Acetate anion: needs neutralization
ACETATE_ANION = "CC(=O)[O-]"


def test_desalt_removes_counterion():
    result = run_script(SCRIPT, ["--smiles", SODIUM_BENZOATE])
    assert result["valid"] is True
    assert "desalted" in result["changes"]
    # Na should be gone
    assert "Na" not in result["output_smiles"]


def test_neutralize_anion():
    result = run_script(SCRIPT, ["--smiles", ACETATE_ANION])
    assert "neutralized" in result["changes"]
    assert result["changed"] is True
    # Neutral acetic acid
    assert "[O-]" not in result["output_smiles"]


def test_salt_and_charge_both_fixed():
    result = run_script(SCRIPT, ["--smiles", SODIUM_BENZOATE])
    assert "desalted" in result["changes"]
    assert "neutralized" in result["changes"]
    assert result["changed"] is True


def test_clean_molecule_unchanged():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["changes"] == []
    assert result["changed"] is False
    assert result["valid"] is True


def test_output_has_required_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for field in ["input_smiles", "output_smiles", "changes", "changed", "valid"]:
        assert field in result, f"Missing field: {field}"


def test_output_smiles_is_valid():
    """Output SMILES must be parseable by RDKit."""
    from rdkit import Chem
    result = run_script(SCRIPT, ["--smiles", SODIUM_BENZOATE])
    mol = Chem.MolFromSmiles(result["output_smiles"])
    assert mol is not None
