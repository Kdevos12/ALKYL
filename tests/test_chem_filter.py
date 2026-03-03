# tests/test_chem_filter.py
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_filter.py"
# Cyclosporin: known Lipinski violator (MW~1202)
BIG_MOL = "CC1NC(=O)C(CC2=CC=CC=C2)N(C)C(=O)C(CC(C)C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)NC(=O)C(C(C)CC)N(C)C1=O"


def test_aspirin_passes_lipinski():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["rules"]["lipinski"]["pass"] is True
    assert result["rules"]["lipinski"]["violations"] == []


def test_aspirin_passes_veber():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["rules"]["veber"]["pass"] is True


def test_aspirin_passes_pains():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["rules"]["pains"]["pass"] is True


def test_aspirin_fails_ghose():
    # Aspirin has 13 heavy atoms < 20 required by Ghose → expected failure
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["rules"]["ghose"]["pass"] is False
    assert result["n_rules_pass"] == 4  # lipinski, veber, egan, pains pass


def test_big_mol_fails_lipinski():
    result = run_script(SCRIPT, ["--smiles", BIG_MOL])
    assert result["rules"]["lipinski"]["pass"] is False
    assert len(result["rules"]["lipinski"]["violations"]) > 0


def test_qed_and_sa_range():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert 0.0 <= result["qed"] <= 1.0
    assert 1.0 <= result["sa_score"] <= 10.0


def test_output_structure():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for f in ["smiles", "descriptors", "qed", "sa_score",
              "rules", "n_rules_pass", "n_rules_total", "overall_pass"]:
        assert f in result
