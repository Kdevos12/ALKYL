# tests/test_chem_analyze.py
from tests.conftest import run_script, ASPIRIN_SMILES, CAFFEINE_SMILES

SCRIPT = "chem_analyze.py"

# Ibuprofen: one stereocenter
IBUPROFEN_SMILES = "CC(Cc1ccc(CC(C)C(=O)O)cc1)C(=O)O"
# (R)-Alanine: one assigned stereocenter
ALANINE_SMILES = "N[C@@H](C)C(=O)O"
# Simple amine: primary amine
METHYLAMINE_SMILES = "CN"


def test_output_structure():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    required = [
        "smiles", "formula", "heavy_atoms", "charge",
        "mw", "exact_mw", "hbd", "hba", "rotbonds",
        "tpsa", "logp", "fsp3",
        "rings", "stereocenters", "functional_groups",
        "qed", "sa_score", "complexity_bertz",
    ]
    for field in required:
        assert field in result, f"Missing field: {field}"


def test_aspirin_formula_and_mw():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["formula"] == "C9H8O4"
    assert abs(result["mw"] - 180.16) < 0.1
    assert result["heavy_atoms"] == 13


def test_aspirin_functional_groups():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    fg = result["functional_groups"]
    assert fg.get("carboxylic_acid", 0) == 1
    assert fg.get("ester", 0) == 1


def test_aspirin_rings():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    rings = result["rings"]
    assert rings["aromatic"] == 1
    assert rings["aliphatic"] == 0
    assert rings["total"] == 1


def test_aspirin_no_stereocenters():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["stereocenters"]["total"] == 0


def test_ibuprofen_stereocenter():
    result = run_script(SCRIPT, ["--smiles", IBUPROFEN_SMILES])
    # Ibuprofen has 2 stereocenters (both unassigned without @)
    assert result["stereocenters"]["total"] >= 1


def test_alanine_assigned_stereocenter():
    result = run_script(SCRIPT, ["--smiles", ALANINE_SMILES])
    assert result["stereocenters"]["assigned"] == 1
    assert result["stereocenters"]["unassigned"] == 0


def test_caffeine_rings():
    result = run_script(SCRIPT, ["--smiles", CAFFEINE_SMILES])
    # Caffeine has fused bicyclic system (pyrimidine + imidazole)
    assert result["rings"]["aromatic"] == 2
    assert result["rings"]["total"] == 2


def test_qed_range():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert 0.0 <= result["qed"] <= 1.0


def test_sa_score_range():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert 1.0 <= result["sa_score"] <= 10.0


def test_aspirin_sa_easy():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    # Aspirin is trivially synthesizable
    assert result["sa_score"] < 3.0


def test_neutral_charge():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["charge"] == 0


def test_charged_molecule():
    result = run_script(SCRIPT, ["--smiles", "CC(=O)[O-]"])
    assert result["charge"] == -1
