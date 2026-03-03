# tests/test_chem_scaffold.py
from tests.conftest import run_script, ASPIRIN_SMILES, CAFFEINE_SMILES

SCRIPT = "chem_scaffold.py"
IBUPROFEN = "CC(Cc1ccc(CC(C)C(=O)O)cc1)C(=O)O"
ACYCLIC   = "CCCCCC"  # no rings


def test_aspirin_murcko():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["has_rings"] is True
    assert result["murcko_scaffold"] != ""
    assert result["n_ring_systems"] == 1


def test_aspirin_generic_scaffold():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    # Generic scaffold: all C, no heteroatoms
    assert "c" not in result["generic_scaffold"]


def test_aspirin_brics():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["n_brics_fragments"] >= 1
    assert isinstance(result["brics_fragments"], list)


def test_caffeine_two_rings():
    result = run_script(SCRIPT, ["--smiles", CAFFEINE_SMILES])
    assert result["n_ring_systems"] == 2


def test_acyclic_no_scaffold():
    result = run_script(SCRIPT, ["--smiles", ACYCLIC])
    assert result["has_rings"] is False
    assert result["murcko_scaffold"] == ""
    assert result["generic_scaffold"] == ""


def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for f in ["smiles", "murcko_scaffold", "generic_scaffold",
              "has_rings", "n_ring_systems", "brics_fragments", "n_brics_fragments"]:
        assert f in result
