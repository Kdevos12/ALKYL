# tests/test_chem_metabolism.py
from tests.conftest import run_script, ASPIRIN_SMILES, CAFFEINE_SMILES

SCRIPT = "chem_metabolism.py"


def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for f in ["smiles", "metabolic_sites", "n_sites", "overall_risk"]:
        assert f in result


def test_aspirin_has_aromatic_site():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    patterns = [s["pattern"] for s in result["metabolic_sites"]]
    assert "aromatic_hydroxylation" in patterns


def test_caffeine_n_dealkylation():
    # Caffeine has N-CH3 groups → N-dealkylation expected
    result = run_script(SCRIPT, ["--smiles", CAFFEINE_SMILES])
    patterns = [s["pattern"] for s in result["metabolic_sites"]]
    assert "n_dealkylation" in patterns


def test_no_sites_simple():
    # Ethane: no metabolic alerts
    result = run_script(SCRIPT, ["--smiles", "CC"])
    assert result["n_sites"] == 0
    assert result["overall_risk"] == "none"


def test_site_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for site in result["metabolic_sites"]:
        for f in ["enzyme", "pattern", "atom_indices", "risk"]:
            assert f in site


def test_overall_risk_is_valid():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["overall_risk"] in {"none", "low", "medium", "high"}


def test_n_sites_matches_list():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    assert result["n_sites"] == len(result["metabolic_sites"])
