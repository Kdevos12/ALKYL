# tests/test_chem_compare.py
from tests.conftest import run_script, ASPIRIN_SMILES, CAFFEINE_SMILES

SCRIPT = "chem_compare.py"
PHENYL_ACETATE = "CC(=O)Oc1ccccc1"


def test_self_comparison():
    result = run_script(SCRIPT, ["--smiles-a", ASPIRIN_SMILES,
                                 "--smiles-b", ASPIRIN_SMILES])
    assert result["tanimoto"] == 1.0
    for k, v in result["delta_b_minus_a"].items():
        assert v == 0, f"delta[{k}] should be 0 for self-comparison"


def test_different_molecules():
    result = run_script(SCRIPT, ["--smiles-a", ASPIRIN_SMILES,
                                 "--smiles-b", CAFFEINE_SMILES])
    assert result["tanimoto"] < 1.0
    assert result["tanimoto"] >= 0.0


def test_mcs_fields():
    result = run_script(SCRIPT, ["--smiles-a", ASPIRIN_SMILES,
                                 "--smiles-b", PHENYL_ACETATE])
    mcs = result["mcs"]
    assert "smarts" in mcs
    assert "n_atoms" in mcs
    assert "n_bonds" in mcs
    assert mcs["n_atoms"] > 0  # phenyl acetate is a subset of aspirin


def test_delta_sign():
    # Caffeine (MW ~194) vs Aspirin (MW ~180): delta should be positive
    result = run_script(SCRIPT, ["--smiles-a", ASPIRIN_SMILES,
                                 "--smiles-b", CAFFEINE_SMILES])
    # delta = B - A; caffeine heavier than aspirin
    assert result["delta_b_minus_a"]["mw"] > 0


def test_properties_present():
    result = run_script(SCRIPT, ["--smiles-a", ASPIRIN_SMILES,
                                 "--smiles-b", CAFFEINE_SMILES])
    for section in ["properties_a", "properties_b", "delta_b_minus_a"]:
        for key in ["mw", "logp", "tpsa", "hbd", "hba", "rotbonds", "heavy_atoms"]:
            assert key in result[section]


def test_maccs_fingerprint():
    result = run_script(SCRIPT, ["--smiles-a", ASPIRIN_SMILES,
                                 "--smiles-b", ASPIRIN_SMILES,
                                 "--fingerprint", "maccs"])
    assert result["tanimoto"] == 1.0
    assert result["fingerprint"] == "maccs"
