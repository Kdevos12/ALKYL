# tests/test_chem_search.py
import json
import subprocess
import sys
from tests.conftest import ASPIRIN_SMILES, CAFFEINE_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_search.py"

# Aspirin fragment: phenyl acetate (shares benzene + ester with aspirin)
PHENYL_ACETATE_SMILES = "CC(=O)Oc1ccccc1"
BENZENE_SMILES = "c1ccccc1"


def run_search(args) -> dict:
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT)] + args,
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, f"Script failed:\n{result.stderr}"
    return json.loads(result.stdout)


def make_library(tmp_path, mols: list[tuple[str, str]]) -> str:
    p = tmp_path / "library.smi"
    p.write_text("\n".join(f"{smi} {name}" for smi, name in mols))
    return str(p)


# ── substructure ──────────────────────────────────────────────────────────────

def test_substructure_hit(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    # Benzene ring present in aspirin, not in caffeine (purine scaffold)
    result = run_search(["--query", BENZENE_SMILES,
                         "--library", lib, "--mode", "substructure"])
    assert result["n_hits"] == 1
    assert result["hits"][0]["name"] == "aspirin"


def test_substructure_no_hits(tmp_path):
    lib = make_library(tmp_path, [(CAFFEINE_SMILES, "caffeine")])
    result = run_search(["--query", BENZENE_SMILES,
                         "--library", lib, "--mode", "substructure"])
    assert result["n_hits"] == 0
    assert result["hits"] == []


def test_substructure_smarts(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    # Ester/carboxylic acid SMARTS — aspirin has both groups
    result = run_search(["--smarts", "C(=O)O",
                         "--library", lib, "--mode", "substructure"])
    names = [h["name"] for h in result["hits"]]
    assert "aspirin" in names


def test_substructure_multiple_hits(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (PHENYL_ACETATE_SMILES, "phenyl_acetate"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_search(["--query", BENZENE_SMILES,
                         "--library", lib, "--mode", "substructure"])
    names = {h["name"] for h in result["hits"]}
    assert names == {"aspirin", "phenyl_acetate"}


# ── exact match ───────────────────────────────────────────────────────────────

def test_exact_match(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_search(["--query", ASPIRIN_SMILES,
                         "--library", lib, "--mode", "exact"])
    assert result["n_hits"] == 1
    assert result["hits"][0]["name"] == "aspirin"


def test_exact_no_match(tmp_path):
    lib = make_library(tmp_path, [(CAFFEINE_SMILES, "caffeine")])
    result = run_search(["--query", ASPIRIN_SMILES,
                         "--library", lib, "--mode", "exact"])
    assert result["n_hits"] == 0


# ── similarity ────────────────────────────────────────────────────────────────

def test_similarity_self(tmp_path):
    lib = make_library(tmp_path, [(ASPIRIN_SMILES, "aspirin")])
    result = run_search(["--query", ASPIRIN_SMILES, "--library", lib,
                         "--mode", "similarity", "--threshold", "0.99"])
    assert result["n_hits"] == 1
    assert result["hits"][0]["tanimoto"] > 0.99


def test_similarity_sorted(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (PHENYL_ACETATE_SMILES, "phenyl_acetate"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_search(["--query", ASPIRIN_SMILES, "--library", lib,
                         "--mode", "similarity", "--threshold", "0.0"])
    scores = [h["tanimoto"] for h in result["hits"]]
    assert scores == sorted(scores, reverse=True)


def test_similarity_threshold_filters(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    # Very high threshold: only aspirin (Tanimoto ≈ 1.0 against itself)
    result = run_search(["--query", ASPIRIN_SMILES, "--library", lib,
                         "--mode", "similarity", "--threshold", "0.99"])
    assert all(h["tanimoto"] >= 0.99 for h in result["hits"])


def test_similarity_top_k(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (PHENYL_ACETATE_SMILES, "phenyl_acetate"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_search(["--query", ASPIRIN_SMILES, "--library", lib,
                         "--mode", "similarity", "--threshold", "0.0",
                         "--top-k", "2"])
    assert len(result["hits"]) <= 2


def test_similarity_maccs(tmp_path):
    lib = make_library(tmp_path, [(ASPIRIN_SMILES, "aspirin")])
    result = run_search(["--query", ASPIRIN_SMILES, "--library", lib,
                         "--mode", "similarity", "--fingerprint", "maccs",
                         "--threshold", "0.99"])
    assert result["n_hits"] == 1


# ── output ────────────────────────────────────────────────────────────────────

def test_output_structure(tmp_path):
    lib = make_library(tmp_path, [(ASPIRIN_SMILES, "aspirin")])
    result = run_search(["--query", BENZENE_SMILES,
                         "--library", lib, "--mode", "substructure"])
    assert "mode" in result
    assert "n_library" in result
    assert "n_hits" in result
    assert "hits" in result
    assert result["hits"][0]["smiles"]  # canonical SMILES present


def test_output_file(tmp_path):
    lib = make_library(tmp_path, [(ASPIRIN_SMILES, "aspirin")])
    out = tmp_path / "hits.json"
    subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--query", BENZENE_SMILES, "--library", lib,
         "--mode", "substructure", "--out", str(out)],
        check=True, timeout=SCRIPT_TIMEOUT,
    )
    with open(str(out)) as f:
        data = json.load(f)
    assert data["n_hits"] == 1
