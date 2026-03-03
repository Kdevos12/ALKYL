# tests/test_chem_batch.py
import json
import subprocess
import sys
from tests.conftest import ASPIRIN_SMILES, CAFFEINE_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_batch.py"


def run_batch(args) -> dict:
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT)] + args,
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, f"Script failed:\n{result.stderr}"
    return json.loads(result.stdout)


def make_smi(tmp_path, lines: list[str]) -> str:
    p = tmp_path / "library.smi"
    p.write_text("\n".join(lines))
    return str(p)


def make_csv(tmp_path, rows: list[tuple[str, str]]) -> str:
    p = tmp_path / "library.csv"
    lines = ["smiles,name"] + [f"{s},{n}" for s, n in rows]
    p.write_text("\n".join(lines))
    return str(p)


# ── descriptors ──────────────────────────────────────────────────────────────

def test_batch_descriptors_basic(tmp_path):
    path = make_smi(tmp_path, [
        f"{ASPIRIN_SMILES} aspirin",
        f"{CAFFEINE_SMILES} caffeine",
    ])
    result = run_batch(["--input", path, "--format", "smi",
                        "--descriptors", "mw,logp"])
    assert result["n_total"] == 2
    assert result["n_processed"] == 2
    by_name = {m["name"]: m for m in result["molecules"]}
    assert abs(by_name["aspirin"]["mw"] - 180.16) < 0.1
    assert "logp" in by_name["caffeine"]


def test_batch_all_descriptors(tmp_path):
    path = make_smi(tmp_path, [f"{ASPIRIN_SMILES} aspirin"])
    result = run_batch(["--input", path, "--format", "smi", "--descriptors", "all"])
    mol = result["molecules"][0]
    for key in ["mw", "logp", "hbd", "hba", "tpsa", "rotbonds", "rings", "fsp3"]:
        assert key in mol, f"Missing descriptor: {key}"


# ── Lipinski / PAINS ─────────────────────────────────────────────────────────

def test_batch_lipinski(tmp_path):
    path = make_smi(tmp_path, [f"{ASPIRIN_SMILES} aspirin"])
    result = run_batch(["--input", path, "--format", "smi", "--lipinski"])
    assert result["molecules"][0]["lipinski"]["pass"] is True
    assert result["molecules"][0]["lipinski"]["violations"] == []


def test_batch_pains(tmp_path):
    path = make_smi(tmp_path, [f"{ASPIRIN_SMILES} aspirin"])
    result = run_batch(["--input", path, "--format", "smi", "--pains"])
    pains = result["molecules"][0]["pains"]
    assert "alerts" in pains
    assert "clean" in pains


# ── invalid molecules ─────────────────────────────────────────────────────────

def test_batch_invalid_mol_included(tmp_path):
    path = make_smi(tmp_path, [
        f"{ASPIRIN_SMILES} aspirin",
        "NOTASMILES bad_mol",
    ])
    result = run_batch(["--input", path, "--format", "smi", "--descriptors", "mw"])
    assert result["n_total"] == 2
    assert result["n_processed"] == 2
    assert result["n_errors"] == 1
    errors = [m for m in result["molecules"] if m.get("error")]
    assert len(errors) == 1
    assert errors[0]["name"] == "bad_mol"


def test_batch_skip_invalid(tmp_path):
    path = make_smi(tmp_path, [
        f"{ASPIRIN_SMILES} aspirin",
        "NOTASMILES bad_mol",
    ])
    result = run_batch(["--input", path, "--format", "smi",
                        "--descriptors", "mw", "--skip-invalid"])
    assert result["n_processed"] == 1
    assert all(m.get("error") is None for m in result["molecules"])


# ── CSV format ────────────────────────────────────────────────────────────────

def test_batch_csv(tmp_path):
    path = make_csv(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_batch(["--input", path, "--format", "csv", "--lipinski"])
    assert result["n_total"] == 2
    by_name = {m["name"]: m for m in result["molecules"]}
    assert by_name["aspirin"]["lipinski"]["pass"] is True


# ── --out file ────────────────────────────────────────────────────────────────

def test_batch_output_file(tmp_path):
    path = make_smi(tmp_path, [f"{ASPIRIN_SMILES} aspirin"])
    out = tmp_path / "results.json"
    summary = run_batch(["--input", path, "--format", "smi",
                         "--descriptors", "mw", "--out", str(out)])
    # stdout = summary (no molecules key)
    assert "n_total" in summary
    assert "molecules" not in summary
    assert summary["out"] == str(out)
    # file = full data
    with open(str(out)) as f:
        full = json.load(f)
    assert "molecules" in full
    assert full["molecules"][0]["name"] == "aspirin"
