# tests/test_chem_qm.py
import json, subprocess, sys
from pathlib import Path
from tests.conftest import run_script, SCRIPTS_DIR, ASPIRIN_SMILES

SCRIPT = "chem_qm.py"

def test_generate_orca_input(tmp_path):
    out_file = tmp_path / "aspirin.inp"
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--smiles", ASPIRIN_SMILES,
         "--engine", "orca", "--task", "opt",
         "--method", "B3LYP", "--basis", "6-31G*",
         "--charge", "0", "--mult", "1",
         "--out", str(out_file)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0
    content = out_file.read_text()
    assert "B3LYP" in content
    assert "6-31G*" in content
    assert "! Opt" in content

def test_generate_gaussian_input(tmp_path):
    out_file = tmp_path / "aspirin.gjf"
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--smiles", ASPIRIN_SMILES,
         "--engine", "gaussian", "--task", "sp",
         "--method", "B3LYP", "--basis", "6-31G*",
         "--out", str(out_file)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0
    content = out_file.read_text()
    assert "#P B3LYP" in content
    assert "6-31G*" in content

def test_stdout_json_when_no_out():
    """When --out is not given, print JSON metadata to stdout."""
    result = run_script(SCRIPT, [
        "--smiles", ASPIRIN_SMILES,
        "--engine", "orca", "--task", "sp",
        "--method", "HF", "--basis", "STO-3G",
    ])
    assert "engine" in result
    assert result["engine"] == "orca"
