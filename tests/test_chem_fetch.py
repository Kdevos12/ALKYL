# tests/test_chem_fetch.py
# Network tests — require internet access to PubChem/ChEMBL REST APIs.
# Run with: pytest -m network
import pytest
from tests.conftest import run_script

SCRIPT = "chem_fetch.py"

@pytest.mark.network
def test_fetch_by_name_pubchem():
    result = run_script(SCRIPT, ["--source", "pubchem",
                                 "--name", "aspirin",
                                 "--properties", "cid,smiles,inchi"])
    assert result.get("cid") == 2244
    assert "smiles" in result
    assert "inchi" in result

@pytest.mark.network
def test_fetch_by_cid():
    result = run_script(SCRIPT, ["--source", "pubchem",
                                 "--cid", "2244"])
    assert "smiles" in result
    assert "iupac_name" in result

@pytest.mark.network
def test_fetch_chembl_by_id():
    result = run_script(SCRIPT, ["--source", "chembl",
                                 "--chembl-id", "CHEMBL25"])
    assert "smiles" in result
    assert result.get("chembl_id") == "CHEMBL25"
