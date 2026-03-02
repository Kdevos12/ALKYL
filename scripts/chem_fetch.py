#!/usr/bin/env python3
"""ALKYL — chem_fetch.py: fetch molecular data from PubChem or ChEMBL."""

import argparse
import json
import sys
import urllib.request
import urllib.parse


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


def _get_json(url: str) -> dict:
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        print(f"HTTP error fetching {url}: {e}", file=sys.stderr)
        sys.exit(1)


def fetch_pubchem(args) -> dict:
    if args.cid:
        namespace, ident = "cid", str(args.cid)
    elif args.name:
        namespace, ident = "name", urllib.parse.quote(args.name)
    elif args.smiles:
        namespace, ident = "smiles", urllib.parse.quote(args.smiles)
    else:
        print("Provide --name, --cid, or --smiles.", file=sys.stderr)
        sys.exit(1)

    props = "CanonicalSMILES,InChI,InChIKey,IUPACName,MolecularWeight,XLogP"
    url = f"{PUBCHEM_BASE}/compound/{namespace}/{ident}/property/{props}/JSON"
    data = _get_json(url)
    props_data = data["PropertyTable"]["Properties"][0]

    result = {
        "cid": props_data.get("CID"),
        "smiles": props_data.get("CanonicalSMILES"),
        "inchi": props_data.get("InChI"),
        "inchikey": props_data.get("InChIKey"),
        "iupac_name": props_data.get("IUPACName"),
        "mw": props_data.get("MolecularWeight"),
        "logp": props_data.get("XLogP"),
    }
    if args.properties:
        keep = {p.strip() for p in args.properties.split(",")}
        # Always keep cid for identification
        result = {k: v for k, v in result.items() if k in keep or k == "cid"}
    return result


def fetch_chembl(args) -> dict:
    if args.chembl_id:
        url = f"{CHEMBL_BASE}/molecule/{args.chembl_id}?format=json"
    elif args.name:
        url = f"{CHEMBL_BASE}/molecule?pref_name={urllib.parse.quote(args.name)}&format=json"
    else:
        print("Provide --chembl-id or --name.", file=sys.stderr)
        sys.exit(1)

    data = _get_json(url)
    if "molecules" in data:
        if not data["molecules"]:
            print("No results found.", file=sys.stderr)
            sys.exit(1)
        mol = data["molecules"][0]
    else:
        mol = data

    struct = mol.get("molecule_structures") or {}
    return {
        "chembl_id": mol.get("molecule_chembl_id"),
        "smiles": struct.get("canonical_smiles"),
        "inchi": struct.get("standard_inchi"),
        "inchikey": struct.get("standard_inchi_key"),
        "name": mol.get("pref_name"),
        "type": mol.get("molecule_type"),
        "mw": mol.get("molecule_properties", {}).get("full_mwt"),
        "logp": mol.get("molecule_properties", {}).get("alogp"),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Fetch molecular data from PubChem or ChEMBL.")
    parser.add_argument("--source", required=True, choices=["pubchem", "chembl"])
    parser.add_argument("--name", help="Molecule name")
    parser.add_argument("--smiles", help="SMILES (PubChem only)")
    parser.add_argument("--cid", type=int, help="PubChem CID")
    parser.add_argument("--chembl-id", dest="chembl_id", help="ChEMBL ID")
    parser.add_argument("--properties", help="Comma-separated property subset")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    result = fetch_pubchem(args) if args.source == "pubchem" else fetch_chembl(args)

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
