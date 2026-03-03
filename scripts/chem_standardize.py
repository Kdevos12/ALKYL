#!/usr/bin/env python3
"""ALKYL — chem_standardize.py: desalt, neutralize, and canonicalize molecules."""

import argparse
import json
import sys


def load_mol(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"Invalid SMILES: {args.smiles}", file=sys.stderr)
            sys.exit(1)
        return mol
    if args.sdf:
        suppl = Chem.SDMolSupplier(args.sdf)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            print(f"No valid molecule in: {args.sdf}", file=sys.stderr)
            sys.exit(1)
        return mol
    print("Provide --smiles or --sdf.", file=sys.stderr)
    sys.exit(1)


def standardize(mol) -> tuple[object, list[str], str]:
    """Run desalt → neutralize → canonicalize pipeline.

    Returns (standardized_mol, list_of_changes, canonical_smiles).
    """
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    changes = []

    # Step 1 — desalt: keep largest organic fragment
    before = Chem.MolToSmiles(mol)
    chooser = rdMolStandardize.LargestFragmentChooser()
    mol = chooser.choose(mol)
    if Chem.MolToSmiles(mol) != before:
        changes.append("desalted")

    # Step 2 — neutralize: remove formal charges where possible
    before = Chem.MolToSmiles(mol)
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    if Chem.MolToSmiles(mol) != before:
        changes.append("neutralized")

    # Step 3 — canonical SMILES (always)
    canonical = Chem.MolToSmiles(mol)

    return mol, changes, canonical


def main():
    parser = argparse.ArgumentParser(
        description="Standardize a molecule: desalt, neutralize, canonicalize.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem

    mol = load_mol(args)
    input_smi = Chem.MolToSmiles(mol)

    _mol_std, changes, output_smi = standardize(mol)

    result = {
        "input_smiles": input_smi,
        "output_smiles": output_smi,
        "changes": changes,
        "changed": len(changes) > 0,
        "valid": True,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
