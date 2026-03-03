#!/usr/bin/env python3
"""ALKYL — chem_enum.py: enumerate stereoisomers of a molecule."""

import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Enumerate stereoisomers of a molecule.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--max-isomers", type=int, default=64,
                        dest="max_isomers",
                        help="Maximum stereoisomers to return (default: 64)")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem.EnumerateStereoisomers import (
        EnumerateStereoisomers, StereoEnumerationOptions)
    from rdkit.Chem import FindMolChiralCenters

    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"Invalid SMILES: {args.smiles}", file=sys.stderr)
            sys.exit(1)
    else:
        suppl = Chem.SDMolSupplier(args.sdf)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            print(f"No valid molecule in: {args.sdf}", file=sys.stderr)
            sys.exit(1)

    opts = StereoEnumerationOptions(unique=True, maxIsomers=args.max_isomers)
    isomers = [Chem.MolToSmiles(iso)
               for iso in EnumerateStereoisomers(mol, options=opts)]

    centers = FindMolChiralCenters(mol, includeUnassigned=True)

    result = {
        "input_smiles":            Chem.MolToSmiles(mol),
        "stereocenters_in_input":  len(centers),
        "n_isomers":               len(isomers),
        "isomers":                 isomers,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
