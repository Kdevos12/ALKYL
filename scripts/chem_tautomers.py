#!/usr/bin/env python3
"""ALKYL — chem_tautomers.py: enumerate tautomers and find the canonical form."""

import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Enumerate tautomers and identify the canonical form.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--max-tautomers", type=int, default=100,
                        dest="max_tautomers",
                        help="Maximum tautomers to enumerate (default: 100)")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

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

    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(args.max_tautomers)

    canonical_smi = Chem.MolToSmiles(enumerator.Canonicalize(mol))
    tautomers = sorted({Chem.MolToSmiles(t) for t in enumerator.Enumerate(mol)})

    result = {
        "input_smiles":       Chem.MolToSmiles(mol),
        "canonical_tautomer": canonical_smi,
        "n_tautomers":        len(tautomers),
        "tautomers":          tautomers,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
