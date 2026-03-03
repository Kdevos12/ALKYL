#!/usr/bin/env python3
"""ALKYL — chem_react.py: apply a reaction SMARTS to a molecule and return products."""

import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Apply a reaction SMARTS to a molecule.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Reactant SMILES")
    inp.add_argument("--sdf", help="Reactant SDF file")
    parser.add_argument("--reaction", required=True,
                        help="Reaction SMARTS (e.g. '[C:1](=O)[OH]>>[C:1](=O)Cl')")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import AllChem

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

    rxn = AllChem.ReactionFromSmarts(args.reaction)
    if rxn is None:
        print(f"Invalid reaction SMARTS: {args.reaction}", file=sys.stderr)
        sys.exit(1)

    product_sets = rxn.RunReactants((mol,))
    products = []
    seen = set()
    for prod_tuple in product_sets:
        for prod_mol in prod_tuple:
            try:
                Chem.SanitizeMol(prod_mol)
                smi = Chem.MolToSmiles(prod_mol)
                if smi not in seen:
                    seen.add(smi)
                    products.append(smi)
            except Exception:
                pass

    result = {
        "reactant":   Chem.MolToSmiles(mol),
        "reaction":   args.reaction,
        "n_products": len(products),
        "products":   products,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
