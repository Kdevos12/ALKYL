#!/usr/bin/env python3
"""ALKYL — chem_scaffold.py: Murcko scaffold and BRICS fragment decomposition."""

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
    suppl = Chem.SDMolSupplier(args.sdf)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        print(f"No valid molecule in: {args.sdf}", file=sys.stderr)
        sys.exit(1)
    return mol


def main():
    parser = argparse.ArgumentParser(
        description="Scaffold decomposition: Murcko scaffold and BRICS fragments.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem import BRICS

    mol = load_mol(args)

    # Murcko scaffold (ring systems + linkers, side chains removed)
    murcko_mol = MurckoScaffold.GetScaffoldForMol(mol)
    murcko_smi = Chem.MolToSmiles(murcko_mol) if murcko_mol.GetNumAtoms() > 0 else ""

    # Generic scaffold (all atoms → C, all bonds → single)
    generic_smi = ""
    if murcko_mol.GetNumAtoms() > 0:
        generic_mol = MurckoScaffold.MakeScaffoldGeneric(murcko_mol)
        generic_smi = Chem.MolToSmiles(generic_mol)

    # BRICS fragments (retrosynthetic disconnection rules)
    brics_frags = sorted(BRICS.BRICSDecompose(mol))

    result = {
        "smiles": Chem.MolToSmiles(mol),
        "murcko_scaffold": murcko_smi,
        "generic_scaffold": generic_smi,
        "has_rings": murcko_mol.GetNumAtoms() > 0,
        "n_ring_systems": murcko_mol.GetRingInfo().NumRings(),
        "brics_fragments": brics_frags,
        "n_brics_fragments": len(brics_frags),
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
