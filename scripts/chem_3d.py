#!/usr/bin/env python3
"""ALKYL — chem_3d.py: 3D conformer generation, minimization, RMSD."""

import argparse
import json
import sys


def load_mol_2d(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print("Invalid SMILES.", file=sys.stderr)
            sys.exit(1)
        mol = Chem.AddHs(mol)
        return mol
    if args.sdf:
        suppl = Chem.SDMolSupplier(args.sdf, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            print("No valid molecule in SDF.", file=sys.stderr)
            sys.exit(1)
        return mol
    print("Provide --smiles or --sdf.", file=sys.stderr)
    sys.exit(1)


def generate_conformers(mol, n: int, method: str) -> list[int]:
    from rdkit.Chem import AllChem
    if method == "etkdg":
        params = AllChem.ETKDGv3()
    else:
        params = AllChem.ETDG()
    params.randomSeed = 42
    params.numThreads = 0
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n, params=params)
    return list(ids)


def minimize_conformers(mol, ff: str) -> list[float | None]:
    from rdkit.Chem import AllChem
    energies = []
    for conf_id in range(mol.GetNumConformers()):
        if ff == "mmff94":
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is None:
                energies.append(None)
                continue
            ff_obj = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
        else:
            ff_obj = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        if ff_obj is None:
            energies.append(None)
            continue
        ff_obj.Minimize()
        energies.append(round(ff_obj.CalcEnergy(), 4))
    return energies


def main():
    parser = argparse.ArgumentParser(description="3D conformer generation.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles")
    inp.add_argument("--sdf")
    parser.add_argument("--conformers", type=int, default=10,
                        help="Number of conformers to generate")
    parser.add_argument("--method", choices=["etkdg", "etdg"], default="etkdg")
    parser.add_argument("--minimize", choices=["mmff94", "uff"],
                        help="Force field for minimization")
    parser.add_argument("--out", help="Write multi-conf SDF to file")
    args = parser.parse_args()

    mol = load_mol_2d(args)
    conf_ids = generate_conformers(mol, args.conformers, args.method)

    energies = None
    if args.minimize:
        energies = minimize_conformers(mol, args.minimize)

    if args.out:
        from rdkit.Chem import SDWriter
        writer = SDWriter(args.out)
        for cid in conf_ids:
            writer.write(mol, confId=cid)
        writer.close()
        # Still print JSON metadata to stdout for Claude to parse
        result = {"n_conformers": len(conf_ids), "out": args.out,
                  "energies_kcal": energies}
    else:
        result = {
            "n_conformers": len(conf_ids),
            "method": args.method,
            "energies_kcal": energies,
        }

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
