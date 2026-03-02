#!/usr/bin/env python3
"""ALKYL — chem_qm.py: QM input generation and output parsing (ORCA, Gaussian)."""

import argparse
import json
import sys
from pathlib import Path


def get_xyz_block(mol) -> str:
    """Return XYZ coordinate block from a 3D RDKit mol."""
    conf = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"  {atom.GetSymbol():2s}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}"
        )
    return "\n".join(lines)


def smiles_to_3d(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES.", file=sys.stderr)
        sys.exit(1)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        print("3D embedding failed for this SMILES.", file=sys.stderr)
        sys.exit(1)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


TASK_ORCA = {"sp": "SP", "opt": "Opt", "freq": "Freq", "scan": "Scan"}
TASK_GAUSSIAN = {"sp": "SP", "opt": "Opt", "freq": "Freq", "scan": "Scan"}


def write_orca(mol, method: str, basis: str, task: str,
               charge: int, mult: int) -> str:
    xyz = get_xyz_block(mol)
    keyword = TASK_ORCA.get(task, task)
    return (
        f"! {method} {basis} {keyword}\n"
        f"\n"
        f"* xyz {charge} {mult}\n"
        f"{xyz}\n"
        f"*\n"
    )


def write_gaussian(mol, method: str, basis: str, task: str,
                   charge: int, mult: int) -> str:
    xyz = get_xyz_block(mol)
    keyword = TASK_GAUSSIAN.get(task, task)
    return (
        f"%mem=4GB\n"
        f"%nprocshared=4\n"
        f"#P {method}/{basis} {keyword}\n"
        f"\n"
        f"ALKYL generated input\n"
        f"\n"
        f"{charge} {mult}\n"
        f"{xyz}\n"
        f"\n"
    )


def parse_orca_output(path: str) -> dict:
    """Parse key results from an ORCA output file."""
    result = {}
    text = Path(path).read_text(errors="replace")
    for line in reversed(text.splitlines()):
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                result["energy_hartree"] = float(line.split()[-1])
            except ValueError:
                pass
            break
    freqs = []
    for line in text.splitlines():
        if ":" in line and "cm**-1" in line:
            try:
                freqs.append(float(line.split(":")[1].split()[0]))
            except (ValueError, IndexError):
                pass
    if freqs:
        result["frequencies_cm1"] = freqs
    if not result:
        result["warning"] = "No recognized patterns found in output file"
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate QM input files or parse QM outputs.")
    inp = parser.add_mutually_exclusive_group()
    inp.add_argument("--smiles", help="Input SMILES (auto 3D via ETKDGv3+MMFF)")
    inp.add_argument("--xyz", help="Input XYZ file (pre-computed geometry)")
    parser.add_argument("--engine", choices=["orca", "gaussian"], default="orca")
    parser.add_argument("--task", choices=["sp", "opt", "freq", "scan"],
                        default="sp")
    parser.add_argument("--method", default="B3LYP")
    parser.add_argument("--basis", default="6-31G*")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--mult", type=int, default=1)
    parser.add_argument("--out", help="Write input file to this path")
    parser.add_argument("--parse", help="Parse existing QM output file -> JSON")
    args = parser.parse_args()

    # Parse mode
    if args.parse:
        if args.engine == "orca":
            result = parse_orca_output(args.parse)
        else:
            print("Gaussian parsing not yet implemented.", file=sys.stderr)
            sys.exit(1)
        print(json.dumps(result, indent=2))
        return

    # Generation mode — require at least one input
    if not args.smiles and not args.xyz:
        parser.error("Provide --smiles or --xyz for input generation.")

    if args.smiles:
        mol = smiles_to_3d(args.smiles)
    else:
        print("--xyz input requires ASE. Use --smiles instead.", file=sys.stderr)
        sys.exit(1)

    if args.engine == "orca":
        content = write_orca(mol, args.method, args.basis,
                             args.task, args.charge, args.mult)
    else:
        content = write_gaussian(mol, args.method, args.basis,
                                 args.task, args.charge, args.mult)

    meta = {"engine": args.engine, "task": args.task,
            "method": args.method, "basis": args.basis}

    if args.out:
        Path(args.out).write_text(content)
        meta["written"] = args.out

    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
