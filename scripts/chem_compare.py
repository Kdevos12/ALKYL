#!/usr/bin/env python3
"""ALKYL — chem_compare.py: compare two molecules by MCS, Tanimoto, and property delta."""

import argparse
import json
import sys


def _load_mol(smiles: str, label: str):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES for {label}: {smiles}", file=sys.stderr)
        sys.exit(1)
    return mol


def _fingerprint(mol, fp_type: str, radius: int, nbits: int):
    if fp_type == "morgan":
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
        return gen.GetFingerprint(mol)
    if fp_type == "maccs":
        from rdkit.Chem import MACCSkeys
        return MACCSkeys.GenMACCSKeys(mol)
    from rdkit.Chem.rdmolops import RDKFingerprint
    return RDKFingerprint(mol, fpSize=nbits)


def main():
    parser = argparse.ArgumentParser(
        description="Compare two molecules: MCS, Tanimoto similarity, property delta.")
    parser.add_argument("--smiles-a", required=True, dest="smiles_a",
                        help="First molecule SMILES")
    parser.add_argument("--smiles-b", required=True, dest="smiles_b",
                        help="Second molecule SMILES")
    parser.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"],
                        default="morgan",
                        help="Fingerprint type for Tanimoto (default: morgan)")
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--nbits", type=int, default=2048)
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFMCS, Descriptors

    mol_a = _load_mol(args.smiles_a, "A")
    mol_b = _load_mol(args.smiles_b, "B")

    # Maximum Common Substructure
    mcs = rdFMCS.FindMCS([mol_a, mol_b], timeout=5,
                         ringMatchesRingOnly=True, completeRingsOnly=True)

    # Tanimoto similarity
    fp_a = _fingerprint(mol_a, args.fingerprint, args.radius, args.nbits)
    fp_b = _fingerprint(mol_b, args.fingerprint, args.radius, args.nbits)
    tanimoto = round(DataStructs.TanimotoSimilarity(fp_a, fp_b), 4)

    # Property delta (B − A)
    fns = {
        "mw":        (Descriptors.MolWt,              True),
        "logp":      (Descriptors.MolLogP,            True),
        "tpsa":      (Descriptors.TPSA,               True),
        "hbd":       (Descriptors.NumHDonors,         False),
        "hba":       (Descriptors.NumHAcceptors,      False),
        "rotbonds":  (Descriptors.NumRotatableBonds,  False),
        "heavy_atoms": (lambda m: m.GetNumHeavyAtoms(), False),
    }
    props_a, props_b, delta = {}, {}, {}
    for name, (fn, is_float) in fns.items():
        va, vb = fn(mol_a), fn(mol_b)
        props_a[name] = round(va, 4) if is_float else int(va)
        props_b[name] = round(vb, 4) if is_float else int(vb)
        diff = vb - va
        delta[name] = round(diff, 4) if is_float else int(diff)

    result = {
        "smiles_a":        Chem.MolToSmiles(mol_a),
        "smiles_b":        Chem.MolToSmiles(mol_b),
        "tanimoto":        tanimoto,
        "fingerprint":     args.fingerprint,
        "mcs": {
            "smarts":   mcs.smartsString,
            "n_atoms":  mcs.numAtoms,
            "n_bonds":  mcs.numBonds,
            "canceled": mcs.canceled,
        },
        "properties_a":    props_a,
        "properties_b":    props_b,
        "delta_b_minus_a": delta,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
