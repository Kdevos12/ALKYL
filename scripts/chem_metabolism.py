#!/usr/bin/env python3
"""ALKYL — chem_metabolism.py: CYP450 metabolic soft spot prediction via SMARTS."""

import argparse
import json
import sys

# (enzyme, pattern_name, SMARTS, risk_level)
METABOLISM_RULES = [
    # CYP3A4
    ("CYP3A4", "aromatic_hydroxylation",  "c1ccccc1",           "medium"),
    ("CYP3A4", "aliphatic_hydroxylation", "[CH2][CH2][CH3]",    "low"),
    ("CYP3A4", "n_dealkylation",          "[N,n][CH3]",         "high"),
    ("CYP3A4", "o_dealkylation",          "[OX2][CH3]",         "medium"),
    # CYP2D6
    ("CYP2D6", "aromatic_hydroxylation",  "c1cccnc1",           "medium"),
    ("CYP2D6", "n_dealkylation",          "[N,n]CC",            "medium"),
    # CYP2C9
    ("CYP2C9", "aromatic_hydroxylation",  "c1ccc(F)cc1",        "low"),
    ("CYP2C9", "sulfoxidation",           "[SX2][CH3]",         "medium"),
    # CYP1A2
    ("CYP1A2", "n_hydroxylation",         "[NX3;H2]c",          "high"),
    ("CYP1A2", "aromatic_hydroxylation",  "c1cccc2ccccc12",     "medium"),
    # Phase II
    ("UGT",    "glucuronidation",         "[OX2H1;!$(OC=O)]",   "medium"),
    ("SULT",   "sulfation",               "[OX2H1;!$(OC=O)]",   "low"),
]

_RISK_RANK = {"none": 0, "low": 1, "medium": 2, "high": 3}


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
        description="Predict CYP450 metabolic soft spots via SMARTS rules.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles")
    inp.add_argument("--sdf")
    parser.add_argument("--out")
    args = parser.parse_args()

    from rdkit import Chem

    mol = load_mol(args)
    canonical = Chem.MolToSmiles(mol)

    sites = []
    overall_risk = "none"

    for (enzyme, pattern_name, smarts, risk) in METABOLISM_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            sites.append({
                "enzyme": enzyme,
                "pattern": pattern_name,
                "atom_indices": list(match),
                "risk": risk,
            })
            if _RISK_RANK[risk] > _RISK_RANK[overall_risk]:
                overall_risk = risk

    result = {
        "smiles": canonical,
        "metabolic_sites": sites,
        "n_sites": len(sites),
        "overall_risk": overall_risk,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
