#!/usr/bin/env python3
"""ALKYL — chem_analyze.py: comprehensive structural analysis of a molecule."""

import argparse
import json
import sys
import os

# SMARTS for functional group detection
# Counts number of matches (not just presence)
FUNCTIONAL_GROUP_SMARTS = {
    "carboxylic_acid":  "[CX3](=O)[OX2H1]",
    "ester":            "[CX3](=O)[OX2][#6]",
    "amide":            "[CX3](=O)[NX3]",
    "amine_primary":    "[NX3;H2;!$(NC=O)]",
    "amine_secondary":  "[NX3;H1;!$(NC=O)]",
    "amine_tertiary":   "[NX3;H0;!$(NC=O);!$([n])]",
    "alcohol":          "[OX2H1;!$(OC=O)]",
    "aldehyde":         "[CX3H1](=O)",
    "ketone":           "[CX3;!H](=O)[#6;!$(C=[OX1])]",
    "nitro":            "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
    "sulfonamide":      "[#16X4](=[OX1])(=[OX1])[NX3]",
    "sulfone":          "[#16X4](=[OX1])(=[OX1])([#6])[#6]",
    "halogen_F":        "[F]",
    "halogen_Cl":       "[Cl]",
    "halogen_Br":       "[Br]",
    "halogen_I":        "[I]",
}


def _load_sascorer():
    """Locate and import sascorer from RDKit Contrib."""
    import rdkit
    contrib_path = os.path.join(os.path.dirname(rdkit.__file__),
                                "Contrib", "SA_Score")
    if contrib_path not in sys.path:
        sys.path.insert(0, contrib_path)
    import sascorer
    return sascorer


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


def analyze_functional_groups(mol) -> dict:
    from rdkit import Chem
    groups = {}
    for name, smarts in FUNCTIONAL_GROUP_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        count = len(mol.GetSubstructMatches(pattern)) if pattern else 0
        if count > 0:
            groups[name] = count
    return groups


def analyze_rings(mol) -> dict:
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    aromatic = sum(
        1 for ring in rings
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
    )
    return {
        "aromatic": aromatic,
        "aliphatic": len(rings) - aromatic,
        "total": len(rings),
    }


def analyze_stereo(mol) -> dict:
    from rdkit.Chem import FindMolChiralCenters
    centers = FindMolChiralCenters(mol, includeUnassigned=True)
    assigned = sum(1 for _, cfg in centers if cfg in ("R", "S"))
    unassigned = len(centers) - assigned
    return {
        "total": len(centers),
        "assigned": assigned,
        "unassigned": unassigned,
    }


def analyze(mol) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from rdkit.Chem.GraphDescriptors import BertzCT

    sascorer = _load_sascorer()

    result = {
        "smiles":           Chem.MolToSmiles(mol),
        "formula":          CalcMolFormula(mol),
        "heavy_atoms":      mol.GetNumHeavyAtoms(),
        "charge":           sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        "mw":               round(Descriptors.MolWt(mol), 4),
        "exact_mw":         round(Descriptors.ExactMolWt(mol), 6),
        "hbd":              Descriptors.NumHDonors(mol),
        "hba":              Descriptors.NumHAcceptors(mol),
        "rotbonds":         Descriptors.NumRotatableBonds(mol),
        "tpsa":             round(Descriptors.TPSA(mol), 2),
        "logp":             round(Descriptors.MolLogP(mol), 4),
        "fsp3":             round(Descriptors.FractionCSP3(mol), 4),
        "rings":            analyze_rings(mol),
        "stereocenters":    analyze_stereo(mol),
        "functional_groups": analyze_functional_groups(mol),
        "qed":              round(QED.qed(mol), 4),
        "sa_score":         round(sascorer.calculateScore(mol), 4),
        "complexity_bertz": round(BertzCT(mol), 1),
    }
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive structural analysis of a molecule.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    mol = load_mol(args)
    result = analyze(mol)

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
