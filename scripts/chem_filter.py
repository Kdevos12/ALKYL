#!/usr/bin/env python3
"""ALKYL — chem_filter.py: multi-rule drug-likeness filtering (Lipinski, Veber, Egan, Ghose)."""

import argparse
import json
import os
import sys


def _load_sascorer():
    import rdkit
    contrib_path = os.path.join(os.path.dirname(rdkit.__file__), "Contrib", "SA_Score")
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
    suppl = Chem.SDMolSupplier(args.sdf)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        print(f"No valid molecule in: {args.sdf}", file=sys.stderr)
        sys.exit(1)
    return mol


def _lipinski(d) -> dict:
    v = []
    if d["mw"] > 500:    v.append("MW > 500")
    if d["logp"] > 5:    v.append("LogP > 5")
    if d["hbd"] > 5:     v.append("HBD > 5")
    if d["hba"] > 10:    v.append("HBA > 10")
    return {"pass": not v, "violations": v, "rule": "Lipinski Ro5"}


def _veber(d) -> dict:
    v = []
    if d["tpsa"] > 140:    v.append("TPSA > 140")
    if d["rotbonds"] > 10: v.append("RotBonds > 10")
    return {"pass": not v, "violations": v, "rule": "Veber (oral bioavailability)"}


def _egan(d) -> dict:
    v = []
    if d["tpsa"] > 132:   v.append("TPSA > 132")
    if d["logp"] > 5.88:  v.append("LogP > 5.88")
    return {"pass": not v, "violations": v, "rule": "Egan (absorption)"}


def _ghose(mol, d) -> dict:
    from rdkit.Chem.Descriptors import MolMR
    mr = round(MolMR(mol), 2)
    heavy = mol.GetNumHeavyAtoms()
    v = []
    if not (-0.4 <= d["logp"] <= 5.6):  v.append(f"LogP {d['logp']:.2f} ∉ [-0.4, 5.6]")
    if not (160 <= d["mw"] <= 480):     v.append(f"MW {d['mw']:.1f} ∉ [160, 480]")
    if not (40 <= mr <= 130):           v.append(f"MR {mr} ∉ [40, 130]")
    if not (20 <= heavy <= 70):         v.append(f"Heavy atoms {heavy} ∉ [20, 70]")
    return {"pass": not v, "violations": v, "mr": mr, "rule": "Ghose"}


def _pains(mol) -> dict:
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    matches = catalog.GetMatches(mol)
    alerts = [m.GetDescription() for m in matches]
    return {"pass": not alerts, "alerts": alerts, "rule": "PAINS"}


def main():
    parser = argparse.ArgumentParser(
        description="Multi-rule drug-likeness filtering.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED

    mol = load_mol(args)
    sascorer = _load_sascorer()

    d = {
        "mw":       round(Descriptors.MolWt(mol), 2),
        "logp":     round(Descriptors.MolLogP(mol), 4),
        "hbd":      Descriptors.NumHDonors(mol),
        "hba":      Descriptors.NumHAcceptors(mol),
        "tpsa":     round(Descriptors.TPSA(mol), 2),
        "rotbonds": Descriptors.NumRotatableBonds(mol),
    }

    rules = {
        "lipinski": _lipinski(d),
        "veber":    _veber(d),
        "egan":     _egan(d),
        "ghose":    _ghose(mol, d),
        "pains":    _pains(mol),
    }
    n_pass = sum(1 for r in rules.values() if r["pass"])

    result = {
        "smiles":        Chem.MolToSmiles(mol),
        "descriptors":   d,
        "qed":           round(QED.qed(mol), 4),
        "sa_score":      round(sascorer.calculateScore(mol), 4),
        "rules":         rules,
        "n_rules_pass":  n_pass,
        "n_rules_total": len(rules),
        "overall_pass":  n_pass == len(rules),
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
