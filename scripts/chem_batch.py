#!/usr/bin/env python3
"""ALKYL — chem_batch.py: compute properties for a library of molecules."""

import argparse
import csv
import json
import sys
from pathlib import Path

DESCRIPTOR_MAP = {
    "mw":       ("rdkit.Chem.Descriptors", "MolWt"),
    "exact_mw": ("rdkit.Chem.Descriptors", "ExactMolWt"),
    "logp":     ("rdkit.Chem.Descriptors", "MolLogP"),
    "hbd":      ("rdkit.Chem.Descriptors", "NumHDonors"),
    "hba":      ("rdkit.Chem.Descriptors", "NumHAcceptors"),
    "tpsa":     ("rdkit.Chem.Descriptors", "TPSA"),
    "rotbonds": ("rdkit.Chem.Descriptors", "NumRotatableBonds"),
    "rings":    ("rdkit.Chem.Descriptors", "RingCount"),
    "fsp3":     ("rdkit.Chem.Descriptors", "FractionCSP3"),
}
ALL_DESCRIPTORS = list(DESCRIPTOR_MAP.keys())


def load_library(path: str, fmt: str) -> list[tuple[str, object]]:
    """Return list of (name, mol) pairs. mol is None for failed parses."""
    from rdkit import Chem
    p = Path(path)
    pairs = []

    if fmt == "sdf":
        suppl = Chem.SDMolSupplier(str(p))
        for i, mol in enumerate(suppl):
            if mol is not None and mol.HasProp("_Name"):
                name = mol.GetProp("_Name").strip() or f"mol_{i}"
            else:
                name = f"mol_{i}"
            pairs.append((name, mol))

    elif fmt == "smi":
        for i, line in enumerate(p.read_text().splitlines()):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{i}"
            pairs.append((name, Chem.MolFromSmiles(smi)))

    elif fmt == "csv":
        with open(str(p), newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames or []
            smiles_col = next(
                (c for c in fieldnames
                 if c.lower() in ("smiles", "smi", "canonical_smiles")), None)
            name_col = next(
                (c for c in fieldnames
                 if c.lower() in ("name", "id", "compound_id")), None)
            if smiles_col is None:
                print("No SMILES column found in CSV. "
                      "Expected: smiles, smi, or canonical_smiles.",
                      file=sys.stderr)
                sys.exit(1)
            for i, row in enumerate(reader):
                smi = row[smiles_col].strip()
                name = row[name_col].strip() if name_col else f"mol_{i}"
                pairs.append((name, Chem.MolFromSmiles(smi)))

    return pairs


def _calc_descriptors(mol, names: list[str]) -> dict:
    import importlib
    result = {}
    for name in names:
        if name not in DESCRIPTOR_MAP:
            continue
        module_path, func_name = DESCRIPTOR_MAP[name]
        module = importlib.import_module(module_path)
        result[name] = round(getattr(module, func_name)(mol), 4)
    return result


def _check_lipinski(mol) -> dict:
    from rdkit.Chem import Descriptors
    violations = []
    if Descriptors.MolWt(mol) > 500:
        violations.append("MW > 500")
    if Descriptors.MolLogP(mol) > 5:
        violations.append("LogP > 5")
    if Descriptors.NumHDonors(mol) > 5:
        violations.append("HBD > 5")
    if Descriptors.NumHAcceptors(mol) > 10:
        violations.append("HBA > 10")
    return {"pass": len(violations) == 0, "violations": violations}


def _check_pains(mol) -> dict:
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    matches = catalog.GetMatches(mol)
    alerts = [m.GetDescription() for m in matches]
    return {"alerts": alerts, "clean": len(alerts) == 0}


def process_mol(name: str, mol, args) -> dict | None:
    """Return a property record for one molecule, or None if skipping invalid."""
    if mol is None:
        if args.skip_invalid:
            return None
        return {"name": name, "error": "invalid_smiles"}

    from rdkit import Chem
    record: dict = {"name": name, "error": None,
                    "smiles": Chem.MolToSmiles(mol)}

    if args.descriptors:
        desc_names = (ALL_DESCRIPTORS if args.descriptors == "all"
                      else [d.strip() for d in args.descriptors.split(",")])
        record.update(_calc_descriptors(mol, desc_names))

    if args.lipinski:
        record["lipinski"] = _check_lipinski(mol)

    if args.pains:
        record["pains"] = _check_pains(mol)

    return record


def main():
    parser = argparse.ArgumentParser(
        description="Compute properties for a library of molecules.")
    parser.add_argument("--input", required=True,
                        help="Input file (SDF, SMILES, or CSV)")
    parser.add_argument("--format", choices=["sdf", "smi", "csv"],
                        default="smi", dest="fmt",
                        help="Input format (default: smi)")
    parser.add_argument("--descriptors",
                        help="Comma-separated descriptors or 'all'")
    parser.add_argument("--lipinski", action="store_true")
    parser.add_argument("--pains", action="store_true")
    parser.add_argument("--skip-invalid", action="store_true",
                        dest="skip_invalid",
                        help="Omit invalid molecules from output")
    parser.add_argument("--out", help="Write full JSON results to file; "
                        "stdout will contain a summary")
    args = parser.parse_args()

    if not any([args.descriptors, args.lipinski, args.pains]):
        parser.error("Specify at least one of: --descriptors, --lipinski, --pains")

    library = load_library(args.input, args.fmt)
    records = []
    for name, mol in library:
        rec = process_mol(name, mol, args)
        if rec is not None:
            records.append(rec)

    full_output = {
        "n_total": len(library),
        "n_processed": len(records),
        "n_errors": sum(1 for r in records if r.get("error")),
        "molecules": records,
    }

    if args.out:
        with open(args.out, "w") as f:
            json.dump(full_output, f, indent=2)
        # Always print a parseable summary to stdout
        summary = {k: v for k, v in full_output.items() if k != "molecules"}
        summary["out"] = args.out
        print(json.dumps(summary, indent=2))
    else:
        print(json.dumps(full_output, indent=2))


if __name__ == "__main__":
    main()
