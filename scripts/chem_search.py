#!/usr/bin/env python3
"""ALKYL — chem_search.py: substructure, similarity, and exact search in a molecular library."""

import argparse
import json
import sys
from pathlib import Path


def load_query(args):
    """Return (mol, kind) where kind is 'smarts' or 'smiles'."""
    from rdkit import Chem
    if args.smarts:
        mol = Chem.MolFromSmarts(args.smarts)
        if mol is None:
            print(f"Invalid SMARTS: {args.smarts}", file=sys.stderr)
            sys.exit(1)
        return mol, "smarts"
    mol = Chem.MolFromSmiles(args.query)
    if mol is None:
        print(f"Invalid query SMILES: {args.query}", file=sys.stderr)
        sys.exit(1)
    return mol, "smiles"


def load_library(path: str) -> list[tuple[int, str, object]]:
    """Load SDF or SMILES library. Returns list of (idx, name, mol)."""
    from rdkit import Chem
    p = Path(path)
    mols = []

    if p.suffix.lower() in (".sdf", ".mol"):
        suppl = Chem.SDMolSupplier(str(p))
        for i, mol in enumerate(suppl):
            if mol is not None and mol.HasProp("_Name"):
                name = mol.GetProp("_Name").strip() or f"mol_{i}"
            else:
                name = f"mol_{i}"
            mols.append((i, name, mol))
    else:
        # Plain SMILES file: one SMILES per line, optional name after whitespace
        for i, line in enumerate(p.read_text().splitlines()):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{i}"
            mols.append((i, name, Chem.MolFromSmiles(smi)))

    return mols


def _to_smiles(mol) -> str:
    from rdkit import Chem
    return Chem.MolToSmiles(mol)


def _fingerprint(mol, fp_type: str, radius: int, nbits: int):
    if fp_type == "morgan":
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
        return gen.GetFingerprint(mol)
    if fp_type == "maccs":
        from rdkit.Chem import MACCSkeys
        return MACCSkeys.GenMACCSKeys(mol)
    # rdkit fingerprint
    from rdkit.Chem.rdmolops import RDKFingerprint
    return RDKFingerprint(mol, fpSize=nbits)


def search_substructure(query, library: list) -> list[dict]:
    hits = []
    for idx, name, mol in library:
        if mol is None:
            continue
        if mol.HasSubstructMatch(query):
            hits.append({"idx": idx, "name": name, "smiles": _to_smiles(mol)})
    return hits


def search_exact(query, library: list) -> list[dict]:
    from rdkit import Chem
    query_smi = Chem.MolToSmiles(query)
    hits = []
    for idx, name, mol in library:
        if mol is None:
            continue
        if Chem.MolToSmiles(mol) == query_smi:
            hits.append({"idx": idx, "name": name, "smiles": _to_smiles(mol)})
    return hits


def search_similarity(query, library: list, fp_type: str,
                       threshold: float, top_k: int | None,
                       radius: int, nbits: int) -> list[dict]:
    from rdkit import DataStructs
    query_fp = _fingerprint(query, fp_type, radius, nbits)
    scored = []
    for idx, name, mol in library:
        if mol is None:
            continue
        mol_fp = _fingerprint(mol, fp_type, radius, nbits)
        sim = DataStructs.TanimotoSimilarity(query_fp, mol_fp)
        if sim >= threshold:
            scored.append({
                "idx": idx,
                "name": name,
                "smiles": _to_smiles(mol),
                "tanimoto": round(sim, 4),
            })
    scored.sort(key=lambda x: x["tanimoto"], reverse=True)
    if top_k is not None:
        scored = scored[:top_k]
    return scored


def main():
    parser = argparse.ArgumentParser(
        description="Search a molecular library by substructure, similarity, or exact match.")

    qgrp = parser.add_mutually_exclusive_group(required=True)
    qgrp.add_argument("--query", help="Query SMILES")
    qgrp.add_argument("--smarts", help="Query SMARTS pattern (substructure mode only)")

    parser.add_argument("--library", required=True,
                        help="Library file (.sdf/.mol or SMILES file)")
    parser.add_argument("--mode", required=True,
                        choices=["substructure", "similarity", "exact"],
                        help="Search mode")

    # Similarity options
    parser.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"],
                        default="morgan",
                        help="Fingerprint type for similarity (default: morgan)")
    parser.add_argument("--radius", type=int, default=2,
                        help="Morgan radius (default: 2)")
    parser.add_argument("--nbits", type=int, default=2048,
                        help="Fingerprint size in bits (default: 2048)")
    parser.add_argument("--threshold", type=float, default=0.7,
                        help="Minimum Tanimoto similarity (default: 0.7)")
    parser.add_argument("--top-k", type=int, dest="top_k",
                        help="Return at most K hits (similarity mode)")

    parser.add_argument("--out", help="Write JSON results to file")
    args = parser.parse_args()

    if args.smarts and args.mode != "substructure":
        parser.error("--smarts is only valid with --mode substructure")

    query, query_kind = load_query(args)
    library = load_library(args.library)

    if args.mode == "substructure":
        hits = search_substructure(query, library)
    elif args.mode == "exact":
        if query_kind == "smarts":
            parser.error("--smarts cannot be used with --mode exact")
        hits = search_exact(query, library)
    else:  # similarity
        if query_kind == "smarts":
            parser.error("--smarts cannot be used with --mode similarity")
        hits = search_similarity(
            query, library, args.fingerprint,
            args.threshold, args.top_k, args.radius, args.nbits)

    result = {
        "mode": args.mode,
        "n_library": len(library),
        "n_hits": len(hits),
        "hits": hits,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
