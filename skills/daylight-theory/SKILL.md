---
name: daylight-theory
description: Use when working with SMILES, SMARTS, SMIRKS, molecular fingerprints, or cheminformatics fundamentals. Covers the complete Daylight theory: molecular graph representation, SMILES specification, SMARTS query language, SMIRKS reaction transforms, and fingerprint-based similarity. Based on the Daylight Theory Manual.
---

# Daylight Theory — Cheminformatics Fundamentals

The Daylight Theory Manual is the canonical reference for the molecular languages underlying modern cheminformatics: SMILES, SMARTS, SMIRKS, and fingerprints. These are not Daylight-proprietary — they are industry-standard formats implemented in RDKit, OpenBabel, CDK, and every major cheminformatics toolkit.

## Router — What to Read

| Topic | Reference |
|-------|-----------|
| Molecular graph model, aromaticity, chirality, SSSR, reaction representation | `references/molecules.md` |
| SMILES syntax: atoms, bonds, branches, rings, stereochemistry, reactions | `references/smiles.md` |
| SMARTS query language: primitives, operators, recursive SMARTS, reaction queries | `references/smarts.md` |
| SMIRKS reaction transforms: atom maps, grammar, stereochemistry | `references/smirks.md` |
| Fingerprints: structural keys, path-based FP, folding, Tanimoto, Tversky, all similarity measures | `references/fingerprints.md` |
| Chemical database concepts: hash tables, identifiers, in-memory search, pools, hitlists | `references/cheminformatics-databases.md` |

## Core Languages at a Glance

| Language | Purpose | Example |
|----------|---------|---------|
| SMILES | Encode a specific molecule | `CC(=O)Oc1ccccc1C(=O)O` |
| SMARTS | Describe a molecular pattern | `[OH]c1ccccc1` (phenol) |
| SMIRKS | Encode a reaction transform | `[C:1][Br:2]>>[C:1][I:2]` |

## Key Principles

- **Molecular graph**: atoms = nodes, bonds = edges; properties are explicit or derived
- **Aromaticity**: Hückel 4N+2 rule, all ring atoms sp² — a derived, not stored, property
- **SMILES organic subset**: B, C, N, O, P, S, F, Cl, Br, I — no brackets needed at normal valence
- **SMARTS unspecified = unrestricted**: `O` in SMARTS matches any aliphatic oxygen; in SMILES it means water
- **Tanimoto**: c/(a+b+c) — the standard fingerprint similarity; double-zero independent

## Related ALKYL Skills

- `rdkit` — implementation of these concepts in Python
- `scientific-skills:matchms` — spectrum similarity (same mathematical ideas)
