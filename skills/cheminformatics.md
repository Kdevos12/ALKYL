# Cheminformatics Workflows

## Molecular Representations

| Format | Use case | Example |
|--------|----------|---------|
| SMILES | Text, databases | `CC(=O)O` |
| InChI  | Unique identifier | `InChI=1S/C2H4O2/...` |
| SDF/MOL | 3D structure storage | File format |
| PDB | Protein structures | File format |

## SMILES Quick Reference
- Single bond: `CC`
- Double bond: `C=C`
- Ring: `c1ccccc1` (benzene)
- Branch: `CC(C)C` (isobutane)
- Charge: `[NH4+]`, `[O-]`

## Key Databases
- **PubChem**: general compounds → `https://pubchem.ncbi.nlm.nih.gov`
- **ChEMBL**: bioactive molecules → `https://www.ebi.ac.uk/chembl`
- **RCSB PDB**: protein structures → `https://www.rcsb.org`
- **ZINC**: purchasable compounds → `https://zinc.docking.org`

## Lipinski's Rule of Five (drug-likeness)
- MW ≤ 500 Da
- LogP ≤ 5
- H-bond donors ≤ 5
- H-bond acceptors ≤ 10

## Common Python Stack
```python
import rdkit          # molecular manipulation
import ase            # atomistic simulations
import openbabel      # format conversion
import py3Dmol        # 3D visualization (Jupyter)
import mdanalysis     # MD trajectory analysis
```
