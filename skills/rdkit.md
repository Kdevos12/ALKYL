# RDKit — Cheminformatics Skill

## Overview
RDKit is the primary Python library for cheminformatics in ALKYL workflows.
Use it for molecular manipulation, descriptor calculation, and visualization.

## Common Patterns

### Loading a molecule
```python
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors

mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')  # aspirin
```

### Calculating descriptors
```python
mw = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
hbd = Descriptors.NumHDonors(mol)
hba = Descriptors.NumHAcceptors(mol)
```

### Substructure search
```python
pattern = Chem.MolFromSmarts('[C:1](=O)[OH]')  # carboxylic acid
matches = mol.GetSubstructMatches(pattern)
```

### 2D Visualization
```python
from rdkit.Chem import Draw
img = Draw.MolToImage(mol, size=(300, 300))
```

## Installation
```bash
pip install rdkit
# or via conda:
conda install -c conda-forge rdkit
```
