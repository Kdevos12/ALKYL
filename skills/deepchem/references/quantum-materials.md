# Quantum Chemistry & Materials Science with DeepChem

## QM9 Dataset — 12 Quantum Properties

QM9 contains 133,885 small organic molecules with 12 DFT-computed properties.

```python
import deepchem as dc

tasks, datasets, transformers = dc.molnet.load_qm9(
    featurizer='CoulombMatrix',
    splitter='random'
)

# Tasks: mu, alpha, homo, lumo, gap, r2, zpve, U0, U, H, G, Cv
# (dipole moment, polarizability, HOMO/LUMO energies, etc.)
print(tasks)

model = dc.models.MPNNModel(
    n_tasks=len(tasks),
    mode='regression',
    learning_rate=1e-3
)
model.fit(datasets[0], nb_epoch=200)
```

### QM9 Properties Reference

| Property | Unit | Description |
|----------|------|-------------|
| `mu` | Debye | Dipole moment |
| `alpha` | a0^3 | Isotropic polarizability |
| `homo` | Ha | HOMO energy |
| `lumo` | Ha | LUMO energy |
| `gap` | Ha | HOMO-LUMO gap |
| `r2` | a0^2 | Electronic spatial extent |
| `zpve` | Ha | Zero point vibrational energy |
| `U0` | Ha | Internal energy at 0K |
| `U` | Ha | Internal energy at 298K |
| `H` | Ha | Enthalpy at 298K |
| `G` | Ha | Free energy at 298K |
| `Cv` | cal/mol/K | Heat capacity |

## Featurizers for Quantum Chemistry

```python
# CoulombMatrix — encodes pairwise nuclear interactions
featurizer = dc.feat.CoulombMatrix(max_atoms=23)

# CoulombMatrixEig — eigenvalue representation (rotation invariant)
featurizer = dc.feat.CoulombMatrixEig(max_atoms=23)

# BPSymmetryFunctionInput — Behler-Parrinello symmetry functions
featurizer = dc.feat.BPSymmetryFunctionInput(max_atoms=23)

# Usage
from rdkit import Chem
mol = Chem.MolFromSmiles('CC(=O)O')
feat = featurizer.featurize([mol])
```

## DeepQMC — Neural Network Potential

```python
# Train neural network potentials on quantum chemistry data
# For learning PES (potential energy surfaces)

from deepchem.models import DeepQMCModel

model = DeepQMCModel(
    n_tasks=1,           # energy prediction
    learning_rate=1e-3
)

# Requires 3D coordinates (not just SMILES)
tasks, datasets, transformers = dc.molnet.load_qm7(
    featurizer='CoulombMatrix',
    move_mean=True
)
model.fit(datasets[0], nb_epoch=100)
```

## Exchange-Correlation Functional with DeepChem

```python
# Train ML exchange-correlation functionals for DFT
# DeepChem interfaces with PySCF for DFT calculations

from deepchem.models import KernelRidgeRegression

# Feature: electron density descriptors
# Target: XC energy
```

## Physics-Informed Neural Networks (PINNs)

```python
# Solve differential equations with neural networks
# Useful for quantum mechanics and fluid dynamics

from deepchem.models import PINNModel, JaxModel

# Solve Schrodinger equation
model = PINNModel(
    n_tasks=1,
    pde_loss_weight=1.0,
    boundary_loss_weight=10.0
)
```

## Neural ODEs

```python
# Torchdiffeq integration for continuous dynamics
from deepchem.models import NeuralODEModel

# Model molecular dynamics trajectories
# Learn continuous-time ODE from observations
```

## Materials Science

```python
# Crystal structure property prediction
# Requires periodic boundary conditions

# CGCNN: Crystal Graph CNN
from deepchem.feat import MaterialStructureFeaturizer

featurizer = MaterialStructureFeaturizer()

# Load materials dataset
tasks, datasets, transformers = dc.molnet.load_band_gap(
    featurizer='cgcnn',
    splitter='random'
)

model = dc.models.CGCNNModel(n_tasks=1, mode='regression')
model.fit(datasets[0], nb_epoch=100)
```

## Best Practices for Quantum Chemistry ML

1. **Always normalize targets** — quantum properties span orders of magnitude
2. **Use CoulombMatrixEig** over CoulombMatrix when rotation invariance matters
3. **max_atoms must be >= largest molecule** in your dataset — pad otherwise
4. **Train/test split**: use random for QM datasets (scaffold less meaningful for small molecules)
5. **Evaluate in physical units** after untransforming predictions
6. **MPNN and AttentiveFP** generally outperform GraphConv on quantum tasks
7. **Baseline**: always compare against DFT-level theory, not just ML baselines

```python
# Proper evaluation flow for quantum properties
metric = dc.metrics.Metric(dc.metrics.mean_absolute_error)
scores = model.evaluate(test_dataset, [metric], transformers)

# Untransform to get MAE in original units (Hartrees, Debye, etc.)
preds = model.predict(test_dataset)
preds_physical = transformers[0].untransform(preds)
true_physical   = transformers[0].untransform(test_dataset.y)
mae = np.mean(np.abs(preds_physical - true_physical))
```
