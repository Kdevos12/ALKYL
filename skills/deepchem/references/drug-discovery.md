# Drug Discovery with DeepChem

## MoleculeNet Datasets

MoleculeNet is DeepChem's curated benchmark suite. All accessible via `dc.molnet.*`.

### Key Datasets

| Dataset | Task | Samples | Use case |
|---------|------|---------|----------|
| `load_delaney` | Solubility (regression) | 1,128 | Baseline benchmark |
| `load_tox21` | Toxicity 12 tasks (clf) | 8,014 | Multitask toxicity |
| `load_bbbp` | Blood-brain barrier (clf) | 2,050 | CNS drug discovery |
| `load_hiv` | HIV replication (clf) | 41,127 | Large-scale clf |
| `load_bace_classification` | BACE inhibition (clf) | 1,522 | Protease inhibition |
| `load_bace_regression` | BACE binding (reg) | 1,522 | Binding affinity |
| `load_clintox` | Clinical toxicity (clf) | 1,491 | DMPK |
| `load_sider` | Side effects 27 tasks (clf) | 1,427 | Side effect prediction |
| `load_qm7` | Electronic energy (reg) | 7,165 | Quantum properties |
| `load_qm9` | 12 quantum props (reg) | 133,885 | Quantum chemistry |
| `load_pdbbind` | Protein-ligand affinity | 11,908 | Structure-based |

```python
# Full list of available loaders
import deepchem as dc
loaders = [m for m in dir(dc.molnet) if 'load_' in m]

# Load with options
tasks, datasets, transformers = dc.molnet.load_tox21(
    featurizer='ECFP',
    splitter='scaffold',
    transformers=['balancing']
)
```

### Splitter Recommendations by Dataset

| Dataset | Recommended splitter | Reason |
|---------|---------------------|--------|
| Tox21, HIV | `scaffold` | Chemotype generalization |
| QM7, QM9 | `random` | Diverse quantum structures |
| PDBBind | `random` | Limited size |
| BBBP, BACE | `scaffold` | Drug-like molecules |

## Protein-Ligand Interaction Modeling

```python
# AtomicConv: 3D structure-based binding affinity
featurizer = dc.feat.AtomicConvFeaturizer(
    frag1_num_atoms=70,
    frag2_num_atoms=634,
    complex_num_atoms=701,
    max_num_neighbors=12,
    neighbor_cutoff=4.0
)

model = dc.models.AtomicConvModel(
    n_tasks=1,
    batch_size=24,
    layer_sizes=[32, 32, 16],
    learning_rate=0.003
)

# Load PDBBind
tasks, datasets, transformers = dc.molnet.load_pdbbind(
    featurizer='atomic',
    splitter='random',
    subset='refined'
)
model.fit(datasets[0], nb_epoch=30)
```

### DeepChem + AlphaFold Integration

```python
# DeepChemXAlphafold: use AlphaFold predicted structures
# for protein-ligand modeling when experimental structure unavailable
from deepchem.models.alphafold import AlphaFoldInteractionModel
```

## Virtual Screening Pipeline

```python
import deepchem as dc
import numpy as np

# 1. Train on known actives/inactives
tasks, datasets, transformers = dc.molnet.load_bace_classification(
    featurizer='GraphConv',
    splitter='scaffold'
)
train, valid, test = datasets

model = dc.models.GraphConvModel(n_tasks=1, mode='classification', dropout=0.2)
model.fit(train, nb_epoch=50)

# Evaluate
metric = dc.metrics.Metric(dc.metrics.roc_auc_score)
print(model.evaluate(test, [metric], transformers))

# 2. Screen a library
library_smiles = load_compound_library()  # your SMILES list
featurizer = dc.feat.MolGraphConvFeaturizer()
library_feats = featurizer.featurize(library_smiles)
library_dataset = dc.data.NumpyDataset(X=library_feats, ids=library_smiles)

# 3. Score all compounds
scores = model.predict(library_dataset)
# scores shape: (n_compounds, n_tasks, n_classes)
activity_scores = scores[:, 0, 1]  # probability of class 1

# 4. Rank and select top hits
top_indices = np.argsort(activity_scores)[::-1][:100]
hits = [library_smiles[i] for i in top_indices]
```

## Conditional Molecular Generation (cGAN)

```python
# Generate molecules with desired properties
from deepchem.models import BasicMolGANModel

model = BasicMolGANModel(
    la=1, edges=5, vertices=9, nodes=5,
    embedding_dim=10, dropout=0.0
)
model.fit(train_dataset, nb_epoch=100)

# Generate molecules
generated = model.predict_generated_molecules(num=100)
# Filter for validity, uniqueness, drug-likeness
```

## Multitask Learning for Drug Discovery

**Key insight from the book:** Multitask models significantly outperform single-task models in drug discovery because biological targets share structural determinants.

```python
# Tox21: 12 toxicity endpoints simultaneously
tasks, datasets, transformers = dc.molnet.load_tox21(featurizer='ECFP')

model = dc.models.MultitaskClassifier(
    n_tasks=12,
    n_features=1024,
    layer_sizes=[1000, 500],
    dropouts=[0.25, 0.25],
    weight_init_stddevs=[0.02, 0.02],
    bias_init_consts=[1.0, 1.0],
    batch_size=100
)

# The model learns shared representations across all 12 toxicity tasks
# Improves performance on low-data tasks via transfer across tasks
```
