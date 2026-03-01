# DeepChem Core Concepts

## Dataset Types

| Class | Use case |
|-------|----------|
| `DiskDataset` | Large datasets saved to disk, lazy loading |
| `NumpyDataset` | Small/medium datasets, all in memory |
| `ImageDataset` | Image-based inputs/outputs |

## Dataset Anatomy

Every DeepChem dataset stores four arrays per sample:

| Field | Name | Description |
|-------|------|-------------|
| `X` | features | Input features (featurized molecules) |
| `y` | labels | Target values (property to predict) |
| `w` | weights | Sample weights (1.0 by default) |
| `ids` | identifiers | Usually SMILES strings |
| `task_names` | task names | Names of prediction targets |

```python
print(test_dataset)
# DiskDataset X.shape: (113,), y.shape: (113, 1), w.shape: (113, 1)
# ids: ['C1c2ccccc...'], task_names: ['measured log solubility']
```

## Accessing Data

```python
# Direct array access (loads all into memory — careful with large datasets)
X = test_dataset.X
y = test_dataset.y
w = test_dataset.w
ids = test_dataset.ids

# Iterate sample by sample (memory-efficient)
for X, y, w, id in test_dataset.itersamples():
    print(y, id)

# Iterate in batches (recommended for training)
for X, y, w, ids in test_dataset.iterbatches(batch_size=50):
    print(y.shape)
# iterbatches(batch_size=100, epochs=10, deterministic=False)

# Convert to Pandas DataFrame (small datasets only)
df = test_dataset.to_dataframe()

# Framework integrations
tf_dataset  = test_dataset.make_tf_dataset()
pt_dataset  = test_dataset.make_pytorch_dataset()
```

## Creating Datasets

```python
import numpy as np
import deepchem as dc

# NumpyDataset — from arrays
X = np.random.random((100, 5))
y = np.random.random((100, 2))
dataset = dc.data.NumpyDataset(X=X, y=y)  # w and ids auto-generated

# DiskDataset — from numpy arrays (persisted)
import tempfile
with tempfile.TemporaryDirectory() as data_dir:
    disk_dataset = dc.data.DiskDataset.from_numpy(X=X, y=y, data_dir=data_dir)
```

## Featurizers

Featurizers convert SMILES strings -> ML-ready features.

```python
# Passed to molnet loaders
tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='GraphConv')
tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='ECFP')
tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='RDKit')

# Used directly
from deepchem.feat import MolGraphConvFeaturizer, CircularFingerprint

featurizer = MolGraphConvFeaturizer(use_edges=True)
mol_graph = featurizer.featurize(['CC(=O)O', 'c1ccccc1'])

featurizer = CircularFingerprint(radius=2, size=2048)  # ECFP4
fps = featurizer.featurize(['CC(=O)O'])
```

| Featurizer | Input | Output | Best for |
|------------|-------|--------|----------|
| `GraphConv` / `MolGraphConvFeaturizer` | SMILES | Molecular graph | GCN models |
| `CircularFingerprint` | SMILES | Bit vector (ECFP) | Random forests, shallow models |
| `RDKitDescriptors` | SMILES | 200 physicochemical descriptors | Interpretable models |
| `CoulombMatrix` | XYZ coordinates | Coulomb matrix | Quantum chemistry |
| `AtomicCoordinates` | PDB/XYZ | 3D coords | Protein-ligand, 3D models |

## Splitters

Splitters divide data into train/valid/test sets.

```python
# Passed to molnet loaders
tasks, datasets, transformers = dc.molnet.load_delaney(splitter='scaffold')

# Available splitters
'random'     # random split — baseline, fast
'scaffold'   # Bemis-Murcko scaffold — tests generalization to new chemotypes
'stratified' # preserves label distribution — for imbalanced classification
'index'      # sequential split by index
'butina'     # Butina clustering — similar to scaffold
'fingerprint'# fingerprint similarity-based split
```

**Rule of thumb:** Use `scaffold` for drug discovery (realistic generalization). Use `random` for quick experiments.

```python
# Manual splitting
splitter = dc.splits.ScaffoldSplitter()
train, valid, test = splitter.train_valid_test_split(dataset, frac_train=0.8, frac_valid=0.1, frac_test=0.1)
```

## Transformers

```python
# Applied automatically when using molnet loaders
# Must be applied to predictions when interpreting results

# Common transformers
dc.trans.NormalizationTransformer   # normalize y values
dc.trans.ClippingTransformer        # clip outliers
dc.trans.BalancingTransformer       # balance class weights

# Untransform predictions
predictions = model.predict(test_dataset)
untransformed = transformers[0].untransform(predictions)
```
