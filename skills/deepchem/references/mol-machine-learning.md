# Molecular Machine Learning with DeepChem

## Molecular Fingerprints

Fingerprints encode molecules as fixed-length bit vectors.

```python
import deepchem as dc
import numpy as np

# Circular (Morgan/ECFP) fingerprints
featurizer = dc.feat.CircularFingerprint(radius=2, size=2048)  # ECFP4
fps = featurizer.featurize(['CC(=O)O', 'c1ccccc1'])  # shape: (2, 2048)

# Extended connectivity variants
dc.feat.CircularFingerprint(radius=1, size=2048)  # ECFP2
dc.feat.CircularFingerprint(radius=3, size=2048)  # ECFP6

# RDKit descriptors (physicochemical)
featurizer = dc.feat.RDKitDescriptors()
descs = featurizer.featurize(['CC(=O)O'])  # ~200 descriptors

# MACCS keys (166-bit structural keys)
featurizer = dc.feat.MACCSKeysFingerprint()
```

**When to use fingerprints:** baseline models, interpretability needed, limited compute, sklearn compatibility.

```python
# Fingerprints work with sklearn models via DeepChem wrappers
from sklearn.ensemble import RandomForestClassifier
from deepchem.models import SklearnModel

sklearn_model = RandomForestClassifier(n_estimators=100)
model = dc.models.SklearnModel(sklearn_model)
model.fit(train_dataset)
```

## Graph Convolutional Networks

GCNs treat molecules as graphs: atoms = nodes, bonds = edges.

```python
# Featurize as molecular graph
featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
mol_graphs = featurizer.featurize(['CC(=O)O'])

# GraphConvModel — classic graph convolution
model = dc.models.GraphConvModel(
    n_tasks=1,
    mode='regression',      # or 'classification'
    dropout=0.2,
    batch_normalize=False,
    graph_conv_layers=[64, 64],
    dense_layer_size=128
)

# AttentiveFP — attention-based, state-of-the-art
model = dc.models.AttentiveFPModel(
    n_tasks=1,
    mode='regression',
    num_layers=2,
    num_timesteps=2,
    graph_feat_size=200
)

# GATModel — Graph Attention Network
model = dc.models.GATModel(n_tasks=1, mode='regression')

# MPNNModel — Message Passing Neural Network
model = dc.models.MPNNModel(n_tasks=1, mode='regression')
```

## Going Deeper: Advanced Featurizations

```python
# Weave featurizer — encodes atom pairs
featurizer = dc.feat.WeaveFeaturizer()

# CoulombMatrix — for quantum chemistry, 3D structure
featurizer = dc.feat.CoulombMatrix(max_atoms=23)

# Atomic contributions (SHAP-compatible)
featurizer = dc.feat.AtomicConvFeaturizer(
    frag1_num_atoms=70,
    frag2_num_atoms=634,
    complex_num_atoms=701
)
```

## Unsupervised Molecular Embeddings

```python
# Mol2Vec embeddings (pretrained)
featurizer = dc.feat.Mol2VecFingerprint()
embeddings = featurizer.featurize(['CC(=O)O'])

# Use embeddings as features for downstream model
model = dc.models.MultitaskRegressor(
    n_tasks=1,
    n_features=300,  # mol2vec dim
    layer_sizes=[512, 256]
)
```

## Transfer Learning with ChemBERTa

```python
# ChemBERTa: BERT pretrained on SMILES strings
model = dc.models.ChemBERTaModel(
    n_tasks=1,
    mode='regression',
    model_type='ChemBERTa-77M-MLM'  # or 'ChemBERTa-77M-MTR'
)

# Fine-tune on your dataset
tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='dummy')
model.fit(datasets[0], nb_epoch=5)
```

## Molecular Attention Transformer (MAT)

```python
# MAT: combines molecular graph + distance matrix + self-attention
featurizer = dc.feat.MATFeaturizer()
model = dc.models.MATModel(n_tasks=1, mode='regression')
```

## Generative Models

```python
# MolGAN: generates valid molecules via GAN
model = dc.models.MolGANModel(vertices=9, edges=5, nodes=5)
model.fit(train_dataset, nb_epoch=100)
generated = model.predict_generated_molecules(num_samples=100)

# Normalizing Flow on QM9
from deepchem.models import NormalizingFlowModel
```

## Model Interpretability / Atomic Contributions

```python
# After training, compute per-atom contributions
atom_shap = dc.utils.rdkit_utils.compute_atom_shap(model, mol_graph)

# Visualize with RDKit
from rdkit.Chem.Draw import SimilarityMaps
```

## Large-Scale Chemical Screening

```python
# For virtual screening of large compound libraries
# Use batch prediction for efficiency

compounds = load_large_library()  # millions of SMILES
featurizer = dc.feat.CircularFingerprint(size=2048)

# Process in chunks
chunk_size = 10000
for i in range(0, len(compounds), chunk_size):
    chunk = compounds[i:i+chunk_size]
    feats = featurizer.featurize(chunk)
    preds = model.predict_on_batch(feats)
    # rank and filter
```
