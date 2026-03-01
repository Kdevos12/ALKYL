---
name: deepchem
description: Use when working with DeepChem for molecular machine learning, drug discovery, quantum chemistry, materials science, or bioinformatics. Handles molecular datasets, featurization strategies, model training/evaluation, and predictions on chemical data.
---

# DeepChem

Deep learning for the life sciences: drug discovery, quantum chemistry, materials science, bioinformatics.

## When to Use This Skill

- Building ML models on molecular datasets (SMILES, graphs, fingerprints)
- Working with MoleculeNet benchmark datasets
- Predicting molecular properties (solubility, toxicity, binding affinity)
- Protein-ligand interaction modeling
- Quantum chemistry property prediction (QM9, GDB datasets)
- Featurizing molecules for downstream ML tasks
- Virtual screening and drug discovery pipelines

## Quick Start — Standard Workflow

```python
import deepchem as dc

# 1. Load dataset with featurizer
tasks, datasets, transformers = dc.molnet.load_delaney(featurizer='GraphConv')
train_dataset, valid_dataset, test_dataset = datasets

# 2. Create model
model = dc.models.GraphConvModel(n_tasks=1, mode='regression', dropout=0.2)

# 3. Train
model.fit(train_dataset, nb_epoch=100)

# 4. Evaluate
metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)
train_score = model.evaluate(train_dataset, [metric], transformers)
test_score  = model.evaluate(test_dataset,  [metric], transformers)

# 5. Predict
predictions = model.predict_on_batch(test_dataset.X[:10])
```

## Router — What to Read

| Task | Reference |
|------|-----------|
| Dataset creation, access, splitters | `references/core-concepts.md` |
| Training workflow, metrics, hyperopt, multitask | `references/model-training.md` |
| Fingerprints, GCN, ChemBERTa, graph models | `references/mol-machine-learning.md` |
| MoleculeNet, protein-ligand, virtual screening | `references/drug-discovery.md` |
| QM9, DeepQMC, materials science | `references/quantum-materials.md` |

## Installation

```bash
pip install --pre deepchem          # with TensorFlow
pip install --pre deepchem[torch]   # with PyTorch
pip install --pre deepchem[jax]     # with JAX
```

```python
import deepchem as dc
dc.__version__   # verify installation
```

## Key Submodules

| Submodule | Role |
|-----------|------|
| `dc.molnet` | MoleculeNet dataset loaders |
| `dc.models` | All model classes |
| `dc.feat` | Featurizers |
| `dc.metrics` | Evaluation metrics |
| `dc.splits` | Dataset splitters |
| `dc.data` | Dataset classes |
| `dc.trans` | Transformers (normalization, etc.) |

## Related Skills

- `rdkit-patterns` - Molecular manipulation before DeepChem ingestion
- `cheminformatics` - SMILES, molecular representations reference
