# Model Training in DeepChem

## Standard 5-Step Workflow

```python
import deepchem as dc

# Step 1: Load dataset
tasks, datasets, transformers = dc.molnet.load_tox21(featurizer='ECFP', splitter='scaffold')
train, valid, test = datasets

# Step 2: Create model
model = dc.models.MultitaskClassifier(
    n_tasks=len(tasks),
    n_features=1024,
    layer_sizes=[1000]
)

# Step 3: Train
model.fit(train, nb_epoch=10)

# Step 4: Evaluate
metric = dc.metrics.Metric(dc.metrics.roc_auc_score, np.mean)
train_scores = model.evaluate(train, [metric], transformers)
test_scores  = model.evaluate(test,  [metric], transformers)
print(f"Train ROC-AUC: {train_scores}, Test ROC-AUC: {test_scores}")

# Step 5: Predict
predictions = model.predict(test)  # on full dataset
batch_preds = model.predict_on_batch(test.X[:10])  # on batch
```

## Core Model Classes

| Model | Mode | Use case |
|-------|------|----------|
| `GraphConvModel` | regression/classification | General molecular property prediction |
| `AttentiveFPModel` | regression/classification | State-of-the-art graph attention |
| `MultitaskClassifier` | classification | Multiple binary classification tasks |
| `MultitaskRegressor` | regression | Multiple regression tasks |
| `ChemBERTaModel` | regression/classification | Transfer learning from SMILES |
| `GATModel` | regression/classification | Graph Attention Networks |
| `MPNNModel` | regression/classification | Message Passing Neural Networks |

```python
# GraphConvModel
model = dc.models.GraphConvModel(n_tasks=1, mode='regression', dropout=0.2, batch_normalize=False)

# AttentiveFP (often best performer)
model = dc.models.AttentiveFPModel(n_tasks=1, mode='regression', num_layers=2, num_timesteps=2)

# MultitaskClassifier (ECFP features)
model = dc.models.MultitaskClassifier(n_tasks=12, n_features=1024, layer_sizes=[1000, 500])
```

## Metrics

```python
# Regression metrics
dc.metrics.pearson_r2_score   # R2 correlation
dc.metrics.mean_absolute_error
dc.metrics.mean_squared_error
dc.metrics.rms_score          # RMSE

# Classification metrics
dc.metrics.roc_auc_score      # ROC-AUC
dc.metrics.accuracy_score
dc.metrics.recall_score
dc.metrics.precision_score
dc.metrics.f1_score

# Usage
metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)
metric = dc.metrics.Metric(dc.metrics.roc_auc_score, np.mean)  # average over tasks

scores = model.evaluate(test_dataset, [metric], transformers)
```

## Callbacks and Training Control

```python
# Early stopping via validation
best_val = -np.inf
for epoch in range(50):
    model.fit(train, nb_epoch=1)
    val_score = model.evaluate(valid, [metric], transformers)
    if val_score['pearson_r2_score'] > best_val:
        best_val = val_score['pearson_r2_score']
        model.save_checkpoint()

model.restore()  # restore best checkpoint
```

## Hyperparameter Optimization

```python
import deepchem as dc

def model_builder(model_dir):
    return dc.models.GraphConvModel(
        n_tasks=1,
        mode='regression',
        model_dir=model_dir
    )

params_dict = {
    'nb_epoch': [50, 100],
    'learning_rate': [0.001, 0.0001],
    'dropout': [0.0, 0.2, 0.5],
}

optimizer = dc.hyper.GridHyperparamOpt(model_builder)
best_model, best_hyperparams, all_results = optimizer.hyperparam_search(
    params_dict, train, valid, transformers,
    metric=dc.metrics.Metric(dc.metrics.pearson_r2_score)
)
```

## Multitask Learning

```python
# Multitask: predict multiple properties simultaneously
# Particularly powerful for drug discovery (Tox21 has 12 tasks)
tasks, datasets, transformers = dc.molnet.load_tox21(featurizer='ECFP')
train, valid, test = datasets

model = dc.models.MultitaskClassifier(
    n_tasks=len(tasks),   # 12 for Tox21
    n_features=1024,
    layer_sizes=[1000, 500]
)
model.fit(train, nb_epoch=10)

# Evaluate per-task
metric = dc.metrics.Metric(dc.metrics.roc_auc_score)
scores = model.evaluate(test, [metric], transformers)
# Returns dict with per-task scores
```

## High-Fidelity Models from Experimental Data

```python
# When you have your own experimental data
import pandas as pd

df = pd.read_csv('my_assay.csv')  # columns: smiles, activity

loader = dc.data.CSVLoader(
    tasks=['activity'],
    feature_field='smiles',
    featurizer=dc.feat.MolGraphConvFeaturizer()
)
dataset = loader.create_dataset('my_assay.csv')

splitter = dc.splits.RandomSplitter()
train, valid, test = splitter.train_valid_test_split(dataset)
```

## PyTorch Lightning Integration

```python
import pytorch_lightning as pl
import deepchem as dc

# DeepChem models can be wrapped for Lightning training
model = dc.models.AttentiveFPModel(n_tasks=1, mode='regression', learning_rate=1e-3)

# Use standard Lightning trainer features: logging, checkpointing, etc.
```
