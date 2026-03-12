---
name: active-learning
description: Use when designing active learning or closed-loop molecular optimization (DMTA cycles). Covers query strategies (UCB/EI/BALD/QBC), batch DPP selection, docking oracles (Vina/Gnina), BEDROC/EF metrics, and Design-Make-Test-Analyze campaign management.
---

# Active Learning for Drug Discovery

## Purpose
Closed-loop molecular optimization: iteratively query the most informative compounds, label with assay/oracle, retrain model. Accelerates hit-to-lead and lead optimization by minimizing wet-lab experiments.

## When to Use This Skill
- Building a surrogate model to replace expensive docking/assay calls
- Running a DMTA (Design-Make-Test-Analyze) loop
- Accelerating virtual screening with AL-ranked acquisition
- Combining QSAR uncertainty with experimental prioritization

## Reference Files
Load specific references on demand:

| File | Content |
|------|---------|
| `references/al-theory.md` | Query strategies, acquisition functions, convergence, pool vs stream |
| `references/molecular-al.md` | Molecular representations, batch AL, diversity-reweighted sampling |
| `references/uncertainty-integration.md` | GP/conformal/ensemble signals → acquisition, calibration |
| `references/docking-al.md` | Surrogate docking oracle, VS acceleration, Logloss/BEDROC metrics |
| `references/dmta-loop.md` | Full DMTA cycle, stopping criteria, experiment prioritization, case studies |

## Quick Routing

**"I want to find actives with fewest assay calls"**
→ `al-theory.md` (query strategy) + `molecular-al.md` (batch AL)

**"I want to accelerate a docking campaign"**
→ `docking-al.md` (surrogate oracle)

**"I have GP/conformal uncertainty, want to plug into AL loop"**
→ `uncertainty-integration.md`

**"I'm running a real DMTA cycle with a CRO"**
→ `dmta-loop.md`

## Core Loop Pattern

```python
# Canonical active learning loop
labeled_pool = initial_dataset          # seed: 50–200 diverse cpds
unlabeled_pool = virtual_library        # 10k–1M candidates

for round in range(n_rounds):
    model.fit(labeled_pool.X, labeled_pool.y)
    scores = acquisition_fn(model, unlabeled_pool.X)  # uncertainty / EI / UCB
    batch = select_batch(unlabeled_pool, scores, k=batch_size)
    labels = oracle(batch)              # assay / docking / human expert
    labeled_pool = labeled_pool + (batch, labels)
    unlabeled_pool = unlabeled_pool - batch
```

## Key Principles
1. **Seed diversity matters** — MaxMin or clustering on initial pool prevents early bias
2. **Batch mode requires diversity** — greedy top-k collapses; use DPP or greedy submodular
3. **Calibrate before querying** — miscalibrated uncertainty → wrong queries (use MAPIE/isotonic)
4. **Track enrichment, not accuracy** — AL goal is finding actives fast, not global R²
5. **Stopping criterion** — plateau in hit rate OR budget exhausted

## Integration with ALKYL Skills
- Uncertainty estimates: `uncertainty-qsar` skill (GP, conformal, deep ensembles)
- Diversity selection: `chem_diversity.py` (MaxMin)
- Docking oracle: `docking` skill (Vina/Gnina)
- Library design: `generative-design` skill (REINVENT + AL reward)
- Property filtering: `chem_filter.py`, `chem_batch.py`
