---
name: fbdd
description: Use for fragment-based drug design (FBDD): Rule of 3 filtering, ligand efficiency metrics (LE/LLE/BEI/LELP), fragment library design, fragment docking (Vina/Gnina), fragment growing/linking/merging strategies, and Abad-Zapatero efficiency plots.
---

# Fragment-Based Drug Design (FBDD)

## Purpose
Design and analyze fragment libraries, compute ligand efficiency metrics,
perform fragment docking, and execute fragment-to-lead elaboration
(growing, linking, merging) with computational support.

## When to Use This Skill
- Building or filtering a fragment library
- Computing LE/LLE/LLEAT efficiency metrics
- Docking fragments into a target (weak binding, requires special settings)
- Growing a fragment hit toward lead-like compounds
- Merging two fragment hits sharing a common substructure
- Analyzing X-ray fragment screening data

## Reference Files

| File | Content |
|------|---------|
| `references/fbdd-theory.md` | Fragment rules (Rule of 3), LE/LLE/LLEAT/BEI/SEI, Hann complexity model, fragment-to-lead strategies (grow/link/merge), success stories |
| `references/fragment-library.md` | Library design: RDKit filters (Ro3/PAINS/flatness/rigidity), 3D sp3 character, commercial sources, diversity selection, quality checks |
| `references/fragment-docking.md` | Low-MW docking pitfalls, Vina fragment settings, ROCS shape screening, Smina fragment mode, pose clustering, hotspot validation |
| `references/fragment-growing.md` | Scaffold growing (R-group enumeration, MMPA vectors), fragment merging (MCS-based), FBDD-aware REINVENT, SynthesizabilityOracle, elaboration scoring |
| `references/efficiency-metrics.md` | LE/LLE/LLEAT/BEI/SEI formulas, efficiency evolution plots, Abad-Zapatero plots, LELP, GE (group efficiency), metric-driven SAR |

## Quick Routing

**"Build a fragment library"** → `fragment-library.md`

**"Dock fragments into my target"** → `fragment-docking.md`

**"I have a fragment hit, want to grow it"** → `fragment-growing.md`

**"Track efficiency as I optimize"** → `efficiency-metrics.md`

**"What makes a good fragment?"** → `fbdd-theory.md`

## Core Concept: Rule of 3

| Property | Fragment (Ro3) | Lead-like | Drug-like (Ro5) |
|----------|---------------|-----------|-----------------|
| MW | ≤ 300 Da | ≤ 400 Da | ≤ 500 Da |
| cLogP | ≤ 3 | ≤ 4 | ≤ 5 |
| HBD | ≤ 3 | ≤ 4 | ≤ 5 |
| HBA | ≤ 3 | ≤ 8 | ≤ 10 |
| PSA | — | ≤ 120 Å² | — |
| Rotatable bonds | ≤ 3 | ≤ 7 | ≤ 10 |

## Minimal LE Calculation

```python
def ligand_efficiency(pIC50, n_heavy_atoms):
    """LE = ΔG / HAC ≈ 1.37 * pIC50 / HAC (kcal/mol per heavy atom)"""
    return 1.37 * pIC50 / n_heavy_atoms

# Good fragment: LE ≥ 0.3 kcal/mol/HA
# Drug-like optimum: LE ≥ 0.3 (maintain or improve during optimization)
```

## Integration with ALKYL Skills
- Fragment docking: `docking` skill (Vina/Gnina, lower exhaustiveness ok)
- Fragment diversity: `chem_diversity.py` (MaxMin)
- Fragment filtering: `chem_filter.py`, `chem_batch.py`
- Growing enumeration: `chem_react.py` (ReactionFromSmarts)
- MMPA elaboration: `mmpa` skill
- Generative growing: `generative-design` skill (REINVENT with fragment constraint)
- 3D visualization: `py3Dmol` skill
