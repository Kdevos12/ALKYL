---
name: binding-kinetics
description: Use when analyzing or predicting drug-target binding kinetics: kon, koff, KD, residence time, SPR data fitting (Langmuir/two-state), ITC thermodynamics, tau-RAMD and funnel metadynamics for unbinding, or kinetic QSAR models.
---

# Binding Kinetics

## Purpose
Analyze, predict, and optimize drug-target binding kinetics:
on-rates (kon), off-rates (koff), residence time (RT = 1/koff),
thermodynamic signatures (ΔH/ΔS), and structure-kinetics relationships (SKR).

## When to Use This Skill
- Analyzing SPR sensorgrams (Biacore/Sierra)
- Fitting ITC thermograms for ΔH/ΔS/ΔG
- Computing residence time from MD simulations (τRAMD, metadynamics)
- Building QSAR models for koff/kon
- Interpreting kinetic selectivity vs equilibrium selectivity
- Prioritizing compounds by residence time, not just KD

## Reference Files

| File | Content |
|------|---------|
| `references/kinetics-theory.md` | kon/koff/KD/RT definitions, kinetic selectivity, two-state binding, conformational selection vs induced fit, thermodynamic signatures |
| `references/spr-analysis.md` | SPR sensorgrams, 1:1 Langmuir fitting, two-state model, Rmax/Rtheor, bulk correction, Biacore data parsing, Python fitting |
| `references/itc-analysis.md` | ITC thermogram integration, n/KD/ΔH/ΔS/ΔG fitting, SEDPHAT equivalents in Python, van't Hoff, enthalpy-entropy compensation |
| `references/residence-time-md.md` | τRAMD (random acceleration MD), funnel metadynamics, WTmetaD koff estimation, HTMD τRAMD Python, unbinding pathway analysis |
| `references/kinetic-qsar.md` | Structure-kinetics relationships (SKR), features for koff/kon models, kinetic maps (LE vs kinetic efficiency), koff cliff detection |

## Quick Routing

**"Fit my SPR data"** → `spr-analysis.md`

**"Fit my ITC experiment"** → `itc-analysis.md`

**"Compute residence time from MD"** → `residence-time-md.md`

**"Build a model to predict koff"** → `kinetic-qsar.md`

**"Why does my drug work despite poor KD?"** → `kinetics-theory.md`

## Key Relationships

```python
# Core kinetic relationships
KD = koff / kon                    # M (equilibrium dissociation constant)
pKD = -log10(KD)                   # analogous to pIC50
RT = 1 / koff                      # seconds (residence time)
t_half = ln(2) / koff              # seconds (half-life of complex)

# Thermodynamics
ΔG = RT_gas * ln(KD)               # kcal/mol (RT_gas = 0.592 at 298K)
ΔG = ΔH - T*ΔS                    # enthalpy-entropy decomposition

# Typical drug ranges
# kon:  10^4 – 10^7 M^-1 s^-1
# koff: 10^-5 – 10^-1 s^-1
# KD:   nM – µM
# RT:   10 s – 10^5 s (hours)
```

## Integration with ALKYL Skills
- Docking poses for kinetic path analysis: `docking` skill
- MD trajectories for τRAMD: `force-fields` skill + MDAnalysis
- SKR feature computation: `chem_props.py`, `chem_analyze.py`
- MMPA for koff SAR: `mmpa` skill
- Uncertainty in kinetic predictions: `uncertainty-qsar` skill
