---
name: coarse-grained
description: Use when running coarse-grained (CG) molecular dynamics with MARTINI 3. Covers protein CG with martinize2, membrane assembly with insane.py (POPC/POPE/CHOL bilayers), GROMACS CG workflows, backmapping CG to all-atom, and CG membrane analysis.
---

# Coarse-Grained Molecular Dynamics

## Purpose
Run µs–ms scale MD simulations using coarse-grained force fields.
Primary use cases: membrane self-assembly, protein-membrane interactions,
lipid nanoparticles, large conformational changes, crowding effects.

## When to Use This Skill
- Simulating lipid bilayers, vesicles, or membrane proteins
- Accessing timescales (µs–ms) beyond all-atom MD reach
- Screening protein-membrane binding or insertion
- Studying large-scale conformational changes (IDPs, domain motion)
- Building membrane systems for subsequent AA MD (backmapping)
- Coarse-grained small molecule parameterization (MARTINI)

## Reference Files

| File | Content |
|------|---------|
| `references/cg-theory.md` | CG resolution levels, mapping schemes, Boltzmann inversion, force matching, MARTINI 3 philosophy, bead types, scaling factors |
| `references/martini-proteins.md` | martinize2, elastic network (ElNeDyn), Go-MARTINI, OpenMM/GROMACS protein CG setup, common pitfalls |
| `references/martini-membranes.md` | Lipid library, insane.py membrane builder, CHARMM-GUI CG, protein-membrane embedding, lipid mixing |
| `references/cgmd-workflows.md` | GROMACS CG workflow (mdp parameters, timestep, thermostat), OpenMM CG, backmapping (backward.py), equilibration protocol |
| `references/cg-analysis.md` | MDAnalysis CG trajectories, membrane thickness/APL/order parameters, lateral diffusion, protein CG RMSD/RMSF, density profiles |

## Quick Routing

**"Set up a lipid bilayer simulation"** → `martini-membranes.md`

**"Convert my protein to MARTINI CG"** → `martini-proteins.md`

**"Run a CG simulation in GROMACS"** → `cgmd-workflows.md`

**"Backmap CG structure to all-atom"** → `cgmd-workflows.md` (backward.py section)

**"Analyze membrane properties from CG trajectory"** → `cg-analysis.md`

**"What resolution should I use?"** → `cg-theory.md`

## Key Numbers (MARTINI 3)

| Property | Value |
|----------|-------|
| Mapping ratio | ~4 heavy atoms per bead |
| Timestep (default) | 20 fs (safe: 10–30 fs) |
| Time scaling factor | ×4 (CG time ≈ 4× real time) |
| vdW cutoff | 1.1 nm |
| Electrostatics cutoff | 1.1 nm |
| Recommended thermostat | v-rescale (τ=1 ps) |
| Recommended barostat | Parrinello-Rahman (τ=12 ps) |
| Effective timestep | 80 fs (20 fs × 4 scaling) |
| Accessible timescale | µs per day (GPU) |

## Integration with ALKYL Skills
- AA structure for CG input: `homology-modeling` or `force-fields` skill
- Post-backmapping refinement: `force-fields` skill (OpenMM minimization)
- Trajectory analysis: `mdanalysis` skill (most tools work on CG trajectories)
- Membrane-protein docking: informed by CG binding mode
