---
name: py3Dmol
description: Use when creating 3D molecular visualizations with py3Dmol (3Dmol.js). Covers PDB/SDF/SMILES loading, cartoon/stick/sphere styles, protein-ligand complexes, docking pose viewers, selection language, PNG/HTML export, and NGLview for MD trajectories.
---

# py3Dmol — Molecular Visualization

## Purpose
Interactive 3D molecular visualization in Jupyter notebooks and scripts.
Wraps 3Dmol.js (WebGL). Used for docking pose inspection, protein-ligand
complexes, conformer overlays, trajectory snapshots.

## When to Use This Skill
- Visualizing docking poses from Vina/Gnina
- Inspecting protein-ligand binding pockets
- Displaying conformer ensembles
- Annotating pharmacophore features on 3D structures
- Quick structure QC after homology modeling or MD prep

## Reference Files

| File | Content |
|------|---------|
| `references/basics.md` | Installation, view creation, loading PDB/SDF/SMILES, stick/sphere/cartoon/surface basics |
| `references/protein-ligand.md` | Protein+ligand display, binding pocket zoom, dual-structure overlay, docking pose batch |
| `references/selections-styles.md` | Selection language (chain/resi/resn/atom), color schemes, surfaces, labels, transparency |
| `references/jupyter-patterns.md` | Jupyter embed, ipywidgets sliders, NGLview alternative, saving PNG, RDKit interop |

## Quick Routing

**"Show me a docking pose"** → `protein-ligand.md`

**"Show all conformers overlaid"** → `jupyter-patterns.md` (animation loop)

**"Highlight binding pocket / surface"** → `selections-styles.md`

**"I just need a quick look at a molecule"** → `basics.md`

## Minimal Pattern

```python
import py3Dmol

view = py3Dmol.view(width=800, height=500)
view.addModel(open('complex.pdb').read(), 'pdb')
view.setStyle({'cartoon': {'color': 'spectrum'}})   # protein
view.setStyle({'resn': 'LIG'}, {'stick': {'colorscheme': 'greenCarbon'}})
view.zoomTo({'resn': 'LIG'})
view.show()
```

## Key Facts
- `py3Dmol` renders via 3Dmol.js in Jupyter — requires a running notebook kernel
- For non-Jupyter contexts: use `view.png()` → base64 PNG, or `view.write_html()`
- NGLview is an alternative with better trajectory support (use for MD)
- `setStyle` is cumulative by default; use `setStyle({}, {})` to reset all
- Ligand residue name varies: 'LIG', 'UNL', 'MOL' — check with grep before scripting
