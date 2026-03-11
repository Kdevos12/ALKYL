# ALKYL

```
 ░▒▓██████▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓████████▓▒░▒▓█▓▒░      ░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓████████▓▒░
```

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![Claude Code](https://img.shields.io/badge/Claude_Code-plugin-blueviolet)](https://claude.ai/code)
[![Skills](https://img.shields.io/badge/skills-23-green)](#skills)
[![Scripts](https://img.shields.io/badge/scripts-23-green)](#scripts)
[![Tests](https://img.shields.io/badge/tests-162-brightgreen)](tests/)

**A [Claude Code](https://claude.ai/code) plugin for computational chemistry and drug discovery.**

ALKYL transforms Claude Code into a specialized computational chemistry assistant. Install once, then work naturally: RDKit cheminformatics, molecular docking, MD simulations, quantum chemistry, free energy calculations, and ML-guided drug design — all through plain conversation, with no wrapper CLI.

> Designed for computational chemists, medicinal chemists, and drug discovery researchers who use Claude Code as their daily driver.

---

## Installation

### Prerequisites

- [Claude Code](https://claude.ai/code) installed and working (`claude --version`)
- Git
- Python ≥ 3.9 (for scripts and tests)

### Step 1 — Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/alkyl
cd alkyl
```

### Step 2 — Install ALKYL

```bash
bash alkyl.sh install
```

This injects a chemistry-specialized context block into `~/.claude/CLAUDE.md` — the global configuration file Claude Code reads at every session start. No daemon, no wrapper, no separate command.

### Step 3 — Set up the Python environment (for scripts)

```bash
bash alkyl.sh venv
```

Creates `.venv/` with RDKit and pytest. Required only if you want to use the standalone scripts or run tests.

### Step 4 — Verify

```bash
bash alkyl.sh status
```

Expected output:
```
✓ Installed — /home/<user>/.claude/CLAUDE.md
  Block size: ~150 lines
  Scripts: /path/to/alkyl/scripts
  Skills:  23 loaded
```

### Step 5 — Open a new Claude Code session

```bash
claude
```

ALKYL is now active. Try:
```
"Compute QED and SA score for aspirin: CC(=O)Oc1ccccc1C(=O)O"
"Set up a virtual screening run against PDB:3HTB"
"Write an ORCA B3LYP-D3BJ/def2-TZVP geometry optimization input"
```

---

## How it works

ALKYL injects a chemistry-specialized context block into `~/.claude/CLAUDE.md`:

```
bash alkyl.sh install
  → appends ALKYL block to ~/.claude/CLAUDE.md
  → idempotent: re-running replaces the old block cleanly

bash alkyl.sh uninstall
  → removes the block via <!-- ALKYL-START/END --> markers
```

**Skills** are Markdown reference files in `skills/`. They are loaded on demand in Claude Code sessions using the built-in `/skill` mechanism — only the relevant skill is loaded, keeping context lean.

### All commands

```bash
bash alkyl.sh install              # install ALKYL context
bash alkyl.sh venv                 # create .venv with RDKit
bash alkyl.sh status               # show installation status and MCP keys
bash alkyl.sh repair               # force re-inject config (fixes corruption)
bash alkyl.sh uninstall            # remove ALKYL from ~/.claude/CLAUDE.md
bash alkyl.sh setup-key perplexity <KEY>   # configure Perplexity API (optional)
```

---

## What you can ask

Once installed, ALKYL responds to natural chemistry requests:

```
"Compute QED, cLogP, and SA score for this SMILES: CC(=O)Oc1ccccc1C(=O)O"
"Set up an AutoDock Vina virtual screening run against PDB:3HTB"
"Write an ORCA input for B3LYP-D3BJ/def2-TZVP geometry optimization"
"Run a MARTINI 3 membrane simulation with POPC bilayer"
"Design a focused library around this fragment using REINVENT 4"
"Estimate RBFE for these two congeneric ligands using OpenMMTools HREX"
"Flag all hERG and PAINS alerts in my SDF library"
"Explain the SN2 mechanism for this substrate using EASE"
```

---

## Skills

23 domain-specific skills, organized by workflow stage. Each skill is a reference file with practical code patterns and theoretical context — loaded only when needed.

### Cheminformatics & molecular tools

| Skill | Description |
|-------|-------------|
| `rdkit` | Molecule I/O, descriptors (MW, cLogP, TPSA, QED), Morgan/MACCS fingerprints, 2D/3D conformer generation, substructure search, SMARTS reactions, SVG/PNG visualization |
| `openbabel` | Format conversion (146 formats), 3D structure generation (MMFF94/UFF/GAFF), conformer search, protonation at pH, FP2/FP3/FP4/MACCS fingerprints, RDKit interoperability |
| `daylight-theory` | Complete SMILES spec, SMARTS query language (all primitives, recursive SMARTS, reaction queries), SMIRKS transforms, path fingerprints, similarity metrics (Tanimoto/Dice/Tversky/Cosine + 15 variants) |
| `chem-brainstorm` | Workflow guide: classify → audit data → map tools → generate directions → sanity checks → literature. 4 rigid protocols (molecule evaluation, SAR hypothesis, reaction design, pipeline). Integrates ALKYL scripts + MCP tools (ChEMBL, OpenTargets, bioRxiv, ClinicalTrials) |

### Molecular dynamics & structure

| Skill | Description |
|-------|-------------|
| `ase` | Atoms objects, geometry optimization (BFGS/LBFGS/FIRE), NVE/NVT/NPT MD (Langevin, Berendsen), NEB/AutoNEB transition states, vibrational analysis, thermochemistry, ORCA/xTB/GPAW/LAMMPS calculators |
| `mdanalysis` | Universe/AtomGroup selection language, RMSD/RMSF/alignment, contact analysis (Q-value), H-bond analysis, Ramachandran/DSSP, PCA free energy landscapes, RDF, MSD diffusion, protein-ligand workflow |
| `force-fields` | AMBER/CHARMM/OPLS-AA/SMIRNOFF families, OpenMM simulation setup (LangevinMiddleIntegrator, NPT barostat, DCD/XTC reporters), OpenFF Sage 2.2, GAFF2 parameterization (antechamber/acpype), AM1-BCC/RESP charges, HMR |
| `coarse-grained` | MARTINI 3 CG simulations: protein CG with martinize2 (ElNeDyn, Go-MARTINI), membrane assembly with insane.py (POPC/POPE/POPS/CHOL/PIP2, asymmetric bilayers), GROMACS workflows, backward.py backmapping, membrane analysis (thickness, APL, Scd, lateral diffusion) |

### Quantum chemistry

| Skill | Description |
|-------|-------------|
| `qm-dft` | DFT functional/basis selection (Jacob's ladder, D3BJ dispersion), ORCA 6.0 (Opt/Freq/TS/TD-DFT/NMR/DLPNO/solvation), xTB/GFN2 (CLI, tblite API, CREST, pKa), PySCF (HF/DFT/MP2/CCSD, GIAO NMR, ESP/CHELPG), standard workflows (opt→freq→SP, barriers, UV-Vis) |
| `organic-mechanisms` | Polar mechanism reasoning via EASE framework (Electrophile/Acid-Base/Sterics/Electron-Flow), SN1/SN2/E1/E2 decision tree, Zaitsev/Hofmann selectivity, HSAB (1,2 vs 1,4), arrow-pushing rules, retrosynthesis (disconnections, synthons, FGI, C–C/C–X toolbox) |

### Drug discovery & docking

| Skill | Description |
|-------|-------------|
| `docking` | Receptor preparation (pdbfixer/propka3), AutoDock Vina Python API + CLI, Gnina CNN rescoring, meeko PDBQT prep, batch parallel docking, ProLIF interaction fingerprints, RMSD pose clustering, ensemble docking on MD snapshots |
| `homology-modeling` | Template search (HHblits/BLAST), BLOSUM62/PIR alignment, MODELLER 10 (automodel/loopmodel/DOPE ranking), AlphaFold2 via ColabFold CLI, ESMFold API, structure quality (pLDDT, Ramachandran, MolProbity), structure prep (HIS tautomers, disulfides, capping) |
| `fbdd` | Rule of 3 filters, Ligand Efficiency metrics (LE/LLE/LLEAT/BEI/SEI/GE/LELP), fragment library design (Ro3+PAINS+reactive+Fsp3), fragment docking (high exhaustiveness, Gnina, RMSD clustering), growing/linking/merging (R-group enumeration, MCS, REINVENT scaffold constraint), Abad-Zapatero plot |
| `free-energy` | Thermodynamic cycles (RBFE/ABFE), FEP/TI/BAR/MBAR estimators, OpenMMTools AlchemicalFactory/MultiStateSampler/HREX, RBFE network design (LOMAP, perses, openfe), ABFE double-decoupling with Boresch restraints, pymbar (overlap matrix, convergence, autocorrelation) |
| `binding-kinetics` | kon/koff/KD/residence time theory (Copeland framework), two-state/induced-fit/conformational-selection models, SPR fitting (Langmuir, Biacore CSV), ITC analysis (Wiseman isotherm, ΔG/ΔH/ΔS/ΔCp), τRAMD (HTMD + GROMACS), funnel metadynamics (PLUMED), kinetic QSAR (RF/GP koff models) |
| `pharmacophore` | Feature types (HBD/HBA/AR/HYD/POS/NEG), FDEF format, RDKit Pharm2D Gobbi fingerprints, Pharm3D 3D matching, structure-based pharmacophore from ProLIF/PLIP interactions, ligand-based alignment (O3A, DBSCAN), full VS pipeline (conformers → scoring → exclusion volumes → EF/ROC) |

### Molecular design & ML

| Skill | Description |
|-------|-------------|
| `generative-design` | SELFIES always-valid grammar, SMILES LMs (LSTM/GPT2/ChemGPT), REINVENT 4 RL (QED/SA/docking oracle/custom scoring), JT-VAE latent space Bayesian optimization (botorch), structure-based generation (DiffSBDD, TargetDiff, DiffLinker), MOSES/GuacaMol evaluation |
| `mmpa` | Matched Molecular Pair Analysis: Hussain-Rea fragmentation, SMIRKS transforms, mmpdb 4 CLI workflow (fragment→index→loadprops→transform→analyze), activity cliff detection (SALI), bioisostere table, focused library generation, REINVENT/docking integration |
| `uncertainty-qsar` | Conformal prediction (MAPIE split/CV+, coverage guarantee, Mondrian), GP with GPyTorch TanimotoKernel, MC Dropout (T=50), deep ensembles (M=5), heteroscedastic head, Laplace approximation, applicability domain (kNN Tanimoto, Williams plot, Mahalanobis), OECD Principle 3 |
| `active-learning` | Query strategies (UCB/EI/BALD/QBC/Core-Set), batch DPP/cluster-then-rank, docking oracle (Vina/Gnina, ~50× screening speedup), BEDROC/EF evaluation, DMTA cycle management (batch composition, stopping criteria, round reports) |

### Visualization & utilities

| Skill | Description |
|-------|-------------|
| `py3Dmol` | 3Dmol.js visualization: PDB/SDF/SMILES loading, cartoon/stick/sphere/surface styles (SES/SAS/VDW), selection language, color schemes (spectrum/b-factor), docking pose batch viewer (ipywidgets), pharmacophore overlay, conformer animation, PNG/HTML export, NGLview for MD |
| `lit-rescue` | Literature search of last resort when hallucination risk is >20%: Perplexity→bioRxiv→PubMed waterfall, 7 query types (METHOD/PARAM/BUG/THEORY/PROTOCOL/BENCHMARK/DOMAIN), confidence reporting (★★★ to ☆☆☆), mandatory negative result block when no source found |

---

## Scripts

23 standalone Python scripts in `scripts/`. Each requires only RDKit (+ stdlib). Run with any Python ≥ 3.9 environment with RDKit.

| Script | Description |
|--------|-------------|
| `chem_convert.py` | Convert molecules between SMILES, SDF, InChI, InChIKey, and SVG. Batch-capable. |
| `chem_props.py` | MW, cLogP, TPSA, HBD, HBA, RotBonds, QED. Lipinski Ro5 + PAINS alerts. Morgan (ECFP4) and MACCS fingerprints. |
| `chem_fetch.py` | Fetch from PubChem (name/CID/InChI) and ChEMBL (name/ID). SMILES, synonyms, properties. stdlib only. |
| `chem_3d.py` | 3D conformer generation (ETKDGv3) + MMFF94/UFF minimization. Outputs SDF. |
| `chem_qm.py` | ORCA/Gaussian input from SMILES (auto 3D embed). Parse ORCA output: energy, frequencies, thermochemistry, IR. |
| `chem_batch.py` | Batch-process SDF/SMI/CSV: descriptors, Lipinski Ro5, PAINS. `--skip-invalid` for robust pipelines. |
| `chem_search.py` | Substructure (SMILES/SMARTS), Tanimoto similarity, or exact match search against SDF/SMI libraries. |
| `chem_standardize.py` | Desalt (largest fragment), neutralize charges, canonicalize SMILES via RDKit MolStandardize. |
| `chem_analyze.py` | Single-molecule deep analysis: formula, 16 functional groups, ring systems, stereocenters, QED, SA score, Bertz complexity. |
| `chem_scaffold.py` | Murcko scaffold, generic scaffold, BRICS fragments. |
| `chem_compare.py` | Two-molecule comparison: MCS (rdFMCS), Tanimoto, Δ properties (MW, cLogP, TPSA, HBD, HBA). |
| `chem_filter.py` | Drug-likeness filters: Lipinski Ro5, Veber, Egan, Ghose, PAINS. |
| `chem_react.py` | Apply SMARTS reaction transforms (RunReactants). Deduplicate and sanitize products. |
| `chem_tautomers.py` | Enumerate tautomers (TautomerEnumerator). Returns canonical + full list with counts. |
| `chem_enum.py` | Enumerate stereoisomers (unique=True). Configurable cap on max_isomers. |
| `chem_pka.py` | SMARTS-based pKa estimation, Henderson-Hasselbalch pH-speciation, dominant protonation state at target pH. |
| `chem_metabolism.py` | CYP450 soft spot prediction: 12 SMARTS rules, five isoforms (CYP3A4/2D6/2C9/1A2/UGT-SULT). |
| `chem_diversity.py` | MaxMin diversity selection (O(k·n)). Morgan (ECFP4) or MACCS. Handles k ≥ library size. |
| `chem_cluster.py` | Butina/Taylor-Butina clustering by Tanimoto distance. Returns cluster IDs, centroids, members. |
| `chem_rgroup.py` | R-group decomposition around a SMARTS core. R1/R2/... table + unmatched count (RGroupDecompose). |
| `chem_admet.py` | Heuristic ADMET: ESOL aqueous solubility (Delaney 2004), BBB score, hERG SMARTS alerts, P-gp substrate, PPB estimate. |
| `chem_highlight.py` | SVG/PNG with SMARTS-highlighted substructure. Stdout = SVG; `--out` = SVG or PNG. |
| `chem_lead.py` | Ligand efficiency metrics (LE/LLE/BEI/LELP) from activity CSV. Tracks evolution across optimization rounds. |

---

## Optional: API keys

ALKYL ships four MCP servers out of the box with no API key required: **bioRxiv**, **ChEMBL**, **ClinicalTrials.gov**, and **PubMed**. They are active immediately after install.

### Perplexity (optional — grounded web search)

For real-time literature search in the `lit-rescue` skill:

```bash
bash alkyl.sh setup-key perplexity pplx-YOUR_KEY_HERE
```

Get a key at [perplexity.ai/settings/api](https://www.perplexity.ai/settings/api). Adds `@perplexity-ai/mcp-server` to your Claude Code MCP settings.

---

## Tests

```bash
# Unit tests (no network)
.venv/bin/python -m pytest tests/ -m "not network" -v

# All tests including network calls
.venv/bin/python -m pytest tests/ -v
```

---

## Project structure

```
alkyl/
├── alkyl.sh                # main management script (install/venv/status/repair/setup-key)
├── install.sh              # shim → alkyl.sh install
├── uninstall.sh            # shim → alkyl.sh uninstall
├── config/
│   └── CLAUDE.md           # ALKYL identity, behavior, and full skill index
├── scripts/
│   ├── chem_convert.py     # format conversion
│   ├── chem_props.py       # molecular properties and fingerprints
│   ├── chem_fetch.py       # PubChem and ChEMBL retrieval
│   ├── chem_3d.py          # ETKDGv3 conformer generation
│   ├── chem_qm.py          # ORCA/Gaussian input + output parsing
│   ├── chem_batch.py       # batch processing
│   ├── chem_search.py      # substructure, similarity, exact search
│   ├── chem_standardize.py # desalting, neutralization
│   ├── chem_analyze.py     # single-molecule deep analysis
│   ├── chem_scaffold.py    # Murcko scaffold and BRICS
│   ├── chem_compare.py     # MCS and property delta
│   ├── chem_filter.py      # drug-likeness filters
│   ├── chem_react.py       # SMARTS reaction application
│   ├── chem_tautomers.py   # tautomer enumeration
│   ├── chem_enum.py        # stereoisomer enumeration
│   ├── chem_pka.py         # pKa estimation and protonation state
│   ├── chem_metabolism.py  # CYP450 soft spot prediction
│   ├── chem_diversity.py   # MaxMin diversity selection
│   ├── chem_cluster.py     # Butina clustering
│   ├── chem_rgroup.py      # R-group decomposition
│   ├── chem_admet.py       # ADMET heuristics
│   ├── chem_highlight.py   # SMARTS-highlighted SVG/PNG
│   └── chem_lead.py        # ligand efficiency metrics per round
└── skills/
    ├── rdkit/              # RDKit cheminformatics
    ├── ase/                # Atomic Simulation Environment
    ├── mdanalysis/         # MD trajectory analysis
    ├── openbabel/          # format conversion and filtering
    ├── deepchem/           # molecular machine learning
    ├── docking/            # virtual screening and docking
    ├── force-fields/       # AMBER/OpenMM/OpenFF/GAFF2
    ├── qm-dft/             # ORCA/xTB/PySCF quantum chemistry
    ├── homology-modeling/  # MODELLER/ColabFold/ESMFold
    ├── free-energy/        # FEP/MBAR/RBFE/ABFE
    ├── pharmacophore/      # pharmacophore modeling and VS
    ├── generative-design/  # de novo molecular generation
    ├── mmpa/               # matched molecular pair analysis
    ├── uncertainty-qsar/   # conformal prediction and GP uncertainty
    ├── active-learning/    # DMTA loop and active screening
    ├── py3Dmol/            # interactive 3D visualization
    ├── coarse-grained/     # MARTINI 3 and membrane simulations
    ├── binding-kinetics/   # SPR/ITC/τRAMD/kinetic QSAR
    ├── fbdd/               # fragment-based drug design
    ├── chem-brainstorm/    # workflow brainstorming guide
    ├── daylight-theory/    # SMILES/SMARTS/SMIRKS/fingerprints theory
    ├── lit-rescue/         # literature search of last resort
    └── organic-mechanisms/ # EASE framework for polar organic mechanisms
```

---

## Requirements

- [Claude Code](https://claude.ai/code) (required)
- Python ≥ 3.9 with RDKit (for scripts and tests — created by `bash alkyl.sh venv`)
- Optional per workflow: ORCA, xTB, GROMACS, OpenMM, MODELLER, AutoDock Vina

---

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgments

The skills in this repository draw on and are informed by the following works and their authors:

- **Daylight Theory Manual** — Daylight Chemical Information Systems (SMILES, SMARTS, SMIRKS, fingerprints)
- **RDKit documentation** — Greg Landrum and RDKit contributors
- **ASE documentation** — Ask Hjorth Larsen, Jens Jørgen Mortensen, and ASE contributors
- **MDAnalysis documentation** — Oliver Beckstein, Richard Gowers, and MDAnalysis contributors
- **MARTINI force field** — Siewert-Jan Marrink, Xavier Periole, D. Peter Tieleman, and CGMD community
- **OpenFF Sage / SMIRNOFF** — Open Force Field Initiative contributors
- **REINVENT** — AstraZeneca Molecular AI team
- **mmpdb** — Andrew Dalke and contributors
- **AlphaFold / ColabFold** — DeepMind, Sergey Ovchinnikov, Martin Steinegger
- **ORCA** — Frank Neese and the ORCA development team
- **EASE organic mechanism framework** — AceOrganicChem.com *Ace Organic Chemistry Mechanisms with E.A.S.E.* (2013); Clayden *Organic Chemistry* (Oxford); March *Advanced Organic Chemistry* (Wiley)
- **Copeland binding kinetics framework** — Robert A. Copeland (*Evaluation of Enzyme Inhibitors in Drug Discovery*, Wiley)
- **Hussain-Rea fragmentation** — Jameed Hussain, Ceara Rea (J. Chem. Inf. Model., 2010)
- **Haussler Tanimoto kernel** — David Haussler (1999)
- **Conformal prediction** — Vladimir Vovk, Alexander Gammerman, Glenn Shafer (*Algorithmic Learning in a Random World*, Springer)
- All open-source tool authors and scientific communities whose work these skills build upon
