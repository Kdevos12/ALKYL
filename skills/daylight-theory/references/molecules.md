# Daylight Theory — Molecular Representation

## Molecular Graph Model

Molecules are represented as **labeled graphs**:
- **Nodes** = atoms (atomic number, charge, isotope, H count, chirality)
- **Edges** = bonds (single, double, triple, aromatic)

Two representations coexist:
- **Hydrogen-suppressed graph** (Hs as atom properties) — default in SMILES
- **Hydrogen-complete graph** (Hs as explicit atoms) — for atom-level tracking

Bond notation: either **Kekulé** (alternating single/double) or **aromatic** — both encode the same ring.

## Explicit vs. Derived Properties

| Type | Meaning | Examples |
|------|---------|---------|
| Explicit | Stored directly in the graph | Atomic number, formal charge, isotope |
| Derived | Computed from the graph | Aromaticity, SSSR, symmetry, canonical numbering |

Derived properties are **not stored** — they are recomputed on demand.

## Ring Analysis

### Ring Bond Detection
A **ring bond** is one whose removal does not disconnect the graph (biconnected property). Atoms connected only via ring bonds are **ring atoms**.

### SSSR — Smallest Set of Smallest Rings
The SSSR is computationally expensive (graph-theoretic). Key facts:
- Not unique: multiple valid SSSRs can describe the same ring system
- Naphthalene has 3 possible paths → 6 valid SSSRs (two 6-membered rings + one 10-membered ring)
- For cheminformatics purposes, consistency matters more than uniqueness

## Bond Order vs. Bond Type

| Property | Nature | Values |
|----------|--------|--------|
| Bond order | Formal (stored) | 1, 2, 3 |
| Bond type | Derived | single, double, triple, **aromatic** |

Aromatic bond type is derived, not stored.

## Aromaticity

**Rule:** A ring is aromatic if and only if:
1. All atoms in the ring are sp² hybridized
2. The number of "excess" p-electrons satisfies **Hückel's 4N+2** rule (N = 0, 1, 2, ...)

Examples:
| Molecule | Electrons | Aromatic? |
|----------|-----------|-----------|
| Benzene | 6 (N=1) | ✓ |
| Naphthalene | 10 (N=2) | ✓ |
| Cyclopentadienyl anion | 6 (N=1) | ✓ |
| Cyclooctatetraene | 8 | ✗ |

Exocyclic double bonds do **not** break aromaticity.

Aromatic atoms: C, N, O, P, S, As, Se (and `*` wildcard in SMARTS).

Note: This is an operational definition for molecular representation, not a prediction of reactivity or NMR.

## Symmetry and Canonical Labeling

**Symmetry detection** identifies equivalent atoms (2D rotational symmetry), enabling:
- Generation of canonical atom numbering
- Chirality classification (chiral vs. pseudo-chiral vs. degenerate-chiral)
- Elimination of redundant computations

**Canonical labeling** assigns a history-independent atom ordering, necessary for:
- Generating a unique SMILES string
- Molecular database hashing and lookup
- Comparison of independently drawn molecules

Without canonical labeling, the same molecule drawn differently yields different strings.

## Chirality

Chirality is both:
- **Explicit input**: specified by the user (`@`, `@@` in SMILES)
- **Derived property**: computed from the graph topology

The system converts any chiral specification to a canonical form during unique SMILES generation. Stereodescriptors (`@`/`@@`) are local — they describe the view from one direction.

## Reaction Representation

Reactions are sets of molecules with assigned roles:

| Role | Meaning |
|------|---------|
| Reactant | Contributes atoms to products (not enforced) |
| Agent | Does not contribute/accept atoms (catalysts, solvents) |
| Product | Final state; atoms should come from reactants (not enforced) |

**Atom maps** track which reactant atoms become which product atoms. They are optional and need not be complete.

Data external to the graph (stoichiometry, yield, conditions, equilibrium) is stored separately — not encoded in the molecular graph.

## Depiction

2D depiction is generated **algorithmically** from scratch using priority ordering:
1. **Correctness** (graph accurately shown)
2. **Comprehensibility** (chemist-readable)
3. **Appearance** (aesthetic quality, lowest priority)

This means any molecule can be depicted without prior visual history.

## Valence Model

The Daylight valence model knows "normal valences of organic compounds":

| Atom | Normal valences |
|------|----------------|
| B | 3 |
| C | 4 |
| N | 3, 5 |
| O | 2 |
| P | 3, 5 |
| S | 2, 4, 6 |
| Halogens | 1 |

Hydrogens are filled automatically by electron counting. Violations trigger warnings but are not rejected.
