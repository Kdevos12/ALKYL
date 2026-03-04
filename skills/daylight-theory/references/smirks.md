# Daylight Theory — SMIRKS Reaction Transforms

SMIRKS is a **reaction transform language** combining SMILES (for products) and SMARTS (for reactant patterns). It encodes **generic reactions** — any set of reactions sharing the same atom and bond changes, regardless of molecular substrate.

> A transform is "a generic reaction within the Daylight system."

---

## 1. Core Concept

A SMIRKS transform captures:
1. **Actual molecular changes** (bonds formed/broken, charges changed)
2. **Activating/deactivating context** (expressed via SMARTS)

Example — Sn2 displacement at a carbon bearing a leaving group:
```
[C:1][Br:2].[I-:3]>>[C:1][I:3].[Br-:2]
```

Bond changes encoded:
- C–Br: single bond → no bond
- C–I: no bond → single bond
- Br: charge 0 → –1
- I: charge –1 → 0

---

## 2. Syntax

**Grammar:**
```
transform    : reactant '>' agent '>' product
             | reactant '>>' product
```

Where:
- `reactant` and `product` are SMARTS patterns with atom maps
- `agent` is optional (catalysts, solvents — not transformed)

**Atom map syntax:** `[expr:N]` where N is a positive integer

---

## 3. Five Core Rules

### Rule 1: Atom Map Conservation
Reactant and product sides must have the **same numbers and types of mapped atoms**. Unmapped atoms may appear (added) or disappear (deleted).

```
✓  [C:1][Cl:2]>>[C:1][I:2]    (same mapped atoms, bond changes)
✓  [C:1][Cl]>>[C:1]           (Cl unmapped → can be dropped)
```

### Rule 2: Stoichiometry is 1:1
Each mapped atom appears exactly once on each side.

### Rule 3: Explicit Hydrogens Must Be Mapped
If explicit H appears on one side, it **must** appear on the other side with the same map number.

```
✓  [C:1][H:2]>>[C:1]          — WRONG: H must appear on product side too
✓  [C:1][H:2]>>[O:3][H:2]     — H migrates from C to O (correctly mapped)
```

### Rule 4: Bond Expressions = SMILES Only (No Queries)
Bond expressions in SMIRKS must be valid SMILES bonds, not SMARTS bond queries.

```
✓  [C:1]=[C:2]>>              (SMILES double bond)
✗  [C:1]~[C:2]>>              (~ is a SMARTS wildcard, not allowed)
```

### Rule 5: Atom Expression Restrictions
- **Atoms with bond changes**: must use SMILES atom expressions
- **Atoms without bond changes**: may use SMARTS expressions (for context restriction)

```
[*:1][N:2](=[O:3])=[O:4]>>[*:1][N+:2](=[O:3])[O-:4]
```
Here `:1` and `:3` don't change bonding → can use SMARTS `[*]`.
`:2` and `:4` change bonding → must use SMILES.

---

## 4. Hydrogen Semantics (post-v4.51)

| Expression | Meaning |
|------------|---------|
| `[H]` | Hydrogen **atom** (explicit, like `[#1]`) |
| `[H1]` | Atom with one attached hydrogen |
| `[#1]` | Hydrogen atom |

This matters for tracking H migration in reactions.

---

## 5. Examples

### 5.1 Nitro Group Interconversion (non-chemical transform)

```smirks
[*:1][N:2](=[O:3])=[O:4]>>[*:1][N+:2](=[O:3])[O-:4]
```

Converts pentavalent N representation to charge-separated trivalent form.

> Note: transforms need not represent real reactions — they are general molecular manipulation tools.

### 5.2 Amide Bond Formation

```smirks
[C:1](=[O:2])[Cl:3].[H:99][N:4]([H:100])[C:5]>>[C:1](=[O:2])[N:4]([H:100])[C:5].[Cl:3][H:99]
```

- H at `:99` changes bonding (N–H breaks, H–Cl forms) → SMILES `[H]`
- H at `:100` retains bonding → can use SMARTS

### 5.3 Tetrahedral Stereochemistry Inversion

```smirks
[*:1][C@:2]([*:3])([*:4])[*:5]>>[*:1][C@@:2]([*:3])([*:4])[*:5]
```

Inverts a tetrahedral center by flipping `@` ↔ `@@`.

### 5.4 Double Bond Stereochemistry Inversion

```smirks
[*:1]/[C:2]([*:3])=[C:4]([*:5])/[*:6]>>[*:1]/[C:2]([*:3])=[C:4]([*:5])\[*:6]
```

Changes E-alkene to Z-alkene (or vice versa). Stereochemistry transforms need "sufficient context."

### 5.5 Halide Exchange (Finkelstein)

```smirks
[C:1][Cl:2].[I-:3]>>[C:1][I:3].[Cl-:2]
```

### 5.6 Generic Esterification

```smirks
[C:1](=[O:2])[OH:3].[O:4][C:5]>>[C:1](=[O:2])[O:4][C:5].[OH2:3]
```

---

## 6. Stereochemistry in SMIRKS

- Stereodescriptors (`@`, `@@`, `/`, `\`) are **local** — referenced to map label ordering
- Transforms **involving** stereocenters should include sufficient context to unambiguously specify the geometry
- Unspecified stereocenters in reactants are not constrained (match any)
- Unspecified in products means stereocenters are not set (retain or destroy)

---

## 7. SMIRKS vs. SMILES vs. SMARTS

| Feature | SMILES | SMARTS | SMIRKS |
|---------|--------|--------|--------|
| Encodes | Specific molecule | Search pattern | Reaction transform |
| Atom expressions | SMILES only | SMARTS primitives | SMILES (bonded) or SMARTS (unbonded) |
| Bond expressions | SMILES | SMARTS (queries allowed) | SMILES only |
| Atom maps | Optional | Optional | **Required** for changed atoms |
| Reaction separator | `>` | `>` | `>` |

---

## 8. Using SMIRKS in RDKit

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Define transform
rxn = AllChem.ReactionFromSmarts('[C:1][Cl:2]>>[C:1][I:2]')

# Apply to molecule
mol = Chem.MolFromSmiles('CCCl')
products = rxn.RunReactants((mol,))

for product_set in products:
    for p in product_set:
        print(Chem.MolToSmiles(p))  # CCI
```

Key distinction: `ReactionFromSmarts()` handles SMIRKS (it interprets the `>>` syntax with atom maps).
