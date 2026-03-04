# Daylight Theory — SMILES Specification

SMILES (Simplified Molecular Input Line Entry System) is a **line notation** — a typographical method using printable ASCII characters — for encoding molecules and reactions.

## SMILES Variants

| Variant | Definition |
|---------|-----------|
| Generic SMILES | Any valid SMILES; isotopes/chirality may be omitted |
| Unique SMILES | Single canonical form; no isotopes or chirality |
| Isomeric SMILES | Includes isotopes and/or chirality |
| Absolute SMILES | Unique isomeric SMILES — fully specified canonical form |

Canonicalization ensures: `OCC`, `C(O)C`, `C-C-O` → all produce `CCO`.

---

## 1. Atoms

### Organic Subset (no brackets needed at normal valence)
B, C, N, O, P, S, F, Cl, Br, I

- Lowercase = aromatic: `c`, `n`, `o`, `p`, `s`, `se`, `as`
- Examples: `C` (methane), `c1ccccc1` (benzene)

### Bracketed Atoms (required when)
- Element not in organic subset: `[Au]`, `[Fe]`, `[Se]`
- Abnormal valence: `[S]` (elemental sulfur, needs brackets)
- Explicit H count or charge needed

**Syntax:** `[isotopeSYMBOLHncharge]`

| Component | Syntax | Example |
|-----------|--------|---------|
| Isotope | integer before symbol | `[13C]`, `[2H]` |
| Symbol | element symbol | `[Fe]` |
| Hydrogens | `H` + optional digit | `[NH3]`, `[CH4]` |
| Charge | `+n` or `-n` | `[Fe+2]`, `[OH-]`, `[NH4+]` |

**Examples:**
- `[H+]` — proton
- `[OH3+]` — hydronium
- `[Fe++]` or `[Fe+2]` — iron(II) (both valid)
- `[13CH4]` — C-13 methane
- `[235U]` — uranium-235

---

## 2. Bonds

| Symbol | Bond |
|--------|------|
| `-` | Single (implicit, usually omitted) |
| `=` | Double |
| `#` | Triple |
| `:` | Aromatic |
| `/` | Directional (stereochemistry) |
| `\` | Directional (stereochemistry) |

Default between adjacent atoms: single or aromatic (auto-detected).

```
CC      ethane
C=O     formaldehyde
C#N     hydrogen cyanide
c:c     aromatic bond (usually implicit in ring)
```

---

## 3. Branches

Branches are enclosed in parentheses. The branch connects back to the atom immediately to the **left**.

```
CCN(CC)CC           triethylamine
CC(C)C(=O)O         isobutyric acid
CC(=O)Oc1ccccc1C(=O)O   aspirin
```

Nesting is allowed: `CC(C)(C)C` (neopentane — 3 methyls on central C).

---

## 4. Ring Closures

Break one bond per ring by inserting the **same digit** at both ends of the broken bond.

```
C1CCCCC1        cyclohexane
c1ccccc1        benzene
C12C3C4C1C5C4C3C25   cubane
```

Rules:
- Digits 0–9; reusable after closure: `O1CCCCC1N1CCCCC1` (two rings both labeled `1`)
- Numbers ≥10: use `%` prefix: `C%13...C%13`
- One atom can bear multiple ring closures: `C12C3C4C1C5C4C3C25` (cubane)

---

## 5. Disconnected Structures

Components separated by `.` (period). No implied charge pairing.

```
[Na+].[Cl-]         sodium chloride
[NH4+].[Cl-]        ammonium chloride
```

Special: matching ring digits across `.` bond the atoms:
```
C1.C1   →   CC  (ethane — unusual but valid)
```

---

## 6. Isomeric SMILES

### 6.1 Isotopic Specification

```
[12C]       carbon-12
[13C]       carbon-13
[C]         unspecified mass carbon
[2H]O[2H]  deuterium oxide (heavy water)
```

### 6.2 Double Bond Stereochemistry

Use `/` and `\` to indicate relative orientation across a double bond:

```
F/C=C/F     (E)-1,2-difluoroethene (same side)
F/C=C\F     (Z)-1,2-difluoroethene (opposite sides)
F\C=C/F     (Z) — alternate notation
```

Rules:
- Symbols indicate **relative** directionality (not absolute)
- Both atoms flanking the double bond must be annotated
- Partial specification is legal: `F/C=C/C=CC` (second DB unspecified)

### 6.3 Tetrahedral Chirality

**Symbols:** `@` (anticlockwise) / `@@` (clockwise) — viewed from first neighbor toward the center.

```
N[C@@H](C)C(=O)O    D-alanine
N[C@H](C)C(=O)O     L-alanine
```

**Reading order of neighbors:**
1. The atom from which we arrive (implicit — previous atom or `[H]` if first)
2. Remaining neighbors in written order

Equivalent representations for D-alanine:
```
N[C@@H](C)C(=O)O
N[C@@](C)(C(=O)O)[H]   (explicit H)
[H][C@](N)(C)C(=O)O
```

### 6.4 Extended Chirality Classes

Full syntax: `@[CLASS][NUMBER]`

| Class | Symbol | Geometry | Values |
|-------|--------|----------|--------|
| Tetrahedral | `@TH` | tetrahedral | 1–2 |
| Allene-like | `@AL` | allenic/spiro | 1–2 |
| Square-planar | `@SP` | SP4 | 1–3 (U, 4, Z shapes) |
| Trigonal-bipyramidal | `@TB` | TBP | 1–20 |
| Octahedral | `@OH` | octahedral | 1–30 |

Shortcuts: `@` = `@TH1`, `@@` = `@TH2`

**Allene example:**
```
OC(Cl)=[C@]=C(C)F     (fully specified allene)
```

**Square-planar (SP1 = U-shape):**
```
F[Po@SP1](Cl)(Br)I
```

---

## 7. Conventions

### 7.1 Hydrogen Treatment

Three modes:
1. **Implicit**: assumed from valence for organic subset atoms
2. **Explicit by count**: `[NH3]`, `[CH4]`
3. **Explicit atoms**: `[H]` as a node

Require explicit H when:
- Charged: `[H+]`
- H–H bond: `[H][H]`
- Isotopic: `[2H]`
- Involved in reaction atom mapping

### 7.2 Aromaticity

Detection: all ring atoms sp², excess electrons = 4N+2.

```
c1ccccc1        benzene (aromatic notation)
C1=CC=CC=C1     benzene (Kekulé — auto-converted to aromatic)
c1cocc1         furan
```

Invalid: `c1cccc1` — only 5 atoms, cannot assign alternating bonds.

### 7.3 Aromatic Nitrogen

All aromatic nitrogens written as lowercase `n`:

```
n1ccccc1        pyridine (lone pair not in ring)
O=n1ccccc1      pyridine-N-oxide
[nH]1cccc1      pyrrole (lone pair in ring, H on N)
```

Pyrrole can also be written `N1C=CC=C1` (Kekulé form).

### 7.4 Bonding Conventions

SMILES does not enforce one valence model — both are valid:

```
CN(=O)=O        nitromethane (pentavalent N)
C[N+](=O)[O-]   nitromethane (charge-separated, preferred)
C=[N+]=[N-]     diazomethane (charge-separated, preferred)
```

### 7.5 Tautomers

Explicit tautomeric form required — no "mobile hydrogen" concept:

```
O=c1[nH]cccc1   2-pyridone form
Oc1ncccc1       2-pyridinol form
```

---

## 8. Reaction SMILES

**Syntax:** `reactant > agent > product` or `reactant >> product`

```
C=CCBr>>C=CCI                              (no agent)
C=CCBr.[Na+].[I-]>CC(=O)C>C=CCI.[Na+].[Br-]  (with acetone as agent)
```

Empty components allowed: `>>` (null reaction).

### 8.1 Atom Maps

Syntax: `[SYMBOL:N]` where N is a non-negative integer map class.

```
[CH3:1][Br:2]>>[CH3:1][I:2]    atom-mapped halide exchange
```

Rules:
- Appear in any SMILES (ignored in non-reaction processing)
- Correlate reactant atoms with product atoms
- Not required to be complete or unique
- Can express mechanistic ambiguity
- Agents lose atom maps during canonicalization
- Mapped hydrogens appear explicitly in absolute SMILES

In unique SMILES: atom maps are **dropped**.
In absolute SMILES: atom maps are **preserved**.
