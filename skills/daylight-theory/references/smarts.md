# Daylight Theory — SMARTS Query Language

SMARTS (SMiles ARbitrary Target Specification) is a **pattern language** for specifying molecular substructures. SMARTS is a superset of SMILES: every valid SMILES is a valid SMARTS, but SMARTS adds query-specific primitives and logic.

**Critical difference from SMILES:**
> In SMARTS, **unspecified properties do not restrict matches**. `O` in SMARTS matches any aliphatic oxygen (not just "water"). In SMILES, `O` is a specific molecule.

---

## 1. Atomic Primitives

### Organic Subset (SMILES-compatible)
`B C N O P S F Cl Br I` — match that specific element at normal valence.
Lowercase (`c n o p s`) — match aromatic atom of that element.

### Query-Specific Primitives

| Primitive | Meaning |
|-----------|---------|
| `*` | Any atom (wildcard) |
| `a` | Any aromatic atom |
| `A` | Any aliphatic atom |
| `#n` | Atomic number n (e.g., `[#6]` = any carbon) |
| `D<n>` | Exactly n explicit connections (degree) |
| `H<n>` | Exactly n total H attached (incl. implicit) |
| `h<n>` | Exactly n implicit H |
| `R<n>` | Atom is in exactly n SSSR rings (or `R` = any ring atom) |
| `r<n>` | Atom is in a SSSR ring of size n |
| `v<n>` | Total bond order (valence) = n |
| `X<n>` | Total connectivity (explicit + implicit H) = n |
| `x<n>` | Ring connectivity = n |
| `+<n>` | Formal charge = +n |
| `-<n>` | Formal charge = -n |
| `@` / `@@` | Chirality: anticlockwise / clockwise |

**Examples:**
```
[#6]        any carbon (aromatic or aliphatic)
[++]        atom with +2 charge
[D3]        atom with exactly 3 explicit bonds
[R]         any ring atom
[r5]        atom in a 5-membered ring
[X4]        atom with connectivity 4 (e.g., quaternary C)
[v4]        atom with total bond order 4
[H1]        atom with exactly one H
[h0]        atom with no implicit H
```

**Practical patterns:**
```
[OH]        sp3 oxygen with one H (hydroxyl)
[OH1]       same (explicit count)
[O;D2]      ether oxygen (2 connections, no H)
[NX3;H2]    primary amine nitrogen
[nH1]       pyrrole-type N (aromatic, 1H)
[CX4]       sp3 carbon (4 connections)
[CX3]       sp2 carbon (3 connections)
[c]         aromatic carbon
[n]         aromatic nitrogen
```

---

## 2. Bond Primitives

| Symbol | Bond type |
|--------|----------|
| `-` | Single (aliphatic) |
| `=` | Double |
| `#` | Triple |
| `:` | Aromatic |
| `~` | Any bond (wildcard) |
| `/` | Directional up |
| `\` | Directional down |
| `@` | Any ring bond |

**Examples:**
```
c:c         aromatic bond between aromatic carbons (benzene ring)
c-c         single bond between two aromatic carbons (biphenyl junction)
*~*         any two atoms connected by any bond
*!@*        any bond that is NOT a ring bond
*@;!:*      non-aromatic ring bond
C=C         aliphatic double bond
[C,c]=,#[C,c]   double or triple bond between carbons (any aromaticity)
```

---

## 3. Logical Operators

Operators combine primitives within `[...]`:

| Operator | Syntax | Precedence | Meaning |
|----------|--------|-----------|---------|
| NOT | `!e` | Highest | Negation |
| AND (high) | `e1&e2` or `e1e2` | High | Both required (tight) |
| AND (low) | `e1;e2` | Low | Both required (loose) |
| OR | `e1,e2` | Lowest | Either |

**Precedence matters for grouping:**
```
[N,O;+1]    → [N,O] AND [+1]  → (N or O) with +1 charge
[N,O+1]     → N or (O+1)     → N, or O with +1 charge
[!C;R]      → (not C) and ring atom
[n,o;H1]    → (n or o) with exactly 1H
```

**Common patterns:**
```
[OH]            hydroxyl oxygen (O + H1)
[O;H1]          same, explicit AND
[c,n;H1]        aromatic C or N with one H
[!#1]           not hydrogen
[!C;!c;!H]      not any carbon, not hydrogen
[F,Cl,Br,I]     any halogen (F, Cl, Br, I)
[CH2]           aliphatic C with 2H (methylene)
[CX4H3]         methyl carbon (4 connections, 3H)
```

---

## 4. Recursive SMARTS

Syntax: `$(SMARTS)` embeds a full SMARTS as an atomic property.

```
$(*C)           atom bonded to a carbon
$(*CC)          atom bonded to carbon-carbon
$([OH]c)        oxygen bonded to aromatic carbon (phenol OH)
```

Recursive SMARTS allow complex environment queries:
```
C[$(aaO)]          carbon adjacent (aromatic path) to oxygen
C[$(aaO);$(aaaN)]  carbon ortho to O AND meta to N
[$(*C);$(*CC)]     atom bonded to both methyl and ethyl
```

Recursive expressions are **atomic properties** and compose with other primitives:
```
[O;$(aaO)]          aromatic oxygen also adjacent via aromatic path to O
[NX3;$([NH2]c)]     sp3 N with 2H bonded to aromatic C (aniline)
```

---

## 5. Component-Level Grouping

Parentheses `()` around dot-disconnected fragments enforce same-component matching.

| Pattern | Meaning |
|---------|---------|
| `C.C` | Two carbons in any two components |
| `(C.C)` | Two carbons in the **same** component |
| `(C).(C)` | Two carbons in **different** components |

```
(C.C)           matches ethane (C-C in same molecule)
(C).(C)         requires two separate molecules, one C each
(C.C).(N)       C-C in one component, N in a different component
```

Critical for specifying **inter-** vs. **intra-molecular** reactions.

---

## 6. Reaction Queries

Format: `reactants>>products` or `reactants>agents>products`

When matching reactions:
- **Molecule SMARTS** matches anywhere in the reaction
- **Atom maps** in SMARTS restrict associations (but never expand them)

Map syntax in SMARTS:
```
[expr:n]    mapped atom, map class n
[expr:?n]   mapped or unmapped atom, class n if mapped
```

Examples:
```
C>>             methyl carbon appearing in reactants
>>C             methyl carbon appearing in products
[C:1]>>[C:1]   C in reactant mapped to C in product (same class)
```

---

## 7. SMARTS vs. SMILES: Key Differences

| Aspect | SMILES | SMARTS |
|--------|--------|--------|
| Purpose | Encode a specific molecule | Specify a search pattern |
| `O` matches | Water (specific molecule) | Any aliphatic oxygen |
| `C1=CC=CC=C1` + benzene | ✓ (Kekulé auto-aromatized) | ✗ (literal non-aromatic pattern) |
| `c1ccccc1` + benzene | ✓ | ✓ |
| `[CH4]` | Methane specifically | Any carbon with 4 H |
| Unspecified property | Defaults to normal valence | **No restriction imposed** |

**Important:** A SMILES string used as a SMARTS query may or may not behave as expected. Kekulé benzene `C1=CC=CC=C1` as SMARTS does **not** match the aromatic form.

---

## 8. Efficiency Guidelines

When atom-order optimization is unavailable, list **least common** features **earliest** in AND-expressions:

```
# Efficient: rare feature first
[N;R;H1]   → try R first (ring N is rarer than any N)

# Efficient OR: common feature last
[O,N,S]    → S least common → list last? No: in OR, list more common first
             Actually: for OR, list more common FIRST (short-circuit)
```

Rules:
- In **AND** expressions: put rare/expensive constraints **first** (fail fast)
- In **OR** expressions: put common alternatives **first** (succeed fast)
- Order of atoms in query SMARTS: put distinctive atoms **early**

---

## 9. Reference Patterns

### Functional Group Matching
```
[OH]c                   phenol
[OH]C                   aliphatic alcohol
[OH]                    any alcohol or phenol
C(=O)[OH]               carboxylic acid
C(=O)[O-]               carboxylate anion
C(=O)N                  amide
C(=O)Cl                 acid chloride
[NH2]c                  primary aromatic amine
[NX3;H2;!$(NC=O)]       primary aliphatic amine (not amide)
[NX3;H1;!$(NC=O)]       secondary amine
[SX2H]                  thiol
[#16X2H]                thiol (sulfur)
[CX3](=O)[OX2H1]        carboxylic acid (precise)
```

### Substructure Matching
```
c1ccccc1                benzene ring (any substituents)
c1ccc(*)cc1             para-substituted benzene
c1cc(*)ccc1*            ortho-disubstituted benzene
[R]                     any ring atom
*!@*                    any rotatable bond
[#6]~[#6]~[#6]          three connected carbons (any bond type)
[F,Cl,Br,I]             halogen
[#7,#8]                 nitrogen or oxygen
```

### Drug-Likeness / Toxicophore Alerts
```
[NH2]c                  aniline (CYP1A2 substrate)
[NX3][CH3]              N-methyl (N-dealkylation)
[SX2][CH3]              S-methyl (sulfoxidation)
c1cccc2ccccc12          naphthalene (CYP1A2 substrate)
[OX2H1;!$(OC=O)]        alcohol (glucuronidation site)
```
