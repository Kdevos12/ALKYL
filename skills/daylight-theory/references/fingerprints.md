# Daylight Theory — Fingerprints & Similarity

## The Screening Problem

Substructure searching is **NP-complete** — worst-case O(K^N). Real molecules rarely hit worst cases (typically O(N²) to O(N³)), but databases of millions of molecules make even polynomial searches slow.

**Core insight:** Detecting that a pattern is **absent** from a molecule can be done in O(N) time. Screening uses this to quickly eliminate non-matches before expensive graph isomorphism.

---

## 1. Molecular Formula Screening

Earliest approach: if the query contains an atom that the candidate lacks, the candidate is eliminated with 100% confidence. Fast but very crude.

---

## 2. Structural Keys

A **structural key** is a fixed-length bit array where each bit represents the presence/absence of a predefined structural feature.

| Property | Value |
|----------|-------|
| Bit meaning | Predefined (each bit = specific substructure) |
| Size | Fixed (tens to thousands of bits) |
| Density | Sparse (most bits 0 for typical molecules) |
| Generation | Requires substructure search for each feature |
| Screening | Molecule passes if all pattern bits are set in molecule |

Examples of features encoded: element presence, ring systems, functional groups, electronic configurations, disjunctions of rare fragments.

**Limitation:** A new feature requires redefining the key for all databases.

---

## 3. Path-Based Fingerprints

Fingerprints solve the generality problem: **no bit has a predefined meaning**.

### Generation Algorithm

1. Enumerate all paths from 0 bonds to 7 bonds through the molecular graph
2. Each path is a sequence of atoms and bonds
3. Each path is **hashed** → pseudo-random set of 4–5 bits
4. All bit sets are ORed together into the final fingerprint

**Path examples:**
- 0-bond: `C`, `O`, `N` (individual atoms)
- 1-bond: `OC`, `C=C`, `CN`
- 2-bond: `OC=C`, `C=CN`
- ...up to 7-bond paths

### Key Property

> If pattern P is a substructure of molecule M, then every bit set in P's fingerprint is also set in M's fingerprint.

This is the **necessary condition** for substructure screening. A 0-bit in the molecule guarantees the pattern is absent; a 1-bit is necessary but not sufficient.

### Advantages Over Structural Keys

| Aspect | Structural Keys | Path Fingerprints |
|--------|----------------|-------------------|
| Universality | Database-specific | Universal (same algorithm everywhere) |
| Generality | Predefined features only | All paths, automatically |
| Bit density | Sparse | 20–40% (better discrimination) |
| New features | Requires schema change | Automatic |
| False positive rate | Higher | Lower |

---

## 4. Folding (Variable-Sized Fingerprints)

Different molecules need different fingerprint sizes for optimal density. **Folding** creates shorter fingerprints while preserving screening validity:

1. Start with a large fingerprint (e.g., 2048 bits)
2. Split into two equal halves
3. OR the halves together
4. Repeat until target density reached

**Critical property:** Folding is **screening-safe** in one direction:
- If a screen was **negative** before folding → it remains negative after (no false negatives introduced)
- A positive screen may become negative after folding (increased false positive rate, acceptable for screening)

This trade-off: smaller storage, higher false positive rate, but zero false negative rate preserved.

---

## 5. In-Memory Screening

Memory is ~10⁵× faster than disk. Modern workstations can load:
- 10–15 million known structures (entire CAS-like databases)
- Corporate libraries of hundreds of thousands: trivially

**Speed:** 100,000 to 1,000,000 structures screened per second.

This makes interactive chemical exploration possible — search time < query formulation time.

---

## 6. Reaction Fingerprints

### 6.1 Structural Reaction Fingerprints

Combines reactant and product fingerprints via OR:
1. Fingerprint of reactant part
2. Fingerprint of product part
3. Bit-shifted fingerprint of product part (to distinguish R vs. P)

Agents are excluded. Standard screening/similarity operations apply.

### 6.2 Reaction Difference Fingerprints

Specialized for identifying **bond changes** in stoichiometric reactions.

Simple XOR fails because fingerprints don't encode path counts — bonds present on both sides would be masked.

**Solution:** Track per-reaction **path count differences**:
- Enumerate paths in reactants (with counts)
- Enumerate paths in products (with counts)
- Set fingerprint bits only for paths where count differs (non-zero Δ)

**Example — Sn2:** `[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI`

| Path length | Changed paths |
|------------|---------------|
| 0-bond | None (same atoms present) |
| 1-bond | C–Br→C–I |
| 2-bond | C–C–Br→C–C–I |
| 3-bond | C=C–C–Br→C=C–C–I |

Only non-zero Δ paths set bits.

**Limitation:** Difference fingerprints **cannot be used as substructure screens** — they violate the subset property required for screening.

---

## 7. Similarity Measures

Similarity is computed from a **2×2 contingency table** of bit comparisons between two fingerprints A and B:

| | B=0 | B=1 |
|-|-----|-----|
| **A=0** | d (both off) | b (B only) |
| **A=1** | c (A only) | a (both on) |

**Notation:**
- `a` = bits on in A only
- `b` = bits on in B only
- `c` = bits on in **both** A and B (common features)
- `d` = bits off in both
- `A` = total bits on in A = a + c
- `B` = total bits on in B = b + c
- `n` = total bits = a + b + c + d

### 7.1 Tanimoto Coefficient (Jaccard)

```
Tc = c / (a + b + c)
```

- Range: [0, 1]
- Interprets as: proportion of "on-bits" that are shared
- **Double-zero independent** (does not use d)
- **Standard choice for path fingerprints**
- Symmetric: Tc(A,B) = Tc(B,A)

```python
from rdkit import DataStructs
sim = DataStructs.TanimotoSimilarity(fp1, fp2)
```

### 7.2 Dice Coefficient

```
Dice = 2c / (2c + a + b)
```

- Range: [0, 1]
- Arithmetic mean version; weights common features more than Tanimoto
- Also double-zero independent

### 7.3 Cosine Coefficient

```
Cosine = c / sqrt(A × B)
```

- Range: [0, 1]
- Geometric mean; from vector algebra

### 7.4 Tversky Index (Asymmetric Similarity)

```
Tversky(α, β) = c / (α·a + β·b + c)
```

Asymmetric — models "prototype vs. variant" comparisons:

| α | β | Result |
|---|---|--------|
| 1 | 0 | "Superstructure-likeness": 1.0 = all prototype features present in variant |
| 0 | 1 | "Substructure-likeness": 1.0 = variant completely embedded in prototype |
| 1 | 1 | Tanimoto index |
| 0.5 | 0.5 | Dice index |

Use case: find all molecules containing a query scaffold (α=0, β=1 → substructure-likeness).

### 7.5 Additional Similarity Measures

All functions of (a, b, c, d):

| Measure | Formula | Range | Notes |
|---------|---------|-------|-------|
| Simpson | c/min(A,B) | [0,1] | Best individual substructure similarity |
| Kulczynski | c/(½(a+b)) | [0,1] | Mean of individual substructure similarities |
| Forbes | cn/[AB+(A-c)(B-c)] | [0,∞] | No upper limit |
| Russell-Rao | c/n | [0,1] | Uses total bits n |
| Rogers-Tanimoto | (c+d)/(a+b+2c+2d) | [0,1] | Double-zero dependent |
| Hamman | (c+d-a-b)/n | [-1,1] | Double-zero dependent |
| Pearson | (cd-ab)/√(ABXY) | [-1,1] | Double-zero dependent |
| Yule | (cd-ab)/(cd+ab) | [-1,1] | Double-zero dependent |
| Manhattan (distance) | (a+b)/n | [0,1] | Distance metric |
| Euclidean | (a+b)/n | [0,1] | Hamming distance form |

### 7.6 Two Coefficient Classes

**Class 1 (double-zero dependent):** Includes `d` (shared off-bits). Problem: for path fingerprints, size can be arbitrarily doubled by adding irrelevant off-bits → `d` is not meaningful.

**Class 2 (double-zero independent):** Ignores `d`. Correct choice for Daylight-style path fingerprints where fingerprint length is variable.

> **Rule:** For path-based fingerprints, use **double-zero independent** measures (Tanimoto, Dice, Cosine, Tversky).

### 7.7 User-Defined Measures

Custom similarity functions f(a, b, c, d) with standard math operations. Requirements:
- Consistent feature inclusion (denominator must not exclude what numerator includes)
- Understand range restrictions
- Daylight fingerprints: avoid d-dependent measures

---

## 8. Practical Guidelines

### Choosing a Similarity Measure

| Use Case | Recommended Measure |
|----------|-------------------|
| General molecular similarity | Tanimoto |
| Finding near-substructures | Tversky (α=0, β=1) |
| Finding molecules containing scaffold | Tversky (α=0, β=1) |
| Comparing feature-rich to feature-poor | Tversky (asymmetric) |
| Clustering | Tanimoto or Dice |

### Thresholds (Rule of Thumb — Morgan/ECFP4)

| Tanimoto | Structural relationship |
|----------|------------------------|
| > 0.85 | Very similar (same scaffold) |
| 0.6–0.85 | Similar |
| 0.4–0.6 | Related |
| < 0.4 | Dissimilar |

Note: thresholds are fingerprint-type dependent.

### Research Note

Studies (Holliday et al., 2002) show that diverse similarity measures tend to produce similar **rankings** — measure choice matters less than consistency. However, different measures are **not monotonic** with each other for absolute values (Hubalek, 1982).

---

## 9. RDKit Implementation

```python
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys

mol1 = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')  # aspirin
mol2 = Chem.MolFromSmiles('Cn1cnc2c1c(=O)n(C)c(=O)n2C')  # caffeine

# Morgan (ECFP4-like)
fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)

# Tanimoto
tc = DataStructs.TanimotoSimilarity(fp1, fp2)

# Dice
dice = DataStructs.DiceSimilarity(fp1, fp2)

# Tversky (substructure-likeness: α=0, β=1)
tv = DataStructs.TverskySimilarity(fp1, fp2, 0.0, 1.0)

# Bulk similarity (against a library)
sims = DataStructs.BulkTanimotoSimilarity(fp1, [fp2, ...])
```
