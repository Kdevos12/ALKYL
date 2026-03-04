# Cheminformatics Database Design Concepts

> **Note:** Chapters 7 (THOR) and 8 (Merlin) of the Daylight Theory Manual describe proprietary software. This reference extracts the **transferable conceptual content** applicable to any chemical database or cheminformatics system.

---

## 1. Hash Tables for Chemical Data

**Problem:** Sequential search through millions of records is O(n) — too slow for interactive use.

**Solution:** Hash tables provide **O(1)** average-case lookup by computing a hash function on the identifier → direct pointer to record location.

Key properties:
- Access time **independent of database size**
- Typically requires 2 disk accesses: one to hash table, one to data file
- Hash table size should ≈ number of records (collision chain length → 1.2–1.5)

**Application in cheminformatics:** Using canonical SMILES or InChI as the hash key enables direct lookup of any molecule by structure, without searching.

---

## 2. Identifier Design Principles

Two classes of database identifiers:

| Class | Nature | Examples |
|-------|--------|---------|
| Arbitrary | Assigned without structural meaning | Registry numbers, catalog IDs, trivial names |
| Structural | Derived from molecular structure itself | IUPAC names, canonical SMILES, InChI |

**Properties of good structural identifiers:**
- Unambiguous and computationally parseable
- Non-arbitrary (derived from structure → same structure, same identifier)
- Universal across databases
- Enable direct lookup without searching
- Canonical (one molecule → one string)

Canonical SMILES and InChI satisfy all these criteria.

---

## 3. Thesaurus / Data Tree Structure

Chemical databases benefit from a **hierarchical thesaurus** model:

```
Root identifier (canonical SMILES / InChI)
├── Names
│   ├── IUPAC name
│   ├── Common name(s)
│   └── Trade name(s)
├── Properties
│   ├── MW, LogP, TPSA, ...
│   └── Experimental data
├── Reactions
│   ├── Reaction as reactant
│   └── Reaction as product
└── References
    └── Literature DOIs
```

This model:
- Groups related information under one key
- Accommodates synonyms without requiring unique names
- Allows "one record, many views"

---

## 4. Datatypes and Standardization

**Datatypes** define the schema for stored data — what each field represents, its format, and how it should be interpreted.

**Standardization** improves retrieval consistency by normalizing data before storage:
- Case normalization
- Punctuation removal
- Whitespace normalization
- Canonical form conversion

Trade-off: minor information loss (e.g., case) for major gain in retrieval consistency. The canonical SMILES stored in a database is a standardization of the input SMILES.

---

## 5. Multi-Database Architecture

Separating concerns into distinct stores:

| Layer | Content |
|-------|---------|
| Schema DB | Datatype definitions, field types |
| Reference DB | Indirect references, cross-links |
| Chemical DB | Structures, properties, reactions |

Benefits:
- Schema can be shared across sites
- Chemical data remains independent
- Modular updates (update properties without re-indexing structures)

---

## 6. In-Memory Search Architecture

### Why In-Memory?

Memory access is ~**100,000× faster** than disk access. This transforms the interaction model:

| Disk-based | In-Memory |
|-----------|-----------|
| Carefully formulate query to minimize wait | Explore interactively |
| Batch queries | Real-time hypothesis testing |
| Data retrieval ("what is X?") | Exploratory analysis ("what looks like X?") |

Modern workstations can load:
- Entire chemical databases (10–15M structures) on large servers
- Corporate libraries (100K–1M) trivially on workstations
- Screening speed: 100,000–1,000,000 structures/second

### Pools (In-Memory Database Views)

A **pool** = entire database loaded into RAM as a "chemical spreadsheet":
- Each row = one compound record
- Enables simultaneous operations across all records
- Multiple pools can coexist (different databases)

In modern systems: equivalent to `pandas.DataFrame` + structure column, fully in memory.

### Columns and Cells

**Columns** provide vertical slices through the pool — one property across all compounds.

When a compound has **multiple values** for a property (e.g., multiple names):

| Selector | Meaning |
|----------|---------|
| First / Last | By order of appearance |
| Longest / Shortest | By string length |
| Least / Greatest | By value |
| Count | Number of values |
| Average | Numerical mean |
| Standard Deviation | Spread |

This separates **data storage** from **presentation**.

### Hitlists (Dynamic Result Sets)

A **hitlist** = ordered subset of pool rows, maintained dynamically:
- Modified by search operations (add/remove rows matching a query)
- Reordered by sort operations
- Multiple concurrent hitlists enable "undo" and saved result preservation

In modern systems: equivalent to filtered and sorted `DataFrame` views or RDKit `SubstructMatchParameters` results stored as index arrays.

---

## 7. Reaction Database Design

Reactions require special handling:

| Aspect | Consideration |
|--------|--------------|
| Atom maps | Track which reactant atoms become which product atoms |
| Agent molecules | Store separately; not part of structural fingerprint |
| Stoichiometry | Stored externally (not in the reaction SMILES/SMARTS) |
| Conditions | Temperature, pressure, yield — stored as annotations |
| Reaction fingerprints | Structural (reactant+product OR) or difference (bond changes) |

Key design principle: **separate the structural graph from the experimental annotation**. The structural reaction (what bonds change) is the key; all other data hangs off it.
