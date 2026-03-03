# ALKYL PhD-Level Scripts — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add 4 scripts that close the gap between ALKYL (M2 level) and a PhD computational chemist.

**Architecture:** Same pattern as existing scripts — standalone Python, lazy RDKit imports, JSON stdout, `--out` flag, pytest tests in `tests/`. Run tests with `/home/de/Bureau/Pepflex/venv/bin/python -m pytest`.

**Tech Stack:** RDKit (already available), rdkit.Contrib.SA_Score, pure Python SMARTS rules. Dimorphite-DL NOT installed — use RDKit MolStandardize ionization instead.

---

## Priority & Rationale

| # | Script | Gap filled | Difficulty |
|---|---|---|---|
| 1 | `chem_pka.py` | Protonation states at pH — foundational for ADMET | Medium |
| 2 | `chem_qm.py` (extend) | Parse IR frequencies + NMR shifts from ORCA output | Medium |
| 3 | `chem_metabolism.py` | CYP450 soft spots via SMARTS | Easy |
| 4 | `chem_diversity.py` | MaxMin diversity selection from a library | Easy |

---

## Task 1: `chem_pka.py` — Protonation states at given pH

**Files:**
- Create: `scripts/chem_pka.py`
- Test: `tests/test_chem_pka.py`

**Interface:**
```
python chem_pka.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --ph 7.4
```

**Output JSON:**
```json
{
  "input_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "ph": 7.4,
  "dominant_form": "CC(=O)Oc1ccccc1C(=O)[O-]",
  "ionization_groups": [
    {"smarts": "C(=O)[OH]", "type": "carboxylic_acid",
     "pka_range": [3.5, 5.0], "state_at_ph": "ionized", "count": 1}
  ],
  "net_charge_at_ph": -1
}
```

**Step 1: Write failing tests**
```python
# tests/test_chem_pka.py
from tests.conftest import run_script

SCRIPT = "chem_pka.py"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"   # pKa ~3.5 → ionized at pH 7.4
ANILINE = "Nc1ccccc1"                 # pKa ~4.6 → neutral at pH 7.4

def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    for f in ["input_smiles", "ph", "dominant_form",
              "ionization_groups", "net_charge_at_ph"]:
        assert f in result

def test_aspirin_ionized_at_physiological_ph():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    assert result["net_charge_at_ph"] < 0   # carboxylate at pH 7.4
    assert "[O-]" in result["dominant_form"]

def test_aspirin_neutral_at_low_ph():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "2.0"])
    assert result["net_charge_at_ph"] == 0

def test_aniline_neutral_at_physiological():
    # aniline pKaH ~4.6 → amine deprotonated at pH 7.4
    result = run_script(SCRIPT, ["--smiles", ANILINE, "--ph", "7.4"])
    assert result["net_charge_at_ph"] == 0

def test_ionization_groups_list():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    assert isinstance(result["ionization_groups"], list)
    types = [g["type"] for g in result["ionization_groups"]]
    assert "carboxylic_acid" in types
```

**Step 2: Run to verify FAIL**
```bash
/home/de/Bureau/Pepflex/venv/bin/python -m pytest tests/test_chem_pka.py -v
# Expected: ERROR (script not found)
```

**Step 3: Implement `scripts/chem_pka.py`**

Key design: use a SMARTS-based rule table with typical pKa ranges.
At a given pH, apply Henderson-Hasselbalch: if pH > pKa_max → ionized, if pH < pKa_min → neutral.

```python
#!/usr/bin/env python3
"""ALKYL — chem_pka.py: estimate protonation states at a given pH."""

import argparse, json, sys

# (type, SMARTS, pKa_min, pKa_max, charge_when_ionized, is_acid)
PKA_RULES = [
    ("carboxylic_acid",  "[CX3](=O)[OX2H1]",      3.5,  5.0,  -1, True),
    ("phenol",           "[OX2H1]c",                8.0, 11.0,  -1, True),
    ("thiol",            "[SX2H1]",                 8.0, 11.0,  -1, True),
    ("sulfonamide_NH",   "[#16X4](=[OX1])(=[OX1])N[H]", 9.0, 11.0, -1, True),
    ("phosphate",        "[P](=O)([OH])([OH])",     1.0,  3.0,  -1, True),
    ("amine_primary",    "[NX3;H2;!$(NC=O)]",       9.0, 11.0,  +1, False),  # pKaH
    ("amine_secondary",  "[NX3;H1;!$(NC=O)]",       9.0, 11.0,  +1, False),
    ("amine_tertiary",   "[NX3;H0;!$(NC=O);!$([n])]", 8.0, 10.5, +1, False),
    ("aniline",          "[NX3;H2]c",               3.5,  5.5,  +1, False),
    ("imidazole",        "c1cnc[nH]1",              6.0,  7.5,  +1, False),
    ("pyridine",         "n1ccccc1",                4.5,  6.5,  +1, False),
    ("guanidine",        "NC(=N)N",                12.0, 14.0,  +1, False),
]

def load_mol(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"Invalid SMILES: {args.smiles}", file=sys.stderr); sys.exit(1)
        return mol
    suppl = Chem.SDMolSupplier(args.sdf)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        print(f"No valid molecule in: {args.sdf}", file=sys.stderr); sys.exit(1)
    return mol

def main():
    parser = argparse.ArgumentParser(
        description="Estimate protonation states at a given pH.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles"); inp.add_argument("--sdf")
    parser.add_argument("--ph", type=float, default=7.4,
                        help="Target pH (default: 7.4)")
    parser.add_argument("--out")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = load_mol(args)
    input_smi = Chem.MolToSmiles(mol)
    ph = args.ph

    groups = []
    net_charge = 0

    for (gtype, smarts, pka_min, pka_max, charge_delta, is_acid) in PKA_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            continue
        # Determine state: acid ionizes when pH > pKa_max (simplified)
        if is_acid:
            state = "ionized" if ph > pka_max else (
                    "mixed"   if ph > pka_min else "neutral")
        else:
            # base: protonated (ionized) when pH < pKa_min
            state = "ionized" if ph < pka_min else (
                    "mixed"   if ph < pka_max else "neutral")
        charge = charge_delta if state == "ionized" else 0
        net_charge += charge * len(matches)
        groups.append({
            "type": gtype,
            "pka_range": [pka_min, pka_max],
            "count": len(matches),
            "state_at_ph": state,
            "charge_contribution": charge * len(matches),
        })

    # Build dominant form using RDKit Uncharger/Reionizer heuristic
    # Apply net charge as modification to canonical SMILES
    # (simplified: use Uncharger for neutral pH, keep input otherwise)
    if net_charge == 0:
        uncharger = rdMolStandardize.Uncharger()
        dominant_mol = uncharger.uncharge(mol)
    else:
        # Re-ionize: add charges via reionizer
        reionizer = rdMolStandardize.Reionizer()
        dominant_mol = reionizer.reionize(mol)
    dominant_smi = Chem.MolToSmiles(dominant_mol)

    result = {
        "input_smiles": input_smi,
        "ph": ph,
        "dominant_form": dominant_smi,
        "ionization_groups": groups,
        "net_charge_at_ph": net_charge,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f: f.write(output)
    else:
        print(output)

if __name__ == "__main__":
    main()
```

**Step 4: Run tests**
```bash
/home/de/Bureau/Pepflex/venv/bin/python -m pytest tests/test_chem_pka.py -v
```

**Step 5: Commit**
```bash
git add scripts/chem_pka.py tests/test_chem_pka.py
git commit -m "feat: add chem_pka.py — protonation state estimation at given pH"
```

---

## Task 2: Extend `chem_qm.py` — Parse IR and NMR from ORCA output

**Files:**
- Modify: `scripts/chem_qm.py` (add `--parse-ir` and `--parse-nmr` flags)
- Test: `tests/test_chem_qm.py` (add new test functions)

**New flags:**
```
python chem_qm.py --parse output.log --parse-ir
python chem_qm.py --parse output.log --parse-nmr
```

**IR output JSON:**
```json
{
  "source": "orca",
  "frequencies_cm1": [
    {"mode": 7, "freq_cm1": 1750.3, "intensity": 234.5, "assignment": "C=O stretch"}
  ],
  "n_imaginary": 0
}
```

**NMR output JSON:**
```json
{
  "source": "orca",
  "nmr_shifts": [
    {"atom_idx": 0, "element": "H", "shielding_ppm": 28.3, "shift_ppm": 3.2}
  ]
}
```

**Step 1: Write failing tests**
```python
# Add to tests/test_chem_qm.py
import tempfile, os

ORCA_IR_OUTPUT = """
$vibrational_frequencies
7
0: 0.0
1: 0.0
2: 0.0
3: 0.0
4: 0.0
5: 0.0
6: 1750.32
$end

$ir_spectrum
7
0: 0.0
1: 0.0
2: 0.0
3: 0.0
4: 0.0
5: 0.0
6: 234.5
$end
"""

def test_parse_ir_frequencies(tmp_path):
    log = tmp_path / "freq.log"
    log.write_text(ORCA_IR_OUTPUT)
    result = run_script("chem_qm.py", ["--parse", str(log), "--parse-ir"])
    assert "frequencies_cm1" in result
    freqs = result["frequencies_cm1"]
    assert any(abs(f["freq_cm1"] - 1750.32) < 0.1 for f in freqs)

def test_parse_ir_no_imaginary(tmp_path):
    log = tmp_path / "freq.log"
    log.write_text(ORCA_IR_OUTPUT)
    result = run_script("chem_qm.py", ["--parse", str(log), "--parse-ir"])
    assert result["n_imaginary"] == 0
```

**Step 2: Run to verify FAIL**
```bash
/home/de/Bureau/Pepflex/venv/bin/python -m pytest tests/test_chem_qm.py::test_parse_ir_frequencies -v
```

**Step 3: Implement in `scripts/chem_qm.py`**

Read the current `chem_qm.py` first. Add new parsing functions:

```python
def parse_orca_ir(text: str) -> dict:
    """Parse ORCA vibrational frequencies and IR intensities."""
    import re
    freqs, intensities = [], []

    # Parse $vibrational_frequencies block
    vib_match = re.search(
        r'\$vibrational_frequencies\s+\d+\s+([\d\s.eE+\-:]+?)\$end',
        text, re.DOTALL)
    if vib_match:
        for line in vib_match.group(1).strip().splitlines():
            parts = line.strip().split(':')
            if len(parts) == 2:
                freqs.append(float(parts[1].strip()))

    # Parse $ir_spectrum block
    ir_match = re.search(
        r'\$ir_spectrum\s+\d+\s+([\d\s.eE+\-:]+?)\$end',
        text, re.DOTALL)
    if ir_match:
        for line in ir_match.group(1).strip().splitlines():
            parts = line.strip().split(':')
            if len(parts) == 2:
                intensities.append(float(parts[1].strip()))

    # Pair freqs with intensities, skip first 6 (translations/rotations)
    modes = []
    n_imaginary = 0
    for i, freq in enumerate(freqs):
        if i < 6:
            continue
        intensity = intensities[i] if i < len(intensities) else 0.0
        if freq < 0:
            n_imaginary += 1
        modes.append({
            "mode": i,
            "freq_cm1": round(freq, 2),
            "intensity": round(intensity, 2),
        })

    return {
        "source": "orca",
        "frequencies_cm1": modes,
        "n_imaginary": n_imaginary,
        "n_modes": len(modes),
    }
```

In `main()`, add to the `--parse` branch:
```python
if args.parse_ir:
    result = parse_orca_ir(text)
elif args.parse_nmr:
    result = parse_orca_nmr(text)
```

Add `--parse-ir` and `--parse-nmr` flags to the parser:
```python
parser.add_argument("--parse-ir", action="store_true", dest="parse_ir")
parser.add_argument("--parse-nmr", action="store_true", dest="parse_nmr")
```

**Step 4: Run tests**
```bash
/home/de/Bureau/Pepflex/venv/bin/python -m pytest tests/test_chem_qm.py -v
```

**Step 5: Commit**
```bash
git add scripts/chem_qm.py tests/test_chem_qm.py
git commit -m "feat: extend chem_qm.py — parse IR frequencies from ORCA output"
```

---

## Task 3: `chem_metabolism.py` — CYP450 metabolic soft spots

**Files:**
- Create: `scripts/chem_metabolism.py`
- Test: `tests/test_chem_metabolism.py`

**Interface:**
```
python chem_metabolism.py --smiles "CC(=O)Oc1ccccc1C(=O)O"
```

**Output:**
```json
{
  "smiles": "...",
  "metabolic_sites": [
    {"enzyme": "CYP3A4", "pattern": "aromatic_hydroxylation",
     "atom_indices": [5, 6, 7, 8, 9, 10], "risk": "medium"}
  ],
  "n_sites": 2,
  "overall_risk": "medium"
}
```

**CYP450 SMARTS rules to implement:**

```python
METABOLISM_RULES = [
    # CYP3A4
    ("CYP3A4", "aromatic_hydroxylation",   "c1ccccc1",              "medium"),
    ("CYP3A4", "aliphatic_hydroxylation",  "[CH2][CH2][CH3]",       "low"),
    ("CYP3A4", "n_dealkylation",           "[NX3][CH3]",            "high"),
    ("CYP3A4", "o_dealkylation",           "[OX2][CH3]",            "medium"),
    # CYP2D6
    ("CYP2D6", "aromatic_hydroxylation",   "c1cccnc1",              "medium"),
    ("CYP2D6", "n_dealkylation",           "[NX3]CC",               "medium"),
    # CYP2C9
    ("CYP2C9", "aromatic_hydroxylation",   "c1ccc(F)cc1",           "low"),
    ("CYP2C9", "sulfoxidation",            "[SX2][CH3]",            "medium"),
    # CYP1A2
    ("CYP1A2", "n_hydroxylation",          "[NX3;H2]c",             "high"),
    ("CYP1A2", "aromatic_hydroxylation",   "c1cccc2ccccc12",        "medium"),
    # Phase II
    ("UGT",    "glucuronidation",          "[OX2H1;!$(OC=O)]",      "medium"),
    ("SULT",   "sulfation",               "[OX2H1;!$(OC=O)]",      "low"),
]
```

**Tests:**
```python
from tests.conftest import run_script, ASPIRIN_SMILES, CAFFEINE_SMILES

SCRIPT = "chem_metabolism.py"

def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    for f in ["smiles", "metabolic_sites", "n_sites", "overall_risk"]:
        assert f in result

def test_aspirin_has_aromatic_site():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES])
    patterns = [s["pattern"] for s in result["metabolic_sites"]]
    assert "aromatic_hydroxylation" in patterns

def test_n_methyl_caffeine():
    # Caffeine has N-CH3 groups → N-dealkylation expected
    result = run_script(SCRIPT, ["--smiles", CAFFEINE_SMILES])
    patterns = [s["pattern"] for s in result["metabolic_sites"]]
    assert "n_dealkylation" in patterns

def test_no_sites_simple():
    # Ethane: no metabolic alerts
    result = run_script(SCRIPT, ["--smiles", "CC"])
    assert result["n_sites"] == 0
    assert result["overall_risk"] == "none"
```

**Step 5: Commit**
```bash
git add scripts/chem_metabolism.py tests/test_chem_metabolism.py
git commit -m "feat: add chem_metabolism.py — CYP450 soft spot prediction"
```

---

## Task 4: `chem_diversity.py` — MaxMin diversity selection

**Files:**
- Create: `scripts/chem_diversity.py`
- Test: `tests/test_chem_diversity.py`

**Interface:**
```
python chem_diversity.py --input library.smi --n 10 [--fingerprint morgan]
```

**Output:**
```json
{
  "n_requested": 10,
  "n_selected": 10,
  "n_library": 100,
  "selected": [
    {"idx": 0, "name": "mol_0", "smiles": "..."},
    ...
  ]
}
```

**Algorithm:** MaxMin — start with the molecule most different from all others, iteratively add the molecule maximally dissimilar to the current selection.

```python
def maxmin_select(mols_fps: list, n: int, seed_idx: int = 0) -> list[int]:
    """Return indices of n maximally diverse molecules."""
    from rdkit import DataStructs
    selected = [seed_idx]
    remaining = list(range(len(mols_fps)))
    remaining.remove(seed_idx)

    while len(selected) < n and remaining:
        max_min_sim = -1
        best_idx = remaining[0]
        for idx in remaining:
            # Min similarity to already selected
            min_sim = min(
                DataStructs.TanimotoSimilarity(mols_fps[idx], mols_fps[s])
                for s in selected)
            if min_sim > max_min_sim:
                max_min_sim = min_sim
                best_idx = idx
        selected.append(best_idx)
        remaining.remove(best_idx)
    return selected
```

**Tests:**
```python
from tests.conftest import ASPIRIN_SMILES, CAFFEINE_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT
import json, subprocess, sys

SCRIPT = "chem_diversity.py"

def make_library(tmp_path, mols):
    p = tmp_path / "lib.smi"
    p.write_text("\n".join(f"{s} {n}" for s, n in mols))
    return str(p)

def run_diversity(args):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT)] + args,
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)

def test_selects_correct_n(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
        ("CC(=O)Oc1ccccc1", "phenyl_acetate"),
        ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", "ibuprofen_like"),
    ])
    result = run_diversity(["--input", lib, "--n", "2"])
    assert result["n_selected"] == 2
    assert len(result["selected"]) == 2

def test_no_duplicates(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"), (CAFFEINE_SMILES, "caffeine"),
        ("CC", "ethane"), ("CCC", "propane"),
    ])
    result = run_diversity(["--input", lib, "--n", "3"])
    indices = [s["idx"] for s in result["selected"]]
    assert len(indices) == len(set(indices))

def test_output_fields(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"), (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_diversity(["--input", lib, "--n", "2"])
    for f in ["n_requested", "n_selected", "n_library", "selected"]:
        assert f in result
    for item in result["selected"]:
        for k in ["idx", "name", "smiles"]:
            assert k in item
```

**Step 5: Commit**
```bash
git add scripts/chem_diversity.py tests/test_chem_diversity.py
git commit -m "feat: add chem_diversity.py — MaxMin diversity selection"
```

---

## Task 5: Update config and memory

**Files:**
- Modify: `config/CLAUDE.md` — add 3 new script rows to table
- Modify: `~/.claude/projects/-home-de-Bureau-ALKYL/memory/MEMORY.md`

**Additions to config table:**
```markdown
| Estimer l'état de protonation à pH donné | `chem_pka.py --smiles SMILES --ph 7.4` |
| Parser fréquences IR depuis output ORCA | `chem_qm.py --parse output.log --parse-ir` |
| Sites de métabolisme CYP450 | `chem_metabolism.py --smiles SMILES` |
| Sélection de diversité MaxMin | `chem_diversity.py --input lib.smi --n 50` |
```

**Commit:**
```bash
git add config/CLAUDE.md
git commit -m "docs: update ALKYL config with PhD-level scripts"
```

---

## Final validation

Run full test suite to confirm nothing broken:
```bash
/home/de/Bureau/Pepflex/venv/bin/python -m pytest tests/ -m "not network" -v
```

Expected: all existing tests pass + new tests pass.

---

## Notes for implementer

- **pKa script**: The `dominant_form` uses RDKit's `Reionizer`/`Uncharger` as a heuristic — it won't be chemically perfect but will be directionally correct. The `ionization_groups` with pKa ranges is the key output for Claude's reasoning.
- **IR parsing**: ORCA format shown uses `$vibrational_frequencies` / `$ir_spectrum` blocks (ORCA 5.x simple_output). Real ORCA `.out` files use a different format — test with `test_data/` fixture files if available, otherwise use synthetic test strings (as shown).
- **Metabolism SMARTS**: Some SMARTS will over-match (e.g., `[OX2H1]` matches all alcohols). This is intentional — false positives are safer than false negatives for metabolic risk assessment.
- **MaxMin performance**: O(n²) — fine for libraries up to ~10k molecules. For larger libraries, consider pre-computing all pairwise similarities.
