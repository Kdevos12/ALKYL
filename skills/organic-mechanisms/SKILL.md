---
name: organic-mechanisms
description: Use when reasoning about polar organic reaction mechanisms, arrow pushing, SN1/SN2/E1/E2, nucleophile/electrophile identification, Zaitsev/Hofmann selectivity, HSAB rules, or retrosynthesis. Applies the EASE framework step-by-step.
---

# Organic Mechanisms — EASE Framework

**Source**: AceOrganicChem.com *Ace Organic Chemistry Mechanisms with E.A.S.E.* (2013) + standard references (Clayden, March).

## When to Use
Any polar organic mechanism or synthesis problem. Apply EASE iteratively — one step at a time, rinse and repeat.

## The Four Steps (Quick Reference)

```
E — Electrophile   Identify E+ and Nu−. Use resonance if unclear.
A — Acid/Base      Strong acid/base present? Move proton FIRST before continuing.
S — Sterics        Bulky group on Nu or E+? Reassess: block / carbocation / overcome.
E — Electron Flow  Nu− → E+. Draw curved arrows. Check valence. Done?
                   If not, repeat from Step 1 with new intermediates.
```

## Reference Files

| File | Content |
|------|---------|
| `ease-framework.md` | Full step-by-step logic, decision rules, worked examples |
| `electrophiles-nucleophiles.md` | Identification rules, pKa table, nucleophilicity scales |
| `sterics-substitution.md` | SN1/SN2/E1/E2 decision tree, steric groups, Zaitsev/Hofmann |
| `electron-flow.md` | Arrow pushing, HSAB (1,2 vs 1,4), intra vs inter, stereochemistry |
| `retrosynthesis.md` | Disconnections, synthons, FGI, toolbox (C–C / C–X / FGI) |
| `method-limits.md` | When EASE fails: RedOx table, radicals, organometallics |

## Quick Routing

**"What is the mechanism?"** → `ease-framework.md`
**"Which product forms?"** → `sterics-substitution.md` (SN/E) or `electron-flow.md` (HSAB)
**"How do I make X from Y?"** → `retrosynthesis.md`
**"Reagent does something weird"** → `method-limits.md`
**"I can't find the nucleophile"** → `electrophiles-nucleophiles.md` (draw resonance structures)

## Core Rule
> Electrons *always* flow from nucleophile to electrophile.
> Arrows point to where electrons are *going*.
> If you cannot identify a nucleophile, draw all resonance structures first.
