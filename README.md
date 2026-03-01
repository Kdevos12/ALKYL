# ALKYL

**Computational Chemistry Assistant** — a Claude Code plugin specialized for molecular modeling, quantum chemistry, and cheminformatics.

## What it does

ALKYL augments Claude Code with:
- A computational chemistry identity and context injected into every session
- Chemistry-specific skills (RDKit, cheminformatics workflows, SMILES reference)
- A dedicated `alkyl` command with splash screen

## Requirements

- [Claude Code](https://claude.ai/code) installed and authenticated
- Bash

## Install

```bash
git clone https://github.com/YOUR_USERNAME/alkyl
cd alkyl
bash install.sh
```

Then run:

```bash
alkyl
```

## Usage

`alkyl` is a drop-in replacement for `claude` with chemistry context pre-loaded:

```bash
alkyl                          # interactive session
alkyl -p "Draw the SMILES for caffeine"   # non-interactive
alkyl --resume <session-id>    # resume a session
```

All Claude Code flags work as-is.

## Project structure

```
ALKYL/
├── alkyl                   # CLI wrapper
├── install.sh              # one-command install
├── config/
│   ├── CLAUDE.md           # chemistry system context
│   └── settings.json       # hooks config
├── hooks/
│   └── banner.sh           # splash screen
└── skills/
    ├── rdkit.md            # RDKit cheminformatics guide
    └── cheminformatics.md  # molecular representations & workflows
```

## Uninstall

```bash
rm -rf ~/.alkyl
rm ~/.local/bin/alkyl
```

---

*Built on [Claude Code](https://claude.ai/code) by Anthropic.*
