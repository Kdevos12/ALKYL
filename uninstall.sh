#!/usr/bin/env bash
# ALKYL — Uninstall script
set -e

CLAUDE_MD="$HOME/.claude/CLAUDE.md"
MARKER_START="<!-- ALKYL-START -->"

GREEN='\033[0;32m'
RESET='\033[0m'

if grep -q "$MARKER_START" "$CLAUDE_MD" 2>/dev/null; then
    sed -i "/<!-- ALKYL-START -->/,/<!-- ALKYL-END -->/d" "$CLAUDE_MD"
    printf "  ${GREEN}✓ ALKYL removed from $CLAUDE_MD${RESET}\n"
else
    echo "  ALKYL was not installed."
fi
