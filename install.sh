#!/usr/bin/env bash
# ALKYL — Install script
# Appends chemistry context to ~/.claude/CLAUDE.md
# Usage: bash install.sh
set -e

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLAUDE_MD="$HOME/.claude/CLAUDE.md"
MARKER_START="<!-- ALKYL-START -->"
MARKER_END="<!-- ALKYL-END -->"

# Neon blue #1F51FF
BLUE='\033[38;2;31;81;255m'
BOLD='\033[1m'
GREEN='\033[0;32m'
DIM='\033[2m'
RESET='\033[0m'

printf "\n${BLUE}${BOLD}"
cat << 'EOF'
 ░▒▓██████▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓████████▓▒░▒▓█▓▒░      ░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓████████▓▒░
EOF
printf "${RESET}"
printf "  ${DIM}Computational Chemistry · Claude Code Plugin${RESET}\n\n"

# --- Inject into ~/.claude/CLAUDE.md ---
mkdir -p "$(dirname "$CLAUDE_MD")"
touch "$CLAUDE_MD"

# Remove existing ALKYL block if present (idempotent install)
if grep -q "$MARKER_START" "$CLAUDE_MD" 2>/dev/null; then
    sed -i "/$MARKER_START/,/$MARKER_END/d" "$CLAUDE_MD"
fi

# Append chemistry context
{
    echo ""
    echo "$MARKER_START"
    cat "$REPO_DIR/config/CLAUDE.md"
    echo "$MARKER_END"
} >> "$CLAUDE_MD"

printf "  ${GREEN}✓ Chemistry context injected into $CLAUDE_MD${RESET}\n"
printf "  ${GREEN}✓ ALKYL active in all future claude sessions${RESET}\n\n"
printf "  ${DIM}Uninstall: bash uninstall.sh${RESET}\n\n"
