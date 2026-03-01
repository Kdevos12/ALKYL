#!/usr/bin/env bash
# ALKYL — Install script
# Usage: bash install.sh
set -e

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALKYL_HOME="$HOME/.alkyl"
BIN_DIR="$HOME/.local/bin"

CYAN='\033[0;36m'
BOLD='\033[1m'
GREEN='\033[0;32m'
DIM='\033[2m'
RESET='\033[0m'

echo ""
printf "${CYAN}${BOLD}"
cat << 'EOF'
    ___   ___       __ __ __  ___
   /   | / / /____ / // // / /  /
  / /| |/ / //_  // // // / /  /
 / ___ / / / / // // // / /  /___
/_/  |_/_/_/ /___//_//_/ /_____/
EOF
printf "${RESET}"
printf "  ${DIM}Installing ALKYL...${RESET}\n\n"

# --- Config directory ---
mkdir -p "$ALKYL_HOME/hooks"
mkdir -p "$ALKYL_HOME/config"
mkdir -p "$ALKYL_HOME/skills"

cp "$REPO_DIR/config/CLAUDE.md"     "$ALKYL_HOME/config/CLAUDE.md"
cp "$REPO_DIR/config/settings.json" "$ALKYL_HOME/config/settings.json"

# Fix absolute path in settings.json for this user
sed -i "s|/home/de/.alkyl|$ALKYL_HOME|g" "$ALKYL_HOME/config/settings.json"

cp "$REPO_DIR/hooks/banner.sh"  "$ALKYL_HOME/hooks/banner.sh"
chmod +x "$ALKYL_HOME/hooks/banner.sh"

cp "$REPO_DIR/skills/"*.md "$ALKYL_HOME/skills/" 2>/dev/null || true

# --- alkyl command ---
mkdir -p "$BIN_DIR"
cp "$REPO_DIR/alkyl" "$BIN_DIR/alkyl"

# Fix paths in alkyl wrapper for this user
sed -i "s|/home/de/.alkyl|$ALKYL_HOME|g" "$BIN_DIR/alkyl"
chmod +x "$BIN_DIR/alkyl"

# --- PATH check ---
if ! echo "$PATH" | grep -q "$BIN_DIR"; then
    echo "  ⚠  Add this to your ~/.bashrc or ~/.zshrc:"
    echo "     export PATH=\"\$HOME/.local/bin:\$PATH\""
    echo ""
fi

printf "  ${GREEN}✓ Installed to $ALKYL_HOME${RESET}\n"
printf "  ${GREEN}✓ Command: $BIN_DIR/alkyl${RESET}\n\n"
printf "  Run ${BOLD}alkyl${RESET} to start.\n\n"
