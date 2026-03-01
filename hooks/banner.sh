#!/usr/bin/env bash
# ALKYL — Splash screen (shown briefly before Claude Code TUI renders)

CYAN='\033[0;36m'
BOLD='\033[1m'
DIM='\033[2m'
RESET='\033[0m'

printf "\n${CYAN}${BOLD}"
cat << 'EOF'
    ___   ___       __ __ __  ___
   /   | / / /____ / // // / /  /
  / /| |/ / //_  // // // / /  /
 / ___ / / / / // // // / /  /___
/_/  |_/_/_/ /___//_//_/ /_____/
EOF
printf "${RESET}"
printf "  ${DIM}Computational Chemistry · Powered by Claude Code${RESET}\n\n"
