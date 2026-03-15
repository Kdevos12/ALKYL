#!/usr/bin/env bash
# ALKYL — Unified management script
# Usage:
#   bash alkyl.sh install
#   bash alkyl.sh repair
#   bash alkyl.sh uninstall
#   bash alkyl.sh setup-key perplexity <API_KEY>
#   bash alkyl.sh status

set -e

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLAUDE_MD="$HOME/.claude/CLAUDE.md"
SETTINGS="$HOME/.claude/settings.json"
MARKER_START="<!-- ALKYL-START -->"
MARKER_END="<!-- ALKYL-END -->"

# ── Colors ─────────────────────────────────────────────────────────────────
BLUE='\033[38;2;31;81;255m'
BOLD='\033[1m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
DIM='\033[2m'
RESET='\033[0m'

# ── Logo ───────────────────────────────────────────────────────────────────
print_logo() {
    printf "\n${BLUE}${BOLD}"
    cat << 'LOGO'
 ░▒▓██████▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓████████▓▒░▒▓█▓▒░      ░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓████████▓▒░
LOGO
    printf "${RESET}"
    printf "  ${DIM}Computational Chemistry · Claude Code Plugin${RESET}\n\n"
}

# ── Helpers ────────────────────────────────────────────────────────────────
is_installed() {
    grep -q "$MARKER_START" "$CLAUDE_MD" 2>/dev/null
}

remove_block() {
    if is_installed; then
        sed -i "/$MARKER_START/,/$MARKER_END/d" "$CLAUDE_MD"
    fi
}

inject_block() {
    mkdir -p "$(dirname "$CLAUDE_MD")"
    touch "$CLAUDE_MD"
    remove_block
    {
        echo ""
        echo "$MARKER_START"
        cat "$REPO_DIR/config/CLAUDE.md"
        echo "$MARKER_END"
    } >> "$CLAUDE_MD"
    # Resolve scripts path placeholder
    SCRIPTS_ABS="$REPO_DIR/scripts"
    sed -i "s|ALKYL_SCRIPTS_PATH|${SCRIPTS_ABS}|g" "$CLAUDE_MD"
}

# ── Commands ───────────────────────────────────────────────────────────────

link_skills() {
    local skills_dst="$HOME/.claude/skills"
    mkdir -p "$skills_dst"
    for skill_dir in "$REPO_DIR/skills"/*/; do
        [ -f "${skill_dir}SKILL.md" ] || continue
        ln -sfn "$skill_dir" "$skills_dst/$(basename "$skill_dir")"
    done
}

unlink_skills() {
    local skills_dst="$HOME/.claude/skills"
    for skill_dir in "$REPO_DIR/skills"/*/; do
        [ -f "${skill_dir}SKILL.md" ] || continue
        local name
        name=$(basename "$skill_dir")
        # Only remove if it's our symlink
        if [ -L "$skills_dst/$name" ] && [ "$(readlink "$skills_dst/$name")" = "$skill_dir" ]; then
            rm "$skills_dst/$name"
        fi
    done
}

register_plugin() {
    mkdir -p "$(dirname "$SETTINGS")"
    [ -f "$SETTINGS" ] || echo '{"plugins":{}}' > "$SETTINGS"
    link_skills
    python3 - "$SETTINGS" "$REPO_DIR" "$HOME/.claude/plugins/installed_plugins.json" <<'PY'
import json, sys, datetime
path, repo, plugins_path = sys.argv[1], sys.argv[2], sys.argv[3]
with open(path) as f:
    d = json.load(f)
now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")
installed = json.load(open(plugins_path))
installed["plugins"]["alkyl@local"] = [{
    "scope": "user",
    "installPath": repo,
    "version": "1.0.0",
    "installedAt": now,
    "lastUpdated": now,
    "gitCommitSha": "local"
}]
with open(plugins_path, "w") as f:
    json.dump(installed, f, indent=2)
# Enable the plugin in settings.json (required for skill discovery)
d.setdefault("enabledPlugins", {})["alkyl@local"] = True
with open(path, "w") as f:
    json.dump(d, f, indent=2, ensure_ascii=False)
PY
}

deregister_plugin() {
    unlink_skills
    python3 - "$HOME/.claude/plugins/installed_plugins.json" "$SETTINGS" <<'PY'
import json, os, sys
plugins_path, settings_path = sys.argv[1], sys.argv[2]
with open(plugins_path) as f:
    d = json.load(f)
d["plugins"].pop("alkyl@local", None)
with open(plugins_path, "w") as f:
    json.dump(d, f, indent=2)
if os.path.exists(settings_path):
    with open(settings_path) as f:
        s = json.load(f)
    s.get("enabledPlugins", {}).pop("alkyl@local", None)
    with open(settings_path, "w") as f:
        json.dump(s, f, indent=2, ensure_ascii=False)
PY
}

cmd_install() {
    print_logo
    if is_installed; then
        printf "  ${YELLOW}⚠ ALKYL already installed. Re-injecting (idempotent)...${RESET}\n"
    fi
    inject_block
    register_plugin
    printf "  ${GREEN}✓ Chemistry context injected into $CLAUDE_MD${RESET}\n"
    printf "  ${GREEN}✓ ALKYL registered as plugin (visible in /plugins)${RESET}\n"
    printf "  ${GREEN}✓ 23 skills discoverable via /skills${RESET}\n\n"
    printf "  ${DIM}Repair:    bash alkyl.sh repair${RESET}\n"
    printf "  ${DIM}Uninstall: bash alkyl.sh uninstall${RESET}\n"
    printf "  ${DIM}API keys:  bash alkyl.sh setup-key perplexity <KEY>${RESET}\n\n"
}

cmd_repair() {
    print_logo
    printf "  ${YELLOW}⟳ Repairing ALKYL — re-injecting config...${RESET}\n"
    inject_block
    printf "  ${GREEN}✓ Config repaired and re-injected${RESET}\n"
    printf "  ${GREEN}✓ Scripts path: $REPO_DIR/scripts${RESET}\n\n"
}

cmd_uninstall() {
    printf "\n  ${RED}Uninstalling ALKYL...${RESET}\n"
    if is_installed; then
        remove_block
        printf "  ${GREEN}✓ ALKYL block removed from $CLAUDE_MD${RESET}\n"
    else
        printf "  ${DIM}ALKYL was not installed in $CLAUDE_MD${RESET}\n"
    fi
    deregister_plugin
    printf "  ${GREEN}✓ ALKYL deregistered from plugin system${RESET}\n"

    # Optional: remove MCP keys from settings.json
    if [ -f "$SETTINGS" ] && grep -q '"perplexity"' "$SETTINGS" 2>/dev/null; then
        printf "  ${YELLOW}Remove Perplexity MCP key from settings.json? (y/N):${RESET} "
        read -r ans
        if [[ "$ans" =~ ^[Yy]$ ]]; then
            python3 - "$SETTINGS" <<'PY'
import json, sys
path = sys.argv[1]
with open(path) as f:
    d = json.load(f)
d.get('mcpServers', {}).pop('perplexity', None)
with open(path, 'w') as f:
    json.dump(d, f, indent=2, ensure_ascii=False)
print("  ✓ Perplexity MCP key removed")
PY
        fi
    fi
    printf "\n"
}

cmd_status() {
    printf "\n  ${BOLD}ALKYL Status${RESET}\n"
    printf "  ─────────────────────────────\n"

    if is_installed; then
        printf "  ${GREEN}✓ Installed${RESET} — $CLAUDE_MD\n"
        # Show block line count
        LINES=$(sed -n "/$MARKER_START/,/$MARKER_END/p" "$CLAUDE_MD" | wc -l)
        printf "  ${DIM}  Block size: $LINES lines${RESET}\n"
    else
        printf "  ${RED}✗ Not installed${RESET}\n"
    fi

    printf "  ${DIM}Repo:    $REPO_DIR${RESET}\n"
    printf "  ${DIM}Scripts: $REPO_DIR/scripts${RESET}\n"

    # Count skills
    N_SKILLS=$(find "$REPO_DIR/skills" -name "SKILL.md" 2>/dev/null | wc -l)
    printf "  ${DIM}Skills:  $N_SKILLS loaded${RESET}\n"

    # MCP keys
    if [ -f "$SETTINGS" ]; then
        printf "\n  ${BOLD}MCP keys in settings.json:${RESET}\n"
        python3 - "$SETTINGS" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    d = json.load(f)
servers = d.get('mcpServers', {})
if not servers:
    print("  (none configured)")
else:
    for name, cfg in servers.items():
        env = cfg.get('env', {})
        keys = [k for k in env if 'KEY' in k or 'TOKEN' in k or 'SECRET' in k]
        masked = {k: env[k][:8] + '...' for k in keys}
        print(f"  • {name}: {masked if masked else '(no key)'}")
PY
    else
        printf "  ${DIM}(settings.json not found)${RESET}\n"
    fi
    printf "\n"
}

cmd_venv() {
    VENV_DIR="$REPO_DIR/.venv"
    printf "\n  ${BOLD}Setting up ALKYL virtual environment${RESET}\n"
    printf "  ${DIM}Location: $VENV_DIR${RESET}\n\n"

    if [ -d "$VENV_DIR" ]; then
        printf "  ${YELLOW}⚠ .venv already exists. Reinstalling dependencies...${RESET}\n"
    else
        python3 -m venv "$VENV_DIR"
        printf "  ${GREEN}✓ Virtual environment created${RESET}\n"
    fi

    "$VENV_DIR/bin/pip" install --quiet --upgrade pip
    "$VENV_DIR/bin/pip" install --quiet "rdkit>=2023.9.1" pytest

    printf "  ${GREEN}✓ Dependencies installed (rdkit, pytest)${RESET}\n\n"
    printf "  ${BOLD}Run scripts:${RESET}\n"
    printf "  ${DIM}  $VENV_DIR/bin/python scripts/chem_props.py --smiles 'CCO'${RESET}\n\n"
    printf "  ${BOLD}Run tests:${RESET}\n"
    printf "  ${DIM}  $VENV_DIR/bin/python -m pytest tests/ -m 'not network' -v${RESET}\n\n"
}

cmd_setup_key() {
    local SERVICE="${1:-}"
    local API_KEY="${2:-}"

    if [ -z "$SERVICE" ] || [ -z "$API_KEY" ]; then
        printf "\n  ${RED}Usage: bash alkyl.sh setup-key <service> <API_KEY>${RESET}\n"
        printf "  Supported services: perplexity\n\n"
        exit 1
    fi

    mkdir -p "$(dirname "$SETTINGS")"
    [ -f "$SETTINGS" ] || echo "{}" > "$SETTINGS"

    case "$SERVICE" in
        perplexity)
            if [[ "$API_KEY" != pplx-* ]]; then
                printf "  ${YELLOW}Warning: key doesn't start with 'pplx-'. Continue? (y/N):${RESET} "
                read -r ans
                [[ "$ans" =~ ^[Yy]$ ]] || exit 1
            fi
            python3 - "$SETTINGS" "$API_KEY" <<'PY'
import json, sys
path, api_key = sys.argv[1], sys.argv[2]
with open(path) as f:
    d = json.load(f)
d.setdefault('mcpServers', {})['perplexity'] = {
    'command': 'npx',
    'args': ['-y', '@perplexity-ai/mcp-server'],
    'env': {'PERPLEXITY_API_KEY': api_key}
}
with open(path, 'w') as f:
    json.dump(d, f, indent=2, ensure_ascii=False)
PY
            printf "  ${GREEN}✓ Perplexity MCP configured in $SETTINGS${RESET}\n"
            printf "  ${DIM}Restart Claude Code for changes to take effect.${RESET}\n\n"
            ;;
        *)
            printf "  ${RED}Unknown service: $SERVICE${RESET}\n"
            printf "  Supported: perplexity\n\n"
            exit 1
            ;;
    esac
}

# ── Dispatch ───────────────────────────────────────────────────────────────
CMD="${1:-}"
case "$CMD" in
    install)   cmd_install ;;
    repair)    cmd_repair ;;
    uninstall) cmd_uninstall ;;
    status)    cmd_status ;;
    venv)      cmd_venv ;;
    setup-key) cmd_setup_key "${2:-}" "${3:-}" ;;
    *)
        printf "\n  ${BOLD}ALKYL — Computational Chemistry Plugin${RESET}\n\n"
        printf "  Usage: bash alkyl.sh <command>\n\n"
        printf "  Commands:\n"
        printf "    ${GREEN}install${RESET}                        Install ALKYL into ~/.claude/CLAUDE.md\n"
        printf "    ${GREEN}venv${RESET}                           Create .venv with RDKit (for scripts & tests)\n"
        printf "    ${GREEN}repair${RESET}                         Force re-inject config (fixes corruption)\n"
        printf "    ${GREEN}uninstall${RESET}                      Remove ALKYL from ~/.claude/CLAUDE.md\n"
        printf "    ${GREEN}status${RESET}                         Show installation status and MCP keys\n"
        printf "    ${GREEN}setup-key perplexity <KEY>${RESET}     Configure Perplexity MCP API key\n\n"
        exit 1
        ;;
esac
