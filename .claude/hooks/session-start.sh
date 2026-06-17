#!/bin/bash
# SessionStart hook: provision R + the refund package's dependencies (and fastFMM)
# for Claude Code on the web. Idempotent and safe to re-run.
set -euo pipefail

# Only run in the remote (web) sandbox; local machines manage their own R setup.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

SUDO=""
if [ "$(id -u)" -ne 0 ]; then
  command -v sudo >/dev/null 2>&1 && SUDO="sudo"
fi

export DEBIAN_FRONTEND=noninteractive

# Refresh apt lists. Some third-party PPAs in the base image are blocked by the
# sandbox network policy; tolerate their failures since the main archive (which
# provides r-base and the r-cran-* binaries) still updates.
$SUDO apt-get update -qq || true

# Install R if it is not already present.
if ! command -v Rscript >/dev/null 2>&1; then
  echo "[session-start] installing R..."
  $SUDO apt-get install -y r-base r-base-dev
fi

# Install refund's dependencies + fastFMM (prefers apt binaries, falls back to
# the GitHub CRAN mirror because CRAN itself is blocked in the web sandbox).
echo "[session-start] installing R package dependencies..."
$SUDO Rscript "${CLAUDE_PROJECT_DIR}/.claude/hooks/install-r-deps.R"

echo "[session-start] done."
