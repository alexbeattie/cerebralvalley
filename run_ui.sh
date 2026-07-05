#!/usr/bin/env bash
#
# run_ui.sh — launch the causality-review web UI.
#
# Finds the Python environment (shared venv from the main worktree if present),
# then starts the server. The server auto-falls back to the next free port if
# the requested one is busy, so "Address already in use" won't stop you.
#
# Usage:
#   ./run_ui.sh              # http://localhost:8000 (or next free port)
#   ./run_ui.sh --port 8080
#   ./run_ui.sh --no-open    # don't auto-open a browser

set -u

if [[ -x ".venv/bin/python" ]]; then
  PY=".venv/bin/python"
elif [[ -x "../biohack/.venv/bin/python" ]]; then
  PY="../biohack/.venv/bin/python"   # shared venv from the main worktree
else
  PY="python3"
fi

exec "$PY" -m causality_review.webapp "$@"
