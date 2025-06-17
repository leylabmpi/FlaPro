#!/bin/bash

# -----------------------------
# Usage:
#   ./run/sync.sh notebooks/ProjectA.ipynb [--update] [--force] [--debug]
# -----------------------------

# Resolve Git root
GIT_ROOT=$(git rev-parse --show-toplevel)

# Use full path to script
SCRIPT_PATH="$GIT_ROOT/notebooks/differential_analysis/notebook_generator/R/reverse_sync.R"

NOTEBOOK_PATH="$1"
shift
EXTRA_ARGS="$@"

if [[ -z "$NOTEBOOK_PATH" ]]; then
  echo "❌ Usage: ./run/sync.sh <notebook.ipynb> [--update] [--force] [--debug]"
  exit 1
fi

# Normalize notebook path to absolute
ABS_NOTEBOOK_PATH=$(realpath --canonicalize-missing "$NOTEBOOK_PATH")
echo "🔄 Normalized notebook path: $ABS_NOTEBOOK_PATH"

cd "$GIT_ROOT" || {
  echo "❌ Failed to cd into Git root at $GIT_ROOT"
  exit 1
}

#  and normalize the notebook path and then back to relative to Git root
REL_NOTEBOOK_PATH="${ABS_NOTEBOOK_PATH#$GIT_ROOT/}"
echo "🔄 Relative notebook path: $REL_NOTEBOOK_PATH"

if [[ ! -f "$REL_NOTEBOOK_PATH" ]]; then
  echo "❌ Notebook not found: $REL_NOTEBOOK_PATH"
  exit 1
fi

echo "🔄 Running reverse sync on: $REL_NOTEBOOK_PATH"
Rscript "$SCRIPT_PATH" "$REL_NOTEBOOK_PATH" $EXTRA_ARGS
