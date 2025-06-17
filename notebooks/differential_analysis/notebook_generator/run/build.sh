#!/bin/bash

# -----------------------------
# Usage:
#   ./run/build.sh                    # uses default config
#   ./run/build.sh configs/dev.yaml   # custom config
# -----------------------------

CONFIG_PATH="$1"
CONFIG_PATH=${CONFIG_PATH:-/notebooks/differential_analysis/notebook_generator/configs/config.yaml}

# Get absolute Git root
GIT_ROOT=$(git rev-parse --show-toplevel)
cd "$GIT_ROOT" || exit 1

# Use full path to script
SCRIPT_PATH="$GIT_ROOT/notebooks/differential_analysis/notebook_generator/R/build_notebook.R"
CONFIG_PATH_ABS="$GIT_ROOT/$CONFIG_PATH"

if [[ ! -f "$SCRIPT_PATH" ]]; then
  echo "❌ build_notebook.R not found at: $SCRIPT_PATH"
  exit 1
fi

if [[ ! -f "$CONFIG_PATH_ABS" ]]; then
  echo "❌ Config not found: $CONFIG_PATH_ABS"
  exit 1
fi

echo "📦 Building notebooks from config: $CONFIG_PATH_ABS"
Rscript "$SCRIPT_PATH" "$CONFIG_PATH_ABS"
