#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUNNER="$SCRIPT_DIR/../run_phase2b.R"

configs=(
  "$SCRIPT_DIR/config_cv_balanced.yaml"
  "$SCRIPT_DIR/config_cv_soncin_1st_only.yaml"
  "$SCRIPT_DIR/config_cv_soncin_2nd_only.yaml"
  "$SCRIPT_DIR/config_cv_mikheev_1st_only.yaml"
  "$SCRIPT_DIR/config_cv_mikheev_2nd_only.yaml"
)

for cfg in "${configs[@]}"; do
  echo "============================================================"
  echo "Running: $(basename "$cfg")"
  echo "============================================================"
  Rscript "$RUNNER" --config="$cfg" --no_archive
  echo ""
done

echo "All cross-validation runs complete."
echo "Run the comparison script:"
echo "  Rscript $SCRIPT_DIR/compare_cv_results.R"
