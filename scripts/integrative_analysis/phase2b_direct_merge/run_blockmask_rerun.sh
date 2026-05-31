#!/bin/bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

Rscript scripts/integrative_analysis/phase2b_direct_merge/run_phase2b.R --config=scripts/integrative_analysis/phase2b_direct_merge/config_combat/config_phase2b_1_2_restoration_blockmask_imputed_0.yaml 2>&1 | tee phase2b_1_2_blockmask_rerun.out

Rscript scripts/integrative_analysis/phase2b_direct_merge/run_phase2b.R --config=scripts/integrative_analysis/phase2b_direct_merge/config_combat/config_phase2b_2_3_restoration_blockmask_imputed_0.yaml 2>&1 | tee phase2b_2_3_blockmask_rerun.out

Rscript scripts/integrative_analysis/phase2b_direct_merge/build_runs_summary_xlsx.R
