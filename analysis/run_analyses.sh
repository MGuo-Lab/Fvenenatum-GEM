#!/usr/bin/env bash
# Run every analysis behind the figures and tables, in dependency order, then draw the figures.
# Run from analysis/ with the environment of the README active:   bash run_analyses.sh
#
# Steps marked "raw" read the raw instrument files (HPLC .lcd data files, PreSens exports), which are not part of
# this repository. They run only when FBAVAL_RAW_DIR points to them; otherwise the committed result files of those
# steps are used by the later steps.
set -euo pipefail
cd "$(dirname "$0")"
export MPLBACKEND=Agg
PY=${PYTHON:-python}

have_raw() { [ -n "${FBAVAL_RAW_DIR:-}" ] && [ -d "${FBAVAL_RAW_DIR}/lcd" ]; }
run() { echo "[$(date -u +%H:%M:%S)] $1"; "$PY" "scripts/$1.py" > /dev/null; }

have_raw && run flask_hplc || echo "skip flask_hplc (raw)"
run flask_products
run model_versions_qc
have_raw && run glycerol_diagnosis || echo "skip glycerol_diagnosis (raw)"
have_raw && run peak_identity || echo "skip peak_identity (raw)"
run yield_sensitivity
run glycerol_oxygen_test
run curation_stage_effects
run energy_rules_effect
run lag_x0_fit
run flask_dynamics
run methods_checks
run uptake_threshold_check
run overprediction_tests
run km_sensitivity
run run_memote
run ph_validation
run paper_figures
run lag_x0_figure
echo "done"
