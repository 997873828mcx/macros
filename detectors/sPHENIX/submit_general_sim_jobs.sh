#!/bin/bash
set -euo pipefail

# Working directory
WORK_DIR="/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX"
CAMPAIGN="${1:-pp_minbias_truth_100}"
HEPMC_FILE="${2:-$WORK_DIR/output/pp_minbias_test.hepmc}"
TOTAL_EVENTS="${3:-100}"
EVENTS_PER_JOB="${4:-10}"

# Create logs directory if it doesn't exist
mkdir -p "$WORK_DIR/log/tmp"
mkdir -p "$WORK_DIR/output"

# Make sure the run script is executable
chmod +x "$WORK_DIR/run_general_sim.sh"

if [ ! -f "$WORK_DIR/Fun4All_G4_sPHENIX.C" ]; then
  echo "Error: Fun4All_G4_sPHENIX.C not found in $WORK_DIR" >&2
  exit 1
fi

if [ ! -f "$WORK_DIR/general_sim.job" ]; then
  echo "Error: general_sim.job not found in $WORK_DIR" >&2
  exit 1
fi

# Submit condor jobs
"$WORK_DIR/submit_general_campaign.sh" "$CAMPAIGN" "$HEPMC_FILE" "$TOTAL_EVENTS" "$EVENTS_PER_JOB"

echo "Jobs submitted. Monitor with: condor_q"
echo "Logs in: $WORK_DIR/log/$CAMPAIGN/"
echo "Outputs in: $WORK_DIR/output/$CAMPAIGN/completed/"
