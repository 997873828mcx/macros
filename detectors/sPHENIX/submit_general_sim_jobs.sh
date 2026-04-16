#!/bin/bash

# Working directory
WORK_DIR="/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX"

# Create logs directory if it doesn't exist
mkdir -p "$WORK_DIR/log/default"
mkdir -p "$WORK_DIR/log/tmp/default"
mkdir -p "$WORK_DIR/output"
mkdir -p "$WORK_DIR/output/default/completed"

# Make sure the run script is executable
chmod +x "$WORK_DIR/run_general_sim.sh"

# Quick sanity checks
if [ ! -f "$WORK_DIR/Fun4All_G4_sPHENIX.C" ]; then
  echo "Error: Fun4All_G4_sPHENIX.C not found in $WORK_DIR" >&2
  exit 1
fi

if [ ! -f "$WORK_DIR/general_sim.job" ]; then
  echo "Error: general_sim.job not found in $WORK_DIR" >&2
  exit 1
fi

# Submit condor jobs
condor_submit "$WORK_DIR/general_sim.job"

echo "Jobs submitted. Monitor with: condor_q"
echo "Logs in: $WORK_DIR/log/"
