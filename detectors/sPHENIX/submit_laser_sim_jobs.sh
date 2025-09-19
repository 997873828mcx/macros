#!/bin/bash

# Working directory
WORK_DIR="/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX"

# Create logs directory if it doesn't exist
mkdir -p "$WORK_DIR/log/tmp"

# Make sure the run script is executable
chmod +x "$WORK_DIR/run_laser_sim.sh"

# Quick sanity checks
if [ ! -f "$WORK_DIR/Fun4All_G4_sPHENIX_LineLaser_sim.C" ]; then
  echo "Error: Fun4All_G4_sPHENIX_LineLaser_sim.C not found in $WORK_DIR" >&2
  exit 1
fi

if [ ! -f "$WORK_DIR/laser_sim.job" ]; then
  echo "Error: laser_sim.job not found in $WORK_DIR" >&2
  exit 1
fi

# Submit condor jobs
condor_submit "$WORK_DIR/laser_sim.job"

echo "Jobs submitted. Monitor with: condor_q"
echo "Logs in: $WORK_DIR/log/"
