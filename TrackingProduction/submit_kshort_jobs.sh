#!/bin/bash

# Working directory - adjust if needed
WORK_DIR="/sphenix/user/dcxchenxi/develope/macros/TrackingProduction"

# Create logs directory if it doesn't exist
mkdir -p $WORK_DIR/log
mkdir -p $WORK_DIR/log/tmp

# Make sure the run script is executable
chmod +x $WORK_DIR/run_kshort_reco.sh

# Check if the macro file exists in the working directory
if [ ! -f "$WORK_DIR/Fun4All_TrackSeeding.C" ]; then
    echo "Error: Fun4All_TrackSeeding.C not found in $WORK_DIR"
    echo "Please ensure the macro is in the correct location"
    exit 1
fi

# Submit the condor job
condor_submit $WORK_DIR/kshort_reco.job

echo "Jobs submitted. Monitor status with: condor_q"
echo "Log files will be in $WORK_DIR/log/"