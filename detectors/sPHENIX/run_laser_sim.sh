#!/bin/bash
set -e  # Exit on error

# Setup environment
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.494
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/

# Args: process index and number of events per job (default 100)
process_id=${1?Error: process ID is not given}
nevents=${2:-100}

# Working directory containing the macro
WORK_DIR="/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX"
cd "$WORK_DIR"

padded_id=$(printf "%05d" "$process_id")
outroot="laser_${padded_id}.root"

echo "Running laser sim: job=$process_id nevents=$nevents out=$outroot"

# Run the laser simulation macro; output file name is used to derive DNL filename
root.exe -l -b -q "Fun4All_G4_sPHENIX_LineLaser_sim.C(${nevents}, \"\", \"${outroot}\")"

echo "Done: ${outroot} and derived DNL file should be created"
