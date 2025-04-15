#!/bin/bash
set -e  # Exit on error

# Setup environment variables
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}

# Source your setup scripts
source /opt/sphenix/core/bin/sphenix_setup.sh new
#source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/

# Print environment for debugging if needed
echo "Running as user: $USER"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

# Get job index from argument
process_id=${1?Error: process ID is not given}
offset=${2?Error: offset is not given}

file_index=$((process_id + offset))
echo "Process ID: $process_id, Offset: $offset, File Index: $file_index"

# Construct filename using pattern
padded_id=$(printf "%05d" $file_index)
input_file="DST_TRKR_CLUSTER_run2pp_ana466_2024p012_v001-00053877-${padded_id}.root"
input_dir="/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana466_2024p012_v001/DST_TRKR_CLUSTER/run_00053800_00053900/dst/"
output_base="kshort_${padded_id}"

# Full path to input file
#full_input_path="${input_dir}/${input_file}"

# Check if input file exists
#if [ ! -f "$full_input_path" ]; then
#    echo "ERROR: Input file not found: $full_input_path"
#    exit 1
#fi

echo "Processing file: $input_file"
echo "Full path: $full_input_path"
echo "Output base: $output_base"

# Run the ROOT macro
root.exe -l -b << EOF
.L Fun4All_TrackSeeding.C
Fun4All_TrackSeeding(
  100,                         // nEvents=0 means process all events
  "${input_file}",           // clusterfilename
  "${input_dir}",            // dir
  "${output_base}",          // outfilename
  false,                     // convertSeeds
  true                       // doKFParticle
);
EOF

echo "Processing complete for $input_file"