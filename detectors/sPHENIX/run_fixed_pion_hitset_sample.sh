#!/bin/bash
set -eo pipefail

export USER="$(id -u -n)"
export LOGNAME="${USER}"
export HOME="/sphenix/u/${LOGNAME}"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/
set -u

sample=${1?Usage: run_fixed_pion_hitset_sample.sh SAMPLE}

case "${sample}" in
  eta_p0p5_phi_0)
    eta="0.5"
    phi="0.0"
    ;;
  eta_p0p5_phi_pi4)
    eta="0.5"
    phi="0.7853981633974483"
    ;;
  eta_m0p5_phi_0)
    eta="-0.5"
    phi="0.0"
    ;;
  eta_m0p5_phi_pi4)
    eta="-0.5"
    phi="0.7853981633974483"
    ;;
  *)
    echo "Error: unknown sample ${sample}" >&2
    exit 2
    ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

output_dir="output/fixed_piplus_pt0p3_z0"
output_file="${output_dir}/DST_TRKR_HITSET_piplus_pt0p3_${sample}_z0_10events.root"
mkdir -p "${output_dir}"

if [[ -e "${output_file}" ]]; then
  echo "Error: refusing to overwrite existing output: ${output_file}" >&2
  exit 3
fi

echo "Starting fixed pi+ sample ${sample}"
echo "  host=$(hostname)"
echo "  events=10"
echo "  pt=0.3 GeV eta=${eta} phi=${phi} vertex=(0,0,0)"
echo "  output=${output_file}"
echo "  start=$(date --iso-8601=seconds)"

# Empty input and truth-tree arguments select the fixed Simple generator and
# write only the filtered tracker-hit DST.
root.exe -l -b -q \
  "Fun4All_G4_sPHENIX.C(10, \"\", \"\", \"\", 0, \".\", false, \"\", \"${output_file}\", true, 211, 0.3, ${eta}, ${phi}, 0.0)"

if [[ ! -s "${output_file}" ]]; then
  echo "Error: nonempty output DST was not produced: ${output_file}" >&2
  exit 4
fi

if ! rootls -t "${output_file}" | /usr/bin/grep -q 'DST#TRKR#TRKR_HITSET'; then
  echo "Error: TRKR_HITSET is missing from ${output_file}" >&2
  exit 5
fi

echo "Completed ${sample} successfully at $(date --iso-8601=seconds)"
ls -lh "${output_file}"

