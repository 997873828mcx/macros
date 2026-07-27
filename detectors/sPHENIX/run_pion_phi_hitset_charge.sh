#!/bin/bash
set -eo pipefail

export USER="$(id -u -n)"
export LOGNAME="${USER}"
export HOME="/sphenix/u/${LOGNAME}"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/
set -u

charge=${1?Usage: run_pion_phi_hitset_charge.sh piplus|piminus}

case "${charge}" in
  piplus)
    input_name="piplus_24_phi_p0p2to1.hepmc"
    output_name="DST_TRKR_HITSET_piplus_24_phi_p0p2to1.root"
    ;;
  piminus)
    input_name="piminus_24_phi_p0p2to1.hepmc"
    output_name="DST_TRKR_HITSET_piminus_24_phi_p0p2to1.root"
    ;;
  *)
    echo "Error: charge must be piplus or piminus, got ${charge}" >&2
    exit 2
    ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

output_dir="output/pion_phi_hitsets"
input_file="${output_dir}/${input_name}"
output_file="${output_dir}/${output_name}"

if [[ ! -s "${input_file}" ]]; then
  echo "Error: missing HepMC input: ${input_file}" >&2
  exit 3
fi

# Never silently replace a completed charge sample.
if [[ -e "${output_file}" ]]; then
  echo "Error: refusing to overwrite existing output: ${output_file}" >&2
  exit 4
fi

echo "Starting ${charge} TPC TRKR_HITSET simulation"
echo "  host=$(hostname)"
echo "  input=${input_file}"
echo "  output=${output_file}"
echo "  start=$(date --iso-8601=seconds)"

root.exe -l -b -q \
  "Fun4All_G4_sPHENIX.C(1, \"${input_file}\", \"\", \"\", 0, \".\", false, \"\", \"${output_file}\")"

# Fun4All can abort a run while ROOT still exits with status zero, so validate
# the actual file and required branch before reporting success.
if [[ ! -s "${output_file}" ]]; then
  echo "Error: nonempty output DST was not produced: ${output_file}" >&2
  exit 5
fi

if ! rootls -t "${output_file}" | /usr/bin/grep -q 'DST#TRKR#TRKR_HITSET'; then
  echo "Error: TRKR_HITSET is missing from ${output_file}" >&2
  exit 6
fi

echo "Completed ${charge} successfully at $(date --iso-8601=seconds)"
ls -lh "${output_file}"

