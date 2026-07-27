#!/bin/bash
set -eo pipefail

export USER="$(id -u -n)"
export LOGNAME="${USER}"
export HOME="/sphenix/u/${LOGNAME}"
source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/
set -u

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_dir}"

output_dir="${1:-output/pion_phi_hitsets}"
mkdir -p "${output_dir}"

plus_hepmc="${output_dir}/piplus_24_phi_p0p2to1.hepmc"
minus_hepmc="${output_dir}/piminus_24_phi_p0p2to1.hepmc"
plus_dst="${output_dir}/DST_TRKR_HITSET_piplus_24_phi_p0p2to1.root"
minus_dst="${output_dir}/DST_TRKR_HITSET_piminus_24_phi_p0p2to1.root"

verify_dst()
{
  local dst_file=$1
  if [[ ! -s "${dst_file}" ]]; then
    echo "Error: expected nonempty DST was not produced: ${dst_file}" >&2
    exit 3
  fi
  if ! rootls -t "${dst_file}" | rg -q 'DST#TRKR#TRKR_HITSET'; then
    echo "Error: TRKR_HITSET is missing from ${dst_file}" >&2
    exit 4
  fi
}

root.exe -l -b -q \
  "GeneratePionPhiScanHepMC.C(\"${plus_hepmc}\", 211, 24, 0.2, 1.0, 0.0, 1)"
root.exe -l -b -q \
  "GeneratePionPhiScanHepMC.C(\"${minus_hepmc}\", -211, 24, 0.2, 1.0, 0.0, 1)"

# Ninth macro argument requests TPC digitization and a filtered DST containing
# TRKR_HITSET and TRKR_HITTRUTHASSOC. Empty outputFile disables the truth tree.
root.exe -l -b -q \
  "Fun4All_G4_sPHENIX.C(1, \"${plus_hepmc}\", \"\", \"\", 0, \".\", false, \"\", \"${plus_dst}\")"
verify_dst "${plus_dst}"
root.exe -l -b -q \
  "Fun4All_G4_sPHENIX.C(1, \"${minus_hepmc}\", \"\", \"\", 0, \".\", false, \"\", \"${minus_dst}\")"
verify_dst "${minus_dst}"

echo "Created:"
echo "  ${plus_dst}"
echo "  ${minus_dst}"
