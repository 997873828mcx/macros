#!/bin/bash
set -eo pipefail

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/
set -u

process_id=${1?Error: process ID is not given}
events_per_job=${2:-1000}
total_events=${3:-0}
seed_base=${4:-100000}

campaign_tag=${CAMPAIGN_TAG:-pp_minbias_hepmc}
output_base_dir=${OUTPUT_BASE_DIR:-output}
pythia_config=${PYTHIA_CONFIG:-}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ "${events_per_job}" -le 0 ]]; then
  echo "Error: events_per_job must be positive, got ${events_per_job}" >&2
  exit 2
fi

skip=$((process_id * events_per_job))
nevents=${events_per_job}
if [[ "${total_events}" -gt 0 ]]; then
  remaining=$((total_events - skip))
  if [[ "${remaining}" -le 0 ]]; then
    echo "Nothing to do: process=${process_id} skip=${skip} total_events=${total_events}"
    exit 0
  fi
  if [[ "${remaining}" -lt "${nevents}" ]]; then
    nevents=${remaining}
  fi
fi

padded_id="$(printf "%05d" "${process_id}")"
job_outdir="${output_base_dir}/${campaign_tag}"
completed_dir="${COMPLETED_OUTPUT_DIR:-${job_outdir}/completed}"
mkdir -p "${job_outdir}" "${completed_dir}"

seed=$((seed_base + process_id))
outfile="${job_outdir}/pp_minbias_${padded_id}_skip${skip}_nev${nevents}.hepmc"

echo "Running Pythia8 HepMC generation:"
echo "  campaign=${campaign_tag}"
echo "  job=${process_id}"
echo "  nevents=${nevents}"
echo "  total_events=${total_events}"
echo "  seed=${seed}"
echo "  config=${pythia_config:-default CALIBRATIONROOT Detroit minbias config}"
echo "  output=${outfile}"

root.exe -l -b -q "GeneratePythia8HepMC.C(${nevents}, \"${outfile}\", \"${pythia_config}\", ${seed}, 0)"

if [[ -f "${outfile}" ]]; then
  base="$(basename "${outfile}")"
  dst="${completed_dir}/${base}"
  if [[ -e "${dst}" ]]; then
    alt="${completed_dir}/${base%.hepmc}__dup_$(date +%s)_${process_id}.hepmc"
    echo "Warning: ${dst} already exists, moving as ${alt}"
    mv -f "${outfile}" "${alt}"
  else
    mv -f "${outfile}" "${dst}"
  fi
  echo "Moved completed HepMC file to ${completed_dir}"
else
  echo "Error: expected output file was not produced: ${outfile}" >&2
  exit 3
fi

echo "Done: ${completed_dir}/${base}"
