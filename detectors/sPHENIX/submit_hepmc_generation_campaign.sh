#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 5 ]]; then
  echo "Usage: $0 <campaign> [total_events] [events_per_job] [seed_base] [pythia_config]"
  echo "Example: $0 pp_minbias_hepmc_100k 100000 1000 100000"
  echo "Example with explicit config: $0 pp_minbias_hepmc_100k 100000 1000 100000 /path/to/phpythia8_detroit_minBias.cfg"
  exit 1
fi

campaign="$1"
total_events="${2:-100000}"
events_per_job="${3:-1000}"
seed_base="${4:-100000}"
pythia_config="${5:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ "${total_events}" -le 0 ]]; then
  echo "Error: total_events must be positive, got ${total_events}" >&2
  exit 2
fi

if [[ "${events_per_job}" -le 0 ]]; then
  echo "Error: events_per_job must be positive, got ${events_per_job}" >&2
  exit 2
fi

if [[ "${seed_base}" -le 0 ]]; then
  echo "Error: seed_base must be positive, got ${seed_base}" >&2
  exit 2
fi

if [[ -n "${pythia_config}" && "${pythia_config}" != /* ]]; then
  pythia_config="${SCRIPT_DIR}/${pythia_config}"
fi

if [[ -n "${pythia_config}" && ! -f "${pythia_config}" ]]; then
  echo "Error: Pythia config file not found: ${pythia_config}" >&2
  exit 2
fi

n_jobs=$(((total_events + events_per_job - 1) / events_per_job))

mkdir -p "log/${campaign}" "log/tmp/${campaign}"
mkdir -p "output/${campaign}" "output/${campaign}/completed"

filelist="output/${campaign}/hepmc_files.list"
: > "${filelist}"
for ((i = 0; i < n_jobs; ++i)); do
  skip=$((i * events_per_job))
  nevents="${events_per_job}"
  remaining=$((total_events - skip))
  if [[ "${remaining}" -lt "${nevents}" ]]; then
    nevents="${remaining}"
  fi
  padded_id="$(printf "%05d" "${i}")"
  echo "${SCRIPT_DIR}/output/${campaign}/completed/pp_minbias_${padded_id}_skip${skip}_nev${nevents}.hepmc" >> "${filelist}"
done

echo "Submitting HepMC generation campaign:"
echo "  campaign=${campaign}"
echo "  total_events=${total_events}"
echo "  events/job=${events_per_job}"
echo "  n_jobs=${n_jobs}"
echo "  seed_base=${seed_base}"
echo "  pythia_config=${pythia_config:-default CALIBRATIONROOT Detroit minbias config}"
echo "  expected_filelist=${SCRIPT_DIR}/${filelist}"

condor_args=(
  -append "campaign = ${campaign}"
  -append "n_jobs = ${n_jobs}"
  -append "total_events = ${total_events}"
  -append "events_per_job = ${events_per_job}"
  -append "seed_base = ${seed_base}"
  -append "output_base_dir = output"
)

if [[ -n "${pythia_config}" ]]; then
  condor_args+=(-append "pythia_config = ${pythia_config}")
fi

condor_submit "${condor_args[@]}" hepmc_generation.job
