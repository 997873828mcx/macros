#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 6 ]]; then
  echo "Usage: $0 <campaign> <hepmc_file_or_list> <total_events> [events_per_job] [events_per_input_file] [use_beam_vertex]"
  echo "Example single file: $0 pp_minbias_truth_100 output/pp_minbias_test.hepmc 100 10"
  echo "Example file list:   $0 pp_minbias_truth_100k output/pp_minbias_hepmc_100k/hepmc_files.list 100000 10 1000 false"
  exit 1
fi

campaign="$1"
hepmc_file="$2"
total_events="$3"
events_per_job="${4:-10}"
events_per_input_file="${5:-1000}"
use_beam_vertex="${6:-false}"
g4_seed_base="${G4_SEED_BASE:-0}"
g4_uniform_bfield_t="${G4_UNIFORM_BFIELD_T:-none}"
write_truth_flat_tree="${WRITE_TRUTH_FLAT_TREE:-true}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ "${hepmc_file}" != /* ]]; then
  hepmc_file="${SCRIPT_DIR}/${hepmc_file}"
fi

if [[ ! -f "${hepmc_file}" ]]; then
  echo "Error: HepMC input file not found: ${hepmc_file}" >&2
  exit 2
fi

if [[ "${total_events}" -le 0 ]]; then
  echo "Error: total_events must be positive, got ${total_events}" >&2
  exit 2
fi

if [[ "${events_per_job}" -le 0 ]]; then
  echo "Error: events_per_job must be positive, got ${events_per_job}" >&2
  exit 2
fi

if [[ "${events_per_input_file}" -le 0 ]]; then
  echo "Error: events_per_input_file must be positive, got ${events_per_input_file}" >&2
  exit 2
fi

if [[ "${use_beam_vertex}" != "true" && "${use_beam_vertex}" != "false" ]]; then
  echo "Error: use_beam_vertex must be true or false, got ${use_beam_vertex}" >&2
  exit 2
fi

if ! [[ "${g4_seed_base}" =~ ^[0-9]+$ ]]; then
  echo "Error: G4_SEED_BASE must be a nonnegative integer, got ${g4_seed_base}" >&2
  exit 2
fi

if [[ "${g4_uniform_bfield_t}" != "none" && "${g4_uniform_bfield_t}" != "NONE" ]] &&
   ! [[ "${g4_uniform_bfield_t}" =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]; then
  echo "Error: G4_UNIFORM_BFIELD_T must be numeric or 'none', got ${g4_uniform_bfield_t}" >&2
  exit 2
fi

if [[ "${write_truth_flat_tree}" != "true" && "${write_truth_flat_tree}" != "false" ]]; then
  echo "Error: WRITE_TRUTH_FLAT_TREE must be true or false, got ${write_truth_flat_tree}" >&2
  exit 2
fi

n_jobs=$(((total_events + events_per_job - 1) / events_per_job))
is_filelist=false
if [[ "${hepmc_file}" == *.list ]]; then
  is_filelist=true
  n_files=$(grep -cv '^[[:space:]]*$' "${hepmc_file}")
  expected_files=$(((total_events + events_per_input_file - 1) / events_per_input_file))
  if [[ "${n_files}" -lt "${expected_files}" ]]; then
    echo "Error: file list has ${n_files} files but ${expected_files} are needed" >&2
    exit 2
  fi
  if ((events_per_input_file % events_per_job != 0)); then
    echo "Error: events_per_job=${events_per_job} must divide events_per_input_file=${events_per_input_file}" >&2
    exit 2
  fi
fi

mkdir -p "log/${campaign}" "log/tmp/${campaign}"
mkdir -p "output/${campaign}" "output/${campaign}/completed"
# safety fallback in case submit variables are not applied for any reason
mkdir -p "log/default" "log/tmp/default"

echo "Submitting campaign:"
echo "  campaign=${campaign}"
if [[ "${is_filelist}" == true ]]; then
  echo "  input_filelist=${hepmc_file}"
  echo "  events/input file=${events_per_input_file}"
else
  echo "  input_hepmc=${hepmc_file}"
fi
echo "  total_events=${total_events}"
echo "  events/job=${events_per_job}"
echo "  n_jobs=${n_jobs}"
echo "  use_beam_vertex=${use_beam_vertex}"
echo "  g4_seed_base=${g4_seed_base} (0 means default random seeding)"
echo "  g4_uniform_bfield_t=${g4_uniform_bfield_t}"
echo "  write_truth_flat_tree=${write_truth_flat_tree}"

condor_args=(
  -append "campaign = ${campaign}" \
  -append "n_jobs = ${n_jobs}" \
  -append "total_events = ${total_events}" \
  -append "events_per_job = ${events_per_job}" \
  -append "events_per_input_file = ${events_per_input_file}" \
  -append "use_beam_vertex = ${use_beam_vertex}" \
  -append "g4_seed_base = ${g4_seed_base}" \
  -append "g4_uniform_bfield_t = ${g4_uniform_bfield_t}" \
  -append "write_truth_flat_tree = ${write_truth_flat_tree}" \
  -append "output_base_dir = output" \
)

if [[ "${is_filelist}" == true ]]; then
  condor_args+=(-append "input_filelist = ${hepmc_file}")
else
  condor_args+=(-append "input_hepmc = ${hepmc_file}")
fi

condor_submit "${condor_args[@]}" general_sim.job
