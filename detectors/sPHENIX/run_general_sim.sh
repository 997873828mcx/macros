#!/bin/bash
set -eo pipefail

# Setup environment
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/
set -u

# Args: process index, events per job, and optional total events.
process_id=${1?Error: process ID is not given}
events_per_job=${2:-10}
total_events=${3:-0}
campaign_tag=${CAMPAIGN_TAG:-default}
output_base_dir=${OUTPUT_BASE_DIR:-output}
hepmc_input=${HEPMC_INPUT:-/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output/pp_minbias_test.hepmc}
hepmc_filelist=${HEPMC_FILELIST:-}
hepmc_events_per_file=${HEPMC_EVENTS_PER_FILE:-1000}
use_beam_vertex=${USE_BEAM_VERTEX:-false}
export SIM_PROCESS_ID=${process_id}

# Guard against malformed CAMPAIGN_TAG values where extra env assignments are
# accidentally appended (e.g. "tag;OUTPUT_BASE_DIR=...").
campaign_tag="${campaign_tag%%;*}"

# Working directory containing the macro (where this script lives)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

padded_id=$(printf "%05d" "$process_id")
job_outdir="${output_base_dir}/${campaign_tag}"
mkdir -p "${job_outdir}"

if [[ "${events_per_job}" -le 0 ]]; then
  echo "Error: events_per_job must be positive, got ${events_per_job}" >&2
  exit 2
fi

global_skip=$((process_id * events_per_job))
skip=${global_skip}
nevents=${events_per_job}
if [[ "${total_events}" -gt 0 ]]; then
  remaining=$((total_events - global_skip))
  if [[ "${remaining}" -le 0 ]]; then
    echo "Nothing to do: process=${process_id} global_skip=${global_skip} total_events=${total_events}"
    exit 0
  fi
  if [[ "${remaining}" -lt "${nevents}" ]]; then
    nevents=${remaining}
  fi
fi

if [[ -n "${hepmc_filelist}" ]]; then
  if [[ ! -f "${hepmc_filelist}" ]]; then
    echo "Error: HepMC file list not found: ${hepmc_filelist}" >&2
    exit 2
  fi
  if [[ "${hepmc_events_per_file}" -le 0 ]]; then
    echo "Error: HEPMC_EVENTS_PER_FILE must be positive, got ${hepmc_events_per_file}" >&2
    exit 2
  fi
  file_index=$((global_skip / hepmc_events_per_file))
  skip=$((global_skip % hepmc_events_per_file))
  hepmc_input="$(sed -n "$((file_index + 1))p" "${hepmc_filelist}")"
  if [[ -z "${hepmc_input}" ]]; then
    echo "Error: no HepMC file for file_index=${file_index} in ${hepmc_filelist}" >&2
    exit 2
  fi
  if ((skip + nevents > hepmc_events_per_file)); then
    echo "Error: job would cross a HepMC file boundary: skip=${skip} nevents=${nevents} events/file=${hepmc_events_per_file}" >&2
    echo "Choose an events_per_job value that divides HEPMC_EVENTS_PER_FILE." >&2
    exit 2
  fi
else
  file_index=0
fi

if [[ ! -f "${hepmc_input}" ]]; then
  echo "Error: HepMC input file not found: ${hepmc_input}" >&2
  exit 2
fi

if [[ "${use_beam_vertex}" != "true" && "${use_beam_vertex}" != "false" ]]; then
  echo "Error: USE_BEAM_VERTEX must be true or false, got ${use_beam_vertex}" >&2
  exit 2
fi

outroot="${job_outdir}/tpc_truthpoints_${padded_id}_gskip${global_skip}_nev${nevents}.root"
outdst="${job_outdir}/tpc_truthpoints_dst_${padded_id}_gskip${global_skip}_nev${nevents}.root"

# If set, completed outputs are moved here at the end of a successful job.
# Set to empty string to disable moving.
completed_dir="${COMPLETED_OUTPUT_DIR:-${job_outdir}/completed}"

echo "Running TPC truth-point simulation:"
echo "  campaign=${campaign_tag}"
echo "  job=${process_id}"
echo "  file_index=${file_index}"
echo "  input_hepmc=${hepmc_input}"
echo "  skip=${skip}"
echo "  global_skip=${global_skip}"
echo "  nevents=${nevents}"
echo "  total_events=${total_events}"
echo "  use_beam_vertex=${use_beam_vertex}"
echo "  flat_tree_out=${outroot}"
echo "  truth_dst_out=${outdst}"

# Run the stripped-down HepMC -> Geant4 -> TPC truth-point macro.
# The last argument writes a DST containing EventHeader, G4TruthInfo, and
# G4HIT_TPC_TRUECLUSTER for second-stage Fun4All V0 processing.
root.exe -l -b -q "Fun4All_G4_sPHENIX.C(${nevents}, \"${hepmc_input}\", \"${outroot}\", \"\", ${skip}, \".\", ${use_beam_vertex}, \"${outdst}\")"

if [[ -n "${completed_dir}" ]]; then
  mkdir -p "${completed_dir}"
  moved=0
  for f in "${outroot}" "${outdst}"; do
    if [[ -f "${f}" ]]; then
      base="$(basename "${f}")"
      dst="${completed_dir}/${base}"
      if [[ -e "${dst}" ]]; then
        if [[ "${base}" == *.* ]]; then
          stem="${base%.*}"
          ext="${base##*.}"
          alt="${completed_dir}/${stem}__dup_$(date +%s)_${process_id}.${ext}"
        else
          alt="${completed_dir}/${base}__dup_$(date +%s)_${process_id}"
        fi
        echo "Warning: ${dst} already exists, moving as ${alt}"
        mv -f "${f}" "${alt}"
      else
        mv -f "${f}" "${dst}"
      fi
      moved=$((moved + 1))
    fi
  done
  echo "Moved ${moved} completed files to ${completed_dir}"
fi

echo "Done: flat_tree_out=${outroot} truth_dst_out=${outdst}"
