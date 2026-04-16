#!/bin/bash
set -e

# Setup environment
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh new

# Args: process index and number of events per job (default 10)
process_id=${1?Error: process ID is not given}
nevents=${2:-10}
campaign_tag=${CAMPAIGN_TAG:-default}
output_base_dir=${OUTPUT_BASE_DIR:-output}
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
completed_dir="${COMPLETED_OUTPUT_DIR:-${job_outdir}/completed}"
mkdir -p "${completed_dir}"

outroot="${job_outdir}/general_${padded_id}.root"
qa_outroot="${outroot%.root}_qa.root"
reco_pt_outroot="${outroot%.root}_reco_pt.root"
svtx_eval_outroot="${outroot%.root}_g4svtx_eval.root"
jet_eval_outroot="${outroot%.root}_g4jet_eval.root"

echo "Running general sim: campaign=${campaign_tag} job=$process_id nevents=$nevents"
echo "  output base=${outroot}"
echo "  expected qa=${qa_outroot}"
echo "  expected reco pt tree=${reco_pt_outroot}"
echo "  expected svtx_eval=${svtx_eval_outroot}"
echo "  completed dir=${completed_dir}"

# Run the default sPHENIX macro. The current branch controls the generator setup,
# and the outputFile path is used as the base for QA/eval file names.
root.exe -l -b -q "Fun4All_G4_sPHENIX.C(${nevents}, \"\", \"${outroot}\", \"\", 0, \".\")"

moved=0
shopt -s nullglob
job_outputs=("${job_outdir}/general_${padded_id}"*)
shopt -u nullglob

for f in "${job_outputs[@]}"; do
  [[ -f "${f}" ]] || continue
  base="$(basename "${f}")"
  dst="${completed_dir}/${base}"
  if [[ -e "${dst}" ]]; then
    if [[ "${base}" == *.* ]]; then
      stem="${base%.*}"
      ext="${base##*.}"
      dst="${completed_dir}/${stem}__dup_$(date +%s)_${process_id}.${ext}"
    else
      dst="${completed_dir}/${base}__dup_$(date +%s)_${process_id}"
    fi
    echo "Warning: destination exists, moving ${base} as $(basename "${dst}")"
  fi
  mv -f "${f}" "${dst}"
  moved=$((moved + 1))
done

echo "Done: moved ${moved} files to ${completed_dir}"
