#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 5 ]]; then
  echo "Usage: $0 <campaign> <dst_filelist> [files_per_job] [events_per_input_file] [total_files]"
  echo "Example: $0 raw_finaltrack_pt_79513 output_dst_79513_files.list 10 500"
  exit 1
fi

campaign="$1"
dst_filelist="$2"
files_per_job="${3:-10}"
events_per_input_file="${4:-500}"
total_files="${5:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ "${dst_filelist}" != /* ]]; then
  dst_filelist="${SCRIPT_DIR}/${dst_filelist}"
fi

if [[ ! -f "${dst_filelist}" ]]; then
  echo "Error: DST file list not found: ${dst_filelist}" >&2
  exit 2
fi

if [[ "${files_per_job}" -le 0 ]]; then
  echo "Error: files_per_job must be positive, got ${files_per_job}" >&2
  exit 2
fi

if [[ "${events_per_input_file}" -le 0 ]]; then
  echo "Error: events_per_input_file must be positive, got ${events_per_input_file}" >&2
  exit 2
fi

available_files=$(grep -cv '^[[:space:]]*$' "${dst_filelist}")
if [[ "${total_files}" -le 0 || "${total_files}" -gt "${available_files}" ]]; then
  total_files=${available_files}
fi

if [[ "${total_files}" -le 0 ]]; then
  echo "Error: no files to process in ${dst_filelist}" >&2
  exit 2
fi

n_jobs=$(((total_files + files_per_job - 1) / files_per_job))
output_base_dir="/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output"
input_path_base="$(dirname "${output_base_dir}")"
node_name="${RAW_FINALTRACK_NODE:-FINALTRACKS}"
pt_min="${RAW_FINALTRACK_PT_MIN:-0}"
pt_max="${RAW_FINALTRACK_PT_MAX:-10}"
pt_bins="${RAW_FINALTRACK_PT_BINS:-250}"

mkdir -p "log/${campaign}" "log/tmp/${campaign}"
mkdir -p "${output_base_dir}/${campaign}" "${output_base_dir}/${campaign}/completed"
chmod +x "${SCRIPT_DIR}/run_raw_finaltrack_pt.sh"

echo "Submitting raw FINALTRACKS pT campaign:"
echo "  campaign=${campaign}"
echo "  dst_filelist=${dst_filelist}"
echo "  total_files=${total_files}"
echo "  files/job=${files_per_job}"
echo "  events/input file=${events_per_input_file}"
echo "  n_jobs=${n_jobs}"
echo "  output=${output_base_dir}/${campaign}/completed"
echo "  node=${node_name}"
echo "  pt range=${pt_min} to ${pt_max} with ${pt_bins} bins"

condor_submit \
  -append "campaign = ${campaign}" \
  -append "input_dst_filelist = ${dst_filelist}" \
  -append "files_per_job = ${files_per_job}" \
  -append "events_per_input_file = ${events_per_input_file}" \
  -append "total_files = ${total_files}" \
  -append "n_jobs = ${n_jobs}" \
  -append "output_base_dir = ${output_base_dir}" \
  -append "input_path_base = ${input_path_base}" \
  -append "node_name = ${node_name}" \
  -append "pt_min = ${pt_min}" \
  -append "pt_max = ${pt_max}" \
  -append "pt_bins = ${pt_bins}" \
  raw_finaltrack_pt.job
