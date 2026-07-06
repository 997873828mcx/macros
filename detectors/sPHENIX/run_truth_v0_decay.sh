#!/bin/bash
set -eo pipefail

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/
set -u

process_id=${1?Error: process ID is not given}
files_per_job=${2:-100}
events_per_input_file=${3:-10}
total_files=${4:-0}

campaign_tag=${CAMPAIGN_TAG:-truth_v0_decays}
output_base_dir=${OUTPUT_BASE_DIR:-/sphenix/user/dcxchenxi/develope/macros/detectors/sPHENIX/output}
input_dst_filelist=${INPUT_DST_FILELIST:?Error: INPUT_DST_FILELIST is not set}
completed_dir=${COMPLETED_OUTPUT_DIR:-${output_base_dir}/${campaign_tag}/completed}
input_path_base=${INPUT_DST_PATH_BASE:-$(dirname "${output_base_dir}")}

campaign_tag="${campaign_tag%%;*}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ "${files_per_job}" -le 0 ]]; then
  echo "Error: files_per_job must be positive, got ${files_per_job}" >&2
  exit 2
fi

if [[ "${events_per_input_file}" -le 0 ]]; then
  echo "Error: events_per_input_file must be positive, got ${events_per_input_file}" >&2
  exit 2
fi

if [[ ! -f "${input_dst_filelist}" ]]; then
  echo "Error: DST file list not found: ${input_dst_filelist}" >&2
  exit 2
fi

available_files=$(grep -cv '^[[:space:]]*$' "${input_dst_filelist}")
if [[ "${total_files}" -le 0 || "${total_files}" -gt "${available_files}" ]]; then
  total_files=${available_files}
fi

file_start=$((process_id * files_per_job))
remaining=$((total_files - file_start))
if [[ "${remaining}" -le 0 ]]; then
  echo "Nothing to do: process=${process_id} file_start=${file_start} total_files=${total_files}"
  exit 0
fi

nfiles=${files_per_job}
if [[ "${remaining}" -lt "${nfiles}" ]]; then
  nfiles=${remaining}
fi

nevents=$((nfiles * events_per_input_file))
padded_id=$(printf "%05d" "${process_id}")
job_outdir="${output_base_dir}/${campaign_tag}"
mkdir -p "${job_outdir}" "${completed_dir}"

chunk_list="${job_outdir}/dst_chunk_${padded_id}_fskip${file_start}_nfiles${nfiles}.list"
grep -v '^[[:space:]]*$' "${input_dst_filelist}" | \
  sed -n "$((file_start + 1)),$((file_start + nfiles))p" | \
  awk -v base="${input_path_base}" '{ if ($0 ~ /^\//) { print $0 } else { print base "/" $0 } }' > "${chunk_list}"

if [[ "$(wc -l < "${chunk_list}")" -ne "${nfiles}" ]]; then
  echo "Error: chunk list has $(wc -l < "${chunk_list}") files, expected ${nfiles}" >&2
  exit 2
fi

missing=0
while IFS= read -r dst_file; do
  if [[ ! -f "${dst_file}" ]]; then
    echo "Error: missing DST input file: ${dst_file}" >&2
    missing=$((missing + 1))
  fi
done < "${chunk_list}"
if [[ "${missing}" -ne 0 ]]; then
  echo "Error: ${missing} missing DST input files in ${chunk_list}" >&2
  exit 2
fi

outroot="${job_outdir}/truth_v0_decays_${padded_id}_fskip${file_start}_nfiles${nfiles}.root"

echo "Running generator-level truth V0 decay tree:"
echo "  campaign=${campaign_tag}"
echo "  job=${process_id}"
echo "  input_dst_filelist=${input_dst_filelist}"
echo "  input_path_base=${input_path_base}"
echo "  chunk_list=${chunk_list}"
echo "  file_start=${file_start}"
echo "  nfiles=${nfiles}"
echo "  events/input file=${events_per_input_file}"
echo "  nevents=${nevents}"
echo "  out=${outroot}"

root.exe -l -b -q "Fun4All_TruthV0Decay.C(${nevents}, \"${chunk_list}\", \"${outroot}\")"

if [[ -f "${outroot}" ]]; then
  base="$(basename "${outroot}")"
  dst="${completed_dir}/${base}"
  if [[ -e "${dst}" ]]; then
    stem="${base%.*}"
    ext="${base##*.}"
    dst="${completed_dir}/${stem}__dup_$(date +%s)_${process_id}.${ext}"
    echo "Warning: output already exists, moving as ${dst}"
  fi
  mv -f "${outroot}" "${dst}"
  echo "Moved completed file to ${dst}"
else
  echo "Error: expected output was not produced: ${outroot}" >&2
  exit 3
fi

echo "Done"
