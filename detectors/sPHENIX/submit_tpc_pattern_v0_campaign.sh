#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 5 ]]; then
  echo "Usage: $0 <campaign> <dst_filelist> [files_per_job] [events_per_input_file] [total_files]"
  echo "Example: $0 tpc_pattern_v0_79513_good /path/to/good_dst_files.list 1 500"
  exit 1
fi

campaign="$1"
dst_filelist="$2"
files_per_job="${3:-1}"
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
pre_track_pt_min="${V0_PRE_TRACK_PT_MIN:-0.05}"
pre_track_dca_xy_min="${V0_PRE_TRACK_DCA_XY_MIN:--1.0}"
pre_pair_dca_max="${V0_PRE_PAIR_DCA_MAX:-10.0}"
pre_lproj_min="${V0_PRE_LPROJ_MIN:--1.0}"
pre_cos_theta_min="${V0_PRE_COS_THETA_MIN:--2.0}"
use_final_track_helix="${V0_USE_FINAL_TRACK_HELIX:-false}"
point_order="${V0_POINT_ORDER:-auto}"

mkdir -p "log/${campaign}" "log/tmp/${campaign}"
mkdir -p "${output_base_dir}/${campaign}" "${output_base_dir}/${campaign}/completed"
chmod +x "${SCRIPT_DIR}/run_tpc_pattern_v0.sh"

echo "Submitting TPC pattern-reco V0 campaign:"
echo "  campaign=${campaign}"
echo "  dst_filelist=${dst_filelist}"
echo "  total_files=${total_files}"
echo "  files/job=${files_per_job}"
echo "  events/input file=${events_per_input_file}"
echo "  n_jobs=${n_jobs}"
echo "  output=${output_base_dir}/${campaign}/completed"
echo "  preselection: pt>${pre_track_pt_min}, dca_xy_min=${pre_track_dca_xy_min}, pairDCA<${pre_pair_dca_max}, Lproj>${pre_lproj_min}, cosTheta>${pre_cos_theta_min}"
echo "  use_final_track_helix=${use_final_track_helix}"
echo "  point_order=${point_order}"

condor_submit \
  -append "campaign = ${campaign}" \
  -append "input_dst_filelist = ${dst_filelist}" \
  -append "files_per_job = ${files_per_job}" \
  -append "events_per_input_file = ${events_per_input_file}" \
  -append "total_files = ${total_files}" \
  -append "n_jobs = ${n_jobs}" \
  -append "output_base_dir = ${output_base_dir}" \
  -append "input_path_base = ${input_path_base}" \
  -append "pre_track_pt_min = ${pre_track_pt_min}" \
  -append "pre_track_dca_xy_min = ${pre_track_dca_xy_min}" \
  -append "pre_pair_dca_max = ${pre_pair_dca_max}" \
  -append "pre_lproj_min = ${pre_lproj_min}" \
  -append "pre_cos_theta_min = ${pre_cos_theta_min}" \
  -append "use_final_track_helix = ${use_final_track_helix}" \
  -append "point_order = ${point_order}" \
  tpc_pattern_v0.job
