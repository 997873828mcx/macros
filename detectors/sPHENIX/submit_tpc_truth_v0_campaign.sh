#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 5 ]]; then
  echo "Usage: $0 <campaign> <dst_filelist> [files_per_job] [events_per_input_file] [total_files]"
  echo "Example: $0 tpc_v0_100k_zero_vtx /path/to/dst_files.list 100 10 10000"
  exit 1
fi

campaign="$1"
dst_filelist="$2"
files_per_job="${3:-100}"
events_per_input_file="${4:-10}"
total_files="${5:-0}"
fit_method="${V0_FIT_METHOD:-helix}"
point_order="${V0_POINT_ORDER:-path}"
kalman_sigma_rphi_cm="${V0_KALMAN_SIGMA_RPHI_CM:-0.03}"
kalman_sigma_r_cm="${V0_KALMAN_SIGMA_R_CM:-0.03}"
kalman_sigma_z_cm="${V0_KALMAN_SIGMA_Z_CM:-0.05}"
write_kalman_innovation_diagnostics="${V0_WRITE_KALMAN_INNOVATION_DIAGNOSTICS:-false}"
kalman_uniform_propagator="${V0_KALMAN_UNIFORM_PROPAGATOR:-analytic}"
bfield_t="${V0_BFIELD_T:-1.4}"
coarse_steps="${V0_COARSE_STEPS:-64}"
pca_candidates="${V0_PCA_CANDIDATES:-32}"

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

case "${fit_method}" in
  helix|HELIX|kalman|KALMAN|none|NONE)
    ;;
  *)
    echo "Error: V0_FIT_METHOD must be helix, kalman, or none, got ${fit_method}" >&2
    exit 2
    ;;
esac

case "${kalman_uniform_propagator}" in
  analytic|ANALYTIC|rk|RK)
    ;;
  *)
    echo "Error: V0_KALMAN_UNIFORM_PROPAGATOR must be analytic or rk, got ${kalman_uniform_propagator}" >&2
    exit 2
    ;;
esac

for boolean_value in "${write_kalman_innovation_diagnostics}"; do
  if [[ "${boolean_value}" != "true" && "${boolean_value}" != "false" ]]; then
    echo "Error: V0_WRITE_KALMAN_INNOVATION_DIAGNOSTICS must be true or false" >&2
    exit 2
  fi
done

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

mkdir -p "log/${campaign}" "log/tmp/${campaign}"
mkdir -p "${output_base_dir}/${campaign}" "${output_base_dir}/${campaign}/completed"
chmod +x "${SCRIPT_DIR}/run_tpc_truth_v0.sh"

echo "Submitting TPC truth V0 campaign:"
echo "  campaign=${campaign}"
echo "  dst_filelist=${dst_filelist}"
echo "  total_files=${total_files}"
echo "  files/job=${files_per_job}"
echo "  events/input file=${events_per_input_file}"
echo "  n_jobs=${n_jobs}"
echo "  output=${output_base_dir}/${campaign}/completed"
echo "  fit_method=${fit_method}"
echo "  point_order=${point_order}"
echo "  bfield_t=${bfield_t}"
echo "  kalman_uniform_propagator=${kalman_uniform_propagator}"
echo "  kalman measurement sigmas: rphi=${kalman_sigma_rphi_cm} cm, r=${kalman_sigma_r_cm} cm, z=${kalman_sigma_z_cm} cm"
echo "  write_kalman_innovation_diagnostics=${write_kalman_innovation_diagnostics}"
echo "  PCA search: coarse_steps=${coarse_steps}, candidates=${pca_candidates}"

condor_submit \
  -append "campaign = ${campaign}" \
  -append "input_dst_filelist = ${dst_filelist}" \
  -append "files_per_job = ${files_per_job}" \
  -append "events_per_input_file = ${events_per_input_file}" \
  -append "total_files = ${total_files}" \
  -append "n_jobs = ${n_jobs}" \
  -append "output_base_dir = ${output_base_dir}" \
  -append "input_path_base = ${input_path_base}" \
  -append "fit_method = ${fit_method}" \
  -append "point_order = ${point_order}" \
  -append "kalman_sigma_rphi_cm = ${kalman_sigma_rphi_cm}" \
  -append "kalman_sigma_r_cm = ${kalman_sigma_r_cm}" \
  -append "kalman_sigma_z_cm = ${kalman_sigma_z_cm}" \
  -append "write_kalman_innovation_diagnostics = ${write_kalman_innovation_diagnostics}" \
  -append "kalman_uniform_propagator = ${kalman_uniform_propagator}" \
  -append "bfield_t = ${bfield_t}" \
  -append "coarse_steps = ${coarse_steps}" \
  -append "pca_candidates = ${pca_candidates}" \
  tpc_truth_v0.job
