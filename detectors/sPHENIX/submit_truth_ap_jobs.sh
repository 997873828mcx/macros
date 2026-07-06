#!/bin/bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
  echo "Usage: $0 <campaign> <truthpoint_root_file_list> [files_per_job] [max_files]"
  echo "Example: $0 pp_minbias_ap_100k_allpairs output/pp_minbias_truth_100k_zero_vtx/root_files.list 10"
  echo
  echo "Optional environment overrides:"
  echo "  AP_PRE_TRACK_PT_MIN=0.2"
  echo "  AP_PRE_TRACK_DCA_XY_MIN=0.03"
  echo "  AP_PRE_PAIR_DCA_MAX=5.0"
  echo "  AP_PRE_LPROJ_MIN=0.2"
  echo "  FIT_HELIX=false FIT_KALMAN=true"
  echo "  KALMAN_MEASUREMENT_PRESET=truth|realistic|custom"
  echo "  POINT_ORDER=path|input|radius|theta-z|auto"
  echo "  COARSE_STEPS=64"
  echo "  PCA_CANDIDATES=32"
  exit 1
fi

campaign="$1"
input_root_list="$2"
files_per_job="${3:-10}"
max_files="${4:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

if [[ "${input_root_list}" != /* ]]; then
  input_root_list="${SCRIPT_DIR}/${input_root_list}"
fi

if [[ ! -f "${input_root_list}" ]]; then
  echo "Error: input ROOT file list not found: ${input_root_list}" >&2
  exit 2
fi
if [[ "${files_per_job}" -le 0 ]]; then
  echo "Error: files_per_job must be positive, got ${files_per_job}" >&2
  exit 2
fi
if [[ "${max_files}" -lt 0 ]]; then
  echo "Error: max_files must be non-negative, got ${max_files}" >&2
  exit 2
fi

available_files=$(grep -cv '^[[:space:]]*$' "${input_root_list}")
if [[ "${available_files}" -le 0 ]]; then
  echo "Error: input ROOT file list is empty: ${input_root_list}" >&2
  exit 2
fi

n_files="${available_files}"
if [[ "${max_files}" -gt 0 && "${max_files}" -lt "${n_files}" ]]; then
  n_files="${max_files}"
fi
n_jobs=$(((n_files + files_per_job - 1) / files_per_job))

ap_script="${AP_SCRIPT:-/sphenix/user/dcxchenxi/develope/macros/TrackingProduction/make_truth_ap.py}"
output_base_dir="${OUTPUT_BASE_DIR:-output}"
use_truth_primary_vertex="${USE_TRUTH_PRIMARY_VERTEX:-true}"
fit_helix="${FIT_HELIX:-true}"
fit_kalman="${FIT_KALMAN:-false}"
theta_extension="${THETA_EXTENSION:-2.0}"
coarse_steps="${COARSE_STEPS:-64}"
pca_search="${PCA_SEARCH:-upstream}"
point_order="${POINT_ORDER:-path}"
downstream_margin="${DOWNSTREAM_MARGIN:-0.2}"
kalman_max_upstream_cm="${KALMAN_MAX_UPSTREAM_CM:-80.0}"
kalman_downstream_margin_cm="${KALMAN_DOWNSTREAM_MARGIN_CM:-5.0}"
fit_first_points="${FIT_FIRST_POINTS:-8}"
pca_candidates="${PCA_CANDIDATES:-32}"
prefer_positive_pointing="${PREFER_POSITIVE_POINTING:-false}"
same_parent_only="${AP_SAME_PARENT_ONLY:-false}"
parent_pids="${AP_PARENT_PIDS:-none}"
write_chunk_size="${WRITE_CHUNK_SIZE:-50000}"
progress_every="${PROGRESS_EVERY:-0}"

kalman_measurement_preset="${KALMAN_MEASUREMENT_PRESET:-custom}"
kalman_meas_sigma_xy="${KALMAN_MEAS_SIGMA_XY:-0.03}"
kalman_meas_sigma_z="${KALMAN_MEAS_SIGMA_Z:-0.05}"
kalman_process_sigma_pos="${KALMAN_PROCESS_SIGMA_POS:-1.0e-4}"
kalman_process_sigma_phi="${KALMAN_PROCESS_SIGMA_PHI:-1.0e-5}"
kalman_process_sigma_qop="${KALMAN_PROCESS_SIGMA_QOP:-1.0e-6}"
kalman_process_sigma_tanl="${KALMAN_PROCESS_SIGMA_TANL:-1.0e-6}"
kalman_material_x0_per_cm="${KALMAN_MATERIAL_X0_PER_CM:-0.0}"
kalman_ms_scale="${KALMAN_MS_SCALE:-1.0}"
kalman_energy_loss_gev_per_cm="${KALMAN_ENERGY_LOSS_GEV_PER_CM:-0.0}"
kalman_energy_loss_sigma_fraction="${KALMAN_ENERGY_LOSS_SIGMA_FRACTION:-0.0}"

pre_track_pt_min="${AP_PRE_TRACK_PT_MIN:-0.2}"
pre_track_dca_xy_min="${AP_PRE_TRACK_DCA_XY_MIN:-0.03}"
pre_track_dca_z_min="${AP_PRE_TRACK_DCA_Z_MIN:-none}"
pre_track_dca_xy_max="${AP_PRE_TRACK_DCA_XY_MAX:-none}"
pre_track_dca_z_max="${AP_PRE_TRACK_DCA_Z_MAX:-none}"
pre_pair_dca_max="${AP_PRE_PAIR_DCA_MAX:-5.0}"
pre_lproj_min="${AP_PRE_LPROJ_MIN:-0.2}"
pre_cos_theta_min="${AP_PRE_COS_THETA_MIN:-none}"

mkdir -p "log/${campaign}" "log/tmp/${campaign}"
mkdir -p "output/${campaign}" "output/${campaign}/completed" "output/${campaign}/splits"
mkdir -p "log/default" "log/tmp/default"
chmod +x "${SCRIPT_DIR}/run_truth_ap.sh"

echo "Submitting truth AP campaign:"
echo "  campaign=${campaign}"
echo "  input_root_list=${input_root_list}"
echo "  available_files=${available_files}"
echo "  n_files=${n_files}"
echo "  files/job=${files_per_job}"
echo "  n_jobs=${n_jobs}"
echo "  output=output/${campaign}/completed/"
echo "  point_order=${point_order}"
echo "  preselection: pt>${pre_track_pt_min}, dca_xy>${pre_track_dca_xy_min}, pairDCA<${pre_pair_dca_max}, Lproj>${pre_lproj_min}"

condor_submit \
  -append "campaign = ${campaign}" \
  -append "input_root_list = ${input_root_list}" \
  -append "files_per_job = ${files_per_job}" \
  -append "n_files = ${n_files}" \
  -append "n_jobs = ${n_jobs}" \
  -append "output_base_dir = ${output_base_dir}" \
  -append "ap_script = ${ap_script}" \
  -append "use_truth_primary_vertex = ${use_truth_primary_vertex}" \
  -append "fit_helix = ${fit_helix}" \
  -append "fit_kalman = ${fit_kalman}" \
  -append "theta_extension = ${theta_extension}" \
  -append "coarse_steps = ${coarse_steps}" \
  -append "pca_search = ${pca_search}" \
  -append "point_order = ${point_order}" \
  -append "downstream_margin = ${downstream_margin}" \
  -append "kalman_max_upstream_cm = ${kalman_max_upstream_cm}" \
  -append "kalman_downstream_margin_cm = ${kalman_downstream_margin_cm}" \
  -append "fit_first_points = ${fit_first_points}" \
  -append "pca_candidates = ${pca_candidates}" \
  -append "prefer_positive_pointing = ${prefer_positive_pointing}" \
  -append "same_parent_only = ${same_parent_only}" \
  -append "parent_pids = ${parent_pids}" \
  -append "write_chunk_size = ${write_chunk_size}" \
  -append "progress_every = ${progress_every}" \
  -append "kalman_measurement_preset = ${kalman_measurement_preset}" \
  -append "kalman_meas_sigma_xy = ${kalman_meas_sigma_xy}" \
  -append "kalman_meas_sigma_z = ${kalman_meas_sigma_z}" \
  -append "kalman_process_sigma_pos = ${kalman_process_sigma_pos}" \
  -append "kalman_process_sigma_phi = ${kalman_process_sigma_phi}" \
  -append "kalman_process_sigma_qop = ${kalman_process_sigma_qop}" \
  -append "kalman_process_sigma_tanl = ${kalman_process_sigma_tanl}" \
  -append "kalman_material_x0_per_cm = ${kalman_material_x0_per_cm}" \
  -append "kalman_ms_scale = ${kalman_ms_scale}" \
  -append "kalman_energy_loss_gev_per_cm = ${kalman_energy_loss_gev_per_cm}" \
  -append "kalman_energy_loss_sigma_fraction = ${kalman_energy_loss_sigma_fraction}" \
  -append "pre_track_pt_min = ${pre_track_pt_min}" \
  -append "pre_track_dca_xy_min = ${pre_track_dca_xy_min}" \
  -append "pre_track_dca_z_min = ${pre_track_dca_z_min}" \
  -append "pre_track_dca_xy_max = ${pre_track_dca_xy_max}" \
  -append "pre_track_dca_z_max = ${pre_track_dca_z_max}" \
  -append "pre_pair_dca_max = ${pre_pair_dca_max}" \
  -append "pre_lproj_min = ${pre_lproj_min}" \
  -append "pre_cos_theta_min = ${pre_cos_theta_min}" \
  truth_ap.job
