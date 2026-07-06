#!/bin/bash
set -eo pipefail

# Setup environment
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
export PYTHONUSERBASE=/sphenix/user/dcxchenxi/scratch/pip
source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dcxchenxi/install/
set -u

python_version="$(python - <<'PY'
import sys
print(f"{sys.version_info.major}.{sys.version_info.minor}")
PY
)"
user_site="${PYTHONUSERBASE}/lib/python${python_version}/site-packages"
if [[ -d "${user_site}" ]]; then
  export PYTHONPATH="${user_site}:${PYTHONPATH:-}"
fi

process_id=${1?Error: process ID is not given}
files_per_job=${2:-10}
total_files=${3:-0}

campaign_tag=${CAMPAIGN_TAG:-ap_truthpoints}
output_base_dir=${OUTPUT_BASE_DIR:-output}
completed_dir=${COMPLETED_OUTPUT_DIR:-}
input_list=${AP_INPUT_LIST?Error: AP_INPUT_LIST is not set}
ap_script=${AP_SCRIPT:-/sphenix/user/dcxchenxi/develope/macros/TrackingProduction/make_truth_ap.py}

use_truth_primary_vertex=${USE_TRUTH_PRIMARY_VERTEX:-true}
fit_helix=${FIT_HELIX:-true}
fit_kalman=${FIT_KALMAN:-false}
theta_extension=${THETA_EXTENSION:-2.0}
coarse_steps=${COARSE_STEPS:-64}
pca_search=${PCA_SEARCH:-upstream}
point_order=${POINT_ORDER:-path}
downstream_margin=${DOWNSTREAM_MARGIN:-0.2}
kalman_max_upstream_cm=${KALMAN_MAX_UPSTREAM_CM:-80.0}
kalman_downstream_margin_cm=${KALMAN_DOWNSTREAM_MARGIN_CM:-5.0}
fit_first_points=${FIT_FIRST_POINTS:-8}
pca_candidates=${PCA_CANDIDATES:-32}
prefer_positive_pointing=${PREFER_POSITIVE_POINTING:-false}
same_parent_only=${AP_SAME_PARENT_ONLY:-false}
parent_pids=${AP_PARENT_PIDS:-none}
write_chunk_size=${WRITE_CHUNK_SIZE:-50000}
progress_every=${PROGRESS_EVERY:-0}

kalman_measurement_preset=${KALMAN_MEASUREMENT_PRESET:-custom}
kalman_meas_sigma_xy=${KALMAN_MEAS_SIGMA_XY:-0.03}
kalman_meas_sigma_z=${KALMAN_MEAS_SIGMA_Z:-0.05}
kalman_process_sigma_pos=${KALMAN_PROCESS_SIGMA_POS:-1.0e-4}
kalman_process_sigma_phi=${KALMAN_PROCESS_SIGMA_PHI:-1.0e-5}
kalman_process_sigma_qop=${KALMAN_PROCESS_SIGMA_QOP:-1.0e-6}
kalman_process_sigma_tanl=${KALMAN_PROCESS_SIGMA_TANL:-1.0e-6}
kalman_material_x0_per_cm=${KALMAN_MATERIAL_X0_PER_CM:-0.0}
kalman_ms_scale=${KALMAN_MS_SCALE:-1.0}
kalman_energy_loss_gev_per_cm=${KALMAN_ENERGY_LOSS_GEV_PER_CM:-0.0}
kalman_energy_loss_sigma_fraction=${KALMAN_ENERGY_LOSS_SIGMA_FRACTION:-0.0}

pre_track_pt_min=${AP_PRE_TRACK_PT_MIN:-0.2}
pre_track_dca_xy_min=${AP_PRE_TRACK_DCA_XY_MIN:-0.03}
pre_track_dca_z_min=${AP_PRE_TRACK_DCA_Z_MIN:-none}
pre_track_dca_xy_max=${AP_PRE_TRACK_DCA_XY_MAX:-none}
pre_track_dca_z_max=${AP_PRE_TRACK_DCA_Z_MAX:-none}
pre_pair_dca_max=${AP_PRE_PAIR_DCA_MAX:-5.0}
pre_lproj_min=${AP_PRE_LPROJ_MIN:-0.2}
pre_cos_theta_min=${AP_PRE_COS_THETA_MIN:-none}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ "${input_list}" != /* ]]; then
  input_list="${SCRIPT_DIR}/${input_list}"
fi
if [[ "${ap_script}" != /* ]]; then
  ap_script="${SCRIPT_DIR}/${ap_script}"
fi
if [[ "${output_base_dir}" != /* ]]; then
  output_base_dir="${SCRIPT_DIR}/${output_base_dir}"
fi

if [[ ! -f "${input_list}" ]]; then
  echo "Error: AP input list not found: ${input_list}" >&2
  exit 2
fi
if [[ ! -f "${ap_script}" ]]; then
  echo "Error: make_truth_ap.py not found: ${ap_script}" >&2
  exit 2
fi
if [[ "${files_per_job}" -le 0 ]]; then
  echo "Error: files_per_job must be positive, got ${files_per_job}" >&2
  exit 2
fi
if [[ "${total_files}" -lt 0 ]]; then
  echo "Error: total_files must be non-negative, got ${total_files}" >&2
  exit 2
fi

start_index=$((process_id * files_per_job))
requested_files=${files_per_job}
if [[ "${total_files}" -gt 0 ]]; then
  remaining=$((total_files - start_index))
  if [[ "${remaining}" -le 0 ]]; then
    echo "Nothing to do: process=${process_id} start_index=${start_index} total_files=${total_files}"
    exit 0
  fi
  if [[ "${remaining}" -lt "${requested_files}" ]]; then
    requested_files=${remaining}
  fi
fi

start_line=$((start_index + 1))
end_line=$((start_index + requested_files))
padded_id=$(printf "%05d" "${process_id}")
job_outdir="${output_base_dir}/${campaign_tag}"
split_dir="${job_outdir}/splits"
mkdir -p "${job_outdir}" "${split_dir}"
job_input_list="${split_dir}/input_${padded_id}_f${start_index}_n${requested_files}.list"

resolve_list_entry()
{
  local entry="$1"
  if [[ "${entry}" == /* ]]; then
    echo "${entry}"
    return 0
  fi

  local list_dir
  list_dir="$(cd "$(dirname "${input_list}")" && pwd -P)"
  local bases=("${SCRIPT_DIR}" "${list_dir}")
  local parent="${list_dir}"
  while [[ "${parent}" != "/" ]]; do
    bases+=("${parent}")
    parent="$(dirname "${parent}")"
  done
  bases+=("/")

  local base
  for base in "${bases[@]}"; do
    if [[ -e "${base}/${entry}" ]]; then
      local base_abs
      base_abs="$(cd "${base}" && pwd -P)"
      echo "${base_abs}/${entry}"
      return 0
    fi
  done

  echo "${list_dir}/${entry}"
}

: > "${job_input_list}"
sed -n "${start_line},${end_line}p" "${input_list}" | while IFS= read -r line; do
  [[ -z "${line}" || "${line}" =~ ^[[:space:]]*# ]] && continue
  resolve_list_entry "${line}" >> "${job_input_list}"
done

n_selected=$(grep -cv '^[[:space:]]*$' "${job_input_list}" || true)
if [[ "${n_selected}" -le 0 ]]; then
  echo "Nothing to do: no input files selected for process=${process_id}"
  exit 0
fi

outroot="${job_outdir}/ap_pairs_${padded_id}_f${start_index}_n${n_selected}.root"
if [[ -z "${completed_dir}" ]]; then
  completed_dir="${job_outdir}/completed"
fi

add_optional_float_arg()
{
  local arg_name="$1"
  local value="$2"
  if [[ -n "${value}" && "${value}" != "none" ]]; then
    cmd+=("${arg_name}" "${value}")
  fi
}

cmd=(python "${ap_script}"
  --input-list "${job_input_list}"
  -o "${outroot}"
  --theta-extension "${theta_extension}"
  --coarse-steps "${coarse_steps}"
  --pca-search "${pca_search}"
  --point-order "${point_order}"
  --downstream-margin "${downstream_margin}"
  --kalman-max-upstream-cm "${kalman_max_upstream_cm}"
  --kalman-downstream-margin-cm "${kalman_downstream_margin_cm}"
  --fit-first-points "${fit_first_points}"
  --pca-candidates "${pca_candidates}"
  --write-chunk-size "${write_chunk_size}"
  --progress-every "${progress_every}"
  --kalman-measurement-preset "${kalman_measurement_preset}"
)

if [[ "${use_truth_primary_vertex}" == "true" ]]; then
  cmd+=(--use-truth-primary-vertex)
elif [[ "${use_truth_primary_vertex}" != "false" ]]; then
  echo "Error: USE_TRUTH_PRIMARY_VERTEX must be true or false, got ${use_truth_primary_vertex}" >&2
  exit 2
fi

if [[ "${fit_helix}" == "true" ]]; then
  cmd+=(--fit-helix)
elif [[ "${fit_helix}" != "false" ]]; then
  echo "Error: FIT_HELIX must be true or false, got ${fit_helix}" >&2
  exit 2
fi

if [[ "${fit_kalman}" == "true" ]]; then
  cmd+=(--fit-kalman)
elif [[ "${fit_kalman}" != "false" ]]; then
  echo "Error: FIT_KALMAN must be true or false, got ${fit_kalman}" >&2
  exit 2
fi

cmd+=(
  --kalman-meas-sigma-xy "${kalman_meas_sigma_xy}"
  --kalman-meas-sigma-z "${kalman_meas_sigma_z}"
  --kalman-process-sigma-pos "${kalman_process_sigma_pos}"
  --kalman-process-sigma-phi "${kalman_process_sigma_phi}"
  --kalman-process-sigma-qop "${kalman_process_sigma_qop}"
  --kalman-process-sigma-tanl "${kalman_process_sigma_tanl}"
  --kalman-material-x0-per-cm "${kalman_material_x0_per_cm}"
  --kalman-ms-scale "${kalman_ms_scale}"
  --kalman-energy-loss-gev-per-cm "${kalman_energy_loss_gev_per_cm}"
  --kalman-energy-loss-sigma-fraction "${kalman_energy_loss_sigma_fraction}"
)

if [[ "${prefer_positive_pointing}" == "true" ]]; then
  cmd+=(--prefer-positive-pointing)
elif [[ "${prefer_positive_pointing}" != "false" ]]; then
  echo "Error: PREFER_POSITIVE_POINTING must be true or false, got ${prefer_positive_pointing}" >&2
  exit 2
fi

if [[ "${same_parent_only}" == "true" ]]; then
  cmd+=(--same-parent-only)
elif [[ "${same_parent_only}" != "false" ]]; then
  echo "Error: AP_SAME_PARENT_ONLY must be true or false, got ${same_parent_only}" >&2
  exit 2
fi

if [[ "${parent_pids}" != "none" && -n "${parent_pids}" ]]; then
  IFS=',' read -r -a parent_pid_array <<< "${parent_pids}"
  cmd+=(--parent-pids "${parent_pid_array[@]}")
fi

add_optional_float_arg --pre-track-pt-min "${pre_track_pt_min}"
add_optional_float_arg --pre-track-dca-xy-min "${pre_track_dca_xy_min}"
add_optional_float_arg --pre-track-dca-z-min "${pre_track_dca_z_min}"
add_optional_float_arg --pre-track-dca-xy-max "${pre_track_dca_xy_max}"
add_optional_float_arg --pre-track-dca-z-max "${pre_track_dca_z_max}"
add_optional_float_arg --pre-pair-dca-max "${pre_pair_dca_max}"
add_optional_float_arg --pre-lproj-min "${pre_lproj_min}"
add_optional_float_arg --pre-cos-theta-min "${pre_cos_theta_min}"

echo "Running truth AP analysis:"
echo "  campaign=${campaign_tag}"
echo "  job=${process_id}"
echo "  input_list=${input_list}"
echo "  selected_files=${n_selected}"
echo "  first_file_index=${start_index}"
echo "  job_input_list=${job_input_list}"
echo "  out=${outroot}"
echo "  point_order=${point_order}"
echo "  preselection: pt>${pre_track_pt_min}, dca_xy>${pre_track_dca_xy_min}, pairDCA<${pre_pair_dca_max}, Lproj>${pre_lproj_min}"
printf '  command:'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

mkdir -p "${completed_dir}"
base="$(basename "${outroot}")"
dst="${completed_dir}/${base}"
if [[ -e "${dst}" ]]; then
  stem="${base%.*}"
  ext="${base##*.}"
  dst="${completed_dir}/${stem}__dup_$(date +%s)_${process_id}.${ext}"
  echo "Warning: completed output exists, moving as ${dst}"
fi
mv -f "${outroot}" "${dst}"

echo "Done: ${dst}"
