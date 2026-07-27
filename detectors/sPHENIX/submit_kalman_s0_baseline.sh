#!/bin/bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <hepmc|g4|fit|status>" >&2
  exit 1
fi

stage="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# These defaults define the frozen S0 sample. Overrides are recorded in the
# manifest, but should not be used when comparing fitter-development stages.
label="${S0_LABEL:-pp_minbias_kalman_s0_10k_new15_uniform1p4}"
total_events="${S0_TOTAL_EVENTS:-10000}"
hepmc_events_per_file="${S0_HEPMC_EVENTS_PER_FILE:-1000}"
g4_events_per_job="${S0_G4_EVENTS_PER_JOB:-10}"
fit_files_per_job="${S0_FIT_FILES_PER_JOB:-10}"
pythia_seed_base="${S0_PYTHIA_SEED_BASE:-730000}"
g4_seed_base="${S0_G4_SEED_BASE:-830000}"
uniform_bfield_t="${S0_UNIFORM_BFIELD_T:-1.4}"

sigma_rphi_cm="${S0_KALMAN_SIGMA_RPHI_CM:-0.08}"
# A cylindrical layer fixes r. Keep a small numerical floor in this 3D
# measurement implementation; S0 applies no random radial displacement.
sigma_r_cm="${S0_KALMAN_SIGMA_R_CM:-1.0e-4}"
sigma_z_cm="${S0_KALMAN_SIGMA_Z_CM:-0.23}"
coarse_steps="${S0_COARSE_STEPS:-64}"
pca_candidates="${S0_PCA_CANDIDATES:-32}"

hepmc_campaign="${label}_hepmc"
g4_campaign="${label}_truthpoints"
fit_campaign="${label}_fit"
metadata_dir="output/${label}_metadata"
manifest="${metadata_dir}/manifest.env"
pythia_config="${metadata_dir}/phpythia8_detroit_minBias.cfg"
hepmc_list="output/${hepmc_campaign}/hepmc_files.list"
dst_list="${metadata_dir}/truth_point_dst_files.list"

for value in "${total_events}" "${hepmc_events_per_file}" "${g4_events_per_job}" \
             "${fit_files_per_job}" "${pythia_seed_base}" "${g4_seed_base}"; do
  if ! [[ "${value}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: event, job, and seed settings must be positive integers; got ${value}" >&2
    exit 2
  fi
done
if ((total_events % hepmc_events_per_file != 0)); then
  echo "Error: S0_TOTAL_EVENTS must be divisible by S0_HEPMC_EVENTS_PER_FILE" >&2
  exit 2
fi
if ((hepmc_events_per_file % g4_events_per_job != 0)); then
  echo "Error: S0_G4_EVENTS_PER_JOB must divide S0_HEPMC_EVENTS_PER_FILE" >&2
  exit 2
fi

mkdir -p "${metadata_dir}"

if [[ ! -f "${pythia_config}" ]]; then
  if [[ -z "${CALIBRATIONROOT:-}" ]]; then
    set +u
    source /opt/sphenix/core/bin/sphenix_setup.sh -n new.15
    set -u
  fi
  source_config="${CALIBRATIONROOT}/Generators/phpythia8_detroit_minBias.cfg"
  if [[ ! -f "${source_config}" ]]; then
    echo "Error: Pythia configuration not found: ${source_config}" >&2
    exit 2
  fi
  cp -p "${source_config}" "${pythia_config}"
fi

write_manifest()
{
  if [[ -f "${manifest}" && "${S0_REFRESH_MANIFEST:-false}" != "true" ]]; then
    return
  fi
  local macros_revision coresoftware_revision
  macros_revision="$(git -C /sphenix/user/dcxchenxi/develope/macros rev-parse HEAD)"
  coresoftware_revision="$(git -C /sphenix/user/dcxchenxi/develope/coresoftware rev-parse HEAD)"
  {
    echo "format_version=1"
    echo "created_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "label=${label}"
    echo "software_release=new.15"
    echo "total_events=${total_events}"
    echo "hepmc_events_per_file=${hepmc_events_per_file}"
    echo "g4_events_per_job=${g4_events_per_job}"
    echo "fit_files_per_job=${fit_files_per_job}"
    echo "pythia_seed_base=${pythia_seed_base}"
    echo "g4_seed_base=${g4_seed_base}"
    echo "use_beam_vertex=false"
    echo "g4_uniform_bfield_t=${uniform_bfield_t}"
    echo "truth_point_node=G4HIT_TPC_TRUECLUSTER"
    echo "truth_point_smearing=none"
    echo "fit_method=kalman"
    echo "point_order=path"
    echo "kalman_uniform_propagator=analytic"
    echo "kalman_sigma_rphi_cm=${sigma_rphi_cm}"
    echo "kalman_sigma_r_cm=${sigma_r_cm}"
    echo "kalman_sigma_z_cm=${sigma_z_cm}"
    echo "kalman_innovation_diagnostics=true"
    echo "coarse_steps=${coarse_steps}"
    echo "pca_candidates=${pca_candidates}"
    echo "hepmc_campaign=${hepmc_campaign}"
    echo "g4_campaign=${g4_campaign}"
    echo "fit_campaign=${fit_campaign}"
    echo "macros_git_revision=${macros_revision}"
    echo "coresoftware_git_revision=${coresoftware_revision}"
    echo "pythia_config=${SCRIPT_DIR}/${pythia_config}"
    echo "pythia_config_sha256=$(sha256sum "${pythia_config}" | awk '{print $1}')"
    echo "generate_macro_sha256=$(sha256sum GeneratePythia8HepMC.C | awk '{print $1}')"
    echo "g4_macro_sha256=$(sha256sum Fun4All_G4_sPHENIX.C | awk '{print $1}')"
    echo "truth_v0_macro_sha256=$(sha256sum Fun4All_TpcTruthV0.C | awk '{print $1}')"
    echo "general_sim_runner_sha256=$(sha256sum run_general_sim.sh | awk '{print $1}')"
    echo "truth_v0_runner_sha256=$(sha256sum run_tpc_truth_v0.sh | awk '{print $1}')"
    echo "s0_driver_sha256=$(sha256sum submit_kalman_s0_baseline.sh | awk '{print $1}')"
    echo "kalman_source_sha256=$(sha256sum /sphenix/user/dcxchenxi/develope/coresoftware/simulation/g4simulation/g4tpc/TpcTrackKalmanFitter.cc | awk '{print $1}')"
  } > "${manifest}"
}

append_history()
{
  printf '%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" >> "${metadata_dir}/history.log"
}

count_files()
{
  local directory="$1"
  local pattern="$2"
  if [[ ! -d "${directory}" ]]; then
    echo 0
    return
  fi
  find "${directory}" -maxdepth 1 -type f -name "${pattern}" | wc -l
}

expected_hepmc=$((total_events / hepmc_events_per_file))
expected_dst=$((total_events / g4_events_per_job))

case "${stage}" in
  hepmc)
    write_manifest
    append_history "submit HepMC generation"
    ./submit_hepmc_generation_campaign.sh \
      "${hepmc_campaign}" \
      "${total_events}" \
      "${hepmc_events_per_file}" \
      "${pythia_seed_base}" \
      "${SCRIPT_DIR}/${pythia_config}"
    ;;

  g4)
    hepmc_done=$(count_files "output/${hepmc_campaign}/completed" 'pp_minbias_*.hepmc')
    if [[ "${hepmc_done}" -ne "${expected_hepmc}" ]]; then
      echo "Error: expected ${expected_hepmc} completed HepMC files, found ${hepmc_done}" >&2
      exit 3
    fi
    write_manifest
    sha256sum "${hepmc_list}" > "${metadata_dir}/hepmc_filelist.sha256"
    append_history "submit deterministic uniform-field G4 truth-point transport"
    G4_SEED_BASE="${g4_seed_base}" \
    G4_UNIFORM_BFIELD_T="${uniform_bfield_t}" \
    WRITE_TRUTH_FLAT_TREE=false \
    ./submit_general_campaign.sh \
      "${g4_campaign}" \
      "${hepmc_list}" \
      "${total_events}" \
      "${g4_events_per_job}" \
      "${hepmc_events_per_file}" \
      false
    ;;

  fit)
    find "output/${g4_campaign}/completed" -maxdepth 1 -type f \
      -name 'tpc_truthpoints_dst_*.root' | LC_ALL=C sort > "${dst_list}"
    dst_done=$(wc -l < "${dst_list}")
    if [[ "${dst_done}" -ne "${expected_dst}" ]]; then
      echo "Error: expected ${expected_dst} completed truth-point DSTs, found ${dst_done}" >&2
      exit 3
    fi
    write_manifest
    sha256sum "${dst_list}" > "${metadata_dir}/truth_point_dst_filelist.sha256"
    append_history "submit unsmeared S0 Kalman diagnostics"
    V0_FIT_METHOD=kalman \
    V0_POINT_ORDER=path \
    V0_KALMAN_SIGMA_RPHI_CM="${sigma_rphi_cm}" \
    V0_KALMAN_SIGMA_R_CM="${sigma_r_cm}" \
    V0_KALMAN_SIGMA_Z_CM="${sigma_z_cm}" \
    V0_WRITE_KALMAN_INNOVATION_DIAGNOSTICS=true \
    V0_KALMAN_UNIFORM_PROPAGATOR=analytic \
    V0_BFIELD_T="${uniform_bfield_t}" \
    V0_COARSE_STEPS="${coarse_steps}" \
    V0_PCA_CANDIDATES="${pca_candidates}" \
    ./submit_tpc_truth_v0_campaign.sh \
      "${fit_campaign}" \
      "${dst_list}" \
      "${fit_files_per_job}" \
      "${g4_events_per_job}" \
      "${expected_dst}"
    ;;

  status)
    write_manifest
    hepmc_done=$(count_files "output/${hepmc_campaign}/completed" 'pp_minbias_*.hepmc')
    dst_done=$(count_files "output/${g4_campaign}/completed" 'tpc_truthpoints_dst_*.root')
    fit_done=$(count_files "output/${fit_campaign}/completed" 'tpc_v0_candidates_*.root')
    expected_fit=$(((expected_dst + fit_files_per_job - 1) / fit_files_per_job))
    echo "S0 baseline: ${label}"
    echo "  HepMC: ${hepmc_done}/${expected_hepmc} files"
    echo "  truth DST: ${dst_done}/${expected_dst} files"
    echo "  Kalman QA: ${fit_done}/${expected_fit} files"
    echo "  manifest: ${SCRIPT_DIR}/${manifest}"
    ;;

  *)
    echo "Error: stage must be hepmc, g4, fit, or status" >&2
    exit 1
    ;;
esac
