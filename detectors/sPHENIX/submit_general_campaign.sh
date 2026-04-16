#!/bin/bash
set -e

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 <campaign> [n_jobs] [nevents_per_job]"
  echo "Example: $0 pionplus_pt10 10000 10"
  exit 1
fi

campaign="$1"
n_jobs="${2:-10000}"
nevents="${3:-10}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

mkdir -p "log/${campaign}" "log/tmp/${campaign}"
mkdir -p "output/${campaign}" "output/${campaign}/completed"
# safety fallback in case submit variables are not applied for any reason
mkdir -p "log/default" "log/tmp/default"

echo "Submitting campaign:"
echo "  campaign=${campaign}"
echo "  n_jobs=${n_jobs}"
echo "  nevents/job=${nevents}"
echo "  macro=Fun4All_G4_sPHENIX.C"
echo "  generator=current branch SIMPLE setup"

condor_submit \
  -append "campaign = ${campaign}" \
  -append "n_jobs = ${n_jobs}" \
  -append "nevents = ${nevents}" \
  -append "output_base_dir = output" \
  general_sim.job
