#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_compare_v0qa.sh HELIX_ROOT KALMAN_ROOT OUT_PREFIX [HELIX_LABEL] [KALMAN_LABEL] [CUT]

Examples:
  ./run_compare_v0qa.sh helix.root kalman.root qa_helix_vs_kalman

  ./run_compare_v0qa.sh helix.root kalman.root qa_cut helix kalman \
    'pairDCA<0.15 && cosThetaReco>0.94'

Outputs:
  OUT_PREFIX.pdf
  OUT_PREFIX.root
  OUT_PREFIX_summary.txt
EOF
}

if [[ $# -lt 3 || $# -gt 6 ]]; then
  usage
  exit 1
fi

file_a=$1
file_b=$2
out_prefix=$3
label_a=${4:-helix}
label_b=${5:-kalman}
cut=${6:-}

macro_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

root_escape() {
  local text=$1
  text=${text//\\/\\\\}
  text=${text//\"/\\\"}
  printf '%s' "$text"
}

root -l -b -q "${macro_dir}/compareV0QA.C(\"$(root_escape "$file_a")\", \"$(root_escape "$file_b")\", \"$(root_escape "$label_a")\", \"$(root_escape "$label_b")\", \"$(root_escape "$out_prefix")\", \"$(root_escape "$cut")\")"
