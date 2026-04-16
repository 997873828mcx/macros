#!/bin/bash
set -e

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <campaign> [merged_output]"
  echo "Example: $0 pionplus_pt10 output/pionplus_pt10_merged_reco_pt.root"
  exit 1
fi

campaign="$1"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

source /opt/sphenix/core/bin/sphenix_setup.sh new

input_dir="output/${campaign}/completed"
output_file="${2:-output/${campaign}_merged_reco_pt.root}"

shopt -s nullglob
tree_files=("${input_dir}"/*_reco_pt.root)
shopt -u nullglob

if [[ ${#tree_files[@]} -eq 0 ]]; then
  echo "No reco pt tree files found in ${input_dir}" >&2
  exit 1
fi

echo "Merging ${#tree_files[@]} reco pt tree files into ${output_file}"
hadd -f "${output_file}" "${tree_files[@]}"
echo "Wrote ${output_file}"
