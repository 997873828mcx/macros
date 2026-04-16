#!/bin/bash
set -e

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <campaign> [merged_output]"
  echo "Example: $0 pionplus_pt10 output/pionplus_pt10_merged_qa.root"
  exit 1
fi

campaign="$1"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

source /opt/sphenix/core/bin/sphenix_setup.sh new

input_dir="output/${campaign}/completed"
output_file="${2:-output/${campaign}_merged_qa.root}"

shopt -s nullglob
qa_files=("${input_dir}"/*_qa.root)
shopt -u nullglob

if [[ ${#qa_files[@]} -eq 0 ]]; then
  echo "No QA files found in ${input_dir}" >&2
  exit 1
fi

echo "Merging ${#qa_files[@]} QA files into ${output_file}"
hadd -f "${output_file}" "${qa_files[@]}"
echo "Wrote ${output_file}"
