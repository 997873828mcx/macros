#!/usr/bin/env python3

"""
Export cluster and hit positions with ADC values from a ROOT file.

This script reads two trees (defaults: clustertree and hittree) that contain
branches gx, gy, gz, and adc, and writes their contents to plain-text files
with one entry per line.
"""

import argparse
import sys
from typing import Iterable, Sequence


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Dump gx/gy/gz/adc from cluster and hit trees into text files."
        )
    )
    parser.add_argument(
        "input_root",
        help="Path to the input ROOT file (e.g. clusters_seeds_53877-0.root_resid.root)",
    )
    parser.add_argument(
        "--cluster-tree",
        default="clustertree",
        help="Name of the TTree holding clusters (default: clustertree)",
    )
    parser.add_argument(
        "--hit-tree",
        default="hittree",
        help="Name of the TTree holding hits (default: hittree)",
    )
    parser.add_argument(
        "--clusters-out",
        default="clusters_xyzadc.txt",
        help="Output text file for cluster entries (default: clusters_xyzadc.txt)",
    )
    parser.add_argument(
        "--hits-out",
        default="hits_xyzadc.txt",
        help="Output text file for hit entries (default: hits_xyzadc.txt)",
    )
    return parser.parse_args()


def ensure_branches(tree, branches: Iterable[str]) -> None:
    missing = [name for name in branches if not tree.GetBranch(name)]
    if missing:
        raise RuntimeError(
            f"TTree '{tree.GetName()}' is missing required branches: {', '.join(missing)}"
        )


def open_tree(root_file, tree_name: str):
    tree = root_file.Get(tree_name)
    if not tree:
        raise RuntimeError(f"Could not find TTree '{tree_name}' in {root_file.GetName()}")
    return tree


def format_value(value) -> str:
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    try:
        float_value = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{float_value:.6f}"


def dump_tree(tree, branches: Sequence[str], output_path: str) -> int:
    entries = tree.GetEntries()
    with open(output_path, "w", encoding="ascii") as output:
        for idx in range(entries):
            tree.GetEntry(idx)
            values = [format_value(getattr(tree, name)) for name in branches]
            output.write(" ".join(values))
            output.write("\n")
    return entries


def main() -> int:
    args = parse_args()

    try:
        import ROOT  # pylint: disable=import-error
    except ImportError as exc:
        raise SystemExit(f"Failed to import ROOT. Ensure PyROOT is available. ({exc})")

    root_file = ROOT.TFile.Open(args.input_root)
    if not root_file or root_file.IsZombie():
        raise SystemExit(f"Could not open ROOT file: {args.input_root}")

    try:
        cluster_tree = open_tree(root_file, args.cluster_tree)
        hit_tree = open_tree(root_file, args.hit_tree)

        required_branches = ["gx", "gy", "gz", "adc"]
        ensure_branches(cluster_tree, required_branches)
        ensure_branches(hit_tree, required_branches)

        cluster_entries = dump_tree(cluster_tree, required_branches, args.clusters_out)
        hit_entries = dump_tree(hit_tree, required_branches, args.hits_out)
    finally:
        root_file.Close()

    print(
        f"Wrote {cluster_entries} cluster entries to '{args.clusters_out}' "
        f"and {hit_entries} hit entries to '{args.hits_out}'."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
