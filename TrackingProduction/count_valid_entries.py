#!/usr/bin/env python
import ROOT, math


def count_valid_entries(
    file_path,
    tree_name,
    t_branch,
    phi_branch,
    check_phi_range=True,
    phi_min=-5,
    phi_max=5,
):
    """
    Open the ROOT file, get the tree, and count entries with valid t and phi.

    Parameters:
      file_path (str): Path to the ROOT file.
      tree_name (str): Name of the tree to check.
      t_branch (str): Name of the branch for time (t).
      phi_branch (str): Name of the branch for phi.
      check_phi_range (bool): If True, check that phi is between phi_min and phi_max.
      phi_min (float): Minimum allowed phi value.
      phi_max (float): Maximum allowed phi value.

    Returns:
      (valid_count, total_entries) tuple
    """
    # Open the file
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print("Error: Could not open file:", file_path)
        return None, 0

    # Get the tree
    tree = f.Get(tree_name)
    if not tree:
        print("Error: Could not find tree", tree_name, "in file", file_path)
        return None, 0

    total_entries = tree.GetEntries()
    valid_count = 0

    # Loop over entries
    for i in range(total_entries):
        tree.GetEntry(i)
        # Get t and phi values from the tree using getattr
        t_val = getattr(tree, t_branch)
        phi_val = getattr(tree, phi_branch)

        valid_t = True
        # Check that t is a float and is not NaN/infinite
        if isinstance(t_val, float):
            if math.isnan(t_val) or math.isinf(t_val):
                valid_t = False

        valid_phi = True
        try:
            # Convert phi to float regardless of its underlying type
            phi_val_float = float(phi_val)
        except Exception:
            valid_phi = False
        else:

            if check_phi_range:
                if phi_val_float < phi_min or phi_val_float > phi_max:
                    valid_phi = False

        if valid_t and valid_phi:
            valid_count += 1

    f.Close()
    return valid_count, total_entries


def main():
    # Set file paths (update these paths as needed)
    clusterizer_file = "/sphenix/tg/tg01/hf/dcxchenxi/kshort_reco/output4/outputFile_53534_0_clusterizer_edgeOff_staticOff_acts_0218.root"
    combined_file = "/sphenix/user/dcxchenxi/develope/macros/TrackingProduction/combined_53534_0_edgeOff_staticOff_acts_0218.root"

    # Tree names
    clusterizer_tree = "hitTree"  # tree name in the clusterizer file
    combined_tree = "combined_hits"  # tree name in the combined file

    # Branch names for t and phi:
    # For the clusterizer file's hit tree, the code uses branch "t" and "iphi".
    # For the combined file's hit tree, the code uses "t" and "tpc_iphi".
    clusterizer_t_branch = "t"
    clusterizer_phi_branch = "phi"
    combined_t_branch = "t"
    combined_phi_branch = "tpc_phi"

    # Count valid entries in the clusterizer file's hit tree
    valid_clusterizer, total_clusterizer = count_valid_entries(
        clusterizer_file,
        clusterizer_tree,
        clusterizer_t_branch,
        clusterizer_phi_branch,
        check_phi_range=False,  # Set to True if you wish to enforce a phi range
    )

    # Count valid entries in the combined file's hit tree
    valid_combined, total_combined = count_valid_entries(
        combined_file,
        combined_tree,
        combined_t_branch,
        combined_phi_branch,
        check_phi_range=False,
    )

    print("Clusterizer file hit tree:")
    print(
        "  Valid entries with t and phi:",
        valid_clusterizer,
        "out of",
        total_clusterizer,
    )
    print("Combined file hit tree:")
    print("  Valid entries with t and phi:", valid_combined, "out of", total_combined)


if __name__ == "__main__":
    main()
