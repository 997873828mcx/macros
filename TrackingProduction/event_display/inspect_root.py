#!/usr/bin/env python3

import uproot
import numpy as np
import sys


def inspect_root_file(filename):
    """
    Inspects the structure of a ROOT file, listing branches in 'combined_clusters' and 'combined_hits' trees.
    Checks for the presence and format of the 'event' branch.

    Parameters:
    -----------
    filename : str
        Path to the ROOT file to inspect.
    """
    try:
        file = uproot.open(filename)
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: An error occurred while opening the file: {e}")
        sys.exit(1)

    # Define the trees to inspect
    trees = ["combined_clusters", "combined_hits"]

    for tree_name in trees:
        print(f"\nInspecting tree: '{tree_name}'")
        try:
            tree = file[tree_name]
        except KeyError:
            print(f"  Error: The tree '{tree_name}' does not exist in the file.")
            continue

        # List all branches in the tree
        branches = tree.keys()
        print(f"  Branches ({len(branches)}):")
        for branch in branches:
            print(f"    - {branch}")

        # Check for 'event' branch
        if "event" in branches:
            print("  'event' branch found.")

            # Attempt to read the 'event' array via tree.arrays(...)
            try:
                arrays_dict = tree.arrays(
                    ["event"], library="np"
                )  # Retrieve only 'event'
                event_data = arrays_dict["event"]

                # Determine if 'event' data is numeric or string-based
                if np.issubdtype(event_data.dtype, np.number):
                    print(f"    'event' data type: Numeric (dtype: {event_data.dtype})")
                    unique_events = np.unique(event_data)
                    print(f"    Unique event IDs: {unique_events}")
                else:
                    # If strings like "Event 3", extract or list them
                    print(
                        f"    'event' data type: Non-numeric (dtype: {event_data.dtype})"
                    )
                    # If it's bytes, decode
                    if event_data.dtype.type is np.bytes_:
                        event_data = np.array([e.decode("utf-8") for e in event_data])
                    unique_events = np.unique(event_data)
                    print(f"    Unique event labels: {unique_events}")

            except Exception as e:
                print(f"    Error accessing 'event' data: {e}")
        else:
            print("  'event' branch NOT found.")

    print("\nInspection complete.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python inspect_root.py <path_to_root_file>")
        sys.exit(1)

    root_filename = sys.argv[1]
    inspect_root_file(root_filename)
