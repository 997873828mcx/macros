import uproot
import numpy as np
import pyvista as pv
from PyQt5.QtWidgets import (
    QMessageBox,
)  # Optional: For error messaging if integrating with GUI


def load_data_from_root(filename):

    try:
        file = uproot.open(filename)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{filename}' was not found.")
    except Exception as e:
        raise IOError(f"An error occurred while opening the file: {e}")

    try:
        cluster_tree = file["combined_clusters"]
    except KeyError:
        raise KeyError("The ROOT file does not contain a 'combined_clusters' tree.")

    try:
        cluster_data = cluster_tree.arrays(library="np")
    except Exception as e:
        raise ValueError(f"Failed to extract arrays from 'combined_clusters' tree: {e}")

    required_cluster_branches = {"gx", "gy", "gz", "side"}
    if not required_cluster_branches.issubset(cluster_data.keys()):
        missing = required_cluster_branches - set(cluster_data.keys())
        raise KeyError(f"Missing branches in 'combined_clusters' tree: {missing}")

    if "event" in cluster_data.keys():
        event_clusters = cluster_data["event"]

        # Determine if 'event' data is numeric or string-based
        if np.issubdtype(event_clusters.dtype, np.number):
            # If numeric, ensure it's integer type
            event_clusters = event_clusters.astype(int)
        else:
            # If strings like "Event 3", extract the numeric part
            event_clusters = np.array(
                [
                    (
                        int(e.split()[1])
                        if isinstance(e, str) and len(e.split()) > 1
                        else 0
                    )
                    for e in event_clusters
                ]
            )

        print(
            f"Loaded 'event' branch for clusters. Unique event IDs: {np.unique(event_clusters)}"
        )
    else:
        # Assign all clusters to event ID 0 if 'event' branch is missing
        print(
            f"Warning: 'event' branch not found in 'combined_clusters' tree. Assigning event=0 to all clusters."
        )
        event_clusters = np.zeros(cluster_data["gx"].shape, dtype=int)
    # Extract cluster coordinates
    cx = np.nan_to_num(cluster_data["gx"], nan=0.0)
    cy = np.nan_to_num(cluster_data["gy"], nan=0.0)
    cz = np.nan_to_num(cluster_data["gz"], nan=0.0)
    cpoints = np.column_stack([cx, cy, cz])

    cluster_polydata = pv.PolyData(cpoints)
    for name in cluster_data.keys():
        if name not in ("gx", "gy", "gz"):
            try:
                if name == "clust_hitkeys":
                    # Convert the vector of hitkeys into a list or object array
                    hitkeys_list = [
                        list(hitkeys) if isinstance(hitkeys, np.ndarray) else []
                        for hitkeys in cluster_data[name]
                    ]
                    cluster_polydata.hitkeys = hitkeys_list
                    print("Hitkeys data verification:")
                    print(f"Number of clusters with hitkeys: {len(hitkeys_list)}")
                    print(f"First 3 clusters' hitkeys: {hitkeys_list[:3]}")
                    print(
                        f"Number of hits in first cluster: {len(hitkeys_list[0]) if hitkeys_list else 0}"
                    )
                else:
                    # Handle other non-vector branches normally
                    cluster_polydata.point_data[name] = cluster_data[name]
            except Exception as e:
                print(f"Warning: Could not add branch {name} to point_data: {e}")

    # Add a 'data_type' attribute to distinguish clusters (0)
    cluster_polydata.point_data["data_type"] = np.zeros(
        cluster_polydata.n_points, dtype=int
    )

    cluster_polydata.point_data["event"] = event_clusters
    try:
        hits_tree = file["combined_hits"]
    except KeyError:
        raise KeyError("The ROOT file does not contain a 'combined_hits' tree.")

    try:
        hits_data = hits_tree.arrays(library="np")
    except Exception as e:
        raise ValueError(f"Failed to extract arrays from 'combined_hits' tree: {e}")

    required_hit_branches = {"gx", "gy", "gz", "side"}
    if not required_hit_branches.issubset(hits_data.keys()):
        missing = required_hit_branches - set(hits_data.keys())
        raise KeyError(f"Missing branches in 'combined_hits' tree: {missing}")

    if "event" in hits_data.keys():
        event_hits = hits_data["event"]

        # Determine if 'event' data is numeric or string-based
        if np.issubdtype(event_hits.dtype, np.number):
            # If numeric, ensure it's integer type
            event_hits = event_hits.astype(int)
        else:
            # If strings like "Event 3", extract the numeric part
            event_hits = np.array(
                [
                    (
                        int(e.split()[1])
                        if isinstance(e, str) and len(e.split()) > 1
                        else 0
                    )
                    for e in event_hits
                ]
            )

        print(
            f"Loaded 'event' branch for hits. Unique event IDs: {np.unique(event_hits)}"
        )
    else:
        # Assign all hits to event ID 0 if 'event' branch is missing
        print(
            f"Warning: 'event' branch not found in 'combined_hits' tree. Assigning event=0 to all hits."
        )
        event_hits = np.zeros(hits_data["gx"].shape, dtype=int)
    hx = np.nan_to_num(hits_data["gx"], nan=0.0)
    hy = np.nan_to_num(hits_data["gy"], nan=0.0)
    hz = np.nan_to_num(hits_data["gz"], nan=0.0)
    hpoints = np.column_stack([hx, hy, hz])
    hit_polydata = pv.PolyData(hpoints)

    for name in hits_data.keys():
        if name not in ("gx", "gy", "gz"):
            hit_polydata.point_data[name] = hits_data[name]

    hit_polydata.point_data["data_type"] = np.ones(hit_polydata.n_points, dtype=int)

    hit_polydata.point_data["event"] = event_hits

    unique_hit_events = np.unique(event_hits)
    print(f"Unique hit event IDs in {filename}: {unique_hit_events}")

    return cluster_polydata, hit_polydata
