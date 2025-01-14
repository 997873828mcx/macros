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

    # Extract cluster coordinates
    cx = np.nan_to_num(cluster_data["gx"], nan=0.0)
    cy = np.nan_to_num(cluster_data["gy"], nan=0.0)
    cz = np.nan_to_num(cluster_data["gz"], nan=0.0)
    cpoints = np.column_stack([cx, cy, cz])

    cluster_polydata = pv.PolyData(cpoints)
    for name in cluster_data.keys():
        if name not in ("gx", "gy", "gz"):
            cluster_polydata.point_data[name] = cluster_data[name]

    # Add a 'data_type' attribute to distinguish clusters (0)
    cluster_polydata.point_data["data_type"] = np.zeros(
        cluster_polydata.n_points, dtype=int
    )

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

    hx = np.nan_to_num(hits_data["gx"], nan=0.0)
    hy = np.nan_to_num(hits_data["gy"], nan=0.0)
    hz = np.nan_to_num(hits_data["gz"], nan=0.0)
    hpoints = np.column_stack([hx, hy, hz])
    hit_polydata = pv.PolyData(hpoints)

    for name in hits_data.keys():
        if name not in ("gx", "gy", "gz"):
            hit_polydata.point_data[name] = hits_data[name]

    hit_polydata.point_data["data_type"] = np.ones(hit_polydata.n_points, dtype=int)

    return cluster_polydata, hit_polydata
