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

    def locate_tree(tree_names):
        for name in tree_names:
            if name in file:
                return file[name], name
        return None, None

    intersections_tree, intersections_name = locate_tree(
        ("truth_intersections", "laser_intersections", "combined_clusters")
    )
    hits_tree, hits_name = locate_tree(
        ("truth_g4hits", "laser_g4hits", "combined_hits")
    )

    if intersections_tree is None:
        raise KeyError(
            "The ROOT file does not contain a 'laser_intersections' (or legacy 'combined_clusters') tree."
        )
    if hits_tree is None:
        raise KeyError(
            "The ROOT file does not contain a 'laser_g4hits' (or legacy 'combined_hits') tree."
        )

    try:
        intersections_data = intersections_tree.arrays(library="np")
    except Exception as e:
        raise ValueError(
            f"Failed to extract arrays from '{intersections_name}' tree: {e}"
        )

    required_intersection_branches = {"gx", "gy", "gz", "side"}
    if not required_intersection_branches.issubset(intersections_data.keys()):
        missing = required_intersection_branches - set(intersections_data.keys())
        raise KeyError(
            f"Missing branches in '{intersections_name}' tree: {missing}"
        )

    if "event" in intersections_data.keys():
        event_intersections = intersections_data["event"]
        if np.issubdtype(event_intersections.dtype, np.number):
            event_intersections = event_intersections.astype(int)
        else:
            event_intersections = np.array(
                [
                    (
                        int(e.split()[1])
                        if isinstance(e, str) and len(e.split()) > 1
                        else 0
                    )
                    for e in event_intersections
                ]
            )
    else:
        print(
            f"Warning: 'event' branch not found in '{intersections_name}' tree. Assigning event=0 to all entries."
        )
        event_intersections = np.zeros(
            intersections_data["gx"].shape, dtype=int
        )

    ix = np.nan_to_num(intersections_data["gx"], nan=0.0)
    iy = np.nan_to_num(intersections_data["gy"], nan=0.0)
    iz = np.nan_to_num(intersections_data["gz"], nan=0.0)
    ipoints = np.column_stack([ix, iy, iz])

    intersections_polydata = pv.PolyData(ipoints)
    for name in intersections_data.keys():
        if name in ("gx", "gy", "gz"):
            continue
        try:
            intersections_polydata.point_data[name] = intersections_data[name]
        except Exception as e:
            print(f"Warning: Could not add branch {name} from '{intersections_name}': {e}")

    intersections_polydata.point_data["data_type"] = np.zeros(
        intersections_polydata.n_points, dtype=int
    )
    intersections_polydata.point_data["event"] = event_intersections

    try:
        hits_data = hits_tree.arrays(library="np")
    except Exception as e:
        raise ValueError(f"Failed to extract arrays from '{hits_name}' tree: {e}")

    required_hits_branches = {"gx", "gy", "gz", "side"}
    if not required_hits_branches.issubset(hits_data.keys()):
        missing = required_hits_branches - set(hits_data.keys())
        raise KeyError(f"Missing branches in '{hits_name}' tree: {missing}")

    if "event" in hits_data.keys():
        event_hits = hits_data["event"]
        if np.issubdtype(event_hits.dtype, np.number):
            event_hits = event_hits.astype(int)
        else:
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
    else:
        print(
            f"Warning: 'event' branch not found in '{hits_name}' tree. Assigning event=0 to all entries."
        )
        event_hits = np.zeros(hits_data["gx"].shape, dtype=int)

    hx = np.nan_to_num(hits_data["gx"], nan=0.0)
    hy = np.nan_to_num(hits_data["gy"], nan=0.0)
    hz = np.nan_to_num(hits_data["gz"], nan=0.0)
    hpoints = np.column_stack([hx, hy, hz])
    hits_polydata = pv.PolyData(hpoints)

    for name in hits_data.keys():
        if name in ("gx", "gy", "gz"):
            continue
        hits_polydata.point_data[name] = hits_data[name]

    hits_polydata.point_data["data_type"] = np.ones(
        hits_polydata.n_points, dtype=int
    )
    hits_polydata.point_data["event"] = event_hits

    return intersections_polydata, hits_polydata
