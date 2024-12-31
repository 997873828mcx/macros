# analysis.py

import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TDirectory
from array import array
import os
import logging

# from helix_fitting import (
#    find_helix_points_at_radius_analytic,
# )  # Ensure helix_fitting.py is in the same directory
from scipy.spatial import cKDTree


# Configure logging
logging.basicConfig(
    filename="analysis.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


track_id_var = array("i", [0])  # Integer
sigma_rphi_var = array("d", [0.0])  # Double
sigma_z_var = array("d", [0.0])  # Double

delta_rphi_vec = ROOT.std.vector("double")()
delta_z_vec = ROOT.std.vector("double")()

tree_initialized = False
tree = None
file = None


def find_helix_points_at_radius_analytic(
    r, helix_params, ref_theta, theta_range=np.pi / 2
):
    """
    Find ALL helix points at a given radius within theta_range of ref_theta.
    Uses analytic solution instead of numerical solver.

    Parameters:
    -----------
    r : float
        Radius where we want to find points (distance from TPC origin)
    helix_params : dict
        Helix parameters including c_x, c_y, r (helix radius), etc.
    ref_theta : float
        Reference theta value (usually from the second point used in fitting)
    theta_range : float
        How far from ref_theta to look for solutions

    Returns:
    --------
    list of tuples (phi, z, theta) for all valid solutions
    """
    # Extract parameters
    cx = helix_params["c_x"]
    cy = helix_params["c_y"]
    cz = helix_params["c_z"]
    R = helix_params["r"]  # helix radius
    alpha = helix_params["alpha"]

    # The equation r² = (cx + R*cos(θ))² + (cy + R*sin(θ))² can be rewritten as:
    # A*cos(θ) + B*sin(θ) = C where:
    A = 2 * R * cx
    B = 2 * R * cy
    C = r**2 - (cx**2 + cy**2 + R**2)

    # This equation has the form: a*cos(θ) + b*sin(θ) = c
    # Solution is: θ = arctan2(b,a) ± arccos(c/sqrt(a² + b²))

    # First check if solution exists
    norm = np.sqrt(A**2 + B**2)
    if abs(C / norm) > 1:
        logging.debug(f"No solutions exist for r={r}, C/norm={C/norm}")
        return []  # No solutions exist

    # Get the two solutions
    beta = np.arctan2(B, A)
    gamma = np.arccos(C / norm)

    theta1 = beta + gamma
    theta2 = beta - gamma

    # Check which solutions are within range of ref_theta
    solutions = []
    for theta in [theta1, theta2]:
        # Normalize theta to be near ref_theta
        while theta < ref_theta - np.pi:
            theta += 2 * np.pi
        while theta > ref_theta + np.pi:
            theta -= 2 * np.pi

        if abs(theta - ref_theta) <= theta_range:
            # Calculate helix point coordinates
            x = cx + R * np.cos(theta)
            y = cy + R * np.sin(theta)
            z = cz + alpha * theta

            # Calculate phi relative to TPC center
            phi = np.arctan2(y, x)

            solutions.append((phi, z, theta))
    logging.debug(f"Found {len(solutions)} solutions for r={r}")
    return solutions


def apply_helix_filter(point, helix_params, rphi_window, z_window):
    """
    Filter a point based on rphi and z windows.
    Returns True if point is within windows, False otherwise.
    """
    x, y, z = point
    r_meas = np.sqrt(x**2 + y**2)
    phi_meas = np.arctan2(y, x)  # phi from TPC center

    # Get all solutions at this radius within ±pi/2 of reference theta
    solutions = find_helix_points_at_radius_analytic(
        r_meas, helix_params, helix_params["ref_theta"], theta_range=np.pi / 2
    )

    if not solutions:
        return False

    # Find solution with phi closest to measured phi
    min_dphi = float("inf")
    best_solution = None

    for solution in solutions:
        phi_sol, z_sol, _ = solution
        dphi = abs(phi_meas - phi_sol)
        dphi = min(dphi, 2 * np.pi - dphi)  # Take smaller angle

        if dphi < min_dphi:
            min_dphi = dphi
            best_solution = solution

    if best_solution is None:
        return False

    phi_sol, z_sol, _ = best_solution

    # Check if point is within windows
    arc_length = r_meas * min_dphi
    dz = abs(z - z_sol)

    within_windows = arc_length <= rphi_window and dz <= z_window
    logging.debug(
        f"Cluster ({x}, {y}, {z}): arc_length={arc_length:.3f}, dz={dz:.3f}, within_windows={within_windows}"
    )

    return within_windows


def calculate_deltas(helix_params, clusters):
    """
    Calculate delta rphi and delta z for a given set of clusters relative to the helix.

    Parameters:
    -----------
    helix_params : dict
        Helix parameters as returned by helix_fitting.py functions.
    clusters : np.ndarray
        Nx3 array of cluster coordinates [x, y, z].

    Returns:
    --------
    delta_rphi : np.ndarray
        Array of delta rphi values.
    delta_z : np.ndarray
        Array of delta z values.
    """

    ref_theta = helix_params.get("ref_theta", 0.0)

    delta_rphi = []
    delta_z = []

    for point in clusters:
        x, y, z = point

        r_meas = np.sqrt(x**2 + y**2)
        phi_meas = np.arctan2(y, x)

        # Find reference points on the helix at this radius
        solutions = find_helix_points_at_radius_analytic(
            r_meas, helix_params, ref_theta, theta_range=np.pi / 2
        )

        if not solutions:
            logging.debug(f"No helix points found for cluster at ({x}, {y}, {z})")
            continue  # No valid reference points found

        # Find the solution with phi closest to measured phi
        min_dphi = float("inf")
        best_solution = None

        for solution in solutions:
            phi_sol, z_sol, theta_sol = solution
            dphi = abs(phi_meas - phi_sol)
            dphi = min(dphi, 2 * np.pi - dphi)  # Normalize to [0, pi]

            if dphi < min_dphi:
                min_dphi = dphi
                best_solution = solution

        if best_solution is None:
            logging.debug(f"No best solution found for cluster at ({x}, {y}, {z})")
            continue  # No valid solution found

        phi_sol, z_sol, theta_sol = best_solution

        # Calculate delta rphi and delta z
        delta_phi = phi_meas - phi_sol
        # Normalize delta_phi to [-pi, pi]
        delta_phi = (delta_phi + np.pi) % (2 * np.pi) - np.pi

        delta_z_val = z - z_sol
        delta_rphi_val = r_meas * delta_phi
        delta_rphi.append(delta_rphi_val)
        delta_z.append(delta_z_val)

    return np.array(delta_rphi), np.array(delta_z)


def save_histograms(root_filename, track_id, delta_rphi, delta_z):
    """
    Save delta rphi and delta z histograms for a track into a ROOT file.

    Parameters:
    -----------
    root_filename : str
        Path to the ROOT file.
    track_id : int
        Unique identifier for the track.
    delta_rphi : np.ndarray
        Array of delta rphi values.
    delta_z : np.ndarray
        Array of delta z values.
    """
    # Open the ROOT file in update mode or create it if it doesn't exist
    file = TFile.Open(root_filename, "UPDATE")
    if not file or file.IsZombie():
        logging.error(f"Cannot open or create ROOT file: {root_filename}")
        raise IOError(f"Cannot open or create ROOT file: {root_filename}")

    # Create a directory for histograms if not exists
    hist_dir = file.Get("Histograms")
    if not hist_dir:
        hist_dir = file.mkdir("Histograms")
    hist_dir.cd()

    # Create histograms
    hist_rphi = TH1D(
        f"delta_rphi_track_{track_id}",
        f"Delta rphi for Track {track_id};Delta rphi (rad);Counts",
        100,
        -1,
        1,
    )
    for drphi in delta_rphi:
        hist_rphi.Fill(drphi)

    hist_z = TH1D(
        f"delta_z_track_{track_id}",
        f"Delta z for Track {track_id};Delta z (cm);Counts",
        100,
        -2,
        2,
    )
    for dz in delta_z:
        hist_z.Fill(dz)

    # Write histograms to file
    hist_rphi.Write()
    hist_z.Write()

    logging.info(f"Saved histograms for Track ID {track_id} to {root_filename}")

    # Close the file
    file.Close()


def save_tree(root_filename, track_id, delta_rphi, delta_z, sigma_rphi, sigma_z):
    """
    Save delta rphi and delta z values along with their sigmas into a TTree.

    Parameters:
    -----------
    root_filename : str
        Path to the ROOT file.
    track_id : int
        Unique identifier for the track.
    delta_rphi : np.ndarray
        Array of delta rphi values.
    delta_z : np.ndarray
        Array of delta z values.
    sigma_rphi : float
        Standard deviation of delta rphi.
    sigma_z : float
        Standard deviation of delta z.
    """
    # Open the ROOT file in update mode or create it if it doesn't exist
    file = TFile.Open(root_filename, "UPDATE")
    if not file or file.IsZombie():
        logging.error(f"Cannot open or create ROOT file: {root_filename}")
        raise IOError(f"Cannot open or create ROOT file: {root_filename}")

    # Get or Create a TTree if it doesn't exist
    tree = file.Get("AnalysisTree")
    if not tree:
        tree = TTree("AnalysisTree", "Helix Fitting Analysis")

        track_id_var = array("i", [0])  # Integer
        sigma_rphi_var = array("d", [0.0])  # Double
        sigma_z_var = array("d", [0.0])  # Double

        tree.Branch("track_id", track_id_var, "track_id/I")
        tree.Branch("sigma_rphi", sigma_rphi_var, "sigma_rphi/D")
        tree.Branch("sigma_z", sigma_z_var, "sigma_z/D")

        # Define vector branches without leaf lists
        delta_rphi_vec = ROOT.std.vector("double")()
        delta_z_vec = ROOT.std.vector("double")()

        tree.Branch("delta_rphi", delta_rphi_vec)
        tree.Branch("delta_z", delta_z_vec)

    # Convert numpy arrays to std::vector<double>
    delta_rphi_vec = ROOT.std.vector("double")(delta_rphi.tolist())
    delta_z_vec = ROOT.std.vector("double")(delta_z.tolist())

    # Set branch values
    tree.track_id = track_id
    tree.delta_rphi = delta_rphi_vec
    tree.delta_z = delta_z_vec
    tree.sigma_rphi = sigma_rphi
    tree.sigma_z = sigma_z

    # Fill the tree
    tree.Fill()

    # Write the tree to file
    tree.Write("", ROOT.TObject.kOverwrite)

    logging.info(f"Saved TTree entry for Track ID {track_id} to {root_filename}")

    # Close the file
    file.Close()
