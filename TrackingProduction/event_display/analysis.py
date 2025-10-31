# analysis.py
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TDirectory
from array import array
import os
import logging
import pyvista as pv

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


def find_line_points_at_radius_analytic(r, line_params):
    """
    Find ALL points on the 3D line that have radius 'r' in the XY-plane.

    Parameters
    ----------
    r : float
        The radius in the XY-plane for which we want to find intersection points.
    line_params : dict
        {
          "x0": float, "y0": float, "z0": float,
          "dir_x": float, "dir_y": float, "dir_z": float
        }

    Returns
    -------
    solutions : list of tuples (phi, z, t)
        Where phi = atan2(y, x), z is the 3D z-coordinate at that point,
        and t is the parameter value on the line.
        The list is empty if no valid intersections are found.
    """
    x0, y0, z0 = line_params["x0"], line_params["y0"], line_params["z0"]
    dx, dy, dz = line_params["dir_x"], line_params["dir_y"], line_params["dir_z"]

    # Quadratic coefficients for the circle radius = r
    A = dx**2 + dy**2
    B = 2.0 * (x0 * dx + y0 * dy)
    C = x0**2 + y0**2 - r**2

    # If the line doesn't vary in the XY-plane, return no solutions
    if abs(A) < 1e-12:
        return []

    # Calculate the discriminant of the quadratic equation
    disc = B**2 - 4 * A * C
    if disc < 0:
        # No real intersection because the line never reaches radius r
        return []

    # Solve for t based on the discriminant
    t_solutions = []
    if abs(disc) < 1e-12:
        # Only one solution (tangent case)
        t_solutions.append(-B / (2 * A))
    else:
        # Two possible solutions for t
        sqrt_disc = np.sqrt(disc)
        t1 = (-B + sqrt_disc) / (2 * A)
        t2 = (-B - sqrt_disc) / (2 * A)
        t_solutions.extend([t1, t2])

    solutions = []
    for t in t_solutions:
        # Compute the 3D point on the line for this t
        X = x0 + dx * t
        Y = y0 + dy * t
        Z = z0 + dz * t

        # Compute the azimuthal angle φ at this point
        phi_line = np.arctan2(Y, X)

        # Collect the solution as (φ, z, t)
        solutions.append((phi_line, Z, t))

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

    ref_theta = helix_params.get("ref_theta_direct", 0.0)

    delta_rphi = []
    delta_z = []
    valid_points = []

    for point in clusters.points:
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
        valid_points.append(point)

    return (
        np.array(delta_rphi),
        np.array(delta_z),
        np.array(valid_points),
    )


def calculate_deltas_line(line_params, clusters):
    """
    Calculate delta rphi and delta z for a set of clusters relative to a fitted line.

    Parameters:
    -----------
    line_params : dict
        Parameters of the fitted line.
    clusters : pv.PolyData
        The set of cluster points.

    Returns:
    --------
    delta_rphi : np.ndarray
        Array of Δrφ values.
    delta_z : np.ndarray
        Array of Δz values.
    """
    delta_rphi = []
    delta_z = []
    valid_points = []

    # Loop over each cluster point
    for point in clusters.points:
        x, y, z = point
        r_meas = np.sqrt(x**2 + y**2)
        phi_meas = np.arctan2(y, x)

        # Find solutions on the line for the radius r_meas
        solutions = find_line_points_at_radius_analytic(r_meas, line_params)
        if not solutions:
            continue  # Skip if no valid line intersections found

        # Select the solution closest in φ to the measured point
        min_dphi = float("inf")
        best_solution = None
        for phi_line, z_line, t_line in solutions:
            dphi = abs(phi_meas - phi_line)
            # Normalize angle difference to [0, π]
            dphi = min(dphi, 2 * np.pi - dphi)
            if dphi < min_dphi:
                min_dphi = dphi
                best_solution = (phi_line, z_line, t_line)

        if best_solution is None:
            continue

        phi_line, z_line, t_line = best_solution

        # Compute differences in φ and z
        delta_phi = phi_meas - phi_line
        # Normalize delta_phi to [-π, π]
        delta_phi = (delta_phi + np.pi) % (2 * np.pi) - np.pi
        # Calculate arc length in the r–φ plane
        delta_rphi_val = r_meas * delta_phi
        delta_z_val = z - z_line

        delta_rphi.append(delta_rphi_val)
        delta_z.append(delta_z_val)
        valid_points.append(point)

    return np.array(delta_rphi), np.array(delta_z), np.array(valid_points)


def find_helix_reference_points(clusters, helix_params, ref_theta):
    """
    Find reference points on the helix for each cluster.
    Returns reference points and their associated deltas.

    Parameters:
    -----------
    clusters : pv.PolyData
        The filtered cluster or hit data as a PyVista PolyData object.
    helix_params : dict
        Dictionary containing helix parameters (e.g., c_x, c_y, r, alpha).
    ref_theta : float
        Reference theta value for helix calculations.

    Returns:
    --------
    reference_points : np.ndarray
        Array of reference point coordinates [x, y, z].
    delta_rphi : np.ndarray
        Array of delta rphi values.
    delta_z : np.ndarray
        Array of delta z values.
    """
    reference_points = []
    delta_rphi = []
    delta_z = []

    # Iterate over the actual points in the PolyData
    for point in clusters.points:
        x, y, z = point
        r_meas = np.sqrt(x**2 + y**2)
        phi_meas = np.arctan2(y, x)

        # Find reference points on the helix at this radius
        solutions = find_helix_points_at_radius_analytic(
            r_meas, helix_params, ref_theta, theta_range=np.pi / 2
        )

        if not solutions:
            continue

        # Find closest solution
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
            continue

        phi_sol, z_sol, theta_sol = best_solution

        # Calculate reference point coordinates
        x_ref = helix_params["c_x"] + helix_params["r"] * np.cos(theta_sol)
        y_ref = helix_params["c_y"] + helix_params["r"] * np.sin(theta_sol)
        z_ref = helix_params["c_z"] + helix_params["alpha"] * theta_sol

        # Calculate deltas
        delta_phi = phi_meas - phi_sol
        delta_phi = (delta_phi + np.pi) % (2 * np.pi) - np.pi  # Normalize to [-pi, pi]

        delta_rphi_val = helix_params["r"] * delta_phi
        delta_z_val = z - z_ref

        reference_points.append([x_ref, y_ref, z_ref])
        delta_rphi.append(delta_rphi_val)
        delta_z.append(delta_z_val)

    return np.array(reference_points), np.array(delta_rphi), np.array(delta_z)


def visualize_reference_points(plotter, clusters, helix_params):
    """
    Add reference points to the visualization.

    Parameters:
    -----------
    plotter : pyvista.Plotter
        The PyVista plotter object for visualization.
    clusters : pv.PolyData
        The filtered cluster or hit data as a PyVista PolyData object.
    helix_params : dict
        Dictionary containing helix parameters.
    """
    # Find reference points
    ref_points, _, _ = find_helix_reference_points(
        clusters, helix_params, helix_params.get("ref_theta", 0.0)
    )

    if len(ref_points) > 0:
        # Create PolyData for reference points
        ref_polydata = pv.PolyData(ref_points)

        # Add reference points to the visualization
        plotter.add_mesh(
            ref_polydata,
            style="points",
            point_size=5,
            color="yellow",  # Different color to distinguish reference points
            pickable=False,
            label="Reference Points",
        )

        # Add lines connecting clusters to their reference points
        for cluster, ref in zip(clusters.points, ref_points):
            line_points = np.array([cluster, ref])
            line = pv.Line(line_points[0], line_points[1])
            plotter.add_mesh(line, color="gray", opacity=0.3)


'''
def calculate_deltas_with_visualization(helix_params, clusters, plotter=None):
    """
    Calculate deltas and optionally visualize reference points.

    Parameters:
    -----------
    helix_params : dict
        Dictionary containing helix parameters.
    clusters : pv.PolyData
        The filtered cluster or hit data as a PyVista PolyData object.
    plotter : pyvista.Plotter, optional
        The PyVista plotter object for visualization.

    Returns:
    --------
    delta_rphi : np.ndarray
        Array of delta rphi values.
    delta_z : np.ndarray
        Array of delta z values.
    """
    # Calculate reference points and deltas
    ref_points, delta_rphi, delta_z = find_helix_reference_points(
        clusters, helix_params, helix_params.get("ref_theta", 0.0)
    )

    if plotter is not None and len(ref_points) > 0:
        # Visualize reference points without recalculating
        visualize_reference_points(plotter, clusters, helix_params)

    return delta_rphi, delta_z

'''


def apply_line_filter(cluster_point, line_params, rphi_window, z_window):
    """
    Check if cluster_point is within rphi_window, z_window of the line, by:
    1) Finding line intersection(s) at radius = r_cluster
    2) Picking the solution with phi closest to cluster phi
    3) Checking arc_length = r_cluster * delta_phi and delta_z
    """
    x, y, z = cluster_point
    r_cluster = np.sqrt(x**2 + y**2)
    phi_meas = np.arctan2(y, x)

    # Find solutions for that radius
    solutions = find_line_points_at_radius_analytic(r_cluster, line_params)

    if not solutions:
        return False

    # Among these, pick the one with phi closest to phi_meas
    best_sol = None
    min_dphi = float("inf")
    for phi_line, z_line, t_line in solutions:
        dphi = abs(phi_meas - phi_line)
        dphi = min(dphi, 2 * np.pi - dphi)
        if dphi < min_dphi:
            min_dphi = dphi
            best_sol = (phi_line, z_line, t_line)

    if not best_sol:
        return False

    phi_line, z_line, t_line = best_sol

    # Now measure arc_length = r_cluster * delta_phi
    delta_phi = phi_meas - phi_line
    delta_phi = (delta_phi + np.pi) % (2 * np.pi) - np.pi
    arc_length = r_cluster * abs(delta_phi)

    dz = abs(z - z_line)

    # Check thresholds
    return (arc_length <= rphi_window) and (dz <= z_window)


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
    Save delta rphi and delta z values along with their sigmas into a TTree in the ROOT file.
    Each call to this function adds a new entry to the TTree.
    """
    delta_rphi_vec = ROOT.std.vector("double")()
    delta_z_vec = ROOT.std.vector("double")()
    # Open the ROOT file in update mode
    file = TFile.Open(root_filename, "UPDATE")
    if not file or file.IsZombie():
        logging.error(f"Cannot open or create ROOT file: {root_filename}")
        raise IOError(f"Cannot open or create ROOT file: {root_filename}")

    tree = file.Get("AnalysisTree")
    if not tree:

        # Create the TTree and define branches

        tree = TTree("AnalysisTree", "Helix Fitting Analysis")

        # Define branches
        tree.Branch("track_id", track_id_var, "track_id/I")
        tree.Branch("sigma_rphi", sigma_rphi_var, "sigma_rphi/D")
        tree.Branch("sigma_z", sigma_z_var, "sigma_z/D")
        tree.Branch("delta_rphi", delta_rphi_vec)
        tree.Branch("delta_z", delta_z_vec)

    else:
        tree.SetBranchAddress("track_id", track_id_var)
        tree.SetBranchAddress("sigma_rphi", sigma_rphi_var)
        tree.SetBranchAddress("sigma_z", sigma_z_var)
        tree.SetBranchAddress("delta_rphi", delta_rphi_vec)
        tree.SetBranchAddress("delta_z", delta_z_vec)

    # Set scalar branch values
    track_id_var[0] = track_id
    sigma_rphi_var[0] = sigma_rphi
    sigma_z_var[0] = sigma_z

    # Clear and fill vector branches
    delta_rphi_vec.clear()
    delta_z_vec.clear()
    for drphi in delta_rphi:
        delta_rphi_vec.push_back(float(drphi))
    for dz in delta_z:
        delta_z_vec.push_back(float(dz))

    # Fill tree and write
    tree.Fill()

    # Write the tree to the file
    tree.Write("", ROOT.TObject.kOverwrite)

    logging.info(f"Saved entry for Track ID {track_id} to '{root_filename}'.")

    # Close the file
    file.Close()
    logging.info(f"Saved entry for Track ID {track_id} to '{root_filename}'.")


def compute_centroid(points):
    """
    Compute the centroid of a set of 3D points.

    Parameters
    ----------
    points : array-like of shape (N, 3)
        Array or list of 3D points.

    Returns
    -------
    np.ndarray of shape (3,)
        The coordinates of the centroid [x_mean, y_mean, z_mean].
    """
    points = np.asarray(points)
    return np.mean(points, axis=0)


def find_pca_two_tracks(a1, b1, a2, b2):

    # Convert inputs to numpy arrays
    a1 = np.asarray(a1, dtype=float)
    b1 = np.asarray(b1, dtype=float)
    a2 = np.asarray(a2, dtype=float)
    b2 = np.asarray(b2, dtype=float)

    # --- 1) Calculate DCA using cross product
    b_cross_b = np.cross(b1, b2)
    mag_b_cross_b = np.linalg.norm(b_cross_b)
    a2_a1 = a2 - a1

    # If cross product is zero, lines are parallel => your snippet returned 999 or skip
    # For directness, replicate the approach:
    if mag_b_cross_b == 0:
        # Lines are parallel or coincide; handle as needed
        dca = 999.0
        pca1 = None
        pca2 = None
        return pca1, pca2, dca

    # The sign of dca is typically made positive by abs
    dca = abs(np.dot(b_cross_b, a2_a1)) / mag_b_cross_b

    # --- 2) Calculate 'c' and 'd' (the scalars) from the snippet
    # Carefully replicate each step:
    #
    # double X = b1.dot(b2) - b1.dot(b1)*b2.dot(b2)/b2.dot(b1);
    # double Y = (a2.dot(b2) - a1.dot(b2)) - (a2.dot(b1) - a1.dot(b1)) * b2.dot(b2)/b2.dot(b1);
    # double c = Y / X;
    #
    # double F = b1.dot(b1)/b2.dot(b1);
    # double G = -( (a2.dot(b1) - a1.dot(b1)) / b2.dot(b1) );
    # double d = c * F + G;

    dot_b1_b2 = np.dot(b1, b2)
    dot_b1_b1 = np.dot(b1, b1)
    dot_b2_b2 = np.dot(b2, b2)
    dot_a2_b2_minus_a1_b2 = np.dot(a2, b2) - np.dot(a1, b2)
    dot_a2_b1_minus_a1_b1 = np.dot(a2, b1) - np.dot(a1, b1)
    dot_b2_b1 = np.dot(b2, b1)  # same as dot_b1_b2 but let's keep naming consistent

    # Check for zero denominators
    # The snippet implicitly assumes b2.dot(b1) != 0
    if dot_b2_b1 == 0:
        # Lines are orthonormal or parallel in some sense => snippet would break
        # Decide how to handle it, for directness we do the same:
        dca = 999.0
        return None, None, dca

    X = dot_b1_b2 - (dot_b1_b1 * dot_b2_b2 / dot_b2_b1)
    Y = (dot_a2_b2_minus_a1_b2) - (dot_a2_b1_minus_a1_b1) * dot_b2_b2 / dot_b2_b1
    c = Y / X

    F = dot_b1_b1 / dot_b2_b1
    G = -(dot_a2_b1_minus_a1_b1 / dot_b2_b1)
    d = c * F + G

    # --- 3) Points of closest approach
    pca1 = a1 + c * b1
    pca2 = a2 + d * b2

    return pca1, pca2, dca


def calculate_dca_to_beam(line_params):
    """
    Calculate DCA and PCA between a fitted line and the beam axis (z-axis).

    Parameters
    ----------
    line_params : dict
        The line parameters from fit_line_initial or fit_line_direct:
        {
            "x0": float,  # point on the line
            "y0": float,
            "z0": float,
            "dir_x": float,  # direction vector
            "dir_y": float,
            "dir_z": float,
        }

    Returns
    -------
    tuple
        (pca_line, pca_beam, dca) where:
        - pca_line is the point of closest approach on the fitted line
        - pca_beam is the point of closest approach on the beam axis
        - dca is the distance of closest approach
    """
    # Point on your fitted line
    a1 = np.array([line_params["x0"], line_params["y0"], line_params["z0"]])

    # Direction vector of your fitted line
    b1 = np.array([line_params["dir_x"], line_params["dir_y"], line_params["dir_z"]])

    # Point on the beam axis (we can use the origin)
    a2 = np.array([0.0, 0.0, 0.0])

    # Direction vector of beam axis (unit vector along z)
    b2 = np.array([0.0, 0.0, 1.0])

    # Calculate PCA and DCA
    pca_line, pca_beam, dca = find_pca_two_tracks(a1, b1, a2, b2)

    return pca_line, pca_beam, dca


def format_dca_info(line_params):
    """
    Format DCA information for a fitted line into a string.

    Parameters
    ----------
    line_params : dict
        The line parameters from fit_line_initial or fit_line_direct

    Returns
    -------
    str
        Formatted string containing DCA and PCA information
    """
    pca_line, pca_beam, dca = calculate_dca_to_beam(line_params)

    info = [f"Distance of Closest Approach (DCA): {dca:.3f} cm"]

    if pca_line is not None:
        info.extend(
            [
                "\nPoint of Closest Approach on fitted line:",
                f"x: {pca_line[0]:.3f} cm",
                f"y: {pca_line[1]:.3f} cm",
                f"z: {pca_line[2]:.3f} cm",
                "\nPoint of Closest Approach on beam axis:",
                f"x: {pca_beam[0]:.3f} cm",
                f"y: {pca_beam[1]:.3f} cm",
                f"z: {pca_beam[2]:.3f} cm",
            ]
        )
    else:
        info.append("Note: PCA points could not be calculated (lines may be parallel)")

    return "\n".join(info)
