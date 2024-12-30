import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from scipy.optimize import least_squares


def fit_circle_2d(p1, p2, p3):
    """
    Fit a circle to three points (2D). Return (center, radius) or (None, None).
    """
    A = 2 * (p2[0] - p1[0])
    B = 2 * (p2[1] - p1[1])
    C = p2[0] ** 2 + p2[1] ** 2 - p1[0] ** 2 - p1[1] ** 2
    D = 2 * (p3[0] - p1[0])
    E = 2 * (p3[1] - p1[1])
    F = p3[0] ** 2 + p3[1] ** 2 - p1[0] ** 2 - p1[1] ** 2

    det = A * E - B * D
    if abs(det) < 1e-9:
        return None, None  # Collinear or insufficient geometry

    cx = (C * E - B * F) / det
    cy = (A * F - C * D) / det
    r = np.sqrt((p1[0] - cx) ** 2 + (p1[1] - cy) ** 2)
    return np.array([cx, cy]), r


def fit_helix_initial(points):
    """
    Fit a helix to exactly three points (initial fitting).
    Returns a dictionary of helix parameters or None on failure.
    """
    if len(points) != 3:
        return None

    p1, p2, p3 = points

    # --- 1. Project onto XY-plane ---
    p1_xy = p1[:2]
    p2_xy = p2[:2]
    p3_xy = p3[:2]

    # --- 2. Fit circle in XY-plane to find center and radius ---
    center_xy, radius = fit_circle_2d(p1_xy, p2_xy, p3_xy)
    if center_xy is None or radius <= 0:
        return None  # Could not fit a circle (collinear or invalid)

    # --- 3. Calculate angular positions (theta) in XY-plane ---
    def calc_theta(p, center):
        return np.arctan2(p[1] - center[1], p[0] - center[0])

    theta1 = calc_theta(p1_xy, center_xy)
    theta2 = calc_theta(p2_xy, center_xy)
    theta3 = calc_theta(p3_xy, center_xy)
    thetas = np.array([theta1, theta2, theta3])
    thetas = np.unwrap(thetas)  # ensure continuous
    ref_theta = thetas[1]

    # --- 4. Check if points span multiple turns ---
    if thetas[-1] - thetas[0] > 2 * np.pi:
        raise ValueError(
            "Selected points span multiple helical turns. Please select three points within the same turn."
        )

    # --- 5. Linear fit between theta and z to find pitch ---
    z_vals = np.array([p1[2], p2[2], p3[2]])
    A = np.vstack([thetas, np.ones_like(thetas)]).T
    try:
        slope, z0 = np.linalg.lstsq(A, z_vals, rcond=None)[0]
    except Exception:
        return None

    # pitch = slope * (2 * np.pi)  # pitch for a 2pi rotation
    t0 = thetas[0]  # phase offset
    center_z = z0

    return {
        "c_x": center_xy[0],
        "c_y": center_xy[1],
        "r": radius,
        "alpha": slope,  # pitch per radian
        "c_z": center_z,
        "t0": t0,
        "ref_theta": ref_theta,
    }


def generate_helix_points_initial(params, num_points=500):
    """
    From initial helix parameters [c_x, c_y, r, alpha, c_z, t0], generate helix points for visualization.
    """
    c_x = params["c_x"]
    c_y = params["c_y"]
    r = params["r"]
    alpha = params["alpha"]
    c_z = params["c_z"]
    t0 = params["ref_theta"]

    # Generate a range of theta values around t0 for visualization
    theta_min = t0 - np.pi / 2
    theta_max = t0 + np.pi / 2  # Adjust as needed for visualization
    t = np.linspace(theta_min, theta_max, num_points)

    x = c_x + r * np.cos(t)
    y = c_y + r * np.sin(t)
    z = c_z + alpha * t

    return np.column_stack((x, y, z))


def generate_helix_points_refined(params, num_points=500):
    """
    From refined helix parameters [c_x, c_y, c_z, r, phi, alpha], generate helix points for visualization.
    """
    c_x = params["c_x"]
    c_y = params["c_y"]
    c_z = params["c_z"]
    r = params["r"]
    phi = params["phi"]
    alpha = params["alpha"]

    # Generate a range of theta values around t0 for visualization
    theta_min = -np.pi
    theta_max = np.pi  # Adjust as needed for visualization
    t = np.linspace(theta_min, theta_max, num_points)
    # Generate a range of theta values for visualization
    # t = np.linspace(0, 4 * np.pi, num_points)  # 2 full turns

    x = c_x + r * np.cos(t + phi)
    y = c_y + r * np.sin(t + phi)
    z = c_z + alpha * t

    return np.column_stack((x, y, z))


def create_helix_tube(helix_points, tube_radius=0.5):
    """
    Create a tubular mesh around the helix points for visualization.
    """
    if helix_points.shape[0] < 2:
        return None

    helix_poly = pv.PolyData(helix_points)
    # Create a single polyline connecting all points
    lines = np.hstack(
        ([helix_points.shape[0]], np.arange(helix_points.shape[0]))
    ).astype(np.int64)
    helix_poly.lines = lines
    try:
        tube = helix_poly.tube(radius=tube_radius)
        return tube
    except Exception as e:
        print(f"Error creating helix tube: {e}")
        return None


def fit_helix_direct(points, initial_params):
    """
    Perform direct helix fitting using nonlinear optimization.
    Parameters:
        - points: Nx3 array of points to fit.
        - initial_params: dict with initial helix parameters.
    Returns:
        - dict with refined helix parameters or None on failure.
    """
    # Extract initial helix parameters
    c_x0 = initial_params["c_x"]
    c_y0 = initial_params["c_y"]
    c_z0 = initial_params["c_z"]
    r0 = initial_params["r"]
    alpha0 = initial_params["alpha"]
    # t0 is not used in direct fitting

    # Initial guess for global parameters
    phi0 = 0.0  # Initial phase offset
    initial_guess = np.array(
        [
            c_x0,  # c_x
            c_y0,  # c_y
            c_z0,  # c_z
            r0,  # r
            phi0,  # phi
            alpha0,  # alpha
        ]
    )

    # Define residuals for least squares
    def residuals(params, points):
        c_x, c_y, c_z, r, phi, alpha = params
        # Calculate theta for each point based on current helix parameters
        theta = np.arctan2(points[:, 1] - c_y, points[:, 0] - c_x) - phi
        # Calculate fitted positions
        x_fit = c_x + r * np.cos(theta + phi)
        y_fit = c_y + r * np.sin(theta + phi)
        z_fit = c_z + alpha * theta
        # Compute residuals as the difference between actual and fitted positions
        residuals = points - np.column_stack((x_fit, y_fit, z_fit))
        return residuals.ravel()

    try:
        # Perform least squares optimization
        result = least_squares(
            residuals,
            initial_guess,
            args=(points,),
            method="lm",  # Levenberg-Marquardt algorithm, need to add huber loss function
            max_nfev=1000,
        )
        if not result.success:
            print("Direct helix fitting did not converge.")
            return None

        fitted = result.x
        c_x, c_y, c_z, r, phi, alpha = fitted[:6]

        return {
            "c_x": c_x,
            "c_y": c_y,
            "c_z": c_z,
            "r": r,
            "phi": phi,
            "alpha": alpha,
        }

    except Exception as e:
        print(f"Error during direct helix fitting: {e}")
        return None


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

    return arc_length <= rphi_window and dz <= z_window
