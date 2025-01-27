import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from scipy.optimize import least_squares
from analysis import find_helix_points_at_radius_analytic


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


def generate_helix_points_initial(params, inner_cut, outer_cut, num_points=500):
    """
    From initial helix parameters [c_x, c_y, r, alpha, c_z, t0], generate helix points for visualization.
    """
    c_x = params["c_x"]
    c_y = params["c_y"]
    r = params["r"]
    alpha = params["alpha"]
    c_z = params["c_z"]
    ref_theta = params["ref_theta"]

    inner_end = None
    outer_end = None
    inner_solutions = find_helix_points_at_radius_analytic(
        inner_cut, params, ref_theta, np.pi / 2
    )
    outer_solutions = find_helix_points_at_radius_analytic(
        outer_cut, params, ref_theta, np.pi / 2
    )

    if not inner_solutions or not outer_solutions:
        return None

    max_dtheta = 0
    for solution in inner_solutions:
        _, _, theta = solution
        dtheta = abs(theta - ref_theta)
        dtheta = min(dtheta, 2 * np.pi - dtheta)  # Take smaller angle

        if dtheta > max_dtheta:
            max_dtheta = dtheta
            inner_end = solution

    if inner_end is None:
        return None

    min_dtheta = float("inf")
    for solution in outer_solutions:
        _, _, theta = solution
        dtheta = abs(theta - ref_theta)
        dtheta = min(dtheta, 2 * np.pi - dtheta)  # Take smaller angle

        if dtheta < min_dtheta:
            min_dtheta = dtheta
            outer_end = solution

    if outer_end is None:
        return None

    _, _, theta_inner = inner_end
    _, _, theta_outer = outer_end

    t = np.linspace(theta_inner, theta_outer, num_points)

    x = c_x + r * np.cos(t)
    y = c_y + r * np.sin(t)
    z = c_z + alpha * t

    return np.column_stack((x, y, z))


def generate_helix_points_refined(params, inner_cut, outer_cut, num_points=500):
    """
    From refined helix parameters [c_x, c_y, c_z, r, phi, alpha], generate helix points for visualization.
    """
    c_x = params["c_x"]
    c_y = params["c_y"]
    c_z = params["c_z"]
    r = params["r"]
    ref_theta_direct = params["ref_theta_direct"]
    alpha = params["alpha"]

    inner_end = None
    outer_end = None

    inner_solutions = find_helix_points_at_radius_analytic(
        inner_cut, params, ref_theta_direct, np.pi / 2
    )
    outer_solutions = find_helix_points_at_radius_analytic(
        outer_cut, params, ref_theta_direct, np.pi / 2
    )

    if not inner_solutions or not outer_solutions:
        return None

    max_dtheta = 0
    for solution in inner_solutions:
        _, _, theta = solution
        dtheta = abs(theta - ref_theta_direct)
        dtheta = min(dtheta, 2 * np.pi - dtheta)  # Take smaller angle

        if dtheta > max_dtheta:
            max_dtheta = dtheta
            inner_end = solution

    if inner_end is None:
        return None

    min_dtheta = float("inf")
    for solution in outer_solutions:
        _, _, theta = solution
        dtheta = abs(theta - ref_theta_direct)
        dtheta = min(dtheta, 2 * np.pi - dtheta)  # Take smaller angle

        if dtheta < min_dtheta:
            min_dtheta = dtheta
            outer_end = solution

    if outer_end is None:
        return None

    _, _, theta_inner = inner_end
    _, _, theta_outer = outer_end

    t = np.linspace(theta_inner, theta_outer, num_points)

    x = c_x + r * np.cos(t)
    y = c_y + r * np.sin(t)
    z = c_z + alpha * t

    return np.column_stack((x, y, z))


def generate_helix_line(helix_points):
    """Create a polyline visualization of the helix"""
    helix_poly = pv.PolyData(helix_points)
    lines = np.hstack(
        ([helix_points.shape[0]], np.arange(helix_points.shape[0]))
    ).astype(np.int64)
    helix_poly.lines = lines
    return helix_poly


def fit_helix_direct(points, initial_params):
    """
    Perform direct helix fitting using nonlinear optimization.
    Parameters:
        - points: Nx3 array of points to fit.
        - initial_params: dict with initial helix parameters.
    Returns:
        - dict with refined helix parameters or None on failure.
    """

    points = np.array(points)
    # Extract initial helix parameters
    c_x0 = initial_params["c_x"]
    c_y0 = initial_params["c_y"]
    c_z0 = initial_params["c_z"]
    r0 = initial_params["r"]
    alpha0 = initial_params["alpha"]
    # t0 is not used in direct fitting

    # Initial guess for global parameters
    # phi0 = 0.0  # Initial phase offset
    initial_guess = np.array(
        [
            c_x0,  # c_x
            c_y0,  # c_y
            c_z0,  # c_z
            r0,  # r
            # phi0,  # phi
            alpha0,  # alpha
        ]
    )

    # theta0=np.arctan2(points[1, 1] - c_y0, points[1, 0] - c_x0)
    # Define residuals for least squares
    def residuals(params, points):
        c_x, c_y, c_z, r, alpha = params
        # Calculate theta for each point based on current helix parameters
        theta = np.arctan2(points[:, 1] - c_y, points[:, 0] - c_x)
        theta = np.unwrap(theta)  # Ensure continuity
        # Calculate theta for each point based on current helix parameters
        # theta = np.arctan2(points[:, 1] - c_y, points[:, 0] - c_x)
        # Calculate fitted positions
        x_fit = c_x + r * np.cos(theta)
        y_fit = c_y + r * np.sin(theta)
        z_fit = c_z + alpha * theta
        # Compute residuals as the difference between actual and fitted positions
        residuals = points - np.column_stack((x_fit, y_fit, z_fit))
        return residuals.ravel()

    try:
        # Perform least squares optimization
        """
        result = least_squares(
            residuals,
            initial_guess,
            args=(points,),
            method="lm",  # Levenberg-Marquardt algorithm, need to add huber loss function
            max_nfev=1000,
        )
        """

        result = least_squares(
            residuals,
            initial_guess,
            args=(points,),
            method="trf",  # Trust Region Reflective algorithm
            loss="huber",  # Specify Huber loss for robustness
            f_scale=0.3,  # Tuning parameter for Huber loss
            max_nfev=3000,
            verbose=2,  # Enable verbosity for debugging; set to 0 for silent
        )

        if not result.success:
            print("Direct helix fitting did not converge.")
            return None

        fitted = result.x
        c_x, c_y, c_z, r, alpha = fitted[:5]

        # Calculate theta0 as the mean of unwrapped theta
        theta = np.arctan2(points[:, 1] - c_y, points[:, 0] - c_x)
        theta = np.unwrap(theta)  # Ensure continuity
        theta0 = np.mean(theta)

        return {
            "c_x": c_x,
            "c_y": c_y,
            "c_z": c_z,
            "r": r,
            # "phi": phi,
            "alpha": alpha,
            "ref_theta_direct": theta0,
        }

    except Exception as e:
        print(f"Error during direct helix fitting: {e}")
        return None
