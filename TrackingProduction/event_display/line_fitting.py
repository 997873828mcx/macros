import numpy as np
import pyvista as pv
from scipy.optimize import least_squares


def fit_line_initial(points):
    """
    Fit a straight line through two 3D points.

    Parameters
    ----------
    points : array-like of shape (2, 3)
        An array containing exactly two distinct 3D points.

    Returns
    -------
    dict or None
        Dictionary with line parameters:
            {
              "x0": float,  # x-coordinate of the midpoint
              "y0": float,  # y-coordinate of the midpoint
              "z0": float,  # z-coordinate of the midpoint
              "dir_x": float,  # x-component of the normalized direction vector
              "dir_y": float,  # y-component of the normalized direction vector
              "dir_z": float,  # z-component of the normalized direction vector
            }
        or None if the points are coincident or too close.
    """
    points = np.asarray(points, dtype=float)
    if points.shape != (2, 3):
        print("fit_line_two_points: Exactly two points of shape (3,) are required.")
        return None

    p1 = points[0]
    p2 = points[1]
    
    # Compute the midpoint of the two points
    centroid = (p1 + p2) / 2.0
    
    # Compute the direction vector from p1 to p2
    direction = p2 - p1
    norm = np.linalg.norm(direction)
    if norm < 1e-9:
        print("fit_line_two_points: Points are coincident or too close.")
        return None
    direction /= norm  # Normalize the direction
    
    return {
        "x0": centroid[0],
        "y0": centroid[1],
        "z0": centroid[2],
        "dir_x": direction[0],
        "dir_y": direction[1],
        "dir_z": direction[2],
    }



def generate_line_points(params, length=300.0, num_points=200):
    """
    Generate line points from the initial line parameters for visualization.

    Parameters
    ----------
    params : dict
        {
          "x0": float,
          "y0": float,
          "z0": float,
          "dir_x": float,
          "dir_y": float,
          "dir_z": float,
        }
    length : float
        Half-length (in both directions) of the generated line segment.
    num_points : int
        Number of points to generate.

    Returns
    -------
    numpy.ndarray of shape (num_points, 3)
        Array of 3D coordinates describing the line segment.
    """
    if params is None:
        return None

    x0 = params["x0"]
    y0 = params["y0"]
    z0 = params["z0"]
    dx = params["dir_x"]
    dy = params["dir_y"]
    dz = params["dir_z"]

    # Parameter t goes from -length to +length
    t_vals = np.linspace(-length, length, num_points)
    # Generate points
    x = x0 + dx * t_vals
    y = y0 + dy * t_vals
    z = z0 + dz * t_vals

    return np.column_stack((x, y, z))


def generate_line_polydata(line_points):
    """
    Create a polyline visualization of the line (similar to generate_helix_line).

    Parameters
    ----------
    line_points : numpy.ndarray
        Nx3 array of points along the line.

    Returns
    -------
    pv.PolyData
        A polydata containing a single polyline.
    """
    if line_points is None or len(line_points) < 2:
        return None

    line_poly = pv.PolyData(line_points)
    # Construct lines array: [num_points, 0,1,2,...,num_points-1]
    lines = np.hstack(([line_points.shape[0]], np.arange(line_points.shape[0])))
    line_poly.lines = lines.astype(np.int64)
    return line_poly


def fit_line_direct(points, initial_params):
    """
    Perform a refined straight-line fitting using nonlinear least squares.
    Similar concept to 'fit_helix_direct', but for a 3D line.

    Parameters
    ----------
    points : array-like of shape (N, 3)
        3D points to be fit (filtered clusters).
    initial_params : dict
        The initial line parameters from 'fit_line_initial':
            {
              "x0": float,  # centroid x
              "y0": float,
              "z0": float,
              "dir_x": float,
              "dir_y": float,
              "dir_z": float,
            }

    Returns
    -------
    dict or None
        Refined line parameters in the same dict format if successful,
        or None if fitting fails.
    """
    points = np.asarray(points)
    if points.shape[0] < 2:
        print("fit_line_direct: Not enough points to refine line fitting.")
        return None

    # Extract initial guess
    x0 = initial_params["x0"]
    y0 = initial_params["y0"]
    z0 = initial_params["z0"]
    dir_x = initial_params["dir_x"]
    dir_y = initial_params["dir_y"]
    dir_z = initial_params["dir_z"]

    # Normalize direction in case it's not normalized
    norm_dir = np.sqrt(dir_x**2 + dir_y**2 + dir_z**2)
    if norm_dir < 1e-9:
        print("fit_line_direct: Initial direction is invalid.")
        return None
    dir_x /= norm_dir
    dir_y /= norm_dir
    dir_z /= norm_dir

    initial_guess = np.array([x0, y0, z0, dir_x, dir_y, dir_z])

    def residuals(params, data):
        x0_, y0_, z0_, dx_, dy_, dz_ = params

        # Normalize direction on the fly to avoid blow-ups
        d_norm = np.sqrt(dx_**2 + dy_**2 + dz_**2)
        if d_norm < 1e-12:
            d_norm = 1e-12  # avoid division by zero
        dx_ /= d_norm
        dy_ /= d_norm
        dz_ /= d_norm

        # For each point, we want the orthogonal distance to the line to be minimized.
        # The vector from the line origin to the point is:
        # v_i = (x_i - x0_, y_i - y0_, z_i - z0_)
        # The projection of v_i onto direction is (v_i . d) d
        # The residual is the difference between v_i and that projection (the orth vector).

        diffs = []
        for x_i, y_i, z_i in data:
            vx = x_i - x0_
            vy = y_i - y0_
            vz = z_i - z0_
            dot = vx * dx_ + vy * dy_ + vz * dz_
            # projected point on the line
            px = x0_ + dot * dx_
            py = y0_ + dot * dy_
            pz = z0_ + dot * dz_
            # residual vector from projected point to actual point
            rx = x_i - px
            ry = y_i - py
            rz = z_i - pz
            diffs.extend([rx, ry, rz])
        return np.array(diffs)

    try:
        result = least_squares(
            residuals,
            initial_guess,
            args=(points,),
            method="trf",
            loss="huber",  # robust to outliers
            f_scale=0.1,
            max_nfev=1000,
            verbose=2,  # set to 0 to silence
        )
        if not result.success:
            print("fit_line_direct: Nonlinear fit did not converge.")
            return None

        x0_, y0_, z0_, dx_, dy_, dz_ = result.x
        # Normalize direction once more
        d_norm = np.sqrt(dx_**2 + dy_**2 + dz_**2)
        if d_norm < 1e-12:
            print("fit_line_direct: Final direction collapsed to zero.")
            return None
        dx_ /= d_norm
        dy_ /= d_norm
        dz_ /= d_norm

        return {
            "x0": x0_,
            "y0": y0_,
            "z0": z0_,
            "dir_x": dx_,
            "dir_y": dy_,
            "dir_z": dz_,
        }

    except Exception as e:
        print(f"fit_line_direct: Exception in fitting: {e}")
        return None
