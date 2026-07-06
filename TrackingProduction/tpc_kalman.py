#!/usr/bin/env python3
"""Small TPC-only Kalman fitter for ordered space points.

This is intentionally lightweight and standalone.  It assumes a uniform
solenoidal magnetic field and that pattern recognition has already grouped the
TPC points belonging to one track.  The state is

    [x, y, z, phi, signed_q_over_pT, tan(lambda)]

with x/y/z in cm, pT in GeV/c, and B in Tesla.  The independent variable is the
signed transverse path length in cm.  This first version is meant for fast
truth-point studies before porting the same idea into a C++ Fun4All helper.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


STATE_DIM = 6
IDX_X = 0
IDX_Y = 1
IDX_Z = 2
IDX_PHI = 3
IDX_QOP_T = 4
IDX_TANL = 5


@dataclass(frozen=True)
class KalmanConfig:
    bfield_t: float = 1.4
    point_sorting: str = "radius"
    meas_sigma_xy_cm: float = 0.03
    meas_sigma_z_cm: float = 0.05
    min_measurement_sigma_cm: float = 1.0e-6
    initial_sigma_pos_cm: float = 0.1
    initial_sigma_phi: float = 0.2
    initial_sigma_qop_t: float = 0.2
    initial_sigma_tanl: float = 0.2
    process_sigma_pos_cm: float = 1.0e-4
    process_sigma_phi: float = 1.0e-5
    process_sigma_qop_t: float = 1.0e-6
    process_sigma_tanl: float = 1.0e-6
    material_x0_per_cm: float = 0.0
    multiple_scattering_scale: float = 1.0
    energy_loss_gev_per_cm: float = 0.0
    energy_loss_sigma_fraction: float = 0.0
    min_pt_gev: float = 0.05


@dataclass
class SeedHelix:
    cx: float
    cy: float
    radius: float
    z0: float
    pitch: float
    direction: float
    theta: np.ndarray
    path_s: np.ndarray


@dataclass
class TpcKalmanFit:
    success: bool
    message: str
    charge: int
    bfield_t: float
    seed: SeedHelix
    point_order: np.ndarray
    path_s: np.ndarray
    states_filtered: np.ndarray
    covs_filtered: np.ndarray
    states_smoothed: np.ndarray
    covs_smoothed: np.ndarray
    states_predicted: np.ndarray
    covs_predicted: np.ndarray
    transport: np.ndarray
    chi2: float
    ndof: int
    config: KalmanConfig
    mass_gev: float


def pdg_mass_gev(pid: int) -> float:
    mass_by_abs_pid = {
        11: 0.00051099895,
        13: 0.1056583755,
        211: 0.13957039,
        321: 0.493677,
        2212: 0.93827208816,
    }
    return mass_by_abs_pid.get(abs(int(pid)), 0.13957039)


def normalize_phi(phi: float) -> float:
    return float(np.arctan2(np.sin(phi), np.cos(phi)))


def state_residual(lhs: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    diff = np.asarray(lhs, dtype=float) - np.asarray(rhs, dtype=float)
    diff[IDX_PHI] = normalize_phi(diff[IDX_PHI])
    return diff


def fit_circle_least_squares(points: np.ndarray):
    xy = np.asarray(points[:, :2], dtype=float)
    if xy.shape[0] < 3:
        return None

    x = xy[:, 0]
    y = xy[:, 1]
    matrix = np.column_stack((x, y, np.ones_like(x)))
    rhs = -(x * x + y * y)
    try:
        a, b, c = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
    except np.linalg.LinAlgError:
        return None

    cx = -0.5 * a
    cy = -0.5 * b
    radius2 = cx * cx + cy * cy - c
    if radius2 <= 0.0 or not np.isfinite(radius2):
        return None

    radius = float(np.sqrt(radius2))
    if radius <= 0.0 or not np.isfinite(radius):
        return None
    return float(cx), float(cy), radius


def order_points_for_fit(points: np.ndarray, config: KalmanConfig) -> tuple[np.ndarray, np.ndarray] | None:
    """Return points in the order used by the fit and their original indices."""

    mode = str(config.point_sorting).lower()
    if mode == "input":
        order = np.arange(points.shape[0], dtype=int)
    elif mode == "radius":
        radii = np.linalg.norm(points[:, :2], axis=1)
        if not np.all(np.isfinite(radii)):
            return None
        order = np.argsort(radii, kind="mergesort")
    else:
        raise ValueError(f"unknown point_sorting mode: {config.point_sorting}")

    return points[order], order


def make_seed_helix(points: np.ndarray) -> SeedHelix | None:
    circle = fit_circle_least_squares(points)
    if circle is None:
        return None
    cx, cy, radius = circle

    theta = np.unwrap(np.arctan2(points[:, 1] - cy, points[:, 0] - cx))
    if theta.size < 2 or np.nanmax(theta) - np.nanmin(theta) < 1.0e-4:
        return None

    try:
        pitch, z0 = np.linalg.lstsq(
            np.column_stack((theta, np.ones_like(theta))),
            points[:, 2],
            rcond=None,
        )[0]
    except np.linalg.LinAlgError:
        return None

    direction = float(np.sign(theta[-1] - theta[0]))
    if direction == 0.0:
        direction = 1.0

    dtheta = np.diff(theta)
    ds = radius * np.abs(dtheta)
    path_s = np.concatenate(([0.0], np.cumsum(ds)))
    if not np.all(np.isfinite(path_s)):
        return None

    return SeedHelix(
        cx=float(cx),
        cy=float(cy),
        radius=float(radius),
        z0=float(z0),
        pitch=float(pitch),
        direction=direction,
        theta=theta,
        path_s=path_s,
    )


def initial_state_from_seed(points: np.ndarray, seed: SeedHelix, config: KalmanConfig) -> np.ndarray | None:
    theta0 = float(seed.theta[0])
    direction = seed.direction
    phi0 = normalize_phi(np.arctan2(direction * np.cos(theta0), direction * -np.sin(theta0)))
    tanl0 = direction * seed.pitch / seed.radius

    denom = 0.003 * config.bfield_t * seed.radius
    if abs(denom) <= 0.0 or not np.isfinite(denom):
        return None
    signed_qop_t0 = direction / denom

    min_abs_qop_t = 1.0 / max(config.min_pt_gev, 1.0e-6)
    if abs(signed_qop_t0) > min_abs_qop_t:
        signed_qop_t0 = np.sign(signed_qop_t0) * min_abs_qop_t

    return np.array(
        [
            points[0, 0],
            points[0, 1],
            points[0, 2],
            phi0,
            signed_qop_t0,
            tanl0,
        ],
        dtype=float,
    )


def omega_from_state(state: np.ndarray, bfield_t: float) -> float:
    return float(0.003 * bfield_t * state[IDX_QOP_T])


def apply_mean_energy_loss(state: np.ndarray, ds_cm: float, config: KalmanConfig, mass_gev: float) -> np.ndarray:
    if config.energy_loss_gev_per_cm <= 0.0 or ds_cm == 0.0:
        return state

    out = np.asarray(state, dtype=float).copy()
    qop_t = float(out[IDX_QOP_T])
    if abs(qop_t) < 1.0e-12:
        return out

    tanl = float(out[IDX_TANL])
    path3d_cm = abs(ds_cm) * np.sqrt(1.0 + tanl * tanl)
    if path3d_cm <= 0.0:
        return out

    pt = 1.0 / abs(qop_t)
    momentum = pt * np.sqrt(1.0 + tanl * tanl)
    energy = np.sqrt(momentum * momentum + mass_gev * mass_gev)
    signed_loss = np.sign(ds_cm) * config.energy_loss_gev_per_cm * path3d_cm
    new_energy = max(mass_gev + 1.0e-9, energy - signed_loss)
    new_momentum = np.sqrt(max(0.0, new_energy * new_energy - mass_gev * mass_gev))
    if new_momentum <= 0.0:
        return out

    new_pt = new_momentum / np.sqrt(1.0 + tanl * tanl)
    min_pt = max(config.min_pt_gev, 1.0e-6)
    new_pt = max(new_pt, min_pt)
    out[IDX_QOP_T] = np.sign(qop_t) / new_pt
    return out


def propagate_state(
    state: np.ndarray,
    ds_cm: float,
    bfield_t: float,
    config: KalmanConfig | None = None,
    mass_gev: float = 0.13957039,
) -> np.ndarray:
    out = np.asarray(state, dtype=float).copy()
    phi = float(state[IDX_PHI])
    omega = omega_from_state(state, bfield_t)

    if abs(omega) < 1.0e-10:
        out[IDX_X] += ds_cm * np.cos(phi)
        out[IDX_Y] += ds_cm * np.sin(phi)
    else:
        phi2 = phi + omega * ds_cm
        out[IDX_X] += (np.sin(phi2) - np.sin(phi)) / omega
        out[IDX_Y] += -(np.cos(phi2) - np.cos(phi)) / omega
        out[IDX_PHI] = phi2

    out[IDX_Z] += state[IDX_TANL] * ds_cm
    out[IDX_PHI] = normalize_phi(out[IDX_PHI])
    if config is not None:
        out = apply_mean_energy_loss(out, ds_cm, config, mass_gev)
    return out


def transport_jacobian(
    state: np.ndarray,
    ds_cm: float,
    bfield_t: float,
    config: KalmanConfig,
    mass_gev: float,
) -> np.ndarray:
    jac = np.eye(STATE_DIM, dtype=float)
    scales = np.array([1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-6], dtype=float)
    for col in range(STATE_DIM):
        step = scales[col] * max(1.0, abs(float(state[col])))
        plus = state.copy()
        minus = state.copy()
        plus[col] += step
        minus[col] -= step
        if col == IDX_PHI:
            plus[col] = normalize_phi(plus[col])
            minus[col] = normalize_phi(minus[col])
        f_plus = propagate_state(plus, ds_cm, bfield_t, config, mass_gev)
        f_minus = propagate_state(minus, ds_cm, bfield_t, config, mass_gev)
        jac[:, col] = state_residual(f_plus, f_minus) / (2.0 * step)
    return jac


def multiple_scattering_theta0(state: np.ndarray, ds_cm: float, config: KalmanConfig, mass_gev: float) -> float:
    if config.material_x0_per_cm <= 0.0 or config.multiple_scattering_scale <= 0.0:
        return 0.0

    tanl = float(state[IDX_TANL])
    path3d_cm = abs(ds_cm) * np.sqrt(1.0 + tanl * tanl)
    x_over_x0 = path3d_cm * config.material_x0_per_cm
    if x_over_x0 <= 0.0:
        return 0.0

    qop_t = float(state[IDX_QOP_T])
    if abs(qop_t) < 1.0e-12:
        return 0.0
    pt = 1.0 / abs(qop_t)
    momentum = pt * np.sqrt(1.0 + tanl * tanl)
    energy = np.sqrt(momentum * momentum + mass_gev * mass_gev)
    beta = max(momentum / energy, 1.0e-6)
    log_term = 1.0 + 0.038 * np.log(max(x_over_x0, 1.0e-12))
    theta0 = 0.0136 / (beta * momentum) * np.sqrt(x_over_x0) * log_term
    return abs(float(config.multiple_scattering_scale * theta0))


def process_noise(state: np.ndarray, ds_cm: float, config: KalmanConfig, mass_gev: float) -> np.ndarray:
    scale = max(1.0, abs(ds_cm))
    sigmas = np.array(
        [
            config.process_sigma_pos_cm,
            config.process_sigma_pos_cm,
            config.process_sigma_pos_cm,
            config.process_sigma_phi,
            config.process_sigma_qop_t,
            config.process_sigma_tanl,
        ],
        dtype=float,
    )
    q = np.diag((sigmas * scale) ** 2)

    theta0 = multiple_scattering_theta0(state, ds_cm, config, mass_gev)
    if theta0 > 0.0:
        tanl = float(state[IDX_TANL])
        cos_lambda = 1.0 / np.sqrt(1.0 + tanl * tanl)
        cos_lambda = max(cos_lambda, 1.0e-4)
        path3d_cm = abs(ds_cm) * np.sqrt(1.0 + tanl * tanl)
        sigma_phi = theta0 / cos_lambda
        sigma_tanl = theta0 / (cos_lambda * cos_lambda)
        sigma_pos = theta0 * path3d_cm / np.sqrt(3.0)
        q[IDX_X, IDX_X] += sigma_pos * sigma_pos
        q[IDX_Y, IDX_Y] += sigma_pos * sigma_pos
        q[IDX_Z, IDX_Z] += sigma_pos * sigma_pos
        q[IDX_PHI, IDX_PHI] += sigma_phi * sigma_phi
        q[IDX_TANL, IDX_TANL] += sigma_tanl * sigma_tanl

    if config.energy_loss_sigma_fraction > 0.0 and config.energy_loss_gev_per_cm > 0.0:
        tanl = float(state[IDX_TANL])
        path3d_cm = abs(ds_cm) * np.sqrt(1.0 + tanl * tanl)
        sigma_e = config.energy_loss_sigma_fraction * config.energy_loss_gev_per_cm * path3d_cm
        qop_t = float(state[IDX_QOP_T])
        if sigma_e > 0.0 and abs(qop_t) > 1.0e-12:
            pt = 1.0 / abs(qop_t)
            momentum = pt * np.sqrt(1.0 + tanl * tanl)
            energy = np.sqrt(momentum * momentum + mass_gev * mass_gev)
            dp_de = energy / max(momentum, 1.0e-9)
            dpt_dp = 1.0 / np.sqrt(1.0 + tanl * tanl)
            dqop_dpt = 1.0 / (pt * pt)
            sigma_qop = dqop_dpt * dpt_dp * dp_de * sigma_e
            q[IDX_QOP_T, IDX_QOP_T] += sigma_qop * sigma_qop

    return q


def measurement_sigma(value_cm: float, config: KalmanConfig) -> float:
    """Return a strictly positive measurement sigma in cm.

    Exact zero covariance makes the Kalman gain numerically fragile.  For
    perfect truth points, use a tiny but nonzero sigma instead of 0.
    """

    value = abs(float(value_cm))
    floor = abs(float(config.min_measurement_sigma_cm))
    return max(value, floor)


def fit_tpc_kalman(
    points: np.ndarray,
    charge: int,
    config: KalmanConfig | None = None,
    pid: int = 0,
) -> TpcKalmanFit | None:
    config = config or KalmanConfig()
    mass_gev = pdg_mass_gev(pid)
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 5:
        return None

    ordered = order_points_for_fit(points, config)
    if ordered is None:
        return None
    points, point_order = ordered

    seed = make_seed_helix(points)
    if seed is None:
        return None

    state0 = initial_state_from_seed(points, seed, config)
    if state0 is None:
        return None

    npoints = points.shape[0]
    states_f = np.zeros((npoints, STATE_DIM), dtype=float)
    covs_f = np.zeros((npoints, STATE_DIM, STATE_DIM), dtype=float)
    states_p = np.zeros_like(states_f)
    covs_p = np.zeros_like(covs_f)
    transports = np.zeros_like(covs_f)

    p0_sigmas = np.array(
        [
            config.initial_sigma_pos_cm,
            config.initial_sigma_pos_cm,
            config.initial_sigma_pos_cm,
            config.initial_sigma_phi,
            config.initial_sigma_qop_t,
            config.initial_sigma_tanl,
        ],
        dtype=float,
    )
    cov = np.diag(p0_sigmas * p0_sigmas)
    state = state0

    hmat = np.zeros((3, STATE_DIM), dtype=float)
    hmat[0, IDX_X] = 1.0
    hmat[1, IDX_Y] = 1.0
    hmat[2, IDX_Z] = 1.0
    meas_sigma_xy = measurement_sigma(config.meas_sigma_xy_cm, config)
    meas_sigma_z = measurement_sigma(config.meas_sigma_z_cm, config)
    meas_cov = np.diag(
        [
            meas_sigma_xy**2,
            meas_sigma_xy**2,
            meas_sigma_z**2,
        ]
    )

    chi2 = 0.0
    ndof = 0
    eye = np.eye(STATE_DIM)

    for index in range(npoints):
        if index == 0:
            pred_state = state
            pred_cov = cov
            fmat = eye.copy()
        else:
            ds = float(seed.path_s[index] - seed.path_s[index - 1])
            fmat = transport_jacobian(state, ds, config.bfield_t, config, mass_gev)
            pred_state = propagate_state(state, ds, config.bfield_t, config, mass_gev)
            pred_cov = fmat @ cov @ fmat.T + process_noise(state, ds, config, mass_gev)
            pred_cov = 0.5 * (pred_cov + pred_cov.T)

        residual = points[index] - pred_state[:3]
        innovation = hmat @ pred_cov @ hmat.T + meas_cov
        try:
            innovation_inv = np.linalg.inv(innovation)
        except np.linalg.LinAlgError:
            innovation_inv = np.linalg.pinv(innovation)

        gain = pred_cov @ hmat.T @ innovation_inv
        update = gain @ residual
        state = pred_state + update
        state[IDX_PHI] = normalize_phi(state[IDX_PHI])
        cov = (eye - gain @ hmat) @ pred_cov @ (eye - gain @ hmat).T + gain @ meas_cov @ gain.T
        cov = 0.5 * (cov + cov.T)

        chi2 += float(residual.T @ innovation_inv @ residual)
        ndof += 3

        states_p[index] = pred_state
        covs_p[index] = pred_cov
        transports[index] = fmat
        states_f[index] = state
        covs_f[index] = cov

    states_s = states_f.copy()
    covs_s = covs_f.copy()
    for index in range(npoints - 2, -1, -1):
        fmat = transports[index + 1]
        try:
            pred_inv = np.linalg.inv(covs_p[index + 1])
        except np.linalg.LinAlgError:
            pred_inv = np.linalg.pinv(covs_p[index + 1])
        smoother_gain = covs_f[index] @ fmat.T @ pred_inv
        residual = state_residual(states_s[index + 1], states_p[index + 1])
        states_s[index] = states_f[index] + smoother_gain @ residual
        states_s[index, IDX_PHI] = normalize_phi(states_s[index, IDX_PHI])
        covs_s[index] = covs_f[index] + smoother_gain @ (covs_s[index + 1] - covs_p[index + 1]) @ smoother_gain.T
        covs_s[index] = 0.5 * (covs_s[index] + covs_s[index].T)

    ndof -= STATE_DIM
    return TpcKalmanFit(
        success=True,
        message="ok",
        charge=int(charge),
        bfield_t=float(config.bfield_t),
        seed=seed,
        point_order=point_order,
        path_s=seed.path_s,
        states_filtered=states_f,
        covs_filtered=covs_f,
        states_smoothed=states_s,
        covs_smoothed=covs_s,
        states_predicted=states_p,
        covs_predicted=covs_p,
        transport=transports,
        chi2=float(chi2),
        ndof=int(ndof),
        config=config,
        mass_gev=float(mass_gev),
    )


def state_position(state: np.ndarray) -> np.ndarray:
    return np.asarray(state[:3], dtype=float).copy()


def state_momentum(state: np.ndarray) -> np.ndarray:
    qop_t = float(state[IDX_QOP_T])
    if abs(qop_t) < 1.0e-12:
        pt = 1.0e12
    else:
        pt = 1.0 / abs(qop_t)
    phi = float(state[IDX_PHI])
    return np.array(
        [
            pt * np.cos(phi),
            pt * np.sin(phi),
            pt * float(state[IDX_TANL]),
        ],
        dtype=float,
    )


def first_smoothed_state(fit: TpcKalmanFit) -> np.ndarray:
    return np.asarray(fit.states_smoothed[0], dtype=float)


def state_at_path(fit: TpcKalmanFit, s_cm: float, use_smoothed: bool = True) -> np.ndarray:
    start = fit.states_smoothed[0] if use_smoothed else fit.states_filtered[0]
    return propagate_state(start, s_cm, fit.bfield_t, fit.config, fit.mass_gev)


def tangent_from_state(state: np.ndarray) -> np.ndarray:
    return np.array(
        [
            np.cos(state[IDX_PHI]),
            np.sin(state[IDX_PHI]),
            state[IDX_TANL],
        ],
        dtype=float,
    )


def line_line_pca_with_steps(pos1: np.ndarray, dir1: np.ndarray, pos2: np.ndarray, dir2: np.ndarray):
    w0 = pos1 - pos2
    a = float(np.dot(dir1, dir1))
    b = float(np.dot(dir1, dir2))
    c = float(np.dot(dir2, dir2))
    d = float(np.dot(dir1, w0))
    e = float(np.dot(dir2, w0))
    denom = a * c - b * b
    if abs(denom) < 1.0e-12:
        return None

    s = (b * e - c * d) / denom
    t = (a * e - b * d) / denom
    pca1 = pos1 + s * dir1
    pca2 = pos2 + t * dir2
    return pca1, pca2, float(np.linalg.norm(pca1 - pca2)), s, t


def kalman_search_range(
    fit: TpcKalmanFit,
    theta_extension: float,
    search_mode: str,
    downstream_margin: float,
    max_upstream_cm: float | None,
    downstream_margin_cm: float | None,
):
    if max_upstream_cm is not None and max_upstream_cm > 0.0:
        extension_cm = float(max_upstream_cm)
    else:
        extension_cm = abs(theta_extension) * fit.seed.radius

    if downstream_margin_cm is not None and downstream_margin_cm >= 0.0:
        downstream_cm = float(downstream_margin_cm)
    else:
        downstream_cm = abs(downstream_margin) * fit.seed.radius

    if search_mode == "full":
        return -extension_cm, float(fit.path_s[-1] + extension_cm)
    if search_mode != "upstream":
        raise ValueError(f"unknown Kalman PCA search mode: {search_mode}")
    return -extension_cm, downstream_cm


def refine_kalman_pair(
    fit1: TpcKalmanFit,
    fit2: TpcKalmanFit,
    s1: float,
    s2: float,
    min1: float,
    max1: float,
    min2: float,
    max2: float,
    max_step: float,
):
    pos1 = state_position(state_at_path(fit1, s1))
    pos2 = state_position(state_at_path(fit2, s2))
    best_dca2 = float(np.sum((pos1 - pos2) ** 2))

    for _ in range(30):
        state1 = state_at_path(fit1, s1)
        state2 = state_at_path(fit2, s2)
        pos1 = state_position(state1)
        pos2 = state_position(state2)
        tan1 = tangent_from_state(state1)
        tan2 = tangent_from_state(state2)
        line_pca = line_line_pca_with_steps(pos1, tan1, pos2, tan2)
        if line_pca is None:
            break

        _, _, _, step1, step2 = line_pca
        step1 = float(np.clip(step1, -max_step, max_step))
        step2 = float(np.clip(step2, -max_step, max_step))
        if abs(step1) < 1.0e-4 and abs(step2) < 1.0e-4:
            break

        candidate1 = float(np.clip(s1 + step1, min1, max1))
        candidate2 = float(np.clip(s2 + step2, min2, max2))
        cand_pos1 = state_position(state_at_path(fit1, candidate1))
        cand_pos2 = state_position(state_at_path(fit2, candidate2))
        candidate_dca2 = float(np.sum((cand_pos1 - cand_pos2) ** 2))
        if candidate_dca2 < best_dca2:
            s1, s2 = candidate1, candidate2
            best_dca2 = candidate_dca2
        else:
            max_step *= 0.5
            if max_step < 1.0e-3:
                break

    state1 = state_at_path(fit1, s1)
    state2 = state_at_path(fit2, s2)
    pca1 = state_position(state1)
    pca2 = state_position(state2)
    return pca1, pca2, float(np.linalg.norm(pca1 - pca2)), s1, s2


def kalman_track_pca_candidates(
    fit1: TpcKalmanFit,
    fit2: TpcKalmanFit,
    theta_extension: float,
    coarse_steps: int,
    search_mode: str,
    downstream_margin: float,
    max_candidates: int,
    max_upstream_cm: float | None = 80.0,
    downstream_margin_cm: float | None = 5.0,
):
    min1, max1 = kalman_search_range(
        fit1,
        theta_extension,
        search_mode,
        downstream_margin,
        max_upstream_cm,
        downstream_margin_cm,
    )
    min2, max2 = kalman_search_range(
        fit2,
        theta_extension,
        search_mode,
        downstream_margin,
        max_upstream_cm,
        downstream_margin_cm,
    )
    s1_values = np.linspace(min1, max1, max(8, coarse_steps))
    s2_values = np.linspace(min2, max2, max(8, coarse_steps))
    points1 = np.vstack([state_position(state_at_path(fit1, s)) for s in s1_values])
    points2 = np.vstack([state_position(state_at_path(fit2, s)) for s in s2_values])
    distances2 = np.sum((points1[:, None, :] - points2[None, :, :]) ** 2, axis=2)
    max_step = max((max1 - min1), (max2 - min2)) / max(8, coarse_steps)

    flat = distances2.ravel()
    n_candidates = min(max_candidates, flat.size)
    if n_candidates <= 0:
        return []

    candidate_indices = np.argpartition(flat, n_candidates - 1)[:n_candidates]
    candidates = []
    seen = set()
    for flat_index in candidate_indices:
        i, j = np.unravel_index(flat_index, distances2.shape)
        if (int(i), int(j)) in seen:
            continue
        seen.add((int(i), int(j)))
        candidates.append(
            refine_kalman_pair(
                fit1,
                fit2,
                float(s1_values[i]),
                float(s2_values[j]),
                min1,
                max1,
                min2,
                max2,
                max_step,
            )
        )

    candidates.sort(key=lambda item: item[2])
    return candidates


def kalman_dca_to_vertex(fit: TpcKalmanFit, vertex: np.ndarray) -> tuple[float, float]:
    state = first_smoothed_state(fit)
    omega = omega_from_state(state, fit.bfield_t)
    if abs(omega) < 1.0e-10:
        pos = state_position(state)
        mom = state_momentum(state)
        dxy = mom[:2]
        rel = pos - vertex
        pt2 = float(np.dot(dxy, dxy))
        if pt2 <= 0.0:
            return np.nan, np.nan
        dca_xy = abs(float(rel[0] * dxy[1] - rel[1] * dxy[0])) / np.sqrt(pt2)
        sxy = -float(np.dot(rel[:2], dxy)) / pt2
        closest = pos + sxy * mom
        return float(dca_xy), abs(float(closest[2] - vertex[2]))

    radius = 1.0 / omega
    center = np.array(
        [
            state[IDX_X] - np.sin(state[IDX_PHI]) / omega,
            state[IDX_Y] + np.cos(state[IDX_PHI]) / omega,
        ],
        dtype=float,
    )
    vec = np.asarray(vertex[:2], dtype=float) - center
    center_distance = float(np.linalg.norm(vec))
    if center_distance <= 0.0:
        theta_closest = np.arctan2(state[IDX_Y] - center[1], state[IDX_X] - center[0])
        dca_xy = abs(radius)
    else:
        closest_xy = center + abs(radius) * vec / center_distance
        theta0 = np.arctan2(state[IDX_Y] - center[1], state[IDX_X] - center[0])
        theta_raw = np.arctan2(closest_xy[1] - center[1], closest_xy[0] - center[0])
        dtheta_raw = normalize_phi(theta_raw - theta0)
        # Pick the nearest turn to the first smoothed state.
        theta_closest = theta0 + dtheta_raw
        dca_xy = abs(center_distance - abs(radius))

    theta0 = np.arctan2(state[IDX_Y] - center[1], state[IDX_X] - center[0])
    dtheta = normalize_phi(theta_closest - theta0)
    s_cm = dtheta / omega
    closest_state = propagate_state(state, s_cm, fit.bfield_t, fit.config, fit.mass_gev)
    return float(dca_xy), abs(float(closest_state[IDX_Z] - vertex[2]))
