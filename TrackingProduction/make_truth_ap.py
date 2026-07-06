#!/usr/bin/env python3
"""Build an Armenteros-Podolanski pairTree from ideal TPC truth points.

The input is the flat ROOT file written by PHG4TpcTruthPointTree.  The output
keeps the same branch names used by KshortReconstruction/drawAP.C so the
existing plotting macro can be reused.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import re

import numpy as np
import uproot

from tpc_kalman import (
    IDX_QOP_T,
    KalmanConfig,
    TpcKalmanFit,
    first_smoothed_state,
    fit_tpc_kalman,
    kalman_dca_to_vertex,
    kalman_track_pca_candidates,
    state_at_path,
    state_momentum,
    state_position,
)


@dataclass
class HelixFit:
    cx: float
    cy: float
    radius: float
    z0: float
    pitch: float
    theta_first: float
    theta_last: float
    theta_min: float
    theta_max: float
    direction: float
    bfield_t: float


@dataclass
class Tracklet:
    event: int
    track_id: int
    pid: int
    parent_id: int
    parent_pid: int
    charge: int
    pos: np.ndarray
    mom: np.ndarray
    truth_mom: np.ndarray
    vertex: np.ndarray
    truth_vertex: np.ndarray
    first_point: np.ndarray
    last_point: np.ndarray
    npoints: int
    helix: HelixFit | None = None
    kalman: TpcKalmanFit | None = None


@dataclass(frozen=True)
class InputSpec:
    path: str
    event_offset: int | None = None


@dataclass(frozen=True)
class PreselectionConfig:
    track_pt_min: float = 0.0
    track_dca_xy_min: float | None = None
    track_dca_z_min: float | None = None
    track_dca_xy_max: float | None = None
    track_dca_z_max: float | None = None
    pair_dca_max: float | None = None
    lproj_min: float | None = None
    cos_theta_min: float | None = None


PAIR_BRANCH_TYPES = {
    "run": np.int32,
    "evt": np.int32,
    "cross1": np.int16,
    "cross2": np.int16,
    "px1": np.float32,
    "py1": np.float32,
    "pz1": np.float32,
    "px2": np.float32,
    "py2": np.float32,
    "pz2": np.float32,
    "dca_xy1": np.float32,
    "dca_z1": np.float32,
    "dca_xy2": np.float32,
    "dca_z2": np.float32,
    "pairDCA": np.float32,
    "alpha": np.float32,
    "qT": np.float32,
    "charge1": np.float32,
    "charge2": np.float32,
    "cosThetaReco": np.float32,
    "Lproj": np.float32,
    "pca_x": np.float32,
    "pca_y": np.float32,
    "pca_z": np.float32,
    "pca1_x": np.float32,
    "pca1_y": np.float32,
    "pca1_z": np.float32,
    "pca2_x": np.float32,
    "pca2_y": np.float32,
    "pca2_z": np.float32,
    "mass_Kshort": np.float32,
    "mass_Lambda": np.float32,
    "mass_AntiLambda": np.float32,
    "true_decay_x": np.float32,
    "true_decay_y": np.float32,
    "true_decay_z": np.float32,
    "pca_to_true_3d": np.float32,
    "pca_to_true_xy": np.float32,
    "pca_to_true_z": np.float32,
    "truth_alpha": np.float32,
    "truth_qT": np.float32,
    "delta_alpha": np.float32,
    "delta_qT": np.float32,
    "truth_px1": np.float32,
    "truth_py1": np.float32,
    "truth_pz1": np.float32,
    "truth_px2": np.float32,
    "truth_py2": np.float32,
    "truth_pz2": np.float32,
    "cos_mom1_truth": np.float32,
    "cos_mom2_truth": np.float32,
    "pca_theta1": np.float32,
    "pca_theta2": np.float32,
    "track_id1": np.int32,
    "track_id2": np.int32,
    "pid1": np.int32,
    "pid2": np.int32,
    "parent_id1": np.int32,
    "parent_id2": np.int32,
    "parent_pid": np.int32,
    "npoints1": np.int16,
    "npoints2": np.int16,
}


TRACK_BRANCH_TYPES = {
    "run": np.int32,
    "evt": np.int32,
    "track_id": np.int32,
    "pid": np.int32,
    "parent_id": np.int32,
    "parent_pid": np.int32,
    "charge": np.int32,
    "npoints": np.int32,
    "has_helix": np.int32,
    "has_kalman": np.int32,
    "fit_method": np.int32,
    "px": np.float32,
    "py": np.float32,
    "pz": np.float32,
    "pt": np.float32,
    "p": np.float32,
    "x": np.float32,
    "y": np.float32,
    "z": np.float32,
    "first_x": np.float32,
    "first_y": np.float32,
    "first_z": np.float32,
    "first_r": np.float32,
    "last_x": np.float32,
    "last_y": np.float32,
    "last_z": np.float32,
    "last_r": np.float32,
    "dca_xy": np.float32,
    "dca_z": np.float32,
    "vertex_x": np.float32,
    "vertex_y": np.float32,
    "vertex_z": np.float32,
    "helix_cx": np.float32,
    "helix_cy": np.float32,
    "helix_radius": np.float32,
    "helix_z0": np.float32,
    "helix_pitch": np.float32,
    "helix_theta_first": np.float32,
    "helix_theta_last": np.float32,
    "helix_direction": np.float32,
    "kalman_seed_direction": np.float32,
    "kalman_expected_direction": np.float32,
    "kalman_charge_consistent": np.int32,
    "kalman_qop_t": np.float32,
    "kalman_chi2": np.float32,
    "kalman_ndof": np.int32,
    "kalman_chi2_ndf": np.float32,
    "truth_px": np.float32,
    "truth_py": np.float32,
    "truth_pz": np.float32,
    "cos_mom_truth": np.float32,
}


GSKIP_PATTERN = re.compile(r"(?:^|_)gskip(?P<gskip>\d+)(?:_|\.|$)")
PION_MASS_GEV = 0.13957039
PROTON_MASS_GEV = 0.93827208816
KALMAN_TRUTH_MEAS_SIGMA_CM = 1.0e-4
KALMAN_REALISTIC_MEAS_SIGMA_XY_CM = 0.03
KALMAN_REALISTIC_MEAS_SIGMA_Z_CM = 0.05


def parse_gskip_from_path(path: str | Path) -> int:
    match = GSKIP_PATTERN.search(Path(path).name)
    if not match:
        raise ValueError(
            f"could not parse gskip from input filename '{path}'. "
            "Expected a name like tpc_truthpoints_00042_gskip420_nev10.root"
        )
    return int(match.group("gskip"))


def resolve_list_entry(entry: str, list_path: Path) -> Path:
    path = Path(entry)
    if path.is_absolute():
        return path

    candidates = [Path.cwd() / path]
    candidates.extend(parent / path for parent in [list_path.parent, *list_path.parents])
    for candidate in candidates:
        if candidate.exists():
            return candidate

    return candidates[0]


def read_input_list(input_list: str) -> list[InputSpec]:
    list_path = Path(input_list).expanduser().resolve()
    specs: list[InputSpec] = []
    with list_path.open("r", encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            path = resolve_list_entry(line, list_path)
            specs.append(InputSpec(str(path), parse_gskip_from_path(path)))

    if not specs:
        raise ValueError(f"input list is empty: {input_list}")
    return specs


def event_ids_from_arrays(arrays: dict[str, np.ndarray], event_offset: int | None) -> np.ndarray:
    if event_offset is None:
        return arrays["event"].astype(np.int64)
    if "event_index" not in arrays:
        raise KeyError("event_index branch is required when using --input-list")
    return arrays["event_index"].astype(np.int64) + int(event_offset)


def rows_to_arrays(
    rows: list[dict[str, object]],
    branch_types: dict[str, object] = PAIR_BRANCH_TYPES,
) -> dict[str, np.ndarray]:
    return {
        branch: np.asarray([row[branch] for row in rows], dtype=dtype)
        for branch, dtype in branch_types.items()
    }


def pdg_charge(pid: int) -> int:
    """Return charge in units of |e| for common stable charged particles."""
    charge_by_abs_pid = {
        11: -1,     # e-
        13: -1,     # mu-
        211: 1,     # pi+
        321: 1,     # K+
        2212: 1,    # proton
        3222: 1,    # Sigma+
        3112: -1,   # Sigma-
        3312: -1,   # Xi-
        3334: -1,   # Omega-
    }
    base = charge_by_abs_pid.get(abs(int(pid)), 0)
    if pid < 0:
        base *= -1
    return base


def unit(vec: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0 or not np.isfinite(norm):
        return None
    return vec / norm


def line_line_pca(pos1: np.ndarray, mom1: np.ndarray, pos2: np.ndarray, mom2: np.ndarray):
    dir1 = unit(mom1)
    dir2 = unit(mom2)
    if dir1 is None or dir2 is None:
        return None

    w0 = pos1 - pos2
    a = float(np.dot(dir1, dir1))
    b = float(np.dot(dir1, dir2))
    c = float(np.dot(dir2, dir2))
    d = float(np.dot(dir1, w0))
    e = float(np.dot(dir2, w0))
    denom = a * c - b * b
    if abs(denom) < 1e-10:
        return None

    s = (b * e - c * d) / denom
    t = (a * e - b * d) / denom
    pca1 = pos1 + s * dir1
    pca2 = pos2 + t * dir2
    dca = float(np.linalg.norm(pca1 - pca2))
    return pca1, pca2, dca


def track_dca_to_vertex(pos: np.ndarray, mom: np.ndarray, vertex: np.ndarray):
    dxy = mom[:2]
    rel = pos - vertex
    pt2 = float(np.dot(dxy, dxy))
    if pt2 <= 0.0:
        return np.nan, np.nan

    dca_xy = abs(float(rel[0] * dxy[1] - rel[1] * dxy[0])) / np.sqrt(pt2)
    sxy = -float(np.dot(rel[:2], dxy)) / pt2
    closest = pos + sxy * mom
    dca_z = abs(float(closest[2] - vertex[2]))
    return dca_xy, dca_z


def track_pt(track: Tracklet) -> float:
    return float(np.linalg.norm(track.mom[:2]))


def approximate_track_dca(track: Tracklet, vertex: np.ndarray) -> tuple[float, float]:
    if track.kalman is not None:
        return kalman_dca_to_vertex(track.kalman, vertex)
    if track.helix is not None:
        return helix_dca_to_vertex(track.helix, vertex)
    return track_dca_to_vertex(track.pos, track.mom, vertex)


def pair_pointing(pair_vertex: np.ndarray, mom1: np.ndarray, mom2: np.ndarray, vertex: np.ndarray):
    flight = pair_vertex - vertex
    total_mom = mom1 + mom2
    flight_norm = float(np.linalg.norm(flight))
    mom_norm = float(np.linalg.norm(total_mom))
    if flight_norm <= 0.0 or mom_norm <= 0.0:
        return np.nan, np.nan
    cos_theta = float(np.dot(flight, total_mom) / (flight_norm * mom_norm))
    return flight_norm, cos_theta


def passes_preselection(
    track1: Tracklet,
    track2: Tracklet,
    primary_vertex: np.ndarray,
    config: PreselectionConfig,
) -> bool:
    if config.track_pt_min > 0.0:
        if track_pt(track1) < config.track_pt_min or track_pt(track2) < config.track_pt_min:
            return False

    dca_xy1, dca_z1 = approximate_track_dca(track1, primary_vertex)
    dca_xy2, dca_z2 = approximate_track_dca(track2, primary_vertex)
    for value in (dca_xy1, dca_z1, dca_xy2, dca_z2):
        if not np.isfinite(value):
            return False

    if config.track_dca_xy_min is not None:
        if dca_xy1 < config.track_dca_xy_min or dca_xy2 < config.track_dca_xy_min:
            return False
    if config.track_dca_z_min is not None:
        if dca_z1 < config.track_dca_z_min or dca_z2 < config.track_dca_z_min:
            return False
    if config.track_dca_xy_max is not None:
        if dca_xy1 > config.track_dca_xy_max or dca_xy2 > config.track_dca_xy_max:
            return False
    if config.track_dca_z_max is not None:
        if dca_z1 > config.track_dca_z_max or dca_z2 > config.track_dca_z_max:
            return False

    if (
        config.pair_dca_max is None
        and config.lproj_min is None
        and config.cos_theta_min is None
    ):
        return True

    rough_pca = line_line_pca(track1.pos, track1.mom, track2.pos, track2.mom)
    if rough_pca is None:
        return False
    rough_pca1, rough_pca2, rough_pair_dca = rough_pca

    if config.pair_dca_max is not None and rough_pair_dca > config.pair_dca_max:
        return False

    rough_vertex = 0.5 * (rough_pca1 + rough_pca2)
    rough_lproj, rough_cos_theta = pair_pointing(
        rough_vertex, track1.mom, track2.mom, primary_vertex
    )
    if config.lproj_min is not None:
        if not np.isfinite(rough_lproj) or rough_lproj < config.lproj_min:
            return False
    if config.cos_theta_min is not None:
        if not np.isfinite(rough_cos_theta) or rough_cos_theta < config.cos_theta_min:
            return False

    return True


def armenteros(pplus: np.ndarray, pminus: np.ndarray):
    v0p = pplus + pminus
    direction = unit(v0p)
    if direction is None:
        return None

    pl_plus = float(np.dot(pplus, direction))
    pl_minus = float(np.dot(pminus, direction))
    denom = pl_plus + pl_minus
    if abs(denom) < 1e-10:
        return None

    alpha = (pl_plus - pl_minus) / denom
    qt = float(np.linalg.norm(pplus - direction * pl_plus))
    return alpha, qt


def invariant_mass(momentum_a: np.ndarray, mass_a: float, momentum_b: np.ndarray, mass_b: float) -> float:
    momentum_a = np.asarray(momentum_a, dtype=float)
    momentum_b = np.asarray(momentum_b, dtype=float)
    energy_a = np.sqrt(float(np.dot(momentum_a, momentum_a)) + mass_a * mass_a)
    energy_b = np.sqrt(float(np.dot(momentum_b, momentum_b)) + mass_b * mass_b)
    total_momentum = momentum_a + momentum_b
    mass2 = (energy_a + energy_b) ** 2 - float(np.dot(total_momentum, total_momentum))
    return float(np.sqrt(max(mass2, 0.0)))


def vector_cosine(lhs: np.ndarray, rhs: np.ndarray) -> float:
    denom = float(np.linalg.norm(lhs) * np.linalg.norm(rhs))
    if denom <= 0.0 or not np.isfinite(denom):
        return np.nan
    return float(np.dot(lhs, rhs) / denom)


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


def order_track_indices(
    arrays: dict[str, np.ndarray],
    indices: list[int],
    point_order: str,
) -> np.ndarray:
    index_array = np.asarray(indices, dtype=np.int64)
    if index_array.size <= 2:
        return index_array

    if point_order == "path":
        return index_array[np.argsort(arrays["path"][index_array], kind="stable")]

    if point_order == "input":
        if "hit_key" in arrays:
            return index_array[np.argsort(arrays["hit_key"][index_array], kind="stable")]
        return index_array

    def inner_first(ordered_indices: np.ndarray) -> np.ndarray:
        if ordered_indices.size < 2:
            return ordered_indices
        first = int(ordered_indices[0])
        last = int(ordered_indices[-1])
        first_r = float(np.hypot(arrays["x"][first], arrays["y"][first]))
        last_r = float(np.hypot(arrays["x"][last], arrays["y"][last]))
        return ordered_indices[::-1] if first_r > last_r else ordered_indices

    def radius_order() -> np.ndarray:
        radius = np.hypot(arrays["x"][index_array], arrays["y"][index_array])
        if "layer" in arrays:
            return index_array[np.lexsort((arrays["layer"][index_array], radius))]
        return index_array[np.argsort(radius, kind="stable")]

    def theta_z_order() -> np.ndarray:
        points = np.column_stack(
            (
                arrays["x"][index_array],
                arrays["y"][index_array],
                arrays["z"][index_array],
            )
        ).astype(float)
        circle = fit_circle_least_squares(points)
        if circle is None:
            return order_track_indices(arrays, indices, "input" if "hit_key" in arrays else "path")

        cx, cy, _ = circle
        raw_theta = np.arctan2(points[:, 1] - cy, points[:, 0] - cx)
        z = points[:, 2]

        if float(np.ptp(z)) > 1.0e-3:
            z_order = np.argsort(z, kind="stable")
            unwrapped_theta = np.empty_like(raw_theta)
            unwrapped_theta[z_order] = np.unwrap(raw_theta[z_order])
            return inner_first(index_array[np.argsort(unwrapped_theta, kind="stable")])

        if "hit_key" in arrays:
            return inner_first(index_array[np.argsort(arrays["hit_key"][index_array], kind="stable")])

        return inner_first(index_array[np.argsort(raw_theta, kind="stable")])

    def looks_like_looper() -> bool:
        if "layer" in arrays:
            layers = arrays["layer"][index_array]
            if np.unique(layers).size < layers.size:
                return True

        sequence = order_track_indices(arrays, indices, "input" if "hit_key" in arrays else "path")
        radius = np.hypot(arrays["x"][sequence], arrays["y"][sequence])
        dr = np.diff(radius)
        dr = dr[np.abs(dr) > 0.5]
        if dr.size >= 2:
            signs = np.sign(dr)
            if np.any(signs[1:] * signs[:-1] < 0.0):
                return True

        points = np.column_stack((arrays["x"][index_array], arrays["y"][index_array], arrays["z"][index_array])).astype(float)
        circle = fit_circle_least_squares(points)
        if circle is not None and points.shape[0] >= 4:
            cx, cy, _ = circle
            theta = np.arctan2(points[:, 1] - cy, points[:, 0] - cx)
            z = points[:, 2]
            if float(np.ptp(z)) > 1.0e-3:
                z_order = np.argsort(z, kind="stable")
                theta_span = float(np.ptp(np.unwrap(theta[z_order])))
                if theta_span > 1.5 * np.pi:
                    return True

        return False

    if point_order == "radius":
        return radius_order()

    if point_order == "auto":
        return theta_z_order() if looks_like_looper() else radius_order()

    if point_order != "theta-z":
        raise ValueError(f"unknown point order mode: {point_order}")

    return theta_z_order()


def fit_helix(points: np.ndarray, bfield_t: float) -> HelixFit | None:
    circle = fit_circle_least_squares(points)
    if circle is None:
        return None
    cx, cy, radius = circle

    theta = np.unwrap(np.arctan2(points[:, 1] - cy, points[:, 0] - cx))
    if np.nanmax(theta) - np.nanmin(theta) < 1e-4:
        return None

    try:
        pitch, z0 = np.linalg.lstsq(
            np.column_stack((theta, np.ones_like(theta))),
            points[:, 2],
            rcond=None,
        )[0]
    except np.linalg.LinAlgError:
        return None

    direction = np.sign(theta[-1] - theta[0])
    if direction == 0.0:
        direction = 1.0

    return HelixFit(
        cx=cx,
        cy=cy,
        radius=radius,
        z0=float(z0),
        pitch=float(pitch),
        theta_first=float(theta[0]),
        theta_last=float(theta[-1]),
        theta_min=float(np.nanmin(theta)),
        theta_max=float(np.nanmax(theta)),
        direction=float(direction),
        bfield_t=float(bfield_t),
    )


def helix_point(helix: HelixFit, theta: float) -> np.ndarray:
    return np.array(
        [
            helix.cx + helix.radius * np.cos(theta),
            helix.cy + helix.radius * np.sin(theta),
            helix.z0 + helix.pitch * theta,
        ],
        dtype=float,
    )


def helix_tangent(helix: HelixFit, theta: float) -> np.ndarray:
    return np.array(
        [
            -helix.radius * np.sin(theta),
            helix.radius * np.cos(theta),
            helix.pitch,
        ],
        dtype=float,
    )


def helix_momentum(helix: HelixFit, theta: float) -> np.ndarray:
    p_t = 0.3 * abs(helix.bfield_t) * (helix.radius / 100.0)
    if p_t <= 0.0 or not np.isfinite(p_t):
        p_t = 1.0

    transverse = helix.direction * p_t * np.array([-np.sin(theta), np.cos(theta)], dtype=float)
    p_z = helix.direction * p_t * helix.pitch / helix.radius
    return np.array([transverse[0], transverse[1], p_z], dtype=float)


def theta_search_range(helix: HelixFit, extension: float, search_mode: str, downstream_margin: float):
    if search_mode == "full":
        return helix.theta_min - extension, helix.theta_max + extension

    if search_mode != "upstream":
        raise ValueError(f"unknown PCA search mode: {search_mode}")

    upstream = helix.theta_first - helix.direction * extension
    downstream = helix.theta_first + helix.direction * downstream_margin
    return min(upstream, downstream), max(upstream, downstream)


def line_line_pca_with_steps(pos1: np.ndarray, dir1: np.ndarray, pos2: np.ndarray, dir2: np.ndarray):
    w0 = pos1 - pos2
    a = float(np.dot(dir1, dir1))
    b = float(np.dot(dir1, dir2))
    c = float(np.dot(dir2, dir2))
    d = float(np.dot(dir1, w0))
    e = float(np.dot(dir2, w0))
    denom = a * c - b * b
    if abs(denom) < 1e-12:
        return None

    s = (b * e - c * d) / denom
    t = (a * e - b * d) / denom
    pca1 = pos1 + s * dir1
    pca2 = pos2 + t * dir2
    return pca1, pca2, float(np.linalg.norm(pca1 - pca2)), s, t


def helix_helix_pca(
    helix1: HelixFit,
    helix2: HelixFit,
    theta_extension: float,
    coarse_steps: int,
    search_mode: str,
    downstream_margin: float,
):
    candidates = helix_helix_pca_candidates(
        helix1,
        helix2,
        theta_extension,
        coarse_steps,
        search_mode,
        downstream_margin,
        max_candidates=1,
    )
    return candidates[0] if candidates else None


def refine_helix_pair(
    helix1: HelixFit,
    helix2: HelixFit,
    theta1: float,
    theta2: float,
    min1: float,
    max1: float,
    min2: float,
    max2: float,
    max_step: float,
):
    best_dca2 = float(np.sum((helix_point(helix1, theta1) - helix_point(helix2, theta2)) ** 2))
    for _ in range(30):
        pos1 = helix_point(helix1, theta1)
        pos2 = helix_point(helix2, theta2)
        tan1 = helix_tangent(helix1, theta1)
        tan2 = helix_tangent(helix2, theta2)
        line_pca = line_line_pca_with_steps(pos1, tan1, pos2, tan2)
        if line_pca is None:
            break

        _, _, _, step1, step2 = line_pca
        step1 = float(np.clip(step1, -max_step, max_step))
        step2 = float(np.clip(step2, -max_step, max_step))
        if abs(step1) < 1e-5 and abs(step2) < 1e-5:
            break

        candidate1 = float(np.clip(theta1 + step1, min1, max1))
        candidate2 = float(np.clip(theta2 + step2, min2, max2))
        candidate_dca2 = float(np.sum((helix_point(helix1, candidate1) - helix_point(helix2, candidate2)) ** 2))
        if candidate_dca2 < best_dca2:
            theta1, theta2 = candidate1, candidate2
            best_dca2 = candidate_dca2
        else:
            max_step *= 0.5
            if max_step < 1e-4:
                break

    pca1 = helix_point(helix1, theta1)
    pca2 = helix_point(helix2, theta2)
    return pca1, pca2, float(np.linalg.norm(pca1 - pca2)), theta1, theta2


def helix_helix_pca_candidates(
    helix1: HelixFit,
    helix2: HelixFit,
    theta_extension: float,
    coarse_steps: int,
    search_mode: str,
    downstream_margin: float,
    max_candidates: int,
):
    min1, max1 = theta_search_range(helix1, theta_extension, search_mode, downstream_margin)
    min2, max2 = theta_search_range(helix2, theta_extension, search_mode, downstream_margin)
    theta1_values = np.linspace(min1, max1, max(8, coarse_steps))
    theta2_values = np.linspace(min2, max2, max(8, coarse_steps))

    points1 = np.vstack([helix_point(helix1, theta) for theta in theta1_values])
    points2 = np.vstack([helix_point(helix2, theta) for theta in theta2_values])
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
        candidate = refine_helix_pair(
            helix1,
            helix2,
            float(theta1_values[i]),
            float(theta2_values[j]),
            min1,
            max1,
            min2,
            max2,
            max_step,
        )
        candidates.append(candidate)

    candidates.sort(key=lambda item: item[2])
    return candidates


def helix_dca_to_vertex(helix: HelixFit, vertex: np.ndarray):
    vec = np.array([vertex[0] - helix.cx, vertex[1] - helix.cy], dtype=float)
    distance_to_center = float(np.linalg.norm(vec))
    if distance_to_center <= 0.0:
        theta_raw = helix.theta_first
        dca_xy = helix.radius
    else:
        theta_raw = float(np.arctan2(vec[1], vec[0]))
        dca_xy = abs(distance_to_center - helix.radius)

    k = np.round((helix.theta_first - theta_raw) / (2.0 * np.pi))
    theta = theta_raw + 2.0 * np.pi * k
    closest = helix_point(helix, theta)
    return float(dca_xy), abs(float(closest[2] - vertex[2]))


def load_event_vertices_from_root(root_file, event_offset: int | None) -> dict[int, np.ndarray]:
    if "truth_particles" not in root_file:
        return {}

    branches = ["event", "is_primary", "vx", "vy", "vz"]
    if event_offset is not None:
        branches.append("event_index")
    arrays = root_file["truth_particles"].arrays(branches, library="np")
    event_ids = event_ids_from_arrays(arrays, event_offset)
    vertices: dict[int, np.ndarray] = {}
    for event in np.unique(event_ids):
        mask = (event_ids == event) & (arrays["is_primary"] != 0)
        if not np.any(mask):
            continue
        first = int(np.flatnonzero(mask)[0])
        vertices[int(event)] = np.array(
            [arrays["vx"][first], arrays["vy"][first], arrays["vz"][first]],
            dtype=float,
        )
    return vertices


def load_parent_pid_map_from_root(root_file, event_offset: int | None) -> dict[tuple[int, int], int]:
    if "truth_particles" not in root_file:
        return {}

    branches = ["event", "track_id", "pid"]
    if event_offset is not None:
        branches.append("event_index")
    arrays = root_file["truth_particles"].arrays(branches, library="np")
    event_ids = event_ids_from_arrays(arrays, event_offset)
    return {
        (int(event), int(track_id)): int(pid)
        for event, track_id, pid in zip(event_ids, arrays["track_id"], arrays["pid"])
    }


def load_tracklets_from_root(
    root_file,
    min_points: int,
    use_truth_kinematics: bool,
    fit_helix_tracks: bool,
    fit_kalman_tracks: bool,
    bfield_t: float,
    kalman_config: KalmanConfig,
    fit_first_points: int,
    event_offset: int | None,
    point_order: str,
) -> dict[int, list[Tracklet]]:
    tree = root_file["tpc_truth_points"]
    parent_pid_map = load_parent_pid_map_from_root(root_file, event_offset)
    branches = [
        "event",
        "track_id",
        "pid",
        "parent_id",
        "x",
        "y",
        "z",
        "px",
        "py",
        "pz",
        "truth_px",
        "truth_py",
        "truth_pz",
        "vx",
        "vy",
        "vz",
        "path",
    ]
    available_branches = set(tree.keys())
    if "hit_key" in available_branches:
        branches.append("hit_key")
    if "layer" in available_branches:
        branches.append("layer")
    if event_offset is not None:
        branches.append("event_index")
    arrays = tree.arrays(
        branches,
        library="np",
    )
    event_ids = event_ids_from_arrays(arrays, event_offset)

    grouped: dict[tuple[int, int], list[int]] = {}
    for index, key in enumerate(zip(event_ids, arrays["track_id"])):
        grouped.setdefault((int(key[0]), int(key[1])), []).append(index)

    by_event: dict[int, list[Tracklet]] = {}
    for (event, track_id), indices in grouped.items():
        if len(indices) < min_points:
            continue

        order = order_track_indices(arrays, indices, point_order)
        first = int(order[0])
        pid = int(arrays["pid"][first])
        parent_id = int(arrays["parent_id"][first])
        charge = pdg_charge(pid)
        if charge == 0:
            continue

        truth_mom = np.array(
            [arrays["truth_px"][first], arrays["truth_py"][first], arrays["truth_pz"][first]],
            dtype=float,
        )
        truth_vertex = np.array([arrays["vx"][first], arrays["vy"][first], arrays["vz"][first]], dtype=float)

        points = np.column_stack(
            (
                arrays["x"][order],
                arrays["y"][order],
                arrays["z"][order],
            )
        ).astype(float)
        fit_points = points
        if fit_first_points > 0:
            fit_points = points[:fit_first_points]
        helix = fit_helix(fit_points, bfield_t) if fit_helix_tracks else None
        if fit_helix_tracks and helix is None:
            continue

        kalman = fit_tpc_kalman(points, charge, kalman_config, pid=pid) if fit_kalman_tracks else None
        if fit_kalman_tracks and kalman is None:
            continue

        if use_truth_kinematics:
            pos = truth_vertex
            mom = truth_mom
        elif fit_kalman_tracks and kalman is not None:
            first_state = first_smoothed_state(kalman)
            pos = state_position(first_state)
            mom = state_momentum(first_state)
        elif fit_helix_tracks and helix is not None:
            pos = helix_point(helix, helix.theta_first)
            mom = helix_momentum(helix, helix.theta_first)
        else:
            pos = np.array([arrays["x"][first], arrays["y"][first], arrays["z"][first]], dtype=float)
            mom = np.array([arrays["px"][first], arrays["py"][first], arrays["pz"][first]], dtype=float)

        tracklet = Tracklet(
            event=event,
            track_id=track_id,
            pid=pid,
            parent_id=parent_id,
            parent_pid=parent_pid_map.get((event, parent_id), 0),
            charge=charge,
            pos=pos,
            mom=mom,
            truth_mom=truth_mom,
            vertex=truth_vertex,
            truth_vertex=truth_vertex,
            first_point=points[0],
            last_point=points[-1],
            npoints=len(indices),
            helix=helix,
            kalman=kalman,
        )
        by_event.setdefault(event, []).append(tracklet)

    return by_event


def make_track_row(track: Tracklet, primary_vertex: np.ndarray, run: int) -> dict[str, object]:
    dca_xy, dca_z = approximate_track_dca(track, primary_vertex)
    fit_method = 2 if track.kalman is not None else 1 if track.helix is not None else 0

    helix_cx = helix_cy = helix_radius = np.nan
    helix_z0 = helix_pitch = np.nan
    helix_theta_first = helix_theta_last = helix_direction = np.nan
    if track.helix is not None:
        helix_cx = track.helix.cx
        helix_cy = track.helix.cy
        helix_radius = track.helix.radius
        helix_z0 = track.helix.z0
        helix_pitch = track.helix.pitch
        helix_theta_first = track.helix.theta_first
        helix_theta_last = track.helix.theta_last
        helix_direction = track.helix.direction

    kalman_chi2 = np.nan
    kalman_ndof = -1
    kalman_chi2_ndf = np.nan
    kalman_seed_direction = np.nan
    kalman_expected_direction = np.nan
    kalman_charge_consistent = -1
    kalman_qop_t = np.nan
    if track.kalman is not None:
        kalman_chi2 = track.kalman.chi2
        kalman_ndof = track.kalman.ndof
        if track.kalman.ndof > 0:
            kalman_chi2_ndf = track.kalman.chi2 / track.kalman.ndof
        kalman_seed_direction = track.kalman.seed.direction
        if track.charge != 0 and track.kalman.bfield_t != 0.0:
            kalman_expected_direction = -float(track.charge) * float(np.sign(track.kalman.bfield_t))
            kalman_charge_consistent = int(kalman_seed_direction == kalman_expected_direction)
        kalman_qop_t = float(first_smoothed_state(track.kalman)[IDX_QOP_T])

    return {
        "run": run,
        "evt": track.event,
        "track_id": track.track_id,
        "pid": track.pid,
        "parent_id": track.parent_id,
        "parent_pid": track.parent_pid,
        "charge": track.charge,
        "npoints": track.npoints,
        "has_helix": 1 if track.helix is not None else 0,
        "has_kalman": 1 if track.kalman is not None else 0,
        "fit_method": fit_method,
        "px": track.mom[0],
        "py": track.mom[1],
        "pz": track.mom[2],
        "pt": float(np.linalg.norm(track.mom[:2])),
        "p": float(np.linalg.norm(track.mom)),
        "x": track.pos[0],
        "y": track.pos[1],
        "z": track.pos[2],
        "first_x": track.first_point[0],
        "first_y": track.first_point[1],
        "first_z": track.first_point[2],
        "first_r": float(np.linalg.norm(track.first_point[:2])),
        "last_x": track.last_point[0],
        "last_y": track.last_point[1],
        "last_z": track.last_point[2],
        "last_r": float(np.linalg.norm(track.last_point[:2])),
        "dca_xy": dca_xy,
        "dca_z": dca_z,
        "vertex_x": primary_vertex[0],
        "vertex_y": primary_vertex[1],
        "vertex_z": primary_vertex[2],
        "helix_cx": helix_cx,
        "helix_cy": helix_cy,
        "helix_radius": helix_radius,
        "helix_z0": helix_z0,
        "helix_pitch": helix_pitch,
        "helix_theta_first": helix_theta_first,
        "helix_theta_last": helix_theta_last,
        "helix_direction": helix_direction,
        "kalman_seed_direction": kalman_seed_direction,
        "kalman_expected_direction": kalman_expected_direction,
        "kalman_charge_consistent": kalman_charge_consistent,
        "kalman_qop_t": kalman_qop_t,
        "kalman_chi2": kalman_chi2,
        "kalman_ndof": kalman_ndof,
        "kalman_chi2_ndf": kalman_chi2_ndf,
        "truth_px": track.truth_mom[0],
        "truth_py": track.truth_mom[1],
        "truth_pz": track.truth_mom[2],
        "cos_mom_truth": vector_cosine(track.mom, track.truth_mom),
    }


def make_pair_row(
    track1: Tracklet,
    track2: Tracklet,
    primary_vertex: np.ndarray,
    run: int,
    same_parent_only: bool,
    parent_pids: set[int],
    fit_helix_tracks: bool,
    fit_kalman_tracks: bool,
    theta_extension: float,
    coarse_steps: int,
    pca_search: str,
    downstream_margin: float,
    kalman_max_upstream_cm: float | None,
    kalman_downstream_margin_cm: float | None,
    pca_candidates: int,
    prefer_positive_pointing: bool,
    preselection: PreselectionConfig,
    counters: dict[str, int] | None = None,
):
    if counters is not None:
        counters["raw_pairs"] = counters.get("raw_pairs", 0) + 1

    if track1.charge == track2.charge:
        if counters is not None:
            counters["reject_charge"] = counters.get("reject_charge", 0) + 1
        return None
    if same_parent_only and (track1.parent_id == 0 or track1.parent_id != track2.parent_id):
        if counters is not None:
            counters["reject_parent"] = counters.get("reject_parent", 0) + 1
        return None
    if parent_pids:
        if track1.parent_id == 0 or track1.parent_id != track2.parent_id:
            if counters is not None:
                counters["reject_parent"] = counters.get("reject_parent", 0) + 1
            return None
        if track1.parent_pid not in parent_pids:
            if counters is not None:
                counters["reject_parent_pid"] = counters.get("reject_parent_pid", 0) + 1
            return None

    if not passes_preselection(track1, track2, primary_vertex, preselection):
        if counters is not None:
            counters["reject_preselection"] = counters.get("reject_preselection", 0) + 1
        return None

    mom1 = track1.mom
    mom2 = track2.mom
    theta1 = np.nan
    theta2 = np.nan
    if fit_kalman_tracks and track1.kalman is not None and track2.kalman is not None:
        candidates = kalman_track_pca_candidates(
            track1.kalman,
            track2.kalman,
            theta_extension,
            coarse_steps,
            pca_search,
            downstream_margin,
            pca_candidates,
            kalman_max_upstream_cm,
            kalman_downstream_margin_cm,
        )
        if not candidates:
            if counters is not None:
                counters["reject_pca"] = counters.get("reject_pca", 0) + 1
            return None

        if prefer_positive_pointing:
            ranked_candidates = []
            for candidate in candidates:
                cand_pca1, cand_pca2, cand_dca, cand_s1, cand_s2 = candidate
                cand_mom1 = state_momentum(state_at_path(track1.kalman, cand_s1))
                cand_mom2 = state_momentum(state_at_path(track2.kalman, cand_s2))
                cand_vertex = 0.5 * (cand_pca1 + cand_pca2)
                cand_flight = cand_vertex - primary_vertex
                cand_total_mom = cand_mom1 + cand_mom2
                if np.linalg.norm(cand_total_mom) <= 0.0 or np.linalg.norm(cand_flight) <= 0.0:
                    cand_cos = -2.0
                else:
                    cand_cos = float(
                        np.dot(cand_flight, cand_total_mom) /
                        (np.linalg.norm(cand_flight) * np.linalg.norm(cand_total_mom))
                    )
                pointing_penalty = 0.0 if cand_cos > 0.0 else 1000.0
                ranked_candidates.append((pointing_penalty, cand_dca, -cand_cos, candidate))

            ranked_candidates.sort(key=lambda item: (item[0], item[1], item[2]))
            pca = ranked_candidates[0][3]
        else:
            pca = min(candidates, key=lambda candidate: candidate[2])
        pca1, pca2, pair_dca, theta1, theta2 = pca
        mom1 = state_momentum(state_at_path(track1.kalman, theta1))
        mom2 = state_momentum(state_at_path(track2.kalman, theta2))
        dca_xy1, dca_z1 = kalman_dca_to_vertex(track1.kalman, primary_vertex)
        dca_xy2, dca_z2 = kalman_dca_to_vertex(track2.kalman, primary_vertex)
    elif fit_helix_tracks and track1.helix is not None and track2.helix is not None:
        candidates = helix_helix_pca_candidates(
            track1.helix,
            track2.helix,
            theta_extension,
            coarse_steps,
            pca_search,
            downstream_margin,
            pca_candidates,
        )
        if not candidates:
            if counters is not None:
                counters["reject_pca"] = counters.get("reject_pca", 0) + 1
            return None

        if prefer_positive_pointing:
            ranked_candidates = []
            for candidate in candidates:
                cand_pca1, cand_pca2, cand_dca, cand_theta1, cand_theta2 = candidate
                cand_mom1 = helix_momentum(track1.helix, cand_theta1)
                cand_mom2 = helix_momentum(track2.helix, cand_theta2)
                cand_vertex = 0.5 * (cand_pca1 + cand_pca2)
                cand_flight = cand_vertex - primary_vertex
                cand_total_mom = cand_mom1 + cand_mom2
                if np.linalg.norm(cand_total_mom) <= 0.0 or np.linalg.norm(cand_flight) <= 0.0:
                    cand_cos = -2.0
                else:
                    cand_cos = float(
                        np.dot(cand_flight, cand_total_mom) /
                        (np.linalg.norm(cand_flight) * np.linalg.norm(cand_total_mom))
                    )
                pointing_penalty = 0.0 if cand_cos > 0.0 else 1000.0
                ranked_candidates.append((pointing_penalty, cand_dca, -cand_cos, candidate))

            ranked_candidates.sort(key=lambda item: (item[0], item[1], item[2]))
            pca = ranked_candidates[0][3]
        else:
            pca = min(candidates, key=lambda candidate: candidate[2])
        pca1, pca2, pair_dca, theta1, theta2 = pca
        mom1 = helix_momentum(track1.helix, theta1)
        mom2 = helix_momentum(track2.helix, theta2)
        dca_xy1, dca_z1 = helix_dca_to_vertex(track1.helix, primary_vertex)
        dca_xy2, dca_z2 = helix_dca_to_vertex(track2.helix, primary_vertex)
    else:
        pca = line_line_pca(track1.pos, track1.mom, track2.pos, track2.mom)
        if pca is None:
            if counters is not None:
                counters["reject_pca"] = counters.get("reject_pca", 0) + 1
            return None
        pca1, pca2, pair_dca = pca
        dca_xy1, dca_z1 = track_dca_to_vertex(track1.pos, track1.mom, primary_vertex)
        dca_xy2, dca_z2 = track_dca_to_vertex(track2.pos, track2.mom, primary_vertex)
    pair_vertex = 0.5 * (pca1 + pca2)

    total_mom = mom1 + mom2
    flight = pair_vertex - primary_vertex
    if np.linalg.norm(total_mom) <= 0.0 or np.linalg.norm(flight) <= 0.0:
        if counters is not None:
            counters["reject_pointing"] = counters.get("reject_pointing", 0) + 1
        return None
    cos_theta = float(np.dot(flight, total_mom) / (np.linalg.norm(flight) * np.linalg.norm(total_mom)))

    if track1.charge > 0:
        pplus, pminus = mom1, mom2
    else:
        pplus, pminus = mom2, mom1
    ap = armenteros(pplus, pminus)
    if ap is None:
        if counters is not None:
            counters["reject_ap"] = counters.get("reject_ap", 0) + 1
        return None
    alpha, qt = ap
    mass_kshort = invariant_mass(pplus, PION_MASS_GEV, pminus, PION_MASS_GEV)
    mass_lambda = invariant_mass(pplus, PROTON_MASS_GEV, pminus, PION_MASS_GEV)
    mass_antilambda = invariant_mass(pplus, PION_MASS_GEV, pminus, PROTON_MASS_GEV)

    if track1.charge > 0:
        truth_pplus, truth_pminus = track1.truth_mom, track2.truth_mom
    else:
        truth_pplus, truth_pminus = track2.truth_mom, track1.truth_mom
    truth_ap = armenteros(truth_pplus, truth_pminus)
    if truth_ap is None:
        truth_alpha = np.nan
        truth_qt = np.nan
    else:
        truth_alpha, truth_qt = truth_ap

    if track1.parent_id != 0 and track1.parent_id == track2.parent_id:
        true_decay = 0.5 * (track1.truth_vertex + track2.truth_vertex)
        pca_delta = pair_vertex - true_decay
        pca_to_true_3d = float(np.linalg.norm(pca_delta))
        pca_to_true_xy = float(np.linalg.norm(pca_delta[:2]))
        pca_to_true_z = float(abs(pca_delta[2]))
    else:
        true_decay = np.array([np.nan, np.nan, np.nan], dtype=float)
        pca_to_true_3d = np.nan
        pca_to_true_xy = np.nan
        pca_to_true_z = np.nan

    return {
        "run": run,
        "evt": track1.event,
        "cross1": 0,
        "cross2": 0,
        "px1": mom1[0],
        "py1": mom1[1],
        "pz1": mom1[2],
        "px2": mom2[0],
        "py2": mom2[1],
        "pz2": mom2[2],
        "dca_xy1": dca_xy1,
        "dca_z1": dca_z1,
        "dca_xy2": dca_xy2,
        "dca_z2": dca_z2,
        "pairDCA": pair_dca,
        "alpha": alpha,
        "qT": qt,
        "charge1": track1.charge,
        "charge2": track2.charge,
        "cosThetaReco": cos_theta,
        "Lproj": float(np.linalg.norm(flight)),
        "pca_x": pair_vertex[0],
        "pca_y": pair_vertex[1],
        "pca_z": pair_vertex[2],
        "pca1_x": pca1[0],
        "pca1_y": pca1[1],
        "pca1_z": pca1[2],
        "pca2_x": pca2[0],
        "pca2_y": pca2[1],
        "pca2_z": pca2[2],
        "mass_Kshort": mass_kshort,
        "mass_Lambda": mass_lambda,
        "mass_AntiLambda": mass_antilambda,
        "true_decay_x": true_decay[0],
        "true_decay_y": true_decay[1],
        "true_decay_z": true_decay[2],
        "pca_to_true_3d": pca_to_true_3d,
        "pca_to_true_xy": pca_to_true_xy,
        "pca_to_true_z": pca_to_true_z,
        "truth_alpha": truth_alpha,
        "truth_qT": truth_qt,
        "delta_alpha": alpha - truth_alpha if np.isfinite(truth_alpha) else np.nan,
        "delta_qT": qt - truth_qt if np.isfinite(truth_qt) else np.nan,
        "truth_px1": track1.truth_mom[0],
        "truth_py1": track1.truth_mom[1],
        "truth_pz1": track1.truth_mom[2],
        "truth_px2": track2.truth_mom[0],
        "truth_py2": track2.truth_mom[1],
        "truth_pz2": track2.truth_mom[2],
        "cos_mom1_truth": vector_cosine(mom1, track1.truth_mom),
        "cos_mom2_truth": vector_cosine(mom2, track2.truth_mom),
        "pca_theta1": theta1,
        "pca_theta2": theta2,
        "track_id1": track1.track_id,
        "track_id2": track2.track_id,
        "pid1": track1.pid,
        "pid2": track2.pid,
        "parent_id1": track1.parent_id,
        "parent_id2": track2.parent_id,
        "parent_pid": track1.parent_pid if track1.parent_id == track2.parent_id else 0,
        "npoints1": track1.npoints,
        "npoints2": track2.npoints,
    }


def build_pair_tree(
    input_specs: list[InputSpec],
    output_file: str,
    primary_vertex: np.ndarray,
    min_points: int,
    run: int,
    same_parent_only: bool,
    parent_pids: set[int],
    use_truth_primary_vertex: bool,
    use_truth_kinematics: bool,
    fit_helix_tracks: bool,
    fit_kalman_tracks: bool,
    bfield_t: float,
    kalman_config: KalmanConfig,
    theta_extension: float,
    coarse_steps: int,
    pca_search: str,
    downstream_margin: float,
    kalman_max_upstream_cm: float | None,
    kalman_downstream_margin_cm: float | None,
    fit_first_points: int,
    pca_candidates: int,
    prefer_positive_pointing: bool,
    write_chunk_size: int,
    progress_every: int,
    preselection: PreselectionConfig,
    point_order: str,
):
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    rows_buffer: list[dict[str, object]] = []
    track_rows_buffer: list[dict[str, object]] = []
    n_rows = 0
    n_track_rows = 0
    visible_events: set[int] = set()
    pca_residuals: list[float] = []
    delta_alpha_values: list[float] = []
    delta_qt_values: list[float] = []
    counters: dict[str, int] = {}

    def record_stats(rows: list[dict[str, object]]) -> None:
        for row in rows:
            pca_residual = float(row["pca_to_true_3d"])
            if np.isfinite(pca_residual):
                pca_residuals.append(pca_residual)
            delta_alpha = float(row["delta_alpha"])
            if np.isfinite(delta_alpha):
                delta_alpha_values.append(delta_alpha)
            delta_qt = float(row["delta_qT"])
            if np.isfinite(delta_qt):
                delta_qt_values.append(delta_qt)

    with uproot.recreate(output_file) as root_file:
        root_file.mktree("pairTree", PAIR_BRANCH_TYPES)
        root_file.mktree("trackTree", TRACK_BRANCH_TYPES)

        def flush_rows() -> None:
            nonlocal n_rows
            if not rows_buffer:
                return
            root_file["pairTree"].extend(rows_to_arrays(rows_buffer))
            record_stats(rows_buffer)
            n_rows += len(rows_buffer)
            rows_buffer.clear()

        def flush_track_rows() -> None:
            nonlocal n_track_rows
            if not track_rows_buffer:
                return
            root_file["trackTree"].extend(rows_to_arrays(track_rows_buffer, TRACK_BRANCH_TYPES))
            n_track_rows += len(track_rows_buffer)
            track_rows_buffer.clear()

        for input_index, spec in enumerate(input_specs, start=1):
            root_input = uproot.open(spec.path)
            by_event = load_tracklets_from_root(
                root_input,
                min_points,
                use_truth_kinematics,
                fit_helix_tracks,
                fit_kalman_tracks,
                bfield_t,
                kalman_config,
                fit_first_points,
                spec.event_offset,
                point_order,
            )
            event_vertices = (
                load_event_vertices_from_root(root_input, spec.event_offset)
                if use_truth_primary_vertex
                else {}
            )
            visible_events.update(by_event.keys())

            for event in sorted(by_event):
                event_vertex = event_vertices.get(event, primary_vertex)
                for track in by_event[event]:
                    track_rows_buffer.append(make_track_row(track, event_vertex, run))
                    if len(track_rows_buffer) >= write_chunk_size:
                        flush_track_rows()

                for track1, track2 in combinations(by_event[event], 2):
                    row = make_pair_row(
                        track1,
                        track2,
                        event_vertex,
                        run,
                        same_parent_only,
                        parent_pids,
                        fit_helix_tracks,
                        fit_kalman_tracks,
                        theta_extension,
                        coarse_steps,
                        pca_search,
                        downstream_margin,
                        kalman_max_upstream_cm,
                        kalman_downstream_margin_cm,
                        pca_candidates,
                        prefer_positive_pointing,
                        preselection,
                        counters,
                    )
                    if row is not None and all(
                        np.isfinite(row[name])
                        for name in ["pairDCA", "alpha", "qT", "cosThetaReco", "Lproj"]
                    ):
                        rows_buffer.append(row)
                        if len(rows_buffer) >= write_chunk_size:
                            flush_rows()

            if progress_every > 0 and (
                input_index == 1
                or input_index == len(input_specs)
                or input_index % progress_every == 0
            ):
                print(
                    "[make_truth_ap] processed "
                    f"{input_index}/{len(input_specs)} input files; "
                    f"rows written/buffered={n_rows + len(rows_buffer)}"
                )

        flush_rows()
        flush_track_rows()

    print(f"[make_truth_ap] wrote {n_rows} pairs to {output_file}")
    print(f"[make_truth_ap] wrote {n_track_rows} track QA rows to trackTree")
    print(f"[make_truth_ap] events with visible charged tracklets: {len(visible_events)}")
    if counters:
        print(
            "[make_truth_ap] pair counters: "
            + ", ".join(f"{key}={value}" for key, value in sorted(counters.items()))
        )
    if n_rows:
        finite_pca = np.asarray(pca_residuals, dtype=float)
        if finite_pca.size:
            print(
                "[make_truth_ap] pca_to_true_3d cm: "
                f"median={np.median(finite_pca):.4g}, "
                f"max={np.max(finite_pca):.4g}, "
                f"n={finite_pca.size}"
            )
        finite_alpha = np.asarray(delta_alpha_values, dtype=float)
        if finite_alpha.size:
            print(
                "[make_truth_ap] delta_alpha: "
                f"median={np.median(finite_alpha):.4g}, "
                f"max_abs={np.max(np.abs(finite_alpha)):.4g}"
            )
        finite_qt = np.asarray(delta_qt_values, dtype=float)
        if finite_qt.size:
            print(
                "[make_truth_ap] delta_qT GeV/c: "
                f"median={np.median(finite_qt):.4g}, "
                f"max_abs={np.max(np.abs(finite_qt)):.4g}"
            )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create KshortReconstruction-style AP pairTree from ideal TPC truth points.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input_file", nargs="?", help="ROOT file from PHG4TpcTruthPointTree")
    parser.add_argument(
        "--input-list",
        help=(
            "Text file containing one PHG4TpcTruthPointTree ROOT file per line. "
            "Each filename must contain gskipNNN; the script uses gskip + event_index "
            "as the global event id."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default="ap_truth_points.root",
        help="Output ROOT file containing pairTree.",
    )
    parser.add_argument(
        "--vertex",
        nargs=3,
        type=float,
        default=(0.0, 0.0, 0.0),
        metavar=("X", "Y", "Z"),
        help="Primary vertex assumption in cm.",
    )
    parser.add_argument("--min-points", type=int, default=5, help="Minimum TPC truth points per tracklet.")
    parser.add_argument(
        "--point-order",
        choices=["path", "input", "radius", "theta-z", "auto"],
        default="path",
        help=(
            "Ordering used for points within each track before helix/Kalman fitting. "
            "'path' uses the input path branch; 'input' uses hit_key/original cluster order; "
            "'radius' sorts inner-to-outer; 'theta-z' unwraps azimuth using z to handle loopers; "
            "'auto' uses radius for normal tracks and theta-z for looper-like tracks."
        ),
    )
    parser.add_argument("--run", type=int, default=1, help="Run number stored in pairTree.")
    parser.add_argument(
        "--same-parent-only",
        action="store_true",
        help="Keep only opposite-sign pairs whose truth particles have the same nonzero parent_id.",
    )
    parser.add_argument(
        "--parent-pids",
        nargs="*",
        type=int,
        default=(),
        help="When set, keep only same-parent pairs with one of these parent PDG IDs, for example 310 3122 -3122.",
    )
    parser.add_argument(
        "--use-truth-primary-vertex",
        action="store_true",
        help="Use the true primary vertex stored in truth_particles for each event instead of the fixed --vertex.",
    )
    parser.add_argument(
        "--use-truth-kinematics",
        action="store_true",
        help="Use each particle's truth production vertex and truth momentum instead of the first TPC point.",
    )
    parser.add_argument(
        "--fit-helix",
        action="store_true",
        help="Fit a helix to each TPC truth-point tracklet and use helix-helix PCA.",
    )
    parser.add_argument(
        "--fit-kalman",
        action="store_true",
        help=(
            "Fit each grouped TPC tracklet with a standalone TPC-only EKF/smoother "
            "and use curved Kalman-track PCA. This does not run pattern recognition."
        ),
    )
    parser.add_argument(
        "--fit-first-points",
        type=int,
        default=0,
        help="If positive, fit each helix using only the first N TPC points along the track.",
    )
    parser.add_argument(
        "--bfield",
        type=float,
        default=1.4,
        help="Approximate solenoidal magnetic field in Tesla used to convert fitted radius to pT.",
    )
    parser.add_argument(
        "--kalman-measurement-preset",
        choices=("custom", "truth", "realistic"),
        default="custom",
        help=(
            "Measurement covariance preset for --fit-kalman. 'truth' uses a tiny "
            "nonzero 1 micron sigma for exact truth points; 'realistic' uses the "
            "rough TPC cluster-resolution defaults; 'custom' uses the explicit "
            "--kalman-meas-sigma-* values."
        ),
    )
    parser.add_argument(
        "--kalman-meas-sigma-xy",
        type=float,
        default=0.03,
        help="TPC-only Kalman measurement sigma for x/y points in cm, used with --kalman-measurement-preset custom.",
    )
    parser.add_argument(
        "--kalman-meas-sigma-z",
        type=float,
        default=0.05,
        help="TPC-only Kalman measurement sigma for z points in cm, used with --kalman-measurement-preset custom.",
    )
    parser.add_argument(
        "--kalman-process-sigma-pos",
        type=float,
        default=1.0e-4,
        help="TPC-only Kalman process noise per cm for x/y/z positions in cm.",
    )
    parser.add_argument(
        "--kalman-process-sigma-phi",
        type=float,
        default=1.0e-5,
        help="TPC-only Kalman process noise per cm for phi.",
    )
    parser.add_argument(
        "--kalman-process-sigma-qop",
        type=float,
        default=1.0e-6,
        help="TPC-only Kalman process noise per cm for signed q/pT.",
    )
    parser.add_argument(
        "--kalman-process-sigma-tanl",
        type=float,
        default=1.0e-6,
        help="TPC-only Kalman process noise per cm for tan(lambda).",
    )
    parser.add_argument(
        "--kalman-material-x0-per-cm",
        type=float,
        default=0.0,
        help=(
            "Material thickness per cm, x/X0/cm, used for Highland multiple-scattering "
            "process noise. Use 0 to disable the physics term."
        ),
    )
    parser.add_argument(
        "--kalman-ms-scale",
        type=float,
        default=1.0,
        help="Scale factor applied to the Highland multiple-scattering angle.",
    )
    parser.add_argument(
        "--kalman-energy-loss-gev-per-cm",
        type=float,
        default=0.0,
        help=(
            "Mean energy loss per 3D path length in GeV/cm used during Kalman propagation. "
            "Use 0 to disable."
        ),
    )
    parser.add_argument(
        "--kalman-energy-loss-sigma-fraction",
        type=float,
        default=0.0,
        help="Optional fractional straggling noise on the configured mean energy loss.",
    )
    parser.add_argument(
        "--theta-extension",
        type=float,
        default=1.5,
        help="Extra helix theta range, in radians, searched beyond the measured TPC points.",
    )
    parser.add_argument(
        "--pca-search",
        choices=("upstream", "full"),
        default="upstream",
        help="Helix PCA search range. Use upstream for V0 candidates before the first TPC point; full searches the full measured helix plus extension.",
    )
    parser.add_argument(
        "--downstream-margin",
        type=float,
        default=0.2,
        help="Small theta margin after the first TPC point used with --pca-search upstream.",
    )
    parser.add_argument(
        "--kalman-max-upstream-cm",
        type=float,
        default=80.0,
        help=(
            "For --fit-kalman, maximum inward propagation distance in cm from the "
            "innermost TPC point when searching the two-track PCA. This replaces "
            "the old theta_extension*radius range for Kalman tracks."
        ),
    )
    parser.add_argument(
        "--kalman-downstream-margin-cm",
        type=float,
        default=5.0,
        help=(
            "For --fit-kalman, small downstream allowance in cm after the innermost "
            "TPC point during PCA search."
        ),
    )
    parser.add_argument(
        "--coarse-steps",
        type=int,
        default=48,
        help="Number of coarse theta samples per track used before PCA refinement.",
    )
    parser.add_argument(
        "--pca-candidates",
        type=int,
        default=32,
        help="Number of coarse helix-helix PCA seeds to refine before choosing a positive-pointing candidate.",
    )
    parser.add_argument(
        "--prefer-positive-pointing",
        action="store_true",
        help="Prefer positive-pointing helix PCA candidates before minimizing pair DCA. Useful for studies, but not the default.",
    )
    parser.add_argument(
        "--pre-track-pt-min",
        type=float,
        default=0.0,
        help="Cheap preselection: require both tracks to have pT above this value in GeV/c before helix PCA.",
    )
    parser.add_argument(
        "--pre-track-dca-xy-min",
        type=float,
        default=None,
        help="Cheap preselection: require both tracks to have |DCA_xy| above this value in cm before helix PCA.",
    )
    parser.add_argument(
        "--pre-track-dca-z-min",
        type=float,
        default=None,
        help="Cheap preselection: require both tracks to have |DCA_z| above this value in cm before helix PCA.",
    )
    parser.add_argument(
        "--pre-track-dca-xy-max",
        type=float,
        default=None,
        help="Cheap preselection: reject tracks with |DCA_xy| above this value in cm before helix PCA.",
    )
    parser.add_argument(
        "--pre-track-dca-z-max",
        type=float,
        default=None,
        help="Cheap preselection: reject tracks with |DCA_z| above this value in cm before helix PCA.",
    )
    parser.add_argument(
        "--pre-pair-dca-max",
        type=float,
        default=None,
        help=(
            "Cheap preselection: reject pairs whose rough straight-line pair DCA exceeds "
            "this value in cm before helix PCA."
        ),
    )
    parser.add_argument(
        "--pre-lproj-min",
        type=float,
        default=None,
        help="Cheap preselection: require rough projected decay length above this value in cm before helix PCA.",
    )
    parser.add_argument(
        "--pre-cos-theta-min",
        type=float,
        default=None,
        help="Cheap preselection: require rough pointing cosine above this value before helix PCA.",
    )
    parser.add_argument(
        "--write-chunk-size",
        type=int,
        default=50000,
        help="Number of selected pairs buffered before extending the output ROOT tree.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=500,
        help="Print progress every N input files when using --input-list. Set to 0 to disable.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.input_list:
        input_specs = read_input_list(args.input_list)
    elif args.input_file:
        input_specs = [InputSpec(args.input_file, None)]
    else:
        raise SystemExit("provide either input_file or --input-list")

    preselection = PreselectionConfig(
        track_pt_min=args.pre_track_pt_min,
        track_dca_xy_min=args.pre_track_dca_xy_min,
        track_dca_z_min=args.pre_track_dca_z_min,
        track_dca_xy_max=args.pre_track_dca_xy_max,
        track_dca_z_max=args.pre_track_dca_z_max,
        pair_dca_max=args.pre_pair_dca_max,
        lproj_min=args.pre_lproj_min,
        cos_theta_min=args.pre_cos_theta_min,
    )
    kalman_meas_sigma_xy = args.kalman_meas_sigma_xy
    kalman_meas_sigma_z = args.kalman_meas_sigma_z
    if args.kalman_measurement_preset == "truth":
        kalman_meas_sigma_xy = KALMAN_TRUTH_MEAS_SIGMA_CM
        kalman_meas_sigma_z = KALMAN_TRUTH_MEAS_SIGMA_CM
    elif args.kalman_measurement_preset == "realistic":
        kalman_meas_sigma_xy = KALMAN_REALISTIC_MEAS_SIGMA_XY_CM
        kalman_meas_sigma_z = KALMAN_REALISTIC_MEAS_SIGMA_Z_CM

    kalman_config = KalmanConfig(
        bfield_t=args.bfield,
        meas_sigma_xy_cm=kalman_meas_sigma_xy,
        meas_sigma_z_cm=kalman_meas_sigma_z,
        process_sigma_pos_cm=args.kalman_process_sigma_pos,
        process_sigma_phi=args.kalman_process_sigma_phi,
        process_sigma_qop_t=args.kalman_process_sigma_qop,
        process_sigma_tanl=args.kalman_process_sigma_tanl,
        material_x0_per_cm=args.kalman_material_x0_per_cm,
        multiple_scattering_scale=args.kalman_ms_scale,
        energy_loss_gev_per_cm=args.kalman_energy_loss_gev_per_cm,
        energy_loss_sigma_fraction=args.kalman_energy_loss_sigma_fraction,
    )

    build_pair_tree(
        input_specs=input_specs,
        output_file=args.output,
        primary_vertex=np.asarray(args.vertex, dtype=float),
        min_points=args.min_points,
        run=args.run,
        same_parent_only=args.same_parent_only,
        parent_pids=set(args.parent_pids),
        use_truth_primary_vertex=args.use_truth_primary_vertex,
        use_truth_kinematics=args.use_truth_kinematics,
        fit_helix_tracks=args.fit_helix,
        fit_kalman_tracks=args.fit_kalman,
        bfield_t=args.bfield,
        kalman_config=kalman_config,
        theta_extension=args.theta_extension,
        coarse_steps=args.coarse_steps,
        pca_search=args.pca_search,
        downstream_margin=args.downstream_margin,
        kalman_max_upstream_cm=args.kalman_max_upstream_cm,
        kalman_downstream_margin_cm=args.kalman_downstream_margin_cm,
        fit_first_points=args.fit_first_points,
        pca_candidates=args.pca_candidates,
        prefer_positive_pointing=args.prefer_positive_pointing,
        write_chunk_size=args.write_chunk_size,
        progress_every=args.progress_every,
        preselection=preselection,
        point_order=args.point_order,
    )


if __name__ == "__main__":
    main()
