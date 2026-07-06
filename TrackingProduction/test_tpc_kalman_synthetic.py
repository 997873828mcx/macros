#!/usr/bin/env python3
"""Synthetic helix checks for the lightweight TPC Kalman fitter.

The goal is to test sign conventions and PCA behavior with tracks whose truth
is exactly known.  This is deliberately independent of ROOT/Geant4 input.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np

from tpc_kalman import (
    IDX_PHI,
    IDX_QOP_T,
    IDX_TANL,
    KalmanConfig,
    fit_tpc_kalman,
    kalman_dca_to_vertex,
    kalman_track_pca_candidates,
    normalize_phi,
    propagate_state,
    state_at_path,
    state_momentum,
    state_position,
)


@dataclass
class CheckResult:
    name: str
    passed: bool
    details: str


def unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError("zero vector")
    return np.asarray(vector, dtype=float) / norm


def armenteros(pplus: np.ndarray, pminus: np.ndarray) -> tuple[float, float]:
    v0p = pplus + pminus
    direction = unit(v0p)
    pl_plus = float(np.dot(pplus, direction))
    pl_minus = float(np.dot(pminus, direction))
    denom = pl_plus + pl_minus
    if abs(denom) < 1.0e-12:
        raise ValueError("zero AP longitudinal momentum")
    alpha = (pl_plus - pl_minus) / denom
    qt = float(np.linalg.norm(pplus - direction * pl_plus))
    return alpha, qt


def truth_config(point_sorting: str = "radius") -> KalmanConfig:
    """Configuration for near-perfect synthetic space points."""

    return KalmanConfig(
        point_sorting=point_sorting,
        meas_sigma_xy_cm=1.0e-5,
        meas_sigma_z_cm=1.0e-5,
        min_measurement_sigma_cm=1.0e-8,
        initial_sigma_pos_cm=0.01,
        initial_sigma_phi=0.02,
        initial_sigma_qop_t=0.02,
        initial_sigma_tanl=0.02,
        process_sigma_pos_cm=1.0e-9,
        process_sigma_phi=1.0e-11,
        process_sigma_qop_t=1.0e-13,
        process_sigma_tanl=1.0e-13,
        material_x0_per_cm=0.0,
        energy_loss_gev_per_cm=0.0,
    )


def state_from_truth(
    position: np.ndarray,
    phi: float,
    pt_gev: float,
    charge: int,
    tanl: float,
    bfield_t: float,
) -> np.ndarray:
    """Build a state using the propagator convention in tpc_kalman.py.

    With this propagator and Bz > 0, physical charge maps to
    state qop_t = -charge / pT.
    """

    if charge not in (-1, 1):
        raise ValueError("charge must be +/-1")
    qop_t = -float(charge) / float(pt_gev)
    return np.array(
        [
            float(position[0]),
            float(position[1]),
            float(position[2]),
            normalize_phi(phi),
            qop_t,
            float(tanl),
        ],
        dtype=float,
    )


def sample_track_points(
    start_state: np.ndarray,
    bfield_t: float,
    s_values: np.ndarray,
    smear_cm: float,
    rng: np.random.Generator,
) -> np.ndarray:
    points = np.vstack(
        [state_position(propagate_state(start_state, float(s), bfield_t)) for s in s_values]
    )
    if smear_cm > 0.0:
        points += rng.normal(0.0, smear_cm, points.shape)
    return points


def prompt_track_check(charge: int, smear_cm: float, rng: np.random.Generator) -> CheckResult:
    config = truth_config("radius")
    vertex = np.array([0.0, 0.0, 0.0], dtype=float)
    pt_truth = 0.8
    tanl_truth = 0.25
    phi_truth = 0.7
    state0 = state_from_truth(vertex, phi_truth, pt_truth, charge, tanl_truth, config.bfield_t)

    s_values = np.linspace(25.0, 82.0, 28)
    points = sample_track_points(state0, config.bfield_t, s_values, smear_cm, rng)
    fit = fit_tpc_kalman(points, charge=charge, config=config, pid=211 * charge)
    if fit is None:
        return CheckResult(f"prompt q={charge:+d}", False, "fit failed")

    dca_xy, dca_z = kalman_dca_to_vertex(fit, vertex)
    first_state = fit.states_smoothed[0]
    pt_fit = 1.0 / abs(float(first_state[IDX_QOP_T]))
    expected_direction = -charge * np.sign(config.bfield_t)
    direction_ok = fit.seed.direction == expected_direction

    passed = (
        dca_xy < 0.02 + 4.0 * smear_cm
        and dca_z < 0.08 + 8.0 * smear_cm
        and abs(pt_fit - pt_truth) / pt_truth < 0.02 + 4.0 * smear_cm
        and direction_ok
    )
    return CheckResult(
        f"prompt q={charge:+d}",
        bool(passed),
        (
            f"dca_xy={dca_xy:.4g} cm, dca_z={dca_z:.4g} cm, "
            f"pt_fit={pt_fit:.5g} truth={pt_truth:.5g}, "
            f"seed_direction={fit.seed.direction:+.0f}, expected={expected_direction:+.0f}"
        ),
    )


def reversed_input_diagnostic(charge: int, smear_cm: float, rng: np.random.Generator) -> CheckResult:
    """Demonstrate why point ordering matters for the sign convention."""

    config = truth_config("input")
    vertex = np.array([0.0, 0.0, 0.0], dtype=float)
    state0 = state_from_truth(vertex, 0.4, 0.8, charge, 0.15, config.bfield_t)
    s_values = np.linspace(25.0, 82.0, 28)
    points = sample_track_points(state0, config.bfield_t, s_values, smear_cm, rng)[::-1]
    fit = fit_tpc_kalman(points, charge=charge, config=config, pid=211 * charge)
    if fit is None:
        return CheckResult(f"reversed input q={charge:+d}", False, "fit failed")

    expected_direction = -charge * np.sign(config.bfield_t)
    direction_ok = fit.seed.direction == expected_direction
    dca_xy, dca_z = kalman_dca_to_vertex(fit, vertex)
    return CheckResult(
        f"reversed input q={charge:+d}",
        bool(not direction_ok),
        (
            "diagnostic expected to expose ordering risk: "
            f"dca_xy={dca_xy:.4g} cm, dca_z={dca_z:.4g} cm, "
            f"seed_direction={fit.seed.direction:+.0f}, expected={expected_direction:+.0f}"
        ),
    )


def conversion_pair_check(smear_cm: float, rng: np.random.Generator) -> CheckResult:
    config = truth_config("radius")
    conv_vertex = np.array([20.0, 0.0, 2.0], dtype=float)
    gamma_phi = 0.03
    opening = 0.004
    tanl = 0.08
    pt_plus = 0.45
    pt_minus = 0.55

    state_plus = state_from_truth(
        conv_vertex, gamma_phi + opening, pt_plus, +1, tanl, config.bfield_t
    )
    state_minus = state_from_truth(
        conv_vertex, gamma_phi - opening, pt_minus, -1, tanl, config.bfield_t
    )

    s_values = np.linspace(5.0, 65.0, 30)
    points_plus = sample_track_points(state_plus, config.bfield_t, s_values, smear_cm, rng)
    points_minus = sample_track_points(state_minus, config.bfield_t, s_values, smear_cm, rng)

    fit_plus = fit_tpc_kalman(points_plus, charge=+1, config=config, pid=-11)
    fit_minus = fit_tpc_kalman(points_minus, charge=-1, config=config, pid=11)
    if fit_plus is None or fit_minus is None:
        return CheckResult("conversion e+e-", False, "one or both fits failed")

    candidates = kalman_track_pca_candidates(
        fit_plus,
        fit_minus,
        theta_extension=2.0,
        coarse_steps=96,
        search_mode="upstream",
        downstream_margin=0.2,
        max_candidates=32,
        max_upstream_cm=80.0,
        downstream_margin_cm=5.0,
    )
    if not candidates:
        return CheckResult("conversion e+e-", False, "PCA search failed")

    pca1, pca2, pair_dca, s1, s2 = min(candidates, key=lambda item: item[2])
    reco_vertex = 0.5 * (pca1 + pca2)
    pca_to_true = float(np.linalg.norm(reco_vertex - conv_vertex))
    pca_to_true_xy = float(np.linalg.norm((reco_vertex - conv_vertex)[:2]))
    mom_plus = state_momentum(state_at_path(fit_plus, s1))
    mom_minus = state_momentum(state_at_path(fit_minus, s2))
    alpha, qt = armenteros(mom_plus, mom_minus)

    truth_mom_plus = state_momentum(state_plus)
    truth_mom_minus = state_momentum(state_minus)
    truth_alpha, truth_qt = armenteros(truth_mom_plus, truth_mom_minus)

    passed = (
        pair_dca < 0.03 + 6.0 * smear_cm
        and pca_to_true_xy < 0.08 + 150.0 * smear_cm
        and pca_to_true < 0.12 + 200.0 * smear_cm
        and abs(qt - truth_qt) < 0.003 + 3.0 * smear_cm
    )
    return CheckResult(
        "conversion e+e-",
        bool(passed),
        (
            f"pairDCA={pair_dca:.4g} cm, pca_to_true={pca_to_true:.4g} cm, "
            f"pca_to_true_xy={pca_to_true_xy:.4g} cm, "
            f"qT={qt:.5g} truth={truth_qt:.5g}, "
            f"alpha={alpha:.5g} truth={truth_alpha:.5g}, "
            f"s_plus={s1:.3g} cm, s_minus={s2:.3g} cm"
        ),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smear-cm", type=float, default=0.0, help="Gaussian point smearing in cm")
    parser.add_argument("--seed", type=int, default=12345, help="random seed")
    parser.add_argument(
        "--skip-ordering-diagnostic",
        action="store_true",
        help="do not run the intentionally reversed-input diagnostic",
    )
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    results = [
        prompt_track_check(+1, args.smear_cm, rng),
        prompt_track_check(-1, args.smear_cm, rng),
        conversion_pair_check(args.smear_cm, rng),
    ]
    if not args.skip_ordering_diagnostic:
        results.extend(
            [
                reversed_input_diagnostic(+1, args.smear_cm, rng),
                reversed_input_diagnostic(-1, args.smear_cm, rng),
            ]
        )

    failed = 0
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(f"[{status}] {result.name}: {result.details}")
        if not result.passed:
            failed += 1

    if failed:
        print(f"[test_tpc_kalman_synthetic] {failed} check(s) failed")
        return 1

    print("[test_tpc_kalman_synthetic] all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
