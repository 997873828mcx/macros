# tests/test_analysis.py

import numpy as np
from analysis import calculate_deltas
import os
from ROOT import TFile


def test_calculate_deltas():
    helix_params = {
        "c_x": 0.0,
        "c_y": 0.0,
        "c_z": 0.0,
        "r": 1.0,
        "alpha": 1.0,
        "ref_theta": 0.0,
    }
    clusters = np.array(
        [
            [1.0, 0.0, 0.0],  # On helix
            [0.0, 1.0, 1.0],  # On helix
            [-1.0, 0.0, 2.0],  # On helix
            [1.0, 0.0, 0.1],  # Slightly off helix
        ]
    )
    delta_rphi, delta_z = calculate_deltas(helix_params, clusters)

    # Expected arc length for on-helix points
    # Since delta_phi = 0 for on-helix points, delta_rphi_val = 0
    assert np.allclose(
        delta_rphi[:3], 0.0, atol=1e-6
    ), "On-helix points should have delta_rphi = 0."

    # For the slightly off helix point
    # Assuming a small angular difference, delta_rphi_val should be small
    assert len(delta_rphi) == 4, "Should have 4 delta_rphi values."
    assert len(delta_z) == 4, "Should have 4 delta_z values."
    assert (
        delta_z[3] == 0.1
    ), "delta_z_val should be correct for the slightly off helix point."
    print("calculate_deltas passed.")


if __name__ == "__main__":
    test_calculate_deltas()
    print("All analysis.py tests passed.")
