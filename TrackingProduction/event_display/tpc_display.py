"""Minimal PyVista scene that shows the sPHENIX TPC envelopes.

This is extracted from the full event-display GUI so collaborators can plug in
new functionality without the Qt machinery.  Running this module directly
launches an interactive PyVista window with the inner/outer TPC cylinders and
the beam axis.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import pyvista as pv


@dataclass(frozen=True)
class TPCGeometry:
    """Convenient container for the key TPC dimensions (in cm)."""

    inner_radius: float = 21.6
    outer_radius: float = 76.4
    full_length: float = 211.0  # total length (±105.5 cm along z)
    beam_extent: float = 200.0  # how far to draw the beam axis line

    @property
    def half_length(self) -> float:
        return self.full_length / 2.0


def add_tpc_geometry(plotter: pv.Plotter, geometry: TPCGeometry) -> None:
    """Create and add the TPC shells + beam axis to the supplied plotter."""

    inner_cylinder = pv.Cylinder(
        center=(0.0, 0.0, 0.0),
        direction=(0.0, 0.0, 1.0),
        radius=geometry.inner_radius,
        height=geometry.full_length,
    )
    outer_cylinder = pv.Cylinder(
        center=(0.0, 0.0, 0.0),
        direction=(0.0, 0.0, 1.0),
        radius=geometry.outer_radius,
        height=geometry.full_length,
    )
    beam_line = pv.Line(
        (0.0, 0.0, -geometry.beam_extent),
        (0.0, 0.0, geometry.beam_extent),
    )

    plotter.add_mesh(
        inner_cylinder,
        color="gray",
        opacity=0.15,
        pickable=False,
        line_width=1,
        style="wireframe",
    )
    plotter.add_mesh(
        outer_cylinder,
        color="gray",
        opacity=0.25,
        pickable=False,
        line_width=1,
        style="wireframe",
    )
    plotter.add_mesh(
        beam_line,
        color="orange",
        line_width=4,
        pickable=False,
    )
    plotter.show_axes()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render the sPHENIX TPC envelopes with PyVista.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--view",
        default="iso",
        choices=["iso", "xy", "xz", "yz"],
        help="Initial camera orientation.",
    )
    return parser


def set_camera(plotter: pv.Plotter, view: str) -> None:
    """Snap the camera to one of a few handy presets."""
    if view == "iso":
        plotter.camera_position = "iso"
    elif view == "xy":
        plotter.camera_position = "xy"
    elif view == "xz":
        plotter.camera_position = "xz"
    elif view == "yz":
        plotter.camera_position = "yz"


def main() -> None:
    args = build_parser().parse_args()
    geometry = TPCGeometry()

    plotter = pv.Plotter(window_size=[900, 700])
    plotter.set_background("black")
    add_tpc_geometry(plotter, geometry)
    set_camera(plotter, args.view)
    plotter.add_text(
        "sPHENIX TPC",
        position="upper_left",
        font_size=14,
        color="white",
    )
    plotter.show(title="TPC visualization")


if __name__ == "__main__":
    main()
