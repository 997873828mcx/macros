from qtpy.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QComboBox, QHBoxLayout
from pyvistaqt import QtInteractor
import pyvista as pv
import numpy as np
from typing import Optional


class GeometricProjectionWindow(QMainWindow):
    """
    A separate window for 2D geometric projections (XY, XZ, YZ) of event display data.
    Allows plotting clusters, hits, and helix lines on a chosen plane
    and overlays TPC boundaries for context.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Geometric Projection View")
        self.resize(800, 800)

        # Create main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)

        # ---------------------- Projection Controls ----------------------
        control_layout = QHBoxLayout()

        self.projection_combo = QComboBox()
        self.projection_combo.addItems(["XY Projection", "ZR Projection"])
        self.projection_combo.currentIndexChanged.connect(self.update_projection)

        control_layout.addWidget(self.projection_combo)
        layout.addLayout(control_layout)

        # ---------------------- PyVista Plotter ----------------------
        self.plotter = QtInteractor()
        layout.addWidget(self.plotter)

        # ---------------------- Data Placeholders ----------------------
        self.current_cluster_data = None  # pv.PolyData for clusters
        self.current_hit_data = None  # pv.PolyData for hits
        self.current_helix_lines = None  # List of VTK actors for helix lines

    # ------------------------------------------------------------------
    #                    Projection / Display Logic
    # ------------------------------------------------------------------

    def project_points(self, points: np.ndarray, projection_type: str) -> np.ndarray:
        """
        Perform geometric projection of points onto a specified plane.

        Parameters
        ----------
        points : np.ndarray
            (N, 3) array of 3D points
        projection_type : str
            One of 'XY Projection', 'XZ Projection', or 'YZ Projection'

        Returns
        -------
        np.ndarray
            (N, 3) array of projected points
        """
        print(f"DEBUG: project_points called with projection_type='{projection_type}'")
        if points is None or len(points) == 0:
            return np.array([])

        projected_points = points.copy()

        if projection_type == "XY Projection":
            # Set Z to 0
            projected_points[:, 2] = 0
        elif projection_type == "ZR Projection":
            # Set Y to 0
            r = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)

            projected_points[:, 0] = r  # R coordinate
            projected_points[:, 1] = points[:, 2]  # Z coordinate
            projected_points[:, 2] = 0  # Set Z to 0 for 2D projection
        print(f"Projection Type: {projection_type}")
        print(f"Original Points Sample:\n{points[:5]}")
        print(f"Projected Points Sample:\n{projected_points[:5]}")

        return projected_points

    def project_polydata(
        self, polydata: pv.DataSet, projection_type: str
    ) -> Optional[pv.PolyData]:
        """
        Project PyVista data (PolyData or similar) onto specified plane.
        """
        if polydata is None or polydata.n_points == 0:
            return None

        # Project points
        projected_points = self.project_points(polydata.points, projection_type)
        if projected_points.size == 0:
            return None

        # Create new PolyData with projected points
        projected_polydata = pv.PolyData(projected_points)

        # Copy over point data from original
        for key in polydata.point_data.keys():
            projected_polydata.point_data[key] = polydata.point_data[key]

        # Copy lines if available and if the input is PolyData
        if isinstance(polydata, pv.PolyData) and polydata.lines is not None:
            projected_polydata.lines = polydata.lines

        return projected_polydata

    def add_projected_circles(self, projection_type: str):
        """
        Add TPC boundary circles or lines representing
        the inner/outer TPC boundaries, depending on projection.
        """
        if projection_type == "XY Projection":
            # For XY, draw actual circles in the XY plane
            theta = np.linspace(0, 2 * np.pi, 100)
            for radius in (21.6, 76.4):  # Inner and outer TPC radii
                x = radius * np.cos(theta)
                y = radius * np.sin(theta)
                z = np.zeros_like(theta)

                circle_points = np.column_stack((x, y, z))
                circle_poly = pv.PolyData(circle_points)

                # Build a closed polyline: need n+1 indices for n points + repeat(0)
                n_points = len(theta)
                line_indices = np.hstack(
                    (
                        n_points + 1,
                        np.arange(n_points),
                        0,  # repeat first point to close
                    )
                )
                circle_poly.lines = line_indices

                self.plotter.add_mesh(
                    circle_poly, color="gray", line_width=1, style="wireframe"
                )

        elif projection_type == "ZR Projection":
            # For XZ, the TPC is a 'vertical' range in Z and a horizontal range in X
            z_range = [-105.5, 105.5]  # half-length of TPC
            # We'll draw lines at x=±21.6 and x=±76.4
            for r_val in (21.6, 76.4):
                line_poly = pv.Line([r_val, z_range[0], 0], [r_val, z_range[1], 0])
                self.plotter.add_mesh(line_poly, color="gray", line_width=1)

    def update_projection(self):
        """
        Called when the user changes projection type in the combo box.
        Re-renders the current data (clusters, hits, lines) to the newly selected projection.
        """
        projection_type = self.projection_combo.currentText()
        print(f"DEBUG: combo selection => {projection_type!r}")

        self.update_display(
            self.current_cluster_data,
            self.current_hit_data,
            self.current_helix_lines,
            projection_type,
        )

    def update_display(
        self,
        clusters_info: Optional[list] = None,
        hits_info: Optional[list] = None,
        helix_lines: Optional[list] = None,
        projection_type: Optional[str] = None,
    ):
        """
        Re-draw the projection window with the given data (clusters, hits, helix lines)
        using the specified plane (XY, XZ, YZ).

        Parameters
        ----------
        cluster_data : pv.PolyData or None
            3D cluster data to project.
        hit_data : pv.PolyData or None
            3D hit data to project.
        helix_lines : list
            List of VTK actors representing the helix lines.
        projection_type : str
            One of 'XY Projection', 'XZ Projection', 'YZ Projection'.
        """
        # Store references for later refresh
        self.current_cluster_data = clusters_info
        self.current_hit_data = hits_info
        self.current_helix_lines = helix_lines

        if projection_type is None:
            projection_type = self.projection_combo.currentText()

        # Clear current view
        self.plotter.clear()

        # ------------- Project Clusters -------------

        if clusters_info:
            for filtered_clusters, color, file_info in clusters_info:
                if filtered_clusters is not None and filtered_clusters.n_points > 0:
                    projected_clusters = self.project_polydata(
                        filtered_clusters, projection_type
                    )
                    if projected_clusters:
                        self.plotter.add_mesh(
                            projected_clusters,
                            style="points",
                            point_size=5,
                            color=color,
                        )

        # ------------- Project Hits -------------
        if hits_info:
            for filtered_hits, color, file_info in hits_info:
                if filtered_hits is not None and filtered_hits.n_points > 0:
                    projected_hits = self.project_polydata(
                        filtered_hits, projection_type
                    )
                    if projected_hits:
                        self.plotter.add_mesh(
                            projected_hits, style="points", point_size=5, color="blue"
                        )

        # ------------- Project Helix Lines -------------
        if helix_lines:
            for actor in helix_lines:
                # actor: a VTK actor (pv.Actor in newer versions) whose mapper data is a PolyData
                polydata = actor.GetMapper().GetInput()
                if polydata and isinstance(polydata, pv.PolyData):
                    projected_line = self.project_polydata(polydata, projection_type)
                    if projected_line:
                        self.plotter.add_mesh(
                            projected_line,
                            color="grey",
                            line_width=3,
                            style="wireframe",
                        )

        # ------------- TPC Boundaries (circles or lines) -------------
        self.add_projected_circles(projection_type)

        # ------------- Camera Setup -------------
        self.setup_camera(projection_type)
        self.plotter.reset_camera()

    def setup_camera(self, projection_type: str):
        """
        Configure a simple orthographic camera for the chosen projection plane.
        You can adjust or remove disable_mouse_movements() if you want to allow
        user panning or zooming in the 2D view.
        """
        # For convenience, start from the built-in 'xy' view preset
        self.plotter.camera_position = "xy"
        # Force orthographic view
        self.plotter.enable_parallel_projection()

        # Minimal transformations for each plane
        if projection_type == "XY Projection":
            # Look down Z-axis
            self.plotter.camera.elevation = 0
            self.plotter.camera.azimuth = 0
            self.plotter.camera.view_up = [0, 1, 0]
        elif projection_type == "ZR Projection":
            # Look down Y-axis
            self.plotter.camera.elevation = 0
            self.plotter.camera.azimuth = 0
            # self.plotter.camera.roll = 0
            self.plotter.camera.view_up = [0, 1, 0]
            """y_distance = 200
            self.plotter.camera_position = [(0, y_distance, 0), (0, 0, 0), (0, 0, 1)]
            self.plotter.camera.view_up = [0, 0, 1]"""

        # Disable all mouse interactions (optional).
        # You can remove or replace this with plotter.disable_rotation_style()
        # if you only want to forbid rotation but allow panning/zooming.
        # self.plotter.disable_mouse_movements()
        self.plotter.reset_camera()
