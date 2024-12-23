import sys
import numpy as np
import uproot
import pyvista as pv
from pyvistaqt import QtInteractor
from qtpy.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QHBoxLayout,
    QVBoxLayout,
    QTextEdit,
    QFileDialog,
    QPushButton,
    QCheckBox,
    QLabel,
    QSplitter,
    QFrame,
    QFormLayout,
    QGroupBox,
    QDoubleSpinBox,
    QMessageBox,
)
from qtpy.QtCore import Qt
from superqt import QRangeSlider
from scipy.spatial import cKDTree


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Display")
        self.resize(1400, 800)

        # Main container widget
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # Horizontal layout for the main interface
        main_layout = QHBoxLayout(main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Use a QSplitter to allow resizing between the 3D view and the sidebar
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)

        # Left side: Vertical layout with top controls and bottom 3D view
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(5, 5, 5, 5)

        # === Control Area ===
        control_area = QWidget()
        control_layout = QVBoxLayout(control_area)
        control_layout.setContentsMargins(0, 0, 0, 0)

        # --- Data Loading Button ---
        btn_load = QPushButton("Load ROOT File")
        btn_load.clicked.connect(self.load_data)
        control_layout.addWidget(btn_load)

        # --- Checkboxes for Clusters and Hits ---
        self.show_clusters = QCheckBox("Show Clusters")
        self.show_clusters.setChecked(True)
        self.show_clusters.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_clusters)

        self.show_hits = QCheckBox("Show Hits")
        self.show_hits.setChecked(True)
        self.show_hits.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_hits)

        # --- Side Selection ---
        side_group = QGroupBox("Select Side")
        side_layout = QVBoxLayout()
        side_group.setLayout(side_layout)

        self.side0_checkbox = QCheckBox("Side 0")
        self.side0_checkbox.setChecked(True)
        self.side0_checkbox.stateChanged.connect(self.update_display)
        side_layout.addWidget(self.side0_checkbox)

        self.side1_checkbox = QCheckBox("Side 1")
        self.side1_checkbox.setChecked(True)
        self.side1_checkbox.stateChanged.connect(self.update_display)
        side_layout.addWidget(self.side1_checkbox)

        control_layout.addWidget(side_group)

        # === Range Selection Area ===
        range_group = QGroupBox("Axis Range Selection")
        range_layout = QFormLayout()
        range_group.setLayout(range_layout)

        slider_width = 300
        label_width = 80
        spinbox_width = 60

        # Define overall ranges for each axis
        overall_ranges = {"X": (-150, 150), "Y": (-150, 150), "Z": (-400, 400)}

        # Function to create axis controls
        def create_axis_controls(
            axis_label, initial_min, initial_max, overall_min, overall_max
        ):
            label = QLabel(f"{axis_label} Range:")
            label.setFixedWidth(label_width)

            min_spin = QDoubleSpinBox()
            min_spin.setRange(overall_min, overall_max)
            min_spin.setDecimals(0)
            min_spin.setValue(initial_min)
            min_spin.setFixedWidth(spinbox_width)

            max_spin = QDoubleSpinBox()
            max_spin.setRange(overall_min, overall_max)
            max_spin.setDecimals(0)
            max_spin.setValue(initial_max)
            max_spin.setFixedWidth(spinbox_width)

            slider = QRangeSlider(Qt.Horizontal)
            slider.setMinimum(overall_min)
            slider.setMaximum(overall_max)
            slider.setValue([initial_min, initial_max])
            slider.setFixedWidth(slider_width)

            # Connect spin boxes to slider
            def on_min_spin_change(value):
                current_max = slider.value()[1]
                if value > current_max:
                    min_spin.blockSignals(True)
                    min_spin.setValue(current_max)
                    min_spin.blockSignals(False)
                    value = current_max
                slider.setValue([value, current_max])  # Only update handle positions
                self.update_display()

            def on_max_spin_change(value):
                current_min = slider.value()[0]
                if value < current_min:
                    max_spin.blockSignals(True)
                    max_spin.setValue(current_min)
                    max_spin.blockSignals(False)
                    value = current_min
                slider.setValue([current_min, value])  # Only update handle positions
                self.update_display()

            min_spin.valueChanged.connect(on_min_spin_change)
            max_spin.valueChanged.connect(on_max_spin_change)

            # Connect slider to spin boxes
            def on_slider_change(values):
                min_val, max_val = values
                min_spin.blockSignals(True)
                max_spin.blockSignals(True)
                min_spin.setValue(min_val)
                max_spin.setValue(max_val)
                min_spin.blockSignals(False)
                max_spin.blockSignals(False)
                self.update_display()

            slider.valueChanged.connect(on_slider_change)

            return label, min_spin, max_spin, slider

        # X Axis
        x_initial_min, x_initial_max = -100, 100
        x_label, self.x_min_spin, self.x_max_spin, self.x_range_slider = (
            create_axis_controls(
                "X",
                x_initial_min,
                x_initial_max,
                overall_ranges["X"][0],
                overall_ranges["X"][1],
            )
        )
        x_hbox = QHBoxLayout()
        x_hbox.addWidget(self.x_min_spin)
        x_hbox.addWidget(self.x_max_spin)
        x_hbox.addWidget(self.x_range_slider)
        range_layout.addRow(x_label, x_hbox)

        # Y Axis
        y_initial_min, y_initial_max = -100, 100
        y_label, self.y_min_spin, self.y_max_spin, self.y_range_slider = (
            create_axis_controls(
                "Y",
                y_initial_min,
                y_initial_max,
                overall_ranges["Y"][0],
                overall_ranges["Y"][1],
            )
        )
        y_hbox = QHBoxLayout()
        y_hbox.addWidget(self.y_min_spin)
        y_hbox.addWidget(self.y_max_spin)
        y_hbox.addWidget(self.y_range_slider)
        range_layout.addRow(y_label, y_hbox)

        # Z Axis
        z_initial_min, z_initial_max = -400, 400
        z_label, self.z_min_spin, self.z_max_spin, self.z_range_slider = (
            create_axis_controls(
                "Z",
                z_initial_min,
                z_initial_max,
                overall_ranges["Z"][0],
                overall_ranges["Z"][1],
            )
        )
        z_hbox = QHBoxLayout()
        z_hbox.addWidget(self.z_min_spin)
        z_hbox.addWidget(self.z_max_spin)
        z_hbox.addWidget(self.z_range_slider)
        range_layout.addRow(z_label, z_hbox)

        control_layout.addWidget(range_group)

        # === Helix Fitting Controls ===
        helix_group = QGroupBox("Helix Fitting")
        helix_layout = QHBoxLayout()
        helix_group.setLayout(helix_layout)

        # Button: Fit Helix
        self.btn_fit_helix = QPushButton("Fit Helix")
        self.btn_fit_helix.clicked.connect(self.initiate_helix_fit)
        helix_layout.addWidget(self.btn_fit_helix)

        # Button: Reset Helix
        self.btn_reset_helix = QPushButton("Reset Helix")
        self.btn_reset_helix.clicked.connect(self.reset_helix)
        helix_layout.addWidget(self.btn_reset_helix)

        control_layout.addWidget(helix_group)

        left_layout.addWidget(control_area)

        # 3D Visualization Area
        self.plotter_widget = QtInteractor()
        left_layout.addWidget(self.plotter_widget)

        # A frame line at the bottom for neatness (optional)
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        left_layout.addWidget(line)

        # Sidebar for info
        sidebar = QWidget()
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(5, 5, 5, 5)

        info_label = QLabel("Cluster/Hit Information")
        sidebar_layout.addWidget(info_label)

        self.info_panel = QTextEdit()
        self.info_panel.setReadOnly(True)
        sidebar_layout.addWidget(self.info_panel)

        splitter.addWidget(left_panel)
        splitter.addWidget(sidebar)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        # Data placeholders
        self.cluster_data = None
        self.hit_data = None
        self.cluster_polydata = None
        self.hit_polydata = None

        self.selected_points = []  # Store 3 picked points
        self.helix_tube = None  # Store the helix tube (for visualization)

        # Show axes
        self.plotter_widget.show_axes()

    def load_data(self):

        filename, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root)"
        )
        if not filename:
            return

        try:

            file = uproot.open(filename)

            cluster_tree = file["combined_clusters"]
            cluster_data = cluster_tree.arrays(library="np")

            if "side" not in cluster_data:
                raise ValueError("The 'side' branch is missing in the cluster data.")

            cx = cluster_data["gx"]
            cy = cluster_data["gy"]
            cz = cluster_data["gz"]
            cpoints = np.column_stack([cx, cy, cz])
            cluster_polydata = pv.PolyData(cpoints)
            for name in cluster_data.keys():
                if name not in ("gx", "gy", "gz"):
                    cluster_polydata.point_data[name] = cluster_data[name]

            cluster_polydata.point_data["data_type"] = np.zeros(
                cluster_polydata.n_points, dtype=int
            )
            self.cluster_data = cluster_polydata

            # Load hits tree
            hits_tree = file["combined_hits"]
            hits_data = hits_tree.arrays(library="np")

            # Ensure 'side' branch exists
            if "side" not in hits_data:
                raise ValueError("The 'side' branch is missing in the hit data.")

            hx = hits_data["gx"]
            hy = hits_data["gy"]
            hz = hits_data["gz"]
            hpoints = np.column_stack([hx, hy, hz])
            hit_polydata = pv.PolyData(hpoints)

            # Add other attributes
            for name in hits_data.keys():
                if name not in ("gx", "gy", "gz"):
                    hit_polydata.point_data[name] = hits_data[name]

            # Add a dataset type attribute (1 for hits)
            hit_polydata.point_data["data_type"] = np.ones(
                hit_polydata.n_points, dtype=int
            )
            self.hit_data = hit_polydata
            self.update_display()

        except Exception as e:
            self.info_panel.setText(f"Error loading file:\n{e}")

    def update_display(self):
        # Save current camera position
        camera_position = self.plotter_widget.camera_position
        # Clear the current plotter
        self.plotter_widget.clear()

        # Retrieve range values from sliders
        x_min, x_max = self.x_range_slider.value()
        y_min, y_max = self.y_range_slider.value()
        z_min, z_max = self.z_range_slider.value()

        # Determine what to show
        show_clusters = self.show_clusters.isChecked()
        show_hits = self.show_hits.isChecked()

        # Retrieve selected sides
        selected_sides = []
        if self.side0_checkbox.isChecked():
            selected_sides.append(0)
        if self.side1_checkbox.isChecked():
            selected_sides.append(1)

        # Handle case when no sides are selected
        if not selected_sides:
            # Optionally, you can decide to show nothing or all sides
            # Here, we'll show nothing
            selected_sides = []

        # Keep track of polydata items added to the scene
        # so we know which is which when picking
        self.cluster_polydata = None
        self.hit_polydata = None

        # --- Clusters ---
        if (
            show_clusters
            and self.cluster_data is not None
            and self.cluster_data.n_points > 0
        ):

            # Apply the slider range to filter points
            points = self.cluster_data.points  # Get cluster points
            side = self.cluster_data.point_data["side"]  # Get side data

            if selected_sides:
                side_mask = np.isin(side, selected_sides)
            else:
                side_mask = False  # No sides selected, no points

            mask = (
                (points[:, 0] >= x_min)
                & (points[:, 0] <= x_max)  # X range
                & (points[:, 1] >= y_min)
                & (points[:, 1] <= y_max)  # Y range
                & (points[:, 2] >= z_min)
                & (points[:, 2] <= z_max)  # Z range
                & side_mask  # Side filter
            )
            # Extract filtered points
            filtered_indices = np.where(mask)[0]
            filtered_points = self.cluster_data.extract_points(filtered_indices)

            if (
                filtered_points.n_points > 0
            ):  # Check if any points are left after filtering
                self.cluster_polydata = filtered_points
                # Add cluster points in one color, e.g., red
                self.plotter_widget.add_mesh(
                    self.cluster_polydata,
                    # render_points_as_spheres=True,
                    style="points",
                    point_size=5,
                    color="red",
                )
        # --- Hits ---
        if show_hits and self.hit_data is not None and self.hit_data.n_points > 0:

            # Apply the slider range to filter points
            points = self.hit_data.points  # Get hit points
            side = self.hit_data.point_data["side"]  # Get side data

            if selected_sides:
                side_mask = np.isin(side, selected_sides)
            else:
                side_mask = False  # No sides selected, no points

            mask = (
                (points[:, 0] >= x_min)
                & (points[:, 0] <= x_max)  # X range
                & (points[:, 1] >= y_min)
                & (points[:, 1] <= y_max)  # Y range
                & (points[:, 2] >= z_min)
                & (points[:, 2] <= z_max)  # Z range
                & side_mask  # Side filter
            )

            # Extract filtered points
            filtered_indices = np.where(mask)[0]
            filtered_points = self.hit_data.extract_points(filtered_indices)

            if (
                filtered_points.n_points > 0
            ):  # Check if any points are left after filtering
                self.hit_polydata = filtered_points
                # Add hit points in another color, e.g., blue
                self.plotter_widget.add_mesh(
                    self.hit_polydata,
                    style="points",
                    point_size=5,
                    color="blue",
                )

        # If a helix tube is already generated, re-add it (ensure it's not pickable)
        if self.helix_tube is not None:
            self.plotter_widget.add_mesh(
                self.helix_tube,
                color="green",
                opacity=0.5,
                pickable=False,  # Prevent picking on the helix
            )
        # Disable and re-enable picking to ensure a fresh start
        self.plotter_widget.disable_picking()
        # Still using use_mesh=True for now
        self.plotter_widget.enable_point_picking(
            callback=self.on_point_picked, show_message=True, use_picker=True
        )
        # Restore the previously saved camera position
        self.plotter_widget.camera_position = camera_position

        self.plotter_widget.show_axes()

    def on_point_picked(self, picked_point, picker):
        """
        Callback function triggered when a point is picked.

        Parameters:
        - picked_point: The coordinates of the picked point.
        - picker: The vtkPointPicker object containing information about the pick.
        """
        # Extract the point ID from the picker
        point_id = picker.GetPointId()

        # Extract the mesh (dataset) from the picker
        mesh = picker.GetDataSet()

        # Validate the picked point
        if point_id < 0:
            return  # No valid point was picked

        if mesh is None:
            return  # No mesh was picked

        # Ensure 'data_type' exists in the mesh's point data
        if "data_type" not in mesh.point_data:
            return  # Ignore picking on meshes without 'data_type'

        # Retrieve the data_type for the picked point
        data_type = mesh.point_data["data_type"][point_id]
        label = "Cluster" if data_type == 0 else "Hit"

        # Display information about the picked point
        info_text = f"{label} (Point ID: {point_id})\n"
        for attr in mesh.point_data.keys():
            value = mesh.point_data[attr][point_id]
            info_text += f"{attr}: {value}\n"

        self.info_panel.setText(info_text)

        # Store the picked point for helix fitting
        # Ensure that 'points' is accessible; if mesh is a PolyData, use mesh.points
        if isinstance(mesh, pv.PolyData) or isinstance(mesh, pv.UnstructuredGrid):
            picked_coordinates = mesh.points[point_id]
        else:
            # If mesh is another type, adjust accordingly
            picked_coordinates = np.array(picked_point)

        self.selected_points.append(picked_coordinates)

        # Inform the user if three points have been selected
        if len(self.selected_points) == 3:
            QMessageBox.information(
                self,
                "Helix Fit",
                "Three points selected. Click 'Fit Helix' to proceed.",
            )

    def initiate_helix_fit(self):
        """Triggered when the user clicks the 'Fit Helix' button."""
        if len(self.selected_points) != 3:
            QMessageBox.warning(
                self,
                "Helix Fit",
                "Please select exactly three points before fitting a helix.",
            )
            return

        helix_params = self.fit_helix(self.selected_points)
        if helix_params is None:
            QMessageBox.warning(self, "Helix Fit", "Helix fitting failed.")
        else:
            # Generate helix points and create a tube
            helix_points = self.generate_helix_points(helix_params)
            self.helix_tube = self.create_helix_tube(helix_points)
            if self.helix_tube is not None:
                self.plotter_widget.add_mesh(
                    self.helix_tube, color="green", opacity=0.5, pickable=False
                )
                QMessageBox.information(
                    self, "Helix Fit", "Helix fitted and displayed successfully."
                )
            else:
                QMessageBox.warning(
                    self, "Helix Fit", "Failed to create helix visualization."
                )

        # Clear selected points for next usage
        self.selected_points.clear()
        self.update_display()

    def reset_helix(self):
        """Clears the fitted helix and resets the selection."""
        self.helix_tube = None
        self.selected_points.clear()
        self.update_display()
        QMessageBox.information(
            self,
            "Helix Reset",
            "Helix has been reset. You can select new points and fit again.",
        )

    def fit_helix(self, points):
        """
        Fit a helix to exactly three points.
        Returns [radius, pitch, t0, center_x, center_y, center_z] or None on failure.
        """
        if len(points) != 3:
            return None

        p1, p2, p3 = points

        # --- 1. Project onto XY-plane ---
        p1_xy = p1[:2]
        p2_xy = p2[:2]
        p3_xy = p3[:2]

        # --- 2. Fit circle in XY-plane to find center and radius ---
        center_xy, radius = self.fit_circle_2d(p1_xy, p2_xy, p3_xy)
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

        # --- 4. Check if points span multiple turns ---
        if thetas[-1] - thetas[0] > 2 * np.pi:
            QMessageBox.warning(
                self,
                "Helix Fit",
                "Selected points span multiple helical turns. Please select three points within the same turn.",
            )
            return None

        # --- 5. Linear fit between theta and z to find pitch ---
        z_vals = np.array([p1[2], p2[2], p3[2]])
        A = np.vstack([thetas, np.ones_like(thetas)]).T
        try:
            slope, z0 = np.linalg.lstsq(A, z_vals, rcond=None)[0]
        except Exception:
            return None

        pitch = slope * (2 * np.pi)  # pitch for a 2pi rotation
        t0 = thetas[0]  # phase offset
        center_z = z0

        return [radius, pitch, t0, center_xy[0], center_xy[1], center_z]

    def fit_circle_2d(self, p1, p2, p3):
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

    def generate_helix_points(self, params, num_points=500):
        """
        From helix parameters [r, pitch, t0, cx, cy, cz], generate helix points.
        """
        r, pitch, t0, cx, cy, cz = params
        # Generate 2 full turns for visualization
        t = np.linspace(t0, t0 + 4 * np.pi, num_points)
        x = cx + r * np.cos(t)
        y = cy + r * np.sin(t)
        z = cz + (pitch / (2 * np.pi)) * t
        return np.column_stack((x, y, z))

    def create_helix_tube(self, helix_points, tube_radius=0.1):
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

    # --------------------------- END HELIX FITTING FUNCTIONS ---------------------------


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
