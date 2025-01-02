import sys
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from qtpy.QtWidgets import (
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
    QGroupBox,
    QDoubleSpinBox,
    QMessageBox,
    QRadioButton,
)
from qtpy.QtCore import Qt
from superqt import QRangeSlider
from scipy.spatial import cKDTree

from helix_fitting import (
    fit_helix_initial,
    fit_helix_direct,
    generate_helix_points_initial,
    generate_helix_points_refined,
    create_helix_tube,
)


from analysis import (
    calculate_deltas,
    save_histograms,
    save_tree,
    apply_helix_filter,
    find_helix_points_at_radius_analytic,
    calculate_deltas_with_visualization,
)

from data_loader import load_data_from_root


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Display")
        self.resize(1600, 900)  # Increased size for better visibility

        # Initialize fitting step (1: initial fitting, 2: direct fitting)
        self.fitting_step = 1

        # Initialize pick mode
        self.pick_mode = "info"  # Default mode

        # Add a new attribute to track the helix tube color
        self.helix_tube_color = "green"  # Default color for initial fitting

        # Initialize tube_radius with a default value (e.g., 0.05 as 5%)
        self.tube_radius_percentage_initial = 0.1  # 10% of the helix radius
        self.tube_radius_percentage_second = 0.01
        self.tube_radius = 0.5
        self.track_counter = 0

        self.rphi_window = 0.5  # Default window size for rphi
        self.z_window = 1.0  # Default window size for z

        # Main container widget
        main_widget = QWidget()
        self.setCentralWidget(main_widget)

        # Horizontal layout for the main interface
        main_layout = QHBoxLayout(main_widget)
        main_layout.setContentsMargins(5, 5, 5, 5)  # Added margins for aesthetics

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

        mode_layout = QHBoxLayout()

        input_group = QGroupBox("Input")
        input_layout = QVBoxLayout()
        input_group.setLayout(input_layout)

        # --- Checkboxes for Clusters and Hits ---
        self.show_clusters = QCheckBox("Clusters")
        self.show_clusters.setChecked(True)
        self.show_clusters.stateChanged.connect(self.update_display)
        input_layout.addWidget(self.show_clusters)

        self.show_hits = QCheckBox("Hits")
        self.show_hits.setChecked(True)
        self.show_hits.stateChanged.connect(self.update_display)
        input_layout.addWidget(self.show_hits)

        mode_layout.addWidget(input_group)

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

        mode_layout.addWidget(side_group)

        # === Pick Mode Selection Area ===
        pick_mode_group = QGroupBox("Pick Mode")
        pick_mode_layout = QVBoxLayout()
        pick_mode_group.setLayout(pick_mode_layout)

        # Radio Button for Viewing Point Info
        self.radio_view_info = QRadioButton("View Point Info")
        self.radio_view_info.setChecked(True)  # Default selection
        self.radio_view_info.toggled.connect(self.on_pick_mode_changed)
        pick_mode_layout.addWidget(self.radio_view_info)

        # Radio Button for Picking Points for Helix Fitting
        self.radio_pick_helix = QRadioButton("Pick for Helix Fitting")
        self.radio_pick_helix.toggled.connect(self.on_pick_mode_changed)
        pick_mode_layout.addWidget(self.radio_pick_helix)

        mode_layout.addWidget(pick_mode_group)
        control_layout.addLayout(mode_layout)

        # === Range Selection Area ===
        range_group = QGroupBox("Axis Range Selection")
        range_layout = QHBoxLayout()
        range_group.setLayout(range_layout)

        slider_width = 200
        label_width = 30
        spinbox_width = 50

        # Define overall ranges for each axis
        overall_ranges = {"X": (-150, 150), "Y": (-150, 150), "Z": (-400, 400)}

        # Function to create axis controls
        def create_axis_controls(
            axis_label, initial_min, initial_max, overall_min, overall_max
        ):
            label = QLabel(f"{axis_label}")
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
            hbox = QHBoxLayout()
            hbox.addWidget(label)
            hbox.addWidget(min_spin)
            hbox.addWidget(max_spin)
            hbox.addWidget(slider)
            container = QWidget()
            container.setLayout(hbox)

            # Assign the slider to the corresponding instance attribute
            if axis_label == "X":
                self.x_range_slider = slider
            elif axis_label == "Y":
                self.y_range_slider = slider
            elif axis_label == "Z":
                self.z_range_slider = slider

            return container

        # Create and add axis controls horizontally
        x_axis_widget = create_axis_controls(
            "X", -100, 100, overall_ranges["X"][0], overall_ranges["X"][1]
        )
        y_axis_widget = create_axis_controls(
            "Y", -100, 100, overall_ranges["Y"][0], overall_ranges["Y"][1]
        )
        z_axis_widget = create_axis_controls(
            "Z", -400, 400, overall_ranges["Z"][0], overall_ranges["Z"][1]
        )
        range_layout.addWidget(x_axis_widget)
        range_layout.addWidget(y_axis_widget)
        range_layout.addWidget(z_axis_widget)
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

        # --- Toggle Checkbox ---
        self.toggle_filter_checkbox = QCheckBox("Show Points Outside Tube")
        self.toggle_filter_checkbox.setChecked(True)  # Default to showing all points
        self.toggle_filter_checkbox.stateChanged.connect(self.update_display)
        helix_layout.addWidget(self.toggle_filter_checkbox)

        # --- Instruction Label (Optional) ---
        self.instruction_label = QLabel(
            "Instruction: Select 3 points for initial helix fitting."
        )
        helix_layout.addWidget(self.instruction_label)

        control_layout.addWidget(helix_group)

        left_layout.addWidget(control_area)

        # 3D Visualization Area
        self.plotter_widget = QtInteractor()
        left_layout.addWidget(self.plotter_widget)

        self.hist_widget = QWidget()
        self.hist_layout = QVBoxLayout(self.hist_widget)
        self.hist_canvas = FigureCanvas(Figure(figsize=(5, 4)))
        self.hist_layout.addWidget(self.hist_canvas)
        self.hist_ax_rphi = self.hist_canvas.figure.add_subplot(121)
        self.hist_ax_z = self.hist_canvas.figure.add_subplot(122)
        self.hist_canvas.draw()

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

        sidebar_layout.addWidget(QLabel("Delta Distributions"))
        sidebar_layout.addWidget(self.hist_widget)

        splitter.addWidget(left_panel)
        splitter.addWidget(sidebar)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        # Data placeholders
        self.cluster_data = None
        self.hit_data = None
        self.cluster_polydata = None
        self.hit_polydata = None

        self.selected_points_first = (
            []
        )  # Store 3 picked points for initial helix fitting
        self.selected_points_second = (
            []
        )  # Store additional picked points for refined helix fitting
        self.helix_tube = None  # Store the helix tube (for visualization)
        self.helix_points = None  # Store helix points for distance calculations.
        self.helix_params_initial = None  # Store initial helix parameters
        self.helix_params_refined = None  # Store refined helix parameters

        self.filtered_clusters = None

        # Show axes
        self.plotter_widget.show_axes()

        # Enable initial point picking for info mode
        self.update_point_picking()

    # --- Method to Handle Mode Changes ---
    def on_pick_mode_changed(self):
        """Handle changes in the pick mode based on radio button selection."""
        if self.radio_view_info.isChecked():
            self.pick_mode = "info"
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )
            self.update_display()
        elif self.radio_pick_helix.isChecked():
            self.pick_mode = "helix"
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )
            self.update_display()

        # Update the point picking callback based on the mode
        self.update_point_picking()

    def update_point_picking(self):
        """Update the point picking callback based on the current mode."""
        # Disable and re-enable picking to ensure a fresh start
        self.plotter_widget.disable_picking()

        if self.pick_mode == "info":
            # Enable point picking for viewing info
            self.plotter_widget.enable_point_picking(
                callback=self.on_point_picked_info, show_message=True, use_picker=True
            )
        elif self.pick_mode == "helix":
            # Enable point picking for helix fitting
            self.plotter_widget.enable_point_picking(
                callback=self.on_point_picked_helix, show_message=True, use_picker=True
            )

    def on_point_picked_info(self, picked_point, picker):
        """
        Callback function for viewing point information.
        Displays info without affecting helix fitting.
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

    def on_point_picked_helix(self, picked_point, picker):
        """
        Callback function for picking points for helix fitting.
        Stores selected points and displays info.
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
        # mesh.points is a numpy array containing the coordinates
        picked_coordinates = mesh.points[point_id]

        if self.fitting_step == 1:
            # First fitting step: store in selected_points_first
            self.selected_points_first.append(picked_coordinates)
            self.instruction_label.setText(
                f"Selected {len(self.selected_points_first)}/3 points for initial helix fitting."
            )
            if len(self.selected_points_first) == 3:
                QMessageBox.information(
                    self,
                    "Helix Fit",
                    "Three points selected. Click 'Fit Helix' to perform initial helix fitting.",
                )
        elif self.fitting_step == 2:
            # Second fitting step: store in selected_points_second
            self.selected_points_second.append(picked_coordinates)
            self.instruction_label.setText(
                f"Selected {len(self.selected_points_second)} points for refined helix fitting."
            )

    def load_data(self):
        """
        Load data from a ROOT file, process it, and display it in the 3D view.
        """
        if self.cluster_data is not None:
            del self.cluster_data
        if self.hit_data is not None:
            del self.hit_data
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root)"
        )
        if not filename:
            return  # User canceled the file dialog

        try:
            # Use the data_loader module to load data
            cluster_polydata, hit_polydata = load_data_from_root(filename)

            # Store the loaded data
            self.cluster_data = cluster_polydata
            self.hit_data = hit_polydata

            # Update the visualization with the loaded data
            self.update_display()

        except FileNotFoundError as fnf_err:
            QMessageBox.critical(self, "Load Error", str(fnf_err))
        except KeyError as key_err:
            QMessageBox.critical(self, "Load Error", str(key_err))
        except ValueError as val_err:
            QMessageBox.critical(self, "Load Error", str(val_err))
        except Exception as e:
            QMessageBox.critical(
                self, "Load Error", f"An unexpected error occurred:\n{e}"
            )

    def update_display(self):
        """
        Update the 3D visualization based on loaded data and current settings.
        Includes optional filtering based on helix tube.
        """
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
            selected_sides = []

        # Define tube radius (should match the helix tube radius)
        tube_radius = 0.5  # Adjust as needed

        # Determine whether to apply tube filtering based on the toggle
        apply_tube_filter = False
        if self.helix_tube is not None and not self.toggle_filter_checkbox.isChecked():
            # If helix is present and the user does not want to show points outside the tube
            apply_tube_filter = True

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

            if filtered_points.n_points > 0:
                self.cluster_polydata = filtered_points

                if apply_tube_filter and (
                    self.helix_params_refined is not None
                    or self.helix_params_initial is not None
                ):
                    # Get current helix parameters
                    current_helix_params = (
                        self.helix_params_refined or self.helix_params_initial
                    )

                    # Apply filter to each point
                    distance_mask = np.array(
                        [
                            apply_helix_filter(
                                point,
                                self.helix_params_initial,
                                self.rphi_window,
                                self.z_window,
                            )
                            for point in filtered_points.points
                        ]
                    )
                    filtered_points = filtered_points.extract_points(distance_mask)

                if filtered_points.n_points > 0:
                    self.filtered_clusters = filtered_points
                    self.plotter_widget.add_mesh(
                        filtered_points,
                        style="points",
                        point_size=5,
                        color="red",
                    )
                    """
                    if apply_tube_filter and (
                        self.helix_params_refined or self.helix_params_initial
                    ):
                        calculate_deltas_with_visualization(
                            self.helix_params_initial,
                            self.filtered_clusters,
                            plotter=self.plotter_widget,
                        )
                    """
                else:
                    self.filtered_clusters = np.array([])

            else:
                self.filtered_clusters = np.array([])  # Empty array if no clusters pass

        else:
            self.filtered_clusters = np.array([])  # Empty array if clusters not shown
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
            """
            # Apply helix tube filtering if required
            if apply_tube_filter and self.helix_points is not None:
                # Build KDTree from helix points
                helix_tree = cKDTree(self.helix_points)

                # Query the nearest distance for each point
                distances, _ = helix_tree.query(filtered_points.points, k=1)

                # Create a mask for points within the tube radius
                distance_mask = distances <= self.tube_radius

                # Apply the distance mask
                filtered_indices = filtered_indices[distance_mask]
                filtered_points = self.hit_data.extract_points(filtered_indices)
            """

            if apply_tube_filter and (
                self.helix_params_refined is not None
                or self.helix_params_initial is not None
            ):
                # Get current helix parameters
                current_helix_params = (
                    self.helix_params_refined or self.helix_params_initial
                )

                # Apply filter to each point
                distance_mask = np.array(
                    [
                        apply_helix_filter(
                            point,
                            self.helix_params_initial,
                            self.rphi_window,
                            self.z_window,
                        )
                        for point in filtered_points.points
                    ]
                )
                # filtered_indices = filtered_indices[distance_mask]
                filtered_points = filtered_points.extract_points(distance_mask)

            if filtered_points.n_points > 0:
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
                color=self.helix_tube_color,
                opacity=0.2,
                pickable=False,  # Prevent picking on the helix
            )

        # Disable and re-enable picking to ensure a fresh start
        self.plotter_widget.disable_picking()
        # Enable point picking based on current mode
        self.update_point_picking()
        # Restore the previously saved camera position
        self.plotter_widget.camera_position = camera_position

        self.plotter_widget.show_axes()

    def initiate_helix_fit(self):
        """Triggered when the user clicks the 'Fit Helix' button."""

        if self.pick_mode != "helix":
            QMessageBox.warning(
                self,
                "Helix Fit",
                "Please switch to 'Pick for Helix Fitting' mode to fit a helix.",
            )
            return

        if self.cluster_data is None:
            QMessageBox.warning(self, "Fit Error", "No cluster data loaded.")
            return

        if self.fitting_step == 1:
            # First fitting step: initial helix fit with three points
            if len(self.selected_points_first) != 3:
                QMessageBox.warning(
                    self,
                    "Helix Fit",
                    "Please select exactly three points before fitting a helix.",
                )
                return

            try:
                helix_params_initial = fit_helix_initial(self.selected_points_first)
            except ValueError as ve:
                QMessageBox.warning(self, "Helix Fit", str(ve))
                return

            # Store initial helix parameters
            self.helix_params_initial = helix_params_initial
            self.tube_radius = (
                self.tube_radius_percentage_initial * helix_params_initial["r"]
            )

            self.helix_tube_color = "green"
            # Generate helix points and create a tube
            helix_points = generate_helix_points_initial(helix_params_initial)
            self.helix_points = (
                helix_points  # Store helix points for distance calculations
            )
            self.helix_tube = create_helix_tube(
                helix_points, tube_radius=self.tube_radius
            )
            if self.helix_tube is not None:
                self.plotter_widget.add_mesh(
                    self.helix_tube,
                    color=self.helix_tube_color,
                    opacity=0.5,
                    pickable=False,
                )
            else:
                QMessageBox.warning(
                    self, "Helix Fit", "Failed to create helix visualization."
                )

            # Clear selected points for second fitting
            self.selected_points_first.clear()
            self.fitting_step = 2  # Move to second fitting step
            self.instruction_label.setText(
                "Instruction: Select additional points within the helix tube for refined fitting."
            )
            QMessageBox.information(
                self,
                "Helix Fit",
                "Initial helix fitted. Now select additional points within the helix tube for refined fitting.",
            )
            self.update_display()

        elif self.fitting_step == 2:
            # Second fitting step: direct helix fit with additional points
            """
            if len(self.selected_points_second) < 3:
                QMessageBox.warning(
                    self,
                    "Helix Fit",
                    "Please select at least three additional points for direct helix fitting.",
                )
                return
            """
            """
            # Fetch all points currently displayed (clusters and hits)
            all_displayed_points = []
            if self.cluster_polydata is not None and self.cluster_polydata.n_points > 0:
                all_displayed_points.append(self.cluster_polydata.points)
            if self.hit_polydata is not None and self.hit_polydata.n_points > 0:
                all_displayed_points.append(self.hit_polydata.points)
            if not all_displayed_points:
                QMessageBox.warning(
                    self,
                    "Helix Fit",
                    "No points available within the helix tube for direct fitting.",
                )
                return
            all_displayed_points = np.vstack(all_displayed_points)
            """

            try:
                
                filtered_points = self.filtered_clusters.points
                #helix_params_refined = fit_helix_direct(
                 #   self.selected_points_second, self.helix_params_initial
                #)
                helix_params_refined = fit_helix_direct(
                    filtered_points, self.helix_params_initial
                )
            except ValueError as ve:
                QMessageBox.warning(self, "Helix Fit", str(ve))
                return

            if helix_params_refined is None:
                QMessageBox.warning(self, "Helix Fit", "Refined helix fitting failed.")
                return

            # Store refined helix parameters
            self.helix_params_refined = helix_params_refined

            self.tube_radius = (
                self.tube_radius_percentage_second * helix_params_refined["r"]
            )

            self.helix_tube_color = "blue"

            # Generate refined helix points and create a tube
            helix_points_refined = generate_helix_points_refined(helix_params_refined)
            self.helix_points = (
                helix_points_refined  # Update helix points for distance calculations
            )

            self.helix_tube = create_helix_tube(
                helix_points_refined, tube_radius=self.tube_radius
            )
            if self.helix_tube is not None:
                self.plotter_widget.add_mesh(
                    self.helix_tube,
                    color=self.helix_tube_color,
                    opacity=0.5,
                    pickable=False,
                )
                QMessageBox.information(
                    self,
                    "Helix Fit",
                    "Refined helix fitted and visualized.",
                )
            else:
                QMessageBox.warning(
                    self, "Helix Fit", "Failed to create refined helix visualization."
                )

            if (
                self.filtered_clusters is not None
                and self.filtered_clusters.n_points > 0
            ):
                delta_rphi, delta_z = calculate_deltas(
                    self.helix_params_refined, self.filtered_clusters
                )
            else:
                QMessageBox.warning(
                    self,
                    "Helix Fit",
                    "No clusters passed the filter. Cannot calculate deltas.",
                )
                return

            # **Calculate sigma (standard deviation)**
            sigma_rphi = np.std(delta_rphi)
            sigma_z = np.std(delta_z)
            self.track_counter += 1
            track_id = self.track_counter

            root_output = "helix_fitting_results.root"
            try:
                save_histograms(root_output, track_id, delta_rphi, delta_z)
                save_tree(
                    root_output, track_id, delta_rphi, delta_z, sigma_rphi, sigma_z
                )
            except Exception as e:
                QMessageBox.critical(
                    self, "Save Error", f"An error occurred while saving results:\n{e}"
                )
                return

            # **Update info panel with sigma values**
            info_text = (
                f"Track ID: {track_id}\n"
                f"Sigma Delta rphi: {sigma_rphi:.4f} cm\n"
                f"Sigma Delta z: {sigma_z:.4f} cm\n"
            )
            self.info_panel.setText(info_text)

            # **Plot histograms within the GUI**
            self.plot_histograms(delta_rphi, delta_z, track_id)

            # **Reset for next fitting**
            self.selected_points_second.clear()
            self.fitting_step = 1  # Reset to initial fitting step
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )
            self.update_display()

    def plot_histograms(self, delta_rphi, delta_z, track_id):
        """
        Plot histograms of delta rphi and delta z for a given track.

        Parameters:
        -----------
        delta_rphi : np.ndarray
            Array of delta rphi (arc length) values.
        delta_z : np.ndarray
            Array of delta z values.
        track_id : int
            Unique identifier for the track.
        """
        # Clear previous histograms
        self.hist_ax_rphi.clear()
        self.hist_ax_z.clear()

        # Plot Delta rphi
        self.hist_ax_rphi.hist(
            delta_rphi, bins=50, range=(-1, 1), color="red", alpha=0.7
        )
        self.hist_ax_rphi.set_title(f"Delta rphi Distribution for Track {track_id}")
        self.hist_ax_rphi.set_xlabel("Delta rphi (cm)")  # Adjust units as needed
        self.hist_ax_rphi.set_ylabel("Counts")

        # Plot Delta z
        self.hist_ax_z.hist(delta_z, bins=50, range=(-2, 2), color="blue", alpha=0.7)
        self.hist_ax_z.set_title(f"Delta z Distribution for Track {track_id}")
        self.hist_ax_z.set_xlabel("Delta z (cm)")
        self.hist_ax_z.set_ylabel("Counts")

        # Refresh the canvas
        self.hist_canvas.draw()

    def reset_helix(self):
        """Clears the fitted helix and resets the selection."""
        self.helix_tube = None
        self.helix_points = None  # Clear helix points
        self.helix_params_initial = None
        self.helix_params_refined = None
        self.selected_points_first.clear()
        self.selected_points_second.clear()
        self.fitting_step = 1  # Reset fitting step
        self.instruction_label.setText(
            "Instruction: Select 3 points for initial helix fitting."
        )
        self.filtered_clusters = None  # Clear filtered clusters
        self.update_display()
        QMessageBox.information(
            self,
            "Helix Reset",
            "Helix has been reset. You can select new points and fit again.",
        )

    # backup
    def export_results(self):
        """
        Export the ROOT file containing all track analyses.
        """
        if self.track_counter == 0:
            QMessageBox.warning(
                self, "Export Error", "No tracks have been analyzed yet."
            )
            return

        # Open a file dialog to choose save location
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Analysis Results",
            "helix_fitting_results.root",
            "ROOT Files (*.root)",
            options=options,
        )
        if filename:
            try:
                # Assuming 'helix_fitting_results.root' already contains all data
                # Copy it to the desired location
                import shutil

                shutil.copy("helix_fitting_results.root", filename)
                QMessageBox.information(
                    self, "Export Success", f"Results saved to {filename}"
                )
            except Exception as e:
                QMessageBox.critical(
                    self, "Export Error", f"Failed to export results:\n{e}"
                )

    def reset_plotter(self):
        """
        Reset the 3D plotter to its initial state.
        """
        self.plotter_widget.clear()
        self.plotter_widget.show_axes()
        self.update_display()
