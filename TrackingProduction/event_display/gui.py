import sys
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from typing import Optional
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
    QComboBox,
)
from qtpy.QtCore import Qt
from superqt import QRangeSlider
from scipy.spatial import cKDTree
from helix_fitting import (
    fit_helix_initial,
    fit_helix_direct,
    generate_helix_points_initial,
    generate_helix_points_refined,
    generate_helix_line,
)



from histogram_window import HistogramWindow
from module_histogram_window import ModuleHistogramWindow
from analysis import (
    calculate_deltas,
    apply_helix_filter,
    apply_line_filter,
    calculate_deltas_line,
    compute_centroid,
)
from data_loader import load_data_from_root
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtCore import Qt


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Display")
        self.resize(1600, 900)  # Increased size for better visibility
        self.file_colors = [
            "red",
            "green",
            "black",
            "cyan",
            "magenta",
            "purple",
        ]
        self.color_index = 0

        self.loaded_files = {}
        # Initialize fitting step (1: initial fitting, 2: direct fitting)
        self.fitting_step = 1
        self.using_line_fitting = True

        # Initialize pick mode
        self.pick_mode = "info"  # Default mode

        self.inner_cut = 21.6
        self.outer_cut = 76.4
        # Add a new attribute to track the helix tube color

        self.helix_line = None
        self.track_lines = []  # List to store all helix lines

        self.helix_line_color = "grey"  # Default color for initial fitting
        # self.helix_colors = []
        # self.available_colors = [
        #   "green", "blue", "red", "yellow", "cyan", "magenta",
        #  "orange", "purple", "brown", "pink"
        # ]
        # self.current_color_index = 0

        self.track_counter = 0

        self.rphi_window = 0.5  # Default window size for rphi
        self.z_window = 1.0  # Default window size for z

        self.histogram_window = HistogramWindow()
        self.module_hist_window = ModuleHistogramWindow()

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

        # --- ADC Threshold Controls ---
        adc_group = QGroupBox("ADC Thresholds")
        adc_layout = QVBoxLayout()
        adc_group.setLayout(adc_layout)

        # Cluster ADC threshold control
        cluster_adc_layout = QHBoxLayout()
        cluster_adc_label = QLabel("Cluster ADC >")
        self.cluster_adc_spinbox = QDoubleSpinBox()
        self.cluster_adc_spinbox.setRange(0, 10000)  # Adjust range as needed
        self.cluster_adc_spinbox.setValue(0)  # Default threshold
        self.cluster_adc_spinbox.setDecimals(0)
        self.cluster_adc_spinbox.valueChanged.connect(self.update_display)
        cluster_adc_layout.addWidget(cluster_adc_label)
        cluster_adc_layout.addWidget(self.cluster_adc_spinbox)
        adc_layout.addLayout(cluster_adc_layout)

        # Hit ADC threshold control
        hit_adc_layout = QHBoxLayout()
        hit_adc_label = QLabel("Hit ADC >")
        self.hit_adc_spinbox = QDoubleSpinBox()
        self.hit_adc_spinbox.setRange(0, 10000)  # Adjust range as needed
        self.hit_adc_spinbox.setValue(0)  # Default threshold
        self.hit_adc_spinbox.setDecimals(0)
        self.hit_adc_spinbox.valueChanged.connect(self.update_display)
        hit_adc_layout.addWidget(hit_adc_label)
        hit_adc_layout.addWidget(self.hit_adc_spinbox)
        adc_layout.addLayout(hit_adc_layout)

        mode_layout.addWidget(adc_group)

        # --- Cluster Filter Controls ---
        cluster_filter_group = QGroupBox("Seed and Track")
        cluster_filter_layout = QVBoxLayout()
        cluster_filter_group.setLayout(cluster_filter_layout)

        # Checkbox for used_in_seed
        self.seed_checkbox = QCheckBox("Seed")
        self.seed_checkbox.setChecked(False)  # Default unchecked
        self.seed_checkbox.stateChanged.connect(self.update_display)
        cluster_filter_layout.addWidget(self.seed_checkbox)

        # Checkbox for used_in_track
        self.track_checkbox = QCheckBox("Track")
        self.track_checkbox.setChecked(False)  # Default unchecked
        self.track_checkbox.stateChanged.connect(self.update_display)
        cluster_filter_layout.addWidget(self.track_checkbox)

        mode_layout.addWidget(cluster_filter_group)
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

        # --- File Selection Dropdown ---
        files_group = QGroupBox("Files")
        files_layout = QVBoxLayout()
        files_group.setLayout(files_layout)

        self.file_combo = QComboBox()
        # Set up a model to allow checkable items in the combo box
        model = QStandardItemModel(self.file_combo)
        self.file_combo.setModel(model)
        # Connect item changes to update_display so changes refresh visualization
        model.itemChanged.connect(lambda item: self.update_display())

        files_layout.addWidget(self.file_combo)
        mode_layout.addWidget(files_group)

        event_info_label = QLabel(
            "<b>Run:</b> 53217<br>"
            "<b>ZDC coincidence:</b> Raw: 1,869750<br>"
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
            "Live: 1,868219"
        )
        event_info_label.setWordWrap(True)
        mode_layout.addWidget(event_info_label)

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
        self.btn_fit_helix.clicked.connect(self.initiate_fit)
        helix_layout.addWidget(self.btn_fit_helix)

        # Button: Reset Helix
        self.btn_reset_helix = QPushButton("Reset Helix")
        self.btn_reset_helix.clicked.connect(self.reset_helix)
        helix_layout.addWidget(self.btn_reset_helix)

        # New button for clearing helix lines
        self.btn_clear_helix = QPushButton("Clear Helix Lines")
        self.btn_clear_helix.clicked.connect(self.clear_helix_lines)
        helix_layout.addWidget(self.btn_clear_helix)

        # Toggle for Module-Based Fitting
        self.module_based_fitting_checkbox = QCheckBox("Module-Based Fitting")
        self.module_based_fitting_checkbox.setChecked(False)  # Default: disabled
        helix_layout.addWidget(self.module_based_fitting_checkbox)
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

        # A frame line at the bottom for neatness (optional)
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        left_layout.addWidget(line)

        # Sidebar for info
        sidebar = QWidget()
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(5, 5, 5, 5)

        # --- Cluster/Hit Information ---
        info_label = QLabel("Cluster/Hit Information")
        sidebar_layout.addWidget(info_label)

        self.info_panel = QTextEdit()
        self.info_panel.setReadOnly(True)
        sidebar_layout.addWidget(self.info_panel)
        """
        # --- Event Information ---
        event_info_label = QLabel("Event Information")
        sidebar_layout.addWidget(event_info_label)

        
        # Replace QTextEdit with QLabel for fixed event information
        self.event_info_panel = QLabel(
            "<b>ZDC coincidence:</b> Raw: 1,869,750 Live: 1,868,219"
        )
        self.event_info_panel.setWordWrap(True)
        sidebar_layout.addWidget(self.event_info_panel)
"""

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
        self.selected_points_second = []
        self.helix_line = None
        self.helix_points = None
        self.helix_params_initial = None
        self.helix_params_refined = None
        self.line_params_initial = None
        self.line_params_refined = None

        self.filtered_clusters = None

        self.plotter_widget.show_axes()

        # Enable initial point picking for info mode
        self.update_point_picking()

    def _filter_data(
        self, data: pv.PolyData, for_fitting: bool = False
    ) -> Optional[pv.PolyData]:
        """
        General-purpose filtering for clusters and hits based on spatial ranges,
        side selection, and optional helix proximity.

        Parameters:
        -----------
        data : pv.PolyData
            The dataset to filter (clusters or hits).

        Returns:
        --------
        Optional[pv.PolyData]
            The filtered data or None if no points pass.
        """
        if data is None:
            return None

        # Retrieve range values from sliders
        x_min, x_max = self.x_range_slider.value()
        y_min, y_max = self.y_range_slider.value()
        z_min, z_max = self.z_range_slider.value()

        # Retrieve selected sides
        selected_sides = []
        if self.side0_checkbox.isChecked():
            selected_sides.append(0)
        if self.side1_checkbox.isChecked():
            selected_sides.append(1)

        if not selected_sides:
            return None  # No sides selected, no points pass

        points = data.points
        side = data.point_data.get("side", np.zeros(data.n_points))

        # Apply spatial and side-based filtering
        side_mask = np.isin(side, selected_sides)
        spatial_mask = (
            (points[:, 0] >= x_min)
            & (points[:, 0] <= x_max)
            & (points[:, 1] >= y_min)
            & (points[:, 1] <= y_max)
            & (points[:, 2] >= z_min)
            & (points[:, 2] <= z_max)
        )
        combined_mask = side_mask & spatial_mask

        # --- ADC Filtering ---
        adc_values = data.point_data.get("adc", None)
        if adc_values is not None:
            # Determine threshold based on data type (cluster or hit)
            if data is self.cluster_data:
                adc_threshold = self.cluster_adc_spinbox.value()
            elif data is self.hit_data:
                adc_threshold = self.hit_adc_spinbox.value()
            else:
                adc_threshold = 0

            adc_mask = adc_values > adc_threshold
            combined_mask = combined_mask & adc_mask

        used_in_seed = data.point_data.get("used_in_seed", None)
        if used_in_seed is not None and self.seed_checkbox.isChecked():
            # Filter: only clusters with used_in_seed == 1
            combined_mask = combined_mask & (used_in_seed == 1)

        used_in_track = data.point_data.get("used_in_track", None)
        if used_in_track is not None and self.track_checkbox.isChecked():
            # Filter: only clusters with used_in_track == 1
            combined_mask = combined_mask & (used_in_track == 1)

        filtered_indices = np.where(combined_mask)[0]
        filtered_points = data.extract_points(filtered_indices)
        if for_fitting:
            if self.helix_params_initial:
                # Apply helix-based filtering using initial helix parameters
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
                helix_filtered = filtered_points.extract_points(distance_mask)
                return helix_filtered if helix_filtered.n_points > 0 else None

        else:
            # If the filter checkbox is checked and initial helix exists, apply helix proximity filtering
            if (
                not self.toggle_filter_checkbox.isChecked()
                and self.helix_params_initial
            ):
                # Apply helix-based filtering using initial helix parameters
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
                helix_filtered = filtered_points.extract_points(distance_mask)
                return helix_filtered if helix_filtered.n_points > 0 else None
            else:
                return filtered_points if filtered_points.n_points > 0 else None

    def _filter_data_line(
        self, data: pv.PolyData, for_fitting: bool = False
    ) -> Optional[pv.PolyData]:
        """
        Example: Filter clusters by distance from the initial line in XY-plane (rphi_window)
                 and by Z-window. This is analogous to 'apply_helix_filter' but for lines.
        """
        if data is None:
            return None

        # 1) Spatial ranges from sliders
        x_min, x_max = self.x_range_slider.value()
        y_min, y_max = self.y_range_slider.value()
        z_min, z_max = self.z_range_slider.value()

        # 2) Side selection
        selected_sides = []
        if self.side0_checkbox.isChecked():
            selected_sides.append(0)
        if self.side1_checkbox.isChecked():
            selected_sides.append(1)
        if not selected_sides:
            return None  # No sides selected => no points

        # 3) Access point arrays
        points = data.points
        side = data.point_data.get("side", np.zeros(data.n_points))

        # 4) Spatial + side masks
        side_mask = np.isin(side, selected_sides)
        spatial_mask = (
            (points[:, 0] >= x_min)
            & (points[:, 0] <= x_max)
            & (points[:, 1] >= y_min)
            & (points[:, 1] <= y_max)
            & (points[:, 2] >= z_min)
            & (points[:, 2] <= z_max)
        )
        combined_mask = side_mask & spatial_mask

        # 5) ADC filtering
        adc_values = data.point_data.get("adc", None)
        if adc_values is not None:
            # Decide threshold based on data type
            if data is self.cluster_data:
                adc_threshold = self.cluster_adc_spinbox.value()
            elif data is self.hit_data:
                adc_threshold = self.hit_adc_spinbox.value()
            else:
                adc_threshold = 0
            adc_mask = adc_values > adc_threshold
            combined_mask &= adc_mask

        # 6) used_in_seed / used_in_track filters
        used_in_seed = data.point_data.get("used_in_seed", None)
        if used_in_seed is not None and self.seed_checkbox.isChecked():
            combined_mask &= used_in_seed == 1

        used_in_track = data.point_data.get("used_in_track", None)
        if used_in_track is not None and self.track_checkbox.isChecked():
            combined_mask &= used_in_track == 1

        # 7) Extract points that pass the above filters
        filtered_indices = np.where(combined_mask)[0]
        filtered_points = data.extract_points(filtered_indices)
        if filtered_points is None or filtered_points.n_points == 0:
            return None

        # 8) Optionally apply line-based filter for fitting or if user toggles
        if for_fitting:
            # We only apply line filtering if we have an initial line fit
            if self.line_params_initial:
                distance_mask = np.array(
                    [
                        apply_line_filter(
                            pt,
                            self.line_params_initial,
                            self.rphi_window,
                            self.z_window,
                        )
                        for pt in filtered_points.points
                    ]
                )
                line_filtered = filtered_points.extract_points(distance_mask)
                return line_filtered if line_filtered.n_points > 0 else None
            else:
                # If we don't have line parameters yet, just return the spatially filtered data
                return filtered_points
        else:
            
            if (
                not self.toggle_filter_checkbox.isChecked()  # analogous to helix usage
                and self.line_params_initial is not None
            ):
                distance_mask = np.array(
                    [
                        apply_line_filter(
                            pt,
                            self.line_params_initial,
                            self.rphi_window,
                            self.z_window,
                        )
                        for pt in filtered_points.points
                    ]
                )
                line_filtered = filtered_points.extract_points(distance_mask)
                return line_filtered if line_filtered.n_points > 0 else None
            else:
                # Return just the spatially filtered data
                return filtered_points if filtered_points.n_points > 0 else None

    def add_file_checkbox(self, filename, cluster_data, hit_data):
        checkbox = QCheckBox(filename)
        checkbox.setChecked(True)  # Display file data by default
        checkbox.stateChanged.connect(self.update_display)
        self.file_toggle_layout.addWidget(checkbox)

        self.loaded_files[filename] = {
            "cluster": cluster_data,
            "hit": hit_data,
            "checkbox": checkbox,
        }

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
        if point_id < 0 or mesh is None:
            return  # No valid point was picked

        # Ensure 'data_type' exists in the mesh's point data
        if "data_type" not in mesh.point_data:
            return  # Ignore picking on meshes without 'data_type'

        # Retrieve the data_type for the picked point
        data_type = mesh.point_data["data_type"][point_id]
        label = "Cluster" if data_type == 0 else "Hit"

        if data_type != 0:
            QMessageBox.warning(
                self,
                "Helix Fit",
                "Please select a Cluster point for helix fitting.",
            )
            return  # Only allow picking clusters for fitting
        # Display information about the picked point
        info_text = f"{label} (Point ID: {point_id})\n"
        for attr in mesh.point_data.keys():
            value = mesh.point_data[attr][point_id]
            info_text += f"{attr}: {value}\n"

        self.info_panel.setText(info_text)

        # Store the picked point for helix fitting
        # mesh.points is a numpy array containing the coordinates
        picked_coordinates = mesh.points[point_id]

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
        """
        # Add a marker for the selected point
        marker = pv.Sphere(radius=2, center=picked_coordinates)
        self.plotter_widget.add_mesh(
            marker,
            color="yellow",
            pickable=False,
        )
        """

    def load_data(self):
        """
        Load data from a ROOT file, process it, and display it in the 3D view.
        """

        filenames, _ = QFileDialog.getOpenFileNames(
            self, "Open ROOT File", "", "ROOT files (*.root)"
        )
        if not filenames:
            return  # User canceled the file dialog

        for filename in filenames:

            try:
                # Use the data_loader module to load data
                cluster_polydata, hit_polydata = load_data_from_root(filename)

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

            self.add_file_to_combo(filename, cluster_polydata, hit_polydata)
        self.update_display()

    def add_file_to_combo(self, filename, cluster_data, hit_data):
        # Create a checkable item for the file
        item = QStandardItem(filename)
        item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
        item.setData(Qt.Checked, Qt.CheckStateRole)  # Default to checked

        # Add item to the combo box model
        self.file_combo.model().appendRow(item)
        assigned_color = self.file_colors[self.color_index % len(self.file_colors)]
        self.color_index += 1

        # Store the file data and its associated item
        self.loaded_files[filename] = {
            "cluster": cluster_data,
            "hit": hit_data,
            "item": item,
            "color": assigned_color,
        }

    def update_display(self):
        """
        Update the 3D visualization based on loaded data and current settings.
        Includes optional filtering based on helix tube.
        """
        # Save current camera position
        camera_position = self.plotter_widget.camera_position
        # Clear the current plotter
        self.plotter_widget.clear()

        # Determine what to show
        show_clusters = self.show_clusters.isChecked()
        show_hits = self.show_hits.isChecked()
        for filename, file_info in self.loaded_files.items():
            item = file_info["item"]
            if item.checkState() != Qt.Checked:
                continue

            # --- Clusters ---
            if (
                show_clusters
                and file_info["cluster"] is not None
                and file_info["cluster"].n_points > 0
            ):
                if not self.using_line_fitting:
                    
                    filtered_clusters_display = self._filter_data(file_info["cluster"])
                else:
                    filtered_clusters_display = self._filter_data_line(file_info["cluster"])
                
                    
                if filtered_clusters_display and filtered_clusters_display.n_points > 0:
                    cluster_color = file_info.get("color", "red")
                    self.plotter_widget.add_mesh(
                        filtered_clusters_display,
                        style="points",
                        point_size=5,
                        # color="red",
                        color=cluster_color,
                    )

            # --- Hits ---
            if (
                show_hits
                and file_info["hit"] is not None
                and file_info["hit"].n_points > 0
            ):
                filtered_hits_display = self._filter_data(file_info["hit"])

                if filtered_hits_display and filtered_hits_display.n_points > 0:
                    self.plotter_widget.add_mesh(
                        filtered_hits_display,
                        style="points",
                        point_size=5,
                        color="blue",
                    )
        """
        for mesh in self.helix_lines:
            self.plotter_widget.add_mesh(
                mesh,  # mesh is a pyvista mesh
                color=self.helix_line_color,
                line_width=3,
                style="wireframe",
                pickable=False,
            )
            
        """

        for actor in self.track_lines:
            # Re-add the stored actor to the renderer
            self.plotter_widget.renderer.add_actor(actor)

        # Disable and re-enable picking to ensure a fresh start
        self.plotter_widget.disable_picking()
        # Enable point picking based on current mode
        self.update_point_picking()
        # Restore the previously saved camera position
        self.plotter_widget.camera_position = camera_position

        self.create_tpc_cylinders()
        # self.plotter_widget.set_scale(xscale=1, yscale=1, zscale=1)

        # Update the camera clipping range to include all data
        # self.plotter_widget.reset_camera_clipping_range()

        self.plotter_widget.show_axes()

    def create_tpc_cylinders(self):
        """
        Create concentric cylinders to represent the TPC detector.
        Inner radius: 21.6 cm
        Outer radius: 76.4 cm
        Half-length: 105.5 cm (total length 211 cm)
        """
        # Create inner cylinder
        inner_cylinder = pv.Cylinder(
            center=(0, 0, 0),
            direction=(0, 0, 1),
            radius=21.6,
            height=211,
        )

        # Create outer cylinder
        outer_cylinder = pv.Cylinder(
            center=(0, 0, 0),
            direction=(0, 0, 1),
            radius=76.4,
            height=211,
        )

        beam_line = pv.Line((0, 0, -200), (0, 0, 200))
        # Add cylinders to the plotter with transparency
        self.plotter_widget.add_mesh(
            inner_cylinder,
            color="gray",
            opacity=0.2,
            pickable=False,
            line_width=1,
            style="wireframe",
        )

        self.plotter_widget.add_mesh(
            outer_cylinder,
            color="gray",
            opacity=0.2,
            pickable=False,
            line_width=1,
            style="wireframe",
        )
        # Add beam axis line
        self.plotter_widget.add_mesh(
            beam_line,
            color="grey",
            line_width=3,
            pickable=False,
        )

    def initiate_fit(self):
        """
        Triggered when the user clicks 'Fit Track'.
        If self.using_line_fitting=True, do line approach.
        Otherwise, do helix approach.
        """
        if self.pick_mode != "helix":
            QMessageBox.warning(
                self,
                "Track Fit",
                "Please switch to 'Pick for Helix/Line Fitting' mode to fit a track.",
            )
            return

       

        # Combine all selected cluster data
        selected_clusters = []
        for filename, file_info in self.loaded_files.items():
            item = file_info["item"]
            if item.checkState() == Qt.Checked and file_info["cluster"] is not None:
                selected_clusters.append(file_info["cluster"])
        if not selected_clusters:
            QMessageBox.warning(self, "Fit Error", "No cluster data loaded.")
            return
        combined_cluster_data = pv.merge(selected_clusters)
        self.cluster_data = combined_cluster_data

        # --------------- If Using Line Fitting ---------------
        if self.using_line_fitting:
            self.do_line_fitting_flow()
        else:
            self.do_helix_fitting_flow()

        self.selected_points_first.clear()
        self.update_display()

    def do_line_fitting_flow(self):
        # Require exactly three points for the initial fit
        if len(self.selected_points_first) != 2:
            QMessageBox.warning(
                self,
                "Track Fit",
                "Please select exactly two cluster points before fitting.",
            )
            return
        """
        Perform the 2-step line fitting approach (initial + direct).
        """
        from line_fitting import (
            fit_line_initial,
            fit_line_direct,
            generate_line_points,
            generate_line_polydata,
        )

        # 1) Initial line fit from 2 picks
        line_params_init = fit_line_initial(self.selected_points_first)
        if not line_params_init:
            QMessageBox.warning(self, "Line Fit", "Initial line fitting failed.")
            return
        
        centroid_initial = compute_centroid(self.selected_points_first)

        self.line_params_initial = line_params_init
        # Visualize the initial line
        line_points_init = generate_line_points(line_params_init, self.inner_cut, self.outer_cut, centroid_initial)
        line_poly_init = generate_line_polydata(line_points_init)
        if line_poly_init is not None:
            actor = self.plotter_widget.add_mesh(
                line_poly_init,
                color="grey",
                line_width=3,
                style="wireframe",
                pickable=False,
            )
            self.track_lines.append(actor)

        # 2) Filter clusters around this line
        filtered_clusters_for_fitting = self._filter_data_line(
            self.cluster_data, for_fitting=True
        )
        if (
            filtered_clusters_for_fitting is None
            or filtered_clusters_for_fitting.n_points == 0
        ):
            QMessageBox.warning(
                self, "Line Fit", "No clusters left after line-based filter."
            )
            return

        # 3) Direct line fit with all filtered points
        line_params_refined = fit_line_direct(
            filtered_clusters_for_fitting.points, line_params_init
        )
        if not line_params_refined:
            QMessageBox.warning(self, "Line Fit", "Refined line fitting failed.")
            return
        self.line_params_refined = line_params_refined
        
        centroid_refined = compute_centroid(filtered_clusters_for_fitting.points)
        

        # Visualize refined line
        line_points_refined = generate_line_points(
            line_params_refined, self.inner_cut, self.outer_cut, centroid_refined
        )
        line_poly_refined = generate_line_polydata(line_points_refined)
        if line_poly_refined is not None:
            actor = self.plotter_widget.add_mesh(
                line_poly_refined,
                color="magenta",
                line_width=3,
                style="wireframe",
                pickable=False,
            )
            self.track_lines.append(actor)
            QMessageBox.information(
                self,
                "Line Fit",
                "Refined line fitted and visualized.",
            )

        # After successful refined line fit and visualization
        if filtered_clusters_for_fitting is not None and filtered_clusters_for_fitting.n_points > 0:
            # Calculate line residuals
            delta_rphi, delta_z = calculate_deltas_line(self.line_params_refined, filtered_clusters_for_fitting)
            
            # Increment track counter and update histograms
            self.track_counter += 1
            if not self.histogram_window.isVisible():
                self.histogram_window.show()
            self.histogram_window.add_histograms(delta_rphi, delta_z, self.track_counter)
        else:
            QMessageBox.warning(
                self,
                "Line Fit",
                "No clusters passed the filter. Cannot calculate deltas.",
            )
            return
        
        
        # Optionally compute residuals, do histograms, etc.
        # self.track_counter += 1
        # ...
        
        
    def do_helix_fitting_flow(self):
        
        
         # Require exactly three points for the initial fit
        if len(self.selected_points_first) != 3:
            QMessageBox.warning(
                self,
                "Track Fit",
                "Please select exactly three cluster points before fitting.",
            )
            return
        try:
            helix_params_initial = fit_helix_initial(self.selected_points_first)
        except ValueError as ve:
            QMessageBox.warning(self, "Helix Fit", str(ve))
            return

        # Store initial helix parameters
        self.helix_params_initial = helix_params_initial
        # self.helix_line_color = "green"
        self.helix_line_color = "grey"
        # Generate helix points and create a tube
        helix_points_initial = generate_helix_points_initial(
            helix_params_initial, self.inner_cut, self.outer_cut
        )
        self.helix_points = (
            helix_points_initial  # Store helix points for distance calculations
        )
        self.helix_line = generate_helix_line(helix_points_initial)
        if self.helix_line is not None:
            self.plotter_widget.add_mesh(
                self.helix_line,
                color=self.helix_line_color,
                line_width=3,  # Make line thicker
                style="wireframe",  # Make line dashed
                pickable=False,
            )
        else:
            QMessageBox.warning(
                self, "Helix Fit", "Failed to create initial helix visualization."
            )

        filtered_clusters_for_fitting = self._filter_data(
            self.cluster_data, for_fitting=True
        )
        if (
            filtered_clusters_for_fitting is None
            or filtered_clusters_for_fitting.n_points == 0
        ):
            QMessageBox.warning(
                self,
                "Helix Fit",
                "No clusters passed the filter. Cannot perform refined fitting.",
            )
            return

        if self.module_based_fitting_checkbox.isChecked():
            # Split clusters by layer into modules
            module_clusters = self.split_clusters_by_layer(
                filtered_clusters_for_fitting
            )
            module_order = {"mod1": 0, "mod2": 1, "mod3": 2}
            module_hist_data = {}

            # For each module, perform a second fitting and generate histograms
            for module_name, module_data in module_clusters.items():
                if module_data is None or module_data.n_points == 0:
                    continue  # Skip empty modules

                try:
                    helix_params_module = fit_helix_direct(
                        module_data.points, self.helix_params_initial
                    )
                except ValueError as ve:
                    QMessageBox.warning(
                        self, "Helix Fit", f"Module {module_name} fit error: {ve}"
                    )
                    continue

                if helix_params_module is None:
                    QMessageBox.warning(
                        self,
                        "Helix Fit",
                        f"Refined helix fitting failed for module {module_name}.",
                    )
                    continue

                helix_points_module = generate_helix_points_refined(
                    helix_params_module, self.inner_cut, self.outer_cut
                )
                helix_line_module = generate_helix_line(helix_points_module)
                if helix_line_module is not None:
                    actor = self.plotter_widget.add_mesh(
                        helix_line_module,
                        color="grey",
                        line_width=3,
                        style="wireframe",
                        pickable=False,
                    )
                    self.track_lines.append(actor)

                delta_rphi, delta_z = calculate_deltas(helix_params_module, module_data)
                self.track_counter += 1

                module_idx = module_order.get(module_name, None)
                if module_idx is not None:
                    module_hist_data[module_idx] = (
                        delta_rphi,
                        delta_z,
                        f"{self.track_counter}_{module_name}",
                    )

                sigma_rphi = np.std(delta_rphi)
                sigma_z = np.std(delta_z)
                print(
                    f"Module {module_name}: sigma_rphi={sigma_rphi}, sigma_z={sigma_z}"
                )
            if module_hist_data:
                if not self.module_hist_window.isVisible():
                    self.module_hist_window.show()
                self.module_hist_window.update_module_histograms(module_hist_data)
        else:
            try:
                # Perform refined helix fitting
                helix_params_refined = fit_helix_direct(
                    filtered_clusters_for_fitting.points, self.helix_params_initial
                )
            except ValueError as ve:
                QMessageBox.warning(self, "Helix Fit", str(ve))
                return

            if helix_params_refined is None:
                QMessageBox.warning(self, "Helix Fit", "Refined helix fitting failed.")
                return

            # Store refined helix parameters
            self.helix_params_refined = helix_params_refined

            # self.helix_line_color = "blue"
            self.helix_line_color = "grey"

            # Generate refined helix points and create a tube
            helix_points_refined = generate_helix_points_refined(
                helix_params_refined, self.inner_cut, self.outer_cut
            )
            self.helix_points = (
                helix_points_refined  # Update helix points for distance calculations
            )

            helix_line_refined = generate_helix_line(helix_points_refined)
            if helix_line_refined is not None:
                actor = self.plotter_widget.add_mesh(
                    helix_line_refined,
                    color=self.helix_line_color,
                    line_width=3,
                    style="wireframe",
                    pickable=False,
                )
                self.track_lines.append(actor)
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
                filtered_clusters_for_fitting is not None
                and filtered_clusters_for_fitting.n_points > 0
            ):
                delta_rphi, delta_z = calculate_deltas(
                    self.helix_params_refined, filtered_clusters_for_fitting
                )
                self.track_counter += 1
                if not self.histogram_window.isVisible():
                    self.histogram_window.show()
                self.histogram_window.add_histograms(
                    delta_rphi, delta_z, self.track_counter
                )

                # **Calculate sigma (standard deviation)**
                sigma_rphi = np.std(delta_rphi)
                sigma_z = np.std(delta_z)

                track_id = self.track_counter
                """
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
                """
                # **Update info panel with sigma values**
                info_text = (
                    f"Track ID: {track_id}\n"
                    f"Sigma Delta rphi: {sigma_rphi:.4f} cm\n"
                    f"Sigma Delta z: {sigma_z:.4f} cm\n"
                )
                self.info_panel.setText(info_text)

                # **Plot histograms within the GUI**
                # self.plot_histograms(delta_rphi, delta_z, track_id)

            else:
                QMessageBox.warning(
                    self,
                    "Helix Fit",
                    "No clusters passed the filter. Cannot calculate deltas.",
                )
                return

        # **Reset for next fitting**
        self.selected_points_first.clear()
        self.instruction_label.setText(
            "Instruction: Select 3 points for initial helix fitting."
        )
        self.update_display()

    '''def plot_histograms(self, delta_rphi, delta_z, track_id):
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
        self.hist_canvas.draw()'''

    def reset_helix(self):
        """Clears the fitted helix and resets the selection."""
        self.helix_line = None
        self.helix_points = None  # Clear helix points
        self.helix_params_initial = None
        self.helix_params_refined = None
        self.line_params_initial = None
        self.line_params_refined = None
        
        self.selected_points_first.clear()
        self.instruction_label.setText(
            "Instruction: Select 3 points for initial helix fitting."
        )
        self.update_display()
        QMessageBox.information(
            self,
            "Helix Reset",
            "Helix has been reset. You can select new points and fit again.",
        )

    def split_clusters_by_layer(self, filtered_clusters):
        """
        Splits clusters into three modules based on layer numbers:
        - Module 1: layers 7 to 22
        - Module 2: layers 23 to 38
        - Module 3: layers 39 to 54
        """
        if filtered_clusters is None or filtered_clusters.n_points == 0:
            return {}

        # Retrieve layer values from cluster point_data
        layers = filtered_clusters.point_data.get("layer", None)
        if layers is None:
            return {}

        # Determine indices for each module based on layer ranges
        mod1_indices = np.where((layers >= 7) & (layers <= 22))[0]
        mod2_indices = np.where((layers >= 23) & (layers <= 38))[0]
        mod3_indices = np.where((layers >= 39) & (layers <= 54))[0]

        mod1 = (
            filtered_clusters.extract_points(mod1_indices)
            if mod1_indices.size > 0
            else None
        )
        mod2 = (
            filtered_clusters.extract_points(mod2_indices)
            if mod2_indices.size > 0
            else None
        )
        mod3 = (
            filtered_clusters.extract_points(mod3_indices)
            if mod3_indices.size > 0
            else None
        )

        return {"mod1": mod1, "mod2": mod2, "mod3": mod3}

    def clear_helix_lines(self):
        """Clears all displayed helix lines from the plot."""
        # Clear the list of helix lines
        self.track_lines = []
        # Optionally update the display to refresh the view
        self.update_display()
