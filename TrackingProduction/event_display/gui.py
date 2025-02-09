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
    QSpinBox,
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
    find_vertex_z,
    find_dca_and_closest_point,
)


from histogram_window import HistogramWindow
from module_histogram_window import ModuleHistogramWindow
from projection_window import GeometricProjectionWindow
from analysis import (
    calculate_deltas,
    apply_helix_filter,
    apply_line_filter,
    calculate_deltas_line,
    compute_centroid,
    calculate_dca_to_beam,
)
from data_loader import load_data_from_root
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtGui import QStandardItem, QStandardItemModel


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Display")
        self.resize(1600, 900)
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
        self.using_line_fitting = False

        # Initialize pick mode
        self.pick_mode = "info"  # Default mode

        self.inner_cut = 21.6
        self.outer_cut = 76.4

        self.helix_line = None
        self.track_lines = []  # List to store all track lines

        self.display_filtered_clusters_info = []
        self.display_filtered_hits_info = []

        self.helix_line_color = "grey"  # Default color for initial fitting

        self.track_counter = 0

        self.rphi_window = 0.5  # Default window size for rphi
        self.z_window = 1  # Default window size for z

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

        # Left side panel
        """left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(5, 5, 5, 5)"""
        left_splitter = QSplitter(Qt.Vertical)
        left_splitter.setContentsMargins(5, 5, 5, 5)

        # === Control Area ===
        control_area = QWidget()
        control_layout = QVBoxLayout(control_area)
        control_layout.setContentsMargins(0, 0, 0, 0)

        # --- Data Loading Button ---
        btn_load = QPushButton("Load ROOT File")
        btn_load.clicked.connect(self.load_data)
        control_layout.addWidget(btn_load)

        top_filter_layout = QHBoxLayout()

        # --------------------------------------------------------------------
        # 2) Silicon group

        silicon_group = QGroupBox("Silicon")
        silicon_layout = QVBoxLayout()
        silicon_group.setLayout(silicon_layout)

        # Create a horizontal layout for the checkboxes
        checkbox_layout = QHBoxLayout()

        # Left column (Seed and Track)
        left_column = QVBoxLayout()
        self.seed_checkbox_silicon = QCheckBox("Seed")
        self.seed_checkbox_silicon.setChecked(False)
        self.seed_checkbox_silicon.stateChanged.connect(self.update_display)
        left_column.addWidget(self.seed_checkbox_silicon)

        self.track_checkbox_silicon = QCheckBox("Track")
        self.track_checkbox_silicon.setChecked(False)
        self.track_checkbox_silicon.stateChanged.connect(self.update_display)
        left_column.addWidget(self.track_checkbox_silicon)

        # Right column (Clusters and Hits)
        right_column = QVBoxLayout()
        self.clusters_checkbox_silicon = QCheckBox("Clusters")
        self.clusters_checkbox_silicon.setChecked(True)
        self.clusters_checkbox_silicon.stateChanged.connect(self.update_display)
        right_column.addWidget(self.clusters_checkbox_silicon)

        self.hits_checkbox_silicon = QCheckBox("Hits")
        self.hits_checkbox_silicon.setChecked(True)
        self.hits_checkbox_silicon.stateChanged.connect(self.update_display)
        right_column.addWidget(self.hits_checkbox_silicon)

        # Add columns to checkbox layout
        checkbox_layout.addLayout(left_column)
        checkbox_layout.addLayout(right_column)
        silicon_layout.addLayout(checkbox_layout)

        # crossing
        crossing_label_silicon = QLabel("Crossing:")
        self.crossing_spin_silicon = QSpinBox()
        self.crossing_spin_silicon.setRange(-2000, 9999)
        self.crossing_spin_silicon.setValue(-2000)
        self.crossing_spin_silicon.valueChanged.connect(self.update_display)
        silicon_layout.addWidget(crossing_label_silicon)
        silicon_layout.addWidget(self.crossing_spin_silicon)

        # T-crossing
        t_crossing_label_silicon = QLabel("T-Crossing:")
        self.t_crossing_spin_silicon = QSpinBox()
        self.t_crossing_spin_silicon.setRange(-2000, 99999)
        self.t_crossing_spin_silicon.setValue(-2000)
        self.t_crossing_spin_silicon.valueChanged.connect(self.update_display)
        silicon_layout.addWidget(t_crossing_label_silicon)
        silicon_layout.addWidget(self.t_crossing_spin_silicon)

        # track ID
        track_id_label_silicon = QLabel("Track ID:")
        self.track_id_spin_silicon = QSpinBox()
        self.track_id_spin_silicon.setRange(-1, 999999)
        self.track_id_spin_silicon.setValue(-1)
        self.track_id_spin_silicon.valueChanged.connect(self.update_display)
        silicon_layout.addWidget(track_id_label_silicon)
        silicon_layout.addWidget(self.track_id_spin_silicon)

        top_filter_layout.addWidget(silicon_group)

        # --------------------------------------------------------------------
        # 3) TPC group

        tpc_group = QGroupBox("TPC")
        tpc_layout = QVBoxLayout()
        tpc_group.setLayout(tpc_layout)

        checkbox_layout = QHBoxLayout()

        left_column = QVBoxLayout()
        self.seed_checkbox_tpc = QCheckBox("Seed")
        self.seed_checkbox_tpc.setChecked(False)
        self.seed_checkbox_tpc.stateChanged.connect(self.update_display)
        left_column.addWidget(self.seed_checkbox_tpc)

        self.track_checkbox_tpc = QCheckBox("Track")
        self.track_checkbox_tpc.setChecked(False)
        self.track_checkbox_tpc.stateChanged.connect(self.update_display)
        left_column.addWidget(self.track_checkbox_tpc)

        right_column = QVBoxLayout()

        self.clusters_checkbox_tpc = QCheckBox("Clusters")
        self.clusters_checkbox_tpc.setChecked(True)
        self.clusters_checkbox_tpc.stateChanged.connect(self.update_display)
        right_column.addWidget(self.clusters_checkbox_tpc)

        self.hits_checkbox_tpc = QCheckBox("Hits")
        self.hits_checkbox_tpc.setChecked(True)
        self.hits_checkbox_tpc.stateChanged.connect(self.update_display)
        right_column.addWidget(self.hits_checkbox_tpc)

        checkbox_layout.addLayout(left_column)
        checkbox_layout.addLayout(right_column)
        tpc_layout.addLayout(checkbox_layout)

        t_crossing_label_tpc = QLabel("T-Crossing:")
        self.t_crossing_spin_tpc = QSpinBox()
        self.t_crossing_spin_tpc.setRange(-2000, 99999)
        self.t_crossing_spin_tpc.setValue(-2000)
        self.t_crossing_spin_tpc.valueChanged.connect(self.update_display)
        tpc_layout.addWidget(t_crossing_label_tpc)
        tpc_layout.addWidget(self.t_crossing_spin_tpc)

        track_id_label_tpc = QLabel("Track ID:")
        self.track_id_spin_tpc = QSpinBox()
        self.track_id_spin_tpc.setRange(-1, 999999)
        self.track_id_spin_tpc.setValue(-1)
        self.track_id_spin_tpc.valueChanged.connect(self.update_display)
        tpc_layout.addWidget(track_id_label_tpc)
        tpc_layout.addWidget(self.track_id_spin_tpc)

        self.side0_checkbox_tpc = QCheckBox("Side 0")
        self.side0_checkbox_tpc.setChecked(True)
        self.side0_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.side0_checkbox_tpc)

        self.side1_checkbox_tpc = QCheckBox("Side 1")
        self.side1_checkbox_tpc.setChecked(True)
        self.side1_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.side1_checkbox_tpc)

        top_filter_layout.addWidget(tpc_group)

        # --------------------------------------------------------------------

        # --- ADC Threshold Controls ---
        adc_group = QGroupBox("Thresholds")
        adc_layout = QVBoxLayout()
        adc_group.setLayout(adc_layout)

        # Cluster ADC threshold control
        cluster_adc_layout = QHBoxLayout()
        cluster_adc_label = QLabel("Cluster ADC ≥")
        self.cluster_adc_spinbox = QDoubleSpinBox()
        self.cluster_adc_spinbox.setRange(0, 10000)
        self.cluster_adc_spinbox.setValue(0)
        self.cluster_adc_spinbox.setDecimals(0)
        self.cluster_adc_spinbox.valueChanged.connect(self.update_display)
        cluster_adc_layout.addWidget(cluster_adc_label)
        cluster_adc_layout.addWidget(self.cluster_adc_spinbox)
        adc_layout.addLayout(cluster_adc_layout)

        # Hit ADC threshold control
        hit_adc_layout = QHBoxLayout()
        hit_adc_label = QLabel("Hit ADC ≥")
        self.hit_adc_spinbox = QDoubleSpinBox()
        self.hit_adc_spinbox.setRange(0, 10000)
        self.hit_adc_spinbox.setValue(0)
        self.hit_adc_spinbox.setDecimals(0)
        self.hit_adc_spinbox.valueChanged.connect(self.update_display)
        hit_adc_layout.addWidget(hit_adc_label)
        hit_adc_layout.addWidget(self.hit_adc_spinbox)
        adc_layout.addLayout(hit_adc_layout)

        nmaps_layout = QHBoxLayout()
        nmaps_label = QLabel("NMAPS ≥")
        self.nmaps_spinbox = QSpinBox()
        self.nmaps_spinbox.setRange(-10, 100)
        self.nmaps_spinbox.setValue(-1)
        self.nmaps_spinbox.valueChanged.connect(self.update_display)
        nmaps_layout.addWidget(nmaps_label)
        nmaps_layout.addWidget(self.nmaps_spinbox)
        adc_layout.addLayout(nmaps_layout)

        top_filter_layout.addWidget(adc_group)

        # === Pick Mode Selection Area ===
        pick_mode_group = QGroupBox("Pick Mode")
        pick_mode_layout = QVBoxLayout()
        pick_mode_group.setLayout(pick_mode_layout)

        # Radio Button for Viewing Point Info
        self.radio_view_info = QRadioButton("View Point Info")
        self.radio_view_info.setChecked(True)
        self.radio_view_info.toggled.connect(self.on_pick_mode_changed)
        pick_mode_layout.addWidget(self.radio_view_info)

        # Radio Button for Picking Points for Helix Fitting
        self.radio_pick_helix = QRadioButton("Pick for Helix Fitting")
        self.radio_pick_helix.toggled.connect(self.on_pick_mode_changed)
        pick_mode_layout.addWidget(self.radio_pick_helix)

        self.radio_pick_vertex = QRadioButton("Pick for Vertex Finding")
        self.radio_pick_vertex.toggled.connect(self.on_pick_mode_changed)
        pick_mode_layout.addWidget(self.radio_pick_vertex)

        top_filter_layout.addWidget(pick_mode_group)

        # --- File Selection Dropdown ---
        files_group = QGroupBox("Files")
        files_layout = QVBoxLayout()
        files_group.setLayout(files_layout)

        self.file_combo = QComboBox()
        # Set up a model to allow checkable items in the combo box
        model = QStandardItemModel(self.file_combo)
        self.file_combo.setModel(model)

        model.itemChanged.connect(lambda item: self.update_display())

        files_layout.addWidget(self.file_combo)

        self.event_combo = QComboBox()
        event_model = QStandardItemModel(self.event_combo)
        self.event_combo.setModel(event_model)

        for event_num in range(31):  # 0 to 30
            item = QStandardItem(f"Event {event_num}")
            item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
            item.setData(Qt.Unchecked, Qt.CheckStateRole)
            event_model.appendRow(item)

        event_model.itemChanged.connect(lambda i: self.update_display())

        files_layout.addWidget(self.event_combo)

        top_filter_layout.addWidget(files_group)

        projection_buttons_layout = QVBoxLayout()
        projection_buttons_layout.setContentsMargins(0, 0, 0, 0)

        self.btn_show_projection = QPushButton("Show 2D Projection")
        self.btn_show_projection.clicked.connect(self.show_projection_window)
        projection_buttons_layout.addWidget(self.btn_show_projection)

        self.btn_toggle_field = QPushButton("Field On")
        self.btn_toggle_field.setCheckable(True)
        self.btn_toggle_field.clicked.connect(self.toggle_field_mode)
        projection_buttons_layout.addWidget(self.btn_toggle_field)

        top_filter_layout.addLayout(projection_buttons_layout)

        self.projection_window = None

        event_info_label = QLabel(
            "<b>Run:</b> 53217<br>"
            "<b>ZDC coincidence:</b> Raw: 1,869750<br>"
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
            "Live: 1,868219"
        )
        event_info_label.setWordWrap(True)
        top_filter_layout.addWidget(event_info_label)

        control_layout.addLayout(top_filter_layout)

        # === Range Selection Area ===
        range_group = QGroupBox("Axis Range Selection")
        range_layout = QHBoxLayout()
        range_group.setLayout(range_layout)

        slider_width = 200
        label_width = 30
        spinbox_width = 50

        overall_ranges = {"X": (-150, 150), "Y": (-150, 150), "Z": (-400, 400)}

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

            def on_min_spin_change(value):
                current_max = slider.value()[1]
                if value > current_max:
                    min_spin.blockSignals(True)
                    min_spin.setValue(current_max)
                    min_spin.blockSignals(False)
                    value = current_max
                slider.setValue([value, current_max])
                self.update_display()

            def on_max_spin_change(value):
                current_min = slider.value()[0]
                if value < current_min:
                    max_spin.blockSignals(True)
                    max_spin.setValue(current_min)
                    max_spin.blockSignals(False)
                    value = current_min
                slider.setValue([current_min, value])
                self.update_display()

            min_spin.valueChanged.connect(on_min_spin_change)
            max_spin.valueChanged.connect(on_max_spin_change)

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

        # === Track Fitting Controls ===
        track_group = QGroupBox("Track Fitting")
        track_layout = QHBoxLayout()
        track_group.setLayout(track_layout)

        self.btn_fit_track = QPushButton("Fit Track")
        self.btn_fit_track.clicked.connect(self.initiate_fit)
        track_layout.addWidget(self.btn_fit_track)

        self.btn_reset_track = QPushButton("Reset Track")
        self.btn_reset_track.clicked.connect(self.reset_track)
        track_layout.addWidget(self.btn_reset_track)

        self.btn_clear_track = QPushButton("Clear Track Lines")
        self.btn_clear_track.clicked.connect(self.clear_track_lines)
        track_layout.addWidget(self.btn_clear_track)

        self.module_based_fitting_checkbox = QCheckBox("Module Fitting")
        self.module_based_fitting_checkbox.setChecked(False)
        track_layout.addWidget(self.module_based_fitting_checkbox)

        self.track_toggle_checkbox = QCheckBox("Show Lines")
        self.track_toggle_checkbox.setChecked(True)
        self.track_toggle_checkbox.stateChanged.connect(self.on_track_toggle)

        track_layout.addWidget(self.track_toggle_checkbox)

        self.toggle_filter_checkbox = QCheckBox("Show Points Outside Tube")
        self.toggle_filter_checkbox.setChecked(True)
        self.toggle_filter_checkbox.stateChanged.connect(self.update_display)
        track_layout.addWidget(self.toggle_filter_checkbox)

        self.instruction_label = QLabel(
            "Instruction: Select 3 points for initial track fitting."
        )
        track_layout.addWidget(self.instruction_label)

        control_layout.addWidget(track_group)

        # left_layout.addWidget(control_area)

        # 3D Visualization Area
        self.plotter_widget = QtInteractor()
        # left_layout.addWidget(self.plotter_widget)

        left_splitter.addWidget(control_area)
        left_splitter.addWidget(self.plotter_widget)

        # line = QFrame()
        # line.setFrameShape(QFrame.HLine)
        # line.setFrameShadow(QFrame.Sunken)
        # left_layout.addWidget(line)

        # Sidebar for info
        sidebar = QWidget()
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(5, 5, 5, 5)

        info_label = QLabel("Cluster/Hit Information")
        sidebar_layout.addWidget(info_label)

        self.info_panel = QTextEdit()
        self.info_panel.setReadOnly(True)
        sidebar_layout.addWidget(self.info_panel)

        splitter.addWidget(left_splitter)
        splitter.addWidget(sidebar)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        # Data placeholders
        self.cluster_data = None
        self.hit_data = None
        self.cluster_polydata = None
        self.hit_polydata = None

        self.selected_points_first = []  # Used for point selection
        self.helix_points = None
        self.helix_params_initial = None
        self.helix_params_refined = None
        self.line_params_initial = None
        self.line_params_refined = None

        self.filtered_clusters = None

        self.plotter_widget.show_axes()

        self.update_point_picking()

    def _filter_data(
        self,
        data: pv.PolyData,
        *,
        for_fitting: bool = False,
    ) -> Optional[pv.PolyData]:

        if data is None or data.n_points == 0:
            # print(f"DEBUG {'(fitting)' if for_fitting else ''}: No data provided.")
            return None

        points = data.points
        # print(
        #   f"DEBUG {'(fitting)' if for_fitting else ''}: Total points before filtering: {len(points)}"
        # )

        event_numbers = data.point_data.get("event", None)
        if event_numbers is None:
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Missing 'event' attribute in data."
            # )
            return None

        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Unique event numbers in data: {np.unique(event_numbers)}"
        # )
        selected_events = []
        for i in range(self.event_combo.model().rowCount()):
            item = self.event_combo.model().item(i)
            if item.checkState() == Qt.Checked:
                event_label = item.text()
                try:

                    event_id = int(event_label.split()[1])
                    selected_events.append(event_id)
                except (ValueError, IndexError):

                    continue

        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Selected events: {selected_events}"
        # )

        """if not selected_events:
            print("No events selected.")
            return None  # No events selected"""

        if not selected_events:
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: No events selected, returning None"
            )
            return None
        event_mask = np.isin(event_numbers, selected_events)
        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after event filter: {np.sum(event_mask)}"
        # )

        x_min, x_max = self.x_range_slider.value()
        y_min, y_max = self.y_range_slider.value()
        z_min, z_max = self.z_range_slider.value()
        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Spatial ranges: X[{x_min}, {x_max}], Y[{y_min}, {y_max}], Z[{z_min}, {z_max}]"
        # )

        spatial_mask = (
            (points[:, 0] >= x_min)
            & (points[:, 0] <= x_max)
            & (points[:, 1] >= y_min)
            & (points[:, 1] <= y_max)
            & (points[:, 2] >= z_min)
            & (points[:, 2] <= z_max)
        )

        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: X range in data: [{np.min(points[:, 0])}, {np.max(points[:, 0])}]"
        # )
        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Y range in data: [{np.min(points[:, 1])}, {np.max(points[:, 1])}]"
        # )
        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Z range in data: [{np.min(points[:, 2])}, {np.max(points[:, 2])}]"
        # )

        # print(f"Points after spatial filter: {np.sum(spatial_mask)}")

        combined_mask = spatial_mask & event_mask
        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after combining masks: {np.sum(combined_mask)}"
        # )

        nmaps_values = data.point_data.get("nmaps", None)

        if nmaps_values is not None:
            nmaps_threshold = self.nmaps_spinbox.value()
            nmaps_mask = nmaps_values >= nmaps_threshold
            combined_mask &= nmaps_mask
            """print(
                f"Points after NMAPS filter (>= {nmaps_threshold}): {np.sum(nmaps_mask)}"
            )"""
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after combining NMAPS filter: {np.sum(combined_mask)}"
            # )
        else:
            print("No 'nmaps' attribute found; skipping NMAPS filter.")

        # --- ADC Filtering ---
        adc_values = data.point_data.get("adc", None)
        data_type_array = data.point_data.get("data_type", np.zeros(data.n_points))

        if adc_values is not None:
            # Determine threshold based on data type (cluster or hit)
            if data is self.cluster_data:
                adc_threshold = self.cluster_adc_spinbox.value()
            elif data is self.hit_data:
                adc_threshold = self.hit_adc_spinbox.value()
            else:
                adc_threshold = 0

            adc_mask = adc_values >= adc_threshold
            combined_mask &= adc_mask
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after combining ADC filter: {np.sum(combined_mask)}"
            # )

        else:
            print("No 'adc' attribute found; skipping ADC filter.")

        layer_array = data.point_data.get("layer", np.full(data.n_points, -1))

        crossing_data = data.point_data.get("crossing", np.full(data.n_points, -2000))
        t_crossing_data = data.point_data.get(
            "t_crossing", np.full(data.n_points, -2000)
        )
        track_id_data = data.point_data.get("trackid", np.full(data.n_points, -1))
        used_in_seed = data.point_data.get("used_in_seed", np.zeros(data.n_points))
        used_in_track = data.point_data.get("used_in_track", np.zeros(data.n_points))

        is_silicon = layer_array < 7
        is_cluster = data_type_array == 0
        is_hit = data_type_array == 1

        # For silicon points only
        silicon_mask = is_silicon
        show_silicon_clusters = self.clusters_checkbox_silicon.isChecked()
        show_silicon_hits = self.hits_checkbox_silicon.isChecked()

        if show_silicon_clusters and show_silicon_hits:
            pass  # Show both
        elif show_silicon_clusters:
            silicon_mask &= is_cluster
        elif show_silicon_hits:
            silicon_mask &= is_hit
        else:
            silicon_mask &= False  # Show neither

        crossing_val_sil = self.crossing_spin_silicon.value()
        t_cross_val_sil = self.t_crossing_spin_silicon.value()
        track_id_val_sil = self.track_id_spin_silicon.value()

        if crossing_val_sil != -2000:
            silicon_mask &= crossing_data == crossing_val_sil
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after crossing_val_sil filter ({crossing_val_sil}): {np.sum(crossing_data == crossing_val_sil)}"
            # )
        if t_cross_val_sil != -2000:
            silicon_mask &= t_crossing_data == t_cross_val_sil
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after t_cross_val_sil filter ({t_cross_val_sil}): {np.sum(t_crossing_data == t_cross_val_sil)}"
            # )
        if track_id_val_sil != -1:
            silicon_mask &= track_id_data == track_id_val_sil
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after track_id_val_sil filter ({track_id_val_sil}): {np.sum(track_id_data == track_id_val_sil)}"
            # )
        if self.seed_checkbox_silicon.isChecked():
            silicon_mask &= used_in_seed == 1
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after seed_checkbox_silicon filter: {np.sum(used_in_seed == 1)}"
            # )
        if self.track_checkbox_silicon.isChecked():
            silicon_mask &= used_in_track == 1
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after track_checkbox_silicon filter: {np.sum(used_in_track == 1)}"
            # )

        is_tpc = layer_array >= 7

        tpc_mask = is_tpc
        show_tpc_clusters = self.clusters_checkbox_tpc.isChecked()
        show_tpc_hits = self.hits_checkbox_tpc.isChecked()

        if show_tpc_clusters and show_tpc_hits:
            pass  # Show both
        elif show_tpc_clusters:
            tpc_mask &= is_cluster
        elif show_tpc_hits:
            tpc_mask &= is_hit
        else:
            tpc_mask &= False  # Show neither

        t_cross_val_tpc = self.t_crossing_spin_tpc.value()
        track_id_val_tpc = self.track_id_spin_tpc.value()
        side_data = data.point_data.get("side", None)

        if t_cross_val_tpc != -2000:
            tpc_mask &= t_crossing_data == t_cross_val_tpc
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after t_cross_val_tpc filter ({t_cross_val_tpc}): {np.sum(t_crossing_data == t_cross_val_tpc)}"
            # )
        if track_id_val_tpc != -1:
            tpc_mask &= track_id_data == track_id_val_tpc
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after track_id_val_tpc filter ({track_id_val_tpc}): {np.sum(track_id_data == track_id_val_tpc)}"
            # )

        # Sides
        if side_data is not None:
            sides_chosen = []
            if self.side0_checkbox_tpc.isChecked():
                sides_chosen.append(0)
            if self.side1_checkbox_tpc.isChecked():
                sides_chosen.append(1)
            if sides_chosen:
                tpc_mask &= np.isin(side_data, sides_chosen)
                # print(
                #    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after after sides_chosen filter ({sides_chosen}): {np.sum(np.isin(side_data, sides_chosen))}"
                # )
            else:
                # If no sides are checked, TPC mask => no points
                tpc_mask &= False
                print("No sides selected for TPC; excluding all TPC points.")

        else:
            print("Missing 'side' attribute in data; skipping side filter.")

        if self.seed_checkbox_tpc.isChecked():
            tpc_mask &= used_in_seed == 1
        if self.track_checkbox_tpc.isChecked():
            tpc_mask &= used_in_track == 1

        combined_mask = combined_mask & (silicon_mask | tpc_mask)
        if not np.any(combined_mask):
            return None

        filtered_indices = np.where(combined_mask)[0]
        filtered_points = data.extract_points(filtered_indices)
        if filtered_points is None or filtered_points.n_points == 0:
            # print(
            #    f"DEBUG {'(fitting)' if for_fitting else ''}: No points passed the filters"
            # )
            return None
        # print(
        #    f"DEBUG {'(fitting)' if for_fitting else ''}: Final point count: {filtered_points.n_points}"
        # )

        if for_fitting:

            if self.using_line_fitting and self.line_params_initial:
                distance_mask = [
                    apply_line_filter(
                        pt,
                        self.line_params_initial,
                        self.rphi_window,
                        self.z_window,
                    )
                    for pt in filtered_points.points
                ]

                line_filtered = filtered_points.extract_points(distance_mask)
                return line_filtered if line_filtered.n_points > 0 else None
            elif not self.using_line_fitting and self.helix_params_initial:
                # Apply helix-based filtering using initial helix parameters
                print("Applying helix-based proximity filtering for fitting.")
                print("Helix initial params:", self.helix_params_initial)
                print("rphi_window:", self.rphi_window, "z_window:", self.z_window)
                print(
                    f"Points before helix proximity filter: {filtered_points.n_points}"
                )
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
                print(f"Points after helix proximity filter: {helix_filtered.n_points}")
                return helix_filtered if helix_filtered.n_points > 0 else None
            else:
                print(
                    "Helix parameters not initialized; skipping helix proximity filter."
                )
                return None

        else:
            # If the filter checkbox is checked and initial helix exists, apply helix proximity filtering
            if not self.toggle_filter_checkbox.isChecked():
                if self.using_line_fitting and self.line_params_initial is not None:
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

                elif (
                    not self.using_line_fitting
                    and self.helix_params_initial is not None
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
        if data is None or data.n_points == 0:
            print(f"DEBUG {'(fitting)' if for_fitting else ''}: No data provided.")
            return None

        points = data.points
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Total points before filtering: {len(points)}"
        )

        event_numbers = data.point_data.get("event", None)

        if event_numbers is None:
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Missing 'event' attribute in data."
            )
            return None

        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Unique event numbers in data: {np.unique(event_numbers)}"
        )
        selected_events = []
        for i in range(self.event_combo.model().rowCount()):
            item = self.event_combo.model().item(i)
            if item.checkState() == Qt.Checked:
                event_label = item.text()  # e.g., "Event 3"
                try:
                    # Extract numeric part from the label
                    event_id = int(event_label.split()[1])
                    selected_events.append(event_id)
                except (ValueError, IndexError):
                    # print(f"Invalid event ID format: {event_label}")
                    continue

        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Selected events: {selected_events}"
        )

        """if not selected_events:
            print("No events selected.")
            return None  # No events selected"""

        if not selected_events:
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: No events selected, returning None"
            )
            return None
        event_mask = np.isin(event_numbers, selected_events)
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Points after event filter: {np.sum(event_mask)}"
        )

        # 1) Spatial ranges from sliders
        x_min, x_max = self.x_range_slider.value()
        y_min, y_max = self.y_range_slider.value()
        z_min, z_max = self.z_range_slider.value()
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Spatial ranges: X[{x_min}, {x_max}], Y[{y_min}, {y_max}], Z[{z_min}, {z_max}]"
        )
        spatial_mask = (
            (points[:, 0] >= x_min)
            & (points[:, 0] <= x_max)
            & (points[:, 1] >= y_min)
            & (points[:, 1] <= y_max)
            & (points[:, 2] >= z_min)
            & (points[:, 2] <= z_max)
        )

        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: X range in data: [{np.min(points[:, 0])}, {np.max(points[:, 0])}]"
        )
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Y range in data: [{np.min(points[:, 1])}, {np.max(points[:, 1])}]"
        )
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Z range in data: [{np.min(points[:, 2])}, {np.max(points[:, 2])}]"
        )

        combined_mask = spatial_mask & event_mask
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Points after combining masks: {np.sum(combined_mask)}"
        )

        nmaps_values = data.point_data.get("nmaps", None)

        if nmaps_values is not None:
            nmaps_threshold = self.nmaps_spinbox.value()
            nmaps_mask = nmaps_values >= nmaps_threshold
            combined_mask &= nmaps_mask
            """print(
                f"Points after NMAPS filter (>= {nmaps_threshold}): {np.sum(nmaps_mask)}"
            )"""
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after combining NMAPS filter: {np.sum(combined_mask)}"
            )
        else:
            print("No 'nmaps' attribute found; skipping NMAPS filter.")

        # --- ADC Filtering ---
        adc_values = data.point_data.get("adc", None)
        data_type_array = data.point_data.get("data_type", np.zeros(data.n_points))
        # We can do a quick rule:

        if adc_values is not None:
            # Determine threshold based on data type (cluster or hit)
            if data is self.cluster_data:
                adc_threshold = self.cluster_adc_spinbox.value()
            elif data is self.hit_data:
                adc_threshold = self.hit_adc_spinbox.value()
            else:
                adc_threshold = 0

            adc_mask = adc_values >= adc_threshold
            combined_mask &= adc_mask
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after combining ADC filter: {np.sum(combined_mask)}"
            )

        else:
            print("No 'adc' attribute found; skipping ADC filter.")

        layer_array = data.point_data.get("layer", np.full(data.n_points, -1))
        crossing_data = data.point_data.get("crossing", np.full(data.n_points, -2000))
        t_crossing_data = data.point_data.get(
            "t_crossing", np.full(data.n_points, -2000)
        )
        track_id_data = data.point_data.get("trackid", np.full(data.n_points, -1))
        used_in_seed = data.point_data.get("used_in_seed", np.zeros(data.n_points))
        used_in_track = data.point_data.get("used_in_track", np.zeros(data.n_points))
        is_silicon = layer_array < 7
        is_cluster = data_type_array == 0
        is_hit = data_type_array == 1
        silicon_mask = is_silicon
        show_silicon_clusters = self.clusters_checkbox_silicon.isChecked()
        show_silicon_hits = self.hits_checkbox_silicon.isChecked()

        if show_silicon_clusters and show_silicon_hits:
            pass  # Show both
        elif show_silicon_clusters:
            silicon_mask &= is_cluster
        elif show_silicon_hits:
            silicon_mask &= is_hit
        else:
            silicon_mask &= False  # Show neither
        crossing_val_sil = self.crossing_spin_silicon.value()
        t_cross_val_sil = self.t_crossing_spin_silicon.value()
        track_id_val_sil = self.track_id_spin_silicon.value()

        # For silicon points only

        if crossing_val_sil != -2000:
            silicon_mask &= crossing_data == crossing_val_sil
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after crossing_val_sil filter ({crossing_val_sil}): {np.sum(crossing_data == crossing_val_sil)}"
            )
        if t_cross_val_sil != -2000:
            silicon_mask &= t_crossing_data == t_cross_val_sil
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after t_cross_val_sil filter ({t_cross_val_sil}): {np.sum(t_crossing_data == t_cross_val_sil)}"
            )
        if track_id_val_sil != -1:
            silicon_mask &= track_id_data == track_id_val_sil
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after track_id_val_sil filter ({track_id_val_sil}): {np.sum(track_id_data == track_id_val_sil)}"
            )
        if self.seed_checkbox_silicon.isChecked():
            silicon_mask &= used_in_seed == 1
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after seed_checkbox_silicon filter: {np.sum(used_in_seed == 1)}"
            )
        if self.track_checkbox_silicon.isChecked():
            silicon_mask &= used_in_track == 1
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after track_checkbox_silicon filter: {np.sum(used_in_track == 1)}"
            )
        is_tpc = layer_array >= 7
        tpc_mask = is_tpc
        show_tpc_clusters = self.clusters_checkbox_tpc.isChecked()
        show_tpc_hits = self.hits_checkbox_tpc.isChecked()

        if show_tpc_clusters and show_tpc_hits:
            pass  # Show both
        elif show_tpc_clusters:
            tpc_mask &= is_cluster
        elif show_tpc_hits:
            tpc_mask &= is_hit
        else:
            tpc_mask &= False  # Show neither

        t_cross_val_tpc = self.t_crossing_spin_tpc.value()
        track_id_val_tpc = self.track_id_spin_tpc.value()
        side_data = data.point_data.get("side", None)

        if t_cross_val_tpc != -2000:
            tpc_mask &= t_crossing_data == t_cross_val_tpc
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after t_cross_val_tpc filter ({t_cross_val_tpc}): {np.sum(t_crossing_data == t_cross_val_tpc)}"
            )
        if track_id_val_tpc != -1:
            tpc_mask &= track_id_data == track_id_val_tpc
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: Points after track_id_val_tpc filter ({track_id_val_tpc}): {np.sum(track_id_data == track_id_val_tpc)}"
            )

        # Sides
        if side_data is not None:
            sides_chosen = []
            if self.side0_checkbox_tpc.isChecked():
                sides_chosen.append(0)
            if self.side1_checkbox_tpc.isChecked():
                sides_chosen.append(1)
            if sides_chosen:
                tpc_mask &= np.isin(side_data, sides_chosen)
                print(
                    f"DEBUG {'(fitting)' if for_fitting else ''}: Points after after sides_chosen filter ({sides_chosen}): {np.sum(np.isin(side_data, sides_chosen))}"
                )
            else:
                # If no sides are checked, TPC mask => no points
                tpc_mask &= False
                print("No sides selected for TPC; excluding all TPC points.")

        else:
            print("Missing 'side' attribute in data; skipping side filter.")

        if self.seed_checkbox_tpc.isChecked():
            tpc_mask &= used_in_seed == 1
        if self.track_checkbox_tpc.isChecked():
            tpc_mask &= used_in_track == 1

        combined_mask = combined_mask & (silicon_mask | tpc_mask)
        if not np.any(combined_mask):
            return None

        # 7) Extract points that pass the above filters
        filtered_indices = np.where(combined_mask)[0]
        filtered_points = data.extract_points(filtered_indices)
        if filtered_points is None or filtered_points.n_points == 0:
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: No points passed the filters"
            )
            return None
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Final point count: {filtered_points.n_points}"
        )

        # 8) Optionally apply line-based filter for fitting or if user toggles
        if for_fitting:
            # We only apply line filtering if we have an initial line fit
            if self.line_params_initial:
                distance_mask = np.array(
                    [
                        apply_line_filter(
                            point,
                            self.line_params_initial,
                            self.rphi_window,
                            self.z_window,
                        )
                        for point in filtered_points.points
                    ]
                )
                line_filtered = filtered_points.extract_points(distance_mask)
                return line_filtered if line_filtered.n_points > 0 else None
            else:
                # If we don't have line parameters yet, just return the spatially filtered data
                return filtered_points if filtered_points.n_points > 0 else None
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

    # --- Method to Handle Mode Changes ---
    def on_pick_mode_changed(self):
        """Handle changes in the pick mode based on radio button selection."""
        if self.radio_view_info.isChecked():
            self.pick_mode = "info"
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )

        elif self.radio_pick_helix.isChecked():
            self.pick_mode = "helix"
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )

        elif self.radio_pick_vertex.isChecked():
            self.pick_mode = "vertex"
            self.instruction_label.setText("Select 2 points for vertex finding.")

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
        elif self.pick_mode in ["helix", "vertex"]:
            # Enable point picking for helix fitting
            self.plotter_widget.enable_point_picking(
                callback=self.on_point_picked_helix, show_message=True, use_picker=True
            )

    def on_point_picked_info(self, picked_point, picker):
        """
        Callback function for viewing point information.
        Displays info without affecting helix fitting.
        """
        actor = picker.GetActor()
        # If the actor is one of the track (helix) actors, ignore the pick.
        if actor in self.track_lines:
            return
        # Extract the point ID from the picker
        point_id = picker.GetPointId()

        # Extract the mesh (dataset) from the picker
        mesh = picker.GetDataSet()

        if point_id < 0 or mesh is None:
            return

        # Ensure 'data_type' exists in the mesh's point data
        if "data_type" not in mesh.point_data:
            return  # Ignore picking on meshes without 'data_type'

        # Retrieve the data_type for the picked point
        data_type = mesh.point_data["data_type"][point_id]

        if data_type not in (0, 1):
            return
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
        required_points = (
            2 if (self.pick_mode == "vertex" or self.using_line_fitting) else 3
        )
        self.instruction_label.setText(
            f"Selected {len(self.selected_points_first)}/{required_points} points for initial fitting."
        )
        if len(self.selected_points_first) == required_points:
            msg = (
                "Click 'Fit Helix' to perform fitting."
                if not self.using_line_fitting
                else "Click 'Fit Helix' to perform line fitting."
            )
            QMessageBox.information(self, "Track Fit", msg)

    def on_track_toggle(self, state):
        visible = state == Qt.Checked
        self.set_track_lines_visibility(visible)

    def set_track_lines_visibility(self, visible):
        for actor in self.track_lines:
            # Set the visibility property on the actor.
            # Using VTK’s method: 1 means visible, 0 means hidden.
            actor.SetVisibility(1 if visible else 0)
        self.plotter_widget.render()

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

                if hasattr(cluster_polydata, "hitkeys"):
                    hitkeys = cluster_polydata.hitkeys
                    print("In GUI - Hitkeys verification:")
                    print(f"First cluster's hitkeys: {hitkeys[0]}")
                    print(f"Type of hitkeys data: {type(hitkeys)}")
                    print(f"Length of hitkeys list: {len(hitkeys)}")
                else:
                    print("Warning: No hitkeys found in cluster_polydata")

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

        self.display_filtered_clusters_info.clear()
        self.display_filtered_hits_info.clear()

        for filename, file_info in self.loaded_files.items():
            item = file_info["item"]
            if item.checkState() != Qt.Checked:
                continue

            # --- Clusters ---
            if file_info["cluster"] is not None and file_info["cluster"].n_points > 0:

                filtered_clusters_display = self._filter_data(file_info["cluster"])

                if filtered_clusters_display and filtered_clusters_display.n_points > 0:
                    cluster_color = file_info.get("color", "red")
                    self.display_filtered_clusters_info.append(
                        (filtered_clusters_display, cluster_color, file_info)
                    )
                    self.plotter_widget.add_mesh(
                        filtered_clusters_display,
                        style="points",
                        point_size=5,
                        color=cluster_color,
                    )

            # --- Hits ---
            if file_info["hit"] is not None and file_info["hit"].n_points > 0:
                filtered_hits_display = self._filter_data(file_info["hit"])

                if filtered_hits_display and filtered_hits_display.n_points > 0:
                    self.display_filtered_hits_info.append(
                        (filtered_hits_display, "blue", file_info)
                    )
                    self.plotter_widget.add_mesh(
                        filtered_hits_display,
                        style="points",
                        point_size=5,
                        color="blue",
                    )

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

        self.plotter_widget.show_axes()
        if self.projection_window is not None and self.projection_window.isVisible():
            self.projection_window.update_display(
                clusters_info=self.display_filtered_clusters_info,
                hits_info=self.display_filtered_hits_info,
                helix_lines=self.track_lines,
                projection_type=self.projection_window.projection_combo.currentText(),
            )

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
        if self.pick_mode not in ["helix", "vertex"]:
            QMessageBox.warning(
                self,
                "Track Fit",
                "Please switch to 'Pick for Helix/Line Fitting' mode to fit a track.",
            )
            return

        if not self.selected_points_first:
            QMessageBox.warning(
                self,
                "Track Fit",
                "Please select points for fitting first.",
            )
            return

        selected_file = None
        selected_points_set = set(map(tuple, self.selected_points_first))

        for filename, file_info in self.loaded_files.items():
            if file_info["cluster"] is not None:
                # Convert points to set of tuples for comparison
                file_points_set = set(map(tuple, file_info["cluster"].points))
                # Check if any of the selected points are in this file
                if any(point in file_points_set for point in selected_points_set):
                    selected_file = filename
                    self.cluster_data = file_info["cluster"]
                    print(f"Using clusters from file: {filename}")
                    print(f"Number of clusters: {self.cluster_data.n_points}")
                    break

        if selected_file is None:
            QMessageBox.warning(
                self,
                "Fit Error",
                "Could not find the file containing the selected points.",
            )
            return
        try:
            if self.pick_mode == "vertex":
                # Perform vertex finding and helix fitting
                vertex_params = find_vertex_z(
                    self.selected_points_first[0], self.selected_points_first[1]
                )

                # Create helix parameters from vertex results
                helix_params = {
                    "c_x": vertex_params["c_x"],
                    "c_y": vertex_params["c_y"],
                    "r": vertex_params["r"],
                    "alpha": vertex_params["alpha"],
                    "c_z": vertex_params["c_z"],
                    "t0": vertex_params["ref_theta"],
                    "ref_theta": vertex_params["ref_theta"],
                }

                # Store parameters and generate visualization
                self.helix_params_initial = helix_params
                helix_points = generate_helix_points_initial(
                    helix_params, self.inner_cut, self.outer_cut
                )

                print("Vertex Params:", vertex_params)
                print("Helix Params:", helix_params)
                print(
                    "Helix Points size:",
                    None if helix_points is None else len(helix_points),
                )

                if helix_points is not None:
                    # helix_line = generate_helix_line(helix_points)
                    actor = self.plotter_widget.add_lines(
                        helix_points,
                        color="grey",
                        width=3,
                        # style="wireframe",
                        # pickable=False,
                    )
                    # actor.GetProperty().SetRepresentationToWireframe()
                    # actor.GetProperty().SetPointSize(0)
                    actor.GetProperty().SetRepresentationToWireframe()
                    actor.GetProperty().SetPointSize(0)
                    actor.PickableOff()
                    self.track_lines.append(actor)

                    vertex_pca = np.array([0.0, 0.0, vertex_params["vertex_z"]])

                    # Display vertex information
                    info_text = "Vertex Finding Results:\n"
                    info_text += f"Vertex Z: {vertex_params['vertex_z']:.3f} cm\n"
                    info_text += f"Helix radius: {vertex_params['r']:.3f} cm\n"
                    info_text += f"Helix center: ({vertex_params['c_x']:.3f}, {vertex_params['c_y']:.3f}) cm\n"
                    info_text += (
                        f"Pitch parameter (alpha): {vertex_params['alpha']:.3f}\n"
                    )
                    self.info_panel.setText(info_text)

                    self.track_counter += 1
                    if not self.histogram_window.isVisible():
                        self.histogram_window.show()
                    self.histogram_window.add_histograms(
                        [],
                        [],
                        self.track_counter,
                        points=None,
                        pca=vertex_pca,
                    )

            elif self.using_line_fitting:
                # Perform line fitting
                self.do_line_fitting_flow()

            else:
                # Perform regular helix fitting
                self.do_helix_fitting_flow()

        except ValueError as e:
            QMessageBox.warning(self, "Fit Error", str(e))
        except Exception as e:
            QMessageBox.critical(
                self, "Error", f"An unexpected error occurred:\n{str(e)}"
            )

        self.selected_points_first.clear()
        self.update_display()
        self.plotter_widget.disable_picking()
        self.update_point_picking()

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

        pca_line, pca_beam, dca = calculate_dca_to_beam(line_params_init)
        info_text = self.info_panel.toPlainText()
        info_text += "\n\nDCA Analysis:"
        info_text += f"\nDCA to beam axis: {dca:.3f} cm"
        if pca_line is not None:
            info_text += f"\nPCA on track: ({pca_line[0]:.3f}, {pca_line[1]:.3f}, {pca_line[2]:.3f}) cm"
            info_text += f"\nPCA on beam: ({pca_beam[0]:.3f}, {pca_beam[1]:.3f}, {pca_beam[2]:.3f}) cm"
        self.info_panel.setText(info_text)
        centroid_initial = compute_centroid(self.selected_points_first)
        line_pca = pca_line if pca_line is not None else None

        self.line_params_initial = line_params_init
        # Visualize the initial line
        line_points_init = generate_line_points(
            line_params_init, self.inner_cut, self.outer_cut, centroid_initial
        )
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
        filtered_clusters_for_fitting = self._filter_data(
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
        if (
            filtered_clusters_for_fitting is not None
            and filtered_clusters_for_fitting.n_points > 0
        ):
            # Calculate line residuals
            delta_rphi, delta_z, valid_points = calculate_deltas_line(
                self.line_params_refined, filtered_clusters_for_fitting
            )

            # Increment track counter and update histograms
            self.track_counter += 1
            if not self.histogram_window.isVisible():
                self.histogram_window.show()
            self.histogram_window.add_histograms(
                delta_rphi,
                delta_z,
                self.track_counter,
                points=valid_points,
                pca=line_pca,
            )
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

            self.helix_line_color = "grey"

            # Generate refined helix points and create a tube
            helix_points_refined = generate_helix_points_refined(
                helix_params_refined, self.inner_cut, self.outer_cut
            )
            print(f"Helix Points Refined: {helix_points_refined}")
            self.helix_points = (
                helix_points_refined  # Update helix points for distance calculations
            )

            helix_line_refined = generate_helix_line(helix_points_refined)
            print(f"Helix Line Refined: {helix_line_refined}")
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

                dca_result = find_dca_and_closest_point(
                    helix_params_refined, helix_params_refined["ref_theta_direct"]
                )
                if dca_result is not None:
                    dca_val = dca_result["dca"]
                    closest_point_3d = dca_result["closest_point_3d"]
                    z_dca = closest_point_3d[2]  # Extract the z-coordinate

                    # Display DCA and Z at DCA in the info panel
                    info_text = self.info_panel.toPlainText()
                    info_text += (
                        f"\nRefined Helix Fit:\n"
                        f"DCA to Beam Axis = {dca_val:.3f} cm\n"
                        f"Z at DCA        = {z_dca:.3f} cm\n"
                    )
                    self.info_panel.setText(info_text)

                else:
                    QMessageBox.warning(
                        self,
                        "DCA Calculation",
                        "Failed to calculate DCA after helix fitting.",
                    )

            else:
                QMessageBox.warning(
                    self, "Helix Fit", "Failed to create refined helix visualization."
                )

            if (
                filtered_clusters_for_fitting is not None
                and filtered_clusters_for_fitting.n_points > 0
            ):
                delta_rphi, delta_z, valid_points = calculate_deltas(
                    self.helix_params_refined, filtered_clusters_for_fitting
                )
                self.track_counter += 1
                if not self.histogram_window.isVisible():
                    self.histogram_window.show()
                self.histogram_window.add_histograms(
                    delta_rphi,
                    delta_z,
                    self.track_counter,
                    points=valid_points,
                    pca=closest_point_3d,
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
                info_text = self.info_panel.toPlainText()
                info_text += (
                    f"\nTrack ID: {track_id}\n"
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

    def reset_track(self):
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

    def clear_track_lines(self):
        """Clears all displayed helix lines from the plot."""
        # Clear the list of helix lines
        self.track_lines = []
        # Optionally update the display to refresh the view
        self.update_display()

    def show_projection_window(self):
        """Handle showing the 2D projection window."""
        if self.projection_window is None:
            self.projection_window = GeometricProjectionWindow(self)

        self.projection_window.update_display(
            clusters_info=self.display_filtered_clusters_info,  # Pass all clusters with their colors
            hits_info=self.display_filtered_hits_info,  # Pass all hits
            helix_lines=self.track_lines,
            projection_type="XY Projection",
        )

        # Show the window
        self.projection_window.raise_()
        self.projection_window.activateWindow()
        self.projection_window.show()

    def toggle_field_mode(self):
        """Toggle between helix and straight line fitting modes."""
        self.using_line_fitting = not self.using_line_fitting

        # Update button text
        if self.using_line_fitting:
            self.btn_toggle_field.setText("Field Off")
            self.instruction_label.setText(
                "Instruction: Select 2 points for initial line fitting."
            )
        else:
            self.btn_toggle_field.setText("Field On")
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )

        # Clear any existing selections
        self.selected_points_first.clear()
        self.update_display()
