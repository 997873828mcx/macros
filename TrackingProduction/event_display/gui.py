import sys
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from pad_adc_window import PadADCWindow
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
    QDockWidget,
    QToolBar,
    QAction,
    QStatusBar,
    QTabWidget,
    QScrollArea,
    QGridLayout,
    QFormLayout,
    QMenu,
    QSizePolicy,
)
from qtpy.QtCore import Qt, QSize
from qtpy.QtGui import QIcon, QStandardItem, QStandardItemModel


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

        self.cluster_actors = []
        self.hit_actors = []
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
        self.pad_adc_window = PadADCWindow()
        self.projection_window = None
        self.avg_position_actor = None

        self.associated_hits_hidden = False  # False means show global hits normally.
        self.hits_override_enabled = False
        self.previous_hits_state = {"tpc_hits": True, "silicon_hits": True}
        self.selected_cluster_index = None
        self.selected_points_first = []
        self.last_picked_mesh = None
        self.setup_ui()

        # Data placeholders
        self.cluster_data = None
        self.hit_data = None
        self.cluster_polydata = None
        self.hit_polydata = None

        # self.selected_points_first = []  # Used for point selection
        self.helix_points = None
        self.helix_params_initial = None
        self.helix_params_refined = None
        self.line_params_initial = None
        self.line_params_refined = None

        self.filtered_clusters = None

        self.plotter_widget.show_axes()

        self.update_point_picking()

    def setup_ui(self):
        """Set up the redesigned UI with dockable panels and maximized 3D view."""
        # Central widget - 3D Visualization Area
        central_widget = QWidget()
        central_layout = QVBoxLayout(central_widget)
        central_layout.setContentsMargins(0, 0, 0, 0)

        # Create 3D visualization widget
        self.plotter_widget = QtInteractor()
        central_layout.addWidget(self.plotter_widget)

        # Add information panel at the bottom
        self.info_panel = QTextEdit()
        self.info_panel.setReadOnly(True)
        self.info_panel.setMaximumHeight(100)
        central_layout.addWidget(self.info_panel)

        self.setCentralWidget(central_widget)

        # Create dock widgets for controls
        self.create_sidebar()  # <-- Change this line from create_dock_widgets() to create_sidebar()

        # Create toolbars
        self.create_toolbars()

        # Create status bar
        self.create_status_bar()

        # Initial status message
        self.instruction_label.setText(
            "Instruction: Select 3 points for initial track fitting."
        )

    def create_toolbars(self):
        """Create main toolbar with commonly used actions."""
        main_toolbar = QToolBar("Main Toolbar")
        main_toolbar.setIconSize(QSize(16, 16))
        self.addToolBar(Qt.TopToolBarArea, main_toolbar)

        # Left side toolbar elements
        self.setWindowTitle("Event Display")
        title_label = QLabel("Event Display")
        title_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        main_toolbar.addWidget(title_label)

        main_toolbar.addSeparator()

        # Load ROOT File button
        load_button = QPushButton("Load ROOT")
        load_button.clicked.connect(self.load_data)
        main_toolbar.addWidget(load_button)

        # Pick Mode dropdown
        pick_mode_button = QPushButton("Pick Mode ▾")
        pick_mode_menu = QMenu(self)

        info_action = QAction("View Point Info", self)
        info_action.setCheckable(True)
        info_action.setChecked(True)
        info_action.triggered.connect(lambda: self.set_pick_mode("info"))

        helix_action = QAction("Pick for Helix Fitting", self)
        helix_action.setCheckable(True)
        helix_action.triggered.connect(lambda: self.set_pick_mode("helix"))

        vertex_action = QAction("Pick for Vertex Finding", self)
        vertex_action.setCheckable(True)
        vertex_action.triggered.connect(lambda: self.set_pick_mode("vertex"))

        pos_action = QAction("Pos Reco", self)
        pos_action.setCheckable(True)
        pos_action.triggered.connect(lambda: self.set_pick_mode("pos"))

        self.pick_mode_actions = [info_action, helix_action, vertex_action, pos_action]

        for action in self.pick_mode_actions:
            pick_mode_menu.addAction(action)

        pick_mode_button.setMenu(pick_mode_menu)
        main_toolbar.addWidget(pick_mode_button)

        # Fit Track button
        fit_track_button = QPushButton("Fit Track")
        fit_track_button.clicked.connect(self.initiate_fit)
        main_toolbar.addWidget(fit_track_button)

        # Settings button
        settings_button = QPushButton("Settings")
        settings_menu = QMenu(self)

        self.toggle_field_action = QAction("Field Mode: On", self)
        self.toggle_field_action.triggered.connect(self.toggle_field_mode)

        reset_track_action = QAction("Reset Track", self)
        reset_track_action.triggered.connect(self.reset_track)

        clear_tracks_action = QAction("Clear Track Lines", self)
        clear_tracks_action.triggered.connect(self.clear_track_lines)

        settings_menu.addAction(self.toggle_field_action)
        settings_menu.addAction(reset_track_action)
        settings_menu.addAction(clear_tracks_action)

        settings_button.setMenu(settings_menu)
        main_toolbar.addWidget(settings_button)

        # Add spacer to push event selector to right
        spacer = QWidget()
        size_policy = QSizePolicy()
        size_policy.setHorizontalPolicy(QSizePolicy.Expanding)
        size_policy.setVerticalPolicy(QSizePolicy.Preferred)
        spacer.setSizePolicy(size_policy)
        main_toolbar.addWidget(spacer)

        file_label = QLabel("File:")
        main_toolbar.addWidget(file_label)

        self.file_combo = QComboBox()
        self.file_combo.setMinimumWidth(120)
        model = QStandardItemModel(self.file_combo)
        self.file_combo.setModel(model)
        model.itemChanged.connect(lambda item: self.update_display())
        main_toolbar.addWidget(self.file_combo)

        # Then add a separator
        main_toolbar.addSeparator()
        # Event selector on right
        event_label = QLabel("Event:")
        main_toolbar.addWidget(event_label)

        self.event_combo = QComboBox()
        self.event_combo.setMinimumWidth(120)
        event_model = QStandardItemModel(self.event_combo)
        self.event_combo.setModel(event_model)

        for event_num in range(31):  # 0 to 30
            item = QStandardItem(f"Event {event_num}")
            item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
            item.setData(Qt.Unchecked, Qt.CheckStateRole)
            event_model.appendRow(item)

        event_model.itemChanged.connect(lambda i: self.update_display())
        main_toolbar.addWidget(self.event_combo)

    def create_sidebar(self):
        """Create left sidebar with tabs for all control groups."""
        # Create dock widget for the sidebar
        self.sidebar_dock = QDockWidget("Controls", self)
        self.sidebar_dock.setAllowedAreas(Qt.LeftDockWidgetArea)
        self.sidebar_dock.setFeatures(QDockWidget.DockWidgetClosable)
        self.sidebar_dock.setFixedWidth(230)  # Set fixed width for the sidebar

        # Create tab widget for the sidebar
        self.sidebar_tabs = QTabWidget()
        self.sidebar_tabs.setTabPosition(QTabWidget.West)  # Tabs on the left

        # Create the four tabs with their respective widgets
        self.create_display_tab()
        self.create_filters_tab()
        self.create_analysis_tab()
        self.create_tracking_tab()

        # Set the tab widget as the dock widget's content
        self.sidebar_dock.setWidget(self.sidebar_tabs)

        # Add dock widget to left area
        self.addDockWidget(Qt.LeftDockWidgetArea, self.sidebar_dock)

    def create_display_tab(self):
        """Create the Display tab content."""
        display_widget = QWidget()
        display_layout = QVBoxLayout(display_widget)
        display_layout.setContentsMargins(8, 8, 8, 8)
        display_layout.setSpacing(8)

        # Silicon section
        silicon_group = QGroupBox("Silicon")
        silicon_layout = QGridLayout()

        # Clusters and hits checkboxes
        self.clusters_checkbox_silicon = QCheckBox("Clusters")
        self.clusters_checkbox_silicon.setChecked(True)
        self.clusters_checkbox_silicon.stateChanged.connect(self.update_display)
        silicon_layout.addWidget(self.clusters_checkbox_silicon, 0, 0)

        self.seed_checkbox_silicon = QCheckBox("Seed")
        self.seed_checkbox_silicon.setChecked(False)
        self.seed_checkbox_silicon.stateChanged.connect(self.update_display)
        silicon_layout.addWidget(self.seed_checkbox_silicon, 0, 1)

        self.hits_checkbox_silicon = QCheckBox("Hits")
        self.hits_checkbox_silicon.setChecked(True)
        self.hits_checkbox_silicon.stateChanged.connect(self.update_display)
        silicon_layout.addWidget(self.hits_checkbox_silicon, 1, 0)

        self.track_checkbox_silicon = QCheckBox("Track")
        self.track_checkbox_silicon.setChecked(False)
        self.track_checkbox_silicon.stateChanged.connect(self.update_display)
        silicon_layout.addWidget(self.track_checkbox_silicon, 1, 1)

        silicon_group.setLayout(silicon_layout)
        display_layout.addWidget(silicon_group)

        # TPC section
        tpc_group = QGroupBox("TPC")
        tpc_layout = QGridLayout()

        # Clusters and hits checkboxes
        self.clusters_checkbox_tpc = QCheckBox("Clusters")
        self.clusters_checkbox_tpc.setChecked(True)
        self.clusters_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.clusters_checkbox_tpc, 0, 0)

        self.seed_checkbox_tpc = QCheckBox("Seed")
        self.seed_checkbox_tpc.setChecked(False)
        self.seed_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.seed_checkbox_tpc, 0, 1)

        self.hits_checkbox_tpc = QCheckBox("Hits")
        self.hits_checkbox_tpc.setChecked(True)
        self.hits_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.hits_checkbox_tpc, 1, 0)

        self.track_checkbox_tpc = QCheckBox("Track")
        self.track_checkbox_tpc.setChecked(False)
        self.track_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.track_checkbox_tpc, 1, 1)

        self.side0_checkbox_tpc = QCheckBox("Side 0")
        self.side0_checkbox_tpc.setChecked(True)
        self.side0_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.side0_checkbox_tpc, 2, 0)

        self.side1_checkbox_tpc = QCheckBox("Side 1")
        self.side1_checkbox_tpc.setChecked(True)
        self.side1_checkbox_tpc.stateChanged.connect(self.update_display)
        tpc_layout.addWidget(self.side1_checkbox_tpc, 2, 1)

        tpc_group.setLayout(tpc_layout)
        display_layout.addWidget(tpc_group)

        # Thresholds
        threshold_group = QGroupBox("Thresholds")
        threshold_layout = QFormLayout()

        self.cluster_adc_spinbox = QDoubleSpinBox()
        self.cluster_adc_spinbox.setRange(0, 10000)
        self.cluster_adc_spinbox.setValue(0)
        self.cluster_adc_spinbox.setDecimals(0)
        self.cluster_adc_spinbox.valueChanged.connect(self.update_display)
        threshold_layout.addRow("Cluster ADC ≥", self.cluster_adc_spinbox)

        self.hit_adc_spinbox = QDoubleSpinBox()
        self.hit_adc_spinbox.setRange(0, 10000)
        self.hit_adc_spinbox.setValue(0)
        self.hit_adc_spinbox.setDecimals(0)
        self.hit_adc_spinbox.valueChanged.connect(self.update_display)
        threshold_layout.addRow("Hit ADC ≥", self.hit_adc_spinbox)

        threshold_group.setLayout(threshold_layout)
        display_layout.addWidget(threshold_group)

        # Axis Range
        axis_group = QGroupBox("Axis Range")
        axis_layout = QVBoxLayout()

        # X Range
        x_layout = QHBoxLayout()  # Change to horizontal layout
        x_label = QLabel("X:")
        x_layout.addWidget(x_label)

        # Create two spinboxes instead of a slider
        self.x_min_spinbox = QSpinBox()
        self.x_min_spinbox.setRange(-150, 150)
        self.x_min_spinbox.setValue(-100)
        self.x_min_spinbox.valueChanged.connect(self.update_display)
        x_layout.addWidget(self.x_min_spinbox)

        x_layout.addWidget(QLabel("to"))  # Add a label between spinboxes

        self.x_max_spinbox = QSpinBox()
        self.x_max_spinbox.setRange(-150, 150)
        self.x_max_spinbox.setValue(100)
        self.x_max_spinbox.valueChanged.connect(self.update_display)
        x_layout.addWidget(self.x_max_spinbox)

        axis_layout.addLayout(x_layout)

        y_layout = QHBoxLayout()  # Change to horizontal layout
        y_label = QLabel("Y:")
        y_layout.addWidget(y_label)

        # Create two spinboxes instead of a slider
        self.y_min_spinbox = QSpinBox()
        self.y_min_spinbox.setRange(-150, 150)
        self.y_min_spinbox.setValue(-100)
        self.y_min_spinbox.valueChanged.connect(self.update_display)
        y_layout.addWidget(self.y_min_spinbox)

        y_layout.addWidget(QLabel("to"))  # Add a label between spinboxes

        self.y_max_spinbox = QSpinBox()
        self.y_max_spinbox.setRange(-150, 150)
        self.y_max_spinbox.setValue(100)
        self.y_max_spinbox.valueChanged.connect(self.update_display)
        y_layout.addWidget(self.y_max_spinbox)

        axis_layout.addLayout(y_layout)

        z_layout = QHBoxLayout()  # Change to horizontal layout
        z_label = QLabel("Z:")
        z_layout.addWidget(z_label)

        # Create two spinboxes instead of a slider
        self.z_min_spinbox = QSpinBox()
        self.z_min_spinbox.setRange(-300, 300)
        self.z_min_spinbox.setValue(-300)
        self.z_min_spinbox.valueChanged.connect(self.update_display)
        z_layout.addWidget(self.z_min_spinbox)

        z_layout.addWidget(QLabel("to"))  # Add a label between spinboxes

        self.z_max_spinbox = QSpinBox()
        self.z_max_spinbox.setRange(-300, 300)
        self.z_max_spinbox.setValue(300)
        self.z_max_spinbox.valueChanged.connect(self.update_display)
        z_layout.addWidget(self.z_max_spinbox)

        axis_layout.addLayout(z_layout)

        axis_group.setLayout(axis_layout)
        display_layout.addWidget(axis_group)

        # Add stretch to push everything to the top
        display_layout.addStretch()

        # Add the tab
        self.sidebar_tabs.addTab(display_widget, "D")

    def create_filters_tab(self):
        """Create the Filters tab content."""
        filters_widget = QWidget()
        filters_layout = QVBoxLayout(filters_widget)
        filters_layout.setContentsMargins(8, 8, 8, 8)
        filters_layout.setSpacing(8)

        # Silicon ID filters
        silicon_group = QGroupBox("Silicon ID Filters")
        silicon_form = QFormLayout()

        # Crossing
        self.crossing_spin_silicon = QSpinBox()
        self.crossing_spin_silicon.setRange(-2000, 9999)
        self.crossing_spin_silicon.setValue(-2000)
        self.crossing_spin_silicon.valueChanged.connect(self.update_display)
        silicon_form.addRow("Crossing:", self.crossing_spin_silicon)

        # T-Crossing
        self.t_crossing_spin_silicon = QSpinBox()
        self.t_crossing_spin_silicon.setRange(-2000, 99999)
        self.t_crossing_spin_silicon.setValue(-2000)
        self.t_crossing_spin_silicon.valueChanged.connect(self.update_display)
        silicon_form.addRow("T-Crossing:", self.t_crossing_spin_silicon)

        # Track ID
        self.track_id_spin_silicon = QSpinBox()
        self.track_id_spin_silicon.setRange(-1, 999999)
        self.track_id_spin_silicon.setValue(-1)
        self.track_id_spin_silicon.valueChanged.connect(self.update_display)
        silicon_form.addRow("Track ID:", self.track_id_spin_silicon)

        silicon_group.setLayout(silicon_form)
        filters_layout.addWidget(silicon_group)

        # TPC ID filters
        tpc_group = QGroupBox("TPC ID Filters")
        tpc_form = QFormLayout()

        # T-Crossing
        self.t_crossing_spin_tpc = QSpinBox()
        self.t_crossing_spin_tpc.setRange(-2000, 99999)
        self.t_crossing_spin_tpc.setValue(-2000)
        self.t_crossing_spin_tpc.valueChanged.connect(self.update_display)
        tpc_form.addRow("T-Crossing:", self.t_crossing_spin_tpc)

        # Track ID
        self.track_id_spin_tpc = QSpinBox()
        self.track_id_spin_tpc.setRange(-1, 999999)
        self.track_id_spin_tpc.setValue(-1)
        self.track_id_spin_tpc.valueChanged.connect(self.update_display)
        tpc_form.addRow("Track ID:", self.track_id_spin_tpc)

        tpc_group.setLayout(tpc_form)
        filters_layout.addWidget(tpc_group)

        # Additional filters
        other_group = QGroupBox("Other Filters")
        other_layout = QVBoxLayout()

        self.nmaps_spinbox = QSpinBox()
        self.nmaps_spinbox.setRange(-10, 100)
        self.nmaps_spinbox.setValue(-1)
        self.nmaps_spinbox.valueChanged.connect(self.update_display)
        nmaps_layout = QFormLayout()
        nmaps_layout.addRow("NMAPS ≥", self.nmaps_spinbox)
        other_layout.addLayout(nmaps_layout)

        # Show points outside tube
        self.toggle_filter_checkbox = QCheckBox("Show Points Outside Tube")
        self.toggle_filter_checkbox.setChecked(True)
        self.toggle_filter_checkbox.stateChanged.connect(self.update_display)
        other_layout.addWidget(self.toggle_filter_checkbox)

        other_group.setLayout(other_layout)
        filters_layout.addWidget(other_group)

        # Add stretch to push everything to the top
        filters_layout.addStretch()

        # Add the tab
        self.sidebar_tabs.addTab(filters_widget, "F")

    def create_analysis_tab(self):
        """Create the Analysis tab content."""
        analysis_widget = QWidget()
        analysis_layout = QVBoxLayout(analysis_widget)
        analysis_layout.setContentsMargins(8, 8, 8, 8)
        analysis_layout.setSpacing(8)

        # Associated hits analysis
        hits_group = QGroupBox("Associated Hits")
        hits_layout = QVBoxLayout()

        self.btn_toggle_associated_hits = QPushButton("Show Associated Hits")
        self.btn_toggle_associated_hits.clicked.connect(self.toggle_associated_hits)
        hits_layout.addWidget(self.btn_toggle_associated_hits)

        self.btn_output_phi_z = QPushButton("ADC Weighted Phi/Z")
        self.btn_output_phi_z.clicked.connect(self.output_adc_averaged_phi_z)
        hits_layout.addWidget(self.btn_output_phi_z)

        self.btn_show_adc_dist = QPushButton("Show ADC Distribution")
        self.btn_show_adc_dist.clicked.connect(self.show_pad_adc_window)
        hits_layout.addWidget(self.btn_show_adc_dist)

        hits_group.setLayout(hits_layout)
        analysis_layout.addWidget(hits_group)

        # Projection controls
        projection_group = QGroupBox("Projection")
        projection_layout = QVBoxLayout()

        self.btn_show_projection = QPushButton("Show 2D Projection")
        self.btn_show_projection.clicked.connect(self.show_projection_window)
        projection_layout.addWidget(self.btn_show_projection)

        projection_group.setLayout(projection_layout)
        analysis_layout.addWidget(projection_group)

        # Add stretch to push everything to the top
        analysis_layout.addStretch()

        # Add the tab
        self.sidebar_tabs.addTab(analysis_widget, "A")

    def create_tracking_tab(self):
        """Create the Tracking tab content."""
        tracking_widget = QWidget()
        tracking_layout = QVBoxLayout(tracking_widget)
        tracking_layout.setContentsMargins(8, 8, 8, 8)
        tracking_layout.setSpacing(8)

        # Track fitting group
        fitting_group = QGroupBox("Track Fitting")
        fitting_layout = QVBoxLayout()

        # Fitting buttons
        button_layout = QVBoxLayout()

        self.btn_fit_track = QPushButton("Fit Track")
        self.btn_fit_track.clicked.connect(self.initiate_fit)
        button_layout.addWidget(self.btn_fit_track)

        self.btn_reset_track = QPushButton("Reset Track")
        self.btn_reset_track.clicked.connect(self.reset_track)
        button_layout.addWidget(self.btn_reset_track)

        self.btn_clear_track = QPushButton("Clear Track Lines")
        self.btn_clear_track.clicked.connect(self.clear_track_lines)
        button_layout.addWidget(self.btn_clear_track)

        fitting_layout.addLayout(button_layout)

        # Track options
        options_layout = QVBoxLayout()

        # Module fitting
        self.module_based_fitting_checkbox = QCheckBox("Module Fitting")
        self.module_based_fitting_checkbox.setChecked(False)
        options_layout.addWidget(self.module_based_fitting_checkbox)

        # Show track lines
        self.track_toggle_checkbox = QCheckBox("Show Lines")
        self.track_toggle_checkbox.setChecked(True)
        self.track_toggle_checkbox.stateChanged.connect(self.on_track_toggle)
        options_layout.addWidget(self.track_toggle_checkbox)

        fitting_layout.addLayout(options_layout)

        fitting_group.setLayout(fitting_layout)
        tracking_layout.addWidget(fitting_group)

        # Picking mode group (moved from toolbar to this tab)
        pick_group = QGroupBox("Pick Mode")
        pick_layout = QVBoxLayout()

        # Radio buttons for pick modes
        self.radio_view_info = QRadioButton("View Point Info")
        self.radio_view_info.setChecked(True)
        self.radio_view_info.toggled.connect(self.on_pick_mode_changed)
        pick_layout.addWidget(self.radio_view_info)

        self.radio_pick_helix = QRadioButton("Pick for Helix Fitting")
        self.radio_pick_helix.toggled.connect(self.on_pick_mode_changed)
        pick_layout.addWidget(self.radio_pick_helix)

        self.radio_pick_vertex = QRadioButton("Pick for Vertex Finding")
        self.radio_pick_vertex.toggled.connect(self.on_pick_mode_changed)
        pick_layout.addWidget(self.radio_pick_vertex)

        self.radio_pick_pos = QRadioButton("Pos Reco")
        self.radio_pick_pos.toggled.connect(self.on_pick_mode_changed)
        pick_layout.addWidget(self.radio_pick_pos)

        pick_group.setLayout(pick_layout)
        tracking_layout.addWidget(pick_group)

        # Add stretch to push everything to the top
        tracking_layout.addStretch()

        # Add the tab
        self.sidebar_tabs.addTab(tracking_widget, "T")

    def create_status_bar(self):
        """Create status bar with run info and instruction label."""
        status_bar = self.statusBar()

        # Run information
        run_info = QLabel(
            "<b>Run:</b> 53217 | <b>ZDC coincidence:</b> Raw: 1,869750 | Live: 1,868219"
        )
        status_bar.addWidget(run_info)

        # Add a spacer
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        status_bar.addWidget(spacer)

        # Instruction label
        self.instruction_label = QLabel()
        status_bar.addWidget(self.instruction_label)

    def set_pick_mode(self, mode):
        """Set the pick mode from toolbar actions."""
        self.pick_mode = mode

        # Update radio buttons on Tracking tab
        if mode == "info":
            self.radio_view_info.setChecked(True)
        elif mode == "helix":
            self.radio_pick_helix.setChecked(True)
        elif mode == "vertex":
            self.radio_pick_vertex.setChecked(True)
        elif mode == "pos":
            self.radio_pick_pos.setChecked(True)

        # Update instruction label
        if mode == "info":
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )
        elif mode == "helix":
            self.instruction_label.setText(
                "Instruction: Select 3 points for initial helix fitting."
            )
        elif mode == "vertex":
            self.instruction_label.setText("Select 2 points for vertex finding.")
        elif mode == "pos":
            self.instruction_label.setText("Select hits for position reconstruction.")

        # Update the point picking callback
        self.update_point_picking()

        # Update toolbar menu actions
        for action in self.pick_mode_actions:
            action.setChecked(action.text().lower().find(mode.lower()) >= 0)

    def _filter_data(
        self,
        data: pv.PolyData,
        *,
        for_fitting: bool = False,
    ) -> Optional[pv.PolyData]:

        if data is None or data.n_points == 0:

            return None

        points = data.points

        event_numbers = data.point_data.get("event", None)
        if event_numbers is None:

            return None

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

        if not selected_events:

            return None
        event_mask = np.isin(event_numbers, selected_events)

        x_min, x_max = self.x_min_spinbox.value(), self.x_max_spinbox.value()
        y_min, y_max = self.y_min_spinbox.value(), self.y_max_spinbox.value()
        z_min, z_max = self.z_min_spinbox.value(), self.z_max_spinbox.value()

        spatial_mask = (
            (points[:, 0] >= x_min)
            & (points[:, 0] <= x_max)
            & (points[:, 1] >= y_min)
            & (points[:, 1] <= y_max)
            & (points[:, 2] >= z_min)
            & (points[:, 2] <= z_max)
        )

        combined_mask = spatial_mask & event_mask

        nmaps_values = data.point_data.get("nmaps", None)

        if nmaps_values is not None:
            nmaps_threshold = self.nmaps_spinbox.value()
            nmaps_mask = nmaps_values >= nmaps_threshold
            combined_mask &= nmaps_mask

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

        if (
            hasattr(self, "hits_override_enabled")
            and self.hits_override_enabled
            and self.associated_hits_hidden
        ):
            # When in associated hits mode, force-enable hits display regardless of checkbox
            show_silicon_hits = True
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

        if t_cross_val_sil != -2000:
            silicon_mask &= t_crossing_data == t_cross_val_sil

        if track_id_val_sil != -1:
            silicon_mask &= track_id_data == track_id_val_sil

        if self.seed_checkbox_silicon.isChecked():
            silicon_mask &= used_in_seed == 1

        if self.track_checkbox_silicon.isChecked():
            silicon_mask &= used_in_track == 1

        is_tpc = layer_array >= 7

        tpc_mask = is_tpc
        show_tpc_clusters = self.clusters_checkbox_tpc.isChecked()
        show_tpc_hits = self.hits_checkbox_tpc.isChecked()

        if (
            hasattr(self, "hits_override_enabled")
            and self.hits_override_enabled
            and self.associated_hits_hidden
        ):
            # When in associated hits mode, force-enable hits display regardless of checkbox
            show_tpc_hits = True

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

        if track_id_val_tpc != -1:
            tpc_mask &= track_id_data == track_id_val_tpc

        # Sides
        if side_data is not None:
            sides_chosen = []
            if self.side0_checkbox_tpc.isChecked():
                sides_chosen.append(0)
            if self.side1_checkbox_tpc.isChecked():
                sides_chosen.append(1)
            if sides_chosen:
                tpc_mask &= np.isin(side_data, sides_chosen)

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
        if hasattr(data, "hitkeys"):
            # Verify that hitkeys list matches the number of points
            if len(data.hitkeys) == data.n_points:
                filtered_points.hitkeys = [data.hitkeys[i] for i in filtered_indices]
            else:
                print(
                    f"Warning: hitkeys size ({len(data.hitkeys)}) doesn't match number of points ({data.n_points})"
                )

        if filtered_points is None or filtered_points.n_points == 0:

            return None

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

                    continue

        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Selected events: {selected_events}"
        )

        if not selected_events:
            print(
                f"DEBUG {'(fitting)' if for_fitting else ''}: No events selected, returning None"
            )
            return None
        event_mask = np.isin(event_numbers, selected_events)
        print(
            f"DEBUG {'(fitting)' if for_fitting else ''}: Points after event filter: {np.sum(event_mask)}"
        )

        x_min, x_max = self.x_min_spinbox.value(), self.x_max_spinbox.value()
        y_min, y_max = self.y_min_spinbox.value(), self.y_max_spinbox.value()
        z_min, z_max = self.z_min_spinbox.value(), self.z_max_spinbox.value()
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

        elif self.radio_pick_pos.isChecked():
            self.pick_mode = "pos"
            self.instruction_label.setText("Select hits for position reconstruction.")

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
        elif self.pick_mode in ["helix", "vertex", "pos"]:
            # Enable point picking for helix fitting and position reconstruction
            if self.pick_mode == "pos":
                callback = self.on_point_picked_pos
            else:
                callback = self.on_point_picked_helix

            self.plotter_widget.enable_point_picking(
                callback=callback, show_message=True, use_picker=True
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

        if data_type == 0:  # Cluster
            self.last_picked_mesh = mesh
            self.selected_cluster_index = point_id
            self.info_panel.setText(f"Selected Cluster (Point ID: {point_id})")

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

    def on_point_picked_pos(self, picked_point, picker):
        """
        Callback function for picking points in position reconstruction mode.
        """
        # Extract the point ID from the picker
        point_id = picker.GetPointId()

        # Extract the mesh (dataset) from the picker
        mesh = picker.GetDataSet()

        if point_id < 0 or mesh is None:
            return  # No valid point was picked

        # Get the picked coordinates
        coordinates = mesh.points[point_id]

        hit_info = {
            "x": coordinates[0],
            "y": coordinates[1],
            "z": coordinates[2],
            "t": mesh.point_data["t"][point_id],
            "adc": mesh.point_data["adc"][point_id],
            "tdriftmax": mesh.point_data["tdriftmax"][point_id],
            "drift_velocity": mesh.point_data["driftVelocity"][point_id],
            "phi": mesh.point_data["phi"][point_id],
            "pad": mesh.point_data["pad"][point_id] if "pad" in mesh.point_data else 0,
            "time": mesh.point_data["tbin"][point_id],
        }
        self.selected_points_first.append(hit_info)

        # Update instruction label with current count
        self.instruction_label.setText(
            f"Selected {len(self.selected_points_first)} hits for position reconstruction."
        )

        if self.pad_adc_window.isVisible():
            self.pad_adc_window.update_distribution(self.selected_points_first)

        # Display point info
        info_text = f"Selected Hit (Point ID: {point_id})\n"
        info_text += f"Position: ({coordinates[0]:.2f}, {coordinates[1]:.2f}, {coordinates[2]:.2f})\n"
        for attr in mesh.point_data.keys():
            value = mesh.point_data[attr][point_id]
            info_text += f"{attr}: {value}\n"
        self.info_panel.setText(info_text)

    def on_track_toggle(self, state):
        visible = state == Qt.Checked
        self.set_track_lines_visibility(visible)

    def set_track_lines_visibility(self, visible):
        for actor in self.track_lines:
            # Set the visibility property on the actor.
            # Using VTK's method: 1 means visible, 0 means hidden.
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
                print("In GUI - Hitkeys verification:")
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
        saved_avg_position_actor = None
        if hasattr(self, "avg_position_actor") and self.avg_position_actor is not None:
            saved_avg_position_actor = self.avg_position_actor
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

                    if (
                        self.associated_hits_hidden
                        and hasattr(self, "last_picked_mesh")
                        and self.last_picked_mesh is not None
                        and self.selected_cluster_index is not None
                    ):
                        try:
                            # Get cluster information
                            selected_cluster_event = self.last_picked_mesh.point_data[
                                "event"
                            ][self.selected_cluster_index]
                            selected_cluster_keys = self.last_picked_mesh.hitkeys[
                                self.selected_cluster_index
                            ]

                            # Get hitsetkey if available
                            if "hitsetkey" in self.last_picked_mesh.point_data:
                                selected_cluster_hitsetkey = (
                                    self.last_picked_mesh.point_data["hitsetkey"][
                                        self.selected_cluster_index
                                    ]
                                )
                            else:
                                selected_cluster_hitsetkey = None

                            # Get hit information
                            all_hitkeys = filtered_hits_display.point_data["hitkeykey"]
                            all_events = filtered_hits_display.point_data["event"]

                            # Get hitsetkey if available
                            if "hitsetkey" in filtered_hits_display.point_data:
                                all_hitsetkeys = filtered_hits_display.point_data[
                                    "hitsetkey"
                                ]
                            else:
                                all_hitsetkeys = [None] * filtered_hits_display.n_points

                            # Create mask for hits that match our criteria
                            mask = np.zeros(len(all_hitkeys), dtype=bool)
                            for i, (hitkey, event, hitsetkey) in enumerate(
                                zip(all_hitkeys, all_events, all_hitsetkeys)
                            ):
                                if (
                                    hitkey in selected_cluster_keys
                                    and event == selected_cluster_event
                                    and (
                                        selected_cluster_hitsetkey is None
                                        or hitsetkey == selected_cluster_hitsetkey
                                    )
                                ):
                                    mask[i] = True

                            # Extract only the matching hits
                            filtered_hits_display = (
                                filtered_hits_display.extract_points(np.where(mask)[0])
                            )

                            # Skip if no matching hits
                            if filtered_hits_display.n_points == 0:
                                continue

                            hit_color = "green"
                        except Exception as e:
                            print(f"Error filtering associated hits: {e}")
                            hit_color = "blue"

                    else:
                        hit_color = "blue"

                    self.display_filtered_hits_info.append(
                        (filtered_hits_display, hit_color, file_info)
                    )

                    self.plotter_widget.add_mesh(
                        filtered_hits_display,
                        style="points",
                        point_size=5,
                        color=hit_color,
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

        if saved_avg_position_actor is not None:
            self.plotter_widget.renderer.add_actor(saved_avg_position_actor)
            self.avg_position_actor = saved_avg_position_actor
            self.plotter_widget.render()

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
                    )

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

            if self.pick_mode == "pos":
                # self.reconstruct_position()
                return

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

        if hasattr(self, "pad_adc_window"):
            self.pad_adc_window.clear_data()
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
        if hasattr(self, "avg_position_actor") and self.avg_position_actor is not None:
            self.plotter_widget.remove_actor(self.avg_position_actor)
            self.avg_position_actor = None

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

    def show_pad_adc_window(self):
        """Show the pad ADC distribution window."""
        if not self.pad_adc_window.isVisible():
            self.pad_adc_window.show()
        self.pad_adc_window.raise_()
        self.pad_adc_window.activateWindow()

        # If there are already selected points, update the distribution
        if self.selected_points_first:
            self.pad_adc_window.update_distribution(self.selected_points_first)

    def toggle_field_mode(self):
        """Toggle between helix and straight line fitting modes."""
        self.using_line_fitting = not self.using_line_fitting

        if (
            hasattr(self, "toggle_field_action")
            and self.toggle_field_action is not None
        ):
            if self.using_line_fitting:
                self.toggle_field_action.setText("Field Off")
                self.instruction_label.setText(
                    "Instruction: Select 2 points for initial line fitting."
                )
            else:
                self.toggle_field_action.setText("Field On")
                self.instruction_label.setText(
                    "Instruction: Select 3 points for initial helix fitting."
                )
        else:
            print("Warning: toggle_field_action not found")

        # Clear any existing selections
        self.selected_points_first.clear()
        self.update_display()

    def toggle_associated_hits(self):
        if self.selected_cluster_index is None or self.last_picked_mesh is None:
            QMessageBox.information(
                self,
                "Toggle Associated Hits",
                "No cluster selected. Please click on a cluster first.",
            )
            return

        print("\nDEBUG DATA ACCESS:")
        print(f"Selected cluster index: {self.selected_cluster_index}")

        """for filename, file_info in self.loaded_files.items():
            cluster_data = file_info["cluster"]
            if cluster_data is not None:
                cluskey = cluster_data.point_data["cluskey"][
                    self.selected_cluster_index
                ]
                print(f"Accessing via loaded_files[{filename}]: {cluskey}")"""

        """cluster_cluskey = self.last_picked_mesh.point_data["cluskey"][
            self.selected_cluster_index
        ]
        selected_cluster_keys = self.last_picked_mesh.hitkeys[
            self.selected_cluster_index
        ]
        cluster_event = self.last_picked_mesh.point_data["event"][
            self.selected_cluster_index
        ]

        print(f"Selected Cluster cluskey: {cluster_cluskey}")
        print(f"Associated hitkeys: {selected_cluster_keys}")
        print(f"Number of associated hitkeys: {len(selected_cluster_keys)}")
        found_match = False
        original_cluskey = None"""

        # Toggle the flag
        self.associated_hits_hidden = not getattr(self, "associated_hits_hidden", False)

        # Update button text
        button_text = (
            "Show All Hits" if self.associated_hits_hidden else "Show Associated Hits"
        )
        self.btn_toggle_associated_hits.setText(button_text)
        if self.associated_hits_hidden:
            # Store current checkbox states for later restoration
            self.previous_hits_state = {
                "tpc_hits": self.hits_checkbox_tpc.isChecked(),
                "silicon_hits": self.hits_checkbox_silicon.isChecked(),
            }

            # Force-enable hits temporarily (without changing checkbox UI)
            self.hits_override_enabled = True

            # Collect associated hits for calculation
            # self.collect_associated_hits_and_calculate_position()
        else:
            # Restore to previous state - no need to use overrides
            self.hits_override_enabled = False

            # Remove the average position marker if it exists
            if (
                hasattr(self, "avg_position_actor")
                and self.avg_position_actor is not None
            ):
                self.plotter_widget.remove_actor(self.avg_position_actor)
                self.avg_position_actor = None

        """if self.associated_hits_hidden:
            try:
                # Get cluster information
                cluster_event = self.last_picked_mesh.point_data["event"][
                    self.selected_cluster_index
                ]
                # Get hitsetkey if available
                if "hitsetkey" in self.last_picked_mesh.point_data:
                    cluster_hitsetkey = self.last_picked_mesh.point_data["hitsetkey"][
                        self.selected_cluster_index
                    ]
                else:
                    cluster_hitsetkey = None
                cluster_cluskey = self.last_picked_mesh.point_data["cluskey"][
                    self.selected_cluster_index
                ]
                selected_cluster_keys = self.last_picked_mesh.hitkeys[
                    self.selected_cluster_index
                ]

                print(f"Selected Cluster cluskey: {cluster_cluskey}")
                print(f"Selected Cluster event: {cluster_event}")
                print(f"Selected Cluster hitsetkey: {cluster_hitsetkey}")
                print(f"Associated hitkeys: {selected_cluster_keys}")
                print(f"Number of associated hitkeys: {len(selected_cluster_keys)}")

                # Collect hit information
                hit_info_list = []
                t_sum = 0
                adc_sum = 0

                for filename, file_info in self.loaded_files.items():
                    hit_data = file_info["hit"]
                    if hit_data is None:
                        continue

                    all_hitkeys = hit_data.point_data.get("hitkeykey", None)
                    all_events = hit_data.point_data.get("event", None)
                    # Get hitsetkey if available
                    if "hitsetkey" in hit_data.point_data:
                        all_hitsetkeys = hit_data.point_data["hitsetkey"]
                    else:
                        all_hitsetkeys = [None] * hit_data.n_points

                    if all_hitkeys is None or all_events is None:
                        continue

                    # Find matching hits with all three identifiers
                    hit_indices = []
                    for i, (key, event, hitsetkey) in enumerate(
                        zip(all_hitkeys, all_events, all_hitsetkeys)
                    ):
                        if (
                            key in selected_cluster_keys
                            and event == cluster_event
                            and (
                                cluster_hitsetkey is None
                                or hitsetkey == cluster_hitsetkey
                            )
                        ):
                            hit_indices.append(i)

                    if not hit_indices:
                        continue

                    print(f"Found {len(hit_indices)} matching hits in file {filename}")

                    # Collect information for these hits
                    for idx in hit_indices:
                        t = hit_data.point_data["t"][idx]
                        adc = hit_data.point_data["adc"][idx]
                        tdriftmax = hit_data.point_data.get(
                            "tdriftmax", [0] * hit_data.n_points
                        )[idx]
                        drift_velocity = hit_data.point_data.get(
                            "driftVelocity", [0] * hit_data.n_points
                        )[idx]
                        hitkey = all_hitkeys[idx]

                        hit_info = {
                            "index": idx,
                            "hitkey": hitkey,
                            "t": t,
                            "adc": adc,
                            "tdriftmax": tdriftmax,
                            "drift_velocity": drift_velocity,
                        }
                        hit_info_list.append(hit_info)

                        # Accumulate for statistical calculations
                        t_sum += t * adc
                        adc_sum += adc

                # Display hit information
                if hit_info_list:
                    info_text = f"Associated Hits for Cluster (cluskey: {cluster_cluskey}, event: {cluster_event}):\n"
                    info_text += f"Total associated hits: {len(hit_info_list)}\n"
                    info_text += "-" * 50 + "\n"

                    # Show a sample of hits (to avoid overwhelming the display)
                    max_hits_to_show = min(10, len(hit_info_list))
                    for i in range(max_hits_to_show):
                        hit = hit_info_list[i]
                        info_text += f"Hit {hit['index']} (hitkey: {hit['hitkey']}):\n"
                        info_text += f"  t: {hit['t']:.2f}\n"
                        info_text += f"  adc: {hit['adc']:.2f}\n"
                        info_text += f"  tdriftmax: {hit['tdriftmax']:.2f}\n"
                        info_text += f"  driftVelocity: {hit['drift_velocity']:.2f}\n"
                        info_text += "-" * 30 + "\n"

                    if len(hit_info_list) > max_hits_to_show:
                        info_text += f"... and {len(hit_info_list) - max_hits_to_show} more hits\n"

                    # Calculate cluster z position if possible
                    if adc_sum > 0:
                        clust = t_sum / adc_sum
                        # Use values from the first hit for tdriftmax and drift_velocity
                        first_hit = hit_info_list[0]
                        tdriftmax = first_hit["tdriftmax"]
                        drift_velocity = first_hit["drift_velocity"]

                        zdriftlength = clust * drift_velocity
                        clustz = tdriftmax * drift_velocity - zdriftlength

                        info_text += "\nCluster Calculations:\n"
                        info_text += f"ADC-weighted average t: {clust:.2f}\n"
                        info_text += f"Calculated clustz: {clustz:.2f} cm\n"

                    self.info_panel.setText(info_text)
                else:
                    self.info_panel.setText(
                        "No associated hits found for the selected cluster."
                    )
            except Exception as e:
                print(f"Error displaying hit information: {e}")
                self.info_panel.setText(f"Error displaying hit information: {e}")"""
        """for filename, file_info in self.loaded_files.items():
            cluster_data = file_info["cluster"]

            if cluster_data is not None:
                try:
                    original_cluskey = cluster_data.point_data["cluskey"][
                        self.selected_cluster_index
                    ]
                    current_cluskey = self.last_picked_mesh.point_data["cluskey"][
                        self.selected_cluster_index
                    ]
                    if original_cluskey == current_cluskey:
                        found_match = True
                        print(
                            f"Verified cluster key matches original data in file: {filename}"
                        )
                        break
                except (IndexError, KeyError):
                    continue

        if not found_match:
            print("Warning: Could not verify cluster key against original data")"""

        """for filename, file_info in self.loaded_files.items():
            hit_data = file_info["hit"]
            if hit_data is not None:
                all_hitkeys = hit_data.point_data.get("hitkeykey", None)
                all_events = hit_data.point_data.get("event", None)
                if all_hitkeys is not None and all_events is not None:
                    # Get indices of associated hits
                    hit_indices = [
                        i
                        for i, (key, event) in enumerate(zip(all_hitkeys, all_events))
                        if key in selected_cluster_keys and event == cluster_event
                    ]

                    if hit_indices:
                        # Initialize sums for clustz calculation

                        print(
                            f"Found {len(hit_indices)} matching hits in file {filename}"
                        )
                        t_sum = 0
                        adc_sum = 0

                        # Display information for each associated hit
                        info_text = f"Associated Hits Information for event {cluster_event} (Total hits: {len(hit_indices)}):\n"
                        info_text += f"Cluster key: {cluster_cluskey}\n"
                        info_text += f"Associated hitkeys: {selected_cluster_keys}\n"
                        info_text += "-" * 50 + "\n"

                        for idx in hit_indices:
                            t = hit_data.point_data["t"][idx]
                            adc = hit_data.point_data["adc"][idx]
                            tdriftmax = hit_data.point_data["tdriftmax"][idx]
                            # zdriftlength = hit_data.point_data["zdriftlength"][idx]
                            drift_velocity = hit_data.point_data["driftVelocity"][idx]
                            hitkey = all_hitkeys[idx]
                            event = all_events[idx]

                            info_text += (
                                f"Hit {idx} (hitkey: {hitkey}, event: {event}):\n"
                            )
                            info_text += f"  t: {t:.2f}\n"
                            info_text += f"  adc: {adc:.2f}\n"
                            info_text += f"  tdriftmax: {tdriftmax:.2f}\n"
                            # info_text += f"  zdriftlength: {zdriftlength:.2f}\n"
                            info_text += f"  driftVelocity: {drift_velocity:.2f}\n"
                            info_text += "-" * 30 + "\n"

                            # Accumulate for clustz calculation
                            t_sum += t * adc
                            adc_sum += adc

                        if adc_sum > 0:
                            clust = t_sum / adc_sum
                            # Use values from the first hit for tdriftmax and drift_velocity
                            first_hit = hit_indices[0]
                            tdriftmax = hit_data.point_data["tdriftmax"][first_hit]
                            drift_velocity = hit_data.point_data["driftVelocity"][
                                first_hit
                            ]

                            zdriftlength = clust * drift_velocity
                            clustz = tdriftmax * drift_velocity - zdriftlength

                            info_text += "\nCluster Calculations:\n"
                            info_text += f"Calculated clustz: {clustz:.2f}\n"

                        self.info_panel.setText(info_text)
"""

        """try:
            # Retrieve the hit keys for the selected cluster.
            selected_cluster_keys = self.cluster_polydata.hitkeys[
                self.selected_cluster_index
            ]
        except Exception as e:
            QMessageBox.warning(
                self,
                "Toggle Associated Hits",
                f"Error retrieving hit keys for selected cluster: {e}",
            )
            return

        try:
            # Assuming the 'cluskey' branch is stored in point_data and is an array.
            cluster_cluskey = self.cluster_polydata.point_data["cluskey"][
                self.selected_cluster_index
            ]
        except Exception as e:
            QMessageBox.warning(
                self,
                "Toggle Associated Hits",
                f"Error retrieving cluskey for selected cluster: {e}",
            )
            return"""

        # print(f"Selected Cluster cluskey: direct {cluster_cluskey}")
        # print(f"Associated hitkeys: {selected_cluster_keys}")

        """if hasattr(self, "last_picked_mesh"):
            mesh_cluskey = self.last_picked_mesh.point_data["cluskey"][
                self.selected_cluster_index
            ]
            print(f"Via last picked mesh: {mesh_cluskey}")"""

        """try:
            # Retrieve all hit keys from the hit polydata.
            all_hitkeys = self.hit_polydata.point_data["hitkeykey"]
        except Exception as e:
            QMessageBox.warning(
                self, "Toggle Associated Hits", f"Hit keys not found in hit data: {e}"
            )
            return

        # Filter the hits: create a mask that is True for hits NOT in the selected cluster.
        mask = ~np.isin(all_hitkeys, selected_cluster_keys)
        associated_indices = np.where(mask)[0]

        if associated_indices.size == 0:
            QMessageBox.information(
                self,
                "Toggle Associated Hits",
                "No associated hits found for the selected cluster.",
            )
            return

        # For a robust solution, use a flag (e.g. self.associated_hits_hidden) and refresh update_display.
        self.associated_hits_hidden = not getattr(self, "associated_hits_hidden", False)
        # Update the toggle button text to reflect the state.
        new_text = (
            "Show Associated Hits"
            if self.associated_hits_hidden
            else "Hide Associated Hits"
        )
        self.btn_toggle_associated_hits.setText(new_text)
        # Call update_display so that when the hits are added, they are filtered accordingly.
        self.update_display()"""

        """new_text = (
            "Show Associated Hits"
            if self.associated_hits_hidden
            else "Hide Associated Hits"
        )"""

        self.update_display()
        if self.associated_hits_hidden:
            self.collect_associated_hits_and_calculate_position()

    '''def reconstruct_position(self):
        """Calculate cluster z-position from selected hits."""
        if not self.selected_points_first:
            QMessageBox.warning(
                self,
                "Position Reconstruction",
                "Please select hits for position reconstruction first.",
            )
            return

        t_sum = 0
        adc_sum = 0

        # Use the stored hit information directly
        for hit_info in self.selected_points_first:
            t_sum += hit_info["t"] * hit_info["adc"]
            adc_sum += hit_info["adc"]

        if adc_sum > 0:
            # Use the first hit's drift parameters
            first_hit = self.selected_points_first[0]
            clust = t_sum / adc_sum
            drift_velocity = first_hit["drift_velocity"]
            tdriftmax = first_hit["tdriftmax"]

            zdriftlength = clust * drift_velocity
            clustz = tdriftmax * drift_velocity - zdriftlength

            info_text = "Position Reconstruction Results:\n"
            info_text += f"Number of hits used: {len(self.selected_points_first)}\n"
            info_text += f"Calculated cluster z: {clustz:.2f} cm\n"
            info_text += f"Average time: {clust:.2f}\n"
            info_text += f"Drift velocity: {drift_velocity:.2f}\n"
            self.info_panel.setText(info_text)
        else:
            QMessageBox.warning(
                self,
                "Position Reconstruction",
                "Could not calculate position. Please ensure valid hits are selected.",
            )

        # Clear selections after reconstruction
        self.selected_points_first.clear()
        self.update_display()'''

    def output_adc_averaged_phi_z(self):
        """Compute and output the ADC-weighted average phi and z position from selected hits."""
        if not self.selected_points_first:
            QMessageBox.warning(
                self,
                "Output Error",
                "Please select hits for position reconstruction first.",
            )
            return

        phi_sum = 0.0
        z_sum = 0.0
        adc_sum = 0.0

        # Calculate weighted sums
        for hit in self.selected_points_first:
            adc = hit.get("adc", 0)
            phi = hit.get("phi", 0)
            z = hit.get("z", 0)  # Get direct z coordinate

            phi_sum += phi * adc
            z_sum += z * adc
            adc_sum += adc

        if adc_sum > 0:
            avg_phi = phi_sum / adc_sum
            avg_z = z_sum / adc_sum

            info_text = "ADC-Weighted Position:\n"
            info_text += f"Phi: {avg_phi:.3f} rad ({np.degrees(avg_phi):.2f}°)\n"
            info_text += f"Z: {avg_z:.2f} cm"
            self.info_panel.setText(info_text)
        else:
            QMessageBox.warning(
                self, "Output Error", "Total ADC is zero; cannot compute averages."
            )

    def collect_associated_hits_and_calculate_position(self):
        """
        Collect information about associated hits, calculate and visualize the ADC-weighted average position.
        """
        try:

            if self.last_picked_mesh is None:
                self.info_panel.setText("Error: No cluster has been selected yet.")
                return None

            if self.selected_cluster_index is None:
                self.info_panel.setText("Error: No cluster index is selected.")
                return None

            if self.selected_cluster_index >= len(self.last_picked_mesh.hitkeys):
                self.info_panel.setText(
                    f"Error: Cluster index {self.selected_cluster_index} is out of bounds."
                )
                return None
            # Get cluster information
            cluster_event = self.last_picked_mesh.point_data["event"][
                self.selected_cluster_index
            ]
            selected_cluster_keys = self.last_picked_mesh.hitkeys[
                self.selected_cluster_index
            ]

            # Get hitsetkey if available
            if "hitsetkey" in self.last_picked_mesh.point_data:
                cluster_hitsetkey = self.last_picked_mesh.point_data["hitsetkey"][
                    self.selected_cluster_index
                ]
            else:
                cluster_hitsetkey = None

            # Variables for ADC-weighted calculations
            x_sum = 0.0
            y_sum = 0.0
            z_sum = 0.0
            phi_sum = 0.0
            t_sum = 0.0
            adc_sum = 0.0
            associated_hit_count = 0

            # Collect hit information from all loaded files
            for filename, file_info in self.loaded_files.items():
                hit_data = file_info["hit"]
                if hit_data is None:
                    continue

                all_hitkeys = hit_data.point_data.get("hitkeykey", None)
                all_events = hit_data.point_data.get("event", None)

                # Get hitsetkey if available
                if "hitsetkey" in hit_data.point_data:
                    all_hitsetkeys = hit_data.point_data["hitsetkey"]
                else:
                    all_hitsetkeys = [None] * hit_data.n_points

                if all_hitkeys is None or all_events is None:
                    continue

                # Find matching hits
                associated_hit_indices = []
                for i, (key, event, hitsetkey) in enumerate(
                    zip(all_hitkeys, all_events, all_hitsetkeys)
                ):
                    if (
                        key in selected_cluster_keys
                        and event == cluster_event
                        and (
                            cluster_hitsetkey is None or hitsetkey == cluster_hitsetkey
                        )
                    ):
                        associated_hit_indices.append(i)

                # Process the associated hits
                for idx in associated_hit_indices:
                    # Get hit position
                    x = hit_data.points[idx][0]
                    y = hit_data.points[idx][1]
                    z = hit_data.points[idx][2]

                    # Get hit properties
                    adc = hit_data.point_data["adc"][idx]
                    t = hit_data.point_data.get("t", [0] * hit_data.n_points)[idx]
                    phi = hit_data.point_data.get("phi", [0] * hit_data.n_points)[idx]

                    # Accumulate weighted values
                    x_sum += x * adc
                    y_sum += y * adc
                    z_sum += z * adc
                    phi_sum += phi * adc if phi != 0 else 0
                    t_sum += t * adc if t != 0 else 0
                    adc_sum += adc
                    associated_hit_count += 1

            # Calculate ADC-weighted averages if we found associated hits
            if associated_hit_count > 0 and adc_sum > 0:
                avg_x = x_sum / adc_sum
                avg_y = y_sum / adc_sum
                avg_z = z_sum / adc_sum
                avg_phi = phi_sum / adc_sum if phi_sum != 0 else 0
                avg_t = t_sum / adc_sum if t_sum != 0 else 0

                # Display the results in the info panel
                info_text = f"Associated Hits for Selected Cluster:\n"
                info_text += f"Total hits: {associated_hit_count}\n"
                info_text += f"Total ADC: {adc_sum:.1f}\n\n"
                info_text += f"ADC-Weighted Average Position:\n"
                info_text += f"X: {avg_x:.2f} cm\n"
                info_text += f"Y: {avg_y:.2f} cm\n"
                info_text += f"Z: {avg_z:.2f} cm\n"
                info_text += f"Phi: {avg_phi:.3f} rad ({np.degrees(avg_phi):.2f}°)\n"
                info_text += f"T: {avg_t:.2f}\n"

                # Display DCA information if available
                if hasattr(self, "helix_params_refined") and self.helix_params_refined:
                    dca_result = find_dca_and_closest_point(
                        self.helix_params_refined,
                        self.helix_params_refined.get("ref_theta_direct", 0),
                    )
                    if dca_result is not None:
                        dca_val = dca_result["dca"]
                        info_text += f"\nDCA to Beam: {dca_val:.3f} cm\n"

                self.info_panel.setText(info_text)

                # Create or update the average position point visualization
                avg_position = np.array([avg_x, avg_y, avg_z])
                self.visualize_average_position(avg_position)

                return {
                    "x": avg_x,
                    "y": avg_y,
                    "z": avg_z,
                    "phi": avg_phi,
                    "t": avg_t,
                    "hit_count": associated_hit_count,
                    "adc_sum": adc_sum,
                }
            else:
                self.info_panel.setText(
                    "No associated hits found for the selected cluster."
                )
                return None

        except Exception as e:
            print(f"Error calculating average position: {e}")
            import traceback

            traceback.print_exc()
            self.info_panel.setText(f"Error calculating average position: {str(e)}")
            return None

    def visualize_average_position(self, position):
        """
        Create a visual marker for the ADC-weighted average position.

        Args:
            position: Numpy array [x, y, z]
        """
        try:
            print(f"Visualizing average position at: {position}")

            # Create a single point at the average position (instead of a sphere)
            # This will make it look like a regular cluster point
            point = pv.PolyData(position.reshape(1, 3))

            # Remove existing actor if it exists
            if (
                hasattr(self, "avg_position_actor")
                and self.avg_position_actor is not None
            ):
                print("Removing existing average position actor")
                self.plotter_widget.remove_actor(self.avg_position_actor)

            # Add the point with the same style as clusters but in black color
            print("Adding new average position actor")
            self.avg_position_actor = self.plotter_widget.add_mesh(
                point,
                style="points",  # Same style as clusters
                point_size=8,  # Slightly larger than regular clusters (which are 5)
                color="black",  # Black color as requested
                render=True,  # Force immediate rendering
                pickable=False,
            )

            # Explicit render call to ensure everything is displayed
            self.plotter_widget.render()

            print("Visualization completed")
            return True
        except Exception as e:
            print(f"Error in visualize_average_position: {e}")
            import traceback

            traceback.print_exc()
            return False
