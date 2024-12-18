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
    QGroupBox,  # Added for grouping side checkboxes
)
from qtpy.QtCore import Qt
from superqt import QRangeSlider


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Display")
        self.resize(1200, 800)

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

        # Controls: Button to Load Data
        btn_load = QPushButton("Load ROOT File")
        btn_load.clicked.connect(self.load_data)
        control_layout.addWidget(btn_load)

        # Controls: Checkboxes for showing clusters and hits
        self.show_clusters = QCheckBox("Show Clusters")
        self.show_clusters.setChecked(True)
        self.show_clusters.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_clusters)

        self.show_hits = QCheckBox("Show Hits")
        self.show_hits.setChecked(True)
        self.show_hits.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_hits)

        # === Side Selection Area ===
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

        # Now let us add range slider to adjust the range
        # Range sliders for X, Y, Z axes
        # range_layout = QFormLayout()
        slider_width = 300
        # self.x_range_slider = QRangeSlider()
        self.x_range_slider = QRangeSlider(Qt.Horizontal)
        self.x_range_slider.setMinimum(-100)
        self.x_range_slider.setMaximum(100)
        self.x_range_slider.setValue([-100, 100])  # Initial range
        # self.x_range_slider.setBarMovesAllHandles(False)
        # self.x_range_slider.setOrientation(Qt.Horizontal)

        self.x_range_slider.setFixedWidth(slider_width)

        self.x_range_slider.sliderReleased.connect(self.update_display)

        self.y_range_slider = QRangeSlider(Qt.Horizontal)
        self.y_range_slider.setMinimum(-100)
        self.y_range_slider.setMaximum(100)
        self.y_range_slider.setValue([-100, 100])
        # self.y_range_slider.setOrientation(Qt.Horizontal)
        self.y_range_slider.setFixedWidth(slider_width)
        self.y_range_slider.sliderReleased.connect(self.update_display)

        self.z_range_slider = QRangeSlider(Qt.Horizontal)
        self.z_range_slider.setMinimum(-400)
        self.z_range_slider.setMaximum(400)
        self.z_range_slider.setValue([-400, 400])
        # self.z_range_slider.setOrientation(Qt.Horizontal)
        self.z_range_slider.setFixedWidth(slider_width)
        self.z_range_slider.sliderReleased.connect(self.update_display)

        # Labels to show the range dynamically
        self.x_label = QLabel("X Range: [-100, 100]")
        self.y_label = QLabel("Y Range: [-100, 100]")
        self.z_label = QLabel("Z Range: [-400, 400]")

        # Set fixed widths for labels so that their size does not change with text length
        label_width = 150
        self.x_label.setFixedWidth(label_width)
        self.y_label.setFixedWidth(label_width)
        self.z_label.setFixedWidth(label_width)

        # Connect range slider signals
        self.x_range_slider.valueChanged.connect(lambda v: self.update_label("X", v))
        self.y_range_slider.valueChanged.connect(lambda v: self.update_label("Y", v))
        self.z_range_slider.valueChanged.connect(lambda v: self.update_label("Z", v))

        range_layout = QFormLayout()

        # Add sliders and labels to form layout
        range_layout.addRow(self.x_label, self.x_range_slider)
        range_layout.addRow(self.y_label, self.y_range_slider)
        range_layout.addRow(self.z_label, self.z_range_slider)
        control_layout.addLayout(range_layout)

        # controls_layout.addWidget()

        # Additional filters if needed (passed, used...)
        # For now, just remove or comment out if not relevant:
        # self.show_passed = QCheckBox("Show Passed Straight")
        # self.show_passed.stateChanged.connect(self.update_display)
        # controls_layout.addWidget(self.show_passed)

        # self.show_used = QCheckBox("Show Used in Seed")
        # self.show_used.stateChanged.connect(self.update_display)
        # controls_layout.addWidget(self.show_used)

        left_layout.addWidget(control_area)

        # 3D Visualization Area
        self.plotter_widget = QtInteractor()
        left_layout.addWidget(self.plotter_widget)

        # Add a frame line at the bottom for neatness (optional)
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

        # Show axes
        self.plotter_widget.show_axes()

    def load_data(self):
        # Prompt the user to choose a ROOT file
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root)"
        )
        if not filename:
            return

        try:
            # Load data using uproot
            file = uproot.open(filename)
            # Load cluster tree
            cluster_tree = file["combined_clusters"]
            cluster_data = cluster_tree.arrays(library="np")

            # Ensure 'side' branch exists
            if "side" not in cluster_data:
                raise ValueError("The 'side' branch is missing in the cluster data.")

            # Create cluster polydata
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

        # Retrieve range values
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

        # Add clusters if requested and we have data
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
            # Add hits if requested and we have data
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
        # Disable and re-enable picking to ensure a fresh start
        self.plotter_widget.disable_picking()
        # Still using use_mesh=True for now
        self.plotter_widget.enable_point_picking(
            callback=self.on_point_picked, show_message=True, use_mesh=True
        )
        # Restore the previously saved camera position
        self.plotter_widget.camera_position = camera_position

        self.plotter_widget.show_axes()

    def update_label(self, axis, value):
        label = f"{axis} Range: [{value[0]}, {value[1]}]"
        if axis == "X":
            self.x_label.setText(label)
        elif axis == "Y":
            self.y_label.setText(label)
        elif axis == "Z":
            self.z_label.setText(label)

    def on_point_picked(self, mesh, point_id):

        if point_id < 0:
            return
        # picked_point = mesh.points[point_id]
        data_type = mesh.point_data["data_type"][point_id]
        # picked_data = None
        # pick_idx = None
        # label = ""
        if data_type == 0:
            label = "Cluster"
        else:
            label = "Hit"

        # Display info
        info_text = f"{label} (Point ID: {point_id})\n"
        for attr in mesh.point_data.keys():
            value = mesh.point_data[attr][point_id]
            info_text += f"{attr}: {value}\n"

        self.info_panel.setText(info_text)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
