import sys
import numpy as np
import uproot
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt5.QtWidgets import (
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
)
from PyQt5.QtCore import Qt

# import qt_material
from qtrangeslider import QRangeSlider


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

        # Controls: Checkboxes
        self.show_clusters = QCheckBox("Show Clusters")
        self.show_clusters.setChecked(True)
        self.show_clusters.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_clusters)

        self.show_hits = QCheckBox("Show Hits")
        self.show_hits.setChecked(True)
        self.show_hits.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_hits)

        # Now let us add range slider to adjust the range
        # Range sliders for X, Y, Z axes
        # range_layout = QFormLayout()
        slider_width = 300
        self.x_range_slider = QRangeSlider()
        self.x_range_slider.setMinimum(-100)
        self.x_range_slider.setMaximum(100)
        self.x_range_slider.setValue([-100, 100])  # Initial range
        self.x_range_slider.setBarMovesAllHandles(False)
        self.x_range_slider.setOrientation(Qt.Horizontal)

        self.x_range_slider.setFixedWidth(slider_width)

        self.x_range_slider.sliderReleased.connect(self.update_display)

        self.y_range_slider = QRangeSlider()
        self.y_range_slider.setMinimum(-100)
        self.y_range_slider.setMaximum(100)
        self.y_range_slider.setValue([-100, 100])
        self.y_range_slider.setOrientation(Qt.Horizontal)
        self.y_range_slider.setFixedWidth(slider_width)
        self.y_range_slider.sliderReleased.connect(self.update_display)

        self.z_range_slider = QRangeSlider()
        self.z_range_slider.setMinimum(-200)
        self.z_range_slider.setMaximum(200)
        self.z_range_slider.setValue([-200, 200])
        self.z_range_slider.setOrientation(Qt.Horizontal)
        self.z_range_slider.setFixedWidth(slider_width)
        self.z_range_slider.sliderReleased.connect(self.update_display)

        # Labels to show the range dynamically
        self.x_label = QLabel("X Range: [-100, 100]")
        self.y_label = QLabel("Y Range: [-100, 100]")
        self.z_label = QLabel("Z Range: [-200, 200]")

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

            # Create cluster polydata
            cx = cluster_data["seed_x"]
            cy = cluster_data["seed_y"]
            cz = cluster_data["seed_z"]
            cpoints = np.column_stack([cx, cy, cz])
            cluster_polydata = pv.PolyData(cpoints)
            for name in cluster_data.keys():
                if name not in ("seed_x", "seed_y", "seed_z"):
                    cluster_polydata.point_data[name] = cluster_data[name]

            self.cluster_data = cluster_polydata

            # Load hits tree
            hits_tree = file["combined_hits"]
            hits_data = hits_tree.arrays(library="np")

            hx = hits_data["gx"]
            hy = hits_data["gy"]
            hz = hits_data["gz"]
            hpoints = np.column_stack([hx, hy, hz])
            hit_polydata = pv.PolyData(hpoints)

            # Add other attributes
            for name in hits_data.keys():
                if name not in ("gx", "gy", "gz"):
                    hit_polydata.point_data[name] = hits_data[name]

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
            mask = (
                (points[:, 0] >= x_min)
                & (points[:, 0] <= x_max)  # X range
                & (points[:, 1] >= y_min)
                & (points[:, 1] <= y_max)  # Y range
                & (points[:, 2] >= z_min)
                & (points[:, 2] <= z_max)  # Z range
            )
            # Extract filtered points
            filtered_points = self.cluster_data.extract_points(np.where(mask)[0])

            if (
                filtered_points.n_points > 0
            ):  # Check if any points are left after filtering
                self.cluster_polydata = filtered_points
                # Add cluster points in one color, e.g., red
                self.plotter_widget.add_points(
                    self.cluster_polydata,
                    render_points_as_spheres=True,
                    point_size=5,
                    color=[1, 0, 0],
                )
            # Add hits if requested and we have data
        if show_hits and self.hit_data is not None and self.hit_data.n_points > 0:

            # Apply the slider range to filter points
            points = self.hit_data.points  # Get hit points
            mask = (
                (points[:, 0] >= x_min)
                & (points[:, 0] <= x_max)  # X range
                & (points[:, 1] >= y_min)
                & (points[:, 1] <= y_max)  # Y range
                & (points[:, 2] >= z_min)
                & (points[:, 2] <= z_max)  # Z range
            )

            # Extract filtered points
            filtered_points = self.hit_data.extract_points(np.where(mask)[0])

            if (
                filtered_points.n_points > 0
            ):  # Check if any points are left after filtering
                self.hit_polydata = filtered_points
                # Add hit points in another color, e.g., blue
                self.plotter_widget.add_points(
                    self.hit_polydata,
                    render_points_as_spheres=True,
                    point_size=5,
                    color=[0, 0, 1],
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
        picked_point = mesh.points[point_id]
        # Determine if the picked mesh corresponds to clusters or hits
        # The 'mesh' returned is a PyVista mesh subset. We need to check which one was picked.
        # A simple approach: compare number of points and point positions or store references.
        # If you have multiple calls to add_points, PyVista returns a composite mesh.
        # In this simplified scenario, you might get the last added mesh. Consider comparing:
        #    - If point_id < self.cluster_polydata.n_points: cluster was picked
        # else it's from hit data if point_id < self.hit_polydata.n_points
        #
        # However, this approach may need refinement depending on how PyVista returns the picked mesh.
        # Another approach: just compare the number of points.
        # If both sets exist, you must determine which set the point belongs to.
        #
        # Let's assume PyVista returns the exact mesh we clicked on. Then:
        picked_data = None
        pick_idx = None
        label = ""
        # Check if mesh is cluster or hit mesh by identity or by comparing with stored polydata
        # PyVista's enable_point_picking creates a new mesh for picked points. Checking identity may be tricky.
        # Instead, check coordinates at point_id and see if they match cluster_polydata or hit_polydata.

        # Extract picked point coordinates
        picked_point = mesh.points[point_id]

        # Check if this point exists in cluster_polydata
        if self.cluster_polydata is not None:
            cluster_coords = self.cluster_polydata.points
            # Find matching point
            c_indices = np.where((cluster_coords == picked_point).all(axis=1))[0]
            if len(c_indices) > 0:
                picked_data = self.cluster_polydata
                pick_idx = c_indices[0]
                label = "Cluster"

        # If not found in clusters, check hits
        if picked_data is None and self.hit_polydata is not None:
            hit_coords = self.hit_polydata.points
            h_indices = np.where((hit_coords == picked_point).all(axis=1))[0]
            if len(h_indices) > 0:
                picked_data = self.hit_polydata
                pick_idx = h_indices[0]
                label = "Hit"

        if picked_data is None:
            # Not found in either dataset; might be an error or the picking mesh is different
            return

        # Display info
        info_text = f"{label} (Point ID: {pick_idx})\n"
        for attr in picked_data.point_data.keys():
            value = picked_data.point_data[attr][pick_idx]
            info_text += f"{attr}: {value}\n"

        self.info_panel.setText(info_text)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
