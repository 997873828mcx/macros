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
)
from PyQt5.QtCore import Qt


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

        # Buttons and Checkboxes at the top
        controls_layout = QHBoxLayout()
        btn_load = QPushButton("Load ROOT File")
        btn_load.clicked.connect(self.load_data)
        controls_layout.addWidget(btn_load)

        self.show_passed = QCheckBox("Show Passed Straight")
        self.show_passed.stateChanged.connect(self.update_display)
        controls_layout.addWidget(self.show_passed)

        self.show_used = QCheckBox("Show Used in Seed")
        self.show_used.stateChanged.connect(self.update_display)
        controls_layout.addWidget(self.show_used)

        left_layout.addLayout(controls_layout)

        # 3D Visualization Area
        self.plotter_widget = QtInteractor()
        left_layout.addWidget(self.plotter_widget)

        # Add a frame line at the bottom for neatness (optional)
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        left_layout.addWidget(line)

        # Sidebar for cluster information
        sidebar = QWidget()
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(5, 5, 5, 5)

        info_label = QLabel("Cluster Information")
        sidebar_layout.addWidget(info_label)

        self.info_panel = QTextEdit()
        self.info_panel.setReadOnly(True)
        sidebar_layout.addWidget(self.info_panel)

        splitter.addWidget(left_panel)
        splitter.addWidget(sidebar)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 1)

        # Data placeholders
        self.polydata = None
        self.full_data = None

        # Show axes
        self.plotter_widget.show_axes()

    def load_data(self):
        # Prompt the user to choose a ROOT file
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root)"
        )
        if not filename:
            return  # User canceled

        try:
            # Load data using uproot
            file = uproot.open(filename)
            tree = file["tracking_clusters"]  # Adjust if your tree name differs
            data = tree.arrays(library="np")

            # Extract coordinates
            x = data["x"]
            y = data["y"]
            z = data["z"]
            points = np.column_stack([x, y, z])

            # Create a PyVista PolyData
            polydata = pv.PolyData(points)

            # Add other attributes
            for name in data.keys():
                if name not in ("x", "y", "z"):
                    polydata.point_data[name] = data[name]

            self.full_data = polydata
            self.update_display()

        except Exception as e:
            self.info_panel.setText(f"Error loading file:\n{e}")

    def update_display(self):
        # Clear the current plotter
        self.plotter_widget.clear()

        # If no data loaded yet, return
        if self.full_data is None:
            return

        # Filter data based on checkboxes
        mask = np.ones(self.full_data.n_points, dtype=bool)
        if (
            "passed_straight" in self.full_data.point_data
            and self.show_passed.isChecked()
        ):
            mask &= self.full_data.point_data["passed_straight"] == 1
        if "used_in_seed" in self.full_data.point_data and self.show_used.isChecked():
            mask &= self.full_data.point_data["used_in_seed"] == 1

        # Create a subset
        subset = self.full_data.extract_points(np.where(mask)[0])

        if subset.n_points > 0:
            # Add points to the plotter_widget directly
            self.plotter_widget.add_points(
                subset, render_points_as_spheres=True, point_size=5
            )
            self.plotter_widget.reset_camera()

            # Enable point picking on plotter_widget directly

            self.plotter_widget.disable_picking()
            self.plotter_widget.enable_point_picking(
                callback=self.on_point_picked, show_message=True, use_picker=True
            )
            self.polydata = subset
        else:
            self.info_panel.setText("No points to display with current filters.")
            self.polydata = None

        self.plotter_widget.show_axes()

    def on_point_picked(self, picker, event):

        if point_id < 0 or self.polydata is None:
            return

        info_text = f"Point ID: {point_id}\n"
        for attr in self.polydata.point_data.keys():
            value = self.polydata.point_data[attr][point_id]
            info_text += f"{attr}: {value}\n"

        self.info_panel.setText(info_text)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
