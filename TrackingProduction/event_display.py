import sys
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QWidget,
    QPushButton,
    QFileDialog,
    QHBoxLayout,
    QCheckBox,
    QTextEdit,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QCursor
import pyqtgraph.opengl as gl
import uproot
import numpy as np
import pandas as pd
import pyqtgraph as pg
from OpenGL.GL import (
    glGetDoublev,
    glGetIntegerv,
    GL_PROJECTION_MATRIX,
    GL_MODELVIEW_MATRIX,
    GL_VIEWPORT,
)
from OpenGL.GLU import gluProject


# Subclass GLViewWidget to handle mouse events properly
class MyGLViewWidget(gl.GLViewWidget):
    def __init__(self, parent=None):
        super().__init__()
        self.parent_widget = parent

    def mousePressEvent(self, event):
        # Call the original method to maintain internal state
        super().mousePressEvent(event)
        # Call the custom method in the parent widget
        if self.parent_widget is not None:
            self.parent_widget.on_mouse_click(event)


class EventDisplay(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Display")
        self.setGeometry(100, 100, 1200, 800)

        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        # Left side: 3D view and controls
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)

        # 3D View widget using the subclassed GLViewWidget
        self.view = MyGLViewWidget(self)  # Use the subclassed widget
        self.view.setBackgroundColor("w")
        self.view.setFocusPolicy(Qt.StrongFocus)
        left_layout.addWidget(self.view)

        # Initialize lists to keep track of items
        self.static_items = []
        self.dynamic_items = []

        # Add coordinate axes and disks
        self.add_coordinate_system()

        # Control panel
        control_layout = QHBoxLayout()

        # Load button
        btn_load = QPushButton("Load ROOT File")
        btn_load.clicked.connect(self.load_data)
        control_layout.addWidget(btn_load)

        # Checkboxes for filtering
        self.show_passed = QCheckBox("Show Passed Straight")
        self.show_passed.setChecked(False)
        self.show_passed.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_passed)

        self.show_used = QCheckBox("Show Used in Final Seed")
        self.show_used.setChecked(False)
        self.show_used.stateChanged.connect(self.update_display)
        control_layout.addWidget(self.show_used)

        left_layout.addLayout(control_layout)

        main_layout.addWidget(left_panel, stretch=7)

        # Right side: Info panel
        self.info_panel = QTextEdit()
        self.info_panel.setReadOnly(True)
        self.info_panel.setMinimumWidth(200)
        self.info_panel.setStyleSheet(
            """
            QTextEdit {
                background-color: white;
                font-family: monospace;
                padding: 10px;
                border: 1px solid gray;
            }
        """
        )
        main_layout.addWidget(self.info_panel, stretch=3)

        # Initialize data storage
        self.data = None
        self.scatter = None
        self.pos_data = None  # Store position data
        self.full_data = None  # Store all branch data

    def add_coordinate_system(self):
        # Add X, Y, Z axes
        axis_length = 500
        axis_width = 5
        x_axis = gl.GLLinePlotItem(
            pos=np.array([[0, 0, 0], [axis_length, 0, 0]]),
            color=(1, 0, 0, 1),
            width=axis_width,
            glOptions="opaque",
        )  # red
        y_axis = gl.GLLinePlotItem(
            pos=np.array([[0, 0, 0], [0, axis_length, 0]]),
            color=(0, 0.5, 0, 1),
            width=axis_width,
            glOptions="opaque",
        )  # green
        z_axis = gl.GLLinePlotItem(
            pos=np.array([[0, 0, 0], [0, 0, axis_length]]),
            color=(0, 0, 1, 1),
            width=axis_width,
            glOptions="opaque",
        )  # blue

        # Add to view and static_items list
        self.view.addItem(x_axis)
        self.view.addItem(y_axis)
        self.view.addItem(z_axis)
        self.static_items.extend([x_axis, y_axis, z_axis])

        # Add TPC disks
        # Create a circular mesh for the disks
        radius = 80  # Adjust based on TPC size
        points = 100  # Number of points to make the circle smooth
        theta = np.linspace(0, 2 * np.pi, points)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        # Create meshes for two disks at different z positions
        z1 = -105  # First disk z position
        z2 = 105  # Second disk z position

        # Create vertices and faces for the disks
        verts1 = np.vstack([x, y, np.ones_like(x) * z1]).T
        verts2 = np.vstack([x, y, np.ones_like(x) * z2]).T

        # Add center point for triangulation
        verts1 = np.vstack([verts1, [0, 0, z1]])
        verts2 = np.vstack([verts2, [0, 0, z2]])

        # Create faces (triangles) for the disks
        faces = []
        for i in range(points - 1):
            faces.append([i, i + 1, points])
        faces.append([points - 1, 0, points])  # Close the circle

        # Convert faces to a NumPy array
        faces = np.array(faces, dtype=np.int32)

        # Create mesh items with semi-transparency
        disk1 = gl.GLMeshItem(
            vertexes=verts1,
            faces=faces,
            smooth=True,
            color=(0.5, 0.5, 0.5, 0.3),  # Gray, semi-transparent
            shader="shaded",
        )
        disk2 = gl.GLMeshItem(
            vertexes=verts2,
            faces=faces,
            smooth=True,
            color=(0.5, 0.5, 0.5, 0.3),  # Gray, semi-transparent
            shader="shaded",
        )

        # Add disks to view and static_items list
        self.view.addItem(disk1)
        self.view.addItem(disk2)
        self.static_items.extend([disk1, disk2])

    def load_data(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root)"
        )
        if filename:
            try:
                # Load ROOT file using uproot
                file = uproot.open(filename)
                tree = file["tracking_clusters"]  # Using your tree name
                self.data = tree.arrays(
                    expressions=["x", "y", "z", "passed_straight", "used_in_seed"],
                    library="pd",
                )
                self.update_display()
            except Exception as e:
                self.info_panel.setText(f"Error loading file:\n{e}")

    def update_display(self):
        if self.data is not None:
            # Remove previous dynamic items
            for item in self.dynamic_items:
                self.view.removeItem(item)
            self.dynamic_items = []

            # Filter data based on checkboxes
            display_data = self.data.copy()

            if self.show_passed.isChecked():
                display_data = display_data[display_data["passed_straight"] == 1]
            if self.show_used.isChecked():
                display_data = display_data[display_data["used_in_seed"] == 1]

            if len(display_data) > 0:
                # Create position array
                pos = np.column_stack(
                    (
                        display_data["x"].to_numpy(),
                        display_data["y"].to_numpy(),
                        display_data["z"].to_numpy(),
                    )
                )

                self.pos_data = pos
                self.full_data = display_data

                # Print data range to debug
                print("X range:", np.min(pos[:, 0]), "to", np.max(pos[:, 0]))
                print("Y range:", np.min(pos[:, 1]), "to", np.max(pos[:, 1]))
                print("Z range:", np.min(pos[:, 2]), "to", np.max(pos[:, 2]))

                # Color points based on their properties
                colors = np.zeros((len(pos), 4))
                colors[:] = [0, 0, 0.7, 1.0]  # Default blue

                # Add points to display
                scatter = gl.GLScatterPlotItem(
                    pos=pos, color=colors, size=5, pxMode=True, glOptions="opaque"
                )
                self.view.addItem(scatter)
                self.dynamic_items.append(scatter)
                self.scatter = scatter

                # Auto-scale view to data
                self.view.opts["center"] = pg.Vector(
                    np.mean(pos[:, 0]), np.mean(pos[:, 1]), np.mean(pos[:, 2])
                )

    def on_mouse_click(self, event):
        if self.scatter is not None and self.pos_data is not None:
            # Ensure OpenGL context is current
            self.view.makeCurrent()

            # Get mouse position
            pos = event.pos()
            mouse_x, mouse_y = pos.x(), pos.y()

            # Retrieve OpenGL matrices and viewport
            projection = glGetDoublev(GL_PROJECTION_MATRIX)
            modelview = glGetDoublev(GL_MODELVIEW_MATRIX)
            viewport = glGetIntegerv(GL_VIEWPORT)

            # Invert y-coordinate because OpenGL's origin is at the bottom-left corner
            win_y = viewport[3] - mouse_y

            # Project each 3D point to 2D screen coordinates
            window_coords = np.array(
                [
                    gluProject(
                        point[0], point[1], point[2], modelview, projection, viewport
                    )
                    for point in self.pos_data
                ]
            )

            # Compute distances between the mouse click and projected points
            distances = np.sqrt(
                (window_coords[:, 0] - mouse_x) ** 2
                + (window_coords[:, 1] - win_y) ** 2
            )
            nearest_idx = np.argmin(distances)
            min_distance = distances[nearest_idx]

            # Threshold to determine if a cluster is close enough to be selected
            if min_distance < 10:  # Adjust the threshold as needed
                point_info = self.full_data.iloc[nearest_idx]
                info_text = "Cluster Information:\n\n"
                for col in point_info.index:
                    value = point_info[col]
                    if isinstance(value, float):
                        info_text += f"{col:<15}: {value:.3f}\n"
                    else:
                        info_text += f"{col:<15}: {value}\n"
                self.info_panel.setText(info_text)

                # Highlight the selected cluster
                colors = np.zeros((len(self.pos_data), 4))
                colors[:] = [0, 0, 0.7, 1.0]  # Default blue color
                colors[nearest_idx] = [
                    1,
                    0,
                    0,
                    1.0,
                ]  # Red color for the selected cluster
                self.scatter.setData(color=colors)


def main():
    app = QApplication(sys.argv)
    window = EventDisplay()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
