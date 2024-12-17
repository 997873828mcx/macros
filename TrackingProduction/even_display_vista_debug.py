import sys
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QWidget,
    QPushButton,
    QFileDialog,
    QHBoxLayout,
    QTextEdit,
)
from PyQt5.QtCore import Qt
import pyqtgraph.opengl as gl
import uproot
import numpy as np
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
        # Custom method in the parent widget (if needed)
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

        # Left side: 3D view
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)

        # 3D View widget using the subclassed GLViewWidget
        self.view = MyGLViewWidget(self)
        self.view.setBackgroundColor("w")
        self.view.setFocusPolicy(Qt.StrongFocus)
        left_layout.addWidget(self.view)

        # Initialize lists to keep track of items
        self.static_items = []
        self.dynamic_items = []

        # Add coordinate axes and TPC disks
        self.add_coordinate_system()

        # Button to load data
        control_layout = QHBoxLayout()
        btn_load = QPushButton("Load ROOT File")
        btn_load.clicked.connect(self.load_data)
        control_layout.addWidget(btn_load)
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
        self.full_data = None  # Store full data (if any attributes beyond coordinates)

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

        self.view.addItem(x_axis)
        self.view.addItem(y_axis)
        self.view.addItem(z_axis)
        self.static_items.extend([x_axis, y_axis, z_axis])

        # Add TPC disks
        radius = 80
        points = 100
        theta = np.linspace(0, 2 * np.pi, points)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        z1 = -105
        z2 = 105

        verts1 = np.vstack([x, y, np.ones_like(x) * z1]).T
        verts2 = np.vstack([x, y, np.ones_like(x) * z2]).T

        verts1 = np.vstack([verts1, [0, 0, z1]])
        verts2 = np.vstack([verts2, [0, 0, z2]])

        faces = []
        for i in range(points - 1):
            faces.append([i, i + 1, points])
        faces.append([points - 1, 0, points])
        faces = np.array(faces, dtype=np.int32)

        disk1 = gl.GLMeshItem(
            vertexes=verts1,
            faces=faces,
            smooth=True,
            color=(0.5, 0.5, 0.5, 0.3),
            shader="shaded",
        )
        disk2 = gl.GLMeshItem(
            vertexes=verts2,
            faces=faces,
            smooth=True,
            color=(0.5, 0.5, 0.5, 0.3),
            shader="shaded",
        )

        self.view.addItem(disk1)
        self.view.addItem(disk2)
        self.static_items.extend([disk1, disk2])

    def load_data(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root)"
        )
        if filename:
            try:
                # Load the ROOT file
                file = uproot.open(filename)

                # Debug: Print available trees
                print("\nAvailable trees in file:")
                for key in file.keys():
                    print(f"- {key}")

                # Try to get the hit tree
                hits_tree = file["hittree"]  # Changed from "combined_hits" to "hittree"
                hits_data = hits_tree.arrays(library="pd")

                # Debug information
                print("\nDebug Information:")
                print("Available branches:", hits_data.columns.tolist())
                print("Data shape:", hits_data.shape)
                print("\nSample of raw data:")
                print(hits_data.head())
                print("\nData types:")
                print(hits_data.dtypes)

                # Check if required columns exist
                required_columns = ["gx", "gy", "gz"]
                if not all(col in hits_data.columns for col in required_columns):
                    missing = [
                        col for col in required_columns if col not in hits_data.columns
                    ]
                    error_msg = f"Missing required columns: {missing}"
                    print(error_msg)
                    self.info_panel.setText(error_msg)
                    return

                # Check for NaN values
                if hits_data[required_columns].isna().any().any():
                    print("\nWarning: NaN values found in coordinates:")
                    print("NaN counts:")
                    print(hits_data[required_columns].isna().sum())

                    # Remove rows with NaN values
                    hits_data = hits_data.dropna(subset=required_columns)
                    print("\nShape after removing NaN values:", hits_data.shape)

                # Verify data ranges
                if len(hits_data) > 0:
                    print("\nCoordinate ranges:")
                    for col in required_columns:
                        print(
                            f"{col} range: {hits_data[col].min():.2f} to {hits_data[col].max():.2f}"
                        )

                    self.data = hits_data
                    self.update_display()
                else:
                    self.info_panel.setText(
                        "No valid data points found after filtering"
                    )

            except KeyError as e:
                error_msg = f"Could not find tree or branch: {str(e)}"
                print(error_msg)
                self.info_panel.setText(error_msg)
            except Exception as e:
                import traceback

                error_msg = f"Error loading file:\n{str(e)}\n\nTraceback:\n{traceback.format_exc()}"
                print(error_msg)
                self.info_panel.setText(error_msg)

    def update_display(self):
        if self.data is not None:
            # Remove previous dynamic items
            for item in self.dynamic_items:
                self.view.removeItem(item)
            self.dynamic_items = []

            if len(self.data) > 0:
                pos = np.column_stack(
                    (
                        self.data["gx"].to_numpy(),
                        self.data["gy"].to_numpy(),
                        self.data["gz"].to_numpy(),
                    )
                )

                self.pos_data = pos
                self.full_data = self.data

                # Print data range for debugging (optional)
                print("X range:", np.min(pos[:, 0]), "to", np.max(pos[:, 0]))
                print("Y range:", np.min(pos[:, 1]), "to", np.max(pos[:, 1]))
                print("Z range:", np.min(pos[:, 2]), "to", np.max(pos[:, 2]))

                # Color all points blue
                colors = np.zeros((len(pos), 4))
                colors[:] = [0, 0, 0.7, 1.0]  # Blue

                scatter = gl.GLScatterPlotItem(
                    pos=pos, color=colors, size=5, pxMode=True, glOptions="opaque"
                )
                self.view.addItem(scatter)
                self.dynamic_items.append(scatter)
                self.scatter = scatter

                # Center view on the data
                self.view.opts["center"] = pg.Vector(
                    np.mean(pos[:, 0]),
                    np.mean(pos[:, 1]),
                    np.mean(pos[:, 2]),
                )

    def on_mouse_click(self, event):
        if self.scatter is not None and self.pos_data is not None:
            # Just an example of picking logic (if needed)
            # Here we can implement the same projection-based picking logic if you want.
            # Since the user asked only for visualization, you can omit this part or leave it as-is.

            self.view.makeCurrent()

            pos = event.pos()
            mouse_x, mouse_y = pos.x(), pos.y()

            projection = glGetDoublev(GL_PROJECTION_MATRIX)
            modelview = glGetDoublev(GL_MODELVIEW_MATRIX)
            viewport = glGetIntegerv(GL_VIEWPORT)
            win_y = viewport[3] - mouse_y

            window_coords = np.array(
                [
                    gluProject(
                        point[0], point[1], point[2], modelview, projection, viewport
                    )
                    for point in self.pos_data
                ]
            )

            distances = np.sqrt(
                (window_coords[:, 0] - mouse_x) ** 2
                + (window_coords[:, 1] - win_y) ** 2
            )
            nearest_idx = np.argmin(distances)
            min_distance = distances[nearest_idx]

            if min_distance < 10:
                point_info = self.full_data.iloc[nearest_idx]
                info_text = "Hit Information:\n\n"
                for col in point_info.index:
                    value = point_info[col]
                    if isinstance(value, float):
                        info_text += f"{col:<15}: {value:.3f}\n"
                    else:
                        info_text += f"{col:<15}: {value}\n"
                self.info_panel.setText(info_text)


def main():
    app = QApplication(sys.argv)
    window = EventDisplay()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
