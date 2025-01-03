from PyQt5.QtWidgets import QWidget, QVBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg


class HistogramWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Track Fitting Histograms")
        self.setGeometry(100, 100, 1200, 800)

        # Create layout
        layout = QVBoxLayout(self)

        # Create figure with subplots that will expand
        self.figure = Figure(figsize=(12, 8))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)

        # Keep track of number of tracks
        self.track_count = 0
        # Store all delta values for possible later use
        self.all_deltas = []

    def add_histograms(self, delta_rphi, delta_z, track_id):
        """Add new histograms for a new track."""
        self.track_count += 1

        # Calculate number of rows and columns needed
        n_rows = (self.track_count + 1) // 2  # 2 tracks per row
        n_cols = 4

        # Clear figure and create new subplots
        self.figure.clear()

        # Create subplots for all tracks
        for i in range(self.track_count):
            # Calculate row and column positions
            row = i // 2
            col = (i % 2) * 2  # Multiply by 2 because each track uses 2 columns
            
            # Create subplot indices
            rphi_index = row * n_cols + col + 1
            z_index = row * n_cols + col + 2
            
            # Create subplots
            ax_rphi = self.figure.add_subplot(n_rows, n_cols, rphi_index)
            ax_z = self.figure.add_subplot(n_rows, n_cols, z_index)

            if i == self.track_count - 1:
                # This is the new track
                self.all_deltas.append((delta_rphi, delta_z))
                # Plot new histograms
                ax_rphi.hist(delta_rphi, bins=50, range=(-1, 1), color="red", alpha=0.7)
                ax_rphi.set_title(f"Track {i+1} Delta rphi")
                ax_rphi.set_xlabel("Delta rphi (cm)")
                ax_rphi.set_ylabel("Counts")

                ax_z.hist(delta_z, bins=50, range=(-2, 2), color="blue", alpha=0.7)
                ax_z.set_title(f"Track {i+1} Delta z")
                ax_z.set_xlabel("Delta z (cm)")
                ax_z.set_ylabel("Counts")
            else:
                # Replot previous histograms
                old_delta_rphi, old_delta_z = self.all_deltas[i]
                ax_rphi.hist(
                    old_delta_rphi, bins=50, range=(-1, 1), color="red", alpha=0.7
                )
                ax_rphi.set_title(f"Track {i+1} Delta rphi")
                ax_rphi.set_xlabel("Delta rphi (cm)")
                ax_rphi.set_ylabel("Counts")

                ax_z.hist(old_delta_z, bins=50, range=(-2, 2), color="blue", alpha=0.7)
                ax_z.set_title(f"Track {i+1} Delta z")
                ax_z.set_xlabel("Delta z (cm)")
                ax_z.set_ylabel("Counts")

        # Adjust layout to prevent overlap
        self.figure.tight_layout()
        self.canvas.draw()
