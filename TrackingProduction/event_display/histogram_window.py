# histogram_window.py

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import numpy as np

class HistogramWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Track Fitting Histograms")
        self.setGeometry(100, 100, 1600, 800)

        # Create layout
        layout = QVBoxLayout(self)

        # Create figure with subplots that will expand
        self.figure = Figure(figsize=(16, 8))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)

        # Initialize subplots
        # Allocate space: individual histograms and accumulated histograms
        # For simplicity, let's stack them vertically
        self.individual_ax_rphi = self.figure.add_subplot(2, 2, 1)
        self.individual_ax_z = self.figure.add_subplot(2, 2, 2)
        self.accumulated_ax_rphi = self.figure.add_subplot(2, 2, 3)
        self.accumulated_ax_z = self.figure.add_subplot(2, 2, 4)

        # Titles for sections
        self.individual_ax_rphi.set_title("Individual Track Delta rphi")
        self.individual_ax_z.set_title("Individual Track Delta z")
        self.accumulated_ax_rphi.set_title("Accumulated Delta rphi")
        self.accumulated_ax_z.set_title("Accumulated Delta z")

        # Labels for axes
        self.individual_ax_rphi.set_xlabel("Delta rphi (cm)")
        self.individual_ax_rphi.set_ylabel("Counts")
        self.individual_ax_z.set_xlabel("Delta z (cm)")
        self.individual_ax_z.set_ylabel("Counts")
        self.accumulated_ax_rphi.set_xlabel("Delta rphi (cm)")
        self.accumulated_ax_rphi.set_ylabel("Counts")
        self.accumulated_ax_z.set_xlabel("Delta z (cm)")
        self.accumulated_ax_z.set_ylabel("Counts")

        # Store all delta values for accumulation
        self.accumulated_delta_rphi = []
        self.accumulated_delta_z = []

    def add_histograms(self, delta_rphi, delta_z, track_id):
        """Add new histograms for a new track and update accumulated histograms."""
        # --- Individual Histograms ---

        # Clear individual histograms
        self.individual_ax_rphi.cla()
        self.individual_ax_z.cla()

        # Plot individual track delta_rphi
        self.individual_ax_rphi.hist(
            delta_rphi, bins=50, range=(-1, 1), color="red", alpha=0.7
        )
        self.individual_ax_rphi.set_title(f"Track {track_id} Delta rphi")
        self.individual_ax_rphi.set_xlabel("Delta rphi (cm)")
        self.individual_ax_rphi.set_ylabel("Counts")
        
        # Calculate and annotate mean and std for delta_rphi
        mean_rphi = np.mean(delta_rphi)
        std_rphi = np.std(delta_rphi)
        self.individual_ax_rphi.text(
            0.95, 0.95, f"mean={mean_rphi:.4f}, std={std_rphi:.4f}", 
            transform=self.individual_ax_rphi.transAxes,
            ha="right", va="top", bbox=dict(facecolor='white', alpha=0.5)
        )

        # Plot individual track delta_z
        self.individual_ax_z.hist(
            delta_z, bins=50, range=(-2, 2), color="blue", alpha=0.7
        )
        self.individual_ax_z.set_title(f"Track {track_id} Delta z")
        self.individual_ax_z.set_xlabel("Delta z (cm)")
        self.individual_ax_z.set_ylabel("Counts")
        
        # Calculate and annotate mean and std for delta_z
        mean_z = np.mean(delta_z)
        std_z = np.std(delta_z)
        self.individual_ax_z.text(
            0.95, 0.95, f"mean={mean_z:.4f}, std={std_z:.4f}", 
            transform=self.individual_ax_z.transAxes,
            ha="right", va="top", bbox=dict(facecolor='white', alpha=0.5)
        )

        # --- Accumulated Histograms ---

        # Accumulate the delta values
        self.accumulated_delta_rphi.extend(delta_rphi)
        self.accumulated_delta_z.extend(delta_z)

        # Clear accumulated histograms
        self.accumulated_ax_rphi.cla()
        self.accumulated_ax_z.cla()

        # Plot accumulated delta_rphi
        self.accumulated_ax_rphi.hist(
            self.accumulated_delta_rphi, bins=100, range=(-1, 1), color="green", alpha=0.7
        )
        self.accumulated_ax_rphi.set_title("Accumulated Delta rphi")
        self.accumulated_ax_rphi.set_xlabel("Delta rphi (cm)")
        self.accumulated_ax_rphi.set_ylabel("Counts")
        
        # Calculate and annotate mean and std for accumulated delta_rphi
        mean_acc_rphi = np.mean(self.accumulated_delta_rphi)
        std_acc_rphi = np.std(self.accumulated_delta_rphi)
        self.accumulated_ax_rphi.text(
            0.95, 0.95, f"mean={mean_acc_rphi:.4f}, std={std_acc_rphi:.4f}", 
            transform=self.accumulated_ax_rphi.transAxes,
            ha="right", va="top", bbox=dict(facecolor='white', alpha=0.5)
        )

        # Plot accumulated delta_z
        self.accumulated_ax_z.hist(
            self.accumulated_delta_z, bins=100, range=(-2, 2), color="purple", alpha=0.7
        )
        self.accumulated_ax_z.set_title("Accumulated Delta z")
        self.accumulated_ax_z.set_xlabel("Delta z (cm)")
        self.accumulated_ax_z.set_ylabel("Counts")
        
        # Calculate and annotate mean and std for accumulated delta_z
        mean_acc_z = np.mean(self.accumulated_delta_z)
        std_acc_z = np.std(self.accumulated_delta_z)
        self.accumulated_ax_z.text(
            0.95, 0.95, f"mean={mean_acc_z:.4f}, std={std_acc_z:.4f}", 
            transform=self.accumulated_ax_z.transAxes,
            ha="right", va="top", bbox=dict(facecolor='white', alpha=0.5)
        )

        # Adjust layout to prevent overlap
        self.figure.tight_layout()
        self.canvas.draw()