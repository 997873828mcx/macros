# histogram_window.py

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import matplotlib.gridspec as gridspec
import numpy as np


class HistogramWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Track Fitting Histograms")
        self.setGeometry(100, 100, 1600, 1000)

        # Create layout
        layout = QVBoxLayout(self)

        # Create figure with subplots that will expand
        self.figure = Figure(figsize=(16, 8))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)

        gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1.2])
        self.individual_ax_rphi = self.figure.add_subplot(gs[0, 0])
        self.individual_ax_z = self.figure.add_subplot(gs[0, 1])
        self.accumulated_ax_rphi = self.figure.add_subplot(gs[1, 0])
        self.accumulated_ax_z = self.figure.add_subplot(gs[1, 1])
        self.rphi_vs_r_ax = self.figure.add_subplot(gs[2, :])

        # Titles for sections
        self.individual_ax_rphi.set_title("Individual Track Delta rphi")
        self.individual_ax_z.set_title("Individual Track Delta z")
        self.accumulated_ax_rphi.set_title("Accumulated Delta rphi")
        self.accumulated_ax_z.set_title("Accumulated Delta z")
        self.rphi_vs_r_ax.set_title("Delta rphi vs. r")

        # Labels for axes
        self.individual_ax_rphi.set_xlabel("Delta rphi (cm)")
        self.individual_ax_rphi.set_ylabel("Counts")
        self.individual_ax_z.set_xlabel("Delta z (cm)")
        self.individual_ax_z.set_ylabel("Counts")
        self.accumulated_ax_rphi.set_xlabel("Delta rphi (cm)")
        self.accumulated_ax_rphi.set_ylabel("Counts")
        self.accumulated_ax_z.set_xlabel("Delta z (cm)")
        self.accumulated_ax_z.set_ylabel("Counts")
        self.rphi_vs_r_ax.set_xlabel("r (cm)")
        self.rphi_vs_r_ax.set_ylabel("Delta rphi (cm)")

        # Store all delta values for accumulation
        self.accumulated_delta_rphi = []
        self.accumulated_delta_z = []
        self.accumulated_r_values = []

    def add_histograms(self, delta_rphi, delta_z, track_id, points=None):
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
            0.95,
            0.95,
            f"mean={mean_rphi:.4f}, std={std_rphi:.4f}",
            transform=self.individual_ax_rphi.transAxes,
            ha="right",
            va="top",
            bbox=dict(facecolor="white", alpha=0.5),
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
            0.95,
            0.95,
            f"mean={mean_z:.4f}, std={std_z:.4f}",
            transform=self.individual_ax_z.transAxes,
            ha="right",
            va="top",
            bbox=dict(facecolor="white", alpha=0.5),
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
            self.accumulated_delta_rphi,
            bins=100,
            range=(-1, 1),
            color="green",
            alpha=0.7,
        )
        self.accumulated_ax_rphi.set_title("Accumulated Delta rphi")
        self.accumulated_ax_rphi.set_xlabel("Delta rphi (cm)")
        self.accumulated_ax_rphi.set_ylabel("Counts")

        # Calculate and annotate mean and std for accumulated delta_rphi
        mean_acc_rphi = np.mean(self.accumulated_delta_rphi)
        std_acc_rphi = np.std(self.accumulated_delta_rphi)
        self.accumulated_ax_rphi.text(
            0.95,
            0.95,
            f"mean={mean_acc_rphi:.4f}, std={std_acc_rphi:.4f}",
            transform=self.accumulated_ax_rphi.transAxes,
            ha="right",
            va="top",
            bbox=dict(facecolor="white", alpha=0.5),
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
            0.95,
            0.95,
            f"mean={mean_acc_z:.4f}, std={std_acc_z:.4f}",
            transform=self.accumulated_ax_z.transAxes,
            ha="right",
            va="top",
            bbox=dict(facecolor="white", alpha=0.5),
        )

        # --- Delta rphi vs. r Plot ---
        if points is not None:
            # Calculate r values from x,y coordinates
            r_values = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)

            # Debugging: Print lengths of r_values and delta_rphi
            print(f"Track {track_id}: Number of r_values = {len(r_values)}")
            print(f"Track {track_id}: Number of delta_rphi = {len(delta_rphi)}")

            if len(r_values) != len(delta_rphi):
                print(f"Error: Mismatch in lengths for Track {track_id}.")
                print(
                    f"r_values length: {len(r_values)}, delta_rphi length: {len(delta_rphi)}"
                )
                # Optionally, raise an exception or handle the mismatch
                # For now, we'll skip plotting this scatter
                self.rphi_vs_r_ax.cla()
                self.rphi_vs_r_ax.set_title(
                    f"Track {track_id} Delta rphi vs. r (Plot Skipped Due to Mismatch)"
                )
            else:
                # Clear the plot
                self.rphi_vs_r_ax.cla()

                # Create scatter plot
                self.rphi_vs_r_ax.scatter(r_values, delta_rphi, alpha=0.5, s=20)
                self.rphi_vs_r_ax.set_title(f"Track {track_id} Delta rphi vs. r")
                self.rphi_vs_r_ax.set_xlabel("r (cm)")
                self.rphi_vs_r_ax.set_ylabel("Delta rphi (cm)")

                # Add grid
                self.rphi_vs_r_ax.grid(True, linestyle="--", alpha=0.7)

                # Set reasonable y-axis limits
                self.rphi_vs_r_ax.set_ylim(-1, 1)

                # Store r values for possible future use
                self.accumulated_r_values.extend(r_values)

        # Adjust layout to prevent overlap
        self.figure.tight_layout()
        self.canvas.draw()
