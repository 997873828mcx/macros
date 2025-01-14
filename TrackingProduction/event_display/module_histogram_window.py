# module_histogram_window.py

from PyQt5.QtWidgets import QWidget, QVBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import numpy as np


class ModuleHistogramWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Module-Based Track Fitting Histograms")
        self.setGeometry(200, 200, 1600, 1200)

        layout = QVBoxLayout(self)
        self.figure = Figure(figsize=(16, 12))
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(self.canvas)

        # Create 3x2 grid of subplots
        self.axes = self.figure.subplots(3, 2)
        module_titles = ["Module 1", "Module 2", "Module 3"]
        for i in range(3):
            self.axes[i, 0].set_title(f"{module_titles[i]} Δrφ")
            self.axes[i, 1].set_title(f"{module_titles[i]} Δz")
            self.axes[i, 0].set_xlabel("Δrφ (cm)")
            self.axes[i, 0].set_ylabel("Counts")
            self.axes[i, 1].set_xlabel("Δz (cm)")
            self.axes[i, 1].set_ylabel("Counts")

    def update_module_histograms(self, module_hist_data):
        """
        module_hist_data: dict mapping module index (0,1,2) to tuple (delta_rphi, delta_z, track_id)
        """
        for module_idx, (delta_rphi, delta_z, track_id) in module_hist_data.items():
            ax_rphi = self.axes[module_idx, 0]
            ax_z = self.axes[module_idx, 1]

            ax_rphi.cla()
            ax_z.cla()

            ax_rphi.hist(delta_rphi, bins=50, range=(-1, 1), color="red", alpha=0.7)
            ax_rphi.set_title(f"Module {module_idx+1} Δrφ\nTrack {track_id}")
            ax_rphi.set_xlabel("Δrφ (cm)")
            ax_rphi.set_ylabel("Counts")
            mean_rphi = np.mean(delta_rphi)
            std_rphi = np.std(delta_rphi)
            ax_rphi.text(
                0.95,
                0.95,
                f"mean={mean_rphi:.4f}, std={std_rphi:.4f}",
                transform=ax_rphi.transAxes,
                ha="right",
                va="top",
                bbox=dict(facecolor="white", alpha=0.5),
            )

            ax_z.hist(delta_z, bins=50, range=(-2, 2), color="blue", alpha=0.7)
            ax_z.set_title(f"Module {module_idx+1} Δz\nTrack {track_id}")
            ax_z.set_xlabel("Δz (cm)")
            ax_z.set_ylabel("Counts")
            mean_z = np.mean(delta_z)
            std_z = np.std(delta_z)
            ax_z.text(
                0.95,
                0.95,
                f"mean={mean_z:.4f}, std={std_z:.4f}",
                transform=ax_z.transAxes,
                ha="right",
                va="top",
                bbox=dict(facecolor="white", alpha=0.5),
            )

        self.figure.tight_layout()
        self.canvas.draw()
