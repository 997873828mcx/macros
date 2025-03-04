# pad_adc_window.py

from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QComboBox,
    QLabel,
    QPushButton,
    QGroupBox,
)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import numpy as np
import matplotlib.gridspec as gridspec


class PadADCWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ADC Distribution Analysis")
        self.setGeometry(100, 100, 1200, 800)

        # Create main layout
        main_layout = QVBoxLayout(self)

        # Create controls layout
        controls_layout = QHBoxLayout()

        # View type selector
        view_group = QGroupBox("View Options")
        view_layout = QHBoxLayout()
        self.view_selector = QComboBox()
        self.view_selector.addItems(
            ["2D Heatmap", "Pad Distribution", "Time Distribution", "3D Surface"]
        )
        self.view_selector.currentIndexChanged.connect(self.update_view)
        view_layout.addWidget(QLabel("View Type:"))
        view_layout.addWidget(self.view_selector)
        view_group.setLayout(view_layout)
        controls_layout.addWidget(view_group)

        # Add clear button
        self.clear_button = QPushButton("Clear Selection")
        self.clear_button.clicked.connect(self.clear_data)
        controls_layout.addWidget(self.clear_button)

        main_layout.addLayout(controls_layout)

        # Create figure with subplots
        self.figure = Figure(figsize=(12, 8))

        self.canvas = FigureCanvasQTAgg(self.figure)
        main_layout.addWidget(self.canvas)
        self.setup_axes()

        # Initialize data storage
        self.pad_values = []
        self.time_values = []
        self.adc_values = []

        # Set initial view
        self.current_view = "2D Heatmap"

    def setup_axes(self):
        """Set up the subplot axes."""
        # Clear existing axes
        self.figure.clear()

        # Create GridSpec
        self.gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], figure=self.figure)

        # Create main plot and projections
        self.ax_main = self.figure.add_subplot(self.gs[0, :])
        self.ax_pad = self.figure.add_subplot(self.gs[1, 0])
        self.ax_time = self.figure.add_subplot(self.gs[1, 1])

        # Adjust layout
        self.figure.tight_layout()

    def update_distribution(self, selected_hits):
        """
        Update the ADC distribution plots with new hit data.

        Args:
            selected_hits: List of dictionaries containing hit information
                         Each dict should have 'pad', 'time', and 'adc' keys
        """
        # Store the new data

        if not selected_hits:
            return
        self.pad_values = [hit.get("pad", 0) for hit in selected_hits]
        self.time_values = [hit.get("time", 0) for hit in selected_hits]
        self.adc_values = [hit.get("adc", 0) for hit in selected_hits]

        self.update_view()

    def update_view(self):
        """Update the visualization based on the selected view type."""
        view_type = self.view_selector.currentText()

        """# Clear all axes
        self.ax_main.cla()
        self.ax_pad.cla()
        self.ax_time.cla()"""

        if not self.pad_values:
            # self.canvas.draw()
            return

        self.setup_axes()

        if view_type == "2D Heatmap":
            self._plot_2d_heatmap()
        elif view_type == "Pad Distribution":
            self._plot_pad_distribution()
        elif view_type == "Time Distribution":
            self._plot_time_distribution()
        elif view_type == "3D Surface":
            self._plot_3d_surface()

        # Update the canvas
        # self.figure.tight_layout()
        self.canvas.draw()

    def _plot_2d_heatmap(self):
        """Create a 2D heatmap of ADC values."""
        # Create 2D histogram
        pad_bins = np.linspace(min(self.pad_values), max(self.pad_values), 30)
        time_bins = np.linspace(min(self.time_values), max(self.time_values), 30)

        hist2d = np.histogram2d(
            self.pad_values,
            self.time_values,
            bins=[pad_bins, time_bins],
            weights=self.adc_values,
        )[0]

        # Plot heatmap
        im = self.ax_main.imshow(
            hist2d.T,
            aspect="auto",
            origin="lower",
            extent=[
                min(self.pad_values),
                max(self.pad_values),
                min(self.time_values),
                max(self.time_values),
            ],
            cmap="viridis",
        )
        self.figure.colorbar(im, ax=self.ax_main, label="ADC Value")

        # Set labels
        self.ax_main.set_xlabel("Pad Number")
        self.ax_main.set_ylabel("Time Bin")
        self.ax_main.set_title("ADC Distribution")

        # Plot projections
        self._plot_projections()

    def _plot_3d_surface(self):
        """Create a 3D surface plot of ADC values."""
        # self.ax_main.remove()
        self.figure.clear()
        self.gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], figure=self.figure)
        self.ax_pad = self.figure.add_subplot(self.gs[1, 0])
        self.ax_time = self.figure.add_subplot(self.gs[1, 1])

        # Create grid for surface plot
        pad_bins = np.linspace(min(self.pad_values), max(self.pad_values), 30)
        time_bins = np.linspace(min(self.time_values), max(self.time_values), 30)

        hist2d = np.histogram2d(
            self.pad_values,
            self.time_values,
            bins=[pad_bins, time_bins],
            weights=self.adc_values,
        )[0]

        pad_centers = (pad_bins[:-1] + pad_bins[1:]) / 2
        time_centers = (time_bins[:-1] + time_bins[1:]) / 2

        PAD, TIME = np.meshgrid(pad_centers, time_centers)

        # Plot surface
        surf = self.ax_main.plot_surface(
            PAD, TIME, hist2d.T, cmap="viridis", edgecolor="none"
        )

        self.figure.colorbar(surf, ax=self.ax_main, label="ADC Value")

        # Set labels
        self.ax_main.set_xlabel("Pad Number")
        self.ax_main.set_ylabel("Time Bin")
        self.ax_main.set_zlabel("ADC Value")
        self.ax_main.set_title("ADC Distribution 3D Surface")

        # Plot projections
        self._plot_projections()

    def _plot_pad_distribution(self):
        """Create a 1D histogram of ADC values vs pad numbers."""
        self.ax_main.hist(self.pad_values, weights=self.adc_values, bins=30, alpha=0.7)
        self.ax_main.set_xlabel("Pad Number")
        self.ax_main.set_ylabel("ADC Value")
        self.ax_main.set_title("ADC Distribution by Pad Number")
        self.ax_main.grid(True, linestyle="--", alpha=0.7)

        # Plot projections
        self._plot_projections()

    def _plot_time_distribution(self):
        """Create a 1D histogram of ADC values vs time bins."""
        self.ax_main.hist(self.time_values, weights=self.adc_values, bins=30, alpha=0.7)
        self.ax_main.set_xlabel("Time Bin")
        self.ax_main.set_ylabel("ADC Value")
        self.ax_main.set_title("ADC Distribution by Time Bin")
        self.ax_main.grid(True, linestyle="--", alpha=0.7)

        # Plot projections
        self._plot_projections()

    def _plot_projections(self):
        """Plot 1D projections for pad and time distributions."""
        # Pad projection
        self.ax_pad.hist(self.pad_values, weights=self.adc_values, bins=30)
        self.ax_pad.set_xlabel("Pad Number")
        self.ax_pad.set_ylabel("Total ADC")
        self.ax_pad.set_title("Pad Projection")
        self.ax_pad.grid(True, linestyle="--", alpha=0.7)

        # Time projection
        self.ax_time.hist(self.time_values, weights=self.adc_values, bins=30)
        self.ax_time.set_xlabel("Time Bin")
        self.ax_time.set_ylabel("Total ADC")
        self.ax_time.set_title("Time Projection")
        self.ax_time.grid(True, linestyle="--", alpha=0.7)

    def clear_data(self):
        """Clear all stored data and reset the plots."""
        self.pad_values = []
        self.time_values = []
        self.adc_values = []
        self.setup_axes()

        self.canvas.draw()
