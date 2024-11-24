import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QFileDialog, QHBoxLayout, QCheckBox
from PyQt5.QtCore import Qt
import pyqtgraph.opengl as gl
import uproot
import numpy as np
import pandas as pd

class EventDisplay(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Event Display")
        self.setGeometry(100, 100, 1200, 800)
        
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)
        
        # 3D View widget
        self.view = gl.GLViewWidget()
        layout.addWidget(self.view)
        
        # Add coordinate axes
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
        
        layout.addLayout(control_layout)
        
        # Initialize data storage
        self.data = None
        
    def add_coordinate_system(self):
        # Add X, Y, Z axes
        x_axis = gl.GLLinePlotItem(pos=np.array([[0,0,0], [100,0,0]]), color=(1,0,0,1), width=2)
        y_axis = gl.GLLinePlotItem(pos=np.array([[0,0,0], [0,100,0]]), color=(0,1,0,1), width=2)
        z_axis = gl.GLLinePlotItem(pos=np.array([[0,0,0], [0,0,100]]), color=(0,0,1,1), width=2)
        self.view.addItem(x_axis)
        self.view.addItem(y_axis)
        self.view.addItem(z_axis)
        
    def load_data(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Open ROOT File", "", "ROOT files (*.root)")
        if filename:
            # Load ROOT file using uproot
            file = uproot.open(filename)
            tree = file["tracking_clusters"]  # Using your tree name
            self.data = tree.arrays(library="pandas", branches=[
                'x', 'y', 'z', 
                'passed_straight',
                'used_in_seed'
            ])
            self.update_display()
            
    def update_display(self):
        if self.data is not None:
            # Clear previous display
            self.view.clear()
            self.add_coordinate_system()
            
            # Filter data based on checkboxes
            display_data = self.data.copy()
            
            if self.show_passed.isChecked():
                display_data = display_data[display_data['passed_straight'] == 1]
            if self.show_used.isChecked():
                display_data = display_data[display_data['used_in_seed'] == 1]
            
            if len(display_data) > 0:
                # Create position array
                pos = np.column_stack((
                    display_data['x'].to_numpy(),
                    display_data['y'].to_numpy(),
                    display_data['z'].to_numpy()
                ))
                
                # Color points based on their properties
                colors = np.zeros((len(pos), 4))
                colors[:] = [0, 0, 1, 0.8]  # Default blue
                
                # Add points to display
                scatter = gl.GLScatterPlotItem(pos=pos, color=colors, size=2)
                self.view.addItem(scatter)

def main():
    app = QApplication(sys.argv)
    window = EventDisplay()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()