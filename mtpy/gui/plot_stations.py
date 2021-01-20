# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 22:47:32 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""

# =============================================================================
# Imports
# =============================================================================
import sys
try:
    from PyQt5 import QtCore, QtWidgets, QtGui
except ImportError:
    raise ImportError("This version needs PyQt5")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
import matplotlib.widgets as mplwidgets

import mtpy.imaging.mtplottools as mtplottools
import mtpy.modeling.modem as modem

# =============================================================================
# Plot stations 
# =============================================================================

class PlotStations(QtWidgets.QWidget):
    """
    plot station locations
    """
    
    def __init__(self, station_locations):
        self.station_locations = station_locations
        self.plot_crs = None
        
        super().__init__()
        self.setup_ui()
        
    def setup_ui(self):
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        self.mpl_widget.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.mpl_widget.setFocus()

        # be able to edit the data
        self.mpl_widget.mpl_connect("pick_event", self.on_pick)
        self.mpl_widget.mpl_connect("axes_enter_event", self.in_axes)
        self.mpl_widget.mpl_connect("button_press_event", self.on_pick)

        # make sure the figure takes up the entire plottable space
        self.mpl_widget.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
        )

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.mpl_toolbar = NavigationToolbar(self.mpl_widget, self)

        # set the layout for the plot
        mpl_vbox = QtWidgets.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.mpl_widget)
        
        self.setLayout(mpl_vbox)
        self.mpl_widget.updateGeometry()
        
    def on_pick(self):
        pass
    
    def in_axes(self):
        pass
 
# ==============================================================================
# Def Main
# ==============================================================================
def main():
    app = QtWidgets.QApplication(sys.argv)
    ui = PlotStations(None)
    ui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()       