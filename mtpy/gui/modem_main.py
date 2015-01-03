# -*- coding: utf-8 -*-
"""
Created on Sun Nov 02 13:47:10 2014

@author: jrpeacock
"""

from PyQt4 import QtCore, QtGui
import mtpy.modeling.modem_new as modem
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import os
import mtpy.analysis.pt as mtpt
import mtpy.utils.exceptions as mtex
from matplotlib.colors import Normalize
import matplotlib.colorbar as mcb
import mtpy.imaging.mtcolors as mtcl
import mtpy.gui.get_edi_files as mt_get_edi_files

class ModEM_Main_Window(QtGui.QMainWindow):
    """
    main window for building a mesh for ModEM
    
    """
    
    def __init__(self):
        super(ModEM_Main_Window, self).__init__()
        
        self.ui_setup()
        
    def ui_setup(self):
        
        self.setWindowTitle("Make ModEM Data and Model Files")
        self.resize(1920, 1080)
        
        #make a widget that will be the station list
        self.list_widget = QtGui.QListWidget()
        self.list_widget.itemClicked.connect(self.get_period)
        self.list_widget.setMaximumWidth(150)
        
        # add a menu bar to the top of the window        
        self.menu_bar = QtGui.QMenuBar(self)
        self.menu_bar.setGeometry(QtCore.QRect(0, 0, 1920, 38))
        self.menu_bar.setObjectName("menu_bar")
        
        #----------> DATA -------------------------------------
        # add a tab for getting data file.
        self.menu_data = QtGui.QMenu(self.menu_bar)
        self.menu_data.setTitle("Data File")
        
        # create an action to the data tab
        self.action_data_open = QtGui.QAction(self)
        self.action_data_open.setText('Open')
        self.action_data_open.triggered.connect(self.get_data_fn)

        # add the action to the data tap
        self.menu_data.addAction(self.action_data_open)
        self.menu_bar.addAction(self.menu_data.menuAction())
        
        #----------> Model -------------------------------------
        self.menu_model = QtGui.QMenu(self.menu_bar)
        self.menu_model.setTitle("Model File")
        
        # create an action to the model tab
        self.action_model_open = QtGui.QAction(self)
        self.action_model_open.setText('Open')
        self.action_model_open.triggered.connect(self.get_model_fn)

        # add the action to the model tap
        self.menu_model.addAction(self.action_model_open)
        self.menu_bar.addAction(self.menu_model.menuAction())
        
        self.setMenuBar(self.menu_bar)
        
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def get_data_fn(self):
        """
        get the filename from a file dialogue
        
        """        

        fn_dialog = QtGui.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                       filter='(*.dat);; (*.data)'))
                                       
        self.modem_data = modem.Data()
        self.modem_data.read_data_file(fn)
        self.modem_data_fn = fn
        
        self.dir_path = os.path.dirname(fn)
        
        self.period_list = sorted(self.modem_data.period_list)
        self.period_dict = dict([('{0:.5f}'.format(key), value) for value, key
                                 in enumerate(self.period_list)])
        
        self.list_widget.clear()
        
        #this will add the station name for each station to the qwidget list
        for period in self.period_list:
            self.list_widget.addItem('{0:.5f}'.format(period))
            
        self.plot_period = self.period_list[0]
        
    def get_period(self, widget_item):
        """
        get the station name from the clicked station 
        """
        self.plot_period = str(widget_item.text()) 
        
#def main():
    
def main():
#if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = ModEM_Main_Window()
    ui.ui_setup()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()
                
        
        
        
    