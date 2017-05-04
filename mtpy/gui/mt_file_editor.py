# -*- coding: utf-8 -*-
"""
Created on Wed May 03 19:01:42 2017

@author: jrpeacock
"""

import sys
from PyQt5 import QtWidgets, QtGui, QtCore
from mtpy.core import mt_new as mt

# header label font
label_font = QtGui.QFont()
label_font.setBold = True
label_font.setPointSize (14)
 
class MTMainWindow(QtWidgets.QMainWindow):
 
    def __init__(self):
        super(MTMainWindow, self).__init__()
        self.setWindowTitle('MT File Editor')
 
        self.setup_ui()
        
    def setup_ui(self):
        """
        setup user interface
        """
        screen_size = QtWidgets.QDesktopWidget().availableGeometry()
        width = screen_size.width()
        
        self.setWindowState(QtCore.Qt.WindowMaximized)
        
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
    
        self.tab_widget = MTTabWidget(self)
        
        self.text_label = QtWidgets.QLabel("File Preview")
        self.text_label.setFont(label_font)
        
        self.text_edit = QtWidgets.QTextEdit()
        self.text_edit.setMaximumWidth(int(width/2.0))

        text_layout = QtWidgets.QVBoxLayout()
        text_layout.addWidget(self.text_label)
        text_layout.addWidget(self.text_edit)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.tab_widget)
        layout.addLayout(text_layout)
        
        self.central_widget.setLayout(layout)
        #self.centeral_widget = self.setCentralWidget(self.tab_widget)
        
 
class MTTabWidget(QtWidgets.QTabWidget):        
 
    def __init__(self, parent=None):   
        super(MTTabWidget, self).__init__(parent)
        
        self.setup_ui()
        
    def setup_ui(self):
 
        self.tab_site = SiteTab(self)
        self.tab_field = QtWidgets.QWidget()
        self.tab_processing = QtWidgets.QWidget()
        
        self.addTab(self.tab_site, "Site")
        self.addTab(self.tab_field, "Field Notes")
        self.addTab(self.tab_processing, "Processing")

#==============================================================================
# Site tab        
#==============================================================================
class SiteTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(SiteTab, self).__init__(parent)
        self.Site = mt.Site()
        
        self.setup_ui()
        
    def setup_ui(self):
        
        self.acq_by_label = QtWidgets.QLabel('Acquired By')
        self.acq_by_label.setFont(label_font)
        self.acq_by_edit = QtWidgets.QLineEdit()
        self.acq_by_edit.setText(self.Site.acquired_by)
        self.acq_by_edit.editingFinished.connect(self.set_acq_by)
        
        self.id_label = QtWidgets.QLabel('Station ID')
        self.id_label.setFont(label_font)
        self.id_edit = QtWidgets.QLineEdit()
        self.id_edit.setText(self.Site.id)
        self.id_edit.editingFinished.connect(self.set_id)
        
        self.lat_label = QtWidgets.QLabel('Station Latitude (decimal or HH:MM:SS.ms)')
        self.lat_label.setFont(label_font)
        self.lat_edit = QtWidgets.QLineEdit()
        self.lat_edit.setText(self.Site.Location.latitude)
        self.lat_edit.editingFinished.connect(self.set_lat)
        
        self.lon_label = QtWidgets.QLabel('Station Longitude (decimal or HH:MM:SS.ms)')
        self.lon_label.setFont(label_font)
        self.lon_edit = QtWidgets.QLineEdit()
        self.lon_edit.setText(self.Site.Location.longitude)
        self.lon_edit.editingFinished.connect(self.set_lon)
        
        self.datum_label = QtWidgets.QLabel('Geographic Datum')
        self.datum_label.setFont(label_font)
        self.datum_combo = QtWidgets.QComboBox()
        self.datum_combo.addItems(['WGS84',
                                   'NAD83',
                                   'AGD84',
                                   'GRS80',
                                   'NAD27',
                                   'ED50',
                                   'JGD2011',
                                   'KGD97',
                                   'GCJ02',
                                   'BD09',
                                   'GTRF',])
        self.datum_combo.activated[str].connect(self.set_datum)
        
        self.elev_label = QtWidgets.QLabel('Station Elevation')
        self.elev_label.setFont(label_font)
        self.elev_edit = QtWidgets.QLineEdit()
        self.elev_edit.setText(self.Site.Location.elevation)
        self.elev_edit.editingFinished.connect(self.set_elev)
        
        self.elev_units_label = QtWidgets.QLabel('Elevation Units')
        self.elev_units_label.setFont(label_font)
        self.elev_units_combo = QtWidgets.QComboBox()
        self.elev_units_combo.addItems(['km', 'm', 'ft', 'miles'])
        self.elev_units_combo.activated[str].connect(self.set_elev_units)
        
        self.coord_label = QtWidgets.QLabel('Z and T Coordinate System')
        self.coord_label.setFont(label_font)
        self.coord_combo = QtWidgets.QComboBox()
        self.coord_combo.addItems(['Geographic North', 'Geomagnetic North'])
        self.coord_combo.activated[str].connect(self.set_coord)
        
        self.dec_label = QtWidgets.QLabel('Geomagnetic Declination (deg)')
        self.dec_label.setFont(label_font)
        self.dec_edit = QtWidgets.QLineEdit()
        self.dec_edit.setText(self.Site.Location.declination)
        self.dec_edit.editingFinished.connect(self.set_dec)
        
        self.proj_label = QtWidgets.QLabel('Project Name')
        self.proj_label.setFont(label_font)
        self.proj_edit = QtWidgets.QLineEdit()
        self.proj_edit.setText(self.Site.project)
        self.proj_edit.editingFinished.connect(self.set_proj)
        
        self.survey_label = QtWidgets.QLabel('Survey Name')
        self.survey_label.setFont(label_font)
        self.survey_edit = QtWidgets.QLineEdit()
        self.survey_edit.setText(self.Site.survey)
        self.survey_edit.editingFinished.connect(self.set_survey)
        
        
        layout = QtWidgets.QFormLayout()
        layout.addRow(self.proj_label, self.proj_edit)
        layout.addRow(self.survey_label, self.survey_edit)
        layout.addRow(self.id_label, self.id_edit)
        layout.addRow(self.lat_label, self.lat_edit)
        layout.addRow(self.lon_label, self.lon_edit)
        layout.addRow(self.datum_label, self.datum_combo)
        layout.addRow(self.coord_label, self.coord_combo)
        layout.addRow(self.elev_label, self.elev_edit)
        layout.addRow(self.elev_units_label, self.elev_units_combo)
        layout.addRow(self.acq_by_label, self.acq_by_edit)
        layout.addRow(self.dec_label, self.dec_edit)
        
        
        self.setLayout(layout)
        
    def set_acq_by(self):
        self.Site.acquired_by = str(self.acq_by_edit.text())
        self.acq_by_edit.setText(self.Site.acquired_by)
        
    def set_id(self):
        self.Site.id = str(self.id_edit.text())
        self.id_edit.setText(self.Site.id)
        
    def set_proj(self):
        self.Site.project = str(self.id_edit.text())
        self.proj_edit.setText(self.Site.project) 
        
    def set_survey(self):
        self.Site.survey = str(self.survey_edit.text())
        self.survey_edit.setText(self.Site.survey)
        
    def set_lat(self):
        self.Site.Location.latitude = str(self.lat_edit.text())
        self.lat_edit.setText('{0:.6f}'.format(self.Site.Location.latitude))
    
    def set_lon(self):
        self.Site.Location.longitude = str(self.lon_edit.text())
        self.lon_edit.setText('{0:.6f}'.format(self.Site.Location.longitude))
    
    def set_elev(self):
        self.Site.Location.elevation = str(self.elev_edit.text())
        self.elev_edit.setText('{0:.3f}'.format(self.Site.Location.elevation))
        
    def set_elev_units(self, text):
        self.Site.Location.elev_units = text
        
    def set_datum(self, text):
        self.Site.Location.datum = text
        
    def set_dec(self):
        try:
            self.Site.Location.declination = float(self.dec_edit.text())
            self.dec_edit.setText('{0:.3f}'.format(self.Site.Location.declination))
        except ValueError:
            self.Site.Location.declination = 0.0
            self.dec_edit.setText('{0:.3f}'.format(self.Site.Location.declination))
        
    def set_coord(self, text):
        self.Site.Location.coordinate_system = text
        
        
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = MTMainWindow()
    ex.show()
    sys.exit(app.exec_())