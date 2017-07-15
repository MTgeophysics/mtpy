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
label_font.setPointSize (13)


 
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
        self.tab_field = FieldNotesTab(self)
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
        
        self.easting_label = QtWidgets.QLabel('Easting (m)')
        self.easting_label.setFont(label_font)
        self.easting_edit = QtWidgets.QLineEdit()
        if self.Site.Location.easting is None:
            self.easting_edit.setText('{0:.3f}'.format(0.0))
        else:
            self.easting_edit.setText('{0:.3f}'.format(self.Site.Location.easting))
        
        self.northing_label = QtWidgets.QLabel('Northing (m)')
        self.northing_label.setFont(label_font)
        self.northing_edit = QtWidgets.QLineEdit()
        if self.Site.Location.northing is None:
            self.northing_edit.setText('{0:.3f}'.format(0.0))
        else:
            self.northing_edit.setText('{0:.3f}'.format(self.Site.Location.northing))
        
        
        self.utm_label = QtWidgets.QLabel('UTM Zone')
        self.utm_label.setFont(label_font)
        self.utm_edit = QtWidgets.QLineEdit()
        if self.Site.Location.northing is None:
            self.utm_edit.setText('00A')
        else:
            self.utm_edit.setText(self.Site.Location.utm_zone)
        
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
        layout.addRow(self.easting_label, self.easting_edit)
        layout.addRow(self.northing_label, self.northing_edit)
        layout.addRow(self.utm_label, self.utm_edit)
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
        self._set_utm()
    
    def set_lon(self):
        self.Site.Location.longitude = str(self.lon_edit.text())
        self.lon_edit.setText('{0:.6f}'.format(self.Site.Location.longitude))
        self._set_utm()
        
    def set_easting(self):
        self.Site.Location.easting = float(self.easting_edit.text())
        self.easting_edit.setText('{0:.3f}'.format(self.Site.Location.easting))
        self._set_ll()
        
    def set_northing(self):
        self.Site.Location.northing = float(self.northing_edit.text())
        self.northing_edit.setText('{0:.3f}'.format(self.Site.Location.northing))
        self._set_ll()
        
    def set_utm(self):
        self.Site.Location.utm_zone = str(self.utm_edit.text())
        self.utm_edit.setText(self.Site.Location.utm_zone)
        self._set_ll()
    
    def _set_ll(self):
        if self.Site.Location.easting is not None and \
           self.Site.Location.northing is not None and \
           self.Site.Location.utm_zone is not None:
            self.Site.Location.project_location2ll()
            self.lat_edit.setText(self.Site.Location.latitude)
            self.lon_edit.setText(self.Site.Location.longitude)
        else:
            print self.Site.Location.easting, self.Site.Location.northing, self.Site.Location.utm_zone
            
    def _set_utm(self):
        if self.Site.Location.latitude is not None and \
           self.Site.Location.longitude is not None:
            self.Site.Location.project_location2utm()
            self.easting_edit.setText('{0:.3f}'.format(self.Site.Location.easting))
            self.northing_edit.setText('{0:.3f}'.format(self.Site.Location.northing))
            self.utm_edit.setText(self.Site.Location.utm_zone)
    
    def set_elev(self):
        self.Site.Location.elevation = str(self.elev_edit.text())
        self.elev_edit.setText('{0:.3f}'.format(self.Site.Location.elevation))
        
    def set_elev_units(self, text):
        self.Site.Location.elev_units = text
        
    def set_datum(self, text):
        self.Site.Location.datum = text
        self._set_utm()
        
    def set_dec(self):
        try:
            self.Site.Location.declination = float(self.dec_edit.text())
            self.dec_edit.setText('{0:.3f}'.format(self.Site.Location.declination))
        except ValueError:
            self.Site.Location.declination = 0.0
            self.dec_edit.setText('{0:.3f}'.format(self.Site.Location.declination))
        
    def set_coord(self, text):
        self.Site.Location.coordinate_system = text
        

#==============================================================================
# Field notes tab
#==============================================================================
class FieldNotesTab(QtWidgets.QWidget):
    """
    Tab to hold field notes
    """      
    def __init__(self, parent=None):
        super(FieldNotesTab, self).__init__(parent)
        self.FieldNotes = mt.FieldNotes()
        self._chn_list = ['{0:d}'.format(ii) for ii in range(1, 7, 1)]
        self._rating_list = ['{0:d}'.format(ii) for ii in range(11)]
        
        self.setup_ui()
        
    def setup_ui(self):
        self.data_logger_label = QtWidgets.QLabel('Data Logger')
        self.data_logger_label.setFont(label_font)
        
        self.data_logger_id_label = QtWidgets.QLabel('ID')
        self.data_logger_id_edit = QtWidgets.QLineEdit(self.FieldNotes.data_logger.id)
        self.data_logger_id_edit.editingFinished.connect(self.set_dl_id)
        
        self.data_logger_man_label = QtWidgets.QLabel('Manufacturer')
        self.data_logger_man_edit = QtWidgets.QLineEdit(self.FieldNotes.data_logger.manufacturer)
        self.data_logger_man_edit.editingFinished.connect(self.set_dl_man)
        
        self.data_logger_type_label = QtWidgets.QLabel('Type')
        self.data_logger_type_edit = QtWidgets.QLineEdit(self.FieldNotes.data_logger.type)
        self.data_logger_type_edit.editingFinished.connect(self.set_dl_type)
        
        
        #--> EX 
        self.ex_label = QtWidgets.QLabel('Electrode EX')
        self.ex_label.setFont(label_font)
        
        self.ex_id_label = QtWidgets.QLabel('ID')
        self.ex_id_edit = QtWidgets.QLineEdit('{0}'.format(self.FieldNotes.electrode_ex.id))
        self.ex_id_edit.editingFinished.connect(self.set_ex_id)
        
        self.ex_man_label = QtWidgets.QLabel('Manufacturer')
        self.ex_man_edit = QtWidgets.QLineEdit(self.FieldNotes.electrode_ex.manufacturer)
        self.ex_man_edit.editingFinished.connect(self.set_ex_man)
        
        self.ex_type_label = QtWidgets.QLabel('Type')
        self.ex_type_edit = QtWidgets.QLineEdit(self.FieldNotes.electrode_ex.type)
        self.ex_type_edit.editingFinished.connect(self.set_ex_type)

        self.ex_x_label = QtWidgets.QLabel("X (m)")
        self.ex_x_edit = QtWidgets.QLineEdit()
        
        self.ex_y_label = QtWidgets.QLabel("Y (m)")
        self.ex_y_edit = QtWidgets.QLineEdit()
        
        self.ex_x2_label = QtWidgets.QLabel("X2 (m)")
        self.ex_x2_edit = QtWidgets.QLineEdit()
        
        self.ex_y2_label = QtWidgets.QLabel("Y2 (m)")
        self.ex_y2_edit = QtWidgets.QLineEdit()
        
        self.ex_acqchn_label = QtWidgets.QLabel("Acq. Channel")
        self.ex_acqchn_combo = QtWidgets.QComboBox()
        self.ex_acqchn_combo.addItems(self._chn_list)
        
        ##--> EY
        self.ey_label = QtWidgets.QLabel('Electrode EY')
        self.ey_label.setFont(label_font)
        
        self.ey_id_label = QtWidgets.QLabel('ID')
        self.ey_id_edit = QtWidgets.QLineEdit('{0}'.format(self.FieldNotes.electrode_ey.id))
        self.ey_id_edit.editingFinished.connect(self.set_ey_id)
        
        self.ey_man_label = QtWidgets.QLabel('Manufacturer')
        self.ey_man_edit = QtWidgets.QLineEdit(self.FieldNotes.electrode_ey.manufacturer)
        self.ey_man_edit.editingFinished.connect(self.set_ey_man)
        
        self.ey_type_label = QtWidgets.QLabel('Type')
        self.ey_type_edit = QtWidgets.QLineEdit(self.FieldNotes.electrode_ey.type)
        self.ey_type_edit.editingFinished.connect(self.set_ey_type)

        self.ey_x_label = QtWidgets.QLabel("X (m)")
        self.ey_x_edit = QtWidgets.QLineEdit()
        
        self.ey_y_label = QtWidgets.QLabel("Y (m)")
        self.ey_y_edit = QtWidgets.QLineEdit()
        
        self.ey_x2_label = QtWidgets.QLabel("X2 (m)")
        self.ey_x2_edit = QtWidgets.QLineEdit()
        
        self.ey_y2_label = QtWidgets.QLabel("Y2 (m)")
        self.ey_y2_edit = QtWidgets.QLineEdit()
        
        self.ey_acqchn_label = QtWidgets.QLabel("Acq. Channel")
        self.ey_acqchn_combo = QtWidgets.QComboBox()
        self.ey_acqchn_combo.addItems(self._chn_list)
        
        ##--> HX
        self.hx_label = QtWidgets.QLabel('Magnetometer HX')
        self.hx_label.setFont(label_font)
        
        self.hx_id_label = QtWidgets.QLabel('ID')
        self.hx_id_edit = QtWidgets.QLineEdit('{0}'.format(self.FieldNotes.magnetometer_hx.id))
        self.hx_id_edit.editingFinished.connect(self.set_hx_id)
        
        self.hx_man_label = QtWidgets.QLabel('Manufacturer')
        self.hx_man_edit = QtWidgets.QLineEdit(self.FieldNotes.magnetometer_hx.manufacturer)
        self.hx_man_edit.editingFinished.connect(self.set_hx_man)
        
        self.hx_type_label = QtWidgets.QLabel('Type')
        self.hx_type_edit = QtWidgets.QLineEdit(self.FieldNotes.magnetometer_hx.type)
        self.hx_type_edit.editingFinished.connect(self.set_hx_type)
        
        self.hx_x_label = QtWidgets.QLabel("X (m)")
        self.hx_x_edit = QtWidgets.QLineEdit()
        
        self.hx_y_label = QtWidgets.QLabel("Y (m)")
        self.hx_y_edit = QtWidgets.QLineEdit()
        
        self.hx_azm_label = QtWidgets.QLabel("Azimtuh (deg)")
        self.hx_azm_edit = QtWidgets.QLineEdit()
        
        self.hx_acqchn_label = QtWidgets.QLabel("Acq. Channel")
        self.hx_acqchn_combo = QtWidgets.QComboBox()
        self.hx_acqchn_combo.addItems(self._chn_list)
        
        ##--> HY
        self.hy_label = QtWidgets.QLabel('Magnetometer HY')
        self.hy_label.setFont(label_font)
        
        self.hy_id_label = QtWidgets.QLabel('ID')
        self.hy_id_edit = QtWidgets.QLineEdit('{0}'.format(self.FieldNotes.magnetometer_hy.id))
        self.hy_id_edit.editingFinished.connect(self.set_hy_id)
        
        self.hy_man_label = QtWidgets.QLabel('Manufacturer')
        self.hy_man_edit = QtWidgets.QLineEdit(self.FieldNotes.magnetometer_hy.manufacturer)
        self.hy_man_edit.editingFinished.connect(self.set_hy_man)
        
        self.hy_type_label = QtWidgets.QLabel('Type')
        self.hy_type_edit = QtWidgets.QLineEdit(self.FieldNotes.magnetometer_hy.type)
        self.hy_type_edit.editingFinished.connect(self.set_hy_type)
        
        self.hy_x_label = QtWidgets.QLabel("X (m)")
        self.hy_x_edit = QtWidgets.QLineEdit()
        
        self.hy_y_label = QtWidgets.QLabel("Y (m)")
        self.hy_y_edit = QtWidgets.QLineEdit()
        
        self.hy_azm_label = QtWidgets.QLabel("Azimtuh (deg)")
        self.hy_azm_edit = QtWidgets.QLineEdit()
        
        self.hy_acqchn_label = QtWidgets.QLabel("Acq. Channel")
        self.hy_acqchn_combo = QtWidgets.QComboBox()
        self.hy_acqchn_combo.addItems(self._chn_list)
        
        ##--> hz
        self.hz_label = QtWidgets.QLabel('Magnetometer HZ')
        self.hz_label.setFont(label_font)
        
        self.hz_id_label = QtWidgets.QLabel('ID')
        self.hz_id_edit = QtWidgets.QLineEdit('{0}'.format(self.FieldNotes.magnetometer_hz.id))
        self.hz_id_edit.editingFinished.connect(self.set_hz_id)
        
        self.hz_man_label = QtWidgets.QLabel('Manufacturer')
        self.hz_man_edit = QtWidgets.QLineEdit(self.FieldNotes.magnetometer_hz.manufacturer)
        self.hz_man_edit.editingFinished.connect(self.set_hz_man)
        
        self.hz_type_label = QtWidgets.QLabel('Type')
        self.hz_type_edit = QtWidgets.QLineEdit(self.FieldNotes.magnetometer_hz.type)
        self.hz_type_edit.editingFinished.connect(self.set_hz_type)
        
        self.hz_x_label = QtWidgets.QLabel("X (m)")
        self.hz_x_edit = QtWidgets.QLineEdit()
        
        self.hz_y_label = QtWidgets.QLabel("Y (m)")
        self.hz_y_edit = QtWidgets.QLineEdit()
        
        self.hz_azm_label = QtWidgets.QLabel("Azimtuh (deg)")
        self.hz_azm_edit = QtWidgets.QLineEdit()
        
        self.hz_acqchn_label = QtWidgets.QLabel("Acq. Channel")
        self.hz_acqchn_combo = QtWidgets.QComboBox()
        self.hz_acqchn_combo.addItems(self._chn_list)
        
        ##--> data quality
        self.dq_label = QtWidgets.QLabel('Data Quality')
        self.dq_label.setFont(label_font)
        self.dq_good_periods = QtWidgets.QLabel('Good Periods (min, max)')
        self.dq_good_periods_min = QtWidgets.QLineEdit()
        self.dq_good_periods_max = QtWidgets.QLineEdit()
        
        self.dq_rating_label = QtWidgets.QLabel('Rating')
        self.dq_rating_combo = QtWidgets.QComboBox()
        self.dq_rating_combo.addItems(self._rating_list)
        
        self.dq_warning_flag_label = QtWidgets.QLabel('Warning Flag')
        self.dq_warning_flag_combo = QtWidgets.QComboBox()
        self.dq_warning_flag_combo.addItems(['False', 'True'])
        
        self.dq_warning_comments_label = QtWidgets.QLabel('Warning Comments')
        self.dq_warning_comments_edit = QtWidgets.QLineEdit()
        
        self.dq_comments = QtWidgets.QTextEdit()
        self.dq_comments.setText('Data Quaility Comments')
        
        
        #-->  layout
        layout = QtWidgets.QFormLayout()
        layout.addRow(self.data_logger_label, None)
        layout.addRow(self.data_logger_id_label, self.data_logger_id_edit)
        layout.addRow(self.data_logger_man_label, self.data_logger_man_edit)
        layout.addRow(self.data_logger_type_label, self.data_logger_type_edit)
        
        dq_layout = QtWidgets.QGridLayout()
        dq_layout.addWidget(self.dq_label, 0, 0)
        dq_layout.addWidget(self.dq_good_periods, 1, 0)
        dq_layout.addWidget(self.dq_good_periods_min, 1, 1)
        dq_layout.addWidget(self.dq_good_periods_max, 1, 2)
        dq_layout.addWidget(self.dq_warning_flag_label, 2, 0)
        dq_layout.addWidget(self.dq_warning_flag_combo, 2, 1)
        dq_layout.addWidget(self.dq_warning_comments_label, 2, 2)
        dq_layout.addWidget(self.dq_warning_comments_edit, 2, 3)
        dq_layout.addWidget(self.dq_comments, 3, 0, 2, 4)
        
        ex_layout = QtWidgets.QGridLayout()
        ex_layout.addWidget(self.ex_label, 0, 0)
        ex_layout.addWidget(self.ex_id_label, 1, 0)
        ex_layout.addWidget(self.ex_id_edit, 1, 1)
        ex_layout.addWidget(self.ex_man_label, 1, 2)
        ex_layout.addWidget(self.ex_man_edit, 1, 3)
        ex_layout.addWidget(self.ex_type_label, 1, 4)
        ex_layout.addWidget(self.ex_type_edit, 1, 5)
        ex_layout.addWidget(self.ex_x_label, 2, 0)
        ex_layout.addWidget(self.ex_x_edit, 2, 1)
        ex_layout.addWidget(self.ex_y_label, 2, 2)
        ex_layout.addWidget(self.ex_y_edit, 2, 3)
        ex_layout.addWidget(self.ex_x2_label, 2, 4)
        ex_layout.addWidget(self.ex_x2_edit, 2, 5)
        ex_layout.addWidget(self.ex_y2_label, 2, 6)
        ex_layout.addWidget(self.ex_y2_edit, 2, 7)
        ex_layout.addWidget(self.ex_acqchn_label, 1, 6)
        ex_layout.addWidget(self.ex_acqchn_combo, 1, 7)
        
        ey_layout = QtWidgets.QGridLayout()
        ey_layout.addWidget(self.ey_label, 0, 0)
        ey_layout.addWidget(self.ey_id_label, 1, 0)
        ey_layout.addWidget(self.ey_id_edit, 1, 1)
        ey_layout.addWidget(self.ey_man_label, 1, 2)
        ey_layout.addWidget(self.ey_man_edit, 1, 3)
        ey_layout.addWidget(self.ey_type_label, 1, 4)
        ey_layout.addWidget(self.ey_type_edit, 1, 5)
        ey_layout.addWidget(self.ey_x_label, 2, 0)
        ey_layout.addWidget(self.ey_x_edit, 2, 1)
        ey_layout.addWidget(self.ey_y_label, 2, 2)
        ey_layout.addWidget(self.ey_y_edit, 2, 3)
        ey_layout.addWidget(self.ey_x2_label, 2, 4)
        ey_layout.addWidget(self.ey_x2_edit, 2, 5)
        ey_layout.addWidget(self.ey_y2_label, 2, 6)
        ey_layout.addWidget(self.ey_y2_edit, 2, 7)
        ey_layout.addWidget(self.ey_acqchn_label, 1, 6)
        ey_layout.addWidget(self.ey_acqchn_combo, 1, 7)
        
        hx_layout = QtWidgets.QGridLayout()
        hx_layout.addWidget(self.hx_label, 0, 0)
        hx_layout.addWidget(self.hx_id_label, 1, 0)
        hx_layout.addWidget(self.hx_id_edit, 1, 1)
        hx_layout.addWidget(self.hx_man_label, 1, 2)
        hx_layout.addWidget(self.hx_man_edit, 1, 3)
        hx_layout.addWidget(self.hx_type_label, 1, 4)
        hx_layout.addWidget(self.hx_type_edit, 1, 5)
        hx_layout.addWidget(self.hx_x_label, 2, 0)
        hx_layout.addWidget(self.hx_x_edit, 2, 1)
        hx_layout.addWidget(self.hx_y_label, 2, 2)
        hx_layout.addWidget(self.hx_y_edit, 2, 3)
        hx_layout.addWidget(self.hx_azm_label, 2, 4)
        hx_layout.addWidget(self.hx_azm_edit, 2, 5)
        hx_layout.addWidget(self.hx_acqchn_label, 1, 6)
        hx_layout.addWidget(self.hx_acqchn_combo, 1, 7)
        
        hy_layout = QtWidgets.QGridLayout()
        hy_layout.addWidget(self.hy_label, 0, 0)
        hy_layout.addWidget(self.hy_id_label, 1, 0)
        hy_layout.addWidget(self.hy_id_edit, 1, 1)
        hy_layout.addWidget(self.hy_man_label, 1, 2)
        hy_layout.addWidget(self.hy_man_edit, 1, 3)
        hy_layout.addWidget(self.hy_type_label, 1, 4)
        hy_layout.addWidget(self.hy_type_edit, 1, 5)
        hy_layout.addWidget(self.hy_x_label, 2, 0)
        hy_layout.addWidget(self.hy_x_edit, 2, 1)
        hy_layout.addWidget(self.hy_y_label, 2, 2)
        hy_layout.addWidget(self.hy_y_edit, 2, 3)
        hy_layout.addWidget(self.hy_azm_label, 2, 4)
        hy_layout.addWidget(self.hy_azm_edit, 2, 5)
        hy_layout.addWidget(self.hy_acqchn_label, 1, 6)
        hy_layout.addWidget(self.hy_acqchn_combo, 1, 7)
        
        hz_layout = QtWidgets.QGridLayout()
        hz_layout.addWidget(self.hz_label, 0, 0)
        hz_layout.addWidget(self.hz_id_label, 1, 0)
        hz_layout.addWidget(self.hz_id_edit, 1, 1)
        hz_layout.addWidget(self.hz_man_label, 1, 2)
        hz_layout.addWidget(self.hz_man_edit, 1, 3)
        hz_layout.addWidget(self.hz_type_label, 1, 4)
        hz_layout.addWidget(self.hz_type_edit, 1, 5)
        hz_layout.addWidget(self.hz_x_label, 2, 0)
        hz_layout.addWidget(self.hz_x_edit, 2, 1)
        hz_layout.addWidget(self.hz_y_label, 2, 2)
        hz_layout.addWidget(self.hz_y_edit, 2, 3)
        hz_layout.addWidget(self.hz_azm_label, 2, 4)
        hz_layout.addWidget(self.hz_azm_edit, 2, 5)
        hz_layout.addWidget(self.hz_acqchn_label, 1, 6)
        hz_layout.addWidget(self.hz_acqchn_combo, 1, 7)
        
        final_layout = QtWidgets.QVBoxLayout()
        final_layout.addLayout(layout)
        final_layout.addLayout(dq_layout)
        final_layout.addLayout(ex_layout)
        final_layout.addLayout(ey_layout)
        final_layout.addLayout(hx_layout)
        final_layout.addLayout(hy_layout)
        final_layout.addLayout(hz_layout)
    
        self.setLayout(final_layout)
        
        
    def set_dl_id(self):
        self.FieldNotes.data_logger.id = self.data_logger_id_edit.text()
  
    def set_dl_man(self):
        self.FieldNotes.data_logger.manufacturer = self.data_logger_man_edit.text()
        
    def set_dl_type(self):
        self.FieldNotes.data_logger.type = self.data_logger_type_edit.text()
        
    def set_ex_id(self):
        self.FieldNotes.electrode_ex.id = self.ex_id_edit.text()
  
    def set_ex_man(self):
        self.FieldNotes.electrode_ex.manufacturer = self.ex_man_edit.text()
        
    def set_ex_type(self):
        self.FieldNotes.electrode_ex.type = self.ex_type_edit.text()
        
    def set_ey_id(self):
        self.FieldNotes.electrode_ey.id = self.ey_id_edit.text()
  
    def set_ey_man(self):
        self.FieldNotes.electrode_ey.manufacturer = self.ey_man_edit.text()
        
    def set_ey_type(self):
        self.FieldNotes.electrode_ey.type = self.ey_type_edit.text()
        
    def set_hx_id(self):
        self.FieldNotes.magnetometer_hx.id = self.hx_id_edit.text()
  
    def set_hx_man(self):
        self.FieldNotes.magnetometer_hx.manufacturer = self.hx_man_edit.text()
        
    def set_hx_type(self):
        self.FieldNotes.magnetometer_hx.type = self.hx_type_edit.text()
        
    def set_hy_id(self):
        self.FieldNotes.magnetometer_hy.id = self.hy_id_edit.text()
  
    def set_hy_man(self):
        self.FieldNotes.magnetometer_hy.manufacturer = self.hy_man_edit.text()
        
    def set_hy_type(self):
        self.FieldNotes.magnetometer_hy.type = self.hy_type_edit.text()
        
    def set_hz_id(self):
        self.FieldNotes.magnetometer_hz.id = self.hz_id_edit.text()
  
    def set_hz_man(self):
        self.FieldNotes.magnetometer_hz.manufacturer = self.hz_man_edit.text()
        
    def set_hz_type(self):
        self.FieldNotes.magnetometer_hz.type = self.hz_type_edit.text()
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = MTMainWindow()
    ex.show()
    sys.exit(app.exec_())