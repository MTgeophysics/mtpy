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
        self.tab_processing = ProcessingTab(self)
        
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
        
        
        #--> Instrument information
        self.ex_widget = Electrode_Widget(self.FieldNotes.electrode_ex,
                                          comp='EX')
        self.ey_widget = Electrode_Widget(self.FieldNotes.electrode_ey,
                                          comp='EY')
    
        self.hx_widget = Magnetometer_Widget(self.FieldNotes.magnetometer_hx,
                                             comp='HX')
        self.hy_widget = Magnetometer_Widget(self.FieldNotes.magnetometer_hy,
                                             comp='HY')
        self.hz_widget = Magnetometer_Widget(self.FieldNotes.magnetometer_hz,
                                             comp='HZ')


        ##--> data quality
        self.dq_label = QtWidgets.QLabel('Data Quality')
        self.dq_label.setFont(label_font)
        self.dq_good_periods = QtWidgets.QLabel('Good Periods (min, max)')
        self.dq_good_periods_min = QtWidgets.QLineEdit()
        self.dq_good_periods_min.editingFinished.connect(self.set_dq_period_min)
        self.dq_good_periods_max = QtWidgets.QLineEdit()
        self.dq_good_periods_max.editingFinished.connect(self.set_dq_period_max)
        
        self.dq_rating_label = QtWidgets.QLabel('Rating')
        self.dq_rating_combo = QtWidgets.QComboBox()
        self.dq_rating_combo.addItems(self._rating_list)
        self.dq_rating_combo.currentIndexChanged.connect(self.set_dq_rating)
        
        self.dq_warning_flag_label = QtWidgets.QLabel('Warning Flag')
        self.dq_warning_flag_combo = QtWidgets.QComboBox()
        self.dq_warning_flag_combo.addItems(['False', 'True'])
        self.dq_warning_flag_combo.currentIndexChanged.connect(self.set_dq_flag)
        
        self.dq_warning_comments_label = QtWidgets.QLabel('Warning Comments')
        self.dq_warning_comments_edit = QtWidgets.QLineEdit()
        self.dq_warning_comments_edit.editingFinished.connect(self.set_dq_warning_comments)
        
        self.dq_comments = QtWidgets.QTextEdit()
        self.dq_comments.setText('Data Quaility Comments')
        self.dq_comments.textChanged.connect(self.set_dq_comments)
        
        
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
        
        final_layout = QtWidgets.QVBoxLayout()
        final_layout.addLayout(layout)
        final_layout.addLayout(dq_layout)
        final_layout.addWidget(self.ex_widget)
        final_layout.addWidget(self.ey_widget)
        final_layout.addWidget(self.hx_widget)
        final_layout.addWidget(self.hy_widget)
        final_layout.addWidget(self.hz_widget)
    
        self.setLayout(final_layout)
        
        
    def set_dl_id(self):
        self.FieldNotes.data_logger.id = self.data_logger_id_edit.text()
  
    def set_dl_man(self):
        self.FieldNotes.data_logger.manufacturer = self.data_logger_man_edit.text()
        
    def set_dl_type(self):
        self.FieldNotes.data_logger.type = self.data_logger_type_edit.text()
        
    def set_dq_period_min(self):
        self.FieldNotes.data_quality.good_from_period = _check_float(self.dq_good_periods_min.text())
        self.dq_good_periods_min.setText('{0:.5g}'.format(self.FieldNotes.data_quality.good_from_period))
        
    def set_dq_period_max(self):
        self.FieldNotes.data_quality.good_to_period = _check_float(self.dq_good_periods_max.text())
        self.dq_good_periods_max.setText('{0:.5g}'.format(self.FieldNotes.data_quality.good_to_period))
        
    def set_dq_rating(self):
        self.FieldNotes.data_quality.rating = self.dq_rating_combo.currentIndex()
        
    def set_dq_flag(self):
        self.FieldNotes.data_quality.warnings_flag = self.dq_warning_flag_combo.currentIndex()
        
    def set_dq_warning_comments(self):
        self.FieldNotes.data_quality.warnings_comments = self.dq_warning_comments_edit.text()
        
    def set_dq_comments(self):
        self.FieldNotes.data_quality.comments = self.dq_comments.toPlainText()
        
#==============================================================================
# Electrode        
#==============================================================================  
class Electrode_Widget(QtWidgets.QWidget):
    """
    class to hold Magnetometer information
    """
    
    def __init__(self, electrode_class, comp='EX', parent=None):
        super(Electrode_Widget, self).__init__(parent)

        self.Electrode = electrode_class
        self.comp = comp
        self._chn_list = ['{0:d}'.format(ii) for ii in range(1, 7, 1)]
        self.setup_ui()
        
    def setup_ui(self):
        
        self.e_label = QtWidgets.QLabel('Electrode {0}'.format(self.comp))
        self.e_label.setFont(label_font)
        
        self.e_id_label = QtWidgets.QLabel('ID')
        self.e_id_edit = QtWidgets.QLineEdit('{0}'.format(self.Electrode.id))
        self.e_id_edit.editingFinished.connect(self.set_e_id)
        
        self.e_man_label = QtWidgets.QLabel('Manufacturer')
        self.e_man_edit = QtWidgets.QLineEdit(self.Electrode.manufacturer)
        self.e_man_edit.editingFinished.connect(self.set_e_man)
        
        self.e_type_label = QtWidgets.QLabel('Type')
        self.e_type_edit = QtWidgets.QLineEdit(self.Electrode.type)
        self.e_type_edit.editingFinished.connect(self.set_e_type)

        self.e_x_label = QtWidgets.QLabel("X (m)")
        self.e_x_edit = QtWidgets.QLineEdit()
        self.e_x_edit.editingFinished.connect(self.set_x)
        
        self.e_y_label = QtWidgets.QLabel("Y (m)")
        self.e_y_edit = QtWidgets.QLineEdit()
        self.e_y_edit.editingFinished.connect(self.set_y)
        
        self.e_x2_label = QtWidgets.QLabel("X2 (m)")
        self.e_x2_edit = QtWidgets.QLineEdit()
        self.e_x2_edit.editingFinished.connect(self.set_x2)
        
        self.e_y2_label = QtWidgets.QLabel("Y2 (m)")
        self.e_y2_edit = QtWidgets.QLineEdit()
        self.e_y2_edit.editingFinished.connect(self.set_y2)
        
        self.e_acqchn_label = QtWidgets.QLabel("Acq. Channel")
        self.e_acqchn_combo = QtWidgets.QComboBox()
        self.e_acqchn_combo.addItems(self._chn_list)
        self.e_acqchn_combo.currentIndexChanged.connect(self.set_chn)

        ##--> set layout
        e_layout = QtWidgets.QGridLayout()
        e_layout.addWidget(self.e_label, 0, 0)
        e_layout.addWidget(self.e_id_label, 1, 0)
        e_layout.addWidget(self.e_id_edit, 1, 1)
        e_layout.addWidget(self.e_man_label, 1, 2)
        e_layout.addWidget(self.e_man_edit, 1, 3)
        e_layout.addWidget(self.e_type_label, 1, 4)
        e_layout.addWidget(self.e_type_edit, 1, 5)
        e_layout.addWidget(self.e_x_label, 2, 0)
        e_layout.addWidget(self.e_x_edit, 2, 1)
        e_layout.addWidget(self.e_y_label, 2, 2)
        e_layout.addWidget(self.e_y_edit, 2, 3)
        e_layout.addWidget(self.e_x2_label, 2, 4)
        e_layout.addWidget(self.e_x2_edit, 2, 5)
        e_layout.addWidget(self.e_y2_label, 2, 6)
        e_layout.addWidget(self.e_y2_edit, 2, 7)
        e_layout.addWidget(self.e_acqchn_label, 1, 6)
        e_layout.addWidget(self.e_acqchn_combo, 1, 7)
        
        self.setLayout(e_layout)
        
        
    def set_e_id(self):
        self.Electrode.id = self.e_id_edit.text()
  
    def set_e_man(self):
        self.Electrode.manufacturer = self.e_man_edit.text()
        
    def set_e_type(self):
        self.Electrode.type = self.e_type_edit.text()
        
    def set_x(self):
        self.Electrode.x = _check_float(self.e_x_edit.text())
        self.e_x_edit.setText('{0:.2f}'.format(self.Electrode.x))

    def set_x2(self):
        self.Electrode.x2 = _check_float(self.e_x2_edit.text())
        self.e_x2_edit.setText('{0:.2f}'.format(self.Electrode.x2))    
    
    def set_y(self):
        self.Electrode.y = _check_float(self.e_y_edit.text())
        self.e_y_edit.setText('{0:.2f}'.format(self.Electrode.y))
    
    def set_y2(self):
        self.Electrode.y2 = _check_float(self.e_y2_edit.text())
        self.e_y2_edit.setText('{0:.2f}'.format(self.Electrode.y2))
        
    def set_chn(self):
        self.Electrode.acqchn = int(self.e_acqchn_combo.currentIndex())
        
#==============================================================================
# Magnetometer        
#==============================================================================
class Magnetometer_Widget(QtWidgets.QWidget):
    """
    class to hold magnetometer information
    """
    
    def __init__(self, magnetometer_class, comp='HX', parent=None):
        super(Magnetometer_Widget, self).__init__(parent)

        self.Magnetometer = magnetometer_class
        self.comp = comp
        self._chn_list = ['{0:d}'.format(ii) for ii in range(1, 7, 1)]
        self.setup_ui()
        
    def setup_ui(self):
        
        self.h_label = QtWidgets.QLabel('Magnetometer {0}'.format(self.comp))
        self.h_label.setFont(label_font)
        
        self.h_id_label = QtWidgets.QLabel('ID')
        self.h_id_edit = QtWidgets.QLineEdit('{0}'.format(self.Magnetometer.id))
        self.h_id_edit.editingFinished.connect(self.set_id)
    
        self.h_man_label = QtWidgets.QLabel('Manufacturer')
        self.h_man_edit = QtWidgets.QLineEdit(self.Magnetometer.manufacturer)
        self.h_man_edit.editingFinished.connect(self.set_man)
    
        self.h_type_label = QtWidgets.QLabel('Type')
        self.h_type_edit = QtWidgets.QLineEdit(self.Magnetometer.type)
        self.h_type_edit.editingFinished.connect(self.set_type)

        self.h_x_label = QtWidgets.QLabel("X (m)")
        self.h_x_edit = QtWidgets.QLineEdit()
        self.h_x_edit.editingFinished.connect(self.set_x)
        
        self.h_y_label = QtWidgets.QLabel("Y (m)")
        self.h_y_edit = QtWidgets.QLineEdit()
        self.h_y_edit.editingFinished.connect(self.set_y)
        
        self.h_azm_label = QtWidgets.QLabel("Azimuth (deg)")
        self.h_azm_edit = QtWidgets.QLineEdit()
        self.h_azm_edit.editingFinished.connect(self.set_azm)
        
        self.h_acqchn_label = QtWidgets.QLabel("Acq. Channel")
        self.h_acqchn_combo = QtWidgets.QComboBox()
        self.h_acqchn_combo.addItems(self._chn_list)
        self.h_acqchn_combo.currentIndexChanged.connect(self.set_chn)

        
        ##--> set layout
        h_layout = QtWidgets.QGridLayout()
        h_layout.addWidget(self.h_label, 0, 0)
        h_layout.addWidget(self.h_id_label, 1, 0)
        h_layout.addWidget(self.h_id_edit, 1, 1)
        h_layout.addWidget(self.h_man_label, 1, 2)
        h_layout.addWidget(self.h_man_edit, 1, 3)
        h_layout.addWidget(self.h_type_label, 1, 4)
        h_layout.addWidget(self.h_type_edit, 1, 5)
        h_layout.addWidget(self.h_x_label, 2, 0)
        h_layout.addWidget(self.h_x_edit, 2, 1)
        h_layout.addWidget(self.h_y_label, 2, 2)
        h_layout.addWidget(self.h_y_edit, 2, 3)
        h_layout.addWidget(self.h_azm_label, 2, 4)
        h_layout.addWidget(self.h_azm_edit, 2, 5)
        h_layout.addWidget(self.h_acqchn_label, 1, 6)
        h_layout.addWidget(self.h_acqchn_combo, 1, 7)
        
        self.setLayout(h_layout)
        
        
    def set_id(self):
        self.Magnetometer.id = self.h_id_edit.text()
  
    def set_man(self):
        self.Magnetometer.manufacturer = self.h_man_edit.text()
        
    def set_type(self):
        self.Magnetometer.type = self.h_type_edit.text()
        
    def set_x(self):
        self.Magnetometer.x = _check_float(self.h_x_edit.text())
        self.h_x_edit.setText('{0:.2f}'.format(self.Magnetometer.x))

    def set_y(self):
        self.Magnetometer.y = _check_float(self.h_y_edit.text())
        self.h_y_edit.setText('{0:.2f}'.format(self.Magnetometer.y))    
    
    def set_azm(self):
        self.Magnetometer.azm = _check_float(self.h_azm_edit.text())
        self.h_azm_edit.setText('{0:.2f}'.format(self.Magnetometer.azm))
        
    def set_chn(self):
        self.Magnetometer.acqchn = int(self.h_acqchn_combo.currentIndex())

#==============================================================================
# Processing
#==============================================================================
class ProcessingTab(QtWidgets.QWidget):
    """
    processing tab
    """

    def __init__(self, parent=None):
        super(ProcessingTab, self).__init__(parent)
        
        self.Processing = mt.Processing()
        
        self.setup_ui()
        
    def setup_ui(self):
        
        self.software_label = QtWidgets.QLabel('Software')
        self.software_label.setFont(label_font)
        
        self.software_name_label = QtWidgets.QLabel('Name')
        self.software_name_edit = QtWidgets.QLineEdit()
        self.software_name_edit.editingFinished.connect(self.set_software_name)
        
        self.software_version_label = QtWidgets.QLabel('Version')
        self.software_version_edit = QtWidgets.QLineEdit()
        self.software_version_edit.editingFinished.connect(self.set_software_version)
        
        
        self.software_author_label = QtWidgets.QLabel('Author')
        self.software_author_edit = QtWidgets.QLineEdit()
        self.software_author_edit.editingFinished.connect(self.set_software_author)
        
        self.software_author_email_label = QtWidgets.QLabel('Author Email')
        self.software_author_email_edit = QtWidgets.QLineEdit()
        self.software_author_email_edit.editingFinished.connect(self.set_software_author_email)
        
        self.software_author_org_label = QtWidgets.QLabel('Author Organization')
        self.software_author_org_edit = QtWidgets.QLineEdit()
        self.software_author_org_edit.editingFinished.connect(self.set_software_author_org)
        
        self.software_author_url_label = QtWidgets.QLabel('URL')
        self.software_author_url_edit = QtWidgets.QLineEdit()
        self.software_author_url_edit.editingFinished.connect(self.set_software_author_url)
        
    
        self.software_date_label = QtWidgets.QLabel('Date (YYYY-MM-DD')
        self.software_date_edit = QtWidgets.QLineEdit()
        self.software_date_edit.editingFinished.connect(self.set_software_date)
        
        self.notes_label = QtWidgets.QLabel('Notes:')
        self.notes_label.setFont(label_font)
        
        self.notes_edit = QtWidgets.QTextEdit()
        self.notes_edit.textChanged.connect(self.set_notes)
        
        h_line_00 = QtWidgets.QFrame(self)
        h_line_00.setFrameShape(QtWidgets.QFrame.HLine)
        h_line_00.setFrameShadow(QtWidgets.QFrame.Sunken)
        
        # layout
        grid_layout = QtWidgets.QGridLayout()
        
        grid_layout.addWidget(self.software_label, 0, 0)
        grid_layout.addWidget(self.software_name_label, 1, 0)
        grid_layout.addWidget(self.software_name_edit, 1, 1)
        grid_layout.addWidget(self.software_version_label, 1, 2)
        grid_layout.addWidget(self.software_version_edit, 1, 3)
        grid_layout.addWidget(self.software_date_label, 1, 4)
        grid_layout.addWidget(self.software_date_edit, 1, 5)
        
        grid_layout.addWidget(self.software_author_label, 2, 0)
        grid_layout.addWidget(self.software_author_edit, 2, 1)
        grid_layout.addWidget(self.software_author_email_label, 2, 2)
        grid_layout.addWidget(self.software_author_email_edit, 2, 3)
        grid_layout.addWidget(self.software_author_org_label, 2, 4)
        grid_layout.addWidget(self.software_author_org_edit, 2, 5)
        grid_layout.addWidget(self.software_author_url_label, 3, 0)
        grid_layout.addWidget(self.software_author_url_edit, 3, 1, 1, 5)
        
        notes_layout = QtWidgets.QVBoxLayout()
        notes_layout.addWidget(self.notes_label)
        notes_layout.addWidget(self.notes_edit)
        
        final_layout = QtWidgets.QVBoxLayout()
        final_layout.addLayout(grid_layout)
        final_layout.addWidget(h_line_00)
        final_layout.addLayout(notes_layout)
        
        self.setLayout(final_layout)
        
    def set_software_name(self):
        pass
    
    def set_software_version(self):
        pass
    
    def set_software_author(self):
        pass
    
    def set_software_date(self):
        pass
    
    def set_software_author_email(self):
        pass

    def set_software_author_org(self):
        pass

    def set_software_author_url(self):
        pass
    
    def set_notes(self):
        pass
    
class ProcessingParameter(QtWidgets.QWidget):
    """
    processing name and value
    """    
    
    def __init__(self, parent=None):
        super(ProcessingParameter, self).__init__(parent)
        
        self.name = None
        self.value = None
        
        self.setup_ui()
    
    def setup_ui(self):
        self.value_edit = QtWidgets.QLineEdit()
        self.value_edit.editingFinished.connect(self.set_value)
        
        self.name_edit = QtWidgets.QLineEdit()
        self.name_edit.editingFinished.connect(self.set_name)
        
        # layout
        self.setLayout(QtWidgets.QFormLayout(self.name_edit, self.value_edit))
        
    def set_name(self):
        self.name = self.name_edit.text()
        
    def set_value(self):
        self.value = self.value_edit.text()
        
#==============================================================================
# Common functions
#==============================================================================
def _check_float(value):
    try:
        return_num = float(value)
    except ValueError:
        return_num = 0.0
        
    return return_num

        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = MTMainWindow()
    ex.show()
    sys.exit(app.exec_())