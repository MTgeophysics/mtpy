# -*- coding: utf-8 -*-
"""
Created on Wed May 03 19:01:42 2017

@author: jrpeacock
"""

import sys
from PyQt5 import QtWidgets, QtGui, QtCore
import mtpy.core.mt as mt

# header label font
label_font = QtGui.QFont()
label_font.setBold = True
label_font.setPointSize (13)


 
class MTMainWindow(QtWidgets.QMainWindow):
 
    def __init__(self):
        super(MTMainWindow, self).__init__()
        self.setWindowTitle('MT File Editor')
        self.mt_obj = mt.MT()
 
        self.setup_ui()
        
    def setup_ui(self):
        """
        setup user interface
        """
        screen_size = QtWidgets.QDesktopWidget().availableGeometry()
        width = screen_size.width()
        
        self.menu_file = self.menuBar().addMenu("&File")
        
        self.action_open_file = self.menu_file.addAction("&Open")
        self.action_open_file.triggered.connect(self.get_mt_file)
        
        self.action_save_file = self.menu_file.addAction("&Save")
        self.action_save_file.triggered.connect(self.save_mt_file)
        
        self.setWindowState(QtCore.Qt.WindowMaximized)
        
        self.central_widget = QtWidgets.QWidget()
        self.setCentralWidget(self.central_widget)
    
        self.tab_widget = MTTabWidget(self, self.mt_obj)
        
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
        
    def get_mt_file(self):
        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose MT file',
                                           directory=self.plot_widget.dir_path,
                                           filter='*.edi;*.xml;*.j')[0])
        self.mt_obj = mt.MT()
        self.mt_obj.read_mt_file(fn)
        self.tab_widget.mt_obj = self.mt_obj
        
 
class MTTabWidget(QtWidgets.QTabWidget):        
 
    def __init__(self, mt_obj, parent=None):   
        super(MTTabWidget, self).__init__(parent)
        self.mt_obj = mt_obj
        
        self.setup_ui()
        
    def setup_ui(self):
 
        self.tab_site = SiteTab(self)
        self.tab_field = FieldNotesTab(self)
        self.tab_processing = ProcessingTab(self)
        self.tab_provenance = ProvenanceTab(self)
        self.tab_copyright = CopyrightTab(self)
        self.tab_data = DataTab(self)
        
        self.addTab(self.tab_site, "Site")
        self.addTab(self.tab_field, "Field Notes")
        self.addTab(self.tab_processing, "Processing")
        self.addTab(self.tab_provenance, "Provenance")
        self.addTab(self.tab_copyright, "Copyright")
        self.addTab(self.tab_data, "Data")

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
            print(self.Site.Location.easting, self.Site.Location.northing, self.Site.Location.utm_zone)
            
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
        self.data_logger_id_edit = QtWidgets.QLineEdit(self.FieldNotes.DataLogger.id)
        self.data_logger_id_edit.editingFinished.connect(self.set_dl_id)
        
        self.data_logger_man_label = QtWidgets.QLabel('Manufacturer')
        self.data_logger_man_edit = QtWidgets.QLineEdit(self.FieldNotes.DataLogger.manufacturer)
        self.data_logger_man_edit.editingFinished.connect(self.set_dl_man)
        
        self.data_logger_type_label = QtWidgets.QLabel('Type')
        self.data_logger_type_edit = QtWidgets.QLineEdit(self.FieldNotes.DataLogger.type)
        self.data_logger_type_edit.editingFinished.connect(self.set_dl_type)
        
        
        #--> Instrument information
        self.ex_widget = Electrode_Widget(self.FieldNotes.Electrode_ex,
                                          comp='EX')
        self.ey_widget = Electrode_Widget(self.FieldNotes.Electrode_ey,
                                          comp='EY')
    
        self.hx_widget = Magnetometer_Widget(self.FieldNotes.Magnetometer_hx,
                                             comp='HX')
        self.hy_widget = Magnetometer_Widget(self.FieldNotes.Magnetometer_hy,
                                             comp='HY')
        self.hz_widget = Magnetometer_Widget(self.FieldNotes.Magnetometer_hz,
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
        self.FieldNotes.DataLogger.id = self.data_logger_id_edit.text()
  
    def set_dl_man(self):
        self.FieldNotes.DataLogger.manufacturer = self.data_logger_man_edit.text()
        
    def set_dl_type(self):
        self.FieldNotes.DataLogger.type = self.data_logger_type_edit.text()
        
    def set_dq_period_min(self):
        self.FieldNotes.DataQuality.good_from_period = _check_float(self.dq_good_periods_min.text())
        self.dq_good_periods_min.setText('{0:.5g}'.format(self.FieldNotes.DataQuality.good_from_period))
        
    def set_dq_period_max(self):
        self.FieldNotes.DataQuality.good_to_period = _check_float(self.dq_good_periods_max.text())
        self.dq_good_periods_max.setText('{0:.5g}'.format(self.FieldNotes.DataQuality.good_to_period))
        
    def set_dq_rating(self):
        self.FieldNotes.DataQuality.rating = self.dq_rating_combo.currentIndex()
        
    def set_dq_flag(self):
        self.FieldNotes.DataQuality.warnings_flag = self.dq_warning_flag_combo.currentIndex()
        
    def set_dq_warning_comments(self):
        self.FieldNotes.DataQuality.warnings_comments = self.dq_warning_comments_edit.text()
        
    def set_dq_comments(self):
        self.FieldNotes.DataQuality.comments = self.dq_comments.toPlainText()
        
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
        
#        self.parameters = []
#        self.add_parameter_button = QtWidgets.QPushButton('Add Parameter')
#        self.add_parameter_button.pressed.connect(self.add_parameter)
        
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
# Provenance
#==============================================================================
class ProvenanceTab(QtWidgets.QWidget):
    """
    Provenance
    """       
    
    def __init__(self, parent=None):
        super(ProvenanceTab, self).__init__(parent)
        
        self.Provenance = mt.Provenance()
        
        self.setup_ui()
        
    def setup_ui(self):
        
        self.creating_app_label = QtWidgets.QLabel('Creating Application')
        self.creating_app_edit = QtWidgets.QLineEdit()
        self.creating_app_edit.editingFinished.connect(self.set_creating_app)
        
        self.creation_time_label = QtWidgets.QLabel('Creation Date')
        self.creation_time_edit = QtWidgets.QDateEdit()
        self.creation_time_edit.setCalendarPopup(True)
        self.creation_time_edit.setDisplayFormat('yyyy-MM-dd')
        self.creation_time_edit.dateChanged.connect(self.set_creation_time)
        
        self.creator_label = QtWidgets.QLabel('Creator')
        self.creator_label.setFont(label_font)
        
        self.creator_name_label = QtWidgets.QLabel('Name')
        self.creator_name_edit = QtWidgets.QLineEdit()
        self.creator_name_edit.editingFinished.connect(self.set_creator_name)

        self.creator_email_label = QtWidgets.QLabel('email')
        self.creator_email_edit = QtWidgets.QLineEdit()
        self.creator_email_edit.editingFinished.connect(self.set_creator_email)
        
        self.creator_org_label = QtWidgets.QLabel('Organization')
        self.creator_org_edit = QtWidgets.QLineEdit()
        self.creator_org_edit.editingFinished.connect(self.set_creator_org)
        
        self.creator_url_label = QtWidgets.QLabel('Organization URL')
        self.creator_url_edit = QtWidgets.QLineEdit()
        self.creator_url_edit.editingFinished.connect(self.set_creator_url)
        
        self.submitter_label = QtWidgets.QLabel('Submitter')
        self.submitter_label.setFont(label_font)
        
        self.submitter_name_label = QtWidgets.QLabel('Name')
        self.submitter_name_edit = QtWidgets.QLineEdit()
        self.submitter_name_edit.editingFinished.connect(self.set_submitter_name)

        self.submitter_email_label = QtWidgets.QLabel('email')
        self.submitter_email_edit = QtWidgets.QLineEdit()
        self.submitter_email_edit.editingFinished.connect(self.set_submitter_email)
        
        self.submitter_org_label = QtWidgets.QLabel('Organization')
        self.submitter_org_edit = QtWidgets.QLineEdit()
        self.submitter_org_edit.editingFinished.connect(self.set_submitter_org)
        
        self.submitter_url_label = QtWidgets.QLabel('Organization URL')
        self.submitter_url_edit = QtWidgets.QLineEdit()
        self.submitter_url_edit.editingFinished.connect(self.set_submitter_url)
        
        ##--> Layout
        creation_layout = QtWidgets.QFormLayout()
        
        creation_layout.addRow(self.creating_app_label, 
                               self.creating_app_edit)
        creation_layout.addRow(self.creation_time_label,
                               self.creation_time_edit)
        creation_layout.setAlignment(QtCore.Qt.AlignTop)
        
        creator_layout = QtWidgets.QGridLayout()
        creator_layout.addWidget(self.creator_label, 0, 0)
        creator_layout.addWidget(self.creator_name_label, 1, 0)
        creator_layout.addWidget(self.creator_name_edit, 1, 1)
        creator_layout.addWidget(self.creator_email_label, 1, 2)
        creator_layout.addWidget(self.creator_email_edit, 1, 3)
        creator_layout.addWidget(self.creator_org_label, 1, 4)
        creator_layout.addWidget(self.creator_org_edit, 1, 5)
        creator_layout.addWidget(self.creator_url_label, 2, 0)
        creator_layout.addWidget(self.creator_url_edit, 2, 1, 1, 5)
        creator_layout.setAlignment(QtCore.Qt.AlignTop)
        
        submitter_layout = QtWidgets.QGridLayout()
        submitter_layout.addWidget(self.submitter_label, 0, 0)
        submitter_layout.addWidget(self.submitter_name_label, 1, 0)
        submitter_layout.addWidget(self.submitter_name_edit, 1, 1)
        submitter_layout.addWidget(self.submitter_email_label, 1, 2)
        submitter_layout.addWidget(self.submitter_email_edit, 1, 3)
        submitter_layout.addWidget(self.submitter_org_label, 1, 4)
        submitter_layout.addWidget(self.submitter_org_edit, 1, 5)
        submitter_layout.addWidget(self.submitter_url_label, 2, 0)
        submitter_layout.addWidget(self.submitter_url_edit, 2, 1, 1, 5)
        submitter_layout.setAlignment(QtCore.Qt.AlignTop)
        
        final_layout = QtWidgets.QVBoxLayout()
        
        final_layout.addLayout(creation_layout)
        final_layout.addLayout(creator_layout)
        final_layout.addLayout(submitter_layout)
        final_layout.addStretch(0)
        final_layout.setAlignment(QtCore.Qt.AlignTop)

        
        self.setLayout(final_layout)
        
    def set_creating_app(self):
        pass
    
    def set_creation_time(self):
        date = self.creation_time_edit.date()
        print(date.toPyDate())

    def set_creator_name(self):
        pass
    
    def set_creator_email(self):
        pass
    
    def set_creator_org(self):
        pass
    
    def set_creator_url(self):
        pass
    
    def set_submitter_name(self):
        pass
    
    def set_submitter_email(self):
        pass
    
    def set_submitter_org(self):
        pass
    
    def set_submitter_url(self):
        pass
    
#==============================================================================
# Copyright
#==============================================================================
class CopyrightTab(QtWidgets.QWidget):
    """
    copyright 
    """

    def __init__(self, parent=None):
        super(CopyrightTab, self).__init__(parent)
        
        self.Copyright = mt.Copyright()
        
        self._release_list = ['Unrestricted Release', 
                              'Academic Use Only',
                              'Restrictions Apply']
        
        self.setup_ui()
        
    def setup_ui(self):
        
        self.citation_label = QtWidgets.QLabel('Citation')
        self.citation_label.setFont(label_font)
        
        self.citation_author_label = QtWidgets.QLabel('Author')
        self.citation_author_edit = QtWidgets.QLineEdit()
        self.citation_author_edit.editingFinished.connect(self.set_author)
        
        self.citation_title_label = QtWidgets.QLabel('Title')
        self.citation_title_edit = QtWidgets.QLineEdit()
        self.citation_title_edit.editingFinished.connect(self.set_title)
        
        self.citation_journal_label = QtWidgets.QLabel('Journal')
        self.citation_journal_edit = QtWidgets.QLineEdit()
        self.citation_journal_edit.editingFinished.connect(self.set_journal)

        self.citation_volume_label = QtWidgets.QLabel('Volume')
        self.citation_volume_edit = QtWidgets.QLineEdit()
        self.citation_volume_edit.editingFinished.connect(self.set_volume)

        self.citation_year_label = QtWidgets.QLabel('Year')
        self.citation_year_edit = QtWidgets.QLineEdit()
        self.citation_year_edit.editingFinished.connect(self.set_year)
        
        self.citation_doi_label = QtWidgets.QLabel('DOI')
        self.citation_doi_edit = QtWidgets.QLineEdit()
        self.citation_doi_edit.editingFinished.connect(self.set_doi)
        
        self.release_status_name = QtWidgets.QLabel('Release Status')
        self.release_status_combo = QtWidgets.QComboBox()
        self.release_status_combo.addItems(self._release_list)
        self.release_status_combo.currentIndexChanged.connect(self.set_release_status)
        
        self.conditions_of_use_label = QtWidgets.QLabel('Conditions of Use')
        self.conditions_of_use_edit = QtWidgets.QTextEdit()
        self.conditions_of_use_edit.setText(self.Copyright.conditions_of_use)
        self.conditions_of_use_edit.textChanged.connect(self.set_conditions)
        
        ##--> layout
        cite_layout = QtWidgets.QGridLayout()
        cite_layout.addWidget(self.citation_label, 0, 0)
        cite_layout.addWidget(self.citation_author_label, 1, 0, 1, 5)
        cite_layout.addWidget(self.citation_author_edit, 1, 1, 1, 5)
        cite_layout.addWidget(self.citation_title_label, 2, 0, 1, 5)
        cite_layout.addWidget(self.citation_title_edit, 2, 1, 1, 5)
        cite_layout.addWidget(self.citation_journal_label, 3, 0)
        cite_layout.addWidget(self.citation_journal_edit, 3, 1)
        cite_layout.addWidget(self.citation_volume_label, 3, 2)
        cite_layout.addWidget(self.citation_volume_edit, 3, 3)
        cite_layout.addWidget(self.citation_year_label, 3, 4)
        cite_layout.addWidget(self.citation_year_edit, 3, 5)
        cite_layout.addWidget(self.citation_doi_label, 4, 0, 1, 5)
        cite_layout.addWidget(self.citation_doi_edit, 4, 1, 1, 5)
        cite_layout.setAlignment(QtCore.Qt.AlignTop)
        
        combo_layout = QtWidgets.QHBoxLayout()
        combo_layout.addWidget(self.release_status_name)
        combo_layout.addWidget(self.release_status_combo)
        
        release_layout = QtWidgets.QVBoxLayout()
        release_layout.addLayout(combo_layout)
        release_layout.addWidget(self.conditions_of_use_label)
        release_layout.addWidget(self.conditions_of_use_edit)
        
        final_layout = QtWidgets.QVBoxLayout()
        final_layout.addLayout(cite_layout)
        final_layout.addLayout(release_layout)
        
        self.setLayout(final_layout)
        
    def set_author(self):
        pass

    def set_title(self):
        pass

    def set_journal(self):
        pass

    def set_volume(self):
        pass

    def set_year(self):
        pass

    def set_doi(self):
        pass

    def set_release_status(self):
        pass
    
    def set_conditions(self):
        pass
    
#==============================================================================
# Data
#==============================================================================
class DataTab(QtWidgets.QWidget):
    """
    hold the data in tabular form
    """
    
    def __init__(self, parent=None):
        super(DataTab, self).__init__(parent)
        
        self.Data = None
        
        self._z_headers = ['Frequency (Hz)', 
                           'Real Zxx', 'Imag Zxx', 'Err Zxx',
                           'Real Zxy', 'Imag Zxy', 'Err Zxy',
                           'Real Zyx', 'Imag Zyx', 'Err Zyx',
                           'Real Zyy', 'Imag Zyy', 'Err Zyy']
        self._t_headers = ['Frequency (Hz)', 
                           'Real Tzx', 'Imag Tzx', 'Err Tzx',
                           'Real Tzy', 'Imag Tzy', 'Err Tzy']
        
        self.setup_ui()
        
    def setup_ui(self):
        
        self.tab = QtWidgets.QTabWidget()
        
        self.data_z_table = QtWidgets.QTableWidget()
        self.data_z_table.setColumnCount(13)
        self.data_z_table.setRowCount(100)
        self.data_z_table.setHorizontalHeaderLabels(self._z_headers)
        #setHorizontalHeaderLabels(headerlist)
        self.tab.addTab(self.data_z_table, 'Impedance')
        
        self.data_t_table = QtWidgets.QTableWidget()
        self.data_t_table.setColumnCount(7)
        self.data_t_table.setRowCount(100)
        self.data_t_table.setHorizontalHeaderLabels(self._t_headers)
        self.tab.addTab(self.data_t_table, 'Tipper')
        
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.tab)
        
        self.setLayout(layout)
        
        

        

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