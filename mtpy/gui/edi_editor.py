# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 15:04:27 2016

@author: jpeacock
"""

#==============================================================================
#  Imports
#==============================================================================

from PyQt4 import QtCore, QtGui
import mtpy.core.mt as mt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import copy
import sys
import os
import mtpy.imaging.plotstrike2d as plotstrike2d

#==============================================================================
# UI
#==============================================================================

class MyStream(QtCore.QObject):
    """
    this class will emit a signal to write standard output to the UI
    """
    message = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(MyStream, self).__init__(parent)

    def write(self, message):
        self.message.emit(str(message))

class EDI_Editor_Window(QtGui.QMainWindow):
    """
    This is the main window for editing an edi file.
    
    Includes:
        * Editing points by removing them
        * Adjusting for static shift
        * Removing distortion
        * Rotations
        * Strike estimations
    """
    
    def __init__(self):
        super(EDI_Editor_Window, self).__init__()
        
        
    def ui_setup(self):
        """
        set up the user interface
        """
        
        self.setWindowTitle('EDI Editor')
#        self.resize(1920, 1080)
        self.setWindowState(QtCore.Qt.WindowMaximized)
        
        self.plot_widget = PlotWidget()
        self.centralWidget = self.setCentralWidget(self.plot_widget)
        
        self.menu_file = self.menuBar().addMenu("&File")
        
        self.action_open_file = self.menu_file.addAction("&Open")
        self.action_open_file.triggered.connect(self.get_edi_file)
        
        self.action_close_file = self.menu_file.addAction("C&lose")
        self.action_close_file.triggered.connect(self.close_edi_file)
        
        self.action_save_file = self.menu_file.addAction("&Save")
        self.action_save_file.triggered.connect(self.save_edi_file)
        
        self.menu_plot_properties = self.menuBar().addMenu("Plot Properties")
        self.action_edit_plot = self.menu_plot_properties.addAction("Edit")
        self.action_edit_plot.triggered.connect(self.edit_plot_properties)
        
        self.menu_metadata = self.menuBar().addMenu("Metadata")
        self.action_edit_metadata = self.menu_metadata.addAction("Edit")
        self.action_edit_metadata.triggered.connect(self.edit_metadata)
        
        
        self.my_stream = MyStream()
        self.my_stream.message.connect(self.plot_widget.normal_output)
        
        sys.stdout = self.my_stream
        
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def get_edi_file(self):
        """
        get edi file
        """
        
        print '='*35
        fn_dialog = QtGui.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose EDI file',
                                           directory=self.plot_widget.dir_path,
                                           filter='*.edi'))
                                           
        self.plot_widget.dir_path = os.path.dirname(fn)
                                           
        self.plot_widget.mt_obj = mt.MT(fn)
        self.plot_widget._mt_obj = copy.deepcopy(self.plot_widget.mt_obj)
        self.plot_widget.fill_metadata()
        self.plot_widget.reset_parameters()
        self.plot_widget.plot_properties = PlotSettings(None)
        self.plot_widget.plot()

    def close_edi_file(self):
        pass

    def save_edi_file(self):
        pass
    
    def edit_plot_properties(self):
        self.plot_widget.plot_properties.setup_ui()
        self.plot_widget.plot_properties.show()
        self.plot_widget.plot_properties.settings_updated.connect(self.update_plot)
        
    def update_plot(self):
        
        self.plot_widget.redraw_plot()
    
    def edit_metadata(self):
        pass
                       
#==============================================================================
# Plot Widget     
#==============================================================================
class PlotWidget(QtGui.QWidget):
    """
    matplotlib plot of the data
    
    """
    
    def __init__(self):
        super(PlotWidget, self).__init__()

        self.mt_obj = mt.MT()
        self._mt_obj = mt.MT()
        self.static_shift_x = 1.0
        self.static_shift_y = 1.0
        self.rotate_z_angle = 0
        self.rotate_tip_angle = 0
        self.plot_properties = PlotSettings(None)
        self._edited_ss = False
        self._edited_dist = False
        self._edited_rot = False
        self._edited_mask = False
        self.dir_path = os.getcwd()
        
        self.setup_ui()
        
        
    def setup_ui(self):
        """
        set up the ui
        """
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        self.mpl_widget.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.mpl_widget.setFocus()
        
        # this will set the minimum width of the mpl plot, important to 
        # resize everything
        screen = QtGui.QDesktopWidget().screenGeometry()
        self.mpl_widget.setMinimumWidth(screen.width()*(1600./1920))
        
        # be able to edit the data
        self.mpl_widget.mpl_connect('pick_event', self.on_pick)
        self.mpl_widget.mpl_connect('axes_enter_event', self.in_axes)
        self.mpl_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                      QtGui.QSizePolicy.Expanding)

        self.mpl_toolbar = NavigationToolbar(self.mpl_widget, self)
        
         # header label font
        header_font = QtGui.QFont()
        header_font.setBold = True
        header_font.setPointSize (16)
        
        # output box for all notes from the program
        self.output_label = QtGui.QLabel("Output")
        self.output_label.setFont(header_font)
        self.output_box = QtGui.QTextEdit()
    
        ## --> SET METADATA
        self.metadata_label = QtGui.QLabel('Metadata')
        self.metadata_label.setFont(header_font)
    
        self.meta_station_name_label = QtGui.QLabel("Station")
        self.meta_station_name_edit = QtGui.QLineEdit(str(self.mt_obj.station))
        self.meta_station_name_edit.editingFinished.connect(self.meta_edit_station)
        
        # sets up with an empty mt object so we can set values to 0
        # once a edi is read in then the metadata will be filled.
        self.meta_lat_label = QtGui.QLabel("Lat (deg)")
        self.meta_lat_edit = QtGui.QLineEdit('{0:.6f}'.format(0.0))
        self.meta_lat_edit.editingFinished.connect(self.meta_edit_lat)
        
        self.meta_lon_label = QtGui.QLabel("Long (deg)")
        self.meta_lon_edit = QtGui.QLineEdit('{0:.6f}'.format(0.0))
        self.meta_lon_edit.editingFinished.connect(self.meta_edit_lon)
        
        self.meta_elev_label = QtGui.QLabel("Elev (m)")
        self.meta_elev_edit = QtGui.QLineEdit('{0:.3f}'.format(0.0))
        self.meta_elev_edit.editingFinished.connect(self.meta_edit_elev)
        
        self.meta_date_label = QtGui.QLabel("Date Acq")
        self.meta_date_edit = QtGui.QLineEdit("YYYY-MM-DD")
        self.meta_date_edit.editingFinished.connect(self.meta_edit_date)
        
        self.meta_acq_label = QtGui.QLabel("Acquired By")
        self.meta_acq_edit = QtGui.QLineEdit("None")
        self.meta_acq_edit.editingFinished.connect(self.meta_edit_acq)
        
        ## Static Shift
        self.static_shift_label = QtGui.QLabel("Static Shift")
        self.static_shift_label.setFont(header_font)
        
        self.static_shift_x_label = QtGui.QLabel("X")
        self.static_shift_x_edit = QtGui.QLineEdit("{0:.4g}".format(self.static_shift_x))
        self.static_shift_x_edit.editingFinished.connect(self.static_shift_set_x)
        
        self.static_shift_y_label = QtGui.QLabel("Y")
        self.static_shift_y_edit = QtGui.QLineEdit("{0:.4g}".format(self.static_shift_y))
        self.static_shift_y_edit.editingFinished.connect(self.static_shift_set_y)
        
        self.static_shift_apply_button = QtGui.QPushButton()
        self.static_shift_apply_button.setText("Apply")
        self.static_shift_apply_button.pressed.connect(self.static_shift_apply)
        
        ## remove distortion 
        self.remove_distortion_label = QtGui.QLabel("Remove Distortion")
        self.remove_distortion_label.setFont(header_font)
        
        self.remove_distortion_button = QtGui.QPushButton()
        self.remove_distortion_button.setText("Remove Distortion [Bibby et al., 2005]")
        self.remove_distortion_button.pressed.connect(self.remove_distortion_apply)
        
        ## rotate data
        self.rotate_data_label = QtGui.QLabel("Rotate")
        self.rotate_data_label.setFont(header_font)

        self.rotate_explanation = QtGui.QLabel("Always rotating original data, assuming: N = 0, E = 90.")
        self.rotate_angle_z_label = QtGui.QLabel("Rotate Z (deg)")
        self.rotate_angle_z_edit = QtGui.QLineEdit("{0:.4g}".format(self.rotate_z_angle))
        self.rotate_angle_z_edit.editingFinished.connect(self.rotate_set_z_angle) 
        
        self.rotate_angle_t_label = QtGui.QLabel("Rotate Tipper (deg)")
        self.rotate_angle_t_edit = QtGui.QLineEdit("{0:.4g}".format(self.rotate_tip_angle))
        self.rotate_angle_t_edit.editingFinished.connect(self.rotate_set_t_angle)        
        
        self.rotate_angle_button = QtGui.QPushButton("Apply")
        self.rotate_angle_button.pressed.connect(self.rotate_data_apply)
        
        self.rotate_estimate_strike_button = QtGui.QPushButton("Estimate Strike")
        self.rotate_estimate_strike_button.pressed.connect(self.rotate_estimate_strike)
        
        ## tools label
        self.tools_label = QtGui.QLabel("Tools")
        self.tools_label.setFont(header_font)
       
        ## apply edits button
        self.edits_apply_button = QtGui.QPushButton()
        self.edits_apply_button.setText("Apply Edits")
        self.edits_apply_button.pressed.connect(self.edits_apply)
       
        ## revert back to original data 
        self.revert_button = QtGui.QPushButton()
        self.revert_button.setText("Revert back to orginal data")
        self.revert_button.pressed.connect(self.revert_back)
        
        ## save edits button
        self.save_edits_button = QtGui.QPushButton()
        self.save_edits_button.setText('Save Edits to new EDI file')
        self.save_edits_button.pressed.connect(self.save_edi_file)
        
        ## horizontal line
        h_line = QtGui.QFrame(self)
        h_line.setFrameShape(QtGui.QFrame.HLine)
        h_line.setFrameShadow(QtGui.QFrame.Sunken)
        
        ## vertical spacer
        v_space = QtGui.QSpacerItem(20, 20, 
                                    QtGui.QSizePolicy.Minimum,
                                    QtGui.QSizePolicy.Maximum)
        ###--> layout ---------------------------------------------
        ## mpl plot --> right panel
        mpl_vbox = QtGui.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.mpl_widget)        
        
        ##--> Left Panel
        ## Metadata grid
        meta_layout = QtGui.QGridLayout()
       
        meta_layout.addWidget(self.metadata_label, 0, 0)
        meta_layout.addWidget(self.meta_station_name_label, 1, 0)
        meta_layout.addWidget(self.meta_station_name_edit, 1, 1)
        meta_layout.addWidget(self.meta_lat_label, 2, 0)
        meta_layout.addWidget(self.meta_lat_edit, 2, 1)
        meta_layout.addWidget(self.meta_lon_label, 3, 0)
        meta_layout.addWidget(self.meta_lon_edit, 3, 1)
        meta_layout.addWidget(self.meta_elev_label, 4, 0)
        meta_layout.addWidget(self.meta_elev_edit, 4, 1)
        meta_layout.addWidget(self.meta_date_label, 5, 0)
        meta_layout.addWidget(self.meta_date_edit, 5, 1)
        meta_layout.addWidget(self.meta_acq_label, 6, 0)
        meta_layout.addWidget(self.meta_acq_edit, 6, 1)
        
        ## static shift
        ss_layout = QtGui.QGridLayout()
        ss_layout.addWidget(self.static_shift_label, 0, 0, 1, 2)
        ss_layout.addWidget(self.static_shift_apply_button, 0, 2, 1, 2)
        ss_layout.addWidget(self.static_shift_x_label, 1, 0)
        ss_layout.addWidget(self.static_shift_x_edit, 1, 1)        
        ss_layout.addWidget(self.static_shift_y_label, 1, 2)
        ss_layout.addWidget(self.static_shift_y_edit, 1, 3) 
        
        ## rotation
        rot_layout = QtGui.QGridLayout()
        rot_layout.addWidget(self.rotate_data_label, 0, 0)
        rot_layout.addWidget(self.rotate_angle_button, 0, 1)
        rot_layout.addWidget(self.rotate_explanation, 1, 0, 1, 2)
        rot_layout.addWidget(self.rotate_angle_z_label, 2, 0)
        rot_layout.addWidget(self.rotate_angle_z_edit, 2, 1)
        rot_layout.addWidget(self.rotate_angle_t_label, 3, 0)
        rot_layout.addWidget(self.rotate_angle_t_edit, 3, 1)
        rot_layout.addWidget(self.rotate_estimate_strike_button, 4, 0, 1, 2)
        
        ## left panel
        info_layout = QtGui.QVBoxLayout()
        info_layout.addLayout(meta_layout)
        info_layout.addWidget(h_line)
        info_layout.addLayout(ss_layout)
        info_layout.addWidget(h_line)
        info_layout.addWidget(self.remove_distortion_label)
        info_layout.addWidget(self.remove_distortion_button)
        info_layout.addLayout(rot_layout)
        info_layout.addWidget(self.tools_label)
        info_layout.addWidget(self.edits_apply_button)
        info_layout.addWidget(self.revert_button)
        info_layout.addWidget(self.save_edits_button)
        info_layout.addItem(v_space)
        info_layout.addWidget(self.output_label)
        info_layout.addWidget(self.output_box)
            
        ## final layout
        final_layout = QtGui.QHBoxLayout()
        final_layout.addLayout(info_layout)
        final_layout.addLayout(mpl_vbox )

        
        self.setLayout(final_layout)
        self.mpl_widget.updateGeometry()
        
    def meta_edit_station(self):
        self.mt_obj.station = str(self.meta_station_name_edit.text())
        
    def meta_edit_lat(self):
        self.mt_obj.lat = float(str(self.meta_lat_edit.text()))
        self.meta_lat_edit.setText('{0:.6f}'.format(self.mt_obj.lat))
        
    def meta_edit_lon(self):
        self.mt_obj.lon = float(str(self.meta_lat_edit.text()))
        self.meta_lon_edit.setText('{0:.6f}'.format(self.mt_obj.lat))
        
    def meta_edit_elev(self):
        self.mt_obj.elev = float(str(self.meta_elev_edit.text()))
        self.meta_elev_edit.setText('{0:.6f}'.format(self.mt_obj.elev))
        
    def meta_edit_date(self):
        pass
    
    def meta_edit_acq(self):
        pass
        
    def fill_metadata(self):
        self.meta_station_name_edit.setText(self.mt_obj.station)
        self.meta_lat_edit.setText('{0:.6f}'.format(self.mt_obj.lat))
        self.meta_lon_edit.setText('{0:.6f}'.format(self.mt_obj.lon))
        self.meta_elev_edit.setText('{0:.6f}'.format(self.mt_obj.elev))
        
    def static_shift_set_x(self):
        self.static_shift_x = float(str(self.static_shift_x_edit.text()))
        self.static_shift_x_edit.setText('{0:.4f}'.format(self.static_shift_x))
        
    def static_shift_set_y(self):
        self.static_shift_y = float(str(self.static_shift_y_edit.text()))
        self.static_shift_y_edit.setText('{0:.4f}'.format(self.static_shift_y))
     
    def static_shift_apply(self):
        """
        shift apparent resistivity up or down
        """
        if self._edited_dist == False and self._edited_rot == False and \
           self._edited_mask == False:
            # be sure to apply the static shift to the original data
            new_z_obj = self._mt_obj.remove_static_shift(ss_x=self.static_shift_x,
                                                         ss_y=self.static_shift_y)
            # print the static shift applied 
            print "\n    Static shift applied to original data:"
        else:
            new_z_obj = self.mt_obj.remove_static_shift(ss_x=self.static_shift_x,
                                                        ss_y=self.static_shift_y)
            # print the static shift applied 
            print "\n    - Static shift applied to edited data:"
        
        # print the static shift applied
        print "        x = {0:.4f}, y = {1:.4f}".format(self.static_shift_x,
                                                        self.static_shift_y)          
        self._edited_ss = True
                                                
        self.mt_obj.Z = new_z_obj
        self.redraw_plot()
        
    def remove_distortion_apply(self):
        """
        remove distortion from the mt repsonse
        """
        if self._edited_dist == False and self._edited_rot == False and \
           self._edited_mask == False:
            # remove distortion from original data
            distortion, new_z_object = self._mt_obj.remove_distortion()
            print '\n    - Removed distortion from original data'
             
        else:
            # remove distortion from edited data
            distortion, new_z_object = self.mt_obj.remove_distortion()
            print '\n    - Removed distortion from edited data'

        self._edited_dist = True
        self.mt_obj.Z = new_z_object

        print '       Distortion matrix:'
        print '          | {0:+.3f}  {1:+.3f} |'.format(distortion[0, 0],
                                                        distortion[0, 1])
        print '          | {0:+.3f}  {1:+.3f} |'.format(distortion[1, 0],
                                                        distortion[1, 1]) 
        
        self.redraw_plot()
        
    def rotate_set_z_angle(self):
        """
        set z rotation angle
        
        Always rotate original data otherwiswe can get messy
        """
        self.rotate_z_angle = float(str(self.rotate_angle_z_edit.text()))
        self.rotate_angle_z_edit.setText("{0:.4f}".format(self.rotate_z_angle))
        
        self.rotate_tip_angle = float(self.rotate_z_angle)
        self.rotate_angle_t_edit.setText("{0:.4f}".format(self.rotate_tip_angle))

                                                         
    def rotate_set_t_angle(self):
        """
        rotate tipper data assuming North is 0 and E is 90
        
        Always rotate original data otherwiswe can get messy
        """
        self.rotate_tip_angle = float(str(self.rotate_angle_t_edit.text()))
        self.rotate_angle_t_edit.setText("{0:.4f}".format(self.rotate_tip_angle))

    def rotate_data_apply(self):
        """
        rotate both Z and tipper original data by the input angles
        """
        
        rot_z_obj = copy.deepcopy(self._mt_obj.Z)
        rot_z_obj.rotate(self.rotate_z_angle)
        
        self.mt_obj.Z = rot_z_obj
        
        rot_t_obj = copy.deepcopy(self._mt_obj.Tipper)
        rot_t_obj.rotate(self.rotate_tip_angle)
        
        self.mt_obj.Tipper = rot_t_obj
        
        if self._edited_ss == True:
            self.static_shift_apply()
        
        if self._edited_dist == True:
            self.remove_distortion_apply()
        
        self.redraw_plot()
        
        self._edited_rot = True
        
        print '\n   Rotated orginal data clockwise by:'
        print '      Z = {0:.3g}'.format(self.rotate_z_angle)
        print '      T = {0:.3g}'.format(self.rotate_tip_angle)
        
    def rotate_estimate_strike(self):
        """
        estimate strike from the invariants, phase tensor, and tipper if 
        measured.
        """
        z_list = [copy.deepcopy(self.mt_obj.Z)]
        t_list = [copy.deepcopy(self.mt_obj.Tipper)]
        strike_plot = plotstrike2d.PlotStrike2D(z_object_list=z_list,
                                                tipper_object_list=t_list,
                                                plot_yn='n')
        strike_plot.plot_type = 1
        if np.any(self.mt_obj.Tipper.tipper == 0) == True:
            strike_plot.plot_tipper = 'n'
        elif np.any(self.mt_obj.Tipper.tipper == 0) == False:
            strike_plot.plot_tipper = 'y'
            
        strike_plot.fold = False
        strike_plot.plot_range = 'data'
        strike_plot.plot()
        
    def edits_apply(self):
        """
        apply edits, well edits are already made, but replot without edited
        points
        """
        
        self.redraw_plot()
           
          
    def revert_back(self):
        """
        revert back to original data
        
        """
        
        self.mt_obj = copy.deepcopy(self._mt_obj)
        self.reset_parameters()
        self.redraw_plot()
        
        print '\n'        
        print '-'*35
        print "Reverted back to original input data."
        print "Reset editing parameters."
        
    def reset_parameters(self):
        self.static_shift_x = 1.0
        self.static_shift_y = 1.0
        self.static_shift_x_edit.setText("{0:.4g}".format(self.static_shift_x))
        self.static_shift_y_edit.setText("{0:.4g}".format(self.static_shift_y))
        self.rotate_z_angle = 0.0
        self.rotate_tip_angle = 0.0
        self.rotate_angle_z_edit.setText("{0:.4g}".format(self.rotate_z_angle))
        self.rotate_angle_t_edit.setText("{0:.4g}".format(self.rotate_tip_angle))
        
        self._edited_ss = False
        self._edited_dist = False
        self._edited_rot = False
        
    def save_edi_file(self):
        """
        save edited edi file to the chosen file name.
        """
        
        save_dialog = QtGui.QFileDialog()
        save_fn = str(save_dialog.getSaveFileName(caption='Choose EDI file',
                                                  directory=self.dir_path,
                                                  filter='*.edi'))
        self.mt_obj.write_edi_file(new_fn=save_fn)
        
    @QtCore.pyqtSlot(str)
    def normal_output(self, message):
        self.output_box.moveCursor(QtGui.QTextCursor.End)
        self.output_box.insertPlainText(message)
        
    def plot(self):
        """
        plot the response
        
        will have original data below the editable data
        """
        
        self.figure.clf()
        
        self.figure.suptitle(self.mt_obj.station, 
                             fontdict={'size':self.plot_properties.fs+4,
                                       'weight':'bold'})
        
        # make some useful internal variabls
        font_dict = {'size':self.plot_properties.fs+2, 'weight':'bold'}
        plot_period = 1./self.mt_obj.Z.freq
            
        if np.all(self.mt_obj.Tipper.tipper == 0) == True:
                print 'No Tipper data for station {0}'.format(self.mt_obj.station)
                self.plot_tipper = False
        else:
            self.plot_tipper = True
            
        #set x-axis limits from short period to long period
        self.plot_properties.xlimits = (10**(np.floor(np.log10(plot_period.min()))),
                                        10**(np.ceil(np.log10((plot_period.max())))))
        
        #--> make key word dictionaries for plotting
        kw_xx = {'color':self.plot_properties.cted,
                 'marker':self.plot_properties.mted,
                 'ms':self.plot_properties.ms,
                 'ls':':',
                 'lw':self.plot_properties.lw,
                 'e_capsize':self.plot_properties.e_capsize,
                 'e_capthick':self.plot_properties.e_capthick,
                 'picker':3}        
       
        kw_yy = {'color':self.plot_properties.ctmd,
                 'marker':self.plot_properties.mtmd,
                 'ms':self.plot_properties.ms,
                 'ls':':',
                 'lw':self.plot_properties.lw,
                 'e_capsize':self.plot_properties.e_capsize,
                 'e_capthick':self.plot_properties.e_capthick,
                 'picker':3} 
                 
        kw_xx_o = {'color':self.plot_properties.cteo,
                   'marker':self.plot_properties.mted,
                   'ms':self.plot_properties.ms,
                   'ls':':',
                   'lw':self.plot_properties.lw,
                   'e_capsize':self.plot_properties.e_capsize,
                   'e_capthick':self.plot_properties.e_capthick,
                   'picker':None}        
       
        kw_yy_o = {'color':self.plot_properties.ctmo,
                   'marker':self.plot_properties.mtmd,
                   'ms':self.plot_properties.ms,
                   'ls':':',
                   'lw':self.plot_properties.lw,
                   'e_capsize':self.plot_properties.e_capsize,
                   'e_capthick':self.plot_properties.e_capthick,
                   'picker':None} 
        
        # create a grid of subplots
        gs = gridspec.GridSpec(3, 2, height_ratios=[2, 1.5, 1])
        
        gs.update(hspace=self.plot_properties.subplot_hspace,
                  wspace=self.plot_properties.subplot_wspace,
                  left=self.plot_properties.subplot_left,
                  right=self.plot_properties.subplot_right,
                  top=self.plot_properties.subplot_top,
                  bottom=self.plot_properties.subplot_bottom)
                  
        #find locations where points have been masked
        nzxx = np.nonzero(self.mt_obj.Z.z[:, 0, 0])[0]
        nzxy = np.nonzero(self.mt_obj.Z.z[:, 0, 1])[0]
        nzyx = np.nonzero(self.mt_obj.Z.z[:, 1, 0])[0]
        nzyy = np.nonzero(self.mt_obj.Z.z[:, 1, 1])[0]
        
        nzxx_o = np.nonzero(self._mt_obj.Z.z[:, 0, 0])[0]
        nzxy_o = np.nonzero(self._mt_obj.Z.z[:, 0, 1])[0]
        nzyx_o = np.nonzero(self._mt_obj.Z.z[:, 1, 0])[0]
        nzyy_o = np.nonzero(self._mt_obj.Z.z[:, 1, 1])[0]
                  
        # make axis od = off-diagonal, d = diagonal, share x axis across plots
        self.ax_res_od = self.figure.add_subplot(gs[0, 0])
        self.ax_res_d = self.figure.add_subplot(gs[0, 1], 
                                                sharex=self.ax_res_od)

        self.ax_phase_od = self.figure.add_subplot(gs[1, 0], 
                                                sharex=self.ax_res_od)
        self.ax_phase_d = self.figure.add_subplot(gs[1, 1], 
                                                sharex=self.ax_res_od)
        
        # include tipper, r = real, i = imaginary
        self.ax_tip_m = self.figure.add_subplot(gs[2, 0], 
                                                sharex=self.ax_res_od)
        self.ax_tip_a = self.figure.add_subplot(gs[2, 1],
                                                sharex=self.ax_res_od)
        
        self.ax_list = [self.ax_res_od, self.ax_res_d,
                        self.ax_phase_od, self.ax_phase_d,
                        self.ax_tip_m, self.ax_tip_a]
        
        ## --> plot apparent resistivity, phase and tipper
        ## plot orginal apparent resistivity
        if self.plot_properties.plot_original_data == True:                                 
            orxx = mtplottools.plot_errorbar(self.ax_res_d, 
                                             plot_period[nzxx_o],
                                             self._mt_obj.Z.resistivity[nzxx_o, 0, 0],
                                             self._mt_obj.Z.resistivity_err[nzxx_o, 0, 0],
                                             **kw_xx_o)
            orxy = mtplottools.plot_errorbar(self.ax_res_od, 
                                             plot_period[nzxy_o],
                                             self._mt_obj.Z.resistivity[nzxy_o, 0, 1],
                                             self._mt_obj.Z.resistivity_err[nzxy_o, 0, 1],
                                             **kw_xx_o)
            oryx = mtplottools.plot_errorbar(self.ax_res_od, 
                                             plot_period[nzyx_o],
                                             self._mt_obj.Z.resistivity[nzyx_o, 1, 0],
                                             self._mt_obj.Z.resistivity_err[nzyx_o, 1, 0],
                                             **kw_yy_o)
            oryy = mtplottools.plot_errorbar(self.ax_res_d, 
                                             plot_period[nzyy_o],
                                             self._mt_obj.Z.resistivity[nzyy_o, 1, 1],
                                             self._mt_obj.Z.resistivity_err[nzyy_o, 1, 1],
                                             **kw_yy_o)
            # plot original phase                                 
            epxx = mtplottools.plot_errorbar(self.ax_phase_d, 
                                             plot_period[nzxx_o],
                                             self._mt_obj.Z.phase[nzxx_o, 0, 0],
                                             self._mt_obj.Z.phase_err[nzxx_o, 0, 0],
                                             **kw_xx_o)
            epxy = mtplottools.plot_errorbar(self.ax_phase_od, 
                                             plot_period[nzxy_o],
                                             self._mt_obj.Z.phase[nzxy_o, 0, 1],
                                             self._mt_obj.Z.phase_err[nzxy_o, 0, 1],
                                             **kw_xx_o)
            epyx = mtplottools.plot_errorbar(self.ax_phase_od, 
                                             plot_period[nzyx_o],
                                             self._mt_obj.Z.phase[nzyx_o, 1, 0]+180,
                                             self._mt_obj.Z.phase_err[nzyx_o, 1, 0],
                                             **kw_yy_o)
            epyy = mtplottools.plot_errorbar(self.ax_phase_d, 
                                             plot_period[nzyy_o],
                                             self._mt_obj.Z.phase[nzyy_o, 1, 1],
                                             self._mt_obj.Z.phase_err[nzyy_o, 1, 1],
                                             **kw_yy_o)
                                         
        # plot manipulated data apparent resistivity
        erxx = mtplottools.plot_errorbar(self.ax_res_d, 
                                         plot_period[nzxx],
                                         self.mt_obj.Z.resistivity[nzxx, 0, 0],
                                         self.mt_obj.Z.resistivity_err[nzxx, 0, 0],
                                         **kw_xx)
        erxy = mtplottools.plot_errorbar(self.ax_res_od, 
                                         plot_period[nzxy],
                                         self.mt_obj.Z.resistivity[nzxy, 0, 1],
                                         self.mt_obj.Z.resistivity_err[nzxy, 0, 1],
                                         **kw_xx)
        eryx = mtplottools.plot_errorbar(self.ax_res_od, 
                                         plot_period[nzyx],
                                         self.mt_obj.Z.resistivity[nzyx, 1, 0],
                                         self.mt_obj.Z.resistivity_err[nzyx, 1, 0],
                                         **kw_yy)
        eryy = mtplottools.plot_errorbar(self.ax_res_d, 
                                         plot_period[nzyy],
                                         self.mt_obj.Z.resistivity[nzyy, 1, 1],
                                         self.mt_obj.Z.resistivity_err[nzyy, 1, 1],
                                         **kw_yy)

                                         
        #--> set axes properties for apparent resistivity
        for aa, ax in enumerate([self.ax_res_od, self.ax_res_d]):
            plt.setp(ax.get_xticklabels(), visible=False)
            if aa == 0:            
                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                              fontdict=font_dict)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_xlim(self.plot_properties.xlimits)

            ax.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)
            # make the plot cleaner by removing the bottom x label
            ylim = ax.get_ylim()
            ylimits = (10**np.floor(np.log10(ylim[0])), 
                      10**np.ceil(np.log10(ylim[1])))
            ax.set_ylim(ylimits)
            ylabels = [' ', ' ']+\
                     [mtplottools.labeldict[ii] for ii 
                     in np.arange(np.log10(ylimits[0])+1, 
                                  np.log10(ylimits[1])+1, 1)]
            ax.set_yticklabels(ylabels)
            
        self.ax_res_od.legend((erxy[0], eryx[0]), 
                              ('$Z_{xy}$', '$Z_{yx}$'),
                                loc=3, 
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02)
        self.ax_res_d.legend((erxx[0], eryy[0]), 
                              ('$Z_{xx}$', '$Z_{yy}$'),
                                loc=3, 
                                markerscale=1, 
                                borderaxespad=.01,
                                labelspacing=.07, 
                                handletextpad=.2, 
                                borderpad=.02)
        self.ax_res_d.set_ylim(self.plot_properties.res_limits_d)
        self.ax_res_od.set_ylim(self.plot_properties.res_limits_od)
        
        ##--> plot phase                           
        # plot manipulated data                         
        epxx = mtplottools.plot_errorbar(self.ax_phase_d, 
                                         plot_period[nzxx],
                                         self.mt_obj.Z.phase[nzxx, 0, 0],
                                         self.mt_obj.Z.phase_err[nzxx, 0, 0],
                                         **kw_xx)
        epxy = mtplottools.plot_errorbar(self.ax_phase_od, 
                                         plot_period[nzxy],
                                         self.mt_obj.Z.phase[nzxy, 0, 1],
                                         self.mt_obj.Z.phase_err[nzxy, 0, 1],
                                         **kw_xx)
        epyx = mtplottools.plot_errorbar(self.ax_phase_od, 
                                         plot_period[nzyx],
                                         self.mt_obj.Z.phase[nzyx, 1, 0]+180,
                                         self.mt_obj.Z.phase_err[nzyx, 1, 0],
                                         **kw_yy)
        epyy = mtplottools.plot_errorbar(self.ax_phase_d, 
                                         plot_period[nzyy],
                                         self.mt_obj.Z.phase[nzyy, 1, 1],
                                         self.mt_obj.Z.phase_err[nzyy, 1, 1],
                                         **kw_yy)
        
                                         
        #--> set axes properties
        for aa, ax in enumerate([self.ax_phase_od, self.ax_phase_d]):
            ax.set_xlabel('Period (s)', font_dict)
            ax.set_xscale('log')
            if aa == 0:
                ax.set_ylabel('Phase (deg)', font_dict)
                #ax.yaxis.set_major_locator(MultipleLocator(15))
                #ax.yaxis.set_minor_locator(MultipleLocator(5))
            ax.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)
        self.ax_phase_od.set_ylim(self.plot_properties.phase_limits_od)
        self.ax_phase_d.set_ylim(self.plot_properties.phase_limits_d)
        
        # set the last label to be an empty string for easier reading
        for ax in [self.ax_phase_od, self.ax_phase_d]:
            y_labels = ax.get_yticks().tolist()
            y_labels[0] = ''
            ax.set_yticklabels(y_labels)

        ## --> plot tipper                                 
        #set th xaxis tick labels to invisible
        plt.setp(self.ax_phase_od.xaxis.get_ticklabels(), visible=False)
        plt.setp(self.ax_phase_d.xaxis.get_ticklabels(), visible=False)
        self.ax_phase_od.set_xlabel('')
        self.ax_phase_d.set_xlabel('')
        
        kw_xx['color'] = self.plot_properties.arrow_color_real
        kw_yy['color'] = self.plot_properties.arrow_color_imag
        
        # plot magnitude of real and imaginary induction vectors
        if self.plot_properties.plot_original_data == True:
            etro = mtplottools.plot_errorbar(self.ax_tip_m,
                                             plot_period,
                                             self._mt_obj.Tipper.amplitude[:, 0, 0],
                                             self._mt_obj.Tipper.amplitude_err[:, 0, 0],
                                             **kw_xx_o)
#            etio = mtplottools.plot_errorbar(self.ax_tip_m,
#                                             plot_period,
#                                             self._mt_obj.Tipper.tipper[:, 0, 0].imag,
#                                             self._mt_obj.Tipper.tippererr[:, 0, 0],
#                                             **kw_yy_o) 
            ## plot angle  of original data                                
            etro = mtplottools.plot_errorbar(self.ax_tip_a,
                                             plot_period,
                                             self._mt_obj.Tipper.amplitude[:, 0, 1],
                                             self._mt_obj.Tipper.amplitude_err[:, 0, 1],
                                             **kw_xx_o)
#            etio = mtplottools.plot_errorbar(self.ax_tip_a,
#                                             plot_period,
#                                             self._mt_obj.Tipper.tipper[:, 0, 1].imag,
#                                             self._mt_obj.Tipper.tippererr[:, 0, 1],
#                                             **kw_yy_o)
#        
        # plot magnitude of edited induction vectors
        etr = mtplottools.plot_errorbar(self.ax_tip_m,
                                        plot_period,
                                        self.mt_obj.Tipper.amplitude[:, 0, 0],
                                        self.mt_obj.Tipper.amplitude_err[:, 0, 0],
                                        **kw_xx) 
#        eti = mtplottools.plot_errorbar(self.ax_tip_m,
#                                        plot_period,
#                                        self.mt_obj.Tipper.tipper[:, 0, 0].imag,
#                                        self.mt_obj.Tipper.tippererr[:, 0, 0],
#                                        **kw_yy) 
#        # plot angle of edited data
        etr = mtplottools.plot_errorbar(self.ax_tip_a,
                                        plot_period,
                                        self.mt_obj.Tipper.amplitude[:, 0, 1],
                                        self.mt_obj.Tipper.amplitude_err[:, 0, 1],
                                        **kw_xx) 
#        eti = mtplottools.plot_errorbar(self.ax_tip_a,
#                                        plot_period,
#                                        self.mt_obj.Tipper.tipper[:, 0, 1].imag,
#                                        self.mt_obj.Tipper.tippererr[:, 0, 1],
#                                        **kw_yy)
                                        
        #--> set axes properties for magnitude and angle of induction vectors
        for aa, ax in enumerate([self.ax_tip_m, self.ax_tip_a]):
            if aa == 0:            
                ax.set_ylabel('Magnitude', fontdict=font_dict)
#            elif aa == 1:
#                ax.set_ylabel('Angle (deg)', fontdict=font_dict)
            
            ax.set_xlim(self.plot_properties.xlimits)
            ax.set_xscale('log')
            ax.grid(True, alpha=.25, which='both', color=(.25, .25, .25),
                          lw=.25)
            
#        self.ax_tip_m.legend((etr[0], eti[0]), 
#                              ('Re{T}', 'Im{T}'),
#                                loc=2,
#                                markerscale=1, 
#                                borderaxespad=.01,
#                                labelspacing=.07, 
#                                handletextpad=.2, 
#                                borderpad=.04)
                                
                                
        self.ax_tip_m.set_ylim(self.plot_properties.tipper_mag_limits)
        self.ax_tip_a.set_ylim(self.plot_properties.tipper_angle_limits)
        
        # set the last label to be an empty string for easier reading
        for ax in [self.ax_tip_m, self.ax_tip_a]:
            y_labels = ax.get_yticks().tolist()
            y_labels[-1] = ''
            ax.set_yticklabels(y_labels)                             
        
        ## --> need to be sure to draw the figure        
        self.mpl_widget.draw()
        
    def redraw_plot(self):
        self.figure.clf()
        self.plot()
        
    def on_pick(self, event):
        """
        mask a data point when it is clicked on.  
        """         
        data_point = event.artist
        data_period = data_point.get_xdata()[event.ind]
        data_value = data_point.get_ydata()[event.ind]
        
        mask_kw = {'color' : self.plot_properties.mask_color,
                   'marker' : self.plot_properties.mask_marker,
                   'ms' : self.plot_properties.mask_ms,
                   'mew' : self.plot_properties.mask_mew}
        
        # modify Z
        if event.mouseevent.button == 1:
            self._edited_mask = True
            if self._ax_index == 0 or self._ax_index == 1:
                d_index = np.where(self.mt_obj.Z.resistivity == data_value)
                comp_ii = d_index[0][0]
                comp_jj = d_index[1][0]
                comp_kk = d_index[2][0]
                
                # mask point in impedance object
                self.mt_obj.Z.z[d_index] = 0.0+0.0*1j            
                
                self._ax.plot(data_period, data_value, **mask_kw)
                
                # mask phase as well
                if self._ax_index == 0:
                    if comp_jj == 1 and comp_kk == 0:
                        self.ax_phase_od.plot(data_period, 
                                              self.mt_obj.Z.phase[d_index]+180,
                                              **mask_kw)
                    else:
                        self.ax_phase_od.plot(data_period, 
                                              self.mt_obj.Z.phase[d_index],
                                              **mask_kw)
                elif self._ax_index == 1:
                    self.ax_phase_d.plot(data_period, 
                                          self.mt_obj.Z.phase[d_index],
                                          **mask_kw)
                self._ax.figure.canvas.draw()
                
            # mask phase points
            elif self._ax_index == 2 or self._ax_index == 3:
                try:
                    d_index = np.where(self.mt_obj.Z.phase == data_value)
                    comp_ii = d_index[0][0]
                    comp_jj = d_index[1][0]
                    comp_kk = d_index[2][0]
                except IndexError:
                    d_index = np.where(self.mt_obj.Z.phase == data_value-180)
                    comp_ii = d_index[0][0]
                    comp_jj = d_index[1][0]
                    comp_kk = d_index[2][0]    
                
                # mask point in impedance object
                self.mt_obj.Z.z[d_index] = 0.0+0.0*1j            
                
                # mask the point in the axis selected
                self._ax.plot(data_period, data_value, **mask_kw)
                
                # mask resistivity as well
                if self._ax_index == 2:
                    self.ax_res_od.plot(data_period, 
                                        self.mt_obj.Z.resistivity[d_index],
                                        **mask_kw)
                elif self._ax_index == 3:
                    self.ax_res_d.plot(data_period, 
                                       self.mt_obj.Z.resistivity[d_index],
                                       **mask_kw)
                self._ax.figure.canvas.draw()
            
#            # mask tipper points
#            elif self._ax_index == 4:
#                data_value = np.round(data_value, 8)
#                # the tipper is only rank on so jj will always be 0
#                comp_jj = 0
#
#                # find the tipper magnitude just picked, need to test for both
#                # real and imaginary, then real again incase imaginary is
#                # picked, probably a better way to do this.
#                try:
#                    d_index = np.where(np.round(self.mt_obj.Tipper.mag_real, 
#                                                8) == data_value)
#                    comp_ii = d_index[0][0]
#                    comp_kk = 0
#                except IndexError:
#                    try:
#                        d_index = np.where(np.round(self.mt_obj.Tipper.mag_imag,
#                                                    8) == data_value)
#                        comp_ii = d_index[0][0]
#                        comp_kk = 1
#                    except IndexError:
#                        d_index = np.where(np.round(self.mt_obj.Tipper.mag_real,
#                                                    8) == data_value)
#                        comp_ii = d_index[0][0]
#                        comp_kk = 0
#                
#            
#                
#                self._ax.plot(data_period, data_value, **mask_kw)
#                
#                #plot in cooresponding direction mode
#                if comp_kk == 0:
#                    plot_angle =  self.mt_obj.Tipper.angle_real[comp_ii]%360
#                elif comp_kk == 1:
#                    plot_angle =  self.mt_obj.Tipper.angle_imag[comp_ii]%360
#                self.ax_tip_a.plot(data_period, plot_angle, **mask_kw)
#                
#                self._ax.figure.canvas.draw()
#                
#                # mask the point
#                self.mt_obj.Tipper.tipper[comp_ii, comp_jj, comp_kk] = 0.0+0.0j                
#                self.mt_obj.Tipper._compute_mag_direction()
                
    def in_axes(self, event):
        """
        check to see which axis the mouse is in
        """
        
        self._ax = event.inaxes
        
        # find the component index so that it can be masked
        for ax_index, ax in enumerate(self.ax_list):
            if ax == event.inaxes:
                self._ax_index = ax_index
        

#==============================================================================
#  Plot setting        
#==============================================================================
class PlotSettings(QtGui.QWidget):
    settings_updated = QtCore.pyqtSignal()
    def __init__(self, parent, **kwargs):
        super(PlotSettings, self).__init__(parent)
        
        self.fs = kwargs.pop('fs', 10)
        self.lw = kwargs.pop('lw', 1.0)
        self.ms = kwargs.pop('ms', 4)
        
        self.e_capthick = kwargs.pop('e_capthick', 1)
        self.e_capsize =  kwargs.pop('e_capsize', 4)

        #color mode
        self.cted = kwargs.pop('cted', (0, 0, .65))
        self.ctmd = kwargs.pop('ctmd', (.65, 0, 0))
        self.cteo = kwargs.pop('cted', (.5, .5, .5))
        self.ctmo = kwargs.pop('ctmd', (.75, .75, .75))
        self.mted = kwargs.pop('mted', 's')
        self.mtmd = kwargs.pop('mtmd', 'o')
        
        #color for occam2d model
        self.ctem = kwargs.pop('ctem', (0, .6, .3))
        self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
        self.mtem = kwargs.pop('mtem', '+')
        self.mtmm = kwargs.pop('mtmm', '+')
         
        self.res_limits_od = kwargs.pop('res_limits_od', None)   
        self.res_limits_d = kwargs.pop('res_limits_d', None)   
        
        self.phase_limits_od = kwargs.pop('phase_limits_od', None)   
        self.phase_limits_d = kwargs.pop('phase_limits_d', None)  
        
        self.tipper_mag_limits = kwargs.pop('tipper_mag_limits', None)
        self.tipper_angle_limits = kwargs.pop('tipper_angle_limits', None)
        
        self.subplot_wspace = kwargs.pop('subplot_wspace', .15)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .00)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .93)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .08)
        
        self.legend_loc = kwargs.pop('legend_loc', 'upper center')
        self.legend_pos = kwargs.pop('legend_pos', (.5, 1.11))
        self.legend_marker_scale = kwargs.pop('legend_marker_scale', 1)
        self.legend_border_axes_pad = kwargs.pop('legend_border_axes_pad', .01)
        self.legend_label_spacing = kwargs.pop('legend_label_spacing', 0.07)
        self.legend_handle_text_pad = kwargs.pop('legend_handle_text_pad', .2)
        self.legend_border_pad = kwargs.pop('legend_border_pad', .15)
        
        self.arrow_direction = kwargs.pop('arrow_direction', 0)
        self.arrow_color_real = kwargs.pop('arrow_color_real', (.4, 0, .2))
        self.arrow_color_imag = kwargs.pop('arrow_color_imag',  (0, .9, .9))
        self.arrow_head_width = kwargs.pop('arrow_head_width', .05)
        self.arrow_head_length = kwargs.pop('arrow_head_length', .05)
        
        self.plot_z = kwargs.pop('plot_z', False)
        self.plot_original_data = kwargs.pop('plot_original_data', True)
        
        self.mask_marker = kwargs.pop('mask_marker', 'x')
        self.mask_ms = kwargs.pop('mask_ms', 10)
        self.mask_mew = kwargs.pop('mask_mew', 3)
        self.mask_color = kwargs.pop('mask_color', 'k')
        
        #self.setup_ui()

    def setup_ui(self):
        #--> line properties
        fs_label = QtGui.QLabel('Font Size')
        fs_edit = QtGui.QLineEdit()
        fs_edit.setText('{0:.1f}'.format(self.fs))
        fs_edit.textChanged[str].connect(self.set_text_fs)
        
        lw_label = QtGui.QLabel('Line Width')
        lw_edit = QtGui.QLineEdit()
        lw_edit.setText('{0:.1f}'.format(self.lw))
        lw_edit.textChanged[str].connect(self.set_text_lw)
        
        e_capthick_label = QtGui.QLabel('Error cap thickness')
        e_capthick_edit = QtGui.QLineEdit()
        e_capthick_edit.setText('{0:.1f}'.format(self.e_capthick))
        e_capthick_edit.textChanged[str].connect(self.set_text_e_capthick)
        
        e_capsize_label = QtGui.QLabel('Error cap size')
        e_capsize_edit = QtGui.QLineEdit()
        e_capsize_edit.setText('{0:.1f}'.format(self.e_capsize))
        e_capsize_edit.textChanged[str].connect(self.set_text_e_capsize)
        
        grid_line = QtGui.QGridLayout()
        grid_line.setSpacing(10)
        
        grid_line.addWidget(fs_label, 1, 0)
        grid_line.addWidget(fs_edit, 1, 1)
        
        grid_line.addWidget(lw_label, 1, 2)
        grid_line.addWidget(lw_edit, 1, 3)
        
        grid_line.addWidget(e_capthick_label, 1, 4)
        grid_line.addWidget(e_capthick_edit, 1, 5)
        
        grid_line.addWidget(e_capsize_label, 1, 6)
        grid_line.addWidget(e_capsize_edit, 1, 7)
        
        #--> marker properties
        ms_label = QtGui.QLabel('Marker Size')
        ms_edit = QtGui.QLineEdit()
        ms_edit.setText('{0:.1f}'.format(self.ms))
        ms_edit.textChanged[str].connect(self.set_text_ms)
        
        dcxy_label = QtGui.QLabel('Data Color xy')
        dcxy_edit = QtGui.QLineEdit()
        dcxy_edit.setText('{0}'.format(self.cted))
        dcxy_edit.textChanged[str].connect(self.set_text_cted)
        
        dcyx_label = QtGui.QLabel('Data Color yx')
        dcyx_edit = QtGui.QLineEdit()
        dcyx_edit.setText('{0}'.format(self.ctmd))
        dcyx_edit.textChanged[str].connect(self.set_text_ctmd)
        
        dmxy_label = QtGui.QLabel('Data Marker xy')
        dmxy_edit = QtGui.QLineEdit()
        dmxy_edit.setText('{0}'.format(self.mted))
        dmxy_edit.textChanged[str].connect(self.set_text_mted)
        
        dmyx_label = QtGui.QLabel('Data Marker yx')
        dmyx_edit = QtGui.QLineEdit()
        dmyx_edit.setText('{0}'.format(self.mtmd))
        dmyx_edit.textChanged[str].connect(self.set_text_mtmd)
        
        mcxy_label = QtGui.QLabel('Model Color xy')
        mcxy_edit = QtGui.QLineEdit()
        mcxy_edit.setText('{0}'.format(self.ctem))
        mcxy_edit.textChanged[str].connect(self.set_text_ctem)
        
        mcyx_label = QtGui.QLabel('Model Color yx')
        mcyx_edit = QtGui.QLineEdit()
        mcyx_edit.setText('{0}'.format(self.ctmm))
        mcyx_edit.textChanged[str].connect(self.set_text_ctmm)
        
        mmxy_label = QtGui.QLabel('Model Marker xy')
        mmxy_edit = QtGui.QLineEdit()
        mmxy_edit.setText('{0}'.format(self.mtem))
        mmxy_edit.textChanged[str].connect(self.set_text_mtem)
    
        mmyx_label = QtGui.QLabel('Model Marker yx')
        mmyx_edit = QtGui.QLineEdit()
        mmyx_edit.setText('{0}'.format(self.mtmm))
        mmyx_edit.textChanged[str].connect(self.set_text_mtmm)

        marker_label = QtGui.QLabel('Maker Properties:')
        
        marker_grid = QtGui.QGridLayout()
        marker_grid.setSpacing(10)
                
        marker_grid.addWidget(marker_label, 1, 0)
        marker_grid.addWidget(ms_label, 1, 2)
        marker_grid.addWidget(ms_edit, 1, 3)
        
        marker_grid.addWidget(dcxy_label, 2, 0)
        marker_grid.addWidget(dcxy_edit, 2, 1)
        
        marker_grid.addWidget(dcyx_label, 2, 2)
        marker_grid.addWidget(dcyx_edit, 2, 3)
        
        marker_grid.addWidget(dmxy_label, 2, 4)
        marker_grid.addWidget(dmxy_edit, 2, 5)
        
        marker_grid.addWidget(dmyx_label, 2, 6)
        marker_grid.addWidget(dmyx_edit, 2, 7)
        
        marker_grid.addWidget(mcxy_label, 3, 0)
        marker_grid.addWidget(mcxy_edit, 3, 1)
        
        marker_grid.addWidget(mcyx_label, 3, 2)
        marker_grid.addWidget(mcyx_edit, 3, 3)
        
        marker_grid.addWidget(mmxy_label, 3, 4)
        marker_grid.addWidget(mmxy_edit, 3, 5)
        
        marker_grid.addWidget(mmyx_label, 3, 6)
        marker_grid.addWidget(mmyx_edit, 3, 7)
        
        #--> plot limits
        ylimr_od_label = QtGui.QLabel('Res_od')
        ylimr_od_edit = QtGui.QLineEdit()
        ylimr_od_edit.setText('{0}'.format(self.res_limits_od))
        ylimr_od_edit.textChanged[str].connect(self.set_text_res_od) 
        
        ylimr_d_label = QtGui.QLabel('Res_d')
        ylimr_d_edit = QtGui.QLineEdit()
        ylimr_d_edit.setText('{0}'.format(self.res_limits_d))
        ylimr_d_edit.textChanged[str].connect(self.set_text_res_d)  
        
        ylimp_od_label = QtGui.QLabel('phase_od')
        ylimp_od_edit = QtGui.QLineEdit()
        ylimp_od_edit.setText('{0}'.format(self.phase_limits_od))
        ylimp_od_edit.textChanged[str].connect(self.set_text_phase_od) 
        
        ylimp_d_label = QtGui.QLabel('phase_d')
        ylimp_d_edit = QtGui.QLineEdit()
        ylimp_d_edit.setText('{0}'.format(self.phase_limits_d))
        ylimp_d_edit.textChanged[str].connect(self.set_text_phase_d) 
        
        limits_grid = QtGui.QGridLayout()
        limits_grid.setSpacing(10)
        
        limits_label = QtGui.QLabel('Plot Limits: (Res=Real, Phase=Imaginary)'
                                    ' --> input on a linear scale')
        
        limits_grid.addWidget(limits_label, 1, 0, 1, 7)
        
        limits_grid.addWidget(ylimr_od_label, 2, 0)
        limits_grid.addWidget(ylimr_od_edit, 2, 1)
        limits_grid.addWidget(ylimr_d_label, 2, 2)
        limits_grid.addWidget(ylimr_d_edit, 2, 3)
        
        limits_grid.addWidget(ylimp_od_label, 3, 0)
        limits_grid.addWidget(ylimp_od_edit, 3, 1)
        limits_grid.addWidget(ylimp_d_label, 3, 2)
        limits_grid.addWidget(ylimp_d_edit, 3, 3)
        
        #--> legend properties
        legend_pos_label = QtGui.QLabel('Legend Position')
        legend_pos_edit = QtGui.QLineEdit()
        legend_pos_edit.setText('{0}'.format(self.legend_pos))
        legend_pos_edit.textChanged[str].connect(self.set_text_legend_pos)
        
        legend_grid = QtGui.QGridLayout()
        legend_grid.setSpacing(10)
        
        legend_grid.addWidget(QtGui.QLabel('Legend Properties:'), 1, 0)
        legend_grid.addWidget(legend_pos_label, 1, 2,)
        legend_grid.addWidget(legend_pos_edit, 1, 3)
        
        update_button = QtGui.QPushButton('Update')
        update_button.clicked.connect(self.update_settings)        
        
        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid_line)
        vbox.addLayout(marker_grid)
        vbox.addLayout(limits_grid)
        vbox.addLayout(legend_grid)
        vbox.addWidget(update_button)
        
        self.setLayout(vbox) 
        
        self.setGeometry(300, 300, 350, 300)
        self.resize(1350, 500)
        self.setWindowTitle('Plot Settings')    
        self.show()

    def set_text_fs(self, text):
        try:
            self.fs = float(text)
        except ValueError:
            print "Enter a float point number"
            
    def set_text_e_capthick(self, text):
        try:
            self.e_capthick = float(text)
        except ValueError:
            print "Enter a float point number"
            
    def set_text_e_capsize(self, text):
        try:
            self.e_capsize = float(text)
        except ValueError:
            print "Enter a float point number"

    
    def set_text_lw(self, text):
        try:
            self.lw = float(text)
        except ValueError:
            print "Enter a float point number"
            
    def set_text_ms(self, text):
        try:
            self.ms = float(text)
        except ValueError:
            print "Enter a float point number"
            
    def set_text_cted(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 3:
            print 'enter as (r, g, b)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 3:
            self.cted = tuple(l_list)
            
    def set_text_ctmd(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 3:
            print 'enter as (r, g, b)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 3:
            self.ctmd = tuple(l_list)
            
    def set_text_mted(self, text):
        try:
            self.mted = str(text)
        except ValueError:
            print "Enter a string"
            
    def set_text_mtmd(self, text):
        try:
            self.mtmd = str(text)
        except ValueError:
            print "Enter a string"
            
    def set_text_ctem(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 3:
            print 'enter as (r, g, b)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 3:
            self.ctem = tuple(l_list)
    
    def set_text_ctmm(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 3:
            print 'enter as (r, g, b)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 3:
            self.ctmm = tuple(l_list)
            
    def set_text_mtem(self, text):
        try:
            self.mtem = str(text)
        except ValueError:
            print "Enter a string"
            
    def set_text_mtmm(self, text):
        try:
            self.mtmm = str(text)
        except ValueError:
            print "Enter a string"
            
    def set_text_res_od(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print 'enter as (min, max)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.res_limits_od = tuple(l_list)
            
    def set_text_res_d(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print 'enter as (min, max)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.res_limits_d = tuple(l_list)
            
    def set_text_phase_od(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print 'enter as (min, max)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.phase_limits_od = tuple(l_list)
            
    def set_text_phase_d(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print 'enter as (min, max)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.phase_limits_d = tuple(l_list)
            
            
    def set_text_legend_pos(self, text):
        if text =='None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print 'enter as (min, max)'
        l_list = []
        for txt in t_list:
            try: 
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.legend_pos = tuple(l_list)
            
    def update_settings(self):
        self.settings_updated.emit()
        
#==============================================================================
# Def Main
#==============================================================================
def main():
    app = QtGui.QApplication(sys.argv)
    ui = EDI_Editor_Window()
    ui.ui_setup()
    ui.show()
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    main()