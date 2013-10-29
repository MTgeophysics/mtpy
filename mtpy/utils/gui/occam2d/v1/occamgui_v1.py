#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys,os
from PyQt4 import QtCore, QtGui
import gui5
reload(gui5)

from gui5 import Ui_occamgui2D as Occam_UI_form

import os.path as op
import shutil

import mtpy.core.edi as MTedi
import mtpy.utils.configfile as MTcf
import mtpy.modeling.occam2d as MTo2



class OccamGui(QtGui.QMainWindow):
    """
    Class for handling OCCAM GUI


    - Connect input fields and buttons.
    - Check for valid input values
    - Generate OCCAM input files
    - Load old/temporary OCCAM run outputs
    - Run OCCAM

    """


    def __init__(self, parent=None):
        
        QtGui.QWidget.__init__(self, parent)
        """
        Call the general QWidget init, then add the connections and slots.
        
        """

        #load the occam gui source module
        self.ui = Occam_UI_form()
        self.ui.setupUi(self)

        #set most basic parameters
        self._load_datafile = None
        self.parameters = {}
        self.ui.wd = op.abspath(op.realpath('.'))

        #Connections
        QtCore.QObject.connect(self.ui.button_browse_wd, QtCore.SIGNAL("clicked()"),  lambda: self.set_path_in_browsefield(self.ui.lineEdit_browse_wd))
        QtCore.QObject.connect(self.ui.button_browse_edis, QtCore.SIGNAL("clicked()"),  lambda: self.set_path_in_browsefield(self.ui.lineEdit_browse_edi))
        QtCore.QObject.connect(self.ui.pushButton_loadstations, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_stations))
        QtCore.QObject.connect(self.ui.button_browse_occam, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_occam))
        #QtCore.QObject.connect(self.ui.button_browse_makemodel, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_makemodel))
        QtCore.QObject.connect(self.ui.pushButton_loaddatafile, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_datafile))
        QtCore.QObject.connect(self.ui.pushButton_loaditerationfile, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_iterationfile))
        
        QtCore.QObject.connect(self.ui.button_browse_configfile, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_configfile))
        QtCore.QObject.connect(self.ui.button_load_configfile, QtCore.SIGNAL("clicked()"),self.load_old_configfile)

        QtCore.QObject.connect(self.ui.pushButton_loaddatafile, QtCore.SIGNAL("clicked()"),self.set_data_filename)
        
        QtCore.QObject.connect(self.ui.pushButton_checkparameter, QtCore.SIGNAL("clicked()"),  self.check_input)
        #QtCore.QObject.connect(self.ui.pushButton_checkparameter, QtCore.SIGNAL("clicked()"),  self.setup_parameter_dict)
        #QtCore.QObject.connect(self.ui.pushButton_runoccam, QtCore.SIGNAL("clicked()"),  self.setup_parameter_dict)
        #QtCore.QObject.connect(self.ui.pushButton_generateinputfile, QtCore.SIGNAL("clicked()"),  self.check_input)
        #QtCore.QObject.connect(self.ui.pushButton_runoccam, QtCore.SIGNAL("clicked()"),  self.check_input)
        
        QtCore.QObject.connect(self.ui.pushButton_generateinputfile, QtCore.SIGNAL("clicked()"),  self.generate_inputfiles)
        #QtCore.QObject.connect(self.ui.pushButton_loadold_meshinmodel, QtCore.SIGNAL("clicked()"),  self.loadold_meshinmodel)
        QtCore.QObject.connect(self.ui.pushButton_runoccam, QtCore.SIGNAL("clicked()"),  self.run_occam)
        QtCore.QObject.connect(self.ui.pushButton_quit, QtCore.SIGNAL("clicked()"), QtCore.QCoreApplication.instance().quit)
        #QtCore.QObject.connect(self.ui.pushButton_generateinputfile, QtCore.SIGNAL("clicked()"),  self.setup_parameter_dict)

        
    #when loading old datafile, copy its name to the filename field:
    def set_data_filename(self):
        self.ui.lineEdit_datafilename.setText(self.ui.lineEdit_browse_datafile.text())
        self._load_datafile = 1
        

    #use standard PyQt tool to browse and set existing directory
    def set_path_in_browsefield(self, browsefield):
        dirname = QtGui.QFileDialog.getExistingDirectory(self, 'Open Directory', '.')
        browsefield.setText(dirname)

    #use standard PyQt tool to browse and set existing file
    def set_filename_in_browsefield(self, browsefield):
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Locate File', '.')
        browsefield.setText(filename)


    def load_old_configfile(self):
        old_cfg_filename =  self.ui.lineEdit_browse_configfile.text()
        #if not a proper file: do nothing
        try:
            if not op.isfile(old_cfg_filename):
                raise
        except:
            messagetext = ''
            messagetext += "<P><FONT COLOR='#000000'>File name: "\
                    "{0}  </FONT></P> \n".format(old_cfg_filename)
            messagetext += "<P><b><FONT COLOR='#800000'>Error:  Not a valid "\
                "configuration file  </FONT></b></P> \n"

            QtGui.QMessageBox.about(self, "Reading configuration file", messagetext)
            return


        #try to read config file into dictionary:
        parameters = {}
        try:
            #to test, if file is readable:
            with open(old_cfg_filename) as F:
                data = F.read()

            temp_dict_outer = MTcf.read_configfile(old_cfg_filename)
            if len(temp_dict_outer) == 0:
                raise 
            
            for k,v in temp_dict_outer.items():
                temp_dict_inner = v
                parameters.update(temp_dict_inner)
        except:
            messagetext = ''
            messagetext += "<P><FONT COLOR='#000000'>File name: "\
                    "{0}  </FONT></P> \n".format(old_cfg_filename)
            messagetext += "<P><b><FONT COLOR='#800000'>Error: File not valid or "\
                "not readable  </FONT></b></P> \n"

            QtGui.QMessageBox.about(self, "Reading configuration file", messagetext)
            return

        #now go through all parameters and see if they are contained in the config file
        #if yes, update the values in the fields

        update_counter = 0

        if 'block_merge_threshold' in parameters:
            try:
                value = float(parameters['block_merge_threshold'])
                self.ui.doubleSpinBox_mergethreshold.setValue(value)
                update_counter += 1
            except:
                pass

        if 'datafile' in parameters:
            try:
                value = str(parameters['datafile'])
                self.ui.lineEdit_browse_datafile.setText(value)
                update_counter += 1
            except:
                pass
        if 'debug_level' in parameters:
            d = {'0':0,'1':1,'2':2 }
            try:
                value = str(int(float((parameters['debug_level'])))).lower()
                self.ui.comboBox_debuglevel.setCurrentIndex(int(d[value]))
                update_counter += 1
            except:
                pass
        if 'edi_directory' in parameters:
            try:
                value = str(parameters['edi_directory'])
                self.ui.lineEdit_browse_edi.setText(value)
                update_counter += 1
            except:
                pass
        if 'edi_type' in parameters:
            d = {'z':0,'resphase':1,'spectra':2 }
            try:
                value = str(parameters['edi_type']).lower()
                self.ui.comboBox_edi_type.setCurrentIndex(int(d[value]))
                update_counter += 1
            except:
                pass
        if 'firstlayer_thickness' in parameters:
            try:
                value = float(parameters['firstlayer_thickness'])
                self.ui.spinBox_firstlayer.setValue(value)
                update_counter += 1
            except:
                pass
        if 'halfspace_resistivity' in parameters:
            try:
                value = float(parameters['halfspace_resistivity'])
                self.ui.doubleSpinBox_rhostart.setValue(value)
                update_counter += 1
            except:
                pass
        if 'max_blockwidth' in parameters:
            try:
                value = float(parameters['max_blockwidth'])
                self.ui.spinBox_maxblockwidth.setValue(value)
                update_counter += 1
            except:
                pass
        if 'max_no_frequencies' in parameters:
            try:
                value = str(parameters['max_no_frequencies'])
                if len(value) == 0 or value.lower().strip() == 'none':
                    self.ui.checkBox_max_no_frequencies.setCheckState(0)
                    self.ui.spinBox_max_no_frequencies.setValue(0)
                else:
                    value = int(float(value))
                    self.ui.checkBox_max_no_frequencies.setCheckState(2)
                    self.ui.spinBox_max_no_frequencies.setValue(value)
                update_counter += 1
            except:
                self.ui.checkBox_max_no_frequencies.setCheckState(0)
        if 'max_frequency' in parameters:
            try:
                value = str(parameters['max_frequency'])
                if len(value) == 0 or value.lower().strip() == 'none':
                    self.ui.checkBox_max_frequency.setCheckState(0)
                    self.ui.doubleSpinBox_max_frequency.setValue(0)
                else:
                    value = int(float(value))
                    self.ui.checkBox_max_frequency.setCheckState(2)
                    self.ui.doubleSpinBox_max_frequency.setValue(value)
                update_counter += 1
            except:
                self.ui.checkBox_max_frequency.setCheckState(0)
        if 'min_frequency' in parameters:
            try:
                value = str(parameters['min_frequency'])
                if len(value) == 0 or value.lower().strip() == 'none':
                    self.ui.checkBox_min_frequency.setCheckState(0)
                    self.ui.doubleSpinBox_min_frequency.setValue(0)
                else:
                    value = int(float(value))
                    self.ui.checkBox_min_frequency.setCheckState(2)
                    self.ui.doubleSpinBox_min_frequency.setValue(value)
                update_counter += 1
            except:
                self.ui.checkBox_min_frequency.setCheckState(0)

        if 'max_no_iterations' in parameters:
            try:
                value = int(float(parameters['max_no_iterations']))
                self.ui.spinBox_max_no_iterations.setValue(value)
                update_counter += 1
            except:
                pass

        if 'mode' in parameters:
            d = {'both':0,'tm':1,'te':2, 'tipper':3, 'all':4 }
            try:
                value = None
                raw_value = str(parameters['mode']).lower()
                if 'te' in raw_value: 
                    value = 'te'
                    if 'tm' in raw_value:
                        value = 'both'
                elif 'tm' in raw_value:
                    value = 'tm'
                if 'both' in raw_value:
                    value = 'both'
                elif 'tipper' in raw_value:
                    value = 'tipper'
                if 'all' in raw_value:
                    value = 'all'
                self.ui.comboBox_mode.setCurrentIndex(int(d[value]))
                update_counter += 1
            except:
                pass

        if 'model_name' in parameters:
            try:
                value = str(parameters['model_name'])
                self.ui.lineEdit_modelname.setText(value)
                update_counter += 1
            except:
                pass
        if 'mu_start' in parameters:
            try:
                value = float(parameters['mu_start'])
                self.ui.doubleSpinBox_lagrange.setValue(value)
                update_counter += 1
            except:
                pass
        if 'no_iteration' in parameters:
            try:
                value = int(float(parameters['no_iteration']))
                self.ui.spinBox_iterationstep.setValue(value)
                update_counter += 1
            except:
                pass
        if 'no_layers' in parameters:
            try:
                value = int(float(parameters['no_layers']))
                self.ui.spinBox_no_layers.setValue(value)
                update_counter += 1
            except:
                pass

        if 'phase_errorfloor' in parameters:
            try:
                value = (parameters['phase_errorfloor'])
                if len(value) == 0 or value.lower().strip() == 'none':
                    self.ui.checkBox_phase_error.setCheckState(0)
                    self.ui.doubleSpinBox_phase_error.setValue(15)
                else:
                    value = float(value)
                    self.ui.checkBox_phase_error.setCheckState(2)
                    self.ui.doubleSpinBox_phase_error.setValue(value)
                update_counter += 1
            except:
                self.ui.checkBox_phase_error.setCheckState(0)
        if 'reached_misfit' in parameters:
            try:
                value = int(float(parameters['reached_misfit']))
                if value == 0:
                    self.ui.checkBox_misfitreached.setCheckState(0)
                else:
                    self.ui.checkBox_misfitreached.setCheckState(2)
                update_counter += 1
            except:
                self.ui.checkBox_misfitreached.setCheckState(0)
                
        if 'rho_errorfloor' in parameters:
            try:
                value = (parameters['rho_errorfloor'])
                if len(value) == 0 or value.lower().strip() == 'none':
                    self.ui.checkBox_rho_error.setCheckState(0)
                    self.ui.doubleSpinBox_rho_error.setValue(10)
                else:
                    value = float(value)
                    self.ui.checkBox_rho_error.setCheckState(2)
                    self.ui.doubleSpinBox_rho_error.setValue(value)
                update_counter += 1
            except:
                self.ui.checkBox_rho_error.setCheckState(0)
        if 'tipper_errorfloor' in parameters:
            try:
                value = (parameters['tipper_errorfloor'])
                if len(value) == 0 or value.lower().strip() == 'none':
                    self.ui.checkBox_tipper_error.setCheckState(0)
                    self.ui.doubleSpinBox_tipper_error.setValue(10)
                else:
                    value = float(value)
                    self.ui.checkBox_tipper_error.setCheckState(2)
                    self.ui.doubleSpinBox_tipper_error.setValue(value)
                update_counter += 1
            except:
                self.ui.checkBox_tipper_error.setCheckState(0)
        if 'strike' in parameters:
            try:
                value = (parameters['strike'])
                if len(value) == 0 or value.lower().strip() == 'none':
                    self.ui.checkBox_strike.setCheckState(0)
                    self.ui.doubleSpinBox_strike.setValue(0)
                else:
                    value = float(value)
                    self.ui.checkBox_strike.setCheckState(2)
                    self.ui.doubleSpinBox_strike.setValue(value)
                update_counter += 1
            except:
                self.ui.checkBox_strike.setCheckState(2)
                self.ui.doubleSpinBox_strike.setValue(0)
        else:
            self.ui.checkBox_strike.setCheckState(0)
            self.ui.doubleSpinBox_strike.setValue(0)

        
        if 'target_rms' in parameters:
            try:
                value = float(parameters['target_rms'])
                self.ui.doubleSpinBox_rms.setValue(value)
                update_counter += 1
            except:
                pass        
        if 'model_depth' in parameters:
            try:
                value = float(parameters['model_depth'])
                self.ui.doubleSpinBox_model_depth.setValue(value)
                update_counter += 1
            except:
                pass        
        if 'wd' in parameters:
            try:
                value = str(parameters['wd'])
                self.ui.lineEdit_browse_wd.setText(value)
                update_counter += 1
            except:
                pass


        
        messagetext = ''
        messagetext += "<P><FONT COLOR='#000000'>Configuration file: "\
                    "{0}  </FONT></P> \n".format(old_cfg_filename)
        messagetext += "<P><b><FONT COLOR='#008080'>Read in {0} parameters"\
                                "</FONT></b></P>".format(update_counter)

        QtGui.QMessageBox.about(self, "Update parameters from file", messagetext )


        

    #check input values for consistency/existence
    def check_input(self):

        #checking files, paths and filenames (including write permission in WD)
    
        messagetext = "<br><P><b><FONT COLOR='#008080'>Parameters are valid !  </FONT></b></P><br>"

        invalid_flag = 0 

        #working directory
        wd_mess = ''
        wd_text = ''
        if str(self.ui.lineEdit_browse_wd.text()):
            wd_text = str(self.ui.lineEdit_browse_wd.text())


        if (wd_text.strip() == '') or (not op.isdir(op.abspath(op.realpath(wd_text)))):
            wd_mess += 'Working directory not existing <br>'
            invalid_flag +=1

        elif not os.access(self.ui.wd, os.W_OK):
            wd_mess += 'Working directory not writable <br>'
            invalid_flag +=1


        #EDI files folder
        edis_mess = ''
        edis_text = ''
        if str(self.ui.lineEdit_browse_edi.text()):
            edis_text = str(self.ui.lineEdit_browse_edi.text())
        if (edis_text.strip() == '') or (not op.isdir(op.abspath(op.realpath(edis_text)))):
            edis_mess += 'EDI files directory not existing <br>'
            invalid_flag +=1

        #OCCAM executable
        occ_exe_mess = ''
        occ_exe_text = ''
        if str(self.ui.lineEdit_browse_occam.text()):
            occ_exe_text = str(self.ui.lineEdit_browse_occam.text())
        if (occ_exe_text.strip() == '') or (not op.isfile(op.abspath(op.realpath(op.join(self.ui.wd,occ_exe_text))))):
            occ_exe_mess += 'Occam executable not existing <br>'
            invalid_flag +=1

        #makemodel_mess = ''
        #makemodel_text = ''
        #if str(self.ui.lineEdit_browse_makemodel.text()):
            #makemodel_text = str(self.ui.lineEdit_browse_makemodel.text())
        #if (makemodel_text.strip() == '') or (not op.isfile(op.abspath(op.realpath(op.join(self.ui.wd,makemodel_text))))):
            #makemodel_mess += 'Make2DModel executable not existing <br>'
            #invalid_flag +=1


        #stations list file
        stations_mess = ''
        stations_text = ''
        if str(self.ui.lineEdit_browse_stationfile .text()):
            stations_text = str(self.ui.lineEdit_browse_stationfile.text())
        if (stations_text.strip() == '') or (not op.isfile(op.abspath(op.realpath(op.join(self.ui.wd,stations_text))))):
            if self.ui.checkBox_usestationlist.checkState():
                stations_mess += 'Stations file not existing <br>'
                invalid_flag +=1

        #OCCAM startup file 
        iteration_mess=''
        
        startup_text=''
        if self.ui.checkBox_useiterationfile.checkState():

            if str(self.ui.lineEdit_browse_iterationfile.text()):
                startup_text = str(self.ui.lineEdit_browse_iterationfile.text())
            else:
                startup_text = None
            if (startup_text.strip() == '') or  (startup_text is None) or \
                    (not op.isfile(op.realpath(op.join(self.ui.wd,startup_text)))):
                iteration_mess += 'startup file not existing <br>'
                invalid_flag += 1
        
        #OCCAM data file (existing file)
        datafile_mess = ''
        datafile_text = ''
        if str(self.ui.lineEdit_browse_datafile.text()):
            datafile_text = str(self.ui.lineEdit_browse_datafile.text())
        if (datafile_text.strip() == '') or (not op.isfile(op.abspath(op.realpath(op.join(self.ui.wd,datafile_text))))):
            if self.ui.checkBox_usedatafile.checkState() and (not self._load_datafile):
                datafile_mess += 'Datafile not existing <br>'
                invalid_flag +=1
        
        

        #OCCAM data file (new file name)
        datafilename_mess = ''
        datafilename_text = ''        
        datafilename_text_raw = str(self.ui.lineEdit_datafilename.text()).strip()
        #check for emtpty slot:
        if not datafilename_text_raw:
            datafilename_mess += 'No datafile given <br>'
            invalid_flag +=1

        elif len(datafilename_text_raw.split()) != 1 :
            datafilename_mess += 'White space found in datafile name <br>'
            invalid_flag +=1
            
        else:
            #check for additional strange characters:
            import re
            m2 = re.match("^[a-zA-Z_\.]+[a-zA-Z0-9_\.]*$", datafilename_text_raw, re.M)
            if not m2 and (not self._load_datafile):
                datafilename_mess += 'Wrong format of datafile (a-z,0-9,_) <br>'
                invalid_flag +=1
                
        #OCCAM model name
        modelname_mess = ''
        modelname_text = ''        
        modelname_text_raw = str(self.ui.lineEdit_modelname.text()).strip()
        #print modelname_text_raw
        if not modelname_text_raw:
            modelname_mess += 'No model name given <br>'
            invalid_flag +=1
        # elif len(modelname_text_raw.split()) != 1 :
        #     modelname_mess += 'White space found in model name <br>'
        #     invalid_flag +=1
        else:
            #check model name for other strange characters:
            import re
            m1 = re.match("\s", modelname_text_raw, re.M)
            m2 = re.match("^[a-zA-Z_]+[a-zA-Z0-9_\.]*$", modelname_text_raw, re.M)
            if m1 or (not m2):
                modelname_mess += 'Wrong format of model name (a-z,0-9,_) <br>'
            

        #set up error message, if problems occurred
        if invalid_flag != 0 :
            conc_mess = wd_mess+occ_exe_mess+iteration_mess+edis_mess+stations_mess+datafile_mess+datafilename_mess+modelname_mess
            messagetext = "<P><b><FONT COLOR='#800000'>Error: %i parameters are invalid !  </FONT></b></P><br> %s"%(invalid_flag,conc_mess)

            
        #define/how message box
        QtGui.QMessageBox.about(self, "Parameter check", messagetext)
        if invalid_flag > 0:
            return 1
        else:
            return 0
        

    def setup_parameter_dict(self):
        """
        Define dictionary, which contains all parameters for OCCAM

        """

        
        D = {}
        D['wd']             = op.abspath( op.realpath( str( self.ui.lineEdit_browse_wd.text() ) ) )
        D['edi_dir']          = op.abspath( op.realpath( str( self.ui.lineEdit_browse_edi.text()) ) )
        D['occam_exe']        = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_occam.text()) ) ) )
        
        D['use_iterationfile'] = self.ui.checkBox_useiterationfile.checkState()
        if D['use_iterationfile']:
            D['iterationfile']  = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_iterationfile.text()) ) ) )
        else:
            D['iterationfile']  = None

        D['use_stationfile'] = self.ui.checkBox_usestationlist.checkState()
        if D['use_stationfile']:
            D['stationlistfile'] = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_stationfile.text()) ) ) )
        else:
            D['stationlistfile'] = None

        D['use_olddatafile'] = self.ui.checkBox_usedatafile.checkState()
        if D['use_olddatafile']:
            D['olddatafile']      = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_datafile.text()) ) ) )
        else:
            D['olddatafile'] = None

        D['datafile']     = str(self.ui.lineEdit_datafilename.text())

        D['modelname']        = str(self.ui.lineEdit_modelname.text())
        D['no_iteration']     = self.ui.spinBox_iterationstep.value()
        D['mu_start']         = self.ui.doubleSpinBox_lagrange.value()
        D['debug_level']      = int(float(self.ui.comboBox_debuglevel.currentText()))
        D['check_misfit_reached']   = self.ui.checkBox_misfitreached.checkState()
        if D['check_misfit_reached'] :
            D['misfit_reached'] = 1
        else:
            D['misfit_reached'] = 0

        D['set_max_no_frequencies'] = self.ui.checkBox_max_no_frequencies.checkState()
        if D['set_max_no_frequencies']:
            D['max_no_frequencies'] = self.ui.spinBox_max_no_frequencies.value()
        else:
            D['max_no_frequencies'] = None

        D['set_min_frequency'] = self.ui.checkBox_min_frequency.checkState()
        if D['set_min_frequency']:
            D['min_frequency'] = self.ui.doubleSpinBox_min_frequency.value()
        else:
            D['min_frequency'] = None

        D['set_max_frequency'] = self.ui.checkBox_max_frequency.checkState()
        if D['set_max_frequency'] :
            D['max_frequency'] = self.ui.doubleSpinBox_max_frequency.value() 
        else:
            D['max_frequency'] = None
        

        D['check_usestationfile'] = self.ui.checkBox_usestationlist.checkState()

        D['check_useolddatafile'] = self.ui.checkBox_usedatafile.checkState()
        
        D['check_useiterationfile'] = self.ui.checkBox_useiterationfile.checkState()
        

        D['mode']             = str(self.ui.comboBox_mode.currentText()).lower()

        D['edi_type']         = str(self.ui.comboBox_edi_type.currentText()).lower()
        if D['edi_type'].startswith('rho'):
            D['edi_type'] = 'resphase'
        #D['freqsteps']        = self.ui.spinBox_freq_steps.value()
        
        D['strikeisknown']    = self.ui.checkBox_strike.checkState()

        if D['strikeisknown'] :
            D['strike']           = float(self.ui.doubleSpinBox_strike.value())
        else:
            D['strike']           = None

        
        D['block_merge_threshold']   = self.ui.doubleSpinBox_mergethreshold.value()
        D['max_no_iterations']     = self.ui.spinBox_max_no_iterations.value()
        D['target_rms']       = self.ui.doubleSpinBox_rms.value()
        D['no_layers']         = self.ui.spinBox_no_layers.value()
        D['firstlayer_thickness']       = self.ui.spinBox_firstlayer.value()
        D['max_blockwidth']   = self.ui.spinBox_maxblockwidth.value()
        D['halfspace_resistivity']  = self.ui.doubleSpinBox_rhostart.value()
        D['model_depth']  = self.ui.doubleSpinBox_model_depth.value()

        D['set_rho_errorfloor']   = self.ui.checkBox_rho_error.checkState()
        if D['set_rho_errorfloor']:
            D['rho_errorfloor']          = self.ui.doubleSpinBox_rho_error.value()
        else:
            D['rho_errorfloor'] = None

        D['set_phase_errorfloor']   = self.ui.checkBox_phase_error.checkState()
        if  D['set_phase_errorfloor']: 
            D['phase_errorfloor']          = self.ui.doubleSpinBox_phase_error.value()
        else:
            D['phase_errorfloor']          = None

        D['set_tipper_errorfloor'] = self.ui.checkBox_tipper_error.checkState()
        if D['set_tipper_errorfloor'] :
            D['tipper_errorfloor']        = self.ui.doubleSpinBox_tipper_error.value()
        else:
            D['tipper_errorfloor']        = None
        
        self.parameters = D



    def _setup_startupfile(self):
        """
        Set up OCCAM startup file, called 'startup' by default
        """
        self.parameters['startupfile'] = 'startup'

        #use old startup file, if box is checked
        if self.parameters['check_usestartupfile']:
            self.parameters['startupfile'] = self.parameters['startupfn']

        
 

    def build_inputfiles(self):
        """
        Set up collection of required input files files for OCCAM.
        
        - data
        - startup
        - model
        - mesh
        """
        returnvalue = 0

        #1. Build OCCAM data file, or use existing one (if field is checked in GUI) 
        

        D = self.parameters

        #print D['check_usedatafile']
        #print D['olddatafile']

        olddatafile = D['olddatafile']

        if olddatafile is not None:
            datafile = op.abspath(op.join(D['wd'],olddatafile))
            messagetext = ''        

            returnvalue = 0 
            try:
                data_object = MTo2.Data()
                data_object.readfile(datafile)

                D['strike'] = data_object.strike
                D['azimuth'] = data_object.azimuth
                D['stationlocations']  = data_object.stationlocations
                
                messagetext += "<P><FONT COLOR='#000000'>Working directory: "\
                    "{0}  </FONT></P> \n".format(data_object.wd)

                messagetext += "<P><b><FONT COLOR='#008080'>Read old data file:</FONT></b></P><br>{0}".format(datafile)
            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error: Cannot read old data file: {0}  </FONT></b></P> ".format(datafile)
                returnvalue = 1

                
            QtGui.QMessageBox.about(self, "Data file generation", messagetext )
            if returnvalue == 1:
                return

        else:
        
            outfilename = D['datafile']
            edidirectory= D['edi_dir']


            def make_stationlist(listfilename):
                """
                Read in stations from file.

                """
                FH = file(listfilename,'r')
                raw_string = FH.read()
                FH.close()
                raw_list1 = raw_string.strip().split()
                raw_list2 = []
                for i in raw_list1:
                    if len(i.split(',')) == 1:
                        raw_list2.append(i)
                    else:
                        for j in i.split(','):
                            raw_list2.append(j)
                return raw_list2

            #define internal station list 
            stationlist  = None
            if D['use_stationfile']:
                stationlist = make_stationlist(D['stationlistfile'])

            D['stationlist'] = stationlist
            
            #make data file  -------------------------------------------
            returnvalue = 0 
            messagetext = ''
            try:
                setup_object = MTo2.Setup(**D)
            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
                "generate setup object - check input parameters!  </FONT></b></P> \n"
                QtGui.QMessageBox.about(self, "Input files generation", messagetext )
                return 1

            try:
                edi_dir = D['edi_dir']
                setup_object.read_edifiles(edi_dir)
                datafile = D['datafile']
                setup_object.datafile = datafile
                try:
                    setup_object.write_datafile()
                except:
                    raise
                datafilename = setup_object.datafile
                self.parameters['stationlocations']  = setup_object.Data.stationlocations
                messagetext += "<P><FONT COLOR='#000000'>Working directory: "\
                        "{0}  </FONT></P> \n".format(setup_object.wd)

                messagetext += "<P><b><FONT COLOR='#008080'>"\
                "Data file: {0}  </FONT></b></P> \n".format(op.split(setup_object.datafile)[1])
            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
                "write data file: {0}  </FONT></b></P> \n".format(setup_object.datafile)
                returnvalue = 1

                    
            QtGui.QMessageBox.about(self, "Data file generation", messagetext )

        #---------------
        #2. other input files:

        D=self.parameters
                        
        #datafile = D['datafile'] 
        messagetext = ''
        # try:
        #     #1. make startup file 
        #     self._setup_startupfile()
        # except:
        #     messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not generate startup file!  </FONT></b></P> \n"

        if olddatafile is not None:       
            try:
                setup_object = MTo2.Setup(**D)
            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
                "generate setup object - check input parameters!  </FONT></b></P> \n"
                QtGui.QMessageBox.about(self, "Input files generation", messagetext )

                return 1

        # edi_dir = D['edi_dir']
        # if not D['check_usedatafile']:
        #     try:
        #         setup_object.read_edifiles(edi_dir)
        #         setup_object.datafile = D['datafilename']
        #         setup_object.write_datafile()
        #         messagetext += "<P><b><FONT COLOR='#008080'>Wrote "\
        #         "data file: {0}  </FONT></b></P> \n".format(setup_object.datafile)
        #     except:
        #         messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
        #         "write data file: {0}  </FONT></b></P> \n".format(setup_object.datafile)
            
        #     QtGui.QMessageBox.about(self, "Data file generation", messagetext )

        try:
            setup_object.setup_mesh_and_model()
        except:
            messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
            "set up mesh and model!  </FONT></b></P> \n"
            returnvalue = 1
        messagetext += "<P><FONT COLOR='#000000'>Working directory: "\
            "{0}  </FONT></P> \n".format(setup_object.wd)
        try:
            setup_object.write_meshfile()
            messagetext += "<P><b><FONT COLOR='#008080'>"\
            "Mesh file: {0}  </FONT></b></P> \n".format(setup_object.meshfile)
        except:
            messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
            "write mesh file: {0}  </FONT></b></P> \n".format(setup_object.meshfile)
            returnvalue = 1

        try:
            setup_object.write_inmodelfile()
            messagetext += "<P><b><FONT COLOR='#008080'>"\
            "Inmodel file: {0}  </FONT></b></P> \n".format(setup_object.inmodelfile)
        except:
            messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
            "write inmodel file: {0}  </FONT></b></P> \n".format(setup_object.inmodelfile)
            returnvalue = 1

        if D['check_useiterationfile']:
            try:
                setup_object.startupfile = D['iterationfile']
                base,short_fn = op.split(setup_object.startupfile)
                if base != setup_object.wd:
                    new_startupfile = op.abspath(op.join(setup_object.wd,short_fn))
                    shutil.copy(setup_object.startupfile,new_startupfile )
                setup_object.startupfile = short_fn
                messagetext += "<P><b><FONT COLOR='#008080'>Using old "\
                "iteration file for startup: {0}  </FONT></b></P> \n".format(setup_object.startupfile)
                D['startupfile'] =  D['iterationfile']
            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error: Could not "\
                "find old iteration file: {0}  </FONT></b></P> \nUsing default 'startup' instead ".format(D['iterationfile'])
                D['startupfile'] ='startup'
                returnvalue = 1


        else:
            try:
                setup_object.write_startupfile()
                messagetext += "<P><b><FONT COLOR='#008080'>"\
                "Startup file: {0}  </FONT></b></P> \n".format(setup_object.startupfile)
                D['startupfile'] = setup_object.startupfile

            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
                "write startup file: {0}  </FONT></b></P> \n".format(setup_object.startupfile)
                D['startupfile'] ='startup'
                returnvalue = 1

        try:
            setup_object.write_configfile()
            messagetext += "<P><b><FONT COLOR='#008080'>"\
            "Configuration file: {0}  </FONT></b></P> \n".format(op.split(setup_object.configfile)[1])
        except:
            messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
            "write configuration file: {0}  </FONT></b></P> \n".format(op.split(setup_object.configfile)[1])
            returnvalue = 1


        QtGui.QMessageBox.about(self, "Input files generation", messagetext )

        return returnvalue


    def generate_inputfiles(self):
        
        if self.check_input() == 1:
            return
        self.setup_parameter_dict()

        if self.build_inputfiles() ==1 :
            return 
        #print 'startupfiles done'
    

#===============================================================================
# problem start
#===============================================================================


    def run_occam(self):
        """Method for calling the external occam code.

        - check input values
        - disable the 'run' button
        - set the executable string and its arguments
        - start occam
        - show message box including 'cancel' button

        - message box should disappear, if occam terminates (how so ever)
        - stdout and stderr of occam should be redirected into files
        - no general blocking of gui

        """

        #import necessary modules
        import threading
        import subprocess
        import Queue


        #check input values and set up startup file
        #self._setup_startupfile()
        
        #deactivate start button to avoid multiple calls
        # re- activation does not work - WHY ??????? - has to stay commented
        #self.ui.pushButton_runoccam.setEnabled(False)

        #define executable and argument(s)
        
        self.generate_inputfiles()
        
        exename    = self.parameters['occam_exe']
        modname    = self.parameters['modelname']
        startfile  = self.parameters['startupfile']

        exec_file      = exename
        exec_args      = [exec_file,startfile,modname]

        #define log and err files
        logfilename    = op.abspath(op.realpath(op.join(self.ui.wd,'occam2d_log.log')))
        errfilename    = op.abspath(op.realpath(op.join(self.ui.wd,'occam2d_err.log')))


        #open log and err files - potentially redirecting sys variables
        
        #save_stdout = sys.stdout
        #save_stderr = sys.stderr
        outfile = open(logfilename,"w")
        errfile = open(errfilename,'w')
        #sys.stdout = outfile
        #sys.stderr = errfile


        #define queue
        Q = Queue.Queue() 


        #===================================================================
        #===================================================================
        # needs 1 extra worker thread for OCCAM
        #===================================================================
        
        
        #1. try: QProcess ....
        
        #occam_process=QtCore.QProcess()#.start(startstring)
        #bb = occam_process.startDetached(startstring)



        #startstring = "%s %s %s"%(exename,startfile,modname)
        #startstring = "nohup %s %s %s >& %s &"%(exename,modname,startfile,logfilename)

        #2. try: os.system call
        #os.system(startstring)


        #3. try pexpect (no windows version!!!)
        
        #import pexpect
        #child1 = pexpect.spawn(startstring)
        #child1.logfile=outfile
        #child1.logfile=sys.stdout
        #child1.expect(pexpect.EOF)


        #4. try: thread setup

        OccamExitFlag = 0
        
        class OccamThread(QtCore.QThread):
            
            def __init__(self, exec_file, exec_args,logfile,errfile):
                #super(OccamThread, self).__init__(self)

                QtCore.QThread.__init__(self)
                self.exiting    =  False
                
                self.exec_file = exec_file
                self.startup = exec_args[1]
                self.outname = exec_args[2]
                self.logfile = logfile
                self.errfile = errfile
                self.running = True

            def __del__(self):

                self.running  =  False
                self.wait()
                

            def run(self):

                #while not self.exiting :
                #self.occam_process =  QtCore.QProcess()
                    #self.occam_process.connect(self.occam_finished)
                #self.occam_process.startDetached(exec_file,exec_args)


#subprocess works, but stderr and stdout redirection fail !!!!!!
#if 'wait' or 'communicate' are included, the gui freezes!!

                 subprocess.Popen([self.exec_file,self.startup,self.outname])#,                stdout=subprocess.PIPE,stderr=subprocess.PIPE)


#                while self.occam_process.poll():
                    #std_out,std_err = self.occam_process.communicate()
                    #print std_out,std_err
                    #self.logfile.write(std_out)
                    #self.logfile.flush()
                    #self.errfile.write(std_err)
                    #self.errfile.flush(


                    
            def occam_finished(self,exitCode, exitStatus):
                
                #self.parent.ui.pushButton_runoccam.setEnabled(True)
                print exitCode, exitStatus
                        
                        
        #starting the worker thread
        ot = OccamThread(exec_file, exec_args,outfile,errfile)
        ot.daemon = True
        try:
            ot.start()
        except (KeyboardInterrupt, SystemExit):
            ot.running = False
            #ot.kill()
            ot.join()
            ot.logfile.close()
            ot.errfile.close()
            print 'OCCAM terminated by keyboard interruption'
            sys.exit()
            


        #starting a message box for start/running of process:
        
        MB_start = QtGui.QMessageBox(self)
        MB_start.setStandardButtons(QtGui.QMessageBox.Cancel)
        MB_start_ret_value = 0
        MB_start.setText('OCCAM is running' )

        try:
            MB_start_ret_value = MB_start.exec_()
        except (KeyboardInterrupt, SystemExit):
            MB_start_ret_value = QtGui.QMessageBox.Cancel
            

        #define message box for end of process:
        MB_end = QtGui.QMessageBox(self)
        MB_end.setStandardButtons(QtGui.QMessageBox.Cancel)
        MB_end_ret_value = 0

        
        
        def _close_occam_mb(self):
            MB_start_ret_value = QtGui.QMessageBox.Cancel

    
        #try to connect a potential termination signal from the OccamThread to the trivial function above
        #Fails, because Thread ot cannot be connected.....WHY???????????

        #self.connect(ot, QtCore.SIGNAL('finished()'), _close_occam_mb)
        #self.connect(ot, QtCore.SIGNAL('terminated()'), _close_occam_mb)
            
            
        #this works - press CANCEL button and the OCCAM run stops
        if MB_start_ret_value == QtGui.QMessageBox.Cancel:
            #print ot.__dict__.items()
            ot.running = False
            #ot.kill()
            ot.join()
            ot.logfile.close()
            ot.errfile.close()
            print 'OCCAM terminated by user input'
            
            endtext = 'OCCAM terminated by user!'

            MB_end.about(self, "OCCAM 2D", endtext )
            #re-activate button
            #does not work - WHY ???????
            
            #self.ui.pushButton_runoccam.setEnabled(False)
            return

        
        
        while ot.isAlive():
            pass
        
        endtext = 'finished'

        MB_end.setText("OCCAM 2D")
        MB_end.setInformativeText(endtext)


        #closing start/running message box

        #send KILL to MB_start widget ....:
        #print MB_start
        #print MB_start.__dir__


        
        MB_end_ret_value = MB_end.exec_()

        
        if (not outfile.closed):
            outfile.close()
        if not errfile.closed():
            errfile.close()


        #closing stdout and err files
        
        #outfil.close()
        #errfil.close()
        #sys.stdout = save_stdout
        #sys.stderr = save_stderr
        
        self.ui.pushButton_runoccam.setEnabled(False)
        
#===============================================================================
# problem end
#===============================================================================


    def add_entry(self):
        self.ui.textEdit_file.selectAll()
        self.ui.textEdit_file.cut()
        self.ui.textEdit.append("")
        self.ui.textEdit.paste()


    def print_entry(self):
        entrylist = str(self.ui.textEdit.toPlainText())
        print entrylist


        
def main():
    #problem with stdout and err files - trying global redirection:

    # stdout and stderr are saved
    #save_stdout = sys.stdout
    #save_stderr = sys.stderr

    #outfile = open("occam2d_gui_log.log","w")
    #errfile = open('occam2d_gui_err.log','w')
    #sys.stdout = outfile
    #sys.stderr = errfile

    app               = QtGui.QApplication(sys.argv)

    occamgui_instance = OccamGui()
    occamgui_instance.show()

    sys.exit(app.exec_())
                
    #outfile.close()
    #errfile.close()
    #sys.stderr = save_stderr
    #sys.stdout = save_stdout
    
if __name__ == "__main__":
    main()                
