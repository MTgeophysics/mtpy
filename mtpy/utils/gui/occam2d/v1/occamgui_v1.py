#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys,os
from PyQt4 import QtCore, QtGui
import gui4
reload(gui4)

from gui4 import Ui_occamgui2D as Occam_UI_form

import os.path as op

import mtpy.core.edi as MTedi
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
        QtCore.QObject.connect(self.ui.button_browse_edis, QtCore.SIGNAL("clicked()"),  lambda: self.set_path_in_browsefield(self.ui.lineEdit_browse_edis))
        QtCore.QObject.connect(self.ui.pushButton_loadstations, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_stations))
        QtCore.QObject.connect(self.ui.button_browse_occam, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_occam))
        #QtCore.QObject.connect(self.ui.button_browse_makemodel, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_makemodel))
        QtCore.QObject.connect(self.ui.pushButton_loaddatafile, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_datafile))
        QtCore.QObject.connect(self.ui.pushButton_loadstartupfile, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_startupfile))       
        QtCore.QObject.connect(self.ui.pushButton_loaddatafile, QtCore.SIGNAL("clicked()"),self.set_data_filename)
        QtCore.QObject.connect(self.ui.pushButton_checkparameter, QtCore.SIGNAL("clicked()"),  self.check_input)
        #QtCore.QObject.connect(self.ui.pushButton_checkparameter, QtCore.SIGNAL("clicked()"),  self.setup_parameter_dict)
        QtCore.QObject.connect(self.ui.pushButton_generateinputfile, QtCore.SIGNAL("clicked()"),  self.setup_parameter_dict)
        QtCore.QObject.connect(self.ui.pushButton_runoccam, QtCore.SIGNAL("clicked()"),  self.setup_parameter_dict)
        QtCore.QObject.connect(self.ui.pushButton_generateinputfile, QtCore.SIGNAL("clicked()"),  self.check_input)
        QtCore.QObject.connect(self.ui.pushButton_runoccam, QtCore.SIGNAL("clicked()"),  self.check_input)
        QtCore.QObject.connect(self.ui.pushButton_generateinputfile, QtCore.SIGNAL("clicked()"),  self.generate_inputfiles)
        #QtCore.QObject.connect(self.ui.pushButton_loadold_meshinmodel, QtCore.SIGNAL("clicked()"),  self.loadold_meshinmodel)
        QtCore.QObject.connect(self.ui.pushButton_runoccam, QtCore.SIGNAL("clicked()"),  self.run_occam)
        QtCore.QObject.connect(self.ui.pushButton_quit, QtCore.SIGNAL("clicked()"), QtCore.QCoreApplication.instance().quit)

    
        
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
        if str(self.ui.lineEdit_browse_edis.text()):
            edis_text = str(self.ui.lineEdit_browse_edis.text())
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
        if str(self.ui.lineEdit_browse_stations.text()):
            stations_text = str(self.ui.lineEdit_browse_stations.text())
        if (stations_text.strip() == '') or (not op.isfile(op.abspath(op.realpath(op.join(self.ui.wd,stations_text))))):
            if self.ui.checkBox_usestationlist.checkState():
                stations_mess += 'Stations file not existing <br>'
                invalid_flag +=1

        #OCCAM startup file 
        startup_mess=''
        startup_text=''
        if str(self.ui.lineEdit_browse_startupfile.text()):
            startup_text = str(self.ui.lineEdit_browse_startupfile.text())
        if (startup_text.strip() == '') or  (not op.isfile(op.realpath(op.join(self.ui.wd,startup_text)))):
            if self.ui.checkBox_usestartupfile.checkState():
                startup_mess += 'Startup file not existing <br>'
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
        elif len(modelname_text_raw.split()) != 1 :
            modelname_mess += 'White space found in model name <br>'
            invalid_flag +=1
        else:
            #check model name for other strange characters:
            import re
            m1 = re.match("\s", modelname_text_raw, re.M)
            m2 = re.match("^[a-zA-Z_]+[a-zA-Z0-9_\.]*$", modelname_text_raw, re.M)
            if m1 or (not m2):
                modelname_mess += 'Wrong format of model name (a-z,0-9,_) <br>'
            

        #set up error message, if problems occurred
        if invalid_flag != 0 :
            conc_mess = wd_mess+occ_exe_mess+startup_mess+edis_mess+stations_mess+datafile_mess+datafilename_mess+modelname_mess
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
        D['wdir']             = op.abspath( op.realpath( str( self.ui.lineEdit_browse_wd.text() ) ) )
        D['edi_dir']         = op.abspath( op.realpath( str( self.ui.lineEdit_browse_edis.text()) ) )
        D['occam_exe']        = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_occam.text()) ) ) )
        D['startupfn']        = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_startupfile.text()) ) ) )
        D['stationlist_file'] = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_stations.text()) ) ) )
        D['olddatafile']      = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_datafile.text()) ) ) )
        D['datafilename']     = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_datafilename.text()) ) ) )

        D['modelname']        = str(self.ui.lineEdit_modelname.text())
        
        D['check_usestations']= self.ui.checkBox_usestationlist.checkState()
        D['check_usedatafile']= self.ui.checkBox_usedatafile.checkState()
        D['check_usestartupfile']= self.ui.checkBox_usestartupfile.checkState()
        

        D['mode']             = self.ui.comboBox_mode.currentIndex()
        D['freqsteps']        = self.ui.spinBox_freq_steps.value()
        D['strike']           = self.ui.doubleSpinBox_strike.value()
        D['orientation']      = self.ui.comboBox_alignment.currentIndex()
        D['mergethreshold']   = self.ui.doubleSpinBox_mergethreshold.value()
        D['n_iterations']     = self.ui.spinBox_iterations.value()
        D['target_rms']       = self.ui.doubleSpinBox_rms.value()
        D['n_layers']         = self.ui.spinBox_layers.value()
        D['thickness1']       = self.ui.spinBox_firstlayer.value()
        D['max_blockwidth']   = self.ui.spinBox_maxblockwidth.value()
        D['decade_layers']    = self.ui.spinBox_layersperdecade.value()
        D['rho0']             = self.ui.doubleSpinBox_rhostart.value()

        D['useDataError_ResXY']   = self.ui.checkBox_usedataerror_resXY.checkState()
        D['Error_ResXY']          = self.ui.doubleSpinBox_errorvalue_resXY.value()
        D['useDataError_ResYX']   = self.ui.checkBox_usedataerror_resYX.checkState()
        D['Error_ResYX']          = self.ui.doubleSpinBox_errorvalue_resYX.value()
        D['useDataError_PhaseXY'] = self.ui.checkBox_usedataerror_phaseXY.checkState()
        D['Error_PhaseXY']        = self.ui.doubleSpinBox_errorvalue_phaseXY.value()
        D['useDataError_PhaseYX'] = self.ui.checkBox_usedataerror_phaseYX.checkState()
        D['Error_PhaseYX']        = self.ui.doubleSpinBox_errorvalue_phaseYX.value()

        D['check_includeTipper']  = self.ui.checkBox_usetipper.checkState()
        D['Error_Tipper']         = self.ui.doubleSpinBox_errorvalue_tipper.value()
        

        self.parameters = D


    def build_datafile(self):
        """
        Build OCCAM data file, or use existing one (if field is checked in GUI) 
        """

        D = self.parameters

        #print D['check_usedatafile']
        #print D['olddatafile']

        checkfilename = D['olddatafile']
        #print op.isfile(checkfilename)


        #Try to load old data file, if box is checked
        # if D['check_usedatafile']:
        #     checkfilename = D['olddatafile']
        #     if op.isfile(checkfilename):
        #         messagetext = "<P><b><FONT COLOR='#008080'>Existing data file loaded successfully:</FONT></b></P><br>%s"%checkfilename
        #         D['datafilename'] = D['olddatafile']
        #         D['datafile']     = D['olddatafile']
        #     else:
        #         messagetext = "<P><b><FONT COLOR='#800000'>Error:  No old data file found!  </FONT></b></P> "
                    
        #     QtGui.QMessageBox.about(self, "Data file", messagetext )
                    
        #     return

        # #Otherwise set up new data filename
        
        outfilename = D['datafilename']
        edidirectory= D['edi_dir']


        def make_stationlist(listfilename):
            """
            Read in stations from file.

            """
            FH=file(listfilename,'r')
            raw_string=FH.read()
            FH.close()
            raw_list1 = raw_string.strip().split()
            raw_list2 = []
            for i in raw_list1:
                if len(i.split(','))==1:
                    raw_list2.append(i)
                else:
                    for j in i.split(','):
                        raw_list2.append(j)
            return raw_list2

        #define internal station list 
        stationlist=None
        if D['check_usestations']:
            stationlist = make_stationlist(D['stationlist_file'])

        #read in error floors-----------------------------------------
        if int(float(D['useDataError_ResXY'])): 
            resxy='data'
        else:
            resxy=float(D['Error_ResXY'] )    


        if int(float(D['useDataError_ResYX'])):
            resyx='data'
        else:
            resyx=float(D['Error_ResYX'] )
    

        if int(float(D['useDataError_PhaseXY'])):
            phasexy='data'
        else:
            phasexy=float(D['Error_PhaseXY'])

        if int(float(D['useDataError_PhaseYX'])):
            phaseyx='data'
        else:
            phaseyx=float(D['Error_PhaseYX'])


        #TIPPER--------------------------------------------------
        tippererror=None
        if D['check_includeTipper']:
            tippererror = float(D['Error_Tipper'])

        #define modes used ------------------------------------------
        _modedict = {'0':'both' , '1':'TM', '2':'TE' }
        mode4modelling = _modedict[str(D['mode'])]
       
        #define orientation -----------------------------------------
        _oridict = {'0':'ns' , '1':'ew' }
        lineorientation = _oridict[str(D['orientation'])]
        
        #make data file  -------------------------------------------
        returnvalue = 0 
        messagetext = ''
        if not D['check_usedatafile']:
            try:
                setup_object = MTo2.Setup(**D)
                edi_dir = D['edi_dir']
                setup_object.read_edifiles(edi_dir)
                setup_object.datafile = D['datafilename']
                setup_object.write_datafile()
                datafilename = setup_object.datafile
                self.parameters['stationlocations']  = setup_object.stationlocations
                messagetext += "<P><b><FONT COLOR='#008080'>Wrote "\
                "data file: {0}  </FONT></b></P> \n".format(setup_object.datafile)
            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
                "write data file: {0}  </FONT></b></P> \n".format(setup_object.datafile)
                returnvalue = 1

        else:
            try:
                data_object = MTo2.Data()
                data_object.readfile(D['datafilename'])
                self.parameters['stationlocations']  = setup_object.stationlocations
                messagetext = "<P><b><FONT COLOR='#008080'>Read old data file:</FONT></b></P><br>{0}".format(D['datafilename'])
            except:
                messagetext = "<P><b><FONT COLOR='#800000'>Error: Cannot read old data file: {0}  </FONT></b></P> ".format(D['datafilename'])
                returnvalue = 1

                
        QtGui.QMessageBox.about(self, "Data file generation", messagetext )

        return returnvalue


            # checkfilename = OCCAMTools.make2DdataFile(edidirectory,\
            #                             savepath=outfilename,\
            #                             mmode= mode4modelling,\
            #                             stationlst=stationlist,\
            #                             title=None,\
            #                             thetar=float(D['strike']),\
            #                             resxyerr=resxy,\
            #                             resyxerr=resyx,\
            #                             phasexyerr=phasexy,\
            #                             phaseyxerr=phaseyx,\
            #                             ss=3*' ',\
            #                             fmt='%2.6f',\
            #                             freqstep=int(float(D['freqsteps'])),\
            #                             plotyn='n',\
            #                             lineori=lineorientation,\
            #                             tippererr=tippererror  \
            #                             )

            #check, if generation was successfull and show message box



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

        - startup
        - model
        - mesh
        """
        returnvalue = 0
        D=self.parameters
                        
        #datafile = D['datafile'] 
        messagetext = ''
        # try:
        #     #1. make startup file 
        #     self._setup_startupfile()
        # except:
        #     messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not generate startup file!  </FONT></b></P> \n"

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

        try:
            setup_object.write_meshfile()
            messagetext += "<P><b><FONT COLOR='#008080'>Wrote "\
            "mesh file: {0}  </FONT></b></P> \n".format(setup_object.meshfile)
        except:
            messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
            "write mesh file: {0}  </FONT></b></P> \n".format(setup_object.meshfile)
            returnvalue = 1

        try:
            setup_object.write_inmodelfile()
            messagetext += "<P><b><FONT COLOR='#008080'>Wrote "\
            " inmodel file: {0}  </FONT></b></P> \n".format(setup_object.inmodelfile)
        except:
            messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
            "write inmodel file: {0}  </FONT></b></P> \n".format(setup_object.inmodelfile)
            returnvalue = 1

        if D['check_usestartupfile']:
            try:
                setup_object.startupfile = D['startupfile']
                messagetext += "<P><b><FONT COLOR='#008080'>Using old "\
                "startup file: {0}  </FONT></b></P> \n".format(setup_object.startupfile)
            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error: Could not "\
                "find old startup file: {0}  </FONT></b></P> \n".format(setup_object.startupfile)
                returnvalue = 1

        else:
            try:
                setup_object.write_startupfile()
                messagetext += "<P><b><FONT COLOR='#008080'>Wrote "\
                " startup file: {0}  </FONT></b></P> \n".format(setup_object.startupfile)

            except:
                messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not "\
                "write startup file: {0}  </FONT></b></P> \n".format(setup_object.startupfile)
                returnvalue = 1



            # #2. make model file
            # m_fn,i_fn,s_fn = OCCAMTools.makeModel(datafile,\
            #                                         niter=int(float(D['n_iterations'])),\
            #                                         targetrms=float(D['target_rms']),\
            #                                         nlayers=int(float(D['n_layers'])),\
            #                                         nlperdec=int(float(D['decade_layers'])),\
            #                                         z1layer=int(float(D['thickness1'])),\
            #                                         bwidth=int(float(D['max_blockwidth'])),\
            #                                         trigger=float(D['mergethreshold']),\
            #                                         #savepath=self.ui.wd,\
            #                                         rhostart=float(D['rho0']),\
            #                                         #occampath=D['makemodel_exe']\
            #                                         cwd=self.ui.wd,\
            #                                         #makemodelexe=D['makemodel_exe'],
            #                                         modelname=D['modelname'],
            #                                         use_existing_startup=D['check_usestartupfile'],\
            #                                         existing_startup_file=D['startupfn']\
            #                                         )
        # except:
        #     messagetext += "<P><b><FONT COLOR='#800000'>Error:  Could not generate model/mesh!  </FONT></b></P> \n"

        # if D['check_usestartupfile']:
        #     messagetext = "<P><b><FONT COLOR='#008080'>Old startup  read in:</FONT></b></P><br>%s"%(s_fn)
        #         QtGui.QMessageBox.about(self, "Startup file", messagetext )
                
        #         messagetext = "<P><b><FONT COLOR='#008080'>Input files generated:</FONT></b></P><br>%s<br>%s"%(m_fn,i_fn)
        #         D['startupfile'] = existing_startup_file
            
        #     else:
        #         messagetext += "<P><b><FONT COLOR='#008080'>Input files generated:</FONT></b></P><br>%s<br>%s<br>%s"%(m_fn,i_fn,s_fn)
        # else:
        #     messagetext += "<P><b><FONT COLOR='#800000'>Error:  No input files generated!  </FONT></b></P> "
        
        QtGui.QMessageBox.about(self, "Input files generation", messagetext )

        return returnvalue


    def generate_inputfiles(self):
        if self.build_datafile() == 1 :
            return

        #print'datafile done'
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
        self._setup_startupfile()
        
        #deactivate start button to avoid multiple calls
        # re- activation does not work - WHY ??????? - has to stay commented
        #self.ui.pushButton_runoccam.setEnabled(False)

        #define executable and argument(s)
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
        ot.setDaemon(True)
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



    def loadold_meshinmodel(self):
        pass



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
