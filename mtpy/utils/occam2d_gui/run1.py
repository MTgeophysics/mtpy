#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys,os
from PyQt4 import QtCore, QtGui
import gui4
reload(gui4)

from gui4 import Ui_occamgui2D as Ui_Form

import MTpy.core.OCCAMTools as OCCAMTools

import os.path as op




class MyForm(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self._load_datafile = None
        self.parameters = {}
        self.ui.wd = op.abspath(op.realpath('.'))
        
        QtCore.QObject.connect(self.ui.button_browse_wd, QtCore.SIGNAL("clicked()"),  lambda: self.set_path_in_browsefield(self.ui.lineEdit_browse_wd))
        
        QtCore.QObject.connect(self.ui.button_browse_edis, QtCore.SIGNAL("clicked()"),  lambda: self.set_path_in_browsefield(self.ui.lineEdit_browse_edis))

        QtCore.QObject.connect(self.ui.pushButton_loadstations, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_stations))

        QtCore.QObject.connect(self.ui.button_browse_occam, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_occam))
        
        #QtCore.QObject.connect(self.ui.button_browse_makemodel, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_makemodel))

        QtCore.QObject.connect(self.ui.pushButton_loaddatafile, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_datafile))

        QtCore.QObject.connect(self.ui.pushButton_loadstartupfile, QtCore.SIGNAL("clicked()"),  lambda: self.set_filename_in_browsefield(self.ui.lineEdit_browse_startupfile))
        
        QtCore.QObject.connect(self.ui.pushButton_loaddatafile, QtCore.SIGNAL("clicked()"),self.set_data_filename)
        
        QtCore.QObject.connect(self.ui.pushButton_checkparameter, QtCore.SIGNAL("clicked()"),  self.check_input)

        QtCore.QObject.connect(self.ui.pushButton_checkparameter, QtCore.SIGNAL("clicked()"),  self.setup_parameter_dict)
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
        

    def set_path_in_browsefield(self, browsefield):
        dirname = QtGui.QFileDialog.getExistingDirectory(self, 'Open Directory', '.')
        browsefield.setText(dirname)

    def set_filename_in_browsefield(self, browsefield):
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Locate File', '.')
        browsefield.setText(filename)


    def check_input(self):

        #checking files, paths and filenames (including write permission in WD)
    
        messagetext = "<br><P><b><FONT COLOR='#008080'>Parameters are valid !  </FONT></b></P><br>"

        invalid_flag = 0 

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


        edis_mess = ''
        edis_text = ''
        if str(self.ui.lineEdit_browse_edis.text()):
            edis_text = str(self.ui.lineEdit_browse_edis.text())
        if (edis_text.strip() == '') or (not op.isdir(op.abspath(op.realpath(edis_text)))):
            edis_mess += 'EDI files directory not existing <br>'
            invalid_flag +=1

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

        stations_mess = ''
        stations_text = ''
        if str(self.ui.lineEdit_browse_stations.text()):
            stations_text = str(self.ui.lineEdit_browse_stations.text())
        if (stations_text.strip() == '') or (not op.isfile(op.abspath(op.realpath(op.join(self.ui.wd,stations_text))))):
            if self.ui.checkBox_usestationlist.checkState():
                stations_mess += 'Stations file not existing <br>'
                invalid_flag +=1

        startup_mess=''
        startup_text=''
        if str(self.ui.lineEdit_browse_startupfile.text()):
            startup_text = str(self.ui.lineEdit_browse_startupfile.text())
        if (startup_text.strip() == '') or  (not op.isfile(op.realpath(op.join(self.ui.wd,startup_text)))):
            if self.ui.checkBox_usestartupfile.checkState():
                startup_mess += 'Startup file not existing <br>'
                invalid_flag += 1
                
         
        datafile_mess = ''
        datafile_text = ''
        if str(self.ui.lineEdit_browse_datafile.text()):
            datafile_text = str(self.ui.lineEdit_browse_datafile.text())
        if (datafile_text.strip() == '') or (not op.isfile(op.abspath(op.realpath(op.join(self.ui.wd,datafile_text))))):
            if self.ui.checkBox_usedatafile.checkState() and (not self._load_datafile):
                datafile_mess += 'Datafile not existing <br>'
                invalid_flag +=1
        
        

    
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
            


        if invalid_flag != 0 :
            conc_mess = wd_mess+occ_exe_mess+startup_mess+edis_mess+stations_mess+datafile_mess+datafilename_mess+modelname_mess
            messagetext = "<P><b><FONT COLOR='#800000'>Error: %i parameters are invalid !  </FONT></b></P><br> %s"%(invalid_flag,conc_mess)

            
        
        QtGui.QMessageBox.about(self, "Parameter check", messagetext)
        

    def setup_parameter_dict(self):
        #set paths
        D = {}
        D['wdir']             = op.abspath( op.realpath( str( self.ui.lineEdit_browse_wd.text() ) ) )
        D['edis_dir']         = op.abspath( op.realpath( str( self.ui.lineEdit_browse_edis.text()) ) )
        D['occam_exe']        = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_occam.text()) ) ) )
        D['startupfile']      = op.abspath( op.realpath( op.join(self.ui.wd, str(self.ui.lineEdit_browse_startupfile.text()) ) ) )
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
        D=self.parameters

        print D['check_usedatafile']
        print D['olddatafile']
        checkfilename = D['olddatafile']
        print op.isfile(checkfilename)
        
        if D['check_usedatafile']:
            checkfilename = D['olddatafile']
            if op.isfile(checkfilename):
                messagetext = "<P><b><FONT COLOR='#008080'>Existing data file loaded successfully:</FONT></b></P><br>%s"%checkfilename
                D['datafilename'] = D['olddatafile']
                D['datafile']     = D['olddatafile']
            else:
                messagetext = "<P><b><FONT COLOR='#800000'>Error: %i No old data file found!  </FONT></b></P> "
                    
            QtGui.QMessageBox.about(self, "Data file", messagetext )
                    
            return


        outfilename = D['datafilename']
        edidirectory= D['edis_dir']

        #define stationlist-----------------------------------------
        def make_stationlist(listfilename):
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

        stationlist=None
        if D['check_usestations']:
            stationlist = make_stationlist(D['stationlist_file'])

        #define error floors-----------------------------------------
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
        
        checkfilename = OCCAMTools.make2DdataFile(edidirectory,\
                                    savepath=outfilename,\
                                    mmode= mode4modelling,\
                                    stationlst=stationlist,\
                                    title=None,\
                                    thetar=float(D['strike']),\
                                    resxyerr=resxy,\
                                    resyxerr=resyx,\
                                    phasexyerr=phasexy,\
                                    phaseyxerr=phaseyx,\
                                    ss=3*' ',\
                                    fmt='%2.6f',\
                                    freqstep=int(float(D['freqsteps'])),\
                                    plotyn='n',\
                                    lineori=lineorientation,\
                                    tippererr=tippererror  \
                                    )
        
        if op.isfile(checkfilename):
            messagetext = "<P><b><FONT COLOR='#008080'>Data file generation successful:</FONT></b></P><br>%s"%checkfilename
            D['datafile'] = checkfilename
        else:
            messagetext = "<P><b><FONT COLOR='#800000'>Error: %i No data file generated!  </FONT></b></P> "
            
        QtGui.QMessageBox.about(self, "Data file", messagetext )

    def _setup_startupfile(self):
        self.parameters['startupfn'] = None    

        if self.parameters['check_usestartupfile']:
            self.parameters['startupfn'] = self.parameters['startupfile']

        
 

    def build_startupfiles(self):
        D=self.parameters
        
        datafile = D['datafile'] 

        self._setup_startupfile()
            
        m_fn,i_fn,s_fn = OCCAMTools.makeModel(datafile,\
                                                niter=int(float(D['n_iterations'])),\
                                                targetrms=float(D['target_rms']),\
                                                nlayers=int(float(D['n_layers'])),\
                                                nlperdec=int(float(D['decade_layers'])),\
                                                z1layer=int(float(D['thickness1'])),\
                                                bwidth=int(float(D['max_blockwidth'])),\
                                                trigger=float(D['mergethreshold']),\
                                                #savepath=self.ui.wd,\
                                                rhostart=float(D['rho0']),\
                                                #occampath=D['makemodel_exe']\
                                                cwd=self.ui.wd,\
                                                #makemodelexe=D['makemodel_exe'],
                                                modelname=D['modelname'],
                                                use_existing_startup=D['check_usestartupfile'],\
                                                existing_startup_file=D['startupfn']\
                                                )

        if op.isfile(m_fn) and op.isfile(i_fn) and op.isfile(s_fn):
            if D['check_usestartupfile']:
                messagetext = "<P><b><FONT COLOR='#008080'>Old startup  read in:</FONT></b></P><br>%s"%(s_fn)
                QtGui.QMessageBox.about(self, "Startup file", messagetext )
                messagetext = "<P><b><FONT COLOR='#008080'>Input files generated:</FONT></b></P><br>%s<br>%s"%(m_fn,i_fn)
            
            else:
                messagetext = "<P><b><FONT COLOR='#008080'>Input files generated:</FONT></b></P><br>%s<br>%s<br>%s"%(m_fn,i_fn,s_fn)
        else:
            messagetext = "<P><b><FONT COLOR='#800000'>Error: %i No startup files generated!  </FONT></b></P> "
        
        QtGui.QMessageBox.about(self, "Startup files generation", messagetext )

        return


    def generate_inputfiles(self):
        self.build_datafile()
        #print'datafile done'
        self.build_startupfiles()
        #print 'startupfiles done'
    

        

    def run_occam(self):
        import subprocess
        exename = self.parameters['occam_exe']
        modname = self.parameters['modelname']
        startfile = self.parameters['startupfile']


        #occam_process=QtCore.QProcess()#.start(startstring)

        #bb = occam_process.startDetached(startstring)

        exec_file = exename
        exec_args = [exec_file,startfile,modname]
        logfilename    = op.abspath(op.realpath(op.join(self.ui.wd,'occam2d_log.log')))
        errfilename    = op.abspath(op.realpath(op.join(self.ui.wd,'occam2d_err.log')))

        #startstring = "%s %s %s"%(exename,startfile,modname)
        #startstring = "nohup %s %s %s >& %s &"%(exename,modname,startfile,logfilename)
        

        #save_stdout = sys.stdout
        #save_stderr = sys.stderr
        #outfile = open(logfilename,"w")
        #errfile = open(errfilename,'w')
        #sys.stdout = outfile
        #sys.stderr = errfile
        
        
        #os.system(startstring)
        #import pexpect
        #child1 = pexpect.spawn(startstring)
        #child1.logfile=outfile
        #child1.logfile=sys.stdout
        #child1.expect(pexpect.EOF)
        
        occam_process = subprocess.Popen(exec_args)#, stdout=logfil,stderr=errfil)
        #while child1.isalive():
        #    outfile.flush()
        
        
        MB_start = QtGui.QMessageBox(self)
        MB_start.setStandardButtons(QtGui.QMessageBox.Cancel)
        MB_start_ret_value = 0
        MB_start.setText('OCCAM is running' )
        MB_start_ret_value = MB_start.exec_()



        MB_end = QtGui.QMessageBox(self)
        
        endtext = 'OCCAM finished!'
        
        if MB_start_ret_value == QtGui.QMessageBox.Cancel:
            #child1.close(force=True)
            occam_process.terminate()
            #child1.expect(pexpect.EOF)
            #outfile.close()
            
            endtext = 'OCCAM terminated by user!'

        #outfil.close()
        #errfil.close()
        #sys.stdout = save_stdout
        #sys.stderr = save_stderr
        
        
        #if not outfile.closed:
            #child1.expect(pexpect.EOF)
            #outfile.close()
            
        

        
        MB_end.about(self, "OCCAM 2D", endtext )
        
        

        
        #occam_process.communicate()

        
        
        #subprocess.os.system("%s %s %s"%(exename,startfile,modname))


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
       




            
if __name__ == "__main__":


    # stdout and stderr are saved
    save_stdout = sys.stdout
    save_stderr = sys.stderr

    outfile = open("occam2d_gui_log.log","w")
    errfile = open('occam2d_gui_err.log','w')
    #sys.stdout = outfile
    #sys.stderr = errfile


    app = QtGui.QApplication(sys.argv)
    occam_vals = MyForm()
    occam_vals.show()

    sys.exit(app.exec_())
                


    outfile.close()
    errfile.close()
    #sys.stderr = save_stderr
    #sys.stdout = save_stdout
    

                
