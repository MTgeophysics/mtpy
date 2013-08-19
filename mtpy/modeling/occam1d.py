# -*- coding: utf-8 -*-
"""
==================
Occam1D
==================

    * Wrapper class to interact with Occam1D written by Kerry Keys at Scripps 
      adapted from the method of Constable et al., [1987].
      
    * This class only deals with the MT functionality of the Fortran code, so 
      it can make the input files for computing the 1D MT response of an input
      model and or data.  It can also read the output and plot them in a 
      useful way.
      
    * Note that when you run the inversion code, the convergence is quite 
      quick, within the first few iterations, so have a look at the L2 cure 
      to decide which iteration to plot, otherwise if you look at iterations
      long after convergence the models will be unreliable.
      
      
      
    * Key, K., 2009, 1D inversion of multicomponent, multi-frequency marine
      CSEM data: Methodology and synthetic studies for resolving thin 
      resistive layers: Geophysics, 74, F9–F20.

    * The original paper describing the Occam's inversion approach is:

    * Constable, S. C., R. L. Parker, and C. G. Constable, 1987, 
      Occam’s inversion –– A practical algorithm for generating smooth 
      models from electromagnetic sounding data, Geophysics, 52 (03), 289–300.
      


@author: JP
"""

import numpy as np
import scipy as sp
import os
import subprocess
import shutil
import fnmatch
from operator import itemgetter
import time
import matplotlib.colorbar as mcb
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import mtpy.core.edi as mtedi
import mtpy.modeling.winglinktools as wlt
import matplotlib.pyplot as plt
import mtpy.utils.latlongutmconversion as utm2ll



occamdict={'1':'resxy','2':'phasexy','3':'realtip','4':'imagtip','5':'resyx',
           '6':'phaseyx'}

class Occam1D:
    """
    ==============================================
    This class will deal with everything occam 1D
    =============================================
    """
    
    def __init__(self, savepath=None):
        
        self.savepath = None             
        self.modelfn = None
        self.inputfn = None
        self.datafn_te = None
        self.datafn_tm = None
        self.iternum = 0
  
    def make1DdataFile(self, station, edipath=None, savepath=None,
                       polarization='both', reserr='data', phaseerr='data',
                       string_fmt='%+.6e', ss=3*' ', thetar=0):
        """
        make1Ddatafile will write a data file for Occam1D
    
        Arguments:
        ---------    
            **station** : string
                          the station name and path if edipath=None
            
            **edipath** : string
                          path to the edi files to be written into a data file,
                          useful for multile data files
                      
            **savepath** : string
                           path to save the file, if None set to dirname of 
                           station if edipath = None.  Otherwise set to 
                           dirname of edipath.
            
            **thetar** : float
                         rotation angle to rotate Z. Clockwise positive and N=0
                         *default* = 0
            
            **polarization** : [ 'both' | 'TE' | 'TM' | 'det']
                              polarization to model can be (*default*='both'):
                
                                - 'both' for TE and TM as separate files
                                - 'TE' for just TE mode
                                - 'TM' for just TM mode
                                - 'det' for the determinant of Z.
                                .. note:: 
                                    
                                    if polarization = 'det' two files 
                                    will be created stationDet_TE.dat and 
                                    stationDet_TM.dat.  These files both use
                                    the determinant however the code for the 
                                    Occam input is different so you can test 
                                    the difference between inverting the 
                                    determinant as TE or TM because there is
                                    no option in Occam for the Det.
                                
                            
            **reserr** : float
                        errorbar for resistivity values.  Can be set to (
                        *default* = 'data'): 
                
                        - 'data' for errorbars from the data
                        - percent number ex. 10 for ten percent
                    
            **phaseerr** : float
                          errorbar for phase values.  Can be set to (
                          *default* = 'data'):
                
                            - 'data' for errorbars from the data
                            - percent number ex. 10 for ten percent
                        
            **string_fmt** : format of the values written to the file. 
                      *default* = %+.6e
            
            **ss** : spacing between values in file.  *default* = ' '*3
            
        Returns:
        --------
            **Occam1D.datafn_te** : full path to data file for TE mode
                
            **Occam1D.datafn_tm** : full path to data file for TM mode
                
        :Example: ::
            
            >>> old = occam.Occam1D()
            >>> old.make1DdataFile('MT01',edipath=r"/home/Line1",
            >>>                    savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                    mode='TE')
            >>> Wrote Data File: /home/Occam1D/Line1/Inv1_TE/MT01TE.dat 
        """    
    
        if os.path.dirname(station)=='':
            if edipath==None:
                raise IOError('Need to input a path for the file.')
            else:
                #find the edifile
                for fn in os.listdir(edipath):
                    if fn.lower().find(station.lower())>=0:
                        edifile=os.path.join(edipath,fn)
        else:
            edifile=station
        
        self.station=os.path.basename(edifile)[:-4]
        
        #raise an error if can't find the edifile        
        if edifile==None:
            raise NameError('No edifile exists, check path and station name')
    
        #read in edifile    
        impz=Z.Z(edifile)    
        
        #make sure the savepath exists, if not create it
        if savepath==None:
            if not self.savepath:
                savepath=os.path.dirname(edifile)
                if not os.path.exists(savepath):
                    os.mkdir(savepath)
                savepath=self.savepath
            else:
                savepath=self.savepath
        elif os.path.basename(savepath).find('.')>0:
            savepath=os.path.dirname(savepath)
            if not os.path.exists(savepath):
                os.mkdir(os.path.dirname(savepath))
            self.savepath=savepath
        else:
            if not os.path.exists(savepath):
                os.mkdir(savepath)
            self.savepath=savepath
        
        #load the edifile and get resistivity and phase
        rp=impz.getResPhase(thetar=thetar)
        freq=impz.frequency
        nf=len(freq)
        returnfn=[]
        
        pdict={'TE':['xy'],
               'TM':['yx'],
               'both':['xy','yx'],
               'det':['det'],
               'all':['xy','yx','det']}
        if polarization=='both' or polarization=='det':
            for pol in ['xy','yx']:
                if pol=='xy':
                    if polarization=='det':
                        dfilesave=os.path.join(self.savepath,
                                               impz.station+'Det_TE.dat')
                    else:
                        dfilesave=os.path.join(self.savepath,
                                               impz.station+'TE.dat')
                elif pol=='yx':
                    if polarization=='det':
                        dfilesave=os.path.join(self.savepath,
                                               impz.station+'Det_TM.dat')
                    else:
                        dfilesave=os.path.join(self.savepath,
                                               impz.station+'TM.dat')
    
                datafid=open(dfilesave,'w')
    
                datafid.write('Format:  EMData_1.1 \n')
                datafid.write('!Polarization:'+ss+pol+'\n')
    
                #needs a transmitter to work so put in a dummy one
                datafid.write('# Transmitters: 1\n')
                datafid.write('0 0 0 0 0 \n')
                
                #write frequencies
                datafid.write('# Frequencies:'+ss+str(nf)+'\n')       
                for ff in freq:
                    datafid.write(ss+'%.6f' % ff+'\n')
                
                #needs a receiver to work so put in a dummy one
                datafid.write('# Receivers: 1 \n')
                datafid.write('0 0 0 0 0 0 \n')
                
                #write data
                datafid.write('# Data:'+2*ss+str(2*nf)+'\n')
                datafid.write('!'+2*ss+'Type'+2*ss+'Freq#'+2*ss+'Tx#'+2*ss+
                             'Rx#'+ 2*ss+'Data'+2*ss+'Std_Error'+'\n')
                              
                #put the yx phase component in the first quadrant as prescribed
                if pol=='yx':
                        rp.phaseyx=rp.phaseyx+180
                        #check if there are any negative phases
                        negphase=np.where(rp.phaseyx>180)
                        if len(negphase)>0:
                            rp.phaseyx[negphase[0]]=rp.phaseyx\
                                                            [negphase[0]]-360
                        
                #write the resistivity and phase components
                for ii in range(nf):
                    #------------write resistivity components------------------
                    if reserr=='data':
                        if pol=='xy':
                            if polarization=='det':
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resdet[ii]+2*ss+
                                              string_fmt % rp.resdeterr[ii]+'\n')
                            else:
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resxy[ii]+2*ss+
                                              string_fmt % rp.resxyerr[ii]+'\n')
                        elif pol=='yx':
                            if polarization=='det':
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resdet[ii]+2*ss+
                                              string_fmt % rp.resdeterr[ii]+'\n')
                            else:
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resyx[ii]+2*ss+
                                              string_fmt % rp.resyxerr[ii]+'\n')
                    #-----------if percent error is given--------------------                          
                    else:
                        if pol=='xy':
                            if polarization=='det':
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resdet[ii]+2*ss+
                                              string_fmt % (rp.resdet[ii]*reserr/100.)+
                                              '\n')
                            else:
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resxy[ii]+2*ss+
                                              string_fmt % (rp.resxy[ii]*reserr/100.)+
                                              '\n')
                        elif pol=='yx':
                            if polarization=='det':
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resdet[ii]+2*ss+
                                              string_fmt % (rp.resdet[ii]*reserr/100.)+
                                              '\n')
                            else:
                                datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.resyx[ii]+2*ss+
                                              string_fmt % (rp.resyx[ii]*reserr/100.)+
                                              '\n')
                    
                    #---------------write phase components--------------------
                    if phaseerr=='data':
                        if pol=='xy':
                            if polarization=='det':
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phasedet[ii]+2*ss+
                                              string_fmt % rp.phasedeterr[ii]+'\n')
                            else:
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phasexy[ii]+2*ss+
                                              string_fmt % rp.phasexyerr[ii]+'\n')
                        if pol=='yx':
                            if polarization=='det':
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phasedet[ii]+2*ss+
                                              string_fmt % rp.phasedeterr[ii]+'\n')
                            else:
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phasexy[ii]+2*ss+
                                              string_fmt % rp.phasexyerr[ii]+'\n')
                    #-----------if percent error is given--------------------                          
                    else:
                        if pol=='xy':
                            if polarization=='det':
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phasedet[ii]+2*ss+
                                              string_fmt % (phaseerr/100.*(180/np.pi))+
                                              '\n')
                            else:
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phasexy[ii]+2*ss+
                                              string_fmt % (phaseerr/100.*(180/np.pi))+
                                              '\n')
                        if pol=='yx':
                            if polarization=='det':
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phasedet[ii]+2*ss+
                                              string_fmt % (phaseerr/100.*(180/np.pi))+
                                              '\n')
                            else:
                                datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+
                                              2*ss+'0'+2*ss+'0'+2*ss+
                                              string_fmt % rp.phaseyx[ii]+2*ss+
                                              string_fmt % (phaseerr/100.*(180/np.pi))+
                                              '\n')
                datafid.write('\n')
                datafid.close()
                print 'Wrote Data File: ',dfilesave
                returnfn.append(dfilesave)
            self.datafn_te=returnfn[0]
            self.datafn_tm=returnfn[1]
        else:
            if polarization=='TE':
                pol='xy'
                dfilesave=os.path.join(savepath,impz.station+'TE.dat')
                self.datafn_te=dfilesave
            elif polarization=='TM':
                pol='yx'
                dfilesave=os.path.join(savepath,impz.station+'TM.dat')
                self.datafn_te=dfilesave

            #open file to write to
            datafid=open(dfilesave,'w')
            datafid.write('Format:  EMData_1.1 \n')
            datafid.write('!Polarization:'+ss+pol+'\n')
    
            #needs a transmitter to work so put in a dummy one
            datafid.write('# Transmitters: 1\n')
            datafid.write('0 0 0 0 0 \n')
            
            #write frequencies
            datafid.write('# Frequencies:'+ss+str(nf)+'\n')       
            for ff in freq:
                datafid.write(ss+'%.6f' % ff+'\n')
            
            #needs a receiver to work so put in a dummy one
            datafid.write('# Receivers: 1 \n')
            datafid.write('0 0 0 0 0 0 \n')
            
            #write header line
            datafid.write('# Data:'+2*ss+str(2*nf)+'\n')
            datafid.write('!'+2*ss+'Type'+2*ss+'Freq#'+2*ss+'Tx#'+2*ss+'Rx#'+
                          2*ss+'Data'+2*ss+'Std_Error'+'\n')
                          
            #put the yx phase component in the first quadrant as prescribed
            if pol=='yx':
                    rp.phaseyx=rp.phaseyx+180
                    #check if there are any negative phases
                    negphase=np.where(rp.phaseyx>180)
                    if len(negphase)>0:
                        rp.phaseyx[negphase[0]]=rp.phaseyx\
                                                        [negphase[0]]-360
                    
            #write the resistivity and phase components
            for ii in range(nf):
                #write resistivity components
                if reserr=='data':
                    if pol=='xy':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.resxy[ii]+2*ss+
                                      string_fmt % rp.resxyerr[ii]+'\n')
                    elif pol=='yx':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.resyx[ii]+2*ss+
                                      string_fmt % rp.resyxerr[ii]+'\n')
                    elif pol=='det':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.resyx[ii]+2*ss+
                                      string_fmt % rp.resyxerr[ii]+'\n')
                else:
                    if pol=='xy':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.resxy[ii]+2*ss+
                                      string_fmt % (rp.resxy[ii]*reserr/100.)+'\n')
                    elif pol=='yx':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.resyx[ii]+2*ss+
                                      string_fmt % (rp.resyx[ii]*reserr/100.)+'\n')
                
                #write phase components
                if phaseerr=='data':
                    if pol=='xy':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.phasexy[ii]+2*ss+
                                      string_fmt % rp.phasexyerr[ii]+'\n')
                    if pol=='yx':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.phaseyx[ii]+2*ss+
                                      string_fmt % rp.phaseyxerr[ii]+'\n')
                else:
                    if pol=='xy':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.phasexy[ii]+2*ss+
                                      string_fmt % (phaseerr/100.*(180/np.pi))+'\n')
                    if pol=='yx':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+string_fmt % rp.phaseyx[ii]+2*ss+
                                      string_fmt % (phaseerr/100.*(180/np.pi))+'\n')
            datafid.close()
            print 'Wrote Data File: ',dfilesave

    def make1DModelFile(self,savepath=None,nlayers=100,bottomlayer=10000,
                        basestep=10,z1layer=10,airlayerheight=10000):
        """
        Makes a 1D model file for Occam1D.  
        
        Arguments:
        ----------
        
            **savepath** :path to save file to, if just path saved as 
                          savepath\model.mod, if None defaults to dirpath
                          
            **nlayers** : number of layers
            
            **bottomlayer** : depth of bottom layer in meters
            
            **basestep** : numerical base of logarithmic depth step 10 or 2 or 
                          1 for linear
                          
            **z1layer** : depth of first layer in meters
            
            **airlayerheight** : height of air layers in meters
              
        Returns:
        --------
        
            **Occam1D.modelfn** = full path to model file
            
        ..Note: This needs to be redone.
        
        :Example: ::
            
            >>> old = occam.Occam1D()
            >>> old.make1DModelFile(savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                     nlayers=50,bottomlayer=10000,z1layer=50)
            >>> Wrote Model file: /home/Occam1D/Line1/Inv1_TE/Model1D 
        """
        
        
        ss='   '
        #if the savepath was not entered test to see if there is one
        if savepath==None:
            if not self.savepath:
                raise IOError('No savepath found.  Please input one.')
            self.modelfn=os.path.join(self.savepath,'Model1D') 
        
        #if the save path was entered as just a path
        elif os.path.basename(savepath).find('.')==-1:
            if not os.path.exists(savepath):
                os.mkdir(savepath)
            self.savepath=savepath
            self.modelfn=os.path.join(self.savepath,'Model1D')
         
        #if the save path was entered as a file with full path
        else:
            self.modelfn=savepath
        
        #---------need to refine this-------------------- 
        
        layers=np.logspace(np.log10(z1layer),np.log10(bottomlayer),num=nlayers)      
        
        #make the model file
        modfid=open(self.modelfn,'w')
        modfid.write('Format: Resistivity1DMod_1.0'+'\n')
        modfid.write('#LAYERS:    '+str(nlayers+2)+'\n')
        modfid.write('!Set free values to -1 or ? \n')
        modfid.write('!penalize between 1 and 0,'+
                     '0 allowing jump between layers and 1 smooth. \n' )
        modfid.write('!preference is the assumed resistivity on linear scale. \n')
        modfid.write('!pref_penalty needs to be put if preference is not 0 [0,1]. \n')
        modfid.write('! top_depth'+ss+'resistivity'+ss+'penalty'+ss+'preference'+ss+
                     'pref_penalty \n')
        modfid.write(ss+'-10000'+ss+'1d12'+ss+'0'+ss+'0'+ss+'0'+ss+'!air layer \n')
        modfid.write(ss+'0'+ss+'-1'+ss+'0'+ss+'0'+ss+'0'+ss+'!first ground layer \n')
        for ll in layers:
            modfid.write(ss+'{0:.2f}'.format(ll)+ss+'-1'+ss+'1'+ss+'0'+ss+'0'+
                         '\n')
        
        modfid.close()
        
        print 'Wrote Model file: ',self.modelfn
        

    def make1DInputFile(self,savepath=None,imode='TE',roughtype=1,
                        maxiter=20,targetrms=1.0,rhostart=100,
                        description='1dInv',lagrange=5.0,roughness=1.0E7,
                        debuglevel=1,iteration=0,misfit=100.0):
        """
        Make a 1D input file for Occam 1D
        
        Arguments:
        ---------
            **savepath** : full path to save input file to, if just path then 
                           saved as savepath/input
                           
            **modelfile** : full path to model file, if None then assumed to be in 
                            savepath/model.mod
                            
            **datafile** : full path to data file, if None then assumed to be 
                            in savepath/TE.dat or TM.dat
                            
            **roughtype** : roughness type. *default* = 0
            
            **maxiter** : maximum number of iterations. *default* = 20 
            
            **targetrms** : target rms value. *default* = 1.0
            
            **rhostart** : starting resistivity value on linear scale. 
                            *default* = 100
            
            **description** : description of the inversion. 
            
            **lagrange** : starting Lagrange multiplier for smoothness.
                           *default* = 5
            
            **roughness** : starting roughness value. *default* = 1E7
            
            **debuglevel** : something to do with how Fortran debuggs the code
                             Almost always leave at *default* = 1
                    
            **iteration** : the starting iteration number, handy if the
                            starting model is from a previous run.
                            *default* = 0
            
            **misfit** : starting misfit value. *default* = 100
                             
        Returns:
        --------
            **Occam1D.inputfn** : full path to input file. 
            
        :Example: ::
            
            >>> old = occam.Occam1D()
            >>> old.make1DdataFile('MT01',edipath=r"/home/Line1",
            >>>                    savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                    mode='TE')
            >>> Wrote Data File: /home/Occam1D/Line1/Inv1_TE/MT01TE.dat
            >>>            
            >>> old.make1DModelFile(savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                     nlayers=50,bottomlayer=10000,z1layer=50)
            >>> Wrote Model file: /home/Occam1D/Line1/Inv1_TE/Model1D
            >>>
            >>> old.make1DInputFile(rhostart=10,targetrms=1.5,maxiter=15)
            >>> Wrote Input File: /home/Occam1D/Line1/Inv1_TE/Input1D  
        """
        
        
        ss='   '
        
        #make input data file name
        #if no savepath is input, test if there is already one
        if savepath==None:
            if not self.savepath:
                raise IOError('No savepath.  Please input one.')
            self.inputfn=os.path.join(self.savepath,'Input1D') 
        
        #if the savepath was input as just a path
        elif os.path.basename(savepath).find('.')==-1:
            if not os.path.exists(savepath):
                os.mkdir(savepath)
            self.inputfn=os.path.join(savepath,'Input1D')
        
        #if the savepath was input as a full path to file
        else:
            self.inputfn=savepath
            
        if not self.modelfn:
            if self.savepath:
                self.modelfn=os.path.join(self.savepath,'Model1D')
            else:
                raise IOError('No savepath.  Please input one.')
                
        #try to get data file name 
        if imode=='TE':
            if not self.datafn_te:
                if self.savepath:
                    try:
                        self.datafn_te=[os.path.join(self.savepath,dd) 
                            for dd in os.listdir(self.savepath) 
                            if dd.find('TE.dat')>0][0]
                    except IndexError:
                        raise IOError('No TE data file found. Please input one.')
                else:
                    raise IOError('No savepth found. Please input one.')
            else:
                pass
        if imode=='TM':
            if not self.datafn_tm:
                if self.savepath:
                    try:
                        self.datafn_tm=[os.path.join(self.savepath,dd) 
                            for dd in os.listdir(self.savepath) 
                            if dd.find('TM.dat')>0][0]
                    except IndexError:
                        raise IOError('No TM data file found. Please input one.') 
            else:
                pass

        #read in the model and get number of parameters
        self.read1DModelFile()
        paramcount=self.mdict['nparam']        
        
        #write input file
        infid=open(self.inputfn,'w')
        infid.write('Format:             OCCAMITER_FLEX      ! Flexible format \n')
        infid.write('Description:        '+description+'     !For your own notes. \n')
        infid.write('Model File:         '+self.modelfn+'       \n')
        if imode=='TE':
            infid.write('Data File:          '+self.datafn_te+'        \n')                                                                     
        if imode=='TM':
            infid.write('Data File:          '+self.datafn_tm+'        \n')                                                                     
        infid.write('Date/Time:          '+time.ctime()+'\n')         
        infid.write('Max Iter:           '+str(maxiter)+'\n')
        infid.write('Target Misfit:      '+str(targetrms)+'\n')
        infid.write('Roughness Type:     '+str(roughtype)+'\n')
        infid.write('!Model Bounds:      min,max             ! Optional, places bounds'+
                    ' on log10(rho) values. \n')
        infid.write('!Model Value Steps: stepsize            ! Optional, forces model'+
                    ' into discrete steps of stepsize. \n')
        infid.write('Debug Level:        '+str(debuglevel)+
                    ' '*19+'! Console output. '+
                    '0: minimal, 1: default, 2: detailed \n')
        infid.write('Iteration:          '+str(iteration)+
                    ' '*19+'! Iteration number,'+
                    ' use 0 for starting from scratch. \n')
        infid.write('Lagrange Value:     '+str(lagrange)+
                    ' '*17+'! log10(largrance '+
                    'multiplier), starting value.\n')
        infid.write('Roughness Value:    '+str(roughness)+
                    ' '*10+'! Roughness of last'+
                    ' model, ignored on startup. \n')
        infid.write('Misfit Value:       '+str(misfit)+
                    ' '*15+'! Misfit of model listed'+
                    'below. Ignored on startup.\n')
        infid.write('Misfit Reached:     0	                ! 0: not reached,'+
                    ' 1: reached.  Useful when restarting.\n')
        infid.write('Param Count:        '+str(paramcount)+
                    ' '*17+'! Number of free' +
                    ' inversion parameters. \n')
        for ii in range(paramcount):
            infid.write(ss+str(np.log10(rhostart))+'\n')
        
        infid.close()
        print 'Wrote Input File: ',self.inputfn
    
    def read1DModelFile(self):
        """
        
        will read in model 1D file
        
        Arguments:
        ----------
            **modelfn** : full path to model file
            
        Returns:
        --------
            **Occam1D.mdict** : dictionary of values with keys: 
                
                *'depth'* : depth of model in meters
                
                *'res'* : value of resisitivity
                
                *'pen'* : penalty
                
                *'pre'* : preference
                
                *'prefpen'* : preference penalty
                
        :Example: ::
            
            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> old.read1DModelFile()
        """
        if not self.modelfn:
            if not self.savepath:
                raise IOError('No model file found.  Please input one.')
            self.modelfn=os.path.join(self.savepath,'Model1D')
            
        mfid=open(self.modelfn,'r')
        mlines=mfid.readlines()
        mfid.close()
        try:
            self.mdict
        except AttributeError:
            mdict={}
            mdict['nparam']=0
            for key in ['depth','res','pen','pref','prefpen']:
                mdict[key]=[]
            
            for mm,mline in enumerate(mlines):
                if mline.find('!')==0:
                    pass
                elif mline.find(':')>=0:
                    mlst=mline.strip().split(':')
                    mdict[mlst[0]]=mlst[1]
                else:
                    mlst=mlst=mline.strip().split()
                    mdict['depth'].append(float(mlst[0]))
                    if mlst[1]=='?':
                        mdict['res'].append(-1)
                    elif mlst[1]=='1d12':
                        mdict['res'].append(1.0E12)
                    else:
                        try:
                            mdict['res'].append(float(mlst[1]))
                        except ValueError:
                            mdict['res'].append(-1)
                    mdict['pen'].append(float(mlst[2]))
                    mdict['pref'].append(float(mlst[3]))
                    mdict['prefpen'].append(float(mlst[4]))
                    if mlst[1]=='-1' or mlst[1]=='?':
                        mdict['nparam']+=1
                        
            #make everything an array
            for key in ['depth','res','pen','pref','prefpen']:
                    mdict[key]=np.array(mdict[key])
                    
            #create an array with empty columns to put the TE and TM models into
            mres=np.zeros((len(mdict['res']),3))
            mres[:,0]=mdict['res']
            mdict['res']=mres
            #make dictionary an attribute of Occam1D class            
            self.mdict=mdict
    
    def read1DInputFile(self):
        """
        reads in a 1D input file
        
        Arguments:
        ---------
            **inputfn** : full path to input file
               
        Returns:
        --------
            **Occam1D.indict** : dictionary with keys following the header and
            
                *'res'* : an array of resistivity values
                
        :Example: ::
            
            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> old.read1DInputFile()
        """
        if not self.inputfn:
            if not self.savepath:
                raise IOError('No input file found.  Please input one.')
            self.inputfn=os.path.join(self.savepath,'Input1D')
            
        infid=open(self.inputfn,'r')
        ilines=infid.readlines()
        infid.close()
    
        self.indict={}
        res=[]
        
        #split the keys and values from the header information
        for iline in ilines:
            if iline.find(':')>=0:
                ikey=iline[0:20].strip()
                ivalue=iline[20:].split('!')[0].strip()
                self.indict[ikey[:-1]]=ivalue
            else:
                try:
                    res.append(float(iline.strip()))
                except ValueError:
                    pass
                
        #make the resistivity array ready for models to be input
        self.indict['res']=np.zeros((len(res),3))
        self.indict['res'][:,0]=res
        
        #get data file
        if self.indict['Data File'].find('TE')>0:
            self.datafn_te=self.indict['Data File']
            
        elif self.indict['Data File'].find('TM')>0:
            self.datafn_tm=self.indict['Data File']
        


    def read1DdataFile(self,imode='TE'):
        """
        reads a 1D data file
        
        Arguments:
        ----------
            **datafile** : full path to data file
            
            **imode** : mode to read from can be TE or TM
        
        Returns:
        --------
            **Occam1D.rpdict** : dictionary with keys:
                
                *'freq'* : an array of frequencies with length nf
                
                *'resxy'* : TE resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr
                         
                *'resyx'* : TM resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr
                         
                *'phasexy'* : TE phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr
                
                *'phaseyx'* : TM phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr
                            
        :Example: ::
            
            >>> old = occam.Occam1d()
            >>> old.datafn_te = r"/home/Occam1D/Line1/Inv1_TE/MT01TE.dat"
            >>> old.read1DdataFile()
        """            
        
        #get the data file for the correct mode
        if imode=='TE':
            if not self.datafn_te:
                raise IOError('No TE data file found.  Please input one.')
                
            dfid=open(self.datafn_te,'r')
            
        elif imode=='TM':
            if not self.datafn_tm:
                raise IOError('No TM data file found.  Please input one.')
                
            dfid=open(self.datafn_te,'r')
        
        #read in lines
        dlines=dfid.readlines()
        dfid.close()
        
        #make a dictionary of all the fields found so can put them into arrays
        finddict={}
        for ii,dline in enumerate(dlines):
            if dline.find('#')<=3:
                fkey=dline[2:].strip().split(':')[0]
                fvalue=ii
                finddict[fkey]=fvalue
                
        #get number of frequencies
        nfreq=int(dlines[finddict['Frequencies']][2:].strip().split(':')[1].strip())
        
        #frequency list
        freq=np.array([float(ff) for ff in dlines[finddict['Frequencies']+1:
                                                            finddict['Receivers']]])
                
        #data dictionary to put things into
        #check to see if there is alread one, if not make a new one
        try:
            self.rpdict
        except NameError:
            self.rpdict={'freq':freq,
                         'resxy':np.zeros((4,nfreq)),
                         'resyx':np.zeros((4,nfreq)),
                         'phasexy':np.zeros((4,nfreq)),
                         'phaseyx':np.zeros((4,nfreq))
                         }
        
        #get data        
        for dline in dlines[finddict['Data']+1:]:
            if dline.find('!')==0:
                pass
            else:
                dlst=dline.strip().split()
                if len(dlst)>4:
                    jj=int(dlst[1])-1
                    dvalue=float(dlst[4])
                    derr=float(dlst[5])
                    if dlst[0]=='RhoZxy' or dlst[0]=='103':
                        self.rpdict['resxy'][0,jj]=dvalue
                        self.rpdict['resxy'][1,jj]=derr
                    if dlst[0]=='PhsZxy' or dlst[0]=='104':
                        self.rpdict['phasexy'][0,jj]=dvalue
                        self.rpdict['phasexy'][1,jj]=derr
                    if dlst[0]=='RhoZyx' or dlst[0]=='105':
                        self.rpdict['resyx'][0,jj]=dvalue
                        self.rpdict['resyx'][1,jj]=derr
                    if dlst[0]=='PhsZyx' or dlst[0]=='106':
                        self.rpdict['phaseyx'][0,jj]=dvalue
                        self.rpdict['phaseyx'][1,jj]=derr
        

    def read1DIterFile(self,iterfn,imode='TE'):
        """
        read an 1D iteration file
        
        Arguments:
        ----------
            **imode** : mode to read from 
        
        Returns:
        --------
            **Occam1D.itdict** : dictionary with keys of the header:
                
            **Occam1D.mdict['res']** : fills this array with the appropriate 
                                        values (0) for data, (1) TE, (2) TM
                                        
        :Example: ::
            
            >>> old = occam.Occam1d()
            >>> old.read1DIterFile(r"/home/Occam1D/Inv1_TE/M01TE_15.iter")
                
        """
        
        if not self.savepath:
            self.savepath=os.path.dirname(iterfn)
            
        self.read1DModelFile()
        
        freeparams=np.where(self.mdict['res']==-1)[0]
        
        if imode=='TE':
            self.iterfn_te=iterfn
            ifid=open(self.iterfn_te,'r')
        elif imode=='TM':
            self.iterfn_tm=iterfn
            ifid=open(self.iterfn_tm,'r')
            
        ilines=ifid.readlines()
        ifid.close()
        
        self.itdict={}
        model=[]    
        for ii,iline in enumerate(ilines):
            if iline.find(':')>=0:
                ikey=iline[0:20].strip()
                ivalue=iline[20:].split('!')[0].strip()
                self.itdict[ikey[:-1]]=ivalue
            else:
                try:
                    ilst=iline.strip().split()
                    for kk in ilst:
                        model.append(float(kk))
                except ValueError:
                    pass
        
        #put the model values into the model dictionary into the res array
        #for easy manipulation and access.  Also so you can compare TE and TM        
        model=np.array(model)
        if imode=='TE':
            self.mdict['res'][:,1]=self.mdict['res'][:,0]
            self.mdict['res'][freeparams,1]=model
        if imode=='TM':
            self.mdict['res'][:,2]=self.mdict['res'][:,0]
            self.mdict['res'][freeparams,2]=model


    def read1DRespFile(self,respfn,imode='TE'):
        """
        read response file
        
        Arguments:
        ---------
            **repsfn** : full path to response file
        
        Returns:
        --------
            *Occam1D.*rpdict** : dictionary with keys:
                
                *freq* : an array of frequencies with length nf
                
                *resxy* : TE resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr
                         
                *resyx* : TM resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr
                         
                *phasexy* : TE phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr
                
                *phaseyx* : TM phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr
       
       :Example: ::
            
            >>> old = occam.Occam1d()
            >>> old.read1DRespFile(r"/home/Occam1D/Inv1_TE/M01TE_15.resp")
        """
        
        if imode=='TE':
            self.respfn_te=respfn
        elif imode=='TM':
            self.respfn_tm=respfn
        
        if not self.savepath:
            self.savepath=os.path.dirname(respfn)
            
        dfid=open(respfn,'r')
        
        dlines=dfid.readlines()
        dfid.close()
    
        finddict={}
        for ii,dline in enumerate(dlines):
            if dline.find('#')<=3:
                fkey=dline[2:].strip().split(':')[0]
                fvalue=ii
                finddict[fkey]=fvalue
        nfreq=int(dlines[finddict['Frequencies']][2:].strip().split(':')[1].strip())
        
        #frequency list
        freq=np.array([float(ff) for ff in dlines[finddict['Frequencies']+1:
                                                            finddict['Receivers']]])
                
        #data dictionary
        try:
            self.rpdict
        except AttributeError:
            self.rpdict={'freq':freq,
                        'resxy':np.zeros((4,nfreq)),
                        'resyx':np.zeros((4,nfreq)),
                        'phasexy':np.zeros((4,nfreq)),
                        'phaseyx':np.zeros((4,nfreq))
                        }
                
        for dline in dlines[finddict['Data']+1:]:
            if dline.find('!')==0:
                pass
            else:
                dlst=dline.strip().split()
                if len(dlst)>4:
                    jj=int(dlst[1])-1
                    dvalue=float(dlst[4])
                    derr=float(dlst[5])
                    rvalue=float(dlst[6])
                    rerr=float(dlst[7])
                    if dlst[0]=='RhoZxy' or dlst[0]=='103':
                        self.rpdict['resxy'][0,jj]=dvalue
                        self.rpdict['resxy'][1,jj]=derr
                        self.rpdict['resxy'][2,jj]=rvalue
                        self.rpdict['resxy'][3,jj]=rerr
                    if dlst[0]=='PhsZxy' or dlst[0]=='104':
                        self.rpdict['phasexy'][0,jj]=dvalue
                        self.rpdict['phasexy'][1,jj]=derr
                        self.rpdict['phasexy'][2,jj]=rvalue
                        self.rpdict['phasexy'][3,jj]=rerr
                    if dlst[0]=='RhoZyx' or dlst[0]=='105':
                        self.rpdict['resyx'][0,jj]=dvalue
                        self.rpdict['resyx'][1,jj]=derr
                        self.rpdict['resyx'][2,jj]=rvalue
                        self.rpdict['resyx'][3,jj]=rerr
                    if dlst[0]=='PhsZyx' or dlst[0]=='106':
                        self.rpdict['phaseyx'][0,jj]=dvalue
                        self.rpdict['phaseyx'][1,jj]=derr
                        self.rpdict['phaseyx'][2,jj]=rvalue
                        self.rpdict['phaseyx'][3,jj]=rerr
    
    def plot1D(self,iternum=10,savepath=None,iterfn=None,respfn=None,
               imode='TE',fignum=1,ms=4,dpi=150,fs=10,lw=2,dlimits=None):
        """
        Plots the results of a 1D inversion.  The left plot is the response
        and the right hand plot is the model as a function of depth.
        
        Arguments:
        ----------
            **iternum** : iteration number to plot. *Default* is 10 
            
            **savepath** : path to the 1D inversion response and iteration 
                           files. *Default* is None.
                           
            **iterfn** : full path to iteration file if savepath is not entered
                         *Default* is None.
                        
            **respfn** : full path to response file if savepath is not entered.
                         *Default* is None.
            
            **imode** : mode to plot. can be input as:
                * 'TE' for TE mode
                * 'TM' for TM mode
                * 'both' for both TE and TM modes
                * Default* is 'TE'
            
            **fignum** : figure number that the plot will be. *Default* is 1
            
            **ms** : marker size. *Default* is 4
            
            **dpi** : dots per inch resolution of the plot. *Default* is 150.
            
            **fs** : font size of the labels. *Default* is 10
            
            **dlimits** : limits on the depth axes. Input as a tuple 
                          (dmin,dmax).  *Default* is None.
            
        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> #Look at only the interval between 5 and 10 kilometers
            >>> old.plot1D(dlimits=(5,10))
        """ 
        
        self.iternum=iternum
        #get files
        try:
            self.modelfn
        except AttributeError:
            if not self.dirpath:
                self.dirpath=os.path.dirname(respfn)
                
            self.modelfn=os.path.join(self.dirpath,'Model1D')
            if os.path.isfile(self.modelfn)==False:
                raise IOError('Could not find '+self.modelfn)
        
        #-------------read in response files---------------------
        if respfn==None:
            if imode=='TE':
                try:
                    self.respfn_te=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TE_{0}.resp'.format(self.iternum))>0][0]
                    
                    self.read1DRespFile(self.respfn_te,imode='TE')
                    
                except IndexError:
                    raise IOError('Could not find response TE file.')
            elif imode=='TM':
                try:
                    self.respfn_tm=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TM_{0}.resp'.format(self.iternum))>0][0]
                    
                    self.read1DRespFile(self.respfn_tm,imode='TM')
                except IndexError:
                    raise IOError('Could not find response TM file.')
            elif imode=='both':
                try:
                    self.respfn_te=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TE_{0}.resp'.format(self.iternum))>0][0]
                    self.read1DRespFile(self.respfn_te,imode='TE')
                except IndexError:
                    raise IOError('Could not find response TE file.')
                try:
                    self.respfn_tm=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TM_{0}.resp'.format(self.iternum))>0][0]
                    self.read1DRespFile(self.respfn_tm,imode='TM')
                except IndexError:
                    raise IOError('Could not find response TM file.')
        
        #if the response files are input read them in 
        else:
            if imode=='TE':
                self.respfn_te=respfn
                self.read1DRespFile(self.respfn_te,imode='TE')
                
            elif imode=='TM':
                self.respfn_tm=respfn
                self.read1DRespFile(self.respfn_tm,imode='TM')
                
            elif imode=='both':
                if type(iterfn) is not list or type(iterfn) is not tuple:
                    raise IOError('Please enter iteration files as a list or tuple.')
                self.respfn_te=respfn[0]
                self.read1DRespFile(self.respfn_te,imode='TE')
                
                self.respfn_tm=respfn[1]
                self.read1DRespFile(self.respfn_tm,imode='TM')
        
        #------Read in iteration files--------------------        
        if iterfn==None:
            if imode=='TE':
                try:
                    self.iterfn_te=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TE_{0}.iter'.format(self.iternum))>0][0]
                    print self.iterfn_te
                    self.read1DIterFile(self.iterfn_te,imode='TE')
                    
                except IndexError:
                    raise IOError('Could not find iteration TE file.')
            elif imode=='TM':
                try:
                    self.iterfn_tm=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TM_{0}.iter'.format(self.iternum))>0][0]
                    
                    self.read1DIterFile(self.iterfn_tm,imode='TM')
                except IndexError:
                    raise IOError('Could not find iteration TM file.')
            elif imode=='both':
                try:
                    self.iterfn_te=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TE_{0}.iter'.format(self.iternum))>0][0]
                    self.read1DIterFile(self.iterfn_te,imode='TE')
                except IndexError:
                    raise IOError('Could not find iteration TE file.')
                try:
                    self.iterfn_tm=[os.path.join(self.savepath,rr) 
                                    for rr in os.listdir(self.savepath)
                                    if rr.find('TM_{0}.iter'.format(self.iternum))>0][0]
                    self.read1DIterFile(self.iterfn_tm,imode='TM')
                except IndexError:
                    raise IOError('Could not find iteration TM file.')
        else:
            if imode=='TE':
                self.iterfn_te=iterfn
                self.read1DIterFile(self.iterfn_te,imode='TE')
                
            elif imode=='TM':
                self.iterfn_tm=iterfn
                self.read1DIterFile(self.iterfn_tm,imode='TM')
                
            elif imode=='both':
                if type(iterfn) is not list or type(iterfn) is not tuple:
                    raise IOError('Please enter iteration files as a list or tuple.')
                self.iterfn_te=iterfn[0]
                self.read1DIterFile(self.iterfn_te,imode='TE')
                
                self.iterfn_tm=iterfn[1]
                self.read1DIterFile(self.iterfn_tm,imode='TM')
                
        period=1/self.rpdict['freq']
        
        #make a grid of subplots
        gs=gridspec.GridSpec(6,5,hspace=.25,wspace=.75)
        
        #make a figure
        fig=plt.figure(fignum,[8,8],dpi=dpi)
        plt.clf()
        
        #set some plot parameters
        plt.rcParams['font.size']=fs-2
        plt.rcParams['figure.subplot.left']=.1
        plt.rcParams['figure.subplot.right']=.93
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.top']=.90
        
        #subplot resistivity
        axr=fig.add_subplot(gs[:4,:4])
        
        #subplot for phase
        axp=fig.add_subplot(gs[4:,:4],sharex=axr)
        
        #check for data in resistivity
        rxy=np.where(self.rpdict['resxy'][0]!=0)[0]
        ryx=np.where(self.rpdict['resyx'][0]!=0)[0]
        
        pxy=np.where(self.rpdict['phasexy'][0]!=0)[0]
        pyx=np.where(self.rpdict['phaseyx'][0]!=0)[0]
        
        #check to make sure a model was read in for resistivity
        rxym=np.where(self.rpdict['resxy'][2]!=0)[0]
        ryxm=np.where(self.rpdict['resyx'][2]!=0)[0]
        
        pxym=np.where(self.rpdict['phasexy'][2]!=0)[0]
        pyxm=np.where(self.rpdict['phaseyx'][2]!=0)[0]
        
        #----------Plot TE mode-------------------
        if imode=='TE':
            titlestr='$Z_{TE}$'
            #plot data resistivity 
            if len(rxy)!=0:
                r1=axr.loglog(period[rxy],self.rpdict['resxy'][0][rxy],
                              ls='None',marker='o',color='k',mfc='k',ms=ms)

            #plot data phase
            if len(pxy)!=0:
                p1=axp.semilogx(period[pxy],self.rpdict['phasexy'][0][pxy],
                          ls='None',marker='o',color='k',mfc='k',ms=ms)
                          
            #plot model resistivity
            if len(rxym)!=0:
                r2=axr.loglog(period[rxym],self.rpdict['resxy'][2][rxym],
                              ls=':',color='b',lw=lw)
            #plot model phase                 
            if len(pxym)!=0:
                p2=axp.semilogx(period[pxym],self.rpdict['phasexy'][2][pxym],
                          ls=':',color='b',lw=lw)
            
            #add legend
            axr.legend([r1[0],r2[0]],['Data','Model'],loc='upper left',
                       markerscale=1,
                       borderaxespad=.15,
                       labelspacing=.18,
                       handletextpad=.15,borderpad=.15)
        
        #--------Plot TM mode-----------------------
        elif imode=='TM':
            titlestr='$Z_{TM}$'
            #plot data resistivity 
            if len(ryx)!=0:
                r1=axr.loglog(period[ryx],self.rpdict['resyx'][0][ryx],
                              ls='None',marker='o',color='k',mfc='k',ms=ms)
            #plot data phase
            if len(pyx)!=0:
                p1=axp.semilogx(period[pyx],self.rpdict['phaseyx'][0][pyx],
                          ls='None',marker='o',color='k',mfc='k',ms=ms)
                          
            #plot model resistivity
            if len(ryxm)!=0:
                r2=axr.loglog(period[ryxm],self.rpdict['resyx'][2][ryxm],
                              ls=':',color='b',lw=lw)
            #plot model phase                 
            if len(pyxm)!=0:
                p2=axp.semilogx(period[pyxm],self.rpdict['phaseyx'][2][pyxm],
                          ls=':',color='b',lw=lw)
                
            axr.legend([r1[0],r2[0]],['Data','Model'],
                       loc='upper left',markerscale=1,
                       borderaxespad=.15,
                       labelspacing=.18,
                       handletextpad=.15,borderpad=.15)
        
        #-------------Plot Both Modes--------------------------------
        elif imode=='both':
            titlestr='$Z_{TE}$ and $Z_{TM}$'
            #plot data resistivity 
            if len(rxy)!=0:
                r1te=axr.loglog(period[rxy],self.rpdict['resxy'][0][rxy],
                              ls='None',marker='s',color='k',mfc='k',ms=ms)
            if len(ryx)!=0:
                r1tm=axr.loglog(period[ryx],self.rpdict['resyx'][0][ryx],
                              ls='None',marker='o',color='k',mfc='k',ms=ms)

            #plot data phase
            if len(pxy)!=0:
                p1te=axp.semilogx(period[pxy],self.rpdict['phasexy'][0][pxy],
                          ls='None',marker='s',color='k',mfc='k',ms=ms)
            
            if len(pyx)!=0:
                p1tm=axp.semilogx(period[pyx],self.rpdict['phaseyx'][0][pyx],
                          ls='None',marker='o',color='k',mfc='k',ms=ms)
                          
            #plot model resistivity
            if len(rxym)!=0:
                r2te=axr.loglog(period[rxym],self.rpdict['resxy'][2][rxym],
                              ls=':',color='b',lw=lw)
        
            if len(ryxm)!=0:
                r2tm=axr.loglog(period[ryxm],self.rpdict['resyx'][2][ryxm],
                              ls=':',color='r',lw=lw)
            #plot model phase                 
            if len(pxym)!=0:
                p2=axp.semilogx(period[pxym],self.rpdict['phasexy'][2][pxym],
                          ls=':',color='b',lw=lw)
            if len(pyxm)!=0:
                p2=axp.semilogx(period[pyxm],self.rpdict['phaseyx'][2][pyxm],
                          ls=':',color='r',lw=lw)
            
            #add legend
            axr.legend([r1te[0],r2te[0],r1tm[0],r2tm[0]],
                       ['Data$_{TE}$','Model$_{TE}$',
                        'Data$_{TM}$','Model$_{TM}$'],
                        loc='upper left',markerscale=1,
                       borderaxespad=.15,
                       labelspacing=.18,
                       handletextpad=.15,borderpad=.15)

                          
        axr.grid(True,alpha=.4,which='both')
        plt.setp(axr.xaxis.get_ticklabels(),visible=False)
        axp.grid(True,alpha=.4,which='both')
        axp.yaxis.set_major_locator(MultipleLocator(10))
        axp.yaxis.set_minor_locator(MultipleLocator(1))
        
        axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                       fontdict={'size':fs,'weight':'bold'})
        axp.set_ylabel('Phase (deg)',
                       fontdict={'size':fs,'weight':'bold'})
        axp.set_xlabel('Period (s)',fontdict={'size':fs,'weight':'bold'})
        plt.suptitle(titlestr,fontsize=fs+2,fontweight='bold')
        
        #------------------plot 1D inversion---------------------------------
        axm=fig.add_subplot(gs[:,4])
        depthp=self.mdict['depth'][1:]
        if imode=='TE':
            modelresp=abs(10**self.mdict['res'][1:,1])
            axm.loglog(modelresp[::-1],depthp[::-1],ls='steps-',color='b',
                       lw=lw)
        elif imode=='TM':
            modelresp=abs(10**self.mdict['res'][1:,2])
            axm.loglog(modelresp[::-1],depthp[::-1],ls='steps-',color='b',
                       lw=lw)
        elif imode=='both':
            modelrespte=abs(10**self.mdict['res'][1:,1])
            axm.loglog(modelrespte[::-1],depthp[::-1],ls='steps-',color='b',
                       lw=lw)
            modelresptm=abs(10**self.mdict['res'][1:,2])
            axm.loglog(modelresptm[::-1],depthp[::-1],ls='steps-',color='r',
                       lw=lw)
        
        if dlimits==None:
            axm.set_ylim(ymin=depthp[-1],ymax=depthp[0])
        else:
            axm.set_ylim(dlimits)
        axm.set_ylabel('Depth (m)',fontdict={'size':fs,'weight':'bold'})
        axm.set_xlabel('Resistivity ($\Omega \cdot m$)',
                       fontdict={'size':fs,'weight':'bold'})
        axm.grid(True,which='both',alpha=.4)
        
        plt.show()
        
    def plotL2Curve(self,savepath=None,imode='TE',fignum=1,dpi=150,fs=10):
        """
        Plot the L curve for RMS vs Iteration and RMS vs Roughness.
        
        Arguments:
        ----------
        
            **savepath** : path to iteration files
            
            **imode** : mode to plot.  Can be:
                    * 'TE' for TE mode. *Default*
                    * 'TM' for TM mode
                    * 'both' for both TE ant TM modes
                    
            **fignum** : figure number for the plot
            
            **dpi** : dots per inch resolution of the plot
            
            **fs** : font size of labels
            
        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> old.plotL2Curve()
        """
        
        if savepath==None:
            if not self.savepath:
                raise IOError('No savepath found, please enter one.')
                
        else:
            self.savepath=savepath
        
        self.rms_te=[]
        self.rms_tm=[]
        self.roughness_te=[]
        self.roughness_tm=[]
        
        #get rms and roughness from each iteration for the different modes
        if imode=='TE':
            #get all iteration files for TE mode
            iterlstte=[os.path.join(self.savepath,itfn) 
                     for itfn in os.listdir(self.savepath)
                     if itfn.find('TE')>0 and itfn.find('iter')>0]
                     
            self.rms_te=np.zeros(len(iterlstte))
            self.roughness_te=np.zeros(len(iterlstte))
            
            #get rms and roughness
            for itfn in iterlstte:
                self.read1DIterFile(itfn,imode='TE')
                
                #get iteration number to make sure the items are in sequence
                ii=int(self.itdict['Iteration'])
                
                #put the values in appropriate place
                self.rms_te[ii]=float(self.itdict['Misfit Value'])
                self.roughness_te[ii]=float(self.itdict['Roughness Value'])
            
        
        elif imode=='TM':
            #get all iteration files for TM mode
            iterlsttm=[os.path.join(self.savepath,itfn) 
                     for itfn in os.listdir(self.savepath)
                     if itfn.find('TM')>0 and itfn.find('iter')>0]
            
            self.rms_tm=np.zeros(len(iterlsttm))
            self.roughness_tm=np.zeros(len(iterlsttm))
            
            #get rms and roughness
            for itfn in iterlsttm:
                self.read1DIterFile(itfn,imode='TM')
                
                #get iteration number to make sure the items are in sequence
                ii=int(self.itdict['Iteration'])
                
                #put the values in appropriate place
                self.rms_tm[ii]=float(self.itdict['Misfit Value'])
                self.roughness_tm[ii]=float(self.itdict['Roughness Value'])
                
        elif imode=='both':
            #get all iteration files for TE mode
            iterlstte=[os.path.join(self.savepath,itfn) 
                     for itfn in os.listdir(self.savepath)
                     if itfn.find('TE')>0 and itfn.find('iter')>0]
                     
            self.rms_te=np.zeros(len(iterlstte))
            self.roughness_te=np.zeros(len(iterlstte))
            
            #get rms and roughness
            for itfn in iterlstte:
                self.read1DIterFile(itfn,imode='TE')
                
                #get iteration number to make sure the items are in sequence
                ii=int(self.itdict['Iteration'])
                
                #put the values in appropriate place
                self.rms_te[ii]=float(self.itdict['Misfit Value'])
                self.roughness_te[ii]=float(self.itdict['Roughness Value'])
                
            #get all iteration files for TM mode
            iterlsttm=[os.path.join(self.savepath,itfn) 
                     for itfn in os.listdir(self.savepath)
                     if itfn.find('TM')>0 and itfn.find('iter')>0]
            
            self.rms_tm=np.zeros(len(iterlsttm))
            self.roughness_tm=np.zeros(len(iterlsttm))
            
            #get rms and roughness
            for itfn in iterlsttm:
                self.read1DIterFile(itfn,imode='TM')
                
                #get iteration number to make sure the items are in sequence
                ii=int(self.itdict['Iteration'])
                
                #put the values in appropriate place
                self.rms_tm[ii]=float(self.itdict['Misfit Value'])
                self.roughness_tm[ii]=float(self.itdict['Roughness Value'])
        
        #plot the rms vs iteration, roughness vs rms
        #---------plot TE mode-------------------
        if imode=='TE':
            fig=plt.figure(fignum,dpi=dpi)
            plt.clf()
            ax1=fig.add_subplot(1,1,1)            
            
            nr=len(self.rms_te)
            #plot the rms vs iteration
            l1,=ax1.plot(np.arange(1,nr,1),self.rms_te[1:],'-k',lw=1,
                         marker='d',ms=5)
            
            #plot the median of the RMS
            medte=np.median(self.rms_te[1:])
            m1,=ax1.plot(np.arange(0,nr,1),
                         np.repeat(medte,nr),
                         '--r',lw=.75)
        
            #make subplot for RMS vs Roughness Plot
            ax2=ax1.twiny()
            
            #plot the rms vs roughness 
            l2,=ax2.plot(self.roughness_te[1:],self.rms_te[1:],
                         '--b',lw=.75,marker='o',ms=7,mfc='white')
            for ii,rms in enumerate(self.rms_te[1:],1):
                ax2.text(self.roughness_te[ii],rms,'{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size':6,'weight':'bold','color':'blue'})
            
            #make a legend
            ax1.legend([l1,l2,m1],['RMS_TE','Roughness_TE',
                       'Median_RMS={0:.2f}'.format(medte)],
                        ncol=4,loc='upper center',columnspacing=.25,
                        markerscale=.75,handletextpad=.15)
            
            ax1.set_ylim(medte-1,medte+1)
            ax1.set_ylabel('RMS',fontdict={'size':fs,'weight':'bold'})                                   
            ax1.set_xlabel('Iteration',fontdict={'size':fs,'weight':'bold'})
            ax1.grid(alpha=.25,which='both')
            ax2.set_xlabel('Roughness',fontdict={'size':fs,'weight':'bold',
                                                 'color':'blue'})
            for t2 in ax2.get_xticklabels():
                t2.set_color('blue')  
                
        #-------Plot TM mode-------------------
        elif imode=='TM':
            fig=plt.figure(fignum,dpi=dpi)
            plt.clf()
            ax1=fig.add_subplot(1,1,1)
            
            nr=len(self.rms_tm)
            #plot the rms vs iteration
            l1,=ax1.plot(np.arange(1,nr,1),self.rms_tm[1:],'-k',lw=1,
                         marker='d',ms=5)
            
            #plot the median of the RMS
            medtm=np.median(self.rms_tm[1:])
            m1,=ax1.plot(np.arange(0,nr,1),
                         np.repeat(medtm,nr),
                         '--r',lw=.75)

        
            #make subplot for RMS vs Roughness Plot
            ax2=ax1.twiny()
            
            #plot the rms vs roughness 
            l2,=ax2.plot(self.roughness_tm[1:],self.rms_tm[1:],
                         '--b',lw=.75,marker='o',ms=7,mfc='white')
            for ii,rms in enumerate(self.rms_tm[1:],1):
                ax2.text(self.roughness_tm[ii],rms,'{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size':fs-2,'weight':'bold','color':'blue'})
            
            #make a legend
            ax1.legend([l1,l2,m1],['RMS_TM','Roughness_TM',
                       'Median_RMS={0:.2f}'.format(medtm)],
                        ncol=4,loc='upper center',columnspacing=.25,
                        markerscale=.75, handletextpad=.15)
            
            ax1.set_ylim(medtm-1,medtm+1)            
            ax1.set_ylabel('RMS',fontdict={'size':fs,'weight':'bold'})                                   
            ax1.set_xlabel('Iteration',fontdict={'size':fs,'weight':'bold'})
            ax1.grid(alpha=.25,which='both')
            ax2.set_xlabel('Roughness',fontdict={'size':fs,'weight':'bold',
                                                 'color':'blue'})
            for t2 in ax2.get_xticklabels():
                t2.set_color('blue')  
                    
        elif imode=='both':
            fig=plt.figure(fignum,dpi=dpi)
            plt.clf()
            ax1=fig.add_subplot(2,1,1)
            ax3=fig.add_subplot(2,1,2,sharex=ax1) 
            
            plt.rcParams['figure.subplot.hspace']=.4
            plt.rcParams['figure.subplot.left']=.1
            plt.rcParams['figure.subplot.right']=.97
            plt.rcParams['figure.subplot.bottom']=.1
            plt.rcParams['figure.subplot.top']=.92
            
            nr=len(self.rms_te)
            #plot the rms vs iteration
            l1,=ax1.plot(np.arange(1,nr,1),self.rms_te[1:],'-k',lw=1,
                         marker='d',ms=5)
            
            #plot the median of the RMS
            medte=np.median(self.rms_te[1:])
            m1,=ax1.plot(np.arange(0,nr,1),
                         np.repeat(medte,nr),
                         '--r',lw=.75)
        
            #make subplot for RMS vs Roughness Plot
            ax2=ax1.twiny()
            
            #plot the rms vs roughness 
            l2,=ax2.plot(self.roughness_te[1:],self.rms_te[1:],
                         '--b',lw=.75,marker='o',ms=7,mfc='white')
            for ii,rms in enumerate(self.rms_te[1:],1):
                ax2.text(self.roughness_te[ii],rms,'{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size':fs-2,'weight':'bold','color':'blue'})
            
            #make a legend
            ax1.legend([l1,l2,m1],['RMS_TE','Roughness_TE',
                       'Median_RMS={0:.2f}'.format(medte)],
                        ncol=4,loc='upper center',columnspacing=.25,
                        markerscale=.75,handletextpad=.15)
            
            ax1.set_ylim(medte-1,medte+1)
            ax1.set_ylabel('RMS',fontdict={'size':fs,'weight':'bold'})                                   
            #ax1.set_xlabel('Iteration',fontdict={'size':8,'weight':'bold'})
            ax1.grid(alpha=.25,which='both')
            ax2.set_xlabel('Roughness',fontdict={'size':fs,'weight':'bold',
                                                 'color':'blue'})
            for t2 in ax2.get_xticklabels():
                t2.set_color('blue') 
            
            #plot TM
            nr=len(self.rms_te)
            #plot the rms vs iteration
            l3,=ax3.plot(np.arange(1,nr,1),self.rms_tm[1:],'-k',lw=1,
                         marker='d',ms=5)
            
            #plot the median of the RMS
            medtm=np.median(self.rms_tm[1:])
            m3,=ax3.plot(np.arange(0,nr,1),
                         np.repeat(medtm,nr),
                         '--r',lw=.75)

        
            #make subplot for RMS vs Roughness Plot
            ax4=ax3.twiny()
            
            #plot the rms vs roughness 
            l4,=ax4.plot(self.roughness_tm[1:],self.rms_tm[1:],
                         '--b',lw=.75,marker='o',ms=7,mfc='white')
            for ii,rms in enumerate(self.rms_tm[1:],1):
                ax4.text(self.roughness_tm[ii],rms,'{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size':6,'weight':'bold','color':'blue'})
            
            #make a legend
            ax3.legend([l1,l2,m1],['RMS_TM','Roughness_TM',
                       'Median_RMS={0:.2f}'.format(medtm)],
                        ncol=4,loc='upper center',columnspacing=.25,
                        markerscale=.75, handletextpad=.15)
            
            ax3.set_ylim(medtm-1,medtm+1)            
            ax3.set_ylabel('RMS',fontdict={'size':fs,'weight':'bold'})                                   
            ax3.set_xlabel('Iteration',fontdict={'size':fs,'weight':'bold'})
            ax3.grid(alpha=.25,which='both')
            ax4.set_xlabel('Roughness',fontdict={'size':fs,'weight':'bold',
                                                 'color':'blue'})
            for t2 in ax4.get_xticklabels():
                t2.set_color('blue') 
        
        plt.show()
