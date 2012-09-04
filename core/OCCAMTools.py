# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 11:54:38 2011

@author: a1185872
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import Z
import MTtools as mt
import fnmatch
import LatLongUTMconversion as utm2ll
from operator import itemgetter
from matplotlib.ticker import MultipleLocator
import WinglinkTools as wlt
import shutil
import subprocess
import time

occamdict={'1':'resxy','2':'phasexy','3':'realtip','4':'imagtip','5':'resyx',
           '6':'phaseyx'}

def make1DdataFile(station,edipath=None,savepath=None,polarization='both',
                   reserr='data',phaseerr='data',fmt='%.6e',ss=3*' ',thetar=0):
    """
    make1Ddatafile will write a data file for Occam1D

    Input:
        station = the station name and path if edipath=None
        edipath = path to the edi files to be written into a data file,
                  useful for multile data files
        savepath = path to save the file
        thetar = rotation angle to rotate Z
    
    
    #===============================================================================
    # Input parameters
    #===============================================================================
    
    edipath=r"C:\Peacock\My Dropbox\Paralana\EDIFilesBaseSurvey\SS" #path to edi's
    station='pb01'  #station to create data file for                   
    savepath=r"C:\Peacock\PHD\OCCAM\1DInversion\Base"  #path to save the data file
    polarization='both'  #choose component either xy or yx or both
    reserr='data'   #either data or set as a percentage
    phaseerr='data' #either data or set as percentage -> 10 for 10 percent error
    ss='   '
    fmt='%.6e'
    #===============================================================================
    # open up edi file and write data file
    #===============================================================================
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
            
    #raise an error if can't find the edifile        
    if edifile==None:
        raise NameError('No edifile exists, check path and station name')

    #read in edifile    
    impz=Z.Z(edifile)    
    
    #make sure the savepath exists, if not create it
    if savepath==None:
        savepath=os.path.dirname(edifile)
        if not os.path.exists(savepath):
            os.mkdir(savepath)
    if savepath.find('.')>0:
        if not os.path.exists(os.path.dirname(savepath)):
            os.mkdir(os.path.dirname(savepath))
    
    #load the edifile and get resistivity and phase
    
    rp=impz.getResPhase(thetar=thetar)
    freq=impz.frequency
    nf=len(freq)
    returnfn=[]
    
    if polarization=='both':
        for pol in ['xy','yx']:
            if savepath.find('.')==0:
                dfilesave=os.path.join(savepath,impz.station+pol.upper()+'.dat')
            else:
                dfilesave=savepath
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
            datafid.write('!'+2*ss+'Type'+2*ss+'Freq#'+2*ss+'Tx#'+2*ss+'Rx#'+
                          2*ss+'Data'+2*ss+'Std_Error'+'\n')
                          
            #put the yx phase component in the first quadrant as prescribed
            if pol=='yx':
                rp.phaseyx=rp.phaseyx+180
                #check if there are any negative phases
#                negphase=np.where(rp.phaseyx>180)
#                if len(negphase)>0:
#                    rp.phaseyx[negphase[0]]=rp.phaseyx\
#                                                    [negphase[0]]-360
            if pol=='xy':
                negphase=np.where(rp.phaseyx<0)
                if len(negphase)>0:
                    rp.phasexy[negphase[0]]=rp.phasexx[negphase[0]]+180
            #write the resistivity and phase components
            for ii in range(nf):
                #write resistivity components
                if reserr=='data':
                    if pol=='xy':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.resxy[ii]+2*ss+
                                      fmt % rp.resxyerr[ii]+'\n')
                    elif pol=='yx':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.resyx[ii]+2*ss+
                                      fmt % rp.resyxerr[ii]+'\n')
                else:
                    if pol=='xy':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.resxy[ii]+2*ss+
                                      fmt % (rp.resxy[ii]*reserr/100.)+'\n')
                    elif pol=='yx':
                        datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.resyx[ii]+2*ss+
                                      fmt % (rp.resyx[ii]*reserr/100.)+'\n')
                
                #write phase components
                if phaseerr=='data':
                    if pol=='xy':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.phasexy[ii]+2*ss+
                                      fmt % rp.phasexyerr[ii]+'\n')
                    if pol=='yx':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.phaseyx[ii]+2*ss+
                                      fmt % rp.phaseyxerr[ii]+'\n')
                else:
                    if pol=='xy':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.phasexy[ii]+2*ss+
                                      fmt % (phaseerr/100.*(180/np.pi))+'\n')
                    if pol=='yx':
                        datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                      2*ss+'1'+2*ss+fmt % rp.phaseyx[ii]+2*ss+
                                      fmt % (phaseerr/100.*(180/np.pi))+'\n')
            datafid.write('\n')
            datafid.close()
            print 'Wrote Data File: ',dfilesave
            returnfn.append(dfilesave)
        return returnfn[0],returnfn[1]
    else:
        if polarization=='TE':
            pol='xy'
        elif polarization=='TM':
            pol='yx'
        
        if savepath.find('.')==-1:
            dfilesave=os.path.join(savepath,impz.station+pol.upper()+'.dat')
        else:
            dfilesave=savepath
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
        datafid.write('!'+2*ss+'Type'+2*ss+'Freq#'+2*ss+'Tx#'+2*ss+'Rx#'+2*ss+
                      'Data'+2*ss+'Std_Error'+'\n')
                      
        #put the yx phase component in the first quadrant as prescribed
        if pol=='TM':
                rp.phaseyx=rp.phaseyx+180
                
        #write the resistivity and phase components
        for ii in range(nf):
            #write resistivity components
            if reserr=='data':
                if pol=='xy':
                    datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.resxy[ii]+2*ss+
                                  fmt % rp.resxyerr[ii]+'\n')
                elif pol=='yx':
                    datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.resyx[ii]+2*ss+
                                  fmt % rp.resyxerr[ii]+'\n')
            else:
                if pol=='xy':
                    datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.resxy[ii]+2*ss+
                                  fmt % (rp.resxy[ii]*reserr/100.)+'\n')
                elif pol=='yx':
                    datafid.write(2*ss+'RhoZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.resyx[ii]+2*ss+
                                  fmt % (rp.resyx[ii]*reserr/100.)+'\n')
            
            #write phase components
            if phaseerr=='data':
                if pol=='xy':
                    datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.phasexy[ii]+2*ss+
                                  fmt % rp.phasexyerr[ii]+'\n')
                if pol=='yx':
                    datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.phaseyx[ii]+2*ss+
                                  fmt % rp.phaseyxerr[ii]+'\n')
            else:
                if pol=='xy':
                    datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.phasexy[ii]+2*ss+
                                  fmt % (phaseerr/100.*(180/np.pi))+'\n')
                if pol=='yx':
                    datafid.write(2*ss+'PhsZ'+pol+2*ss+str(ii+1)+2*ss+'0'+
                                  2*ss+'1'+2*ss+fmt % rp.phaseyx[ii]+2*ss+
                                  fmt % (phaseerr/100.*(180/np.pi))+'\n')
        datafid.close()
        print 'Wrote Data File: ',dfilesave
        return dfilesave

def make1DModelFile(savepath,nlayers=100,bottomlayer=10000,basestep=10,
                    z1layer=50,airlayerheight=10000):
    """
    Makes a 1D model file
    
    Input:
        savepath =path to save file to, if just path saved as savepath\model.mod
        nlayers = number of layers
        bottomlayer = depth of bottom layer in meters
        basestep = numerical base of logarithmic depth step 10 or 2 or 1 for 
                    linear
        z1layer = depth of first layer in meters
        airlayerheight = height of air layers in meters
    Output:
        modelfilename = full path to model file
    """
    
    #===========================================================================
    # Make model file
    #===========================================================================
    #make sure the savepath exists, if not create it
    
    ss='   '
    if savepath.find('.')==-1:
        if not os.path.exists(savepath):
            os.mkdir(savepath)
        modfn=os.path.join(savepath,'model.txt')
    else:
        modfn=savepath
    
    #---------need to refine this-------------------- 
    #make an array of layer depths on a scale specified by base
    layers=np.zeros(nlayers+1)
    layers[0]=z1layer
    
    if basestep==10:
    #    l=np.log10(np.logspace(0,np.log10(bottomlayer),num=nlayers))
        l=np.linspace(1,np.log10(bottomlayer)+1,num=nlayers+1)
    
    elif basestep==2:
    #    l=np.log2(np.logspace(0,np.log2(bottomlayer),num=nlayers)
        l=np.linspace(1,np.log2(bottomlayer)+1,num=nlayers+1)
    elif basestep==1:
        l=np.linspace(1,np.log10(bottomlayer)+1,num=nlayers+1)
    else:
        raise ValueError('Base '+str(basestep)+' not supported')
    
    for ll in range(1,nlayers+1):
        layers[ll]=layers[ll-1]+z1layer*np.floor(l[ll])
    
    
    #make the model file
    modfid=open(modfn,'w')
    modfid.write('Format: Resistivity1DMod_1.0'+'\n')
    modfid.write('#LAYERS:    '+str(nlayers+3)+'\n')
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
        modfid.write(ss+str(ll)+ss+'-1'+ss+'1'+ss+'0'+ss+'0'+'\n')
    
    modfid.close()
    print 'Wrote Model file: ',modfn
    
    return modfn

def make1DInputFile(savepath,modelfile=None,datafile=None,roughtype=1,
                    maxiter=100,targetrms=1.0,rhostart=100,description='1dInv',
                    lagrange=5.0,roughness=1.0E7,debuglevel=1,iteration=0,
                    misfit=100.0):
    """
    Make a 1D input file for occam
    
    Input:
        savepath = full path to save input file to, if just path then saved as
                   savepath/input.input
        modelfile = full path to model file, if None then assumed to be in 
                    savepath/model.mod
        datafile = full path to data file, if None then assumed to be in 
                    savepath/data.data
        roughtype = roughness type
        maxiter = maximum number of iterations
        targetrms = target rms value
        rhostart = starting resistivity value on linear scale
        paramcount = 
    """
    
    
    ss='   '
    
    #make input data file name
    if savepath.find('.')==-1:
        if not os.path.exists(savepath):
            os.mkdir(savepath)
        inputfn=os.path.join(savepath,'input.txt')
        ipath=savepath
    else:
        inputfn=savepath
        ipath=os.path.dirname(savepath)
        
    if modelfile==None:
        modelfile=os.path.join(ipath,'model.txt')
    if datafile==None:
        datafile=os.path.join(ipath,'Data.txt')
    
    mdict=read1DModelFile(modelfile)
    paramcount=mdict['nparam']        
    
    infid=open(os.path.join(savepath,inputfn),'w')
    infid.write('Format:             OCCAMITER_FLEX      ! Flexible format \n')
    infid.write('Description:        '+description+'     !For your own notes. \n')
    infid.write('Model File:         '+modelfile+'       \n')
    infid.write('Data File:          '+datafile+'        \n')                                                                     
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
    print 'Wrote Input File: ',os.path.join(savepath,inputfn)
    
def read1DModelFile(modelfile):
    """
    will read in model file
    
    Inputs:
        modelfile = full path to model file
        
    Outputs:
        mdict = dictionary of values with keys:
            depth = depth of model in meters
            res = value of resisitivity
            pen = penalty
            pre = preference
            prefpen = preference penalty
            
    """

    mfid=open(modelfile,'r')
    mlines=mfid.readlines()
    mfid.close()
    
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
    for key in ['depth','res','pen','pref','prefpen']:
            mdict[key]=np.array(mdict[key])
                
    return mdict
    
def read1DInputFile(inputfile):
    """
    reads in a 1D input file
    """

    infid=open(inputfile,'r')
    ilines=infid.readlines()
    infid.close()

    indict={}
    indict['res']=[]
    
    for iline in ilines:
        if iline.find(':')>=0:
            ikey=iline[0:20].strip()
            ivalue=iline[20:].split('!')[0].strip()
            indict[ikey]=ivalue
        else:
            try:
                indict['res'].append(float(iline.strip()))
            except ValueError:
                pass
    return indict

def read1DdataFile(datafile):
    """
    reads a 1D data file
    
    Iputs:
        datafile = full path to data file
    """            
    
    dfid=open(datafile,'r')
    
    dlines=dfid.readlines()
    dfid.close()
    
    #get format of input data
    fmt=dlines[0].strip().split(':')[1].strip()
    
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
    rpdict={'freq':freq,
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
                if dlst[0]=='RhoZxy' or dlst[0]=='103':
                    rpdict['resxy'][0,jj]=dvalue
                    rpdict['resxy'][1,jj]=derr
                if dlst[0]=='PhsZxy' or dlst[0]=='104':
                    rpdict['phasexy'][0,jj]=dvalue
                    rpdict['phasexy'][1,jj]=derr
                if dlst[0]=='RhoZyx' or dlst[0]=='105':
                    rpdict['resyx'][0,jj]=dvalue
                    rpdict['resyx'][1,jj]=derr
                if dlst[0]=='PhsZyx' or dlst[0]=='106':
                    rpdict['phaseyx'][0,jj]=dvalue
                    rpdict['phaseyx'][1,jj]=derr
    
    return rpdict

def read1DIterFile(iterfile,modelfile):
    """
    read an iteration file
    
    Input:
        iterfile = full path to iteration file
        modelfile = full path to model file
    
    Output:
        idict = dictionary with keys:
            depth = depth of model in m
            modelres = resistivity of model, any fixed parameters will be put in
                        correct place.
            
    """
    
    mdict=read1DModelFile(modelfile)
    
    freeparams=np.where(mdict['res']==-1)[0]
    
    ifid=open(iterfile,'r')
    ilines=ifid.readlines()
    ifid.close()
    
    idict={}
    model=[]    
    for ii,iline in enumerate(ilines):
        if iline.find(':')>=0:
            ikey=iline[0:20].strip()
            ivalue=iline[20:].split('!')[0].strip()
            idict[ikey]=ivalue
        else:
            try:
                ilst=iline.strip().split()
                for kk in ilst:
                    model.append(float(kk))
            except ValueError:
                pass
            
    model=np.array(model)
    idict['modelres']=mdict['res'].copy()
    idict['modelres'][freeparams]=model
    idict['depth']=mdict['depth'].copy()
    
    return idict
    

def read1DRespFile(respfile):
    """
    read response file
    
    Input:
        repsfile = full path to response file
        
    Outputs:
        rpdict = dictionary with keys:
            freq = frequency array
            resxy,resyx,phasexy,phaseyx = array for corresponding component with
                index as:
                    0 = input data
                    1 = input error
                    2 = model data
                    3 = model error
    """
           
    
    dfid=open(respfile,'r')
    
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
    rpdict={'freq':freq,
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
                    rpdict['resxy'][0,jj]=dvalue
                    rpdict['resxy'][1,jj]=derr
                    rpdict['resxy'][2,jj]=rvalue
                    rpdict['resxy'][3,jj]=rerr
                if dlst[0]=='PhsZxy' or dlst[0]=='104':
                    rpdict['phasexy'][0,jj]=dvalue
                    rpdict['phasexy'][1,jj]=derr
                    rpdict['phasexy'][2,jj]=rvalue
                    rpdict['phasexy'][3,jj]=rerr
                if dlst[0]=='RhoZyx' or dlst[0]=='105':
                    rpdict['resyx'][0,jj]=dvalue
                    rpdict['resyx'][1,jj]=derr
                    rpdict['resyx'][2,jj]=rvalue
                    rpdict['resyx'][3,jj]=rerr
                if dlst[0]=='PhsZyx' or dlst[0]=='106':
                    rpdict['phaseyx'][0,jj]=dvalue
                    rpdict['phaseyx'][1,jj]=derr
                    rpdict['phaseyx'][2,jj]=rvalue
                    rpdict['phaseyx'][3,jj]=rerr
    
    return rpdict
    
def plot1D(respfile,iterfile,modelfile,fignum=1,ms=4):
    """
    
    """
    
    #read in data
    rpdict=read1DRespFile(respfile)
    mdict=read1DIterFile(iterfile,modelfile)

    period=1/rpdict['freq']
    
    #make a grid of subplots
    gs=gridspec.GridSpec(6,5,hspace=.15,wspace=.35)
    
    #make a figure
    fig=plt.figure(fignum,[8,8])
    
    #plot resistivity
    axr=fig.add_subplot(gs[:4,:4])
    rxy=np.where(rpdict['resxy'][0]!=0)[0]
    ryx=np.where(rpdict['resyx'][0]!=0)[0]
    if len(rxy)!=0:
        r1=axr.loglog(period[rxy],rpdict['resxy'][0][rxy],
                      ls='None',marker='o',color='k',mfc='k',ms=ms)
        titlestr='$Z_{xy}$'
    if len(ryx)!=0:
        r1=axr.loglog(period[ryx],rpdict['resyx'][0][ryx],
                      ls='None',marker='o',color='k',mfc='k',ms=ms)
        titlestr='$Z_{yx}$'
                      
    rxym=np.where(rpdict['resxy'][2]!=0)[0]
    ryxm=np.where(rpdict['resyx'][2]!=0)[0]
    if len(rxym)!=0:
        r2=axr.loglog(period[rxym],rpdict['resxy'][2][rxym],
                      ls=':',color='b',lw=2)
    if len(ryxm)!=0:
        r2=axr.loglog(period[ryxm],rpdict['resyx'][2][ryxm],
                      ls=':',color='b',lw=2)
                      
    axr.legend([r1,r2],['Data','Model'],loc='upper left',markerscale=2,
               borderaxespad=.05,
               labelspacing=.08,
               handletextpad=.15,borderpad=.05)
    
    #plot Phase
    axp=fig.add_subplot(gs[4:,:4],sharex=axr)
    pxy=np.where(rpdict['phasexy'][0]!=0)[0]
    pyx=np.where(rpdict['phaseyx'][0]!=0)[0]
    if len(pxy)!=0:
        p1=axp.semilogx(period[pxy],rpdict['phasexy'][0][pxy],
                      ls='None',marker='o',color='k',mfc='k',ms=ms)
    if len(pyx)!=0:
        p1=axp.semilogx(period[pyx],rpdict['phaseyx'][0][pyx],
                      ls='None',marker='o',color='k',mfc='k',ms=ms)
                      
    pxym=np.where(rpdict['phasexy'][2]!=0)[0]
    pyxm=np.where(rpdict['phaseyx'][2]!=0)[0]
    if len(pxym)!=0:
        p2=axp.semilogx(period[pxym],rpdict['phasexy'][2][pxym],
                      ls=':',color='b')
    if len(pyxm)!=0:
        p2=axp.semilogx(period[pyxm],rpdict['phaseyx'][2][pyxm],
                      ls=':',color='b')
                      
    axr.grid(True)
    axr.set_xticklabels(['' for ii in range(10)])
    axp.grid(True)
    axp.yaxis.set_major_locator(MultipleLocator(10))
    axp.yaxis.set_minor_locator(MultipleLocator(1))
    
    axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                   fontdict={'size':12,'weight':'bold'})
    axp.set_ylabel('Phase (deg)',
                   fontdict={'size':12,'weight':'bold'})
    axp.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})
    axr.yaxis.set_label_coords(-.15,.5)
    axp.yaxis.set_label_coords(-.15,.5)
    plt.suptitle(titlestr,fontsize=14,fontweight='bold')
    
    #plot 1D inversion
    axm=fig.add_subplot(gs[:,4])
    modelresp=abs(mdict['modelres'][1:])
    depthp=mdict['depth'][1:]
    axm.loglog(modelresp[::-1],depthp[::-1])
    axm.set_ylim(ymin=depthp[-1],ymax=depthp[1])
#    axm.set_xlim(xmax=10**3)
    axm.set_ylabel('Depth (m)',fontdict={'size':12,'weight':'bold'})
    axm.set_xlabel('Resistivity ($\Omega \cdot m$)',
                   fontdict={'size':12,'weight':'bold'})
#    axm.yaxis.set_minor_locator(MultipleLocator(200))
#    axm.yaxis.set_major_locator(MultipleLocator(1000))
    axm.grid(True,which='both')
    
    plt.show()
    

def make2DdataFile(edipath,mmode='both',savepath=None,stationlst=None,title=None,
                   thetar=0,resxyerr=10,resyxerr=10,phasexyerr=10,phaseyxerr=10,
                   ss=3*' ',fmt='%2.6f',freqstep=1,plotyn='y',lineori='ew',
                   tippererr=None):
    """
    make2DdataFile will make a data file for occam2D.  
    
    Input:
        edipath = path to edifiles
        mmode = modes to invert for.  Can be: 
                'both' -> will model both TE and TM modes
                'TM'   -> will model just TM mode
                'TE'   -> will model just TE mode
        savepath = path to save the data file to, this can include the name of
                   the data file, if not the file will be named:
                       savepath\Data.dat or edipath\Data.dat if savepath=None
        stationlst = list of stations to put in the data file, doesn't need to
                     be in order, the relative distance will be calculated
                     internally.  If stationlst=None, it will be assumed all the
                     files in edipath will be input into the data file
        title = title input into the data file
        thetar = rotation angle (deg) of the edifiles if you want to align the
                 components with the profile.  Angle is on the unit circle with 
                 an orientation that north is 0 degree, east -90.
        resxyerr = percent error in the res_xy component (TE), 
                  can be entered as 'data' where the errors from the data are
                  used.  
        resyxerr = percent error in the res_yx component (TM), 
                  can be entered as 'data' where the errors from the data are
                  used.  
        phasexyerr = percent error in the phase_xy component (TE), 
                  can be entered as 'data' where the errors from the data are
                  used.  
        phaseyxerr = percent error in the phase_yx component (TM), 
                  can be entered as 'data' where the errors from the data are
                  used.  
        ss = is the spacing parameter for the data file
        fmt = format of the numbers for the data file, see string formats for 
              a full description
        freqstep = take frequencies at this step, so if you want to take every
                   third frequency enter 3.  
                   Can input as a list of specific frequencies.  Note that the
                   frequencies must match the frequencies in the EDI files,
                   otherwise they will not be input.  
        plotyn = y or n to plot the stations on the profile line.
        lineori = predominant line orientation with respect to geographic north
                 ew for east-west line-> will orientate so first station is 
                                         farthest to the west
                 ns for north-south line-> will orientate so first station is 
                                         farthest to the south
        tippererr = error for tipper in percent.  If this value is entered than
                    the tipper will be included in the inversion, if the value
                    is None than the tipper will not be included.
              
    Output:
        datfilename = full path of data file
                 
    """
    
    if abs(thetar)>2*np.pi:
        thetar=thetar*(np.pi/180)
    #create rotation matrix
    rotmatrix=np.array([[np.cos(thetar), np.sin(thetar)],
                         [-np.sin(thetar), np.cos(thetar)]])
    
    #-----------------------Station Locations-----------------------------------    
    #create a list to put all the station dictionaries into
    surveylst=[]
    eastlst=[]
    northlst=[]
    pstationlst=[]
    freqlst=[]
    
    if stationlst==None:
        stationlst=[edifile[:-4] 
            for edifile in os.listdir(edipath) if edifile.find('.edi')]
    
    for kk,station in enumerate(stationlst):
        #search for filenames in the given directory and match to station name
        for filename in os.listdir(edipath):
            if fnmatch.fnmatch(filename,station+'*.edi'):
                print 'Found station edifile: ', filename
                surveydict={} #create a dictionary for the station data and info
                edifile=os.path.join(edipath,filename) #create filename path
                edidict=mt.readedi(edifile) #read in edifile as a dictionary
                freq=edidict['frequency']
                #check to see if the frequency is in descending order
                if freq[0]<freq[-1]:
                    freq=freq[::-1]
                    z=edidict['z'][::-1,:,:]
                    zvar=edidict['zvar'][::-1,:,:]
                    tip=edidict['tipper'][::-1,:,:]
                    tipvar=edidict['tippervar'][::-1,:,:]
                    print 'Flipped to descending frequency for station '+station
                else:
                    z=edidict['z']
                    zvar=edidict['zvar']
                    tip=edidict['tipper']
                    tipvar=edidict['tippervar']
                #rotate matrices if angle is greater than 0
                if thetar!=0:
                    for rr in range(len(z)):
                        z[rr,:,:]=np.dot(rotmatrix,np.dot(z[rr],rotmatrix.T))
                        zvar[rr,:,:]=np.dot(rotmatrix,np.dot(zvar[rr],
                                                             rotmatrix.T))
                else:
                    pass
                        
                zone,east,north=utm2ll.LLtoUTM(23,edidict['lat'],edidict['lon'])
                #put things into a dictionary to sort out order of stations
                surveydict['station']=station
                surveydict['east']=east
                surveydict['north']=north
                surveydict['zone']=zone
                surveydict['z']=z
                surveydict['zvar']=zvar
                surveydict['freq']=freq
                surveydict['tipper']=tip
                surveydict['tippervar']=tipvar
                surveydict['lat']=edidict['lat']
                surveydict['lon']=edidict['lon']
                freqlst.append(freq)
                eastlst.append(east)
                northlst.append(north)
                pstationlst.append(station)
                surveylst.append(surveydict)
                
    #project stations onto a best fitting line
    #plot a bestfitting line
    p=sp.polyfit(eastlst,northlst,1)
    theta=np.arctan(p[0])
    print 'Profile Line Angle is: {0:.4g} (E=0,N=90)'.format(theta*180/np.pi)
    
    #plot stations on profile line
    if plotyn=='y':
        plt.figure(4)
        plt.title('Projected Stations')
        plt.plot(eastlst,sp.polyval(p,eastlst),'-b',lw=2)
        
    for ii in range(len(surveylst)):
        d=(northlst[ii]-sp.polyval(p,eastlst[ii]))*np.cos(theta)
        x0=eastlst[ii]+d*np.sin(theta)
        y0=northlst[ii]-d*np.cos(theta)
        surveylst[ii]['east']=x0
        surveylst[ii]['north']=y0
        if plotyn=='y':
            plt.plot(x0,y0,'v',color='k',ms=8,mew=3)
            plt.text(x0,y0+y0*.0005,pstationlst[ii],horizontalalignment='center',
                 verticalalignment='baseline',fontdict={'size':12,
                                                        'weight':'bold'})
        
        #need to figure out a way to account for zone changes
        
        if lineori=='ew': 
            if surveylst[0]['east']<surveylst[ii]['east']:
                surveylst[ii]['offset']=np.sqrt((surveylst[0]['east']-
                                                surveylst[ii]['east'])**2+
                                                (surveylst[0]['north']-
                                                surveylst[ii]['north'])**2)
            elif surveylst[0]['east']>surveylst[ii]['east']:
                surveylst[ii]['offset']=-1*np.sqrt((surveylst[0]['east']-
                                                surveylst[ii]['east'])**2+
                                                (surveylst[0]['north']-
                                                surveylst[ii]['north'])**2)
            else:
                surveylst[ii]['offset']=0
        elif lineori=='ns': 
            if surveylst[0]['north']<surveylst[ii]['north']:
                surveylst[ii]['offset']=np.sqrt((surveylst[0]['east']-
                                                surveylst[ii]['east'])**2+
                                                (surveylst[0]['north']-
                                                surveylst[ii]['north'])**2)
            elif surveylst[0]['north']>surveylst[ii]['north']:
                surveylst[ii]['offset']=-1*np.sqrt((surveylst[0]['east']-
                                                surveylst[ii]['east'])**2+
                                                (surveylst[0]['north']-
                                                surveylst[ii]['north'])**2)
            else:
                surveylst[ii]['offset']=0
    
    #sort by ascending order of distance from first station
    surveylst=sorted(surveylst,key=itemgetter('offset'))
    
    #number of stations read    
    nstat=len(surveylst)    
    
    #--------------------------Match Frequencies--------------------------------
    #a dictionary is created with the frequency as the key and the value is the
    #frequency number in the list. Each edi file is iterated over extracting
    #only the matched frequencies.  This makes it necessary to have the same
    #frequency content in each edifile.    
    
    #make a list to iterate over frequencies
    if type(freqstep) is list or type(freqstep) is not int:
        if type(freqstep[0]) is int:
            #find the median frequency list
            maxflen=max([len(ff) for ff in freqlst])
            farray=np.zeros((nstat,maxflen))
            for ii in range(nstat):
                farray[ii,0:len(freqlst[ii])]=freqlst[ii]
        
            mfreq=np.median(farray,axis=0)
            print len(mfreq),len(freqstep)
            fdict=dict([('%.6g' % mfreq[ff],ii) 
                            for ii,ff in enumerate(freqstep,1) if mfreq[ff]!=0])
        else:
            fdict=dict([('%.6g' % ff,ii) for ii,ff in enumerate(freqstep,1)])
    else:
        #find the median frequency list
        maxflen=max([len(ff) for ff in freqlst])
        farray=np.zeros((nstat,maxflen))
        for ii in range(nstat):
            farray[ii,0:len(freqlst[ii])]=freqlst[ii]
        
        mfreq=np.median(farray,axis=0)
    
        #make a dictionary of values        
        fdict=dict([('%.6g' % ff,ii) for ii,ff in 
                    enumerate(mfreq[range(0,maxflen,freqstep)],1) if ff!=0])

    #print the frequencies to look for to make sure its what the user wants
    #make a list of keys that is sorted in descending order
    klst=[float(dd) for dd in fdict.keys()]
    klst.sort(reverse=True)
    klst=['%.6g' % dd for dd in klst]    
    
    print 'Frequencies to look for are: (# freq(Hz) Period(s)) '
    for key in klst:
        print fdict[key],key, 1./float(key)
    
    #make lists of parameters to write to file    
    reslst=[]
    offsetlst=[]
    stationlstsort=[]
    for kk in range(nstat):
        z=surveylst[kk]['z']
        zvar=surveylst[kk]['zvar']
        freq=surveylst[kk]['freq']
        offsetlst.append(surveylst[kk]['offset'])  
        stationlstsort.append(surveylst[kk]['station'])
        tip=surveylst[kk]['tipper']
        tipvar=surveylst[kk]['tippervar']
        #loop over frequencies to pick out the ones desired
        for jj,ff in enumerate(freq):
            #jj is the index of edi file frequency list, this index corresponds
            #to the impedance tensor component index
            #ff is the frequency from the edi file frequency list
            try:
                #nn is the frequency number out of extracted frequency list
                nn=fdict['%.6g' % ff]
                
                #calculate resistivity 
                wt=.2/(ff)
                resxy=wt*abs(z[jj,0,1])**2
                resyx=wt*abs(z[jj,1,0])**2
        
                #calculate the phase putting the yx in the 1st quadrant        
                phasexy=np.arctan2(z[jj,0,1].imag,z[jj,0,1].real)*(180/np.pi)
                phaseyx=np.arctan2(z[jj,1,0].imag,z[jj,1,0].real)*(180/np.pi)+\
                        180
                #put phases in correct quadrant if should be negative
                if phaseyx>180:
                    phaseyx=phaseyx-360
                    print 'Found Negative Phase',surveylst[kk]['station'],ff    
                
                #calculate errors
                #res_xy (TE)
                if resxyerr=='data':
                    dresxyerr=wt*(abs(z[jj,0,1])+zvar[jj,0,1])**2-resxy
                    lresxyerr=(dresxyerr/resxy)/np.log(10)
                
                else:
                    lresxyerr=(resxyerr/100.)/np.log(10)
                
                #Res_yx(TM)
                if resyxerr=='data':
                    dresyxerr=wt*(abs(z[jj,1,0])+zvar[jj,1,0])**2-resyx
                    lresyxerr=(dresyxerr/resyx)/np.log(10)
                else:
                    lresyxerr=(resyxerr/100.)/np.log(10)
                
                #phase_xy(TE)
                if phasexyerr=='data':
                    dphasexyerr=np.arcsin(zvar[jj,0,1]/abs(z[jj,0,1]))*\
                                (180/np.pi)
                else:
                    dphasexyerr=(phasexyerr/100.)*57/2.
                    
                #phase_yx (TM)
                if phaseyxerr=='data':
                    dphaseyxerr=np.arcsin(zvar[jj,1,0]/abs(z[jj,1,0]))*\
                                (180/np.pi)
                else:
                    dphaseyxerr=(phaseyxerr/100.)*57/2.
                
                #calculate log10 of resistivity as prescribed by OCCAM
                lresyx=np.log10(resyx)
                lresxy=np.log10(resxy)
                
                #if include the tipper
                if tippererr!=None:
                    if tip[jj,0].real==0.0 or tip[jj,1]==0.0:
                        tipyn='n'
                    else:
                        #calculate the projection angle for real and imaginary
                        tipphir=np.arctan(tip[jj,0].real/tip[jj,1].real)-theta
                        tipphii=np.arctan(tip[jj,0].imag/tip[jj,1].imag)-theta
                        
                        #project the tipper onto the profile line
                        projtipr=np.sqrt(tip[jj,0].real**2+tip[jj,1].real**2)*\
                                  np.cos(tipphir)
                        projtipi=np.sqrt(tip[jj,0].imag**2+tip[jj,1].imag**2)*\
                                  np.cos(tipphii)
                                  
                        #error of tipper is a decimal percentage
                        projtiperr=tippererr/100.
                        
                        tipyn='y'
                    
                
                #make a list of lines to write to the data file
                if mmode=='both':
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'1'+ss+
                                    fmt % lresxy +ss+fmt % lresxyerr+'\n')
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'2'+ss+
                                    fmt % phasexy +ss+fmt % dphasexyerr+'\n')
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'5'+ss+
                                    fmt % lresyx+ss+fmt % lresyxerr+'\n')
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'6'+ss+
                                    fmt % phaseyx +ss+fmt % dphaseyxerr+'\n')
                    if tippererr!=None and tipyn=='y':
                        reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'3'+ss+
                                    fmt % projtipr +ss+fmt % projtiperr+'\n')
                        reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'4'+ss+
                                    fmt % projtipi +ss+fmt % projtiperr+'\n')
                elif mmode=='TM':
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'5'+ss+
                                    fmt % lresyx +ss+fmt % lresyxerr+'\n')
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'6'+ss+
                                    fmt % phaseyx +ss+fmt % dphaseyxerr+'\n')
                    if tippererr!=None and tipyn=='y':
                        reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'3'+ss+
                                    fmt % projtipr +ss+fmt % projtiperr+'\n')
                        reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'4'+ss+
                                    fmt % projtipi +ss+fmt % projtiperr+'\n')
                elif mmode=='TE':
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'1'+ss+
                                    fmt % lresxy+ss+fmt % lresxyerr+'\n')
                    reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'2'+ss+
                                    fmt % phasexy+ss+fmt % dphasexyerr+'\n')
                    if tippererr!=None and tipyn=='y':
                        reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'3'+ss+
                                    fmt % projtipr +ss+fmt % projtiperr+'\n')
                        reslst.append(ss+str(kk+1)+ss+str(nn)+ss+'4'+ss+
                                    fmt % projtipi +ss+fmt % projtiperr+'\n')
                else:
                    raise NameError('mmode' +mmode+' not defined')
            except KeyError:
                pass
    
    #===========================================================================
    #                             write dat file
    #===========================================================================
    
    if savepath!=None:
        if savepath.find('.')>0:
            datfilename=savepath
        else:
            if not os.path.exists(savepath):
                os.mkdir(savepath)
            datfilename=os.path.join(savepath,'Data.dat')
    else:
        datfilename=os.path.join(edipath,'Data.dat')
        
    if title==None:
        title='Occam Inversion'
        
    datfid=open(datfilename,'w')
    datfid.write('FORMAT:'+' '*11+'OCCAM2MTDATA_1.0'+'\n')
    datfid.write('TITLE:'+' '*12+'{0:.4g}--'.format(theta*180/np.pi)+' '+title+'\n')
    
    #write station sites
    datfid.write('SITES:'+' '*12+str(nstat)+'\n')
    for station in stationlstsort:
        datfid.write(ss+station+'\n')
    
    #write offsets
    datfid.write('OFFSETS (M):'+'\n')
    for offset in offsetlst:
        datfid.write(ss+fmt % offset+'\n')
    
    #write frequencies
    #writefreq=[freq[ff] for ff in range(0,len(freq),freqstep)]
    datfid.write('FREQUENCIES:'+' '*8+str(len(fdict))+'\n')
    for fkey in klst:
        datfid.write(ss+fmt % float(fkey) +'\n')
    
    #write data block
    datfid.write('DATA BLOCKS:'+' '*10+str(len(reslst))+'\n')
    datfid.write('SITE'+ss+'FREQ'+ss+'TYPE'+ss+'DATUM'+ss+'ERROR'+'\n')
    for ll,datline in enumerate(reslst):
        if datline.find('#IND')>=0:
            print 'Found #IND on line ',ll
            ndline=datline.replace('#IND','00')
            print 'Replaced with 00'
            datfid.write(ndline)
        else:
            datfid.write(datline)
    datfid.close()
    
    print 'Wrote Occam2D data file to: ',datfilename
    
    return datfilename
    
def makeModel(datafilename,niter=20,targetrms=1.0,nlayers=100,nlperdec=30,
              z1layer=50,bwidth=200,trigger=.75,savepath=None,rhostart=100,
              occampath=r"c:\Peacock\PHD\OCCAM\MakeFiles"):
    """
    makeModel will make an the input files for occam using Steve Constable's
    MakeModel2DMT.f code.
    
    Inputs:
        datafilename = full path to data file
        niter = maximum number of iterations
        targetrms = target root mean square error
        nlayers = total number of layers in mesh
        nlperdec = number of layers per decade
        z1layer = thickness of the first layer in meters
        bwidth = maximum block width for regularization grid in meters
        trigger = triger point to amalgamate blocks
        savepath = path to save files to
        rhostart = starting resistivity for homogeneous half space in ohm-m
        occampath = path to MakeModel2DMT.exe
        
    Outputs:
        meshfn = mesh file for finite element grid saved ats MESH
        inmodelfn = input model, starting model with rhostart as starting value
                    saved as INMODEL
        startupfn = start up filepath, saved as startup
    """
    #get the base name of data file    
    dfnb=os.path.basename(datafilename)
    
    #put data file into the same directory as MakeModel2DMT
    if os.path.dirname(datafilename)!=occampath:
        shutil.copy(datafilename,os.path.join(occampath,dfnb))
    
    #write input file for MakeModel2DMT
    mmfid=open(os.path.join(occampath,'inputMakeModel.txt'),'w')
    mmfid.write(dfnb+'\n')
    mmfid.write(str(niter)+'\n')    
    mmfid.write(str(targetrms)+'\n')    
    mmfid.write(str(nlayers)+'\n')
    mmfid.write(str(nlperdec)+'\n')
    mmfid.write(str(z1layer)+'\n')
    mmfid.write(str(bwidth)+'\n')
    mmfid.write(str(trigger)+'\n')
    mmfid.write('\n')
    mmfid.close()
    
    shutil.copy(os.path.join(occampath,'inputMakeModel.txt'),
                os.path.join(os.path.dirname(datafilename),
                             'inputMakeModel.txt'))
    #get current working directory
    cdir=os.getcwd() 
    
    #change directory path to occam path
    os.chdir(occampath) 
    
    #---call MakeModel2DMT---
    subprocess.os.system("MakeModel2DMT < inputMakeModel.txt")
    
    #change back to original working directory    
    os.chdir(cdir)
    
    if savepath==None:
        savepath=os.path.dirname(datafilename)
    
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    
    meshfn=os.path.join(savepath,'MESH')    
    inmodelfn=os.path.join(savepath,'INMODEL')    
    startupfn=os.path.join(savepath,'startup')    
    
    #copy ouput files to savepath
    shutil.copy(os.path.join(occampath,'MESH'),meshfn)
    shutil.copy(os.path.join(occampath,'INMODEL'),inmodelfn)
    shutil.copy(os.path.join(occampath,'startup'),startupfn)
    if not os.path.exists(os.path.join(savepath,dfnb)):
        shutil.copy(datafilename,os.path.join(savepath,dfnb))
    if os.path.getctime(os.path.join(savepath,dfnb))<\
        os.path.getctime(datafilename):
        shutil.copy(datafilename,os.path.join(savepath,dfnb))
    
    #rewrite mesh so it contains the right number of columns and rows
    rewriteMesh(meshfn)
    
    #write startup file to have the starting desired starting rho value
    ifid=open(startupfn,'r')
    ilines=ifid.readlines()
    ifid.close()
    
    if rhostart!=100:
        #make startup model a homogeneous half space of rhostart
        rhostart=np.log10(rhostart)
        ifid=open(startupfn,'w')
        for line in ilines:
            if line.find('2.000000')>=0:
                line=line.replace('2.000000','%.6f' % rhostart)
            ifid.write(line)
    ifid.close()
    
    print 'Be sure to check the INMODEL file for clumped numbers near the bottom.'
    print 'Also, check the MESH and startup files to make sure they are correct.'
    
    return meshfn,inmodelfn,startupfn


def rewriteMesh(meshfn):
    """
    checkMesh will check to see if the number of lines are correct and the 
    length of the line is correct
    """
    
    #check the mesh
    mfid=open(meshfn,'r')
    mlines=mfid.readlines()
    mfid.close()
    
    #get parameters for number of lines and number of unknowns per line    
    pstr=mlines[1].strip().split()
    nhn=int(pstr[1])
    nvn=(int(pstr[2])-1)*4
    
    print 'Number of horizontal nodes: ',nhn
    print 'Number of vertical nodes: ',nvn
    #find the first line
    for ii,mline in enumerate(mlines[2:],2):
        if mline.find('?')==0:
            qspot=ii
            break

    #rewrite the file to have the proper amount of stuff
    mfid=open(meshfn,'w')
    for line in mlines[0:qspot]:
        mfid.write(line)
    
    for kk in range(qspot,qspot+nvn):
        mfid.write('?'*(nhn-1)+'\n')
    mfid.close()
 

def read2Dmesh(meshfile,ncol=8):
    """
    read a 2D meshfile
    """
    
    mfid=file(meshfile,'r')
    
    mlines=mfid.readlines()
    
    nm=len(mlines)
    
    nh=int(mlines[1].strip().split()[1])-1
    nv=int(mlines[1].strip().split()[2])-1
    
    hnodes=np.zeros(nh)
    vnodes=np.zeros(nv)
    mdata=np.zeros((nh,nv,4),dtype=str)    
    
    
    for ii,mm in enumerate(mlines[2:nh/ncol+2]):
        hline=mm.strip().split()
        for jj,hh in enumerate(hline):
            try:
                hnodes[ii*8+jj]=float(hh)
            except IndexError:
                pass
    
    for ii,mm in enumerate(mlines[nh/ncol+2:nv/ncol+nh/ncol+3]):
        vline=mm.strip().split()
        for jj,hh in enumerate(vline):
            try:
                vnodes[ii*8+jj]=float(hh)
            except IndexError:
                pass
            
    for ii,mm in enumerate(mlines[nv/ncol+nh/ncol+4:]):
        kk=0
        while kk<4:        
            mline=mm.rstrip()
            for jj in range(nh):
                try:
                    mdata[jj,ii,kk]=mline[jj]
                except IndexError:
                    pass
            kk+=1
    
    return hnodes,vnodes,mdata
    
    

def read2DdataFile(datafilename):
    """
    read2DdataFile will read in data from a 2D occam data file.  Only supports
    the first 6 data types of occam2D
    
    Input: 
        datafilename = full path to data file
    
    Output:
        rplst = list of dictionaries for each station with keywords:
            'station' = station name
            'offset' = relative offset,
            'resxy' = TE resistivity and error as row 0 and 1 ressectively,
            'resyx'= TM resistivity and error as row 0 and 1 respectively,
            'phasexy'= TE phase and error as row 0 and 1 respectively,
            'phaseyx'= Tm phase and error as row 0 and 1 respectively,
            'realtip'= Real Tipper and error as row 0 and 1 respectively,
            'imagtip'= Imaginary Tipper and error as row 0 and 1 respectively
            
            Note: that the resistivity will be in log10 space.  Also, there are
            2 extra rows in the data arrays, this is to put the response from
            the inversion. 
        
        stationlst = list of stations in order from one side of the profile
                     to the other.
        freq = list of frequencies used in the inversion
        title = title, could be useful for plotting.
    """
    
    dfid=open(datafilename,'r')
    
    dlines=dfid.readlines()
    #get format of input data
    fmt=dlines[0].strip().split(':')[1].strip()
    
    #get title
    titlestr=dlines[1].strip().split(':')[1].strip()

    if titlestr.find('--')>0:
        tstr=titlestr.split('--')
        theta=tstr[0]
        title=tstr[1]
    else:
        title=titlestr
        theta=0
        print 'Need to figure out angle of profile line'
    #get number of sits
    nsites=int(dlines[2].strip().split(':')[1].strip())
    
    #get station names
    stationlst=[dlines[ii].strip() for ii in range(3,nsites+3)]
    
    #get offsets in meters
    offsets=[float(dlines[ii].strip()) for ii in range(4+nsites,4+2*nsites)]
    
    #get number of frequencies
    nfreq=int(dlines[4+2*nsites].strip().split(':')[1].strip())

    #get frequencies
    freq=[float(dlines[ii].strip()) for ii in range(5+2*nsites,
                                                      5+2*nsites+nfreq)]
                                                      

    #-----------get data-------------------
    #set zero array size the first row will be the data and second the error
    asize=(4,nfreq)
    #make a list of dictionaries for each station.
    rplst=[{'station':station,'offset':offsets[ii],
            'resxy':np.zeros(asize),
            'resyx':np.zeros(asize),
            'phasexy':np.zeros(asize),
            'phaseyx':np.zeros(asize),
            'realtip':np.zeros(asize),
            'imagtip':np.zeros(asize),
            } for ii,station in enumerate(stationlst)]
    for line in dlines[7+2*nsites+nfreq:]:
        ls=line.split()
        #station index
        ss=int(float(ls[0]))-1
        #component key
        comp=str(int(float(ls[2])))
        #frequency index        
        ff=int(float(ls[1]))-1
        #print ls,ss,comp,ff
        #put into array
        #input data
        rplst[ss][occamdict[comp]][0,ff]=float(ls[3]) 
        #error       
        rplst[ss][occamdict[comp]][1,ff]=float(ls[4])
    
    return rplst,stationlst,np.array(freq),title,theta
    
def rewrite2DdataFile(datafilename,edipath=None,thetar=0,resxyerr=None,
                      resyxerr=None,phasexyerr=None,phaseyxerr=None,
                      tippererr=None,mmode='both',flst=None,removestation=None):
    """
    rewrite2DDataFile will rewrite an existing data file so you can redefine 
    some of the parameters, such as rotation angle, or errors for the different
    components or only invert for one mode or add one or add tipper or remove
    tipper.
    
    Inputs:
        datafilename = full path to data file to rewrite
        thetar = rotation angle with positive clockwise
        resxyerr = error for TE mode resistivity (percent) or 'data' for data
        resyxerr = error for TM mode resistivity (percent) or 'data' for data
        phasexyerr = error for TE mode phase (percent) or 'data' for data
        phaseyxerr = error for TM mode phase (percent) or 'data' for data
        tippererr = error for tipper (percent) input only if you want to invert
                    for the tipper or 'data' for data errors
        mmodes = 'both' for both TE and TM
                 'TE' for TE
                 'TM' for TM
        flst = frequency list in Hz to rewrite, needs to be similar to the 
                datafile, cannot add frequencies
        removestation = list of stations to remove if desired
    """
    ss=3*' '
    fmt='%2.6f'
    
    #load the data for the data file    
    rplst,stationlst,freq,title,theta=read2DdataFile(datafilename)
    
    #make a dictionary of rplst for easier extraction of data
    rpdict=dict([(station,rplst[ii]) for ii,station in enumerate(stationlst)])

    #remove stations from rplst and stationlst if desired
    if removestation!=None:
        #if removestation is not a list make it one
        if type(removestation) is not list:
            removestation=[removestation]
        
#        #remove station dictionary from rplst
#        for rstation in removestation:
#            for hh,hdict in enumerate(rplst):
#                if hdict['station']==rstation:
#                    rplst.remove(rplst[hh])
        
        #remove station from station list           
        for rstation in removestation:        
            try:
                stationlst.remove(rstation)
            except ValueError:
                print 'Did not find '+rstation
    
    #if flst is not the same as freq make freq=flst
    if flst!=None:
        freq=flst
    
    
    #if the rotation angle is not 0 than need to read the original data in
    if thetar!=0:
        if edipath==None:
            raise IOError('Need to input the edipath to original edifiles to'+
                           ' get rotations correct')
        
        #get list of edifiles already in data file
        edilst=[os.path.join(edipath,edi) for stat in stationlst 
                for edi in os.listdir(edipath) if edi[0:len(stat)]==stat]
        reslst=[]
        for kk,edifn in enumerate(edilst,1):
            imp1=Z.Z(edifn)
            rp=imp1.getResPhase(thetar=thetar)
            imptip=imp1.getTipper()
            tip=imptip.tipper
            station=stationlst[kk-1]
            fdict=dict([('{0:.6g}'.format(fr),ii) for ii,fr in enumerate(imp1.frequency)])
            #loop over frequencies to pick out the ones desired
            for jj,ff in enumerate(freq,1):
                #jj is the index of edi file frequency list, this index corresponds
                #to the impedance tensor component index
                #ff is the frequency from the edi file frequency list
                try:
                    #nn is the frequency number out of extracted frequency list
                    nn=fdict['%.6g' % ff]
                    
                    #calculate resistivity
                    resxy=rp.resxy[nn]
                    resyx=rp.resyx[nn]
            
                    #calculate the phase putting the yx in the 1st quadrant
                    phasexy=rp.phasexy[nn]
                    phaseyx=rp.phaseyx[nn]+180
                    #put phases in correct quadrant if should be negative
                    if phaseyx>180:
                        phaseyx=phaseyx-360
                        print 'Found Negative Phase at',imp1.station,ff    
                    
                    #calculate errors
                    #res_xy (TE)
                    if resxyerr=='data':
                        lresxyerr=(rp.resxyerr[nn]/resxy)/np.log(10)
                    else:
                        lresxyerr=(resxyerr/100.)/np.log(10)
                    
                    #Res_yx(TM)
                    if resyxerr=='data':
                        lresyxerr=(rp.resyxerr[nn]/resyx)/np.log(10)
                    else:
                        lresyxerr=(resyxerr/100.)/np.log(10)
                    
                    #phase_xy(TE)
                    if phasexyerr=='data':
                        dphasexyerr=rp.phasexyerr[nn]
                    else:
                        dphasexyerr=(phasexyerr/100.)*57/2.
                        
                    #phase_yx (TM)
                    if phaseyxerr=='data':
                        dphaseyxerr=rp.phaseyxerr[nn]
                    else:
                        dphaseyxerr=(phaseyxerr/100.)*57/2.
                    
                    #calculate log10 of resistivity as prescribed by OCCAM
                    lresyx=np.log10(resyx)
                    lresxy=np.log10(resxy)
                    
                    #if include the tipper
                    if tippererr!=None:
                        if tip.tipper[nn,0]==0.0 or tip[nn,1]==0.0:
                            tipyn='n'
                        else:
                            #calculate the projection angle for real and imaginary
                            tipphir=np.arctan(tip[nn,0].real/tip[nn,1].real)-\
                                    theta
                            tipphii=np.arctan(tip[nn,0].imag/tip[nn,1].imag)-\
                                    theta
                            
                            #project the tipper onto the profile line
                            projtipr=np.sqrt(tip[nn,0].real**2+tip[nn,1].real**2)*\
                                      np.cos(tipphir)
                            projtipi=np.sqrt(tip[nn,0].imag**2+tip[nn,1].imag**2)*\
                                      np.cos(tipphii)
                                      
                            #error of tipper is a decimal percentage
                            projtiperr=tippererr/100.
                            
                            tipyn='y'
                        
                    
                    #make a list of lines to write to the data file
                    if mmode=='both':
                        if rpdict[station]['resxy'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'1'+ss+
                                        fmt % lresxy +ss+fmt % lresxyerr+'\n')
                        if rpdict[station]['phasexy'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'2'+ss+
                                        fmt % phasexy+ss+fmt % dphasexyerr+'\n')
                        if rpdict[station]['resyx'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'5'+ss+
                                        fmt % lresyx+ss+fmt % lresyxerr+'\n')
                        if rpdict[station]['phaseyx'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'6'+ss+
                                        fmt % phaseyx+ss+fmt % dphaseyxerr+'\n')
                        if tippererr!=None and tipyn=='y':
                            if rpdict[station]['realtip'][0,jj-1]!=0.0:
                                reslst.append(ss+str(kk)+ss+str(jj)+ss+'3'+ss+
                                            fmt % projtipr+ss+fmt % projtiperr+
                                            '\n')
                            if rpdict[station]['imagtip'][0,jj-1]!=0.0:
                                reslst.append(ss+str(kk)+ss+str(jj)+ss+'4'+ss+
                                            fmt % projtipi+ss+fmt % projtiperr+
                                            '\n')
                    elif mmode=='TM':
                        if rpdict[station]['resyx'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'5'+ss+
                                        fmt % lresyx +ss+fmt % lresyxerr+'\n')
                        if rpdict[station]['phaseyx'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'6'+ss+
                                        fmt % phaseyx+ss+fmt % dphaseyxerr+'\n')
                        if tippererr!=None and tipyn=='y':
                            if rpdict[station]['realtip'][0,jj-1]!=0.0:
                                reslst.append(ss+str(kk)+ss+str(jj)+ss+'3'+ss+
                                            fmt % projtipr+ss+fmt % projtiperr+
                                            '\n')
                            if rpdict[station]['imagtip'][0,jj-1]!=0.0:
                                reslst.append(ss+str(kk)+ss+str(jj)+ss+'4'+ss+
                                            fmt % projtipi+ss+fmt % projtiperr+
                                            '\n')
                    elif mmode=='TE':
                        if rpdict[station]['resxy'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'1'+ss+
                                        fmt % lresxy +ss+fmt % lresxyerr+'\n')
                        if rpdict[station]['phasexy'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'2'+ss+
                                        fmt % phasexy+ss+fmt % dphasexyerr+'\n')
                        if tippererr!=None and tipyn=='y':
                            if rpdict[station]['realtip'][0,jj-1]!=0.0:
                                reslst.append(ss+str(kk)+ss+str(jj)+ss+'3'+ss+
                                            fmt % projtipr+ss+fmt % projtiperr+
                                            '\n')
                            if rpdict[station]['imagtip'][0,jj-1]!=0.0:
                                reslst.append(ss+str(kk)+ss+str(jj)+ss+'4'+ss+
                                            fmt % projtipi+ss+fmt % projtiperr+
                                            '\n')
                    else:
                        raise NameError('mmode' +mmode+' not defined')
                except KeyError:
                    pass
    
    #If no rotation is desired but error bars are than...
    else:
        reslst=[]
        for kk,station in enumerate(stationlst,1):
            srp=rpdict[station]
            nr=srp['resxy'].shape[1]
            #calculate errors and rewrite
            #res_xy (TE)
            if resxyerr!=None:
                lresxyerr=(resxyerr/100.)/np.log(10)
                srp['resxy'][1,:]=np.repeat(lresxyerr,nr)
            
            #Res_yx(TM)
            if resyxerr!=None:
                lresyxerr=(resyxerr/100.)/np.log(10)
                srp['resyx'][1,:]=np.repeat(lresyxerr,nr)
            
            #phase_xy(TE)
            if phasexyerr!=None:
                dphasexyerr=(phasexyerr/100.)*57/2.
                srp['phasexy'][1,:]=np.repeat(dphasexyerr,nr)
                
            #phase_yx (TM)
            if phaseyxerr!=None:
                dphaseyxerr=(phaseyxerr/100.)*57/2.
                srp['phaseyx'][1,:]=np.repeat(dphaseyxerr,nr)
            
            if tippererr!=None:
                #error of tipper is a decimal percentage
                projtiperr=tippererr/100.
                srp['realtip'][1,:]=np.repeat(projtiperr,nr)
                srp['imagtip'][1,:]=np.repeat(projtiperr,nr)
            
            for jj,ff in enumerate(freq,1):
                #make a list of lines to write to the data file
                if mmode=='both':
                    if srp['resxy'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'1'+ss+
                                    fmt % srp['resxy'][0,jj-1]+ss+
                                    fmt % srp['resxy'][1,jj-1]+'\n')
                    if srp['phasexy'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'2'+ss+
                                    fmt % srp['phasexy'][0,jj-1]+ss+
                                    fmt % srp['phasexy'][1,jj-1]+'\n')
                    if srp['resyx'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'5'+ss+
                                    fmt % srp['resyx'][0,jj-1]+ss+
                                    fmt % srp['resyx'][1,jj-1]+'\n')
                    if srp['phaseyx'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'6'+ss+
                                    fmt % srp['phaseyx'][0,jj-1]+ss+
                                    fmt % srp['phaseyx'][1,jj-1]+'\n')
                    if tippererr!=None and tipyn=='y':
                        if srp['realtip'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'3'+ss+
                                    fmt % srp['realtip'][0,jj-1]+ss+
                                    fmt % srp['realtip'][1,jj-1]+'\n')
                        if srp['imagtip'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'4'+ss+
                                    fmt % srp['imagtip'][0,jj-1]+ss+
                                    fmt % srp['imagtip'][1,jj-1]+'\n')
                elif mmode=='TM':
                    if srp['resyx'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'5'+ss+
                                    fmt % srp['resyx'][0,jj-1]+ss+
                                    fmt % srp['resyx'][1,jj-1]+'\n')
                    if srp['phaseyx'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'6'+ss+
                                    fmt % srp['phaseyx'][0,jj-1]+ss+
                                    fmt % srp['phaseyx'][1,jj-1]+'\n')
                    if tippererr!=None and tipyn=='y':
                        if srp['realtip'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'3'+ss+
                                    fmt % srp['realtip'][0,jj-1]+ss+
                                    fmt % srp['realtip'][1,jj-1]+'\n')
                        if srp['imagtip'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'4'+ss+
                                    fmt % srp['imagtip'][0,jj-1]+ss+
                                    fmt % srp['imagtip'][1,jj-1]+'\n')
                elif mmode=='TE':
                    if srp['resxy'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'1'+ss+
                                    fmt % srp['resxy'][0,jj-1]+ss+
                                    fmt % srp['resxy'][1,jj-1]+'\n')
                    if srp['phasexy'][0,jj-1]!=0.0:
                        reslst.append(ss+str(kk)+ss+str(jj)+ss+'2'+ss+
                                    fmt % srp['phasexy'][0,jj-1]+ss+
                                    fmt % srp['phasexy'][1,jj-1]+'\n')
                    if tippererr!=None and tipyn=='y':
                        if srp['realtip'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'3'+ss+
                                    fmt % srp['realtip'][0,jj-1]+ss+
                                    fmt % srp['realtip'][1,jj-1]+'\n')
                        if srp['imagtip'][0,jj-1]!=0.0:
                            reslst.append(ss+str(kk)+ss+str(jj)+ss+'4'+ss+
                                    fmt % srp['imagtip'][0,jj-1]+ss+
                                    fmt % srp['imagtip'][1,jj-1]+'\n')

    #===========================================================================
    #                             write dat file
    #===========================================================================
    
    dfn=datafilename[:-4]+'RW.dat'
    nstat=len(stationlst)
        
    if title==None:
        title='Occam Inversion'
        
    datfid=open(dfn,'w')
    datfid.write('FORMAT:'+' '*11+'OCCAM2MTDATA_1.0'+'\n')
    datfid.write('TITLE:'+' '*12+title+'\n')
    
    #write station sites
    datfid.write('SITES:'+' '*12+str(nstat)+'\n')
    for station in stationlst:
        datfid.write(ss+station+'\n')
    
    #write offsets
    datfid.write('OFFSETS (M):'+'\n')
    for station in stationlst:
        datfid.write(ss+fmt % rpdict[station]['offset']+'\n')
    
    #write frequencies
    #writefreq=[freq[ff] for ff in range(0,len(freq),freqstep)]
    datfid.write('FREQUENCIES:'+' '*8+str(len(freq))+'\n')
    for ff in freq:
        datfid.write(ss+fmt % ff +'\n')
    
    #write data block
    datfid.write('DATA BLOCKS:'+' '*10+str(len(reslst))+'\n')
    datfid.write('SITE'+ss+'FREQ'+ss+'TYPE'+ss+'DATUM'+ss+'ERROR'+'\n')
    for ll,datline in enumerate(reslst):
        if datline.find('#IND')>=0:
            print 'Found #IND on line ',ll
            ndline=datline.replace('#IND','00')
            print 'Replaced with 00'
            datfid.write(ndline)
        else:
            datfid.write(datline)
    datfid.close()
    
    print 'Wrote Occam2D data file to: ',dfn
    
    return dfn
            
    
    
                                
def read2DRespFile(respfilename,datafilename):
    """
    read2DRespFile will read in a response file and combine the data with info 
    from the data file.

    Input:
        respfilename = full path to the response file
        datafilename = full path to data file

    Outputs:
        for each data array, the rows are ordered as:
            0 -> input data
            1 -> input error
            2 -> model output
            3 -> relative error (data-model)/(input error)
            
        rplst = list of dictionaries for each station with keywords:
            'station' = station name
            'offset' = relative offset,
            'resxy' = TE resistivity 
            'resyx'= TM resistivity 
            'phasexy'= TE phase 
            'phaseyx'= TM phase a
            'realtip'= Real Tipper 
            'imagtip'= Imaginary Tipper 
            
            Note: that the resistivity will be in log10 space.  Also, there are
            2 extra rows in the data arrays, this is to put the response from
            the inversion. 
        
        stationlst = list of stations in order from one side of the profile
                     to the other.
        freq = list of frequencies used in the inversion
        title = title, could be useful for plotting.
        
    """
    
    rplst,stationlst,freq,title,theta=read2DdataFile(datafilename)
    
    rfid=open(respfilename,'r')
    
    rlines=rfid.readlines()
    for line in rlines:
        ls=line.split()
        #station index
        ss=int(float(ls[0]))-1
        #component key
        comp=str(int(float(ls[2])))
        #frequency index        
        ff=int(float(ls[1]))-1
        #put into array
        #model response
        rplst[ss][occamdict[comp]][2,ff]=float(ls[5]) 
        #relative error        
        rplst[ss][occamdict[comp]][3,ff]=float(ls[6]) 
        
    return rplst,stationlst,np.array(freq),title

def read2DIterFile(iterfilename,datafilename,iterpath=None):
    """
    read2DIterFile will read an iteration file and combine that info from the 
    datafilename and return a dictionary of variables.
    
    Inputs:
        iterfilename = full path to iteration file if iterpath=None.  If 
                       iterpath is input then iterfilename is just the name
                       of the file without the full path.
        datafilename = full path to data file if iterpath=None.  If 
                       iterpath is input then datafilename is just the name
                       of the file without the full path.
    
    Outputs:
        
    """

    #get full paths if not already input
    if iterpath!=None and os.path.dirname(iterfilename)=='':
        ifn=os.path.join(iterpath,iterfilename)
    else:
        ifn=iterfilename
    if os.path.exists(ifn)==False:
        raise IOError('File: '+ifn+' does not exist, check name and path')
        
    if iterpath!=None and os.path.dirname(datafilename)=='':
        dfn=os.path.join(iterpath,datafilename)
    else:
        dfn=datafilename
    if os.path.exists(dfn)==False:
        raise IOError('File: '+dfn+' does not exist, check name and path')
        
    pass


def plot2DResponses(datafilename,respfilename=None,winglinkdatafile=None,
                    maxcol=8,plottype='1',ms=4,**kwargs):
    """
    plotResponse will plot the responses modeled from winglink against the 
    observed data.
    
    Inputs:
        respfilename = full path to response file
        datafilename = full path to data file
        winglinkdatfile = full path to a winglink data file used for a similar
                          inversion.  This will be plotted on the response
                          plots for comparison of fits.
        maxcol = maximum number of columns for the plot
        plottype = 'all' to plot all on the same plot
                   '1' to plot each respones in a different figure
                   station to plot a single station or enter as a list of 
                   stations to plot a few stations [station1,station2].  Does
                   not have to be verbatim but should have similar unique 
                   characters input pb01 for pb01cs in outputfile
    Outputs:
        None
    """
    
    if respfilename!=None:
        #read in the data    
        rplst,stationlst,freq,title=read2DRespFile(respfilename,datafilename)
        #make a legend list for plotting
        legendlst=['$Obs_{xy}$','$Obs_{yx}$','$Mod_{xy}$','$Mod_{yx}$']
        plotresp=True
    else:
        rplst,stationlst,freq,title,theta=read2DdataFile(datafilename)
        #make a legend list for plotting
        legendlst=['$Obs_{xy}$','$Obs_{yx}$']
        plotresp=False
    
    #boolean for adding winglink output to the plots 0 for no, 1 for yes
    addwl=0
    hspace=.15
    #read in winglink data file
    if winglinkdatafile!=None:
        addwl=1
        hspace=.25
        wld,wlrplst,plst,slst,tlst=wlt.readOutputFile(winglinkdatafile)
        legendlst=['$Obs_{xy}$','$Obs_{yx}$','Occ_$Mod_{xy}$','Occ_$Mod_{yx}$',
                   'Wl_$Mod_{xy}$','Wl_$Mod_{yx}$']
        sdict=dict([(ostation,wlstation) for wlstation in slst 
                    for ostation in stationlst if wlstation.find(ostation)>=0])
    period=1./freq
    nf=len(period)
    
    nstations=len(stationlst)
    
    #plot all responses onto one plot
    if plottype=='all':
        maxcol=8         
        nrows=int(np.ceil(nstations/float(maxcol)))
        
        fig=plt.figure(1,[14,10])
        gs=gridspec.GridSpec(nrows,1,hspace=hspace,left=.05,right=.98)
        count=0
        for rr in range(nrows):
            g1=gridspec.GridSpecFromSubplotSpec(6,maxcol,subplot_spec=gs[rr],
                                                    hspace=.15,wspace=.05)
            count=rr*(maxcol)
            for cc in range(maxcol):
                rlst=[]
                try:
                    ii=count+cc
                    stationlst[ii]
                except IndexError:
                    break
                rmslst=np.hstack((rplst[ii]['resxy'][3],
                                       rplst[ii]['resyx'][3],
                                        rplst[ii]['phasexy'][3],
                                        rplst[ii]['phaseyx'][3]))
                rms=np.sqrt(np.sum(ms**2 for ms in rmslst)/len(rmslst))
                #plot resistivity
                axr=plt.Subplot(fig,g1[:4,cc])
                fig.add_subplot(axr)
                #cut out missing data points first
                rxy=np.where(rplst[ii]['resxy'][0]!=0)[0]
                ryx=np.where(rplst[ii]['resyx'][0]!=0)[0]
                r1=axr.loglog(period[rxy],10**rplst[ii]['resxy'][0][rxy],
                              ls=':',marker='s',ms=ms,color='b',mfc='b')
                r2=axr.loglog(period[ryx],10**rplst[ii]['resyx'][0][ryx],
                              ls=':',marker='o',ms=ms,color='r',mfc='r')
                if plotresp==True:
                    mrxy=[np.where(rplst[ii]['resxy'][2]!=0)[0]]
                    mryx=[np.where(rplst[ii]['resyx'][2]!=0)[0]]
                    r3=axr.loglog(period[mrxy],10**rplst[ii]['resxy'][2][mrxy],
                                  ls='--',marker='+', ms=2*ms,color='g',mfc='g')
                    r4=axr.loglog(period[mryx],10**rplst[ii]['resyx'][2][mryx],
                                  ls='--',marker='+',ms=2*ms,color='m',mfc='m')
                
                    rlst=[r1,r2,r3,r4]
                else:
                    rlst=[r1,r2]
                #plot phase
                axp=plt.Subplot(fig,g1[-2:,cc])
                fig.add_subplot(axp)
                #cut out missing data points first
                pxy=[np.where(rplst[ii]['phasexy'][0]!=0)[0]]
                pyx=[np.where(rplst[ii]['phaseyx'][0]!=0)[0]]
                axp.semilogx(period[pxy],rplst[ii]['phasexy'][0][pxy],
                             ls=':',marker='s',ms=ms,color='b',mfc='b')
                axp.semilogx(period[pyx],rplst[ii]['phaseyx'][0][pyx],
                             ls=':',marker='o',ms=ms,color='r',mfc='r')
                if plotresp==True:
                    mpxy=[np.where(rplst[ii]['phasexy'][2]!=0)[0]]
                    mpyx=[np.where(rplst[ii]['phaseyx'][2]!=0)[0]]
                    axp.semilogx(period[mpxy],rplst[ii]['phasexy'][2][mpxy],
                                 ls='--',marker='+',ms=2*ms,color='g',mfc='g')
                    axp.semilogx(period[mpyx],rplst[ii]['phaseyx'][2][mpyx],
                                 ls='--',marker='+',ms=2*ms,color='m',mfc='m')
                
                #add in winglink responses
                if addwl==1:
                    try:
                        wlrms=wld[sdict[stationlst[ii]]]['rms']
                        axr.set_title(stationlst[ii]+'\n'+\
                                    'rms_occ, rms_wl= %.2f, %.2f' % (rms,wlrms),
                                     fontdict={'size':12,'weight':'bold'})
                        r5=axr.loglog(wld[sdict[stationlst[ii]]]['period'],
                                   wld[sdict[stationlst[ii]]]['modresxy'],
                                   ls='-.',marker='x',ms=2*ms,color='cyan',
                                   mfc='cyan')
                        r6=axr.loglog(wld[sdict[stationlst[ii]]]['period'],
                                   wld[sdict[stationlst[ii]]]['modresyx'],
                                   ls='-.',marker='x',ms=2*ms,color='orange',
                                   mfc='orange')
                        axp.semilogx(wld[sdict[stationlst[ii]]]['period'],
                                     wld[sdict[stationlst[ii]]]['modphasexy'],
                                     ls='-.',marker='x',ms=2*ms,color='cyan',
                                     mfc='cyan')
                        axp.semilogx(wld[sdict[stationlst[ii]]]['period'],
                                     wld[sdict[stationlst[ii]]]['modphaseyx'],
                                     ls='-.',marker='x',ms=2*ms,color='orange',
                                     mfc='orange')
                        rlst.append(r5)
                        rlst.append(r6)
                    except IndexError:
                        print 'Station not present'
                else:
                    if plotresp==True:
                        axr.set_title(stationlst[ii]+'; rms= %.2f' % rms,
                                      fontdict={'size':12,'weight':'bold'})
                    else:
                        axr.set_title(stationlst[ii],
                                      fontdict={'size':12,'weight':'bold'})
                
                #make plot nice with labels
                if cc==0 and rr==0:
                    fig.legend(rlst,legendlst,
                                loc='upper center',markerscale=2,
                                borderaxespad=.35,
                                labelspacing=.08,
                                handletextpad=.15,borderpad=.1,
                                ncol=len(rlst))
                axr.grid(True)
                axr.set_xticklabels(['' for ii in range(10)])
                if cc>0:
                    axr.set_yticklabels(['' for ii in range(6)])
                    
                axp.set_ylim(-10,120)
                axp.grid(True)
                axp.yaxis.set_major_locator(MultipleLocator(30))
                axp.yaxis.set_minor_locator(MultipleLocator(5))
                
                if cc==0:
                    axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                   fontdict={'size':12,'weight':'bold'})
                    axp.set_ylabel('Phase (deg)',
                                   fontdict={'size':12,'weight':'bold'})
                    axr.yaxis.set_label_coords(-.15,.5)
                    axp.yaxis.set_label_coords(-.15,.5)
        
                if cc>0:
                    axr.set_yticklabels(['' for ii in range(6)])
                    axp.set_yticklabels(['' for ii in range(6)])
                if rr==nrows-1:
                    axp.set_xlabel('Period (s)',
                                   fontdict={'size':12,'weight':'bold'})
                                   
    #plot each respones in a different figure
    elif plottype=='1':
        gs=gridspec.GridSpec(6,2,wspace=.05)
        for ii,station in enumerate(stationlst):
            
            rlst=[]
            llst=[]
            
            rmslst=np.hstack((rplst[ii]['resxy'][3],
                                       rplst[ii]['resyx'][3],
                                        rplst[ii]['phasexy'][3],
                                        rplst[ii]['phaseyx'][3]))
            rms=np.sqrt(np.sum(ms**2 for ms in rmslst)/len(rmslst))
            fig=plt.figure(ii+1,[7,8])
            
            #plot resistivity
            axr=fig.add_subplot(gs[:4,:])
            #cut out missing data points first
            rxy=np.where(rplst[ii]['resxy'][0]!=0)[0]
            ryx=np.where(rplst[ii]['resyx'][0]!=0)[0]
            
            #check to see if there is a yx component
            if len(rxy)>0:
                r1=axr.errorbar(period[rxy],10**rplst[ii]['resxy'][0][rxy],
                                ls=':',marker='s',ms=ms,mfc='b',mec='b',
                                color='b',yerr=rplst[ii]['resxy'][1][rxy],
                                ecolor='b')
                rlst.append(r1[0])
                llst.append('$Obs_{xy}$')
            else:
                pass
            
            #check to see if there is a yx component
            if len(ryx)>0:
                r2=axr.errorbar(period[ryx],10**rplst[ii]['resyx'][0][ryx],
                                ls=':',marker='o',ms=ms,mfc='r',mec='r',
                                color='r',yerr=rplst[ii]['resyx'][1][ryx],
                                ecolor='r')
                rlst.append(r2[0])
                llst.append('$Obs_{yx}$')
            else:
                pass                                
            
            if plotresp==True:
                mrxy=np.where(rplst[ii]['resxy'][2]!=0)[0]
                mryx=np.where(rplst[ii]['resyx'][2]!=0)[0]
                
                #check for the xy of model component
                if len(mrxy)>0:
                    r3=axr.errorbar(period[mrxy],10**rplst[ii]['resxy'][2][mrxy],
                                    ls='--',marker='+',ms=ms,mfc='g',mec='g',
                                    color='g',yerr=rplst[ii]['resxy'][3][mrxy],
                                    ecolor='g')
                    rlst.append(r3[0])
                    llst.append('$Mod_{xy}$')
                else:
                    pass
                
                #check for the yx model component  of resisitivity
                if len(mryx)>0:
                    r4=axr.errorbar(period[mryx],10**rplst[ii]['resyx'][2][mryx],
                                    ls='--',marker='+',ms=ms,mfc='m',mec='m',
                                    color='m',yerr=rplst[ii]['resyx'][3][mryx],
                                    ecolor='m')
                    rlst.append(r4[0])
                    llst.append('$Mod_{yx}$')
                                
            #plot phase
            axp=fig.add_subplot(gs[-2:,:],sharex=axr)
            
            #cut out missing data points first
            pxy=np.where(rplst[ii]['phasexy'][0]!=0)[0]
            pyx=np.where(rplst[ii]['phaseyx'][0]!=0)[0]

            if len(pxy)>0:
                axp.errorbar(period[pxy],rplst[ii]['phasexy'][0][pxy],
                             ls=':',marker='s',ms=ms,mfc='b',mec='b',color='b',
                             yerr=rplst[ii]['phasexy'][1][pxy],ecolor='b')
            else:
                pass
            
            if len(pyx)>0:
                axp.errorbar(period[pyx],rplst[ii]['phaseyx'][0][pyx],
                             ls=':',marker='o',ms=ms,mfc='r',mec='r',color='r',
                             yerr=rplst[ii]['phaseyx'][1][pyx],ecolor='r')
            else:
                pass
            
            if plotresp==True:
                mpxy=np.where(rplst[ii]['phasexy'][2]!=0)[0]
                mpyx=np.where(rplst[ii]['phaseyx'][2]!=0)[0]
                
                if len(mpxy)>0:
                    axp.errorbar(period[mpxy],rplst[ii]['phasexy'][2][mpxy],
                                 ls='--',marker='+',ms=ms,mfc='g',mec='g',
                                 color='g',yerr=rplst[ii]['phasexy'][3][mpxy],
                                 ecolor='g')
                else:
                    pass
                
                if len(mpyx)>0:
                    axp.errorbar(period[mpyx],rplst[ii]['phaseyx'][2][mpyx],
                                 ls='--',marker='+',ms=ms,mfc='m',mec='m',
                                 color='m',yerr=rplst[ii]['phaseyx'][3][mpyx],
                                 ecolor='m')
                else:
                    pass
#                axp.semilogx(period[mpxy],rplst[ii]['phasexy'][2][mpxy],
#                             ls='--',marker='+',ms=2*ms,color='g',mfc='g')
#                axp.semilogx(period[mpyx],rplst[ii]['phaseyx'][2][mpyx],
#                             ls='--',marker='+',ms=2*ms,color='m',mfc='m')
                         
            #add in winglink responses
            if addwl==1:
                try:
                    wlrms=wld[sdict[station]]['rms']
                    axr.set_title(stationlst[ii]+'\n'+\
                                'rms_occ, rms_wl= %.2f, %.2f' % (rms,wlrms),
                                 fontdict={'size':12,'weight':'bold'})
                    
                                 
                    r5=axr.loglog(wld[sdict[station]]['period'],
                               wld[sdict[station]]['modresxy'],
                               ls='-.',marker='x',ms=2*ms,color='cyan',
                               mfc='cyan')
                    r6=axr.loglog(wld[sdict[station]]['period'],
                               wld[sdict[stationlst[ii]]]['modresyx'],
                               ls='-.',marker='x',ms=2*ms,color='orange',
                               mfc='orange')
                    axp.semilogx(wld[sdict[station]]['period'],
                                 wld[sdict[station]]['modphasexy'],
                                 ls='-.',marker='x',ms=2*ms,color='cyan',
                                 mfc='cyan')
                    axp.semilogx(wld[sdict[station]]['period'],
                                 wld[sdict[station]]['modphaseyx'],
                                 ls='-.',marker='x',ms=2*ms,color='orange',
                                 mfc='orange')
                except IndexError:
                    print 'Station not present'
            else:
                axr.set_title(stationlst[ii]+'; rms= %.2f' % rms,
                              fontdict={'size':12,'weight':'bold'})
                              
            axr.set_xscale('log')
            axp.set_xscale('log')
            axr.set_yscale('log')
            axr.grid(True)
            axr.set_xticklabels(['' for ii in range(10)])
            axp.set_ylim(-10,120)
            axp.grid(True)
            axp.yaxis.set_major_locator(MultipleLocator(10))
            axp.yaxis.set_minor_locator(MultipleLocator(1))
            
            axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_ylabel('Phase (deg)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})
            axr.legend(rlst,llst,
                       loc=2,markerscale=1,borderaxespad=.05,
                       labelspacing=.08,
                       handletextpad=.15,borderpad=.05)
            axr.yaxis.set_label_coords(-.05,.5)
            axp.yaxis.set_label_coords(-.05,.5)
    
    else:
        pstationlst=[]
        if type(plottype) is not list:
            plottype=[plottype]
        for ii,station in enumerate(stationlst):
            for pstation in plottype:
                if station.find(pstation)>=0:
                    print 'plotting ',station
                    pstationlst.append(ii)

        gs=gridspec.GridSpec(6,2,wspace=.05,left=.1,top=.93,bottom=.07)
        for ii in pstationlst:
            rlst=[]
            pstation=stationlst[ii]
            rmslst=np.hstack((rplst[ii]['resxy'][3],
                                       rplst[ii]['resyx'][3],
                                        rplst[ii]['phasexy'][3],
                                        rplst[ii]['phaseyx'][3]))
            rms=np.sqrt(np.sum(ms**2 for ms in rmslst)/len(rmslst))
            fig=plt.figure(ii+1,[7,7])
            
            #plot resistivity
            #cut out missing data points first
            axr=fig.add_subplot(gs[:4,:])
            rxy=np.where(rplst[ii]['resxy'][0]!=0)[0]
            ryx=np.where(rplst[ii]['resyx'][0]!=0)[0]
            r1=axr.errorbar(period[rxy],10**rplst[ii]['resxy'][0][rxy],
                            ls=':',marker='s',ms=ms,mfc='b',mec='b',color='b',
                            yerr=rplst[ii]['resxy'][1][rxy],ecolor='b')
            r2=axr.errorbar(period[ryx],10**rplst[ii]['resyx'][0][ryx],
                            ls=':',marker='o',ms=ms,mfc='r',mec='r',color='r',
                            yerr=rplst[ii]['resyx'][1][ryx],ecolor='r')
#            r1=axr.loglog(period[rxy],10**rplst[ii]['resxy'][0][rxy],
#                          ls=':',marker='s',ms=ms,color='b',mfc='b')
#            r2=axr.loglog(period[ryx],10**rplst[ii]['resyx'][0][ryx],
#                          ls=':',marker='o',ms=ms,color='r',mfc='r')
            if plotresp==True:
                mrxy=[np.where(rplst[ii]['resxy'][2]!=0)[0]]
                mryx=[np.where(rplst[ii]['resyx'][2]!=0)[0]]
                r3=axr.errorbar(period[mrxy],10**rplst[ii]['resxy'][2][mrxy],
                                ls='--',marker='+',ms=ms,mfc='g',mec='g',color='g',
                                yerr=rplst[ii]['resxy'][3][mrxy],ecolor='g')
                r4=axr.errorbar(period[mryx],10**rplst[ii]['resyx'][2][mryx],
                                ls='--',marker='+',ms=ms,mfc='m',mec='m',color='m',
                                yerr=rplst[ii]['resyx'][3][mryx],ecolor='m')
#                r3=axr.loglog(period[mrxy],10**rplst[ii]['resxy'][2][mrxy],
#                              ls='--',marker='+', ms=2*ms,color='g',mfc='g')
#                r4=axr.loglog(period[mryx],10**rplst[ii]['resyx'][2][mryx],
#                              ls='--',marker='+',ms=2*ms,color='m',mfc='m')
            
                rlst=[r1[0],r2[0],r3[0],r4[0]]
            else:
                rlst=[r1[0],r2[0]]
                                
            #plot phase
            axp=fig.add_subplot(gs[-2:,:],sharex=axr)
            #cut out missing data points first
            pxy=[np.where(rplst[ii]['phasexy'][0]!=0)[0]]
            pyx=[np.where(rplst[ii]['phaseyx'][0]!=0)[0]]
            axp.errorbar(period[pxy],rplst[ii]['phasexy'][0][pxy],
                            ls=':',marker='s',ms=ms,mfc='b',mec='b',color='b',
                            yerr=rplst[ii]['phasexy'][1][pxy],ecolor='b')
            axp.errorbar(period[pyx],rplst[ii]['phaseyx'][0][pyx],
                            ls=':',marker='o',ms=ms,mfc='r',mec='r',color='r',
                            yerr=rplst[ii]['phaseyx'][1][pyx],ecolor='r')
#            axp.semilogx(period[pxy],rplst[ii]['phasexy'][0][pxy],
#                         ls=':',marker='s',ms=ms,color='b',mfc='b')
#            axp.semilogx(period[pyx],rplst[ii]['phaseyx'][0][pyx],
#                         ls=':',marker='o',ms=ms,color='r',mfc='r')
            if plotresp==True:
                mpxy=[np.where(rplst[ii]['phasexy'][2]!=0)[0]]
                mpyx=[np.where(rplst[ii]['phaseyx'][2]!=0)[0]]
                
                axp.errorbar(period[mpxy],rplst[ii]['phasexy'][2][mpxy],
                            ls='--',marker='+',ms=ms,mfc='g',mec='g',color='g',
                            yerr=rplst[ii]['phasexy'][3][mpxy],ecolor='g')
                axp.errorbar(period[mpyx],rplst[ii]['phaseyx'][2][mpyx],
                            ls='--',marker='+',ms=ms,mfc='m',mec='m',color='m',
                            yerr=rplst[ii]['phaseyx'][3][mpyx],ecolor='m')
#                axp.semilogx(period[mpxy],rplst[ii]['phasexy'][2][mpxy],
#                             ls='--',marker='+',ms=2*ms,color='g',mfc='g')
#                axp.semilogx(period[mpyx],rplst[ii]['phaseyx'][2][mpyx],
#                             ls='--',marker='+',ms=2*ms,color='m',mfc='m')
                         
            #add in winglink responses
            if addwl==1:
                try:
                    wlrms=wld[sdict[station]]['rms']
                    print 'plotting WL station: ', sdict[pstation]
                    axr.set_title(pstation+'\n'+\
                                'rms_occ, rms_wl= %.2f, %.2f' % (rms,wlrms),
                                 fontdict={'size':12,'weight':'bold'})
                    
                    r5=axr.loglog(wld[sdict[pstation]]['period'],
                               wld[sdict[pstation]]['modresxy'],
                               ls='-.',marker='x',ms=2*ms,color='cyan',
                               mfc='cyan')
                    r6=axr.loglog(wld[sdict[pstation]]['period'],
                               wld[sdict[pstation]]['modresyx'],
                               ls='-.',marker='x',ms=2*ms,color='orange',
                               mfc='orange')
                    axp.semilogx(wld[sdict[pstation]]['period'],
                                 wld[sdict[pstation]]['modphasexy'],
                                 ls='-.',marker='x',ms=2*ms,color='cyan',
                                 mfc='cyan')
                    axp.semilogx(wld[sdict[pstation]]['period'],
                                 wld[sdict[pstation]]['modphaseyx'],
                                 ls='-.',marker='x',ms=2*ms,color='orange',
                                 mfc='orange')
                    rlst.append(r5)
                    rlst.append(r6)
                except IndexError:
                    print 'Station not present'
            else:
                if plotresp==True:
                    axr.set_title(pstation+'; rms= %.2f' % rms,
                                  fontdict={'size':12,'weight':'bold'})
                else:
                    axr.set_title(pstation,
                                  fontdict={'size':12,'weight':'bold'})
            
            axr.set_xscale('log')
            axr.set_yscale('log')
            axp.set_xscale('log')                  
            axr.grid(True)
            axr.set_xticklabels(['' for ii in range(10)])
            axp.set_ylim(-10,120)
            axp.grid(True)
            axp.yaxis.set_major_locator(MultipleLocator(10))
            axp.yaxis.set_minor_locator(MultipleLocator(1))
            
            axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_ylabel('Phase (deg)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})
            axr.legend(rlst,legendlst,
                       loc=2,markerscale=2,borderaxespad=.05,
                       labelspacing=.08,
                       handletextpad=.15,borderpad=.05)
            axr.yaxis.set_label_coords(-.05,.5)
            axp.yaxis.set_label_coords(-.05,.5)
            
def plotTipper(datafile,):
    pass
    
def plotAllResponses(datafile,station,fignum=1):
    """
    Plot all the responses of occam inversion from data file.  This assumes
    the response curves are in the same folder as the datafile.

    Input:
        datafile = full path to occam data file
        
    Output:
        Plot
    
    """    
    
    rpath=os.path.dirname(datafile)
    
    gs=gridspec.GridSpec(6,2,wspace=.15)


    rlst=[os.path.join(rpath,rfile) for rfile in os.listdir(rpath) 
            if rfile.find('.resp')>0]
    
    nresp=len(rlst)
    
    colorlst=[(cc,0,1-cc) for cc in np.arange(0,1,1./nresp)]
    fig=plt.figure(fignum,[7,8])
    axrte=fig.add_subplot(gs[:4,0])
    axrtm=fig.add_subplot(gs[:4,1])
    axpte=fig.add_subplot(gs[-2:,0])
    axptm=fig.add_subplot(gs[-2:,1])
    rmstelst=[]
    rmstmlst=[]
    #read responses
    for jj,rfile in enumerate(rlst):
        rplst,stationlst,freq,title=read2DRespFile(os.path.join(rpath,rfile),
                                                 datafile)
        ii=np.where(np.array(stationlst)==station)[0][0]
        
        period=1./freq
        
        rmslstte=np.hstack((rplst[ii]['resxy'][3],
                            rplst[ii]['phasexy'][3]))
        rmslsttm=np.hstack((rplst[ii]['resyx'][3],
                            rplst[ii]['phaseyx'][3]))
        rmste=np.sqrt(np.sum(ms**2 for ms in rmslstte)/len(rmslstte))
        rmstm=np.sqrt(np.sum(ms**2 for ms in rmslsttm)/len(rmslsttm))
        rmstelst.append('%d rms=%.3f ' % (jj,rmste))
        rmstmlst.append('%d rms=%.3f ' % (jj,rmstm))
        #plot resistivity
        
        
        if jj==0:
            #cut out missing data points first
            rxy=np.where(rplst[ii]['resxy'][0]!=0)[0]
            ryx=np.where(rplst[ii]['resyx'][0]!=0)[0]
            r1=axrte.loglog(period[rxy],10**rplst[ii]['resxy'][0][rxy],
                          ls=':',marker='s',ms=4,color='k',mfc='k')
            r2=axrtm.loglog(period[ryx],10**rplst[ii]['resyx'][0][ryx],
                          ls=':',marker='o',ms=4,color='k',mfc='k')
            rlstte=[r1]
            rlsttm=[r2]
    
        mrxy=[np.where(rplst[ii]['resxy'][2]!=0)[0]]
        mryx=[np.where(rplst[ii]['resyx'][2]!=0)[0]]
        r3=axrte.loglog(period[mrxy],10**rplst[ii]['resxy'][2][mrxy],
                        ls='-',color=colorlst[jj])
        r4=axrtm.loglog(period[mryx],10**rplst[ii]['resyx'][2][mryx],
                        ls='-',color=colorlst[jj])
    
        rlstte.append(r3)
        rlsttm.append(r4)
                            
        #plot phase
        #cut out missing data points first
        pxy=[np.where(rplst[ii]['phasexy'][0]!=0)[0]]
        pyx=[np.where(rplst[ii]['phaseyx'][0]!=0)[0]]
        
        if jj==0:            
            axpte.semilogx(period[pxy],rplst[ii]['phasexy'][0][pxy],
                         ls=':',marker='s',ms=4,color='k',mfc='k')
            axptm.semilogx(period[pyx],rplst[ii]['phaseyx'][0][pyx],
                         ls=':',marker='o',ms=4,color='k',mfc='k')
                         
        mpxy=[np.where(rplst[ii]['phasexy'][2]!=0)[0]]
        mpyx=[np.where(rplst[ii]['phaseyx'][2]!=0)[0]]
        axpte.semilogx(period[mpxy],rplst[ii]['phasexy'][2][mpxy],
                     ls='-',color=colorlst[jj])
        axptm.semilogx(period[mpyx],rplst[ii]['phaseyx'][2][mpyx],
                     ls='-',color=colorlst[jj])
                     
    axrte.legend(rlstte,rmstelst,loc=2,markerscale=2,borderaxespad=.05,
               labelspacing=.08,
               handletextpad=.15,borderpad=.05)
    
    axrtm.legend(rlsttm,rmstmlst,loc=2,markerscale=2,borderaxespad=.05,
               labelspacing=.08,
               handletextpad=.15,borderpad=.05)
                   
    axrte.grid(True)
    axrtm.grid(True)
    
    
    axrtm.set_xticklabels(['' for ii in range(10)])
    axrte.set_xticklabels(['' for ii in range(10)])
    #axpte.set_ylim(-10,120)
    
    axrte.set_title('TE')
    axrtm.set_title('TM')
    
    axpte.grid(True)
    axpte.yaxis.set_major_locator(MultipleLocator(10))
    axpte.yaxis.set_minor_locator(MultipleLocator(1))
    
    axrte.set_ylabel('App. Res. ($\Omega \cdot m$)',
                   fontdict={'size':12,'weight':'bold'})
    axpte.set_ylabel('Phase (deg)',
                   fontdict={'size':12,'weight':'bold'})
    axpte.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})

    axrte.yaxis.set_label_coords(-.05,.5)
    axpte.yaxis.set_label_coords(-.05,.5)
    
    axrtm.set_xticklabels(['' for ii in range(10)])
    axptm.set_ylim(-10,120)
    axptm.grid(True)
    axptm.yaxis.set_major_locator(MultipleLocator(10))
    axptm.yaxis.set_minor_locator(MultipleLocator(1))
    
    axrtm.set_ylabel('App. Res. ($\Omega \cdot m$)',
                   fontdict={'size':12,'weight':'bold'})
    axptm.set_ylabel('Phase (deg)',
                   fontdict={'size':12,'weight':'bold'})
    axptm.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})

    axrtm.yaxis.set_label_coords(-.05,.5)
    axptm.yaxis.set_label_coords(-.05,.5)
    plt.suptitle(station,fontsize=14,fontweight='bold')
    plt.show()
    
    
def plotModel(iterfile,meshfile):
    """
    plot the mode as an image
    """