#!/usr/bin/env python

"""
This modules contains helper functions for the calibration of raw time series. 

The various functions deal with the calibration of data from fluxgates, coils, dipoles,...
The calibration depends on  the instrument as well as on the respective data logger. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import calendar
import time


from mtpy.utils.exceptions import *
#=================================================================

list_of_loggers = ['edl','elogger']
list_of_instruments = ['dipole','fluxgate', 'coils', 'bartington']
dict_of_EDL_gain_factors = {'high': 0.1, 'low': 1. , 'verylow': 2.5}


#=================================================================
# section for amplification factors

# list of gain levels

# z channel scaling factor



#=================================================================


def calibrate(raw_data, field, dipole_length=1., instrument, amplification, logger, gain):

    """
    Convert a given time series from raw data (instrument counts) 
    into field strength amplitude values. 
    The dipole length factor must be set to 1 for the calibration of B-field data.

    The output is a time series of field strngth values in basic units: 
    Tesla for the B-field and V/m for the E-field.

    input:
    - 1D time series (list or numpy array)
    - type of measured field (E/B)
    - instrument (fluxgate, coils, dipole, bartington)
    - amplification factor for the instrument
    - logger (EDL, elogger, other)
    - gain factor of the logger

    output:
    - 1D time series (numpy array)

    """

    raw_data = np.array(raw_data)

    _data_instrument_consitency_check(raw_data,field, dipole_length, instrument, amplification, logger, gain )
    

    data = raw_data / dipole_length * amplification * gain


    return data

#=================================================================




def _data_instrument_consitency_check(data,field, dipole, instrument, amplification, logger, gain):
    """
    Check, if the input values given  for the calibration make any sense at all.

    """

    if len(data) == 0 :
        raise MTpyError_ts_data( 'no data provided for calibration' )

    if not field.lower() in ['e','b']:
        raise MTpyError_inputarguments( 'Field must be E or B' )

    if float(dipole) <= 0:
        raise MTpyError_inputarguments( 'Dipole length must be positive' )

    if float(amplification) <= 0:
        raise MTpyError_inputarguments( 'Amplification factor must be positive' )

    if float(gain) <= 0:
        raise MTpyError_inputarguments( 'Instrument gain must be positive' )

    try:
        if not logger.lower() in list_of_loggers:
            raise
    except:
        raise MTpyError_inputarguments( 'wrong choice of logger')

    try:
        if not instrument.lower() in list_of_instruments:
            raise
    except:
        raise MTpyError_inputarguments( 'wrong choice of instrument')

    if field.lower == 'b':
        if logger.lower() == 'elogger':
            raise MTpyError_inputarguments( 'wrong choice of logger')
        if instrument.lower() == 'dipole'
            raise MTpyError_inputarguments( 'wrong choice of instrument')
        if not float(dipole) == 1:
            raise MTpyError_inputarguments( 'Dipole length must be 1 for B-field calibration')
           

    if field.lower == 'e':
        if not instrument.lower() == 'dipole'
            raise MTpyError_inputarguments( 'wrong choice of instrument')



def e_field(raw_data):



    return data



def b_field(raw_data):


    return data

#-=---------------------------------------------------------
# old functions





def convertfiles(dirpath,folder,infodict,fmt='%.6g'):

    """
    convertfiles will convert data of counts from data logger to units.
    """
    aconvstr=' has already been converted check data file'+'\n'
    delemptyfile=' has been deleted because the file was empty'+'\n'
    clines=[]
    clines.append('======'+folder+'======'+'\n')
    for dayfolder in os.listdir(os.path.join(dirpath,folder)):
        if dayfolder.find('.')==-1:
            clines.append('---'+dayfolder+'---'+'\n')
            for filename in os.listdir(os.path.join(dirpath,folder,dayfolder)):
                if filename.find('.')>=0:
                    if fnmatch.fnmatch(filename,'*.EX'):
                        exfid=file(os.path.join(dirpath,folder,dayfolder,
                                                filename),'r')
                        exlines=exfid.readlines()
                        if len(exlines)==0:
                            #os.remove(os.path.join(dirpath,folder,dayfolder,filename))
                            clines.append(filename+delemptyfile)                        
                        
                        elif exlines[0].find('.')>=0:
                            exfid.close()
                            clines.append(filename+aconvstr)
                        
                        else:
                            exconv=convertE(exlines,infodict['dlgain'],
                                                    infodict['egain'],
                                                    infodict['ex'])
                            exfid.close()
                            exconvlst=[fmt % exconv[ii] +'\n' for ii in range(len(exconv))]
                            exfidn=file(os.path.join(dirpath,folder,dayfolder,
                                                     filename),'w')
                            exfidn.writelines(exconvlst)
                            exfidn.close()
                            clines.append(filename+'\n')
                    elif fnmatch.fnmatch(filename,'*.EY'):
                        eyfid=file(os.path.join(dirpath,folder,dayfolder,
                                                filename),'r')
                        eylines=eyfid.readlines()
                        if len(eylines)==0:
#                            os.remove(os.path.join(dirpath,folder,dayfolder,
#                                                   filename))
                            clines.append(filename+delemptyfile)
                        elif eylines[0].find('.')>=0:
                            eyfid.close()
                            #clines.append(filename+aconvstr)
                        
                        else:
                            eyconv=convertE(eylines,infodict['dlgain'],
                                                    infodict['egain'],
                                                    infodict['ey'])
                            eyfid.close()
                            eyconvlst=[fmt % eyconv[ii] +'\n' for ii in range(len(eyconv))]
                            eyfidn=file(os.path.join(dirpath,folder,dayfolder,
                                                     filename),'w')
                            eyfidn.writelines(eyconvlst)
                            eyfidn.close()
                            clines.append(filename+'\n')
                else:
                    clines.append('Found Folder: '+filename+'\n')
            if infodict['magtype']=='lp':
                magoristr=infodict['magori'].replace('"','')
                magorilst=magoristr.split(',')
                for filename in os.listdir(os.path.join(dirpath,folder,dayfolder)):
                    if filename.find('.')>=0:
                        if fnmatch.fnmatch(filename,'*.'+magorilst[0]):
                            bxfid=file(os.path.join(dirpath,folder,dayfolder,
                                                    filename),'r')
                            bxlines=bxfid.readlines()
                            if len(bxlines)==0:
#                                os.remove(os.path.join(dirpath,folder,dayfolder,
#                                                       filename))
                                clines.append(filename+delemptyfile)
                            elif bxlines[0].find('.')>=0:
                                bxfid.close()
                                clines.append(filename+aconvstr)                                
                            else:
                                bxconv=convertlpB(bxlines,infodict['dlgain'])
                                bxfid.close()
                                bxconvlst=[fmt % bxconv[ii] + '\n' for ii in range(len(bxconv))]
                                bxfidn=file(os.path.join(dirpath,folder,
                                                         dayfolder,filename),'w')
                                bxfidn.writelines(bxconvlst)
                                bxfidn.close()
                                clines.append(filename+' as BX'+'\n')
                        elif fnmatch.fnmatch(filename,'*.'+magorilst[1]):
                            byfid=file(os.path.join(dirpath,folder,dayfolder,
                                                    filename),'r')
                            bylines=byfid.readlines()
                            
                            if len(bylines)==0:
#                                os.remove(os.path.join(dirpath,folder,dayfolder,
#                                                       filename))
                                clines.append(filename+delemptyfile)
                            elif bylines[0].find('.')>=0:
                                byfid.close()
                                clines.append(filename+aconvstr)
                            
                            else:
                                byconv=convertlpB(bylines,
                                                          infodict['dlgain'])
                                byfid.close()
                                byconvlst=[fmt % byconv[ii] +'\n' for ii in range(len(byconv))]
                                byfidn=file(os.path.join(dirpath,folder,
                                                         dayfolder,filename),'w')
                                byfidn.writelines(byconvlst)
                                byfidn.close()
                                clines.append(filename+' as BY'+'\n')
                        elif fnmatch.fnmatch(filename,'*.'+magorilst[2]):
                            bzfid=file(os.path.join(dirpath,folder,dayfolder,
                                                    filename),'r')
                            bzlines=bzfid.readlines()
                            
                            if len(bzlines)==0:
#                                os.remove(os.path.join(dirpath,folder,dayfolder,
#                                                       filename))
                                clines.append(filename+delemptyfile)
                            elif bzlines[0].find('.')>=0:
                                bzfid.close()
                                clines.append(filename+aconvstr)
                            
                            else:
                                bzconv=convertlpB(bzlines,
                                                          infodict['dlgain'],
                                                          zadj=infodict['lpbzcor'])
                                bzfid.close()
                                bzconvlst=[fmt % bzconv[ii]+'\n' for ii in range(len(bzconv))]
                                bzfidn=file(os.path.join(dirpath,folder,
                                                         dayfolder,filename),'w')
                                bzfidn.writelines(bzconvlst)
                                bzfidn.close()
                                clines.append(filename+' as BZ'+'\n')
                        else:
                            pass
                    else:
                        clines.append('Found Folder: '+filename+'\n')
    return clines



def convertlpB(bfield,dlgain=1,zadj=1):
    """
    Convert the magnetic field from counts to units of microV/nT.
    bfield is a list of numbers. dlain is amount of gain applied
    by data logger(verylow=2.5,low=1, high=.1)
    
    Inputs:
        bfield = 1D array of long period data
        dlgain = data logger gain (very low= 2.5,low = 1, high = .1)
        zadj = Bz adjustment if using corrected Bartingtons
    
    Outputs:
        bfieldc = scaled bfield 1D array
    """
    bfieldc=np.array(bfield,dtype='float')/10.E7*70000.*float(dlgain)*float(zadj)

    return bfieldc

def convertE(efield,dlgain,egain,dlength):
    """
    Convert the electric field from counts to units of microV/m.
    efield is a list, dlgain is gain applied by data logger(verylow=2.5,
    low=1, high=.1), egain is interface box gain (10,100),dlength is
    the dipole length of that component.

    Inputs:
        efield = 1d array of electric field measurements
        dlgain = data logger gain (very low= 2.5,low = 1, high = .1)
        egain = gain from interface box
        dlength = length of dipole in (m)
    
    Outputs:
        efieldc = scaled electric field 1D array
    """
    efieldc=np.array(efield,dtype='float')*float(dlgain)/(float(dlength)*\
             float(egain))

    return efieldc



def convertCounts2Units(filenames,eyn='n',lpyn='n',egain=1.0,dlgain=1.0,exlen=100.,eylen=100.,magtype='lp',zadj=2):
    """convertCounts2Units(filenames,eyn='n',lpyn='n',egain=1.0,dlgain=1.0,exlen=100.,eylen=100.,
    magtype='lp',zadj=2) will convert a set of files given by filenames with
    parameters defined:

    filenames,
    eyn => electric channels converted y or n
    lpyn => long period magnetic channels converted y or n
    egain => electric gain
    dlgain => data logger gain
    exlen => ex length (m)
    eylen => ey length (m)
    magtype => magnetic sensor type lp for longperiod or bb for broadband
    zadj => bz adjusting parameter for bartington sensor"""
    
    for ii in range(len(filenames)):
        if eyn=='n':
            #convert EX chanel
            if fnmatch.fnmatch(filenames[ii],'*.EX'):
                exfid=file(filenames[ii],'r')
                exlines=exfid.readlines()
                if exlines[0].find('.')>=0:
                        print 'Found decimal point in '+filenames[ii]+'. Check File'
                        exyn=input('Still convert? (y/n) as a string')
                        if exyn=='n':
                            exfid.close()
                        else:
                            exconv=convertE(exlines,dlgain,egain,exlen)
                            exfid.close()
                            exconvlst=[str(exconv[ii])+'\n' for ii in range(len(exconv))]
                            exfidn=file(filenames[ii],'w')
                            exfidn.writelines(exconvlst)
                            exfidn.close()
                else:
                    exconv=convertE(exlines,dlgain,egain,exlen)
                    exfid.close()
                    exconvlst=[str(exconv[ii])+'\n' for ii in range(len(exconv))]
                    exfidn=file(filenames[ii],'w')
                    exfidn.writelines(exconvlst)
                    exfidn.close()
            elif fnmatch.fnmatch(filenames[ii],'*.EY'):
                eyfid=file(filenames[ii],'r')
                eylines=eyfid.readlines()
                if eylines[0].find('.')>=0:
                        print 'Found decimal point '+filenames[ii]+'. Check File.'
                        eyyn=input('Still convert? (y/n) as a string')
                        if eyyn=='n':
                            eyfid.close()
                        else:
                            eyconv=convertE(eylines,dlgain,egain,eylen)
                            eyfid.close()
                            eyconvlst=[str(eyconv[ii])+'\n' for ii in range(len(eyconv))]
                            eyfidn=file(filenames[ii],'w')
                            eyfidn.writelines(eyconvlst)
                            eyfidn.close()
                else:
                    eyconv=convertE(eylines,dlgain,egain,eylen)
                    eyfid.close()
                    eyconvlst=[str(eyconv[ii])+'\n' for ii in range(len(eyconv))]
                    eyfidn=file(filenames[ii],'w')
                    eyfidn.writelines(eyconvlst)
                    eyfidn.close()
        else:
            pass
            
        #convert Magnetic Channels for long period surveys
        if magtype=='lp' and lpyn=='n':
            #Convert BX
            if fnmatch.fnmatch(filenames[ii],'*.BX'):
                bxfid=file(filenames[ii],'r')
                bxlines=bxfid.readlines()
                if bxlines[0].find('.')>=0:
                    print 'Found decimal point '+filenames[ii]+'. Check File.'
                    bxyn=input('Still convert? (y/n) as a string')
                    if bxyn=='n':
                            bxfid.close()
                    else:
                            bxconv=convertlpB(bxlines,dlgain)
                            bxfid.close()
                            bxconvlst=[str(bxconv[ii])+'\n' for ii in range(len(bxconv))]
                            bxfidn=file(filenames[ii],'w')
                            bxfidn.writelines(bxconvlst)
                            bxfidn.close()
            #convert BY
            elif fnmatch.fnmatch(filenames[ii],'*.BY'):
                byfid=file(filenames[ii],'r')
                bylines=byfid.readlines()
                if bylines[0].find('.')>=0:
                    print 'Found decimal point '+filenames[ii]+'. Check File.'
                    byyn=input('Still convert? (y/n) as a string')
                    if byyn=='n':
                            byfid.close()
                    else:
                            byconv=convertlpB(bylines,dlgain)
                            byfid.close()
                            byconvlst=[str(byconv[ii])+'\n' for ii in range(len(byconv))]
                            byfidn=file(filenames[ii],'w')
                            byfidn.writelines(byconvlst)
                            byfidn.close()
                else:
                    byconv=convertlpB(bylines,dlgain)
                    byfid.close()
                    byconvlst=[str(byconv[ii])+'\n' for ii in range(len(byconv))]
                    byfidn=file(filenames[ii],'w')
                    byfidn.writelines(byconvlst)
                    byfidn.close()
            #convert BZ
            elif fnmatch.fnmatch(filenames[ii],'*.BZ'):
                bzfid=file(filenames[ii],'r')
                bzlines=bzfid.readlines()
                if bzlines[0].find('.')>=0:
                    print 'Found decimal point '+filenames[ii]+'. Check File.'
                    bzyn=input('Still convert? (y/n) as a string')
                    if bzyn=='n':
                            bzfid.close()
                    else:
                            bzconv=convertlpB(bzlines,dlgain,zadj=zadj)
                            bzfid.close()
                            bzconvlst=[str(bzconv[ii])+'\n' for ii in range(len(bzconv))]
                            bzfidn=file(filenames[ii],'w')
                            bzfidn.writelines(bzconvlst)
                            bzfidn.close()
                else:
                    bzconv=convertlpB(bzlines,dlgain,zadj=zadj)
                    bzfid.close()
                    bzconvlst=[str(bzconv[ii])+'\n' for ii in range(len(bzconv))]
                    bzfidn=file(filenames[ii],'w')
                    bzfidn.writelines(bzconvlst)
                    bzfidn.close()
            else:
                pass
        else:
            pass
    
    







