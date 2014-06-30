#!/usr/bin/env python

"""
mtpy/processing/calibration.py

This modules contains functions for the calibration of raw time series. 

The various functions deal with the calibration of data from 
fluxgates, coil, dipoles,... The calibration depends on  the instrument 
as well as on the respective data logger. 

All information needed for the calibration must be provided by a configuration
file. This has to contain station names as section headers. Each section 
must contain a set of mandatory keywords. The keyword list is defined in the 
function mtpy.processing.filehandling.read_survey_configfile()

For the calibration of one or more station data, use the appropriate
scripts from the mtpy.utils subpackage. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import sys, os
import os.path as op
import copy

import  mtpy.utils.exceptions as MTex
import mtpy.utils.configfile as MTcf

#=================================================================

list_of_channels = ['ex','ey','bx','by','bz']


list_of_bfield_loggers = MTcf.dict_of_allowed_values_bfield['B_logger_type']
list_of_bfield_instruments = MTcf.dict_of_allowed_values_bfield['B_instrument_type']
list_of_efield_loggers = MTcf.dict_of_allowed_values_efield['E_logger_type']
list_of_efield_instruments = MTcf.dict_of_allowed_values_efield['E_instrument_type']

list_of_loggers = list(set(list_of_bfield_loggers + list_of_efield_loggers))
list_of_instruments = list(set(list_of_bfield_instruments+list_of_efield_instruments))

# section for amplification and scaling factors:


dict_of_calibration_factors_volt2nanotesla = {'fluxgate': 70000/0.1, 'coil': 1.}

#...dict_of_instrument_amplification = {'electrodes' :10. , 'fluxgate' = 1., 'coil': 1.}
#dict_of_channel_amplification = {'ex':1, 'ey':1.,'bx':1. ,'by':1., 'bz': 0.5}

dict_of_bz_instrument_amplification = {'edl': 0.5, 'elogger': 1.}

dict_of_EDL_gain_factors = {'high': 10., 'low': 1., 'verylow': 0.4,
                             str(10): 10., str(1): 1. , str(0.4): 0.4}

list_of_elogger_gain_factors = [11.,1]

#dict_of_efield_amplification = {'edl': 10., 'elogger': 1.}

#=================================================================




#=================================================================


def calibrate(raw_data, field, instrument, logger,dipole_length=1.,
                calibration_factor=1., amplification=1., gain=1., offset = 0.):

    """
    Convert a given time series from raw data (voltage) 
    into field strength amplitude values. 
    The dipole length factor must be set to 1 for the calibration of B-field data.

    The output is a time series of field strength values in basic units: 
    Tesla for the B-field and V/m for the E-field.

    input:
    - 1D time series (list or numpy array)
    - type of measured field (E/B)
    - dipole length (if value is negative, the data will change sign )
    - actual calibration factor ('Volts/Tesla per count')
    - instrument (fluxgate, coil, electrodes, bartington)
    - amplification factor for the instrument
    - logger (EDL, elogger, other)
    - gain factor of the logger
    - offset of the counts (will be subtracted)

    output:
    - 1D time series (numpy array)

    """

    raw_data = np.array(raw_data) - offset

    units_per_count = calibration_factor

    _data_instrument_consistency_check(raw_data, field, dipole_length,
                                     instrument, amplification, logger, gain )
    

    # converting counts into units, including  
    # - transistion from voltage to field
    # - correcting for instrument amplification
    # - correcting for logger amplification (gain)
    # - raw conversion factor counts2unit (Volt or Tesla)
    data = raw_data * units_per_count / dipole_length / amplification / gain 


    return data

#=================================================================

def EDL_e_field(data, edl_gain, dipole, instrument_amplification):
    """
    Convert EDL output (channels MTex and EY) into E field values.

    input:
    - time series of field values in microvolt (standard EDL output)
    - gain factor set at the EDL 
    - dipole length in meters
    - instrument amplification

    output:
    - time series of E field values in V/m

    """

    # Since the conversion is straight from volt into V/m, no further calibration factor is needed
    e_field = calibrate(data, 'e', 'electrodes', 'edl', dipole_length = dipole, 
                        calibration_factor = 1., 
                        amplification = instrument_amplification, 
                        gain = edl_gain, offset = 0.)

    return e_field

#=================================================================

def EDL_b_field(data, edl_gain, instrument , instrument_amplification):
    """
    Convert EDL output (channels BX, BY and BZ) into B field values.

    input:
    - time series of field values in microvolt (standard EDL output)
    - gain factor set at the EDL 
    - type of instrument
    - instrument amplification

    output:
    - time series of B field values in nanotesla (!!)

    """

    #setting the calibration factor for converting voltage into electric field strength:
    #nanotesla_per_volt = dict_of_calibration_factors_volt2nanotesla[instrument]

    #nanotesla_per_microvolt = nanotesla_per_volt / (10 ** 6)

    b_field = calibrate(data, 'b', instrument, 'edl', dipole_length = 1., 
                        calibration_factor = 1. , 
                        amplification = instrument_amplification,
                        gain = edl_gain, offset = 0.)

    return b_field

#=================================================================

def elogger_e_field(data, elogger_gain, dipole, instrument_amplification):
    """
    Convert elogger output (channels MTex and EY) into E field values.

    input:
    - time series of field values in microvolt (standard EDL output)
    - gain factor set at the elogger 
    - dipole length in meters
    - instrument amplification

    output:
    - time series of E field values in microvolt/meter

    """
    # Since the conversion is straight from volt into V/m, no further calibration factor is needed
    e_field = calibrate(data, 'e', 'electrodes', 'elogger',dipole_length = dipole, 
                        calibration_factor = 1., 
                        amplification = instrument_amplification, 
                        gain = elogger_gain, offset = 0.)

    return e_field

#=================================================================

def calibrate_file(filename, outdir, instrument, instrument_amplification, 
                    logger, gain, dipole, stationname, channel, latitude, 
                    longitude, elevation, offset = 0 ):
    """
    Calibrate data from one given file and store the output to another file.
    If the channel is not given explicitly, it's taken from the filename suffix.

    E field values will be present in microvolt/meter
    B fields are given in nanotesla 

    input:
    - data file name
    - foldername for saving the output
    - instrument type
    - instrument amplification factor
    - data logger type
    - logger gain factor
    - station name
    - channel  

    """
    time_axis = None

    if not instrument.lower() in list_of_instruments:
        raise MTex.MTpyError_inputarguments('instrument type not known')

    if not logger.lower() in list_of_loggers:
        raise MTex.MTpyError_inputarguments('data logger type not known')



    if not op.isfile(filename):
        raise MTex.MTpyError_inputarguments('data file not existing')

    infile_base = op.basename(filename)

    try:
        data_in = np.loadtxt(filename)
    except:
        raise MTex.MTpyError_inputarguments('cannot read data file')


    data_out = copy.copy(data_in)

    #read in first line of input file, checking, if header line exists
    FH = open(filename,'r')
    firstline = FH.readline().strip()
    FH.close()


    if np.size(data_in.shape) > 1:
        if data_in.shape[1] > 1:
            #at least 2 columns - assume, first is time, second data - ignore, if there are more
            time_axis = data_in[:,0]
            data_in = data_in[:,1]



    if not op.isdir(outdir):
        try:
            os.makedirs(outdir)
        except:
            raise MTex.MTpyError_inputarguments('output directory is not '
                                            'existing and cannot be generated')

    if channel == None:
        channel = filename[-2:].lower()
        
    if not channel in list_of_channels:
        raise MTex.MTpyError_inputarguments('wrong channel specification')
   
    field = channel[0]
    
    #print 'channel:...........',channel, field
    #print 'read file',filename ,'wrote file...'
    #return

    #separate way for B and E fields here:
    if field == 'e':

        if dipole <= 1:
            print 'Check dipole length value ! - It is highly improbable to '\
                    'have a 1 meter dipole!!'
            
            answer = raw_input('\t\tContinue anyway? [y/N] \n')
            
            if not answer[0].lower() == 'y':
                sys.exit('Calibration process interrupted by user input!')


        instrument = 'electrodes'

        logger = logger.lower()

        if logger == 'elogger':

            if not type(gain) in [float, int]:#list_of_elogger_gain_factors:
                raise MTex.MTpyError_inputarguments('invalid gain for elogger:'
                                                    ' {0}'.format(gain))

            #instrument_amplification = dict_of_efield_amplification[logger]

            outfile_data = elogger_e_field(data_in, gain, dipole,
                                             instrument_amplification)
        
        elif logger == 'edl':

            if not type(gain) in [float, int, str]:
                raise MTex.MTpyError_inputarguments('invalid gain for EDL: '
                                                    '{0}'.format(gain))

            #instrument_amplification = dict_of_efield_amplification[logger]

            if type(gain) == str:
                EDLgain = dict_of_EDL_gain_factors[gain] 
            else:
                EDLgain = float(gain)

            outfile_data = EDL_e_field(data_in, EDLgain, dipole,
                                         instrument_amplification)
    
        dataunit ='microvoltpermeter'

    #B-field part 
    elif field == 'b':
        instrument = instrument.lower()
        if not instrument in list_of_bfield_instruments:
            raise MTex.MTpyError_inputarguments('invalid instrument for B-'
                                                'field measurements')


        logger = logger.lower()

        if not logger in list_of_bfield_loggers:
            raise MTex.MTpyError_inputarguments('invalid logger for B-field'
                                                ' measurements')


        #instrument_amplification = 1.
 
        #calibration_factor = dict_of_calibration_factors_volt2nanotesla[instrument]


        if logger == 'edl':

            if not type(gain) in [float,int,str]:
                raise MTex.MTpyError_inputarguments('invalid gain: '
                                                    '{0}'.format(gain))

            if type(gain) == str:
                EDLgain = dict_of_EDL_gain_factors[gain] 
            else:
                EDLgain = float(gain)


            if instrument == 'fluxgate' and channel == 'bz':
                instrument_amplification *= dict_of_bz_instrument_amplification[logger]


            outfile_data = EDL_b_field(data_in, EDLgain, instrument, 
                                        instrument_amplification)

        dataunit = 'nanotesla'



    newbasename = '{0}_{1}.{2}'.format(op.splitext(infile_base)[0], dataunit, 
                                infile_base.split('.')[-1].lower())


    #set up output file
    outfile = op.join(outdir, newbasename)
    
    additional_header_info = ' {0} {1:02.5f} {2:03.5f} {3:.1f} \n'.format(
                                    dataunit, latitude, longitude, elevation)

    if firstline[0][0] == '#':
        newfirstline = firstline + additional_header_info

    else:
        newfirstline = '# {0} {1} {2}'.format(stationname, channel, 
                                                additional_header_info)


    if time_axis != None:
        data_out[:,1] = outfile_data
    else:
        data_out = outfile_data

    Fout = open(outfile,'w')

    Fout.write(newfirstline)
    np.savetxt(Fout,data_out,fmt='%.8e')
    Fout.close()

    print 'read file',filename ,'  ->  wrote file %s'%(outfile)
    

#=================================================================

def _data_instrument_consistency_check(data, field, dipole_length, instrument, 
                                        amplification, logger, gain):
    """
    Check, if the input values given  for the calibration make any sense at all.

    """

    if len(data) == 0 :
        raise MTex.MTpyError_ts_data( 'no data provided for calibration' )

    if not field.lower() in ['e','b']:
        raise MTex.MTpyError_inputarguments( 'Field must be E or B' )

#    if float(dipole_length) <= 0:
#        raise MTpyError_inputarguments( 'Dipole length must be positive' )

    if float(amplification) <= 0:
        raise MTex.MTpyError_inputarguments('Amplification factor must be positive')

    if float(gain) <= 0:
        raise MTex.MTpyError_inputarguments('Instrument gain must be positive')

    try:
        if not logger.lower() in list_of_loggers:
            raise
    except:
        raise MTex.MTpyError_inputarguments( 'wrong choice of logger')

    try:
        if not instrument.lower() in list_of_instruments:
            raise
    except:
        raise MTex.MTpyError_inputarguments( 'wrong choice of instrument')

    if field.lower == 'b':
        if logger.lower() == 'elogger':
            raise MTex.MTpyError_inputarguments( 'wrong choice of logger')
        if instrument.lower() == 'electrodes':
            raise MTex.MTpyError_inputarguments( 'wrong choice of instrument')
        if not float(dipole_length) == 1:
            raise MTex.MTpyError_inputarguments( 'Dipole length must be "1" for'
                                                ' B-field calibration')
           

    if field.lower == 'e':
        if not instrument.lower() == 'electrodes':
            raise MTex.MTpyError_inputarguments( 'wrong choice of instrument')






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
                    if fnmatch.fnmatch(filename,'*.MTex'):
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



def convertCounts2Units(filenames,eyn='n',lpyn='n',egain=1.0,dlgain=1.0,
                        exlen=100.,eylen=100.,magtype='lp',zadj=2):

    """
    convertCounts2Units(filenames,eyn='n',lpyn='n',egain=1.0,dlgain=1.0,
    exlen=100.,eylen=100., magtype='lp',zadj=2) 
    will convert a set of files given by filenames with parameters defined:

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
            #convert MTex chanel
            if fnmatch.fnmatch(filenames[ii],'*.MTex'):
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
    
    







