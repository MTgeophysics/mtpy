#!/usr/bin/env python

"""
mtpy/processing/calibration.py

This modules contains functions for the calibration of raw time series. 

The various functions deal with the calibration of data from fluxgates, coil, dipoles,...
The calibration depends on  the instrument as well as on the respective data logger. 

All information needed for the calibration must be provided by a configuration file. This has to 
consist of station names as section headers. Each section must contain a set of mandatory keywords.
The keyword list is defined in the function mtpy.processing.filehandling.read_survey_configfile()

For calling a batch calibration rather than just one file, use the appropriate scripts from the mtpy.utils subpackage. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import glob
import os.path as op
import calendar
import time
import copy


import  mtpy.utils.exceptions as MTex

#=================================================================

list_of_channels = ['ex','ey','bx','by','bz']

list_of_loggers = ['edl','elogger']
list_of_instruments = ['electrodes','fluxgate', 'coil']
list_of_bfield_loggers = ['edl']
list_of_bfield_instruments = ['fluxgate', 'coil']

# section for amplification and scaling factors:


dict_of_calibration_factors_volt2nanotesla = {'fluxgate': 70000/0.1, 'coil': 1.}

#...dict_of_instrument_amplification = {'electrodes' :10. , 'fluxgate' = 1., 'coil': 1.}
#dict_of_channel_amplification = {'ex':1, 'ey':1.,'bx':1. ,'by':1., 'bz': 0.5}

dict_of_bz_instrument_amplification = {'edl': 0.5, 'elogger': 1.}

dict_of_EDL_gain_factors = {'high': 10., 'low': 1., 'verylow': 0.4, str(10): 10., str(1): 1. , str(0.4): 0.4}

list_of_elogger_gain_factors = [11.,1]

dict_of_efield_amplification = {'edl': 1., 'elogger': 10.}

#=================================================================




#=================================================================


def calibrate(raw_data, field, instrument, logger,dipole_length=1.,
                calibration_factor = 1.,  amplification = 1.,  gain = 1., offset = 0.):

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

    _data_instrument_consistency_check(raw_data, field, dipole_length, instrument, amplification, logger, gain )
    

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
                        calibration_factor = 1., amplification = instrument_amplification, 
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
    nanotesla_per_volt = dict_of_calibration_factors_volt2nanotesla[instrument]

    nanotesla_per_microvolt = nanotesla_per_volt / (10 ** 6)

    b_field = calibrate(data, 'b', instrument, 'edl', dipole_length = 1., 
                        calibration_factor = 1. , amplification = instrument_amplification,
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
                        calibration_factor = 1., amplification = instrument_amplification, 
                        gain = elogger_gain, offset = 0.)

    return e_field

#=================================================================

def calibrate_file(filename, outdir, instrument, logger, gain, dipole, stationname, channel, latitude, longitude, elevation, offset = 0 ):
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
            raise MTex.MTpyError_inputarguments('output directory is not existing and cannot be generated')

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
            print 'Check dipole length value ! - It is highly improbable to have a 1 meter dipole!!'
            answer = raw_input('\t\tContinue anyway? [y/N] \n')
            if not answer[0].lower() == 'y':
                sys.exit('Calibration process interrupted by user input!')


        instrument = 'electrodes'

        logger = logger.lower()

        if logger == 'elogger':

            if not type(gain) in [float, int]:#list_of_elogger_gain_factors:
                raise MTex.MTpyError_inputarguments('invalid gain for elogger: {0}'.format(gain))

            instrument_amplification = dict_of_efield_amplification[logger]

            outfile_data = elogger_e_field(data_in, gain, dipole, instrument_amplification)
        
        elif logger == 'edl':

            if not type(gain) in [float, int, str]:
                raise MTex.MTpyError_inputarguments('invalid gain for EDL: {0}'.format(gain))

            instrument_amplification = dict_of_efield_amplification[logger]

            if type(gain) == str:
                EDLgain = dict_of_EDL_gain_factors[gain] 
            else:
                EDLgain = float(gain)

            outfile_data = EDL_e_field(data_in, EDLgain, dipole, instrument_amplification)
    
        dataunit ='microvoltpermeter'

    #B-field part 
    elif field == 'b':
        instrument = instrument.lower()
        if not instrument in list_of_bfield_instruments:
            raise MTex.MTpyError_inputarguments('invalid instrument for B field measurements')


        logger = logger.lower()

        if not logger in list_of_bfield_loggers:
            raise MTex.MTpyError_inputarguments('invalid logger for B field measurements')


        instrument_amplification = 1.
 
        calibration_factor = dict_of_calibration_factors_volt2nanotesla[instrument]


        if logger == 'edl':

            if not type(gain) in [float,int,str]:
                raise MTex.MTpyError_inputarguments('invalid gain: {0}'.format(gain))

            if type(gain) == str:
                EDLgain = dict_of_EDL_gain_factors[gain] 
            else:
                EDLgain = float(gain)


            if instrument == 'fluxgate' and channel == 'bz':
                instrument_amplification = dict_of_bz_instrument_amplification[logger]


            outfile_data = EDL_b_field(data_in, EDLgain, instrument, instrument_amplification)

        dataunit = 'nanotesla'



    newbasename = '%s_%s.%s'%(op.splitext(infile_base)[0], dataunit, infile_base.split('.')[-1].lower())


    #set up output file
    outfile = op.join(outdir, newbasename)
    
    additional_header_info = ' %s %02.5f %03.5f %.1f \n'%(dataunit, latitude, longitude, elevation)

    if firstline[0][0] == '#':
        newfirstline = firstline + additional_header_info

    else:
        newfirstline = '# %s %s %s'%(stationname, channel, additional_header_info)


    if time_axis != None:
        data_out[:,1] = outfile_data
    else:
        data_out = outfile_data

    Fout = open(outfile,'w')

    Fout.write(newfirstline)
    np.savetxt(Fout,data_out)
    Fout.close()

    print 'read file',filename ,'  ->  wrote file %s'%(outfile)
    

#=================================================================

def _data_instrument_consistency_check(data, field, dipole_length, instrument, amplification, logger, gain):
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
        raise MTex.MTpyError_inputarguments( 'Amplification factor must be positive' )

    if float(gain) <= 0:
        raise MTex.MTpyError_inputarguments( 'Instrument gain must be positive' )

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
            raise MTex.MTpyError_inputarguments( 'Dipole length must be 1 for B-field calibration')
           

    if field.lower == 'e':
        if not instrument.lower() == 'electrodes':
            raise MTex.MTpyError_inputarguments( 'wrong choice of instrument')


#----------------------------------------------------------
#calibration/instrument correction factors for coils from old fortran code:
"""

      data fcal/
      -1.698970004,-1.143414469,-1,-1,-0.950000176,-0.918367218,&
     -0.900000142,-0.84999986,-0.836734836,-0.800000053,-0.755101948,-0.7500001,&
     -0.700000069,-0.673469365,-0.650000027,-0.600000075,-0.591836722,-0.587858949,&
     -0.549999989,-0.510204025,-0.499999953,-0.449999987,-0.428571396,-0.399999968,&
     -0.349999992,-0.346938744,-0.299999942,-0.265306158,-0.250000019,-0.199999962,&
     -0.183673507,-0.14999999,-0.102040773,-0.099999964,-0.04999997,-0.032303369,&
     -0.020408153,0,0.049999824,0.061224339,0.099999858,0.142856988,0.15000014,&
     0.199999947,0.224489812,0.2499999,0.299999931,0.306122525,0.349999973,0.387755085,&
     0.399999925,0.450000011,0.469387799,0.500000047,0.523252208,0.550000013,0.551020371,&
     0.600000032,0.632653035,0.650000008,0.699999971,0.714285741,0.749999981,0.795918373,&
     0.800000038,0.85000001,0.877551017,0.900000036,0.95000003,0.959183684,1,1.04081627,&
     1.049999824,1.078807954,1.099999858,1.12244886,1.15000014,1.199999947,1.204081709,&
     1.2499999,1.285714347,1.299999931,1.349999973,1.367346974,1.399999925,1.448979638,&
     1.450000011,1.500000047,1.530612274,1.550000013,1.600000032,1.612244891,1.634363291,&
     1.650000008,1.693877519,1.699999971,1.749999981,1.775510181,1.799999969,1.85000001,&
     1.857142873,1.899999981,1.938775523,1.949999981,2,2.020408108,2.049999824,2.099999858,&
     2.102040742,2.15000014,2.183673479,2.189918782,2.199999947,2.2499999,2.26530613,&
     2.299999931,2.346938877,2.349999973,2.399999925,2.428571462,2.450000011,2.500000047,&
     2.510204143,2.550000013,2.591836742,2.600000032,2.650000008,2.673469354,2.699999971,&
     2.745474455,2.749999981,2.755102039,2.799999969,2.836734729,2.85000001,2.899999981,&
     2.918367361,2.949999981,3,3,3.049999824,3.099999858,3.15000014,3.199999947,3.2499999,&
     3.299999931,3.301029996,3.349999973,3.399999925,3.450000011,3.500000047,3.550000013,&
     3.600000032,3.650000008,3.700000058,3.749999981,3.799999969,3.85000001,3.899999981,&
     3.95000003,4./

      data bbreal/0.08751027,0.991550192,1.932381149,1.912288443,2.411190921,2.754305947,&
    2.988249574,3.739075611,3.992621642,4.702044979,5.762355864,5.904239902,&
    7.366392861,8.35493777,9.167182541,11.37253817,11.83207151,12.0305058,14.12676304,&
    16.79331044,17.43937534,21.45092302,23.43415045,26.23873616,31.92253211,32.37847604,&
    38.58290002,43.87056582,46.17393689,54.79515782,57.83382015,64.27013599,74.0963827,&
    74.53149823,85.32784198,89.34121468,92.11352358,96.50259747,107.7184107,110.4354569,&
    118.6086542,127.5768849,129.0051705,138.649938,143.1004377,147.3925129,155.2149482,&
    156.1498695,161.9971905,166.5446629,167.8300754,172.7686301,174.4063644,176.8930595,&
    178.5916332,180.3623584,180.4559098,183.1665899,184.7034306,185.4936039,187.3538583,&
    187.6872816,188.8728864,189.946305,190.075909,190.9736596,191.4714411,191.7815779,  &
    192.4726385,192.560534,192.9090899,193.1581461,193.2561537,193.1082516,193.2773939, &
    193.4418833,193.9895741,193.6051332,193.4131278,193.6173632,193.5920319,193.8826849,&
    193.4712867,193.4396524,193.4834029,192.4605288,193.2654432,193.4109519,192.2431426,&
    192.777544,191.0945088,189.9369073,190.2595909,192.8336427,199.0939097,212.0531716, &
    189.9708944,187.0299701,185.7458902,184.1742584,184.375308,180.5125671,178.6752186, &
    180.6873016,173.8307254,170.3666283,170.488904,162.1373034,162.9634387,156.2315774, &
    150.1511192,146.7016962,144.2943709,132.6565939,127.6704154,123.3076614,101.2930001,&
    101.0248442,78.11257604,63.50840658,52.1427182,21.54794583,15.20055558,-13.54964865,&
    -43.52213147,-50.04713873,-90.58974186,-109.4390568,-129.6500829,-158.6440058,      &
    -160.4292032,-161.5623377,-175.6708413,-171.6665928,-165.8740514,-120.0742394,      &
    -92.9019606,-36.20424614,59.81240942,59.42844342,88.50219043,44.38480271,7.867028419,&
    -3.01091948,-2.528702377,-1.91966993,-1.836025186,-0.90110632,-0.365153733,-0.12709547,&
    0.018406311,-0.243037017,0.001443097,-0.10361397,0.020712788,0.012486207,0.001846097,&
    -0.069616301,0.120219932,-0.02063392,0.030414613/
 
      data bbimag/3.848183107,13.75145843,19.18515322,19.0494526,21.29797644,22.89311286,&
     23.82129191,26.62749805,27.25472718,29.73288341,32.76611495,33.15811302,          &
     36.91149455,39.16080683,41.00270904,45.44852767,46.38882319,46.61936714,          &
     50.24291696,54.41279185,55.32325763,60.65569969,62.97454721,66.14951417,          &
     71.71683201,72.14938494,77.20477585,80.72751174,82.35644874,87.03971189,          &
     88.45874146,90.99444582,93.839937,94.02914439,95.94383058,96.27605347,            &
     96.32786209,96.5822972,95.95040301,95.53757009,94.01759997,91.40092668,           &
     90.95964536,86.88860593,84.67137233,82.13590066,76.80060822,76.15185696,          &
     71.12511551,66.67172106,65.26922257,59.45221423,57.08444434,53.74853734,          &
     51.10310655,48.22060587,48.0925069,42.98022352,39.73423938,38.0106161,            &
     33.43178049,31.89999825,28.91367384,25.1806512,24.75078285,20.9578769,            &
     18.90476185,17.30739285,13.94779668,13.29422881,10.66639781,8.141509683,          &
     7.521084244,5.815522234,4.444986422,2.891989598,1.249788826,-1.72299622,          &
     -1.191850094,-4.572580358,-6.802236755,-8.085339895,-10.90637053,-12.14682078,    &
     -14.3996928,-18.08595291,-18.44309206,-21.83746531,-25.27133007,-25.96587104,     &
     -30.08704364,-30.52136955,-31.75743372,-37.11418625,-62.25344718,-30.94718289,    &
     -47.52625026,-50.39690761,-52.23813562,-59.79692721,-59.473502,-66.1242099,       &
     -72.90920343,-75.19773176,-65.09391549,-88.50581562,-92.84724945,-104.2954989,    &
     -105.6268134,-115.5907706,-117.8152969,-121.3235606,-127.352269,-140.4479552,     &
     -142.7664632,-151.6617244,-162.6776787,-164.3943341,-175.6993301,-180.3320503,    &
     -183.7880545,-189.1793755,-189.268938,-189.6716401,-182.525685,-179.457892,       &
     -162.8238423,-149.6895162,-131.5341539,-89.95650577,-85.07450402,-79.57454501,    &
     -25.2006442,26.8307038,45.60399439,112.6015387,131.3457142,148.0418325,114.237942,&
     115.1202867,24.79314617,-25.84466571,-23.34388932,-8.99707201,-4.323179975,&
     -0.988029522,-1.074438893,-0.486058304,0.055081554,0.078239812,-0.011510172,&
     -0.045552074,0.015427292,0.091207003,0.116829919,0.010301625,-0.000215012,&
     -0.020632919,0.026161701,0.005712443,0.026240797/

"""

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
    
    







