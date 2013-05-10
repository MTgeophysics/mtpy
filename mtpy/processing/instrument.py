
"""
mtpy/processing/instrument.py

Functions for the removal of instrument response from time series data.

Either working on ASCII data or on miniSeed



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import os.path as op

import copy

import mtpy.processing.filter as MTfi
import mtpy.utils.exceptions as MTex

#=================================================================

def correct_for_instrument_response(data, samplingrate, responsedata):
    """Correct input time series for instrument response.
        Instr.Resp. is given as 3 column array: frequency, real part, imaginary part

        The given section is demeaned, window tapered, zero padded, (potentially bandpassed with the extreme frequencies of the response frequency axis), FFT-ed, "deconvolved" (straight division by the array values or interpolated values inbetween - in frequency domain), re-transformed, mean-re-added and returned.

    """

    datamean = np.mean(data)
    data -= datamean
    
    freqmin = responsedata[0,0]
    freqmax = responsedata[-1,0]
    #print freqmin,freqmax
    
    N = len(data)
    if N < 1:
        raise MTex.MTpyError_ts_data('Error - Length of TS to correct is zero!')

    #use double sided cosine taper function
    window = MTfi.tukey(N, 0.2)

    tapered_data = data * window


    #check, if length is directly equal to power of 2:
    if np.log2(N)%1 == 0:
        next2power = int(np.log2(N))
    else:
        next2power = int(np.log2(N)) + 1

    #zero pad data for significantly faster fft - NO, Numpy does not do that automatically
    padded_data = np.zeros((2**next2power))
    padded_data[:len(tapered_data)] = tapered_data


    #bandpass data for excluding all frequencies that are not covered by the known instrument response
    bp_data = padded_data #butter_bandpass_filter(padded_data, freqmin, freqmax, samplingrate)

    #get the spectrum of the data 
    data_spectrum = np.fft.rfft(bp_data)
    data_freqs = np.fft.fftfreq(len(bp_data),1./samplingrate)

    #and the same for the instrument
    instr_spectrum = responsedata[:,1] + np.complex(0,1) * responsedata[:,2]
    instr_freqs = responsedata[:,0]

    lo_mags = []
    lo_freqs = []

    #return data_freqs,data_spectrum,instr_freqs,instr_spectrum

    #now correct for all frequencies on the data_freqs-axis:
    corrected_spectrum = np.zeros((len(data_spectrum)),'complex')
    for i in range(len(data_spectrum)):
        freq = data_freqs[i]
        spec = data_spectrum[i]

        #this is effectively a boxcar window - maybe to be replaced by proper windowing function ?
        if not (freqmin <= np.abs(freq) <= freqmax):
            print 'no instrument response in this frequency range - spectrum set to zero here: ', freq
            corrected_spectrum[i] = 0#spec
            #lo_mags.append([np.abs(spec),np.abs(corrected_spectrum[i]), 0 ])
            #lo_freqs.append([freq,0,0,0,0,0,0])
            continue

        #find the appropriate frequencies ( since the current freq-value is most likely inbetween two values on the instr_freqs-axis)
        #get the respective value by linear interpolation

        #find the value closest to the current freq, assume it's lower
        closest_lower = np.abs(freq-instr_freqs).argmin()
        
        #if it coincides with the highest frequency/last entry:
        if closest_lower == len(instr_freqs)-1:
            factor = data_spectrum[closest_lower]
        # or the lowest
        elif closest_lower == 0:
            factor = data_spectrum[0]
        else:
            #in case the closest frequency value is not lower but higher, take the freq value below as lower bound for the interval:        
            if instr_freqs[closest_lower] > freq:
                closest_lower -= 1
            
            #define the interval:
            instrfreq1 = instr_freqs[closest_lower]
            instrfreq2 = instr_freqs[closest_lower+1]
            instrfactor1 = instr_spectrum[closest_lower]
            instrfactor2 = instr_spectrum[closest_lower+1]
            intervallength = instrfreq2 - instrfreq1
            weight = (freq-instrfreq1)/intervallength

            factor = weight * instrfactor2 + (1-weight) * instrfactor1
            #lo_freqs.append([freq,instrfreq1,instrfreq2,weight,instrfactor1,instrfactor2,factor])

        
        #finally correct the data for the instrument influence by dividing the complex data spectral value by the instrument response value: 
        corrected_spectrum[i] = spec / factor
        #lo_mags.append([np.abs(spec),np.abs(corrected_spectrum[i]), factor ])
    #invert into time domain
    correctedTS = np.fft.irfft(corrected_spectrum)

    #cut the zero padding
    correctedTS = correctedTS[:N]

    #re-attach the mean
    correctedTS += datamean

    return correctedTS#, lo_mags,lo_freqs