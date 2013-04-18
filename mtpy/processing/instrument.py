
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
from scipy.signal import butter, lfilter


import mtpy.utils.mseed as MTmseed
#import obspy.mseed as Omseed
import mtpy.utils.exceptions as MTexceptions

#=================================================================


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def correct_for_instrument_resonse(data, samplingrate, responsedata):
    """Correct input time series for instrument response.
        Instr.Resp. is given as 3 column array: frequency, real part, imaginary part

        The given section is demeaned, window tapered, zero padded, bandpassed(with the extreme frequencies of the response frequency axis), deconvolved (by straight division by the array values or interpolated values inbetween), re-transformed, mean-re-added and returned.

    """

    datamean = np.mean(data)
    data -= datamean
    
    freqmin = responsedata[0,0]
    freqmax = responsedata[-1,0]
    
    N = len(data)
    
    window = np.blackman(N)

    tapered_data = window * data

    #check, if length is directly equal to power of 2:
    if np.log2(N)%1 == 0:
        next2power = int(np.log2(N))
    else:
        next2power = int(np.log2(N)) + 1

    #zero pad data for significantly faster fft - YES, Numpy does not do that automatically
    padded_data = zeros((2**next2power))
    padded_data[:len(tapered_data)] = tapered_data

    #bandpass data for excluding all frequencies that are not covered by the known instrument response
    bp_data = butter_bandpass_filter(padded_data, freqmin, freqmax, samplingrate)

    #get the spectrum of the data 
    data_spectrum = np.fft.rfft(bp_data)
    data_freqs = np.fft.fftfreq(len(bp_data),1./samplingrate)

    #and the same for the instrument
    instr_spectrum = responsedata[:,1] + np.complex(0,1) * responsedata[:,2]
    instr_freqs = responsedata[:,0]

    #now correct for all frequencies on the data_freqs-axis:
    corrected_spectrum = zeros((len(data_spectrum)),'complex')
    for i in range(len(data_spectrum)):
        freq = data_freqs[i]
        spec = data_spectrum[i]

        factor = 1.

        #find the appropriate frequency most likely inbetween two values on the instr_freqs-axis
        #get the respective value by linear interpolation
        closest_lower = np.abs(freq-instr_freqs).argmin()
        
#Todo: catch frequencies outside the range

        if closest_lower == len(instr_freqs)-1:
            factor = data_spectrum[closest_lower]
        elif closest_lower == 0:
            factor = data_spectrum[0]
        else:            
            if instr_freqs[closest] > freq:
                closest_lower -= 1

            instrfreq1 = instr_freqs[closest_lower]
            instrfreq2 = instr_freqs[closest_lower+1]
            instrfactor1 = data_spectrum[closest_lower]
            instrfactor2 = data_spectrum[closest_lower+1]


    #invert into time domain
    correctedTS = fft.irfft(corrected_spectrum)

    #cut the zero padding
    correctedTS = correctedTS[:N]

    #re-attach the mean
    correctedTS += datamean

    return correctedTS