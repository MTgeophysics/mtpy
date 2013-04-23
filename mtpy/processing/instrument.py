
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
from scipy.signal import butter, lfilter, buttord


import mtpy.utils.mseed as MTmseed
#import obspy.mseed as Omseed
import mtpy.utils.exceptions as MTexceptions

#=================================================================


def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    print high, low
    if high >=1. and low == 0.:
        b = np.array([1.])
        a = np.array([1.])

    elif high < 0.95 and low > 0. :
        wp = [1.05*low,high-0.05]
        ws = [0.95*low,high+0.05]
 
        print wp,ws
        order,wn = buttord(wp,ws,0., 30.)
        b, a = butter(order, wn, btype='band')
    
    elif high>= 0.95:
        print 'highpass',low,1.2*low,0.8*low
        order,wn = buttord( 15*low,0.05*low,gpass=0.0, gstop=10.0)
        print order,wn
        b, a = butter(order, wn, btype='high')
    elif low <= 0.05:
        print 'lowpass',high
        order,wn = buttord( high-0.05,high+0.05,gpass=0.0, gstop=10.0)
        b, a = butter(order, wn, btype='low')

    return b, a



def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    print b,a
    y = lfilter(b, a, data)
    return y


def tukey(window_length, alpha=0.5):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width \alpha * N / 2
    that is convolved with a rectangle window of width (1 - \alpha / 2). At \alpha = 1 it becomes rectangular, and
    at \alpha = 0 it becomes a Hann window.
 
    We use the same reference as MATLAB to provide the same results in case users compare a MATLAB output to this function
    output
 
    Reference
    ---------
 
    http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html
 
    '''
    # Special cases
    if alpha <= 0:
        return np.ones(window_length) #rectangular window
    elif alpha >= 1:
        return np.hanning(window_length)
 
    # Normal case
    x = np.linspace(0, 1, window_length)
    w = np.ones(x.shape)
 
    # first condition 0 <= x < alpha/2
    first_condition = x<alpha/2
    w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))
 
    # second condition already taken care of
 
    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x>=(1 - alpha/2)
    w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))
 
    return w    

def correct_for_instrument_response(data, samplingrate, responsedata):
    """Correct input time series for instrument response.
        Instr.Resp. is given as 3 column array: frequency, real part, imaginary part

        The given section is demeaned, window tapered, zero padded, bandpassed(with the extreme frequencies of the response frequency axis), deconvolved (by straight division by the array values or interpolated values inbetween), re-transformed, mean-re-added and returned.

    """

    datamean = np.mean(data)
    data -= datamean
    
    freqmin = responsedata[0,0]
    freqmax = responsedata[-1,0]
    print freqmin,freqmax
    
    N = len(data)
    
    window = tukey(N,0.2)#np.blackman(N)

    tapered_data = data * window

    #check, if length is directly equal to power of 2:
    if np.log2(N)%1 == 0:
        next2power = int(np.log2(N))
    else:
        next2power = int(np.log2(N)) + 1

    #zero pad data for significantly faster fft - YES, Numpy does not do that automatically
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


    #return data_freqs,data_spectrum,instr_freqs,instr_spectrum

    #now correct for all frequencies on the data_freqs-axis:
    corrected_spectrum = np.zeros((len(data_spectrum)),'complex')
    for i in range(len(data_spectrum)):
        freq = data_freqs[i]
        spec = data_spectrum[i]

        #this is effectively a boxcar window - maybe to be replaced by proper windowing function ?
        if not (freqmin <= np.abs(freq) <= freqmax):
            print freq, i
            corrected_spectrum[i] = 0#spec
            continue

        factor = 1.

        #find the appropriate frequency most likely inbetween two values on the instr_freqs-axis
        #get the respective value by linear interpolation
        closest_lower = np.abs(freq-instr_freqs).argmin()
        
        #if it coincides with the highest frequency:
        if closest_lower == len(instr_freqs)-1:
            factor = data_spectrum[closest_lower]
        # or the lowest (here: zero)
        elif closest_lower == 0:
            factor = data_spectrum[0]
        #in case the closest frequency value is not lower but higher, take the freq value below as lower bound for the interval:
        else:            
            if instr_freqs[closest_lower] > freq:
                closest_lower -= 1

            instrfreq1 = instr_freqs[closest_lower]
            instrfreq2 = instr_freqs[closest_lower+1]
            instrfactor1 = data_spectrum[closest_lower]
            instrfactor2 = data_spectrum[closest_lower+1]
            interval = instrfreq2 - instrfreq1
            weight = (freq-instrfreq1)/interval

            factor = weight * instrfactor2 + (1-weight) * instrfactor1
        
        corrected_spectrum[i] = spec #/ factor


    #invert into time domain
    correctedTS = np.fft.irfft(corrected_spectrum)

    #cut the zero padding
    correctedTS = correctedTS[:N]

    #re-attach the mean
    correctedTS += datamean

    return correctedTS