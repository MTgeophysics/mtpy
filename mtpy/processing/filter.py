#!/usr/bin/env python

"""
mtpy/processing/filter.py

Functions for the frequency filtering of raw time series. 

@UofA, 2013
(LK)

Revised 2017
JP

"""

#=================================================================
import numpy as np
import  os

import scipy.signal as signal

#=================================================================


def butter_bandpass(lowcut, highcut, samplingrate, order=4):
    nyq = 0.5 * samplingrate
    low = lowcut / nyq
    high = highcut / nyq

    if high >= 1. and low == 0.:
        b = np.array([1.])
        a = np.array([1.])

    elif high < 0.95 and low > 0. :
        wp = [1.05*low,high-0.05]
        ws = [0.95*low,high+0.05]

        order,wn = signal.buttord(wp, ws, 3., 40.)
        b, a = signal.butter(order, wn, btype='band')
    
    elif high >= 0.95:
        print('highpass', low, 1.2*low, 0.8*low)
        order,wn = signal.buttord(15*low, 0.05*low, gpass=0.0, gstop=10.0)
        print(order,wn)
        b, a = signal.butter(order, wn, btype='high')
    
    elif low <= 0.05:
        print('lowpass', high)
        order,wn = signal.buttord(high-0.05, high+0.05, gpass=0.0, gstop=10.0)
        b, a = signal.butter(order, wn, btype='low')

    return b, a

def butter_bandpass_filter(data, lowcut, highcut, samplingrate, order=4):
    b, a = butter_bandpass(lowcut, highcut, samplingrate, order=order)
    y = signal.lfilter(b, a, data)
    return y
    
def low_pass(f, low_pass_freq, cutoff_freq, sampling_rate):
    nyq = .5*sampling_rate
    filt_order, wn = signal.buttord(low_pass_freq/nyq, 
                                    cutoff_freq/nyq, 
                                    3, 40)
                                    
    b, a = signal.butter(filt_order, wn, btype='low')
    f_filt = signal.filtfilt(b, a, f)
    
    return f_filt


def tukey(window_length, alpha=0.2):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width alpha * N / 2  that is convolved with a rectangle window of width (1 - alpha / 2). At alpha = 0 it becomes rectangular, and at alpha = 1 it becomes a Hann window.
 
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
    
def zero_pad(input_array, power=2, pad_fill=0):
    """
    pad the input array with pad_fill to the next power of power.  
    
    For faster fft computation pad the array to the next power of 2 with zeros
    
    Arguments:
    -----------
        **input_array** : np.ndarray (only 1-d arrays are supported at the 
                                      moment)
        
        **power** : [ 2 | 10 ]
                    power look for
        
        **pad_fill** : float or int
                       pad the array with this
                       
    Output:
    --------
        **pad_array** : np.ndarray padded with pad_fill
        
    """

    len_array = input_array.shape[0]
    if power == 2:
        npow = int(np.ceil(np.log2(len_array)))
    if power == 10:
        npow = int(np.ceil(np.log10(len_array)))
    
    if npow > 32:
        print('Exceeding memory allocation inherent in your computer 2**32')
        print('Limiting the zero pad to 2**32')
        
    
    pad_array = np.zeros(power**npow)
    if pad_fill is not 0:
        pad_array[:] = pad_fill
        
    pad_array[0:len_array] = input_array
    
    return pad_array
    
    

def adaptive_notch_filter(bx, df=100, notches=[50, 100], notchradius=.5, 
                          freqrad=.9, rp=.1, dbstop_limit=5.0):
    """
    adaptive_notch_filter(bx, df, notches=[50,100], notchradius=.3, freqrad=.9)
    will apply a notch filter to the array bx by finding the nearest peak 
    around the supplied notch locations.  The filter is a zero-phase 
    Chebyshev type 1 bandstop filter with minimal ripples.
    
    Arguments:
    -----------
        **bx** : np.ndarray(len_time_series)
                 time series to filter
                 
        **df** : float
                 sampling frequency in Hz
                 
        **notches** : list of frequencies (Hz) to filter
                      
        **notchradius** : float
                          radius of the notch in frequency domain (Hz)
        
        **freqrad** : float
                      radius to searching for peak about notch from notches
                      
        **rp** : float
                 ripple of Chebyshev type 1 filter, lower numbers means less
                 ripples
                 
        **dbstop_limit** : float (in decibels)
                           limits the difference between the peak at the 
                           notch and surrounding spectra.  Any difference 
                           above dbstop_limit will be filtered, anything
                           less will not

    Outputs:
    ---------
        
        **bx** : np.ndarray(len_time_series) 
                 filtered array 
                 
        **filtlst** : list
                      location of notches and power difference between peak of
                      notch and average power.
                      
    ..Example: ::
        
        >>> import RemovePeriodicNoise_Kate as rmp
        >>> # make a variable for the file to load in
        >>> fn = r"/home/MT/mt01_20130101_000000.BX"
        >>> # load in file, if the time series is not an ascii file
        >>> # might need to add keywords to np.loadtxt or use another 
        >>> # method to read in the file
        >>> bx = np.loadtxt(fn)
        >>> # create a list of frequencies to filter out
        >>> freq_notches = [50, 150, 200]
        >>> # filter data
        >>> bx_filt, filt_lst = rmp.adaptiveNotchFilter(bx, df=100. 
        >>> ...                                         notches=freq_notches)
        >>> #save the filtered data into a file
        >>> np.savetxt(r"/home/MT/Filtered/mt01_20130101_000000.BX", bx_filt)
    
    Notes:
    -------
        Most of the time the default parameters work well, the only thing
        you need to change is the notches and perhaps the radius.  I would
        test it out with a few time series to find the optimum parameters.
        Then make a loop over all you time series data. Something like
        
        >>> import os
        >>> dirpath = r"/home/MT"
        >>> #make a director to save filtered time series
        >>> save_path = r"/home/MT/Filtered"
        >>> if not os.path.exists(save_path):
        >>>     os.mkdir(save_path)
        >>> for fn in os.listdir(dirpath):
        >>>     bx = np.loadtxt(os.path.join(dirpath, fn)
        >>>     bx_filt, filt_lst = rmp.adaptiveNotchFilter(bx, df=100. 
        >>>     ...                                         notches=freq_notches)
        >>>     np.savetxt(os.path.join(save_path, fn), bx_filt)
         
    """
    
    bx = np.array(bx)
    
    if type(notches) is list:
        notches = np.array(notches)
    elif type(notches) in [float, int]:
        notches = np.array([notches], dtype=np.float)
    
    df = float(df)         #make sure df is a float
    dt = 1./df             #sampling rate


    # transform data into frequency domain to find notches
    BX = np.fft.fft(zero_pad(bx))
    n = len(BX)            #length of array
    dfn = df/n             #frequency step
    dfnn = int(freqrad/dfn)          #radius of frequency search
    fn = notchradius               #filter radius
    freq = np.fft.fftfreq(n,dt)
    
    filtlst = []
    for notch in notches:
        if notch > freq.max():
            break
            #print 'Frequency too high, skipping {0}'.format(notch)
        else:
            fspot = int(round(notch/dfn))
            nspot = np.where(abs(BX)== max(abs(BX[max([fspot-dfnn, 0]):\
                                                 min([fspot+dfnn, n])])))[0][0]   

            med_bx =np.median(abs(BX[max([nspot-dfnn*10, 0]):\
                                     min([nspot+dfnn*10, n])])**2)
            
            #calculate difference between peak and surrounding spectra in dB
            dbstop = 10*np.log10(abs(BX[nspot])**2/med_bx) 
            if np.nan_to_num(dbstop) == 0.0 or dbstop < dbstop_limit:
                filtlst.append('No need to filter \n')
                pass
            else:
                filtlst.append([freq[nspot], dbstop])
                ws = 2*np.array([freq[nspot]-fn, freq[nspot]+fn])/df
                wp = 2*np.array([freq[nspot]-2*fn, freq[nspot]+2*fn])/df
                ford, wn = signal.cheb1ord(wp, ws, 1, dbstop)
                b, a = signal.cheby1(1, .5, wn, btype='bandstop')
                bx = signal.filtfilt(b, a, bx)
    
    return bx, filtlst

def remove_periodic_noise(filename, dt, noiseperiods, save='n'):
    """
    removePeriodicNoise will take a window of length noise period and 
    compute the median of signal for as many windows that can fit within the 
    data.  This median window is convolved with a series of delta functions at
    each window location to create a noise time series. This is then 
    subtracted from the data to get a 'noise free' time series.
    
    Arguments:
    ----------
        **filename** : string (full path to file) or array
                      name of file to have periodic noise removed from
                      can be an array
                      
        **dt** : float
                 time sample rate (s)
                 
        **noiseperiods** : list
                           a list of estimated periods with a range of values
                           to look around [[noiseperiod1,df1]...] where df1 is
                           a fraction value find the peak about noiseperiod1 
                           must be less than 1. (0 is a good start, but
                           if you're periodic noise drifts, might need to
                           adjust df1 to .2 or something)
        **save** : [ 'y' | 'n' ]
                    * 'y' to save file to:
                        os.path.join(os.path.dirname(filename), 'Filtered', fn)
                    * 'n' to return the filtered time series
    
    Outputs:
    --------
    
        **bxnf** : np.ndarray
                   filtered time series
                   
        **pn** : np.ndarray
                periodic noise time series
                
        **fitlst** : list
                     list of peaks found in time series
                     
    ..Example: ::
        
        >>> import RemovePeriodicNoise_Kate as rmp
        >>> # make a variable for the file to load in
        >>> fn = r"/home/MT/mt01_20130101_000000.BX"
        >>> # filter data assuming a 12 second period in noise and save data
        >>> rmp.remove_periodic_noise(fn, 100., [[12,0]], save='y')
    
    Notes:
    -------
        Test out the periodic noise period at first to see if it drifts.  Then
        loop over files
        
        >>> import os
        >>> dirpath = r"/home/MT"
        >>> for fn in os.listdir(dirpath):
        >>>     rmp.remove_periodic_noise(fn, 100., [[12,0]], save='y') 
        
    """
    
    if type(noiseperiods) != list:
        noiseperiods=[noiseperiods]
    
    dt = float(dt)    
    #sampling frequency    
    df = 1./dt
 
    
    filtlst = []
    pnlst = []
    for kk, nperiod in enumerate(noiseperiods):
        #if the nperiod is the first one load file or make an array of the input
        if kk == 0:
            #load file
            if type(filename) is str:
                bx = np.loadtxt(filename)
                m = len(bx)
            else:
                bx = np.array(filename)
                m = len(bx)
        #else copy the already filtered array
        else:
            bx = bxnf.copy()
            m = len(bx)
        
        #get length of array
        T = len(bx)
        
        #frequency step
        dfn = df/T
        #make a frequency array that describes BX
        pfreq = np.fft.fftfreq(int(T), dt)
        
        #get noise period in points along frequency axis
        nperiodnn = round((1./nperiod[0])/dfn)
        
        #get region to look around to find exact peak
        try:
            dfnn = nperiodnn*nperiod[1]
        except IndexError:
            dfnn = .2*nperiodnn
        
        
        #comput FFT of input to find peak value
        BX = np.fft.fft(bx)
        #if dfnn is not 0 then look for max with in region nperiod+-dfnn
        if dfnn != 0:
            nspot = np.where(abs(BX)==
                            max(abs(BX[nperiodnn-dfnn:nperiodnn+dfnn])))[0][0]
        else:
            nspot = nperiodnn
        #output the peak frequency found
        filtlst.append('Found peak at : '+str(pfreq[nspot])+' Hz \n')
        
        #make nperiod the peak period in data points
        nperiod = (1./pfreq[nspot])/dt

        #create list of time instances for windowing
        #nlst=np.arange(start=nperiod,stop=T-nperiod,step=nperiod,dtype='int')
        nlst = np.arange(start=0, stop=m, step=nperiod, dtype='int')
        
        #convolve a series of delta functions with average of periodic window
        dlst = np.zeros(T)               #delta function list
        dlst[0] = 1
        winlst = np.zeros((len(nlst),int(nperiod)))
        for nn,ii in enumerate(nlst):
            if T-ii<nperiod:
                dlst[ii] = 1
            else:
                winlst[nn] = bx[ii:ii+int(nperiod)]
                dlst[ii] = 1
                
        #compute median window to remove any influence of outliers
        medwin = np.median(winlst, axis=0)
        
        #make a time series by convolving
        pn = np.convolve(medwin, dlst)[0:T]

        #remove noise from data
        bxnf = (bx-pn)
        
        pnlst.append(pn)
    if len(pnlst) > 1:
        pn = np.sum(pnlst,axis=0)
    else:
        pn = np.array(pn)
    if save == 'y':
        savepath = os.path.join(os.path.dirname(filename),'Filtered')
        if not os.path.exists(savepath):
            os.mkdir(savepath)
        #savepathCN=os.path.join(savepath,'CN')
        np.savetxt(os.path.join(savepath, filename), bxnf, fmt='%.7g')
        print('Saved filtered file to {0}'.format(os.path.join(savepath, 
                                                               filename)))
    else:
        return bxnf, pn, filtlst