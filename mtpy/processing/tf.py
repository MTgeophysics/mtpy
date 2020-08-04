# -*- coding: utf-8 -*-
"""
mtpy/processing/tf.py

Functions for the time-frequency analysis of time series data.

Output can be visualised with the help of mtpy/imaging/spectrogram.py

JP, 2013

"""

# =================================================================

import numpy as np
import scipy.signal as sps


# =================================================================


def padzeros(f, npad=None, pad_pattern=None):
    """
    padzeros(f) returns a function that is padded with zeros to the next
    power of 2 for faster processing for fft or to length npad if given.

    Arguments:
    ----------
        **f** :  np.ndarray(m, n)
                 array to pad

        **npad** : int
                  length to pad to
                  if None finds next power of 2

        **pad_pattern** : int or float
                          pattern to pad with
                          if None set zeros

    Returns:
    --------
        **fpad** : np.ndarray(m, npad)
                    array f padded to length npad with pad_pattern

    :Example: ::
        >>> x_array = np.sin(np.arange(0, 2, .01)*np.pi/3)
        >>> print len(x_array)
        >>> x_array_pad = padzeros(x_array)
        >>> print len(x_array_pad)

    """

    # make f an array
    f = np.array(f)
    # check dimensions of f
    try:
        n, m = f.shape
    except ValueError:
        n = f.shape[0]
        m = 0
    if npad is None:
        power = np.log2(n)
        fpow = np.floor(power)
        if power != fpow:
            npad = 2 ** (fpow + 1)
        else:
            npad = 2 ** power

    else:
        pass
    if m != 0:
        fpad = np.zeros((npad, m), dtype=type(f[0, 0]))
        fpad[0:n, m - 1] = f[0:n, m - 1]
        if pad_pattern is not None:
            fpad[n:npad, m - 1] = pad_pattern
    else:
        fpad = np.zeros(npad, dtype=type(f[0]))
        fpad[0:n] = f[0:n]
        if pad_pattern is not None:
            fpad[n:npad] = pad_pattern

    return fpad


def sinc_filter(f, fcutoff=10., w=10.0, dt=.001):
    """
    Applies a sinc filter of width w to the function f by multipling in
    the frequency domain.

    Arguments:
    ----------
        **f** : np.ndarray()
                array to filter
        **fcuttoff** : float
                       cutoff frequency of sinc filter
        **w** : float
                length of filter
        **dt** : float
                 sampling rate in time (s)

    Returns:
    ---------
        **f_filt** : np.ndarray()
                     f with sinc filter applied
    """
    # time shift for the peak of the sinc function
    tshift = float(w) / 2.

    # pad the input array with zeros for faster fft
    fpad = padzeros(f)
    Fpad = np.fft.fft(fpad)

    fc = fcutoff

    # make a time series for sinc function
    t = np.arange(start=-tshift, stop=tshift, step=dt)
    # make it the same length as the input array
    filt = np.zeros(len(fpad))

    # calculate the sinc function
    fs = 2 * fc * np.sinc(2 * t * fc)

    # be sure to normalize it so there is no scaling applied to input array
    norm = sum(fs)
    filt[0:len(t)] = fs / norm

    # transform to Fourier domain
    Filt = np.fft.fft(filt)

    # mulitiply the input array and sinc function in the frequency domain
    Filt_F = Fpad * Filt

    # tranform back to time domain and cull edges so the length of the returned
    # array is the same as the input array
    filt_f = np.fft.ifft(Filt_F)
    filt_f = filt_f[len(t) / 2:len(f) + len(t) / 2]

    return filt_f


def dctrend(f):
    """
    dctrend(f) will remove a dc trend from the function f.

    Arguments:
    -----------
        **f** : np.ndarray()
                array to remove dc trend from

    Returns:
    --------
        **fdc** : np.ndarray()
                 array f with dc component removed
    """

    fdc = sps.detrend(f)

    return fdc


def normalize_L2(f):
    """
    normalize_L2(f) returns the array f normalized by the L2 norm ->
    f/(sqrt(sum(abs(x_i)^2))).

    Arguments
        **f** : np.ndarray()
                array to be normalized

    Returns:
    --------
        **fnorm** : np.ndarray()
                    array f normalized in L2 sense
    """

    f = np.array(f)
    fsum = np.sum(np.abs(f))
    if fsum == 0:
        fnorm = f
    else:
        fnorm = f / np.sqrt(np.sum(np.abs(f) ** 2))

    return fnorm


def decimate(f, m, window_function='hanning'):
    """
    resamples the data at the interval m so that the returned array is
    len(f)/m samples long

    Arguments:
    -----------
        **f** : np.ndarray
                array to be decimated

        **m** : int
                decimation factor

        **window_function** : windowing function to apply to the data
                              to make sure the is no Gibbs ringing or aliasing
                              see scipy.signal.window for all the options

    Returns:
    ----------
        **fdec** : np.ndarray()
                    array f decimated by factor m
    """

    n = len(f)
    fdec = sps.resample(f, n / m, window=window_function)

    return fdec


def dwindow(window):
    """
    Calculates the derivative of the given window. Used for reassignment
    methods

    Arguments:
    ----------
        **window** : np.ndarray
                     some sort of windowed array

    Returns:
    --------
        **dwin** : np.ndarray
                   derivative of window
    """

    h = window
    nh = len(h)
    lh = (nh - 1) / 2
    stepheight = (h[0] + h[-1]) / 2.
    ramp = float((h[-1] - h[0])) / nh
    h2 = np.zeros(nh + 2)
    h2[1:nh + 1] = h - stepheight - ramp * \
        np.arange(start=-lh, stop=lh + 1, step=1)

    dwin = (h2[2:nh + 2] - h2[0:nh]) / 2. + ramp
    dwin[0] = dwin[0] + stepheight
    dwin[-1] = dwin[-1] - stepheight

    return dwin


def gausswin(window_len, alpha=2.5):
    """
    gausswin will compute a gaussian window of length winlen with a variance of
    alpha

    Arguments:
    ----------
        **window_len**: int
                        length of desired window
        **alpha** : float
                   1/standard deviation of window,
                   ie full width half max of window

    Returns:
    --------
        **gauss_window** : np.array
                           gaussian window
    """
    lh = (window_len - 1) / 2 + 1 - np.remainder(window_len, 2)
    gt = np.arange(start=-lh, stop=lh + 1, step=1)

    gauss_window = np.exp(-.5 * (alpha * gt / float(lh)) ** 2)

    return gauss_window


def wvd_analytic_signal(fx):
    """
    Computes the analytic signal for WVVD as defined by
    J. M. O' Toole, M. Mesbah, and B. Boashash, (2008),
    "A New Discrete Analytic Signal for Reducing Aliasing in the
     Discrete Wigner-Ville Distribution", IEEE Trans.  on Signal Processing,

    Argument:
    ---------
        **fx** : np.ndarray()
                 signal to compute anlytic signal for with length N

    Returns:
    --------
        **fxa** : np.ndarray()
                  analytic signal of fx with length 2*N
    """

    n = len(fx)

    # pad the time series with zeros
    fxp = padzeros(fx, npad=2 * n)

    # compute the fourier transform
    FX = np.fft.fft(fxp)
    # apply analytic signal
    FX[1:n - 1] = 2 * FX[1:n - 1]
    FX[n:] = 0

    # inverse fourier transform and set anything outside of length n to zero
    fxa = np.fft.ifft(FX)
    fxa[n:] = 0

    return fxa


def stft(fx, nh=2 ** 8, tstep=2 ** 7, ng=1, df=1.0, nfbins=2 ** 10):
    """
    calculate the spectrogam of the given function by calculating the fft of
    a window of length nh at each time instance with an interval of tstep.
    The frequency resolution is nfbins.

    Can compute the cross STFT by inputting fx as [fx1, fx2]

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **nh** : int (should be power of 2)
                 window length for each time step
                 *default* is 2**8 = 256

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **ng** : int (should be odd)
                 length of smoothing window along frequency plane

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated
       """

    # get length of input time series if there is two columns
    if isinstance(fx, list):
        fx = np.array(fx)
    try:
        fn, fm = fx.shape
        if fm < fn:
            fm, fn = fx.shape
    except ValueError:
        fn = fx.shape[0]
        fm = 1
    if fm > 1:
        fx = fx.reshape(fn)
    else:
        fx = fx.reshape(fn)

    # make a hanning window to minimize aliazing and Gibbs effect of short time
    # windows
    h = normalize_L2(np.hanning(nh))

    # make a hanning window to smooth in frequency domain
    if ng != 1:
        if np.remainder(ng, 2) != 1:
            ng = ng - 1
            print('ng forced to be odd as ng-1')
        else:
            pass
        g = normalize_L2(np.hanning(ng))
    else:
        pass

    # make time step list
    tlst = np.arange(start=0, stop=fn - nh + 1, step=tstep)

    df = float(df)

    # get only positive frequencies
    flst = np.fft.fftfreq(nfbins, 1 / df)[0:int(nfbins / 2)]

    # initialize the TFD array
    tfarray = np.zeros((int(nfbins / 2), len(tlst)), dtype='complex128')

    # calculate the analytic signal to fold negative frequencies onto the
    # positive ones
    fa = sps.hilbert(dctrend(fx))

    # compute the fft at each window instance
    for place, ii in enumerate(tlst):
        fxwin = fa[ii:ii + nh] * h

        # get only positive frequencies
        FXwin = np.fft.fft(padzeros(fxwin, npad=nfbins))[:int(nfbins / 2)]

        # smooth in frequency plane
        if ng != 1:
            FXwin = np.convolve(padzeros(FXwin, npad=len(FXwin) + ng - 1), g,
                                'valid')
        else:
            pass

        # pull out only positive quadrant, flip array for plotting
        tfarray[:, place] = FXwin[::-1]

    return tfarray, tlst, flst


def reassigned_stft(fx, nh=2 ** 6 - 1, tstep=2 ** 5, nfbins=2 ** 10, df=1.0, alpha=4,
                    threshold=None):
    """
    Computes the reassigned spectrogram by estimating the center of gravity of
    the signal and condensing dispersed energy back to that location.  Works
    well for data with minimal noise and strong spectral structure.

    Arguments:
    ----------
        **fx** : np.ndarray
                 time series to be analyzed

        **nh** : int(should be odd)
                 length of gaussian window that is applied to the short
                 time intervals
                 *default* is 127

        **tstep** : int
                    time step for each window calculation
                    *default* is 64

        **nfbins** : int (should be a power of 2 and larger or equal to nh
                     number of frequency bins to calculate, note result will be
                     length nfbins/2
                     *default* is 1024

        **df** : float or int
                 sampling frequency (Hz)

        **alpha** : float
                   reciprocal of full width half max of gaussian window
                   *default* is 4

        **threshold** : float
                        threshold value for reassignment
                        If None the threshold is automatically calculated
                        *default* is None

    Returns:
        **rtfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                       reassigned spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated

        **stft** : np.ndarray(nfbins/2, len(fx)/tstep)
                   standard spectrogram calculated from stft
                   in units of amplitude
    """

    # make sure fx is type array
    fx = np.array(fx)

    # compute length of fx
    nx = len(fx)

    # make sure window length is odd
    if np.remainder(nh, 2) == 0:
        nh = nh + 1

    # compute gaussian window
    h = gausswin(nh, alpha=alpha)
    lh = (nh - 1) / 2

    # compute ramp window
    th = h * np.arange(start=-lh, stop=lh + 1, step=1)

    # compute derivative of window
    dh = dwindow(h)

    # make a time list of indexes
    tlst = np.arange(start=0, stop=nx, step=tstep)
    nt = len(tlst)

    # make a frequency list
    flst = np.fft.fftfreq(nfbins, 1. / df)[nfbins / 2:]
    return_flst = np.fft.fftfreq(nfbins, 1. / df)[0:nfbins / 2]

    # initialize some time-frequency arrays
    tfr = np.zeros((nfbins, nt), dtype='complex')
    tf2 = np.zeros((nfbins, nt), dtype='complex')
    tf3 = np.zeros((nfbins, nt), dtype='complex')

    # compute components for reassignment
    for ii, tt in enumerate(tlst):
        # create a time shift list
        tau = np.arange(start=-min([np.round(nx / 2.), lh, tt - 1]),
                        stop=min([np.round(nx / 2.), lh, nx - tt - 1]) + 1)
        # compute the frequency spots to be calculated
        ff = np.remainder(nfbins + tau, nfbins)
        xlst = tt + tau
        hlst = lh + tau
        normh = np.sqrt(np.sum(abs(h[hlst]) ** 2))
        tfr[ff, ii] = fx[xlst] * h[hlst].conj() / normh
        tf2[ff, ii] = fx[xlst] * th[hlst].conj() / normh
        tf3[ff, ii] = fx[xlst] * dh[hlst].conj() / normh

    # compute Fourier Transform
    spec = np.fft.fft(tfr, axis=0)
    spect = np.fft.fft(tf2, axis=0)
    specd = np.fft.fft(tf3, axis=0)

    # get only positive frequencies
    spec = spec[nfbins / 2:, :]
    spect = spect[nfbins / 2:, :]
    specd = specd[nfbins / 2:, :]

    # check to make sure no spurious zeros floating around
    spec[np.where(abs(spec) < 1.E-6)] = 0.0
    zerofind = np.nonzero(abs(spec))
    twspec = np.zeros((nfbins / 2, nt), dtype='float')
    dwspec = np.zeros((nfbins / 2, nt), dtype='float')
    twspec[zerofind] = np.round(np.real(spect[zerofind] / spec[zerofind]) / 1)
    dwspec[zerofind] = np.round(np.imag((nfbins / 2.) *
                                        specd[zerofind] / spec[zerofind]) / (np.pi))

    # compute reassignment
    rtfarray = np.zeros_like(spec)

    if threshold is None:
        threshold = 1.E-4 * np.mean(fx[tlst])

    for nn in range(nt):
        for kk in range(nfbins / 2):
            if abs(spec[kk, nn]) > threshold:
                # get center of gravity index in time direction
                nhat = int(nn + twspec[kk, nn])
                nhat = int(min([max([nhat, 1]), nt - 1]))
                # get center of gravity index in frequency direction
                khat = int(kk - dwspec[kk, nn])
                khat = int(np.remainder(np.remainder(khat - 1, nfbins / 2) + nfbins / 2,
                                        nfbins / 2))
                # reassign energy
                rtfarray[khat, nhat] = rtfarray[khat, nhat] + spec[kk, nn]
                spect[kk, nn] = khat + 1j * nhat
            else:
                spect[kk, nn] = np.inf * (1 + 1j)
                rtfarray[kk, nn] = rtfarray[kk, nn] + spec[kk, nn]

    return rtfarray, tlst, return_flst, spec


def wvd(fx, nh=2 ** 8 - 1, tstep=2 ** 5, nfbins=2 ** 10, df=1.0):
    """
    calculates the Wigner-Ville distribution of f.

    Can compute the cross spectra by inputting fx as [fx1,fx2]

    Arguments:
    ----------
                **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **nh** : int (should be odd)
                 window length for each time step
                 *default* is 2**8-1 = 255

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated

    """
    # check to see if computing auto or cross-spectra
    if isinstance(fx, list):
        fx = np.array(fx)
    try:
        fn, fm = fx.shape
        if fm > fn:
            fm, fn = fx.shape
    except ValueError:
        fn = len(fx)
        fm = 1

    if fm > 1:
        fn = fn[0]
        print('computing cross spectra')
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx[0])
        fb = wvd_analytic_signal(fx[1])
    else:
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx)
        fa = sps.hilbert(dctrend(fx))
        fb = fa.copy()

    fn = len(fa)
    # sampling period
    df = float(df)
    dt = 1. / df
    # time shift
    tau = (nh - 1) / 2

    # create a time array such that the first point is centered on time window
    tlst = np.arange(start=0, stop=fn - 1, step=tstep, dtype='int')

    # create an empty array to put the tf in
    tfarray = np.zeros((nfbins, len(tlst)), dtype='complex')

    # create a frequency array with just positive frequencies
    flst = np.fft.fftfreq(nfbins, dt)[0:nfbins / 2]

    # calculate pseudo WV
    for point, nn in enumerate(tlst):
        # calculate the smallest timeshift possible
        tau_min = min(nn, tau, fn - nn - 1)
        # make a timeshift array
        tau_lst = np.arange(start=-tau_min, stop=tau_min + 1, step=1,
                            dtype='int')
        # calculate rectangular windowed correlation function of analytic
        # signal
        Rnn = 4 * np.conjugate(fa[nn - tau_lst]) * fb[nn + tau_lst]
        # calculate fft of windowed correlation function
        # put into tfarray
        tfarray[:, point] = padzeros(Rnn, npad=nfbins)[::-1]

    # compute Fourier Transform of array along the time axis
    tfarray = np.fft.fft(tfarray, axis=0)
    # normalize
    tfarray = tfarray / nh

    return tfarray, tlst, flst


def spwvd(fx, tstep=2 ** 5, nfbins=2 ** 10, df=1.0, nh=None, ng=None, sigmat=None,
          sigmaf=None):
    """
    Calculates the smoothed pseudo Wigner-Ville distribution for an array
    fx. Smoothed with Gaussians windows to get best localization.

    Can be input as [fx1, fx2] to compute cross spectra.

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **nh** : int (should be odd)
                 window length for each time step
                 *default* is None and window is calculated automatically

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **ng** : int (should be odd)
                 length of smoothing window along frequency plane

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

        **sigmat** : float
                    std of window h, ie full width half max of gaussian
                    *default* is None and sigmat is calculate automatically

        **sigmaf** : float
                     std of window g, ie full width half max of gaussian
                     *default* is None and sigmaf is calculate automatically

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      SPWVD spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated
    """
    # check to see if calculating the auto or cross spectra
    if isinstance(fx, list):
        fx = np.array(fx)

    try:
        fn, fm = fx.shape
        if fm > fn:
            fm, fn = fx.shape
    except ValueError:
        fn = len(fx)
        fm = 1
    if fm > 1:
        print('computing cross spectra')
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx[0])
        fb = wvd_analytic_signal(fx[1])

    else:
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx)
        fa = sps.hilbert(dctrend(fx))
        fb = fa.copy()
        print('Computed Analytic signal')

    # sampling period
    df = float(df)
    dt = 1 / df

    # create normalize windows in time (g) and frequency (h)
    # note window length should be odd so that h,g[0]=1,nh>ng
    if nh is None:
        nh = np.floor(fn / 2.)
    # make sure the window length is odd
    if np.remainder(nh, 2) == 0:
        nh += 1
    # calculate length for time smoothing window
    if ng is None:
        ng = np.floor(fn / 5.)
    if np.remainder(ng, 2) == 0:
        ng += 1
    # calculate standard deviations for gaussian windows
    if sigmat is None:
        sigmah = nh / (6 * np.sqrt(2 * np.log(2)))
    else:
        sigmah = sigmat

    if sigmaf is None:
        sigmag = ng / (6 * np.sqrt(2 * np.log(2)))
    else:
        sigmag = sigmaf
    nh = int(nh)
    ng = int(ng)
    print('nh=' + str(nh) + '; ng=' + str(ng))
    # calculate windows and normalize
    h = sps.gaussian(nh, sigmah)
    h /= sum(h)

    g = sps.gaussian(ng, sigmag)
    g /= sum(g)

    Lh = (nh - 1) / 2  # midpoint index of window h
    Lg = (ng - 1) / 2  # midpoint index of window g

    # create a time array such that the first point is centered on time window
    tlst = np.arange(start=0, stop=fn + 1, step=tstep, dtype='int')

    # create an empty array to put the tf in
    # make sure data type is complex
    tfarray = np.zeros((nfbins, len(tlst)), dtype='complex')

    # create a frequency array with just positive frequencies
    flst = np.fft.fftfreq(nfbins, dt)[0:nfbins / 2]

    # calculate pseudo WV
    for point, t in enumerate(tlst):
        # find the smallest possible time shift
        tau_max = min(t + Lg - 1, fn - t + Lg, round(nfbins / 2), Lh)
        # create time lag list
        taulst = np.arange(start=-min(Lg, fn - t), stop=min(Lg, t - 1) + 1, step=1,
                           dtype='int')
        # calculate windowed correlation function of analytic function for
        # zero frequency
        tfarray[0, point] = sum(2 * (g[Lg + taulst] / sum(g[Lg + taulst])) *
                                fa[t - taulst - 1] * np.conjugate(fb[t - taulst - 1]))
        # calculate tfd by calculating convolution of window and correlation
        # function as sum of correlation function over the lag period times the
        # window at that point. Calculate symmetrical segments for FFT later
        for mm in range(tau_max):
            taulst = np.arange(start=-min(Lg, fn - t - mm - 1), stop=min(Lg, t - mm - 1) + 1,
                               step=1, dtype='int')
            # compute positive half
            gm = 2 * (g[Lg + taulst] / sum(g[Lg + taulst]))
            Rmm = sum(gm * fa[t + mm - taulst - 1] *
                      np.conjugate(fb[t - mm - taulst]))
            tfarray[mm, point] = h[Lh + mm - 1] * Rmm
            # compute negative half
            Rmm = sum(gm * fa[t - mm - taulst] *
                      np.conjugate(fb[t + mm - taulst - 1]))
            tfarray[nfbins - mm - 1, point] = h[Lh - mm] * Rmm
        mm = round(nfbins / 2)

        if t <= fn - mm and t >= mm and mm <= Lh:
            print('doing weird thing')
            taulst = np.arange(start=-min(Lg, fn - t - mm), stop=min(Lg, fn - t, mm) + 1,
                               step=1, dtype='int')
            gm = g[Lg + taulst] / sum(g[Lg + taulst])
            tfarray[mm - 1, point] = .5 * \
                (sum(h[Lh + mm] * (gm * fa[t + mm - taulst - 1] *
                                   np.conjugate(fb[t - mm - taulst]))) +
                 sum(h[Lh - mm] * (gm * fa[t - mm - taulst] *
                                   np.conjugate(fb[t + mm - taulst - 1]))))

    tfarray = np.fft.fft(tfarray, axis=0)
    # rotate for plotting purposes so that (t=0,f=0) is at the lower left
    tfarray = np.rot90(tfarray.T, 1)

    return tfarray, tlst, flst


def robust_wvd(fx, nh=2 ** 7 - 1, ng=2 ** 4 - 1, tstep=2 ** 4, nfbins=2 ** 8, df=1.0,
               sigmat=None, sigmaf=None):
    """
    Calculate the robust Wigner-Ville distribution for an array
    fx. Smoothed with Gaussians windows to get best localization.

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **nh** : int (should be power of 2)
                 window length for each time step
                 *default* is 2**8 = 256

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **ng** : int (should be odd)
                 length of smoothing window along frequency plane
                 *default* is 2**4-1 = 15

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

        **sigmat** : float
                    std of window h, ie full width half max of gaussian
                    *default* is None and sigmat is calculate automatically

        **sigmaf** : float
                     std of window g, ie full width half max of gaussian
                     *default* is None and sigmaf is calculate automatically

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated
    """
    # check to see if computing the auto or cross spectra
    if isinstance(fx, list):
        fx = np.array(fx)

    try:
        fn, fm = fx.shape
        if fm > fn:
            fm, fn = fx.shape
    except ValueError:
        fn = len(fx)
        fm = 1
    if fm > 1:
        print('computing cross spectra')
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx[0])
        fb = wvd_analytic_signal(fx[1])

    else:
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx)
        fa = sps.hilbert(dctrend(fx))
        fb = fa.copy()
        print('Computed Analytic signal')

        # make sure window length is odd
    if nh is None:
        nh = np.floor(fn / 2.)
    # make sure the window length is odd
    if np.remainder(nh, 2) == 0:
        nh += 1
    # calculate length for time smoothing window
    if ng is None:
        ng = np.floor(fn / 5.)
    if np.remainder(ng, 2) == 0:
        ng += 1
    nh = int(nh)
    ng = int(ng)
    print('nh = {0}'.format(nh))
    print('ng = {0}'.format(ng))

    dt = 1. / (df * 2.)

    # get length of input time series
    nfx = len(fa)

    # make frequency smoothing window
    if sigmat is None:
        sigmat = nh / (5 * np.sqrt(2 * np.log(2)))
    h = sps.gaussian(nh, sigmat)
    h /= sum(h)

    # make a time smoothing window
    if sigmaf is None:
        sigmaf = ng / (5 * np.sqrt(2 * np.log(2)))
    g = sps.gaussian(ng, sigmaf)
    g /= sum(g)

    mlst = np.arange(start=-nh / 2 + 1, stop=nh / 2 + 1, step=1, dtype='int')
    tlst = np.arange(start=nh / 2, stop=nfx - nh / 2, step=tstep)
    # make a frequency list for plotting exporting only positive frequencies
    # get only positive frequencies
    flst = np.fft.fftfreq(nfbins, dt)[nfbins / 2:]
    flst[-1] = 0
    flstp = np.fft.fftfreq(nfbins, 2 * dt)[0:nfbins / 2]

    # create an empty array to put the tf in
    tfarray = np.zeros((nfbins / 2, len(tlst)), dtype='complex')

    for tpoint, nn in enumerate(tlst):
        # calculate windowed correlation function of analytic function
        fxwin = h * fa[nn + mlst] * fb[nn - mlst].conj()
        for fpoint, mm in enumerate(flst):
            fxmed = np.convolve(g, fxwin * np.exp(1j * 4 * np.pi * mlst * mm * dt),
                                mode='same') / (nh * ng)
            fxmedpoint = np.median(fxmed.real)
            if fxmedpoint == 0.0:
                tfarray[fpoint, tpoint] = 1E-10
            else:
                tfarray[fpoint, tpoint] = fxmedpoint

    tfarray = (4. * nh / dt) * tfarray

    return tfarray, tlst, flstp


def specwv(fx, tstep=2 ** 5, nfbins=2 ** 10, nhs=2 ** 8, nhwv=2 ** 9 - 1, ngwv=2 ** 3 - 1,
           df=1.0):
    """
    Calculates the Wigner-Ville distribution mulitplied by the STFT windowed
    by the common gaussian window h for an array f.  Handy for removing cross
    terms in the wvd.

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for

       **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **nhs** : int (should be power of 2)
                 window length for each time step to caclulate STFT
                 *default* is 2**8 = 256 and window is calculated automatically


        **nhwv** : int (should be odd)
                 length of smoothing window for each time step to calculate
                 WVD. *default* is 2**9-1 = 511

        **ngwv** : int (should be odd)
                 length of frequency smoothing window for each time step to
                 calculate WVD. *default* is 2**3-1 = 7

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins


    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      SPWVD spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated
    """

    # calculate stft
    pst, tlst, flst = stft(fx, nh=nhs, tstep=tstep, nfbins=nfbins, df=df)

    # calculate new time step so WVD and STFT will align
    ntstep = len(fx) / (len(tlst) * 2.)

    # calculate spwvd
    pwv, twv, fwv = spwvd(fx, tstep=ntstep, nfbins=nfbins, df=df,
                          nh=nhwv, ng=ngwv)

    # multiply the two together normalize
    tfarray = pst / pst.max() * pwv / pwv.max()

    return tfarray, tlst, flst


def modifiedb(fx, tstep=2 ** 5, nfbins=2 ** 10,
              df=1.0, nh=2 ** 8 - 1, beta=.2):
    """
    Calculates the modified b distribution as defined by cosh(n)^-2 beta
    for an array fx.  Supposed to remove cross terms in the WVD.

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **nh** : int (should be odd)
                 window length for each time step
                 *default* is None and window is calculated automatically

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**5 = 32

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins
        **beta** : float
                   smoothing coefficient ussully between [0, 1]

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      SPWVD spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated
    """
    # check to see if calculating the auto or cross spectra
    if isinstance(fx, list):
        fx = np.array(fx)
    try:
        fn, fm = fx.shape
        if fm > fn:
            fm, fn = fx.shape
    except ValueError:
        fn = len(fx)
        fm = 1

    if fm > 1:
        fn = fn[0]
        print('computing cross spectra')
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx[0])
        fb = wvd_analytic_signal(fx[1])
    else:
        # compute the analytic signal of function f and dctrend
        fa = wvd_analytic_signal(fx)
        fa = sps.hilbert(dctrend(fx))
        fb = fa.copy()

    # sampling period
    df = float(df)
    dt = 1. / df

    tau = (nh - 1) / 2  # midpoint index of window h

    # create a time array such that the first point is centered on time window
    tlst = np.arange(start=0, stop=fn - 1, step=tstep, dtype='int')

    # create an empty array to put the tf in
    tfarray = np.zeros((nfbins, len(tlst)), dtype='complex')

    # create a frequency array with just positive frequencies
    flst = np.fft.fftfreq(nfbins, dt)[0:nfbins / 2]

    # calculate pseudo WV
    for point, nn in enumerate(tlst):
        # calculate the smallest timeshift possible
        tau_min = min(nn, tau, fn - nn - 1)
        # make a timeshift array
        taulst = np.arange(start=-tau_min, stop=tau_min +
                           1, step=1, dtype='int')
        # create modified b window
        mbwin = np.cosh(taulst) ** (-2 * beta)
        mbwin = mbwin / sum(mbwin)
        MBwin = np.fft.fft(padzeros(mbwin, npad=nfbins))
        # calculate windowed correlation function of analytic function
        Rnn = np.conjugate(fa[nn - taulst]) * fb[nn + taulst]
        # calculate fft of windowed correlation function
        FTRnn = MBwin * np.fft.fft(padzeros(Rnn, npad=nfbins))
        # put into tfarray
        tfarray[:, point] = FTRnn[::-1]

    # need to cut the time frequency array in half due to the WVD assuming
    # time series sampled at twice nyquist.
    # tfarray = tfarray

    return tfarray, tlst, flst


def robust_stft_median(fx, nh=2 ** 8, tstep=2 ** 5, df=1.0, nfbins=2 ** 10):
    """
    Calculates the robust spectrogram using the vector median simplification.

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **nh** : int (should be power of 2)
                 window length for each time step
                 *default* is 2**8 = 256

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **df** : float
                 sampling frequency (Hz)

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated
    """

    # get length of input time series
    nfx = len(fx)

    # compute time shift list
    mlst = np.arange(start=-nh / 2 + 1, stop=nh / 2 + 1, step=1, dtype='int')
    # compute time locations to take STFT
    tlst = np.arange(start=0, stop=nfx - nh + 1, step=tstep)

    # make a frequency list for plotting exporting only positive frequencies
    flst = np.fft.fftfreq(nfbins, 1 / df)
    flstc = flst[nfbins / 2:]
    # Note: these are actually the negative frequencies but works better for
    # calculations
    flstp = flst[0:nfbins / 2]

    # make time window and normalize
    sigmanh = nh / (6 * np.sqrt(2 * np.log(2)))
    h = sps.gaussian(nh, sigmanh)
    h = h / sum(h)

    # create an empty array to put the tf in and initialize a complex value
    tfarray = np.zeros((nfbins / 2, len(tlst)), dtype='complex')

    # take the hilbert transform of the signal to make complex and remove
    # negative frequencies
    fa = sps.hilbert(dctrend(fx))
    fa = fa / fa.std()

    # make a frequency list for plotting exporting only positive frequencies
    # get only positive frequencies
    flst = np.fft.fftfreq(nfbins, 1 / df)[nfbins / 2:]

    for tpoint, nn in enumerate(tlst):
        # calculate windowed correlation function of analytic function
        fxwin = h * fa[nn:nn + nh]
        for fpoint, mm in enumerate(flstc):
            fxmed = fxwin * np.exp(1j * 2 * np.pi * mlst * mm / df)
            fxmedreal = np.median(fxmed.real)
            fxmedimag = np.median(fxmed.imag)
            if fxmedreal + 1j * fxmedimag == 0.0:
                tfarray[fpoint, tpoint] = 1E-10
            else:
                tfarray[fpoint, tpoint] = fxmedreal + 1j * fxmedimag
    # normalize tfarray
    tfarray = (4. * nh * df) * tfarray

    return tfarray, tlst, flstp


def robust_stft_L(fx, alpha=.325, nh=2 ** 8, tstep=2 **
                  5, df=1.0, nfbins=2 ** 10):
    """
    Calculates the robust spectrogram by estimating the vector median and
    summing terms estimated by alpha coefficients.

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **alpha** : float
                    robust parameter [0,.5] -> 0 gives spectrogram,
                    0.5 gives median stft
                    *default* is 0.325

        **nh** : int (should be power of 2)
                 window length for each time step
                 *default* is 2**8 = 256

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated

    """

    # get length of input time series
    nfx = len(fx)

    # compute time shift list
    mlst = np.arange(start=-nh / 2 + 1, stop=nh / 2 + 1, step=1, dtype='int')
    # compute time locations to take STFT
    tlst = np.arange(start=0, stop=nfx - nh + 1, step=tstep)

    # make a frequency list for plotting exporting only positive frequencies
    flst = np.fft.fftfreq(nfbins, 1 / df)
    flstc = flst[nfbins / 2:]
    # Note: these are actually the negative frequencies but works better for
    # calculations
    flstp = flst[0:nfbins / 2]

    # make time window and normalize
    sigmanh = nh / (6 * np.sqrt(2 * np.log(2)))
    h = sps.gaussian(nh, sigmanh)
    h /= sum(h)

    # create an empty array to put the tf in and initialize a complex value
    tfarray = np.zeros((nfbins / 2, len(tlst)), dtype='complex')

    # take the hilbert transform of the signal to make complex and remove
    # negative frequencies
    fa = sps.hilbert(dctrend(fx))
    fa /= fa.std()

    # make a frequency list for plotting exporting only positive frequencies
    # get only positive frequencies
    flst = np.fft.fftfreq(nfbins, 1 / df)[nfbins / 2:]

    # create list of coefficients
    a = np.zeros(nh)
    a[(nh - 2) * alpha:alpha * (2 - nh) + nh - 1] = 1. / \
        (nh * (1 - 2 * alpha) + 4 * alpha)

    for tpoint, nn in enumerate(tlst):
        # calculate windowed correlation function of analytic function
        fxwin = h * fa[nn:nn + nh]
        for fpoint, mm in enumerate(flstc):
            fxelement = fxwin * np.exp(1j * 2 * np.pi * mlst * mm / df)
            fxreal = np.sort(fxelement.real)[::-1]
            fximag = np.sort(fxelement.imag)[::-1]
            tfpoint = sum(a * (fxreal + 1j * fximag))
            if tfpoint == 0.0:
                tfarray[fpoint, tpoint] = 1E-10
            else:
                tfarray[fpoint, tpoint] = tfpoint
    # normalize tfarray
    tfarray = (4. * nh * df) * tfarray

    return tfarray, tlst, flstp


def smethod(fx, L=11, nh=2 ** 8, tstep=2 ** 7, ng=1, df=1.0, nfbins=2 ** 10,
            sigmaL=None):
    """
    Calculates the smethod by estimating the STFT first and computing the WV
    of window length L in the frequency domain.

    For larger L more of WV estimation, if L=0 get back STFT

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **L** : int (should be odd)
                length of window for S-method calculation, higher numbers tend
                toward WVD

        **nh** : int (should be power of 2)
                 window length for each time step
                 *default* is 2**8 = 256

        **ng** : int (should be odd)
                 length of smoothing window along frequency plane

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      S-method spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated

        **pxx** : np.ndarray(nfbins/2, len(fx)/tstep)
                  STFT spectrogram in units of amplitude

    """

    df = float(df)

    # check to see if computing auto or cross spectra
    if isinstance(fx, list):
        fx = np.array(fx)
    try:
        fn, fm = fx.shape
        if fm > fn:
            fm, fn = fx.shape
    except ValueError:
        fn = len(fx)
        fm = 1
    if fm > 1:
        print('computing cross spectra')
        # compute the analytic signal of function f and dctrend
        # fa=sps.hilbert(dctrend(fx[0]))
        # fb=sps.hilbert(dctrend(fx[1]))
        fa = fx[0]
        fb = fx[1]
        fa = fa.reshape(fn)
        fb = fb.reshape(fn)
        pxa, tlst, flst = stft(fa, nh=nh, tstep=tstep, ng=ng, df=df,
                               nfbins=nfbins)
        pxb, tlst, flst = stft(fb, nh=nh, tstep=tstep, ng=ng, df=df,
                               nfbins=nfbins)
        pxx = pxa * pxb.conj()
    else:
        # compute the analytic signal of function f and dctrend
        # fa=sps.hilbert(dctrend(fx))
        fa = fx
        fa = fa.reshape(fn)
        fb = fa
        pxx, tlst, flst = stft(fa, nh=nh, tstep=tstep, ng=ng, df=df,
                               nfbins=nfbins)

    # make an new array to put the new tfd in
    tfarray = abs(pxx) ** 2
    # get shape of spectrogram
    nf, nt = tfarray.shape
    # create a list of frequency shifts
    Llst = np.arange(start=-L / 2 + 1, stop=L / 2 + 1, step=1, dtype='int')
    # create a frequency gaussian window
    if sigmaL is None:
        sigmaL = L / (1 * np.sqrt(2 * np.log(2)))
    p = sps.gaussian(L, sigmaL)
    # make a matrix of windows
    pm = np.zeros((L, nt))
    for kk in range(nt):
        pm[:, kk] = p

    # loop over frequency and calculate the s-method
    for ff in range(int(L / 2), nf - int(L / 2) - 1):
        tfarray[ff, :] = tfarray[ff, :] + 2 * np.real(np.sum(pm * pxx[ff + Llst, :] *
                                                             pxx[ff - Llst, :].conj(), axis=0))
    # normalize
    tfarray[int(L / 2):int(-L / 2)] /= L

    return tfarray, tlst, flst, pxx


def robust_smethod(fx, L=5, nh=2 ** 7, tstep=2 ** 5, nfbins=2 ** 10, df=1.0,
                   robusttype='median', sigmaL=None, alpha=.325):
    """
    Computes the robust Smethod via the robust spectrogram.

    Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **L** : int (should be odd)
                length of window for S-method calculation, higher numbers tend
                toward WVD

        **nh** : int (should be power of 2)
                 window length for each time step
                 *default* is 2**8 = 256

        **ng** : int (should be odd)
                 length of smoothing window along frequency plane

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

        **robusttype** : [ 'median' | 'L' ]
                         type of robust STFT to compute. *default* is 'median'

        **simgaL** : float
                    full-width half max of gaussian window applied in frequency

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      S-method spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated

        **pxx** : np.ndarray(nfbins/2, len(fx)/tstep)
                  STFT spectrogram in units of amplitude
    """
    # check to see if computing auto or cross spectra
    if isinstance(fx, list):
        fx = np.array(fx)
    try:
        fn, fm = fx.shape
        if fm > fn:
            fm, fn = fx.shape
    except ValueError:
        fn = len(fx)
        fm = 1
    if fm > 1:
        print('computing cross spectra')
        # compute the analytic signal of function f and dctrend
        fa = fx[0].reshape(fn)
        fb = fx[1].reshape(fn)
        if robusttype == 'median':
            pxa, tlst, flst = robust_stft_median(fa, nh=nh, tstep=tstep, df=df,
                                                 nfbins=nfbins)
            pxb, tlst, flst = robust_stft_median(fb, nh=nh, tstep=tstep, df=df,
                                                 nfbins=nfbins)
        elif robusttype == 'L':
            pxa, tlst, flst = robust_stft_L(fa, nh=nh, tstep=tstep, df=df,
                                            nfbins=nfbins, alpha=alpha)
            pxb, tlst, flst = robust_stft_L(fb, nh=nh, tstep=tstep, df=df,
                                            nfbins=nfbins, alpha=alpha)
        else:
            raise NameError('robusttype {0} undefined'.format(robusttype))
        pxx = pxa * pxb.conj()
    else:
        fa = fx.reshape(fn)
        if robusttype == 'median':
            pxx, tlst, flst = robust_stft_median(fa, nh=nh, tstep=tstep, df=df,
                                                 nfbins=nfbins)
        elif robusttype == 'L':
            pxx, tlst, flst = robust_stft_L(fa, nh=nh, tstep=tstep, df=df,
                                            nfbins=nfbins, alpha=alpha)
        else:
            raise NameError('robusttype {0} undefined'.format(robusttype))

    # compute frequency shift list
    Llst = np.arange(start=-L / 2 + 1, stop=L / 2 + 1, step=1, dtype='int')

    # compute the frequency window of length L
    if sigmaL is None:
        sigmaL = L / 3 * (np.sqrt(2 * np.log(2)))
    lwin = gausswin(L, sigmaL)
    lwin /= sum(lwin)
    pm = np.zeros((L, len(tlst)))
    for kk in range(len(tlst)):
        pm[:, kk] = lwin

    smarray = pxx.copy()
    # compute S-method
    for ff in range(L / 2, nfbins / 2 - L / 2):
        smarray[ff, :] = smarray[ff, :] + 2 * np.real(np.sum(pm * pxx[ff + Llst, :] *
                                                             pxx[ff - Llst,
                                                                 :].conj(),
                                                             axis=0))
    # normalize
    smarray = (2. / (L * nh)) * smarray

    return smarray, tlst, flst, pxx


def reassigned_smethod(fx, nh=2 ** 7 - 1, tstep=2 ** 4, nfbins=2 ** 9, df=1.0, alpha=4,
                       thresh=.01, L=5):
    """
    Calulates the reassigned S-method as described by Djurovic[1999] by
    using the spectrogram to estimate the reassignment.

        Arguments:
    -----------
        **fx** : list or np.ndarray
                 the function to have a spectrogram computed for
                 for cross-correlation input as [fx1, fx2]

        **L** : int (should be odd)
                length of window for S-method calculation, higher numbers tend
                toward WVD

        **nh** : int (should be power of 2)
                 window length for each time step
                 *default* is 2**8 = 256

        **alpha** : float
                    inverse of full-width half max of gaussian window, smaller
                    numbers mean broader windows.

        **thresh** : float
                     threshold for reassignment, lower numbers more points
                     reassigned, higer numbers less points reassigned

        **tstep** : int
                    number of sample between short windows
                    *default* is 2**7 = 128

        **df** : float
                 sampling frequency

        **nfbins** : int (should be power of 2 and equal or larger than nh)
                     number of frequency bins

    Returns:
    --------
        **tfarray** : np.ndarray(nfbins/2, len(fx)/tstep)
                      S-method spectrogram in units of amplitude

        **tlst** : np.array()
                   array of time instances for each window calculated

        **flst** : np.ndarray(nfbins/2)
                   frequency array containing only positive frequencies where
                   the Fourier coeffients were calculated

        **sm** : np.ndarray(nfbins/2, len(fx)/tstep)
                  S-method spectrogram in units of amplitude

    """
    # check to see if computing auto or cross spectra
    if isinstance(fx, list):
        fx = np.array(fx)
    try:
        fn, fm = fx.shape
        if fm > fn:
            fm, fn = fx.shape
    except ValueError:
        fn = len(fx)
        fm = 1
    if fm > 1:
        print('computing cross spectra')
        fa = fx[0]
        fb = fx[1]
        fa = fa.reshape(fn)
        fb = fb.reshape(fn)
    else:
        fa = fx
        fa = fa.reshape(fn)
        fb = fa.copy()

    nx = len(fx)

    # compute gaussian window
    h = gausswin(nh, alpha=alpha)
    lh = (nh - 1) / 2

    # compute ramp window
    th = h * np.arange(start=-lh, stop=lh + 1, step=1)

    # compute derivative of window
    dh = dwindow(h)

    # make a time list of indexes
    tlst = np.arange(start=0, stop=nx, step=tstep)
    nt = len(tlst)

    # make frequency list for plotting
    flst = np.fft.fftfreq(nfbins, 1. / df)[:nfbins / 2]

    # initialize some time-frequency arrays
    tfh = np.zeros((nfbins, nt), dtype='complex')
    tfth = np.zeros((nfbins, nt), dtype='complex')
    tfdh = np.zeros((nfbins, nt), dtype='complex')

    # compute components for reassignment
    for ii, tt in enumerate(tlst):
        # create a time shift list
        tau = np.arange(start=-min([np.round(nx / 2.), lh, tt - 1]),
                        stop=min([np.round(nx / 2.), lh, nx - tt - 1]) + 1)
        # compute the frequency spots to be calculated
        ff = np.remainder(nfbins + tau, nfbins)
        # make lists of data points for each window calculation
        xlst = tt + tau
        hlst = lh + tau
        normh = np.sqrt(np.sum(abs(h[hlst]) ** 2))
        tfh[ff, ii] = fx[xlst] * h[hlst].conj() / normh
        tfth[ff, ii] = fx[xlst] * th[hlst].conj() / normh
        tfdh[ff, ii] = fx[xlst] * dh[hlst].conj() / normh

    # compute Fourier Transform
    spech = np.fft.fft(tfh, axis=0)
    specth = np.fft.fft(tfth, axis=0)
    specdh = np.fft.fft(tfdh, axis=0)

    # get only positive frequencies
    spech = spech[nfbins / 2:, :]
    specth = specth[nfbins / 2:, :]
    specdh = specdh[nfbins / 2:, :]

    # check to make sure no spurious zeros floating around
    szf = np.where(abs(spech) < 1.E-6)
    spech[szf] = 0.0 + 0.0j
    zerofind = np.nonzero(abs(spech))
    twspec = np.zeros((nfbins / 2, nt), dtype='float')
    dwspec = np.zeros((nfbins / 2, nt), dtype='float')
    twspec[zerofind] = np.round(np.real(specth[zerofind] / spech[zerofind]))
    dwspec[zerofind] = np.round(np.imag((nfbins / 2.) * specdh[zerofind] /
                                        spech[zerofind]) / (np.pi))

    # get shape of spectrogram
    nf, nt = spech.shape

    # -----calculate s-method-----
    Llst = np.arange(start=-L / 2 + 1, stop=L / 2 + 1, step=1, dtype='int')

    # make and empty array of zeros
    sm = np.zeros_like(spech)

    # put values where L cannot be value of L, near top and bottom
    sm[0:L / 2, :] = abs(spech[0:L / 2, :]) ** 2
    sm[-L / 2:, :] = abs(spech[-L / 2:, :]) ** 2

    # calculate s-method
    for ff in range(L / 2, nf - L / 2 - 1):
        sm[ff, :] = 2 * np.real(np.sum(spech[ff + Llst, :] * spech[ff - Llst, :].conj(),
                                       axis=0)) / L

    # ------compute reassignment-----
    rtfarray = np.zeros((nfbins / 2, nt))

    threshold = thresh * np.max(abs(sm))

    for nn in range(nt):
        for kk in range(nf):
            if abs(spech[kk, nn]) > threshold:
                # get center of gravity index in time direction from
                # spectrogram
                nhat = int(nn + twspec[kk, nn])
                nhat = int(min([max([nhat, 1]), nt - 1]))
                # get center of gravity index in frequency direction from spec
                khat = int(kk - dwspec[kk, nn])
                khat = int(np.remainder(np.remainder(khat - 1, nfbins / 2) + nfbins / 2,
                                        nfbins / 2))
                rtfarray[khat, nhat] = rtfarray[khat, nhat] + abs(sm[kk, nn])
            else:
                rtfarray[kk, nn] = rtfarray[kk, nn] + sm[kk, nn]

    # place values where L cannot be L
    rtfarray[:L / 2, :] = abs(sm[:L / 2, :])
    rtfarray[-L / 2:, :] = abs(sm[-L / 2:, :])

    # for plotting purposes set anything that is 0 to 1, so that log(1) will
    # plot as zero
    rtfarray[np.where(rtfarray == 0)] = 1.0

    sm[np.where(sm == 0.0)] = 1.0

    # scale
    rtfarray = abs(rtfarray)

    return rtfarray, tlst, flst, sm


def stfbss(X, nsources=5, ng=2 ** 5 - 1, nh=2 ** 9 - 1, tstep=2 ** 6 - 1, df=1.0, nfbins=2 ** 10,
           tftol=1.E-8, L=7, normalize=True, tftype='spwvd', alpha=.38):
    """
    btfssX,nsources=5,ng=2**5-1,nh=2**9-1,tstep=2**6-1,df=1.0,nfbins=2**10,
          tftol=1.E-8,normalize=True)
    estimates sources using a blind source algorithm based on spatial
    time-frequency distributions.  At the moment this algorithm uses the SPWVD
    to estimate TF distributions.

    Arguments
        X = m x n array of time series, where m is number of time series and n
           is length of each time series
        nsources = number of estimated sources
        ng = frequency window length
        nh = time window length
        tstep = time step increment
        df = sampling frequency (Hz)
        nfbins = number of frequencies
        tftol = tolerance for a time-frequency point to be estimated as a cross
                term or as an auto term, the higher the number the more auto
                terms.
        normalization = True or False, True to normalize, False if already
                        normalized

    Returns:

        Se = estimated individual signals up to a permutation and scale
        Ae = estimated mixing matrix as X=A*S
    """

    # get shape of timeseries matrix,
    # m=number of channels
    # tlen=length of timeseries
    m, maxn = X.shape

    n = nsources
    # get number of time bins
    ntbins = int(float(maxn) / tstep)

    if tftype == 'spwvd':
        tfkwargs = {'ng': ng, 'nh': nh, 'df': df,
                    'nfbins': nfbins, 'tstep': tstep}

    elif tftype == 'smethod':
        tfkwargs = {'ng': ng, 'nh': nh, 'df': df,
                    'nfbins': nfbins, 'tstep': tstep, 'L': L}

    elif tftype == 'Lestimate':
        tfkwargs = {'nh': nh, 'df': df, 'nfbins': nfbins, 'tstep': tstep, 'L': L,
                    'alpha': alpha, 'robusttype': 'L'}
    # remove dc component from time series and normalize
    if normalize == True:
        for ii in range(m):
            X[ii, :] = X[ii, :] - np.mean(X[ii, :])
            X[ii, :] = X[ii, :] / X[ii, :].std()

            # ===============================================================================
    # Whiten data and Compute Whitening matrix
    # ===============================================================================
    # whiten data to get a unitary matrix with unit variance and zero mean
    # compute covariance matrix
    Rxx = (np.dot(X, X.T)) / float(maxn)

    # calculate eigen decomposition
    [l, u] = np.linalg.eig(Rxx)

    # sort eigenvalues from smallest to largest assuming largest are sources and
    # smallest are noise
    lspot = l.argsort()
    eigval = l[lspot]
    eigvec = u[:, lspot]

    # calculate the noise variance as mean of non-principal components
    sigman = np.mean(eigval[0:m - n])

    # compute scaling factor for whitening matrix
    wscale = 1 / np.sqrt(eigval[m - n:m] - sigman)

    # compute whitening matrix
    W = np.zeros((m, n))
    for kk in range(n):
        W[:, kk] = wscale[kk] * eigvec[:, m - n + kk].T
    W = W.T
    # compute whitened signal vector. Note the dimensionality is reduced from [mxn]
    # to [nxn] making the computation simpler.
    Z = np.dot(W, X)

    # ===============================================================================
    #  Compute Spatial Time Frequency Distribution
    # ===============================================================================

    Za = np.array(Z.copy())

    if tftype == 'spwvd':
        stfd = np.zeros((n, n, nfbins, ntbins + 1), dtype='complex128')
        # compute auto terms
        for ii in range(n):
            pswvd, tswvd, fswvd = spwvd(Za[ii].reshape(maxn), **tfkwargs)
            stfd[ii, ii, :, :] = pswvd

        # compute cross terms
        for jj in range(n):
            for kk in range(jj, n):
                pswvd, tswvd, fswvd = spwvd([Za[jj].reshape(maxn),
                                             Za[kk].reshape(maxn)],
                                            **tfkwargs)
                stfd[jj, kk, :, :] = pswvd
                stfd[kk, jj, :, :] = pswvd.conj()

    elif tftype == 'smethod':
        nfbins = nfbins / 2
        stfd = np.zeros((n, n, nfbins, ntbins), dtype='complex128')
        # compute auto terms
        for ii in range(n):
            psm, tsm, fsm, pst = smethod(Za[ii].reshape(maxn), **tfkwargs)
            stfd[ii, ii, :, :] = psm

        # compute cross terms
        for jj in range(n):
            for kk in range(jj, n):
                psm, tsm, fsm, pst = smethod([Za[jj].reshape(maxn),
                                              Za[kk].reshape(maxn)],
                                             **tfkwargs)
                stfd[jj, kk, :, :] = psm
                stfd[kk, jj, :, :] = psm.conj()

    elif tftype == 'Lestimate':
        nfbins = nfbins / 2
        stfd = np.zeros((n, n, nfbins, ntbins), dtype='complex128')
        # compute auto terms
        for ii in range(n):
            psm, tsm, fsm, pst = robust_smethod(
                Za[ii].reshape(maxn), **tfkwargs)
            stfd[ii, ii, :, :psm.shape[1]] = psm

        # compute cross terms
        for jj in range(n):
            for kk in range(jj, n):
                psm, tsm, fsm, pst = robust_smethod([Za[jj].reshape(maxn),
                                                     Za[kk].reshape(maxn)],
                                                    **tfkwargs)
                stfd[jj, kk, :, :psm.shape[1]] = psm
                stfd[kk, jj, :, :psm.shape[1]] = psm.conj()

    # ===============================================================================
    # Compute criteria for cross terms
    # ===============================================================================

    stfdTr = np.zeros((nfbins, ntbins))
    C = np.zeros((nfbins, ntbins))

    for ff in range(nfbins):
        for tt in range(ntbins):
            # compensate for noise
            stfd[:, :, ff, tt] = stfd[:, :, ff, tt] - \
                sigman * np.matrix(W) * np.matrix(W.T)
            # compute the trace
            stfdTr[ff, tt] = abs(np.trace(stfd[:, :, ff, tt]))
    # compute mean over entire t-f plane
    trmean = stfdTr.mean()

    # find t-f points that meet the criteria
    fspot, tspot = np.nonzero(stfdTr > trmean)

    for ll in range(len(fspot)):
        treig = abs(np.linalg.eig(stfd[:, :, fspot[ll], tspot[ll]])[0])
        if sum(treig) != 0 and sum(treig) > tftol:
            C[fspot[ll], tspot[ll]] = max(treig) / sum(treig)
        else:
            C[fspot[ll], tspot[ll]] = 0

    # compute gradients and jacobi matrices
    negjacobi = np.zeros((nfbins, ntbins))
    smallgrad = np.zeros((nfbins, ntbins))
    maxpoints = np.zeros((nfbins, ntbins))

    gradt, gradf = np.gradient(C)
    Jtt, Jtf = np.gradient(gradt)
    Jft, Jff = np.gradient(gradf)

    # get points when L2 of gradient is smaller than tolerance level
    smallgrad = np.where(np.sqrt(gradt ** 2 + gradf ** 2) < tftol, 1, 0)

    # get points where the Jacobi is negative definite
    detjacobi = Jtt * Jff - Jtf * Jft
    negjacobi = np.where(detjacobi > 0, 1, 0) * np.where(Jtt < 0, 1, 0) \
        * np.where((Jtt + Jff) < 0, 1, 0)

    maxpoints = smallgrad * negjacobi

    gfspot, gtspot = np.nonzero(maxpoints)
    ntfpoints = len(gfspot)

    if ntfpoints == 0:
        raise ValueError('Found no tf points, relax tolerance')
    else:
        print('Found ' + str(ntfpoints) + ' t-f points')

    for rr in range(ntfpoints):
        if rr == 0:
            Rjd = stfd[:, :, gfspot[rr], gtspot[rr]]
        else:
            Rjd = np.concatenate(
                (Rjd, stfd[:, :, gfspot[rr], gtspot[rr]]), axis=1)
    Rjd = np.array(Rjd)

    # ===============================================================================
    # Calculate Joint Diagonalization
    # ===============================================================================
    # get size of array of matrices to be diagonalized
    mtf, nm = Rjd.shape  # mtf is number of t-f points, nm is number of matrices
    # set up some initial parameters
    V = np.eye(mtf)

    # update boolean
    encore = True
    # Total number of rotations
    updates = 0
    sweep = 0
    # print 'Computing Joint Diagonalization'
    # Joint diagonalization proper
    # ============================
    while encore:
        # reset some parameters
        encore = False
        sweep += 1
        upds = 0
        Vkeep = V
        for p in range(mtf):

            for q in range(p + 1, mtf):
                # set up indices
                qi = np.arange(start=q, stop=nm, step=mtf)
                pi = np.arange(start=p, stop=nm, step=mtf)
                # computation of Givens angle
                g = np.array([Rjd[p, pi] - Rjd[q, qi], Rjd[p, qi], Rjd[q, pi]])
                gg = np.real(np.dot(g, g.T))
                ton = gg[0, 0] - gg[1, 1]
                toff = gg[0, 1] + gg[1, 0]
                theta = 0.5 * np.arctan2(toff, ton +
                                         np.sqrt(ton ** 2 + toff ** 2))
                # Givens update
                if abs(theta) > tftol:
                    encore = True
                    upds += 1
                    c = np.cos(theta)
                    s = np.sin(theta)
                    G = np.matrix([[c, -s], [s, c]])
                    pair = np.array([p, q])
                    V[:, pair] = V[:, pair] * G
                    Rjd[pair, :] = G.T * Rjd[pair, :]
                    Rjd[:, np.concatenate([pi, qi])] = np.append(
                        c * Rjd[:, pi] + s * Rjd[:, qi],
                        -s * Rjd[:, pi] + c * Rjd[:, qi], axis=1)
        updates += upds
    print('Updated ' + str(updates) + ' times.')

    # compute estimated signal matrix
    Se = np.dot(V.T, Z)
    # compute estimated mixing matrix
    Ae = np.dot(np.linalg.pinv(W), V)

    return Se, Ae
