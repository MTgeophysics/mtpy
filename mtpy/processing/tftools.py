# -*- coding: utf-8 -*-
"""
Created on Mon May 03 14:53:54 2010

@author: a1185872
"""

import numpy as np
import scipy.signal as sps
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def padzeros(f,npad=None,padpattern=None):
    """
    padzeros(f) will return a function that is padded with zeros to the next
    power of 2 for faster processing for fft or to length npad if given.
    
    Inputs:
        f = array to pad 
        npad = length to pad to defaults to next power of two
        padpattern = pattern to pad with default is zero
    
    Outputs:
        fpad = array f padded to length npad with padpattern
    """
    
    #make f an array
    f=np.array(f)
    #check dimensions of f
    try:
        n,m=f.shape
    except ValueError:
        n=f.shape[0]
        m=0
    if npad==None:
        power=np.log2(n)
        fpow=np.floor(power)
        if power!=fpow:
            npad=2**(fpow+1)
        else:
            npad=2**power

    else:
        pass
    if m!=0:
        fpad=np.zeros((npad,m),dtype=type(f[0,0]))
        fpad[0:n,m-1]=f[0:n,m-1]
        if padpattern!=None:
            fpad[n:npad,m-1]=padpattern
    else:
        fpad=np.zeros(npad,dtype=type(f[0]))
        fpad[0:n]=f[0:n]
        if padpattern!=None:
            fpad[n:npad]=padpattern
        
    return fpad
    
def sfilter(f,fcutoff=10.,w=10.0,dt=.001):
    """
    Will apply a sinc filter of width w to the function f by multipling in
    the frequency domain. Returns filtered function
    
    Inputs:
        f = array to filter
        fcuttoff = cutoff frequency
        w = length of filter
        dt = sampling time (s)
    
    Outputs:
        filtfunc = filtered function
    """

    tshift=float(w)/2.
    
    fpad=padzeros(f)
    Fpad=np.fft.fft(fpad)
    fc=fcutoff
    
    t=np.arange(start=-tshift,stop=tshift,step=dt)
    filt=np.zeros(len(fpad))
    fs=2*fc*np.sinc(2*t*fc)
    norm=sum(fs)
    filt[0:len(t)]=fs/norm
    Filt=np.fft.fft(filt)
    
    Filtfunc=Fpad*Filt
    filtfunc=np.fft.ifft(Filtfunc)
    filtfunc=filtfunc[len(t)/2:len(f)+len(t)/2]
    
    return filtfunc

def dctrend(f):
    """
    dctrend(f) will remove a dc trend from the function f.
    
    Inputs:
        f = array to dctrend
        
    Outputs:
        fdc = array f with dc component removed
    """
    
    fdc=sps.detrend(f)
    
    return fdc

def normalizeL2(f):
    """
    normalizeL2(f) will return the function f normalized by the L2 norm ->
    f/(sqrt(sum(abs(x_i)^2))).
    
    Inputs:
        f = array to be normalized
    
    Outputs:
        fnorm = array f normalized in L2 sense
    """
    
    f=np.array(f)
    fsum=np.sum(np.abs(f))
    if fsum==0:
        fnorm=f
    else:
        fnorm=f/np.sqrt(np.sum(np.abs(f)**2))
    
    return fnorm

def decimatef(f,m):
    """
    Will decimate a function by the factor m. First an 8th order Cheybechev 
    type I filter with a cuttoff frequency of .8/m  is applied in both 
    directions to minimize any phase distortion and remove any aliasing. Note 
    decimation values above 10 will typically result in bad coefficients, 
    therefore if you decimation is more than 10 just repeat the decimation until
    the desired decimation is reached.
    
    Inputs:
        f = array to be decimated
        m = decimation factor
    
    Outputs:
        fdec = array f decimated by factor m
    """
    
    n=len(f)
    fdec=sps.resample(f,n/m,window='hanning')
#    n=len(f)
#    nout=np.ceil(n/m)
#    nfilt=8
#    rip=.05
#    
#    #make a cheybeshev1 zero-phase filter with cuttoff frequency of .8/m
#    b,a=sps.iirfilter(nfilt,.8/m,rp=rip,btype='low',ftype='cheby1',output='ba')
#    ffilt=sps.filtfilt(b,a,f)
#    nbeg=n-m*nout
#    fdec=np.array([ffilt[ii] for ii in np.arange(start=nbeg,stop=int(n),step=m)])
    
    return fdec

def dwindow(window):
    """
    Calculates the derivative of the given window
    
    Input:
        window = some sort of window function
    
    Output:
        dwin = derivative of window
    """
    
    h=window
    nh=len(h)
    lh=(nh-1)/2
    stepheight=(h[0]+h[-1])/2.
    ramp=float((h[-1]-h[0]))/nh
    h2=np.zeros(nh+2)
    h2[1:nh+1]=h-stepheight-ramp*np.arange(start=-lh,stop=lh+1,step=1)
    
    dwin=(h2[2:nh+2]-h2[0:nh])/2.+ramp
    dwin[0]=dwin[0]+stepheight
    dwin[-1]=dwin[-1]-stepheight
    
    return dwin

def gausswin(winlen,alpha=2.5):
    """
    gausswin will compute a gaussian window of length winlen with a variance of
    alpha
    
    Inputs:
        winlen = length of desired window
        alpha = 1/standard deviation of window, ie full width half max of window 
    
    Outputs:
        gwin = gaussian window
    """
    lh=(winlen-1)/2+1-np.remainder(winlen,2)
    gt=np.arange(start=-lh,stop=lh+1,step=1)
    
    gwin=np.exp(-.5*(alpha*gt/float(lh))**2)
    
    return gwin
    
def wvdas(fx):
    """
    wvdas(fx) will compute the analytic signal for WVVD as defined by \
    J. M. O' Toole, M. Mesbah, and B. Boashash, (2008), "A New Discrete Analytic\
    Signal for Reducing Aliasing in the Discrete Wigner-Ville Distribution", \
    IEEE Trans.  on Signal Processing,
    
    Inputs:
        fx = signal to compute anlytic signal for with length N
    
    Outputs:
        fxa = analytic signal of fx with length 2*N
    """
    
    n=len(fx)
    
    #pad the time series with zeros
    fxp=padzeros(fx,npad=2*n)
    
    #compute the fourier transform
    FX=np.fft.fft(fxp)
    #apply analytic signal    
    FX[1:n-1]=2*FX[1:n-1]
    FX[n:]=0
    
    #inverse fourier transform and set anything outside of length n to zero 
    fxa=np.fft.ifft(FX)
    fxa[n:]=0
    
    return fxa
    
    


def stft(fx,nh=2**8,tstep=2**7,ng=1,df=1.0,nfbins=2**10):
    """stft(fx,nh=2**8,tstep=2**7,ng=1,df=1.0) will calculate the spectrogam of
    the given function by calculating the fft of a window of length nh at each
    time instance with an interval of tstep.  The frequency resolution is nfbins
    Can compute the cross STFT by inputting fx as [fx1,fx2]
    
    Inputs:
        fx = the function to have a spectrogram computed for can be two functions
            input as [fx1,fx2] 
        nh = window length for each time step 
        tstep = time step between short windows 
        ng = smoothing window along frequency plane should be odd
        df = sampling frequency 
        nfbins = number of frequency bins
    
    Outputs:
        tfarray = spectrogram in units of amplitude 
        tlst = time instance array where each window was calculated
        flst = frequency array containing only positive frequencies
       """
    
    #get length of input time series if there is two columns
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm<fn:
            fm,fn=fx.shape
    except ValueError:
        fn=fx.shape[0]
        fm=1
    if fm>1:
        fx=fx.reshape(fn)
    else:
        fx=fx.reshape(fn)
    #make a hanning window to minimize aliazing and Gibbs effect of short time 
    #windows
    h=normalizeL2(np.hanning(nh))
    #make a hanning window to smooth in frequency domain
    if ng!=1:
        if np.remainder(ng,2)!=1:
            ng=ng-1
            print 'ng forced to be odd as ng-1'
        else:
            pass
        g=normalizeL2(np.hanning(ng))
    else:
        pass
    #make time step list
    tlst=np.arange(start=0,stop=fn-nh+1,step=tstep)
    #make a frequency list for plotting exporting only positive frequencies
    df=float(df)
    flst=np.fft.fftfreq(nfbins,1/df)[0:nfbins/2] #get only positive frequencies
    #initialize the TFD array
    tfarray=np.zeros((nfbins/2,len(tlst)),dtype='complex128')
    
    fa=sps.hilbert(dctrend(fx))
    
    for place,ii in enumerate(tlst):
        fxwin=fa[ii:ii+nh]*h
        #get only positive frequencies
        FXwin=np.fft.fft(padzeros(fxwin,npad=nfbins))[:nfbins/2]
        #smooth in frequency plane
        if ng!=1:
            FXwin=np.convolve(padzeros(FXwin,npad=len(FXwin)+ng-1),g,'valid')
        else:
            pass
        #pull out only positive quadrant, flip array for plotting
        tfarray[:,place]=FXwin[::-1]
        
    return tfarray,tlst,flst
    
def reassignedstft(fx,nh=2**6-1,tstep=2**5,nfbins=2**10,df=1.0,alpha=4,
                   threshold=None):
    """
    reassignedstft(fx,nh=2**5-1,tstep=2**8,nfbins=2**10,df=1.0,alpha=20) will 
    compute the reassigned spectrogram by estimating the center of gravity of 
    the signal and condensing dispersed energy back to that location.
   
    Inputs:
        fx = time series to be analyzed
        nh = length of gaussian window, should be odd
        tstep = time step for each window calculation
        nfbins = number of frequency bins to calculate, note result will be 
                 length nfbins/2
        df = sampling frequency (Hz)
        alpha = reciprocal of full width half max of gaussian window
        threshold = threshold value for reassignment
        
    Outputs:
        rtfarray = reassigned spectrogram in units of amplitude
        tlst = array of time instances where windows were calculated for ploting
        flst = array of frequencies for plotting
        stft = standard spectrogram in units of amplitude
    """
    #make sure fx is type array
    fx=np.array(fx)
    #compute length of fx
    nx=len(fx)
    
    #make sure window length is odd
    if np.remainder(nh,2)==0:
        nh=nh+1
    
    #compute gaussian window
    h=gausswin(nh,alpha=alpha)
    #h=np.hanning(nh)
    lh=(nh-1)/2
    
    #compute ramp window
    th=h*np.arange(start=-lh,stop=lh+1,step=1)
    
    #compute derivative of window
    dh=dwindow(h)
    
    #make a time list of indexes
    tlst=np.arange(start=0,stop=nx,step=tstep)
    nt=len(tlst)
    
    #make a frequency list
    flst=np.fft.fftfreq(nfbins,1./df)[nfbins/2:]
    
    #initialize some time-frequency arrays
    tfr=np.zeros((nfbins,nt),dtype='complex128')
    tf2=np.zeros((nfbins,nt),dtype='complex128')
    tf3=np.zeros((nfbins,nt),dtype='complex128')
    
    #compute components for reassignment
    for ii,tt in enumerate(tlst):
    	#create a time shift list
        tau=np.arange(start=-min([np.round(nx/2.),lh,tt-1]),
                   stop=min([np.round(nx/2.),lh,nx-tt-1])+1)
    	#compute the frequency spots to be calculated
        ff=np.remainder(nfbins+tau,nfbins)
        xlst=tt+tau
        hlst=lh+tau
        normh=np.sqrt(np.sum(abs(h[hlst])**2))
        tfr[ff,ii]=fx[xlst]*h[hlst].conj()/normh
        tf2[ff,ii]=fx[xlst]*th[hlst].conj()/normh
        tf3[ff,ii]=fx[xlst]*dh[hlst].conj()/normh
    
    #compute Fourier Transform
    spec=np.fft.fft(tfr,axis=0)
    spect=np.fft.fft(tf2,axis=0)
    specd=np.fft.fft(tf3,axis=0)
    
    #get only positive frequencies
    spec=spec[nfbins/2:,:]
    spect=spect[nfbins/2:,:]
    specd=specd[nfbins/2:,:]
    
    #check to make sure no spurious zeros floating around
    szf=np.where(abs(spec)<1.E-6)
    spec[szf]=0.0
    zerofind=np.nonzero(abs(spec))
    twspec=np.zeros((nfbins/2,nt),dtype='float')
    dwspec=np.zeros((nfbins/2,nt),dtype='float')
    twspec[zerofind]=np.round(np.real(spect[zerofind]/spec[zerofind])/1)
    dwspec[zerofind]=np.round(np.imag((nfbins/2.)*specd[zerofind]/spec[zerofind])/
                                    (np.pi))
    
    #compute reassignment
    rtfarray=np.zeros_like(spec)
    
    if threshold==None:
        threshold=1.E-4*np.mean(fx[tlst])

    for nn in range(nt):
        for kk in range(nfbins/2):
            if abs(spec[kk,nn])>threshold:
                #get center of gravity index in time direction
                nhat=int(nn+twspec[kk,nn])
                nhat=int(min([max([nhat,1]),nt-1]))
                #get center of gravity index in frequency direction
                khat=int(kk-dwspec[kk,nn])
                khat=int(np.remainder(np.remainder(khat-1,nfbins/2)+nfbins/2,
                                      nfbins/2))
                #reassign energy
                rtfarray[khat,nhat]=rtfarray[khat,nhat]+spec[kk,nn]
                #rtfarray[kk,nn]=spec[khat,nhat]
                spect[kk,nn]=khat+1j*nhat
            else:
                spect[kk,nn]=np.inf*(1+1j)
                rtfarray[kk,nn]=rtfarray[kk,nn]+spec[kk,nn]
        
    return rtfarray,tlst,flst,spec
    
    
def wvd(fx,nh=2**8-1,tstep=2**5,nfbins=2**10,df=1.0):
    """
    wvd(f,nh=2**8-1,tstep=2**5,nfbins=2**10,df=1.0) will calculate the 
    Wigner-Ville distribution for a function f. Can compute the cross spectra
    by inputting fx as [fx1,fx2] 
    
    Inputs:
        fx = array for which WVD will be calculated, input as [fx1,fx2] for 
            cross-spectra calculation
        nh = window length, needs to be odd so centered on zero
        tstep = time step between windows
        nfbins = number of frequencies
        df = sampling frequency (Hz)
    
    Outputs:
        tfarray = WVD estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm>fn:
            fm,fn=fx.shape
    except ValueError:
        fn=len(fx)
        fm=1
    if fm>1:
        fn=fn[0]
        print 'computing cross spectra'
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx[0])
        fb=wvdas(fx[1])

    else:
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx)        
        fa=sps.hilbert(dctrend(fx))
        fb=fa.copy()
    fn=len(fa)
    #sampling period
    df=float(df)
    dt=1./df
    tau=(nh-1)/2
    
    #create a time array such that the first point is centered on time window
    tlst=np.arange(start=0,stop=fn-1,step=tstep,dtype='int')
      
    #create an empty array to put the tf in 
    tfarray=np.zeros((nfbins,len(tlst)),dtype='complex128')
    
    #create a frequency array with just positive frequencies
    flst=np.fft.fftfreq(nfbins,dt)[0:nfbins/2]
    
    #calculate pseudo WV
    for point,nn in enumerate(tlst):
        #calculate the smallest timeshift possible
        taun=min(nn,tau,fn-nn-1)
        #make a timeshift array
        taulst=np.arange(start=-taun,stop=taun+1,step=1,dtype='int')
        #calculate rectangular windowed correlation function of analytic signal
        Rnn=4*np.conjugate(fa[nn-taulst])*fb[nn+taulst]  
        #calculate fft of windowed correlation function
        #FTRnn=np.fft.fft(padzeros(Rnn,npad=nfbins))
        #put into tfarray
        tfarray[:,point]=padzeros(Rnn,npad=nfbins)[::-1]
        
    #normalize
    tfarray=np.fft.fft(tfarray,axis=0)
    tfarray=tfarray/nh
    
    return tfarray,tlst,flst

def spwvd(fx,tstep=2**5,nfbins=2**10,df=1.0,nh=None,ng=None,sigmat=None,
          sigmaf=None):
    """
    spwvd(fx,tstep=2**5,nfbins=2**10,df=1.0,nh=2**8-1,ng=2**5-1,sigmat=None,
          sigmaf=None)
    will calculate the smoothed pseudo Wigner-Ville distribution for a function
    fx. smoothed with Gaussians windows to get best localization. 
    
    Inputs:
        fx = array to estimate spwvd, input as [fx1,fx2] if computing cross
            spectra
        tstep = time step between windows
        nfbins = number of frequencies 
        df = sampling frequency (Hz)
        ng = length of time-domain smoothing window (needs to be odd)
        nh = length of frequency-domain smoothing window (needs to be odd) 
        sigmat = std of window h, ie full width half max of gaussian 
        sigmaf = std of window g, ie full width half max of gaussian
    
    Outputs:
        tfarray = SPWVD estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm>fn:
            fm,fn=fx.shape
    except ValueError:
        fn=len(fx)
        fm=1
    if fm>1:
        print 'computing cross spectra'
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx[0])
        fb=wvdas(fx[1])

    else:
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx)        
        fa=sps.hilbert(dctrend(fx))
        fb=fa.copy()
        print 'Computed Analytic signal'
    
    #sampling period
    df=float(df)
    dt=1/df
    
    #create normalize windows in time (g) and frequency (h)
    #note window length should be odd so that h,g[0]=1,nh>ng
    if nh==None:
        nh=np.floor(fn/2.)
    #make sure the window length is odd
    if np.remainder(nh,2)==0:
        nh=nh+1
    #calculate length for time smoothing window
    if ng==None:
        ng=np.floor(fn/5.)
    if np.remainder(ng,2)==0:
        ng=ng+1
    #calculate standard deviations for gaussian windows    
    if sigmat==None:
        sigmah=nh/(6*np.sqrt(2*np.log(2)))
    else:
        sigmah=sigmat
    
    if sigmaf==None:
        sigmag=ng/(6*np.sqrt(2*np.log(2)))
    else:
        sigmag=sigmaf
    nh=int(nh)
    ng=int(ng)
    print 'nh='+str(nh)+'; ng='+str(ng)
    #calculate windows and normalize
    h=sps.gaussian(nh,sigmah)
    h=h/sum(h)
    
    g=sps.gaussian(ng,sigmag)
    g=g/sum(g)
    
    Lh=(nh-1)/2  #midpoint index of window h
    Lg=(ng-1)/2   #midpoint index of window g
    
    #create a time array such that the first point is centered on time window
    tlst=np.arange(start=0,stop=fn+1,step=tstep,dtype='int')
    
    #create an empty array to put the tf in 
    #make sure data type is complex 
    tfarray=np.zeros((nfbins,len(tlst)),dtype='complex128')
    
    #create a frequency array with just positive frequencies
    flst=np.fft.fftfreq(nfbins,dt)[0:nfbins/2]
    
    #calculate pseudo WV
    for point,t in enumerate(tlst):
        #find the smallest possible time shift
        maxtau=min(t+Lg-1,fn-t+Lg,round(nfbins/2),Lh)
        #create time lag list
        taulst=np.arange(start=-min(Lg,fn-t),stop=min(Lg,t-1)+1,step=1,
                         dtype='int')
        #calculate windowed correlation function of analytic function for
        #zero frequency 
        tfarray[0,point]=sum(2*(g[Lg+taulst]/sum(g[Lg+taulst]))*fa[t-taulst-1]*
                            np.conjugate(fb[t-taulst-1]))
        #calculate tfd by calculating convolution of window and correlation 
        #function as sum of correlation function over the lag period times the
        #window at that point. Calculate symmetrical segments for FFT later
        for mm in range(maxtau):
            taulst=np.arange(start=-min(Lg,fn-t-mm-1),stop=min(Lg,t-mm-1)+1,
                             step=1,dtype='int')
            #compute positive half
            gm=2*(g[Lg+taulst]/sum(g[Lg+taulst]))
            Rmm=sum(gm*fa[t+mm-taulst-1]*np.conjugate(fb[t-mm-taulst]))
            tfarray[mm,point]=h[Lh+mm-1]*Rmm
            #compute negative half 
            Rmm=sum(gm*fa[t-mm-taulst]*np.conjugate(fb[t+mm-taulst-1]))
            tfarray[nfbins-mm-1,point]=h[Lh-mm]*Rmm
        mm=round(nfbins/2)
        
        if t<=fn-mm and t>=mm and mm<=Lh:
            print 'doing weird thing'
            taulst=np.arange(start=-min(Lg,fn-t-mm),stop=min(Lg,fn-t,mm)+1,step=1,
                             dtype='int')
            gm=g[Lg+taulst]/sum(g[Lg+taulst])
            tfarray[mm-1,point]=.5*\
                (sum(h[Lh+mm]*(gm*fa[t+mm-taulst-1]*
                np.conjugate(fb[t-mm-taulst])))+\
                sum(h[Lh-mm]*(gm*fa[t-mm-taulst]*
                np.conjugate(fb[t+mm-taulst-1]))))
    
    tfarray=np.fft.fft(tfarray,axis=0)
    #rotate for plotting purposes so that (t=0,f=0) is at the lower left
    tfarray=np.rot90(tfarray.T,1)
    
    return tfarray,tlst,flst
    
def robustwvd(fx,nh=2**7-1,ng=2**4-1,tstep=2**4,nfbins=2**8,df=1.0,
              sigmanh=None,sigmang=None):
    """
    robustwvd(fx,tstep=2**5,nfbins=2**10,df=1.0,nh=2**8-1,ng=2**5-1,
              sigmanh=None,sigmang=None)
    will calculate the smoothed pseudo Wigner-Ville distribution for a function
    fx. smoothed with Gaussians windows to get best localization. 
    
    Inputs:
        fx = array to estimate spwvd, input as [fx1,fx2] if computing cross
            spectra
        tstep = time step between windows
        nfbins = number of frequencies 
        df = sampling frequency (Hz)
        ng = length of time-domain smoothing window (needs to be odd)
        nh = length of frequency-domain smoothing window (needs to be odd) 
        sigmanh = std of window h, ie full width half max of gaussian 
        sigmang = std of window g, ie full width half max of gaussian
    
    Outputs:
        tfarray = SPWVD estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm>fn:
            fm,fn=fx.shape
    except ValueError:
        fn=len(fx)
        fm=1
    if fm>1:
        print 'computing cross spectra'
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx[0])
        fb=wvdas(fx[1])

    else:
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx)        
        fa=sps.hilbert(dctrend(fx))
        fb=fa.copy()
        print 'Computed Analytic signal'    
    
    #make sure window length is odd
    if nh==None:
        nh=np.floor(fn/2.)
    #make sure the window length is odd
    if np.remainder(nh,2)==0:
        nh=nh+1
    #calculate length for time smoothing window
    if ng==None:
        ng=np.floor(fn/5.)
    if np.remainder(ng,2)==0:
        ng=ng+1
    nh=int(nh)
    ng=int(ng)
    print 'nh= ',nh
    print 'ng= ',ng    
    
    dt=1./(df*2.)
    
    
    #get length of input time series 
    nfx=len(fa)
    
    #make frequency smoothing window
    if sigmanh==None:
        sigmanh=nh/(5*np.sqrt(2*np.log(2)))
    h=sps.gaussian(nh,sigmanh)
    h=h/sum(h)
    
    #make a time smoothing window
    if sigmang==None:
        sigmang=ng/(5*np.sqrt(2*np.log(2)))
    g=sps.gaussian(ng,sigmang)
    
    mlst=np.arange(start=-nh/2+1,stop=nh/2+1,step=1,dtype='int')
    #mlst=np.arange(nh,dtype='int')
    tlst=np.arange(start=nh/2,stop=nfx-nh/2,step=tstep)
    #make a frequency list for plotting exporting only positive frequencies
    flst=np.fft.fftfreq(nfbins,dt)[nfbins/2:]#get only positive frequencies
    flst[-1]=0
    flstp=np.fft.fftfreq(nfbins,2*dt)[0:nfbins/2]
    
    #create an empty array to put the tf in 
    tfarray=np.zeros((nfbins/2,len(tlst)),dtype='complex128')

    
    for tpoint,nn in enumerate(tlst):
        #calculate windowed correlation function of analytic function
        fxwin=h*fa[nn+mlst]*fb[nn-mlst].conj()
        for fpoint,mm in enumerate(flst):
            fxmed=np.convolve(g,fxwin*np.exp(1j*4*np.pi*mlst*mm*dt),
                              mode='same')/(nh*ng)
            fxmedpoint=np.median(fxmed.real)
            if fxmedpoint==0.0:
                tfarray[fpoint,tpoint]=1E-10
            else:
                tfarray[fpoint,tpoint]=fxmedpoint
    
    tfarray=(4.*nh/dt)*tfarray


    return tfarray,tlst,flstp

def specwv(fx,tstep=2**5,nfbins=2**10,nhs=2**8,nhwv=2**9-1,ngwv=2**3-1,df=1.0):
    """
    specwv(f,tstep=2**5,nfbins=2**10,nh=2**8-1,ng=1,df=1.0) will calculate 
    the Wigner-Ville distribution mulitplied by the STFT windowed by the common 
    gaussian window h for a function f. 

    Inputs:
        fx = array to compute the specwv
        tstep = time step between windows 
        nfbins = number of frequencies 
        nhs = length of time-domain smoothing window for STFT should be even
        nhwv = length of time-domain smoothing window for WV (needs to be odd) 
        ngwv = lenght of frequency-domain smoothing window (needs to be odd) 
        df = sampling frequency (Hz)
    
    Outputs:
        tfarray = SPECWV estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    #calculate stft
    pst,tlst,flst=stft(fx,nh=nhs,tstep=tstep,nfbins=nfbins,df=df)
    
    #calculate new time step so WVD and STFT will align
    ntstep=len(fx)/(len(tlst)*2.)
    
    #calculate spwvd
    pwv,twv,fwv=spwvd(fx,tstep=ntstep,nfbins=nfbins,df=df,nh=nhwv,ng=ngwv)
    
    #multiply the two together normalize
    tfarray=pst/pst.max()*pwv/pwv.max()
        
    return tfarray,tlst,flst
    
def modifiedb(fx,tstep=2**5,nfbins=2**10,df=1.0,nh=2**8-1,beta=.2):
    """modifiedb(fx,tstep=2**5,nfbins=2**10,df=1.0,nh=2**8-1,beta=.2)
    will calculate the modified b distribution as defined by cosh(n)^-2 beta 
    for a function fx. 

    Inputs:
        fx = array from which modifiedb will be calculated if computing cross
            spectra input as [fx1,fx2]
        tstep = time step between windows 
        nfbins = number of frequencies 
        df = sampling frequency (Hz)
        nh = length of time-domain smoothing window (needs to be odd)
        beta = smoothing coefficient
        
    Outputs:
        tfarray = modifiedB estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm>fn:
            fm,fn=fx.shape
    except ValueError:
        fn=len(fx)
        fm=1
    if fm>1:
        fn=fn[0]
        print 'computing cross spectra'
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx[0])
        fb=wvdas(fx[1])

    else:
        #compute the analytic signal of function f and dctrend
        fa=wvdas(fx)        
        fa=sps.hilbert(dctrend(fx))
        fb=fa.copy()
    
    #sampling period
    df=float(df)
    dt=1./df
    
    tau=(nh-1)/2  #midpoint index of window h
    
    #create a time array such that the first point is centered on time window
    tlst=np.arange(start=0,stop=fn-1,step=tstep,dtype='int')
	
    #create an empty array to put the tf in 
    tfarray=np.zeros((nfbins,len(tlst)),dtype='complex')
    
    #create a frequency array with just positive frequencies
    flst=np.fft.fftfreq(nfbins,dt)[0:nfbins/2]
    
    #calculate pseudo WV
    for point,nn in enumerate(tlst):
        #calculate the smallest timeshift possible
        taun=min(nn,tau,fn-nn-1)
        #make a timeshift array
        taulst=np.arange(start=-taun,stop=taun+1,step=1,dtype='int')
        #create modified b window
        mbwin=np.cosh(taulst)**(-2*beta)
        mbwin=mbwin/sum(mbwin)
        MBwin=np.fft.fft(padzeros(mbwin,npad=nfbins))
        #calculate windowed correlation function of analytic function
        Rnn=np.conjugate(fa[nn-taulst])*fb[nn+taulst]  
        #calculate fft of windowed correlation function
        FTRnn=MBwin*np.fft.fft(padzeros(Rnn,npad=nfbins))
        #put into tfarray
        tfarray[:,point]=FTRnn[::-1]
        
    #need to cut the time frequency array in half due to the WVD assuming 
    #time series sampled at twice nyquist.
    tfarray=tfarray
        
    return tfarray,tlst,flst

def robuststftMedian(fx,nh=2**8,tstep=2**5,df=1.0,nfbins=2**10):
    """
    robuststftMedian(fx,nh=2**8,tstep=2**5,ng=1,df=1.0) will output an array 
    of the time-frequency robust spectrogram calculated using the vector median
    simplification.
     
    Inputs:
        fx = the function to have a spectrogram computed for 
        nh = window length for each time step 
        tstep = time step between short windows 
        df = sampling frequency
        nfbins = number of frequency bins
    
    Outputs:
        tfarray = WVD estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    #get length of input time series 
    nfx=len(fx)
     
    #compute time shift list
    mlst=np.arange(start=-nh/2+1,stop=nh/2+1,step=1,dtype='int')
    #compute time locations to take STFT
    tlst=np.arange(start=0,stop=nfx-nh+1,step=tstep)
    
    #make a frequency list for plotting exporting only positive frequencies
    flst=np.fft.fftfreq(nfbins,1/df)
    flstc=flst[nfbins/2:]
    #Note: these are actually the negative frequencies but works better for
    #calculations
    flstp=flst[0:nfbins/2]
    
    #make time window and normalize
    sigmanh=nh/(6*np.sqrt(2*np.log(2)))
    h=sps.gaussian(nh,sigmanh)
    h=h/sum(h)
    
    #create an empty array to put the tf in and initialize a complex value
    tfarray=np.zeros((nfbins/2,len(tlst)),dtype='complex128')
    
    #take the hilbert transform of the signal to make complex and remove
    #negative frequencies
    fa=sps.hilbert(dctrend(fx))
    fa=fa/fa.std()
    
    #make a frequency list for plotting exporting only positive frequencies
    flst=np.fft.fftfreq(nfbins,1/df)[nfbins/2:]#get only positive frequencies
    
    for tpoint,nn in enumerate(tlst):
        #calculate windowed correlation function of analytic function
        fxwin=h*fa[nn:nn+nh]
        for fpoint,mm in enumerate(flstc):
            fxmed=fxwin*np.exp(1j*2*np.pi*mlst*mm/df)
            fxmedreal=np.median(fxmed.real)
            fxmedimag=np.median(fxmed.imag)
            if fxmedreal+1j*fxmedimag==0.0:
                tfarray[fpoint,tpoint]=1E-10
            else:
                tfarray[fpoint,tpoint]=fxmedreal+1j*fxmedimag
    #normalize tfarray
    tfarray=(4.*nh*df)*tfarray
        
    return tfarray,tlst,flstp

def robuststftL(fx,alpha=.325, nh=2**8,tstep=2**5,df=1.0,nfbins=2**10):
    """
    robuststftL(fx,nh=2**8,tstep=2**5,ng=1,df=1.0) will output an array of the
    time-frequency robust spectrogram by estimating the vector median and 
    summing terms estimated by alpha coefficients.
    
    Inputs:
        fx = the function to have a spectrogram computed for 
        alpha = robust parameter [0,.5] -> 0 gives spectrogram, .5 gives median stft
        nh = window length for each time step 
        tstep = time step between short windows 
        df = sampling frequency
        nfbins = number of frequency bins
    
    Outputs:
        tfarray = robust L-estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    #get length of input time series 
    nfx=len(fx)
     
    #compute time shift list
    mlst=np.arange(start=-nh/2+1,stop=nh/2+1,step=1,dtype='int')
    #compute time locations to take STFT
    tlst=np.arange(start=0,stop=nfx-nh+1,step=tstep)
    
    #make a frequency list for plotting exporting only positive frequencies
    flst=np.fft.fftfreq(nfbins,1/df)
    flstc=flst[nfbins/2:]
    #Note: these are actually the negative frequencies but works better for
    #calculations
    flstp=flst[0:nfbins/2]
    
    #make time window and normalize
    sigmanh=nh/(6*np.sqrt(2*np.log(2)))
    h=sps.gaussian(nh,sigmanh)
    h=h/sum(h)
    
    #create an empty array to put the tf in and initialize a complex value
    tfarray=np.zeros((nfbins/2,len(tlst)),dtype='complex128')
    
    #take the hilbert transform of the signal to make complex and remove
    #negative frequencies
    fa=sps.hilbert(dctrend(fx))
    fa=fa/fa.std()
    
    #make a frequency list for plotting exporting only positive frequencies
    flst=np.fft.fftfreq(nfbins,1/df)[nfbins/2:]#get only positive frequencies
    
    #create list of coefficients
    a=np.zeros(nh)
    a[(nh-2)*alpha:alpha*(2-nh)+nh-1]=1./(nh*(1-2*alpha)+4*alpha)
    
    for tpoint,nn in enumerate(tlst):
        #calculate windowed correlation function of analytic function
        fxwin=h*fa[nn:nn+nh]
        for fpoint,mm in enumerate(flstc):
            fxelement=fxwin*np.exp(1j*2*np.pi*mlst*mm/df)
            fxreal=np.sort(fxelement.real)[::-1]
            fximag=np.sort(fxelement.imag)[::-1]
            tfpoint=sum(a*(fxreal+1j*fximag))
            if tfpoint==0.0:
                tfarray[fpoint,tpoint]=1E-10
            else:
                tfarray[fpoint,tpoint]=tfpoint
    #normalize tfarray
    tfarray=(4.*nh*df)*tfarray
        
    return tfarray,tlst,flstp

def smethod(fx,L=11,nh=2**8,tstep=2**7,ng=1,df=1.0,nfbins=2**10,sigmaL=None):
    """
    smethod(fx,L=11,nh=2**8,tstep=2**7,ng=1,df=1.0,nfbins=2**10) will calculate
    the smethod by estimating the STFT first and computing the WV of window 
    length L in the frequency domain.  For larger L more of WV estimation, if 
    L=0 get back STFT
    
    Inputs:
        fx = the function to have a S-methoc computed for, if computing cross
           spectra input as [fx1,fx2] 
        L = window length in frequency domain
        nh = window length for each time step 
        tstep = time step between short windows 
        ng = smoothing window along frequency plane should be odd
        df = sampling frequency 
        nfbins = number of frequency bins
    
    Outputs:
        tfarray = S-method estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    
    """
    	
    df=float(df)
    
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm>fn:
            fm,fn=fx.shape
    except ValueError:
        fn=len(fx)
        fm=1
    if fm>1:
        print 'computing cross spectra'
        #compute the analytic signal of function f and dctrend
        #fa=sps.hilbert(dctrend(fx[0]))
        #fb=sps.hilbert(dctrend(fx[1]))
        fa=fx[0]
        fb=fx[1]
        fa=fa.reshape(fn)
        fb=fb.reshape(fn)
        pxa,tlst,flst=stft(fa,nh=nh,tstep=tstep,ng=ng,df=df,nfbins=nfbins)
        pxb,tlst,flst=stft(fb,nh=nh,tstep=tstep,ng=ng,df=df,nfbins=nfbins)
        pxx=pxa*pxb.conj()
    else:
        #compute the analytic signal of function f and dctrend
        #fa=sps.hilbert(dctrend(fx))
        fa=fx
        fa=fa.reshape(fn)
        fb=fa
        pxx,tlst,flst=stft(fa,nh=nh,tstep=tstep,ng=ng,df=df,nfbins=nfbins)
#        pxb=pxa

    #make an new array to put the new tfd in
    tfarray=abs(pxx)**2
    #get shape of spectrogram
    nf,nt=tfarray.shape
    #create a list of frequency shifts
    Llst=np.arange(start=-L/2+1,stop=L/2+1,step=1,dtype='int')
    #create a frequency gaussian window
    if sigmaL==None:
        sigmaL=L/(1*np.sqrt(2*np.log(2)))
    p=sps.gaussian(L,sigmaL)
    #make a matrix of windows
    pm=np.zeros((L,nt))
    for kk in range(nt):
        pm[:,kk]=p
    
    #loop over frequency and calculate the s-method    
    for ff in range(L/2,nf-L/2):
        tfarray[ff,:]=tfarray[ff,:]+2*np.real(np.sum(pm*pxx[ff+Llst,:]*
                                             pxx[ff-Llst,:].conj(),axis=0))
    tfarray[L/2:-L/2]=tfarray[L/2:-L/2]/L
    
    return tfarray,tlst,flst,pxx
    
def robustSmethod(fx,L=5,nh=2**7,tstep=2**5,nfbins=2**10,df=1.0,
                  robusttype='median',sigmal=None):
    """
    robustSmethod(fx,L=15,nh=2**7,tstep=2**5,nfbins=2**10,df=1.0) computes the 
    robust Smethod via the robust spectrogram.
    
    Inputs:
        fx = array of data, if computing cross-spectra input as [fa,fb]
        L = frequency smoothing window if robusttype='median'
        nh = window length for STFT
        tstep = time step for each STFT to be computed
        nfbins = number of frequency bins to be calculate
        df = sampling frequency
        robusttype = type of robust STFT calculation can be 'median' or 'L'
        simgal = full-width half max of gaussian window applied in frequency
    
    Outputs:
        tfarray = robust S-method estimation of array fx
        tlst = time instances of each calculation
        flst = array of positive frequencies
    """
    
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm>fn:
            fm,fn=fx.shape
    except ValueError:
        fn=len(fx)
        fm=1
    if fm>1:
        print 'computing cross spectra'
        #compute the analytic signal of function f and dctrend
        fa=fx[0].reshape(fn)
        fb=fx[1].reshape(fn)
        if robusttype=='median':
            pxa,tlst,flst=robuststftMedian(fa,nh=nh,tstep=tstep,df=df,
                                           nfbins=nfbins)
            pxb,tlst,flst=robuststftMedian(fb,nh=nh,tstep=tstep,df=df,
                                           nfbins=nfbins)
        elif robusttype=='L':
            pxa,tlst,flst=robuststftL(fa,nh=nh,tstep=tstep,df=df,nfbins=nfbins)
            pxb,tlst,flst=robuststftL(fb,nh=nh,tstep=tstep,df=df,nfbins=nfbins)
        else:
            raise ValueError('robusttype undefined')
        pxx=pxa*pxb.conj()
    else:
        fa=fx.reshape(fn)
        if robusttype=='median':
            pxx,tlst,flst=robuststftMedian(fa,nh=nh,tstep=tstep,df=df,
                                           nfbins=nfbins)
        elif robusttype=='L':
            pxx,tlst,flst=robuststftL(fa,nh=nh,tstep=tstep,df=df,nfbins=nfbins)
        else:
            raise ValueError('robusttype undefined')
    
    #compute frequency shift list
    Llst=np.arange(start=-L/2+1,stop=L/2+1,step=1,dtype='int')

    #compute the frequency window of length L
    if sigmal==None:    
        sigmal=L/3*(np.sqrt(2*np.log(2)))        
    lwin=gausswin(L,sigmal)
    lwin=lwin/sum(lwin)
    pm=np.zeros((L,len(tlst)))
    for kk in range(len(tlst)):
        pm[:,kk]=lwin
    
    smarray=pxx.copy()
    #compute S-method
    for ff in range(L/2,nfbins/2-L/2):
        smarray[ff,:]=smarray[ff,:]+2*np.real(np.sum(pm*pxx[ff+Llst,:]*
                                            pxx[ff-Llst,:].conj(),axis=0))
    
#    for tt in range(len(tlst)):
#        for kk in range((L-1)/2,len(flst)-(L-1)/2):
#            smarray[kk,tt]=abs(pxx[kk,tt])+np.sqrt(abs(2*sum(lwin*
#   pxx[kk+Llst,tt]*pxx[kk-Llst,tt].conj())))
    smarray=(2./(L*nh))*smarray
    
    return smarray,tlst,flst,pxx
    
def reassignedSmethod(fx,nh=2**7-1,tstep=2**4,nfbins=2**9,df=1.0,alpha=4,
                      thresh=.01,L=5):
    """
    reassignedSmethod(fx,nh=2**7-2,tstep=2**4,nfbins=2**9,df=1.0,alpha=4,
                      thresh=.05,L=5) 
    will calulate the reassigned S-method as described by Djurovic[1999] by 
    using the spectrogram to estimate the reassignment
    
    Inputs:
        fx = 1-d array to be processed
        nh = window length for each time instance
        tstep = step between time instances
        nfbins = number of frequency bins, note output will be nfbins/2 due to 
                symmetry of the FFT
        df = sampling rate (Hz)
        alpha = inverse of full-width half max of gaussian window, smaller 
                numbers mean broader windows
        thresh = threshold for reassignment, lower numbers more points
                reassigned, higer numbers less points reassigned
        L = length of window for S-method calculation, higher numbers tend 
            tend toward WVD
    Outputs:
        rtfarray = reassigned S-method shape of (nfbins/2,len(fx)/tstep)
        tlst = list of time instances where rtfarray was calculated 
        flst = positive frequencies
        sm = S-method array
        
    """  
    
    if type(fx) is list:
        fx=np.array(fx)
    try:
        fn,fm=fx.shape
        if fm>fn:
            fm,fn=fx.shape
    except ValueError:
        fn=len(fx)
        fm=1
    if fm>1:
        print 'computing cross spectra'
        #compute the analytic signal of function f and dctrend
        #fa=sps.hilbert(dctrend(fx[0]))
        #fb=sps.hilbert(dctrend(fx[1]))
        fa=fx[0]
        fb=fx[1]
        fa=fa.reshape(fn)
        fb=fb.reshape(fn)
    else:
        fa=fx
        fa=fa.reshape(fn)
        fb=fa.copy()

    
    nx=len(fx)    
    
    #compute gaussian window
    h=gausswin(nh,alpha=alpha)
    #h=np.hanning(nh)
    lh=(nh-1)/2
    
    #compute ramp window
    th=h*np.arange(start=-lh,stop=lh+1,step=1)
    
    #compute derivative of window
    dh=dwindow(h)
    
    #make a time list of indexes
    tlst=np.arange(start=0,stop=nx,step=tstep)
    nt=len(tlst)
    
    #make frequency list for plotting
    flst=np.fft.fftfreq(nfbins,1./df)[:nfbins/2]
    
    #initialize some time-frequency arrays
    tfh=np.zeros((nfbins,nt),dtype='complex128')
    tfth=np.zeros((nfbins,nt),dtype='complex128')
    tfdh=np.zeros((nfbins,nt),dtype='complex128')
    
    #compute components for reassignment
    for ii,tt in enumerate(tlst):
        #create a time shift list
        tau=np.arange(start=-min([np.round(nx/2.),lh,tt-1]),
                   stop=min([np.round(nx/2.),lh,nx-tt-1])+1)
        #compute the frequency spots to be calculated
        ff=np.remainder(nfbins+tau,nfbins)
        #make lists of data points for each window calculation
        xlst=tt+tau
        hlst=lh+tau
        normh=np.sqrt(np.sum(abs(h[hlst])**2))
        tfh[ff,ii]=fx[xlst]*h[hlst].conj()/normh
        tfth[ff,ii]=fx[xlst]*th[hlst].conj()/normh
        tfdh[ff,ii]=fx[xlst]*dh[hlst].conj()/normh
    
    #compute Fourier Transform
    spech=np.fft.fft(tfh,axis=0)
    specth=np.fft.fft(tfth,axis=0)
    specdh=np.fft.fft(tfdh,axis=0)
    
    #get only positive frequencies
    spech=spech[nfbins/2:,:]
    specth=specth[nfbins/2:,:]
    specdh=specdh[nfbins/2:,:]
    
    #check to make sure no spurious zeros floating around
    szf=np.where(abs(spech)<1.E-6)
    spech[szf]=0.0+0.0j
    zerofind=np.nonzero(abs(spech))
    twspec=np.zeros((nfbins/2,nt),dtype='float')
    dwspec=np.zeros((nfbins/2,nt),dtype='float')
    twspec[zerofind]=np.round(np.real(specth[zerofind]/spech[zerofind]))
    dwspec[zerofind]=np.round(np.imag((nfbins/2.)*specdh[zerofind]/
                                spech[zerofind])/(np.pi))
    
    #get shape of spectrogram
    nf,nt=spech.shape
    
    #-----calculate s-method-----
    Llst=np.arange(start=-L/2+1,stop=L/2+1,step=1,dtype='int')

    #make and empty array of zeros
    sm=np.zeros_like(spech)
    
    #put values where L cannot be value of L, near top and bottom
    sm[0:L/2,:]=abs(spech[0:L/2,:])**2
    sm[-L/2:,:]=abs(spech[-L/2:,:])**2

    #calculate s-method
    for ff in range(L/2,nf-L/2-1):
        sm[ff,:]=2*np.real(np.sum(spech[ff+Llst,:]*spech[ff-Llst,:].conj(),
                            axis=0))/L
    
    #------compute reassignment-----    

    
    rtfarray=np.zeros((nfbins/2,nt))
    
    threshold=thresh*np.max(abs(sm))
    
    for nn in range(nt):
        for kk in range(nf):
            if abs(spech[kk,nn])>threshold:
                #get center of gravity index in time direction from spectrogram           
                nhat=int(nn+twspec[kk,nn])
                nhat=int(min([max([nhat,1]),nt-1]))
                #get center of gravity index in frequency direction from spec
                khat=int(kk-dwspec[kk,nn])
                khat=int(np.remainder(np.remainder(khat-1,nfbins/2)+nfbins/2,
                                      nfbins/2))
                rtfarray[khat,nhat]=rtfarray[khat,nhat]+abs(sm[kk,nn])
            else:
                rtfarray[kk,nn]=rtfarray[kk,nn]+sm[kk,nn]

    #place values where L cannot be L    
    rtfarray[:L/2,:]=abs(sm[:L/2,:])
    rtfarray[-L/2:,:]=abs(sm[-L/2:,:])
    
    tz=np.where(rtfarray==0)
    rtfarray[tz]=1.0
    
    tz=np.where(sm==0.0)
    sm[tz]=1.0    
    
    #scale
    rtfarray=abs(rtfarray)
    
    return rtfarray,tlst,flst,sm
    
    
def plottf(tfarray,tlst,flst,fignum=1,starttime=0,timeinc='hrs',
           dt=1.0,title=None,vmm=None,cmap=None,aspect=None,interpolation=None,
           cbori=None,cbshrink=None,cbaspect=None,cbpad=None,powscale='log',
           normalize='n',yscale='log',period='n'):
    """plottf(tfarray,tlst,flst,fignum=1) will plot a calculated tfarray with 
    limits corresponding to tlst and flst. 
   
    Inputs:
    
        starttime = starttime measured in timeincrement 
        tinc = 'hrs','min' or 'sec' 
        vmm = [vmin,vmax] a list for min and max 
        title = title string 
        cmap = colormap scheme default is jet, type help on matplotlib.cm 
        aspect = aspect of plot, default is auto, can be 'equal' or a scalar 
        interpolation = type of color interpolation, type help on 
            matplotlib.pyplot.imshow 
        cbori = colorbar orientation 'horizontal' or 'vertical' 
        cbshrink = percentage of 1 for shrinking colorbar 
        cbaspect = aspect ratio of long to short dimensions 
        cbpad = pad between colorbar and axis
        powscale = linear or log for power
        normalize = y or n, yes for normalization, n for no
        yscale = linear or log plot yscale
        period = 'y' or 'n' to plot in period instead of frequency
    Outputs:
        plot
     """
    
    #time increment
    if timeinc=='hrs':
        tinc=3600/dt
    elif timeinc=='min':
        tinc=60/dt
    elif timeinc=='sec':
        tinc=1/dt
    else:
        raise ValueError(timeinc+'is not defined')
    #colormap
    if cmap==None:
        cmap='jet'
    else:
        cmap=cmap
    #aspect ratio
    if aspect==None:
        aspect='auto'
    else:
        aspect=aspect
    #interpolation
    if interpolation==None:
        interpolation='gaussian'
    else:
        interpolation=interpolation
    #colorbar orientation
    if cbori==None:
        cbori='vertical'
    else:
        cbori=cbori
    #colorbar shinkage
    if cbshrink==None:
        cbshrink=.8
    else:
        cbshrink=cbshrink
    #colorbar aspect
    if cbaspect==None:
        cbaspect=20
    else:
        cbaspect=cbaspect
    #colorbar pad
    if cbpad==None:
        cbpad=.05
    else:
        cbpad=cbpad
    #scale
    if powscale=='log':
        zerofind=np.where(abs(tfarray)==0)
        tfarray[zerofind]=1.0
        if normalize=='y':
            plottfarray=10*np.log10(abs(tfarray/np.max(abs(tfarray))))
        else:
            plottfarray=10*np.log10(abs(tfarray))
    elif powscale=='linear':
        if normalize=='y':
            plottfarray=abs(tfarray/np.max(abs(tfarray)))
        else:
            plottfarray=abs(tfarray)

    #period or frequency
    if period=='y':
        flst[1:]=1./flst[1:]
        flst[0]=2*flst[1]
    elif period=='n':
        pass
    
    #set properties for the plot         
    plt.rcParams['font.size']=9
    plt.rcParams['figure.subplot.left']=.12
    plt.rcParams['figure.subplot.right']=.99
    plt.rcParams['figure.subplot.bottom']=.12
    plt.rcParams['figure.subplot.top']=.96
    plt.rcParams['figure.subplot.wspace']=.25
    plt.rcParams['figure.subplot.hspace']=.20
    
    #set the font dictionary
    fdict={'size':10,'weight':'bold'}    
    
    #make a meshgrid if yscale is logarithmic    
    if yscale=='log':
        logt,logf=np.meshgrid(tlst,flst)
    
    #make figure    
    fig1=plt.figure(fignum,[10,10],dpi=300)
    ax=fig1.add_subplot(1,1,1)
    if vmm!=None:
        vmin=vmm[0]
        vmax=vmm[1]
        #add in log yscale
        if yscale=='log':
            #need to flip the matrix so that origin is bottom right
            cbp=ax.pcolormesh(logt,logf,np.flipud(plottfarray),
                          cmap=cmap,vmin=vmin,vmax=vmax)
            ax.semilogy()
            ax.set_ylim(flst[1],flst[-1])
            ax.set_xlim(tlst[0],tlst[-1])
            cb=plt.colorbar(cbp,orientation=cbori,shrink=cbshrink,pad=cbpad,
                            aspect=cbaspect)
        else:    
            plt.imshow(plottfarray,extent=(tlst[0]/tinc+starttime,
                       tlst[-1]/tinc+starttime,flst[1],flst[-1]),aspect=aspect,
                       vmin=vmin,vmax=vmax,cmap=cmap,
                       interpolation=interpolation)
            cb=plt.colorbar(orientation=cbori,shrink=cbshrink,pad=cbpad,
                            aspect=cbaspect)
    else:
        if yscale=='log':
            cbp=ax.pcolormesh(logt,logf,np.flipud(plottfarray),
                          cmap=cmap)
            ax.semilogy()
            ax.set_ylim(flst[1],flst[-1])
            ax.set_xlim(tlst[0],tlst[-1])
            cb=plt.colorbar(cbp,orientation=cbori,shrink=cbshrink,pad=cbpad,
                            aspect=cbaspect)
        else:
            
            plt.imshow(plottfarray,extent=(tlst[0]/tinc+starttime,
                       tlst[-1]/tinc+starttime,flst[1],flst[-1]),aspect=aspect,
                       cmap=cmap,interpolation=interpolation)
            cb=plt.colorbar(orientation=cbori,shrink=cbshrink,pad=cbpad,
                            aspect=cbaspect)
    ax.set_xlabel('time('+timeinc+')',fontdict=fdict)
    if period=='y':
        ax.set_ylabel('period (s)',fontdict=fdict)
    else:
        ax.set_ylabel('frequency (Hz)',fontdict=fdict)
    if title!=None:
        ax.set_title(title,fontdict=fdict)
    
    plt.show()
    
    
def plotAll(fx,tfarray,tlst,flst,fignum=1,starttime=0,timeinc='hrs',
           dt=1.0,title=None,vmm=None,cmap=None,aspect=None,interpolation=None,
           cbori=None,cbshrink=None,cbaspect=None,cbpad=None,normalize='n',
           scale='log'):
    """plottf(tfarray,tlst,flst,fignum=1) will plot a calculated tfarray with 
    limits corresponding to tlst and flst. Can have: 
    
    Inputs:
        starttime = starttime measured in timeincrement 
        timeincrement = 'hrs','min' or 'sec' 
        vmm = [vmin,vmax] a list for min and max 
        title = title string 
        cmap = colormap scheme default is jet, type help on matplotlib.cm 
        aspect = aspect of plot, default is auto, can be 'equal' or a scalar 
        interpolation = type of color interpolation, type help on 
                     matplotlib.pyplot.imshow 
        cbori = colorbar orientation 'horizontal' or 'vertical' 
        cbshrink = percentage of 1 for shrinking colorbar 
        cbaspect = aspect ratio of long to short dimensions 
        cbpad = pad between colorbar and axis
        normalization = y or n, y for normalization n for none
    
    Outputs:
        plot
    """
    
    #time increment
    if timeinc=='hrs':
        tinc=3600/dt
    elif timeinc=='min':
        tinc=60/dt
    elif timeinc=='sec':
        tinc=1/dt
    else:
        raise ValueError(timeinc+'is not defined')
    #colormap
    if cmap==None:
        cmap='jet'
    else:
        cmap=cmap
    #aspect ratio
    if aspect==None:
        aspect='auto'
    else:
        aspect=aspect
    #interpolation
    if interpolation==None:
        interpolation='gaussian'
    else:
        interpolation=interpolation
    #colorbar orientation
    if cbori==None:
        cbori='vertical'
    else:
        cbori=cbori
    #colorbar shinkage
    if cbshrink==None:
        cbshrink=.99
    else:
        cbshrink=cbshrink
    #colorbar aspect
    if cbaspect==None:
        cbaspect=20
    else:
        cbaspect=cbaspect
    #colorbar pad
    if cbpad==None:
        cbpad=.1
    else:
        cbpad=cbpad
        
    #scale
    if scale=='log':
        zerofind=np.where(abs(tfarray)==0)
        tfarray[zerofind]=1.0
        if normalize=='y':
            plottfarray=20*np.log10(abs(tfarray/np.max(abs(tfarray))))
        else:
            plottfarray=20*np.log10(abs(tfarray))
    elif scale=='linear':
        if normalize=='y':
            plottfarray=abs(plottfarray/np.max(abs(plottfarray)))**2
        else:
            plottfarray=abs(tfarray)**2
        
    t=np.arange(len(fx))*dt+starttime*dt
    FX=np.fft.fft(padzeros(fx))
    FXfreq=np.fft.fftfreq(len(FX),dt)
        
    #set some plot parameters
    plt.rcParams['font.size']=10
    plt.rcParams['figure.subplot.left']=.13
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.07
    plt.rcParams['figure.subplot.top']=.96
    plt.rcParams['figure.subplot.wspace']=.25
    plt.rcParams['figure.subplot.hspace']=.20
    #plt.rcParams['font.family']='helvetica'
    
    fig=plt.figure(fignum)
    plt.clf()
    
    #plot FFT of fx
    fax=fig.add_axes([.05,.25,.1,.7])
    plt.semilogx(abs(FX[0:len(FX)/2]/max(abs(FX))),FXfreq[0:len(FX)/2],'-k')
    plt.axis('tight')
    plt.ylim(0,FXfreq[len(FX)/2-1])
#    fax.xaxis.set_major_locator(MultipleLocator(.5))
    
    #plot TFD
    pax=fig.add_axes([.25,.25,.75,.7])
    if vmm!=None:
        vmin=vmm[0]
        vmax=vmm[1]
        plt.imshow(plottfarray,extent=(tlst[0]/tinc,tlst[-1]/tinc,
               flst[0],flst[-1]),aspect=aspect,vmin=vmin,vmax=vmax,cmap=cmap,
               interpolation=interpolation)
    else:
        plt.imshow(plottfarray,extent=(tlst[0]/tinc,tlst[-1]/tinc,
               flst[0],flst[-1]),aspect=aspect,cmap=cmap,
               interpolation=interpolation)
    plt.xlabel('Time('+timeinc+')',fontsize=12,fontweight='bold')
    plt.ylabel('Frequency (Hz)',fontsize=12,fontweight='bold')
    if title!=None:
        plt.title(title,fontsize=14,fontweight='bold')
    plt.colorbar(orientation=cbori,shrink=cbshrink,pad=cbpad,aspect=cbaspect)
    
    #plot timeseries
    tax=fig.add_axes([.25,.05,.60,.1])
    plt.plot(t,fx,'-k')
    plt.axis('tight')
    plt.show()
    

def stfbss(X,nsources=5,ng=2**5-1,nh=2**9-1,tstep=2**6-1,df=1.0,nfbins=2**10,
          tftol=1.E-8,normalize=True):
    """
    btfssX,nsources=5,ng=2**5-1,nh=2**9-1,tstep=2**6-1,df=1.0,nfbins=2**10,
          tftol=1.E-8,normalize=True) 
    estimates sources using a blind source algorithm based on spatial 
    time-frequency distributions.  At the moment this algorithm uses the SPWVD 
    to estimate TF distributions.  
    
    Inputs:
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
    
    Outputs:
        
        Se = estimated individual signals up to a permutation and scale
        Ae = estimated mixing matrix as X=A*S
    """
    
    
    #get shape of timeseries matrix, 
    #m=number of channels
    #tlen=length of timeseries
    m,maxn=X.shape
    
    n=nsources
    #get number of time bins
    ntbins=int(float(maxn)/tstep)
    
    tfkwargs={'ng':ng,'nh':nh,'df':df,'nfbins':nfbins,'tstep':tstep}
    
    
    #remove dc component from time series and normalize
    if normalize==True:
        for ii in range(m):
            X[ii,:]=X[ii,:]-np.mean(X[ii,:])
            X[ii,:]=X[ii,:]/X[ii,:].std()   
    
    #===============================================================================
    # Whiten data and Compute Whitening matrix
    #===============================================================================
    #whiten data to get a unitary matrix with unit variance and zero mean
    #compute covariance matrix
    Rxx=(np.dot(X,X.T))/float(maxn)
    
    #calculate eigen decomposition
    [l,u]=np.linalg.eig(Rxx)
    
    #sort eigenvalues from smallest to largest assuming largest are sources and 
    #smallest are noise
    lspot=l.argsort()
    eigval=l[lspot]
    eigvec=u[:,lspot]
    
    #calculate the noise variance as mean of non-principal components
    sigman=np.mean(eigval[0:m-n])
    
    #compute scaling factor for whitening matrix
    wscale=1/np.sqrt(eigval[m-n:m]-sigman)
    
    #compute whitening matrix
    W=np.zeros((m,n))
    for kk in range(n):
        W[:,kk]=wscale[kk]*eigvec[:,m-n+kk].T
    W=W.T
    #compute whitened signal vector. Note the dimensionality is reduced from [mxn] 
    #to [nxn] making the computation simpler.
    Z=np.dot(W,X)
    
    #===============================================================================
    #  Compute Spatial Time Frequency Distribution
    #===============================================================================
    
    stfd=np.zeros((n,n,nfbins,ntbins+1),dtype='complex128')
    Za=np.array(Z.copy())
    #compute auto terms
    for ii in range(n):
        pswvd,tswvd,fswvd=spwvd(Za[ii].reshape(maxn),**tfkwargs)
        stfd[ii,ii,:,:]=pswvd
    
    #compute cross terms
    for jj in range(n):
        for kk in range(jj,n):
            pswvd,tswvd,fswvd=spwvd([Za[jj].reshape(maxn),Za[kk].reshape(maxn)],
                                        **tfkwargs)
            stfd[jj,kk,:,:]=pswvd
            stfd[kk,jj,:,:]=pswvd.conj()
    
    #===============================================================================
    # Compute criteria for cross terms 
    #===============================================================================
    
    stfdTr=np.zeros((nfbins,ntbins))
    C=np.zeros((nfbins,ntbins))
    
    for ff in range(nfbins):
        for tt in range(ntbins):
            #compensate for noise
            stfd[:,:,ff,tt]=stfd[:,:,ff,tt]-sigman*np.matrix(W)*np.matrix(W.T)
            #compute the trace
            stfdTr[ff,tt]=abs(np.trace(stfd[:,:,ff,tt]))
    #compute mean over entire t-f plane
    trmean=stfdTr.mean()
    
    #find t-f points that meet the criteria
    fspot,tspot=np.nonzero(stfdTr>trmean)
    
    for ll in range(len(fspot)):
        treig=abs(np.linalg.eig(stfd[:,:,fspot[ll],tspot[ll]])[0])
        if sum(treig)!=0 and sum(treig)>tftol:
            C[fspot[ll],tspot[ll]]=max(treig)/sum(treig)
        else:
            C[fspot[ll],tspot[ll]]=0
    
    #compute gradients and jacobi matrices
    negjacobi=np.zeros((nfbins,ntbins))
    smallgrad=np.zeros((nfbins,ntbins))
    maxpoints=np.zeros((nfbins,ntbins))
    
    gradt,gradf=np.gradient(C)
    Jtt,Jtf=np.gradient(gradt)
    Jft,Jff=np.gradient(gradf)
    
    #get points when L2 of gradient is smaller than tolerance level
    smallgrad=np.where(np.sqrt(gradt**2+gradf**2)<tftol,1,0)
    
    #get points where the Jacobi is negative definite
    detjacobi=Jtt*Jff-Jtf*Jft
    negjacobi=np.where(detjacobi>0,1,0)*np.where(Jtt<0,1,0)\
                        *np.where((Jtt+Jff)<0,1,0)
    
    maxpoints=smallgrad*negjacobi
    
    gfspot,gtspot=np.nonzero(maxpoints)
    ntfpoints=len(gfspot)
    
    if ntfpoints==0:
        raise ValueError('Found no tf points, relax tolerance')
    else:
        print 'Found '+str(ntfpoints)+' t-f points'
    
    for rr in range(ntfpoints):
        if rr==0:
            Rjd=stfd[:,:,gfspot[rr],gtspot[rr]]
        else:
            Rjd=np.concatenate((Rjd,stfd[:,:,gfspot[rr],gtspot[rr]]),axis=1)
    Rjd=np.array(Rjd)
    
    #===============================================================================
    # Calculate Joint Diagonalization
    #===============================================================================
    #get size of array of matrices to be diagonalized
    mtf,nm=Rjd.shape #mtf is number of t-f points, nm is number of matrices
    #set up some initial parameters
    V=np.eye(mtf)
    
    #update boolean
    encore=True 
    #Total number of rotations
    updates=0
    sweep=0
    #print 'Computing Joint Diagonalization'
    # Joint diagonalization proper
    # ============================
    while encore:
        #reset some parameters
        encore=False
        sweep+=1 
        upds=0
        Vkeep=V
        for p in range(mtf):
            
            for q in range(p+1,mtf):
                #set up indices
                qi=np.arange(start=q,stop=nm,step=mtf)
                pi=np.arange(start=p,stop=nm,step=mtf)
                # computation of Givens angle
                g=np.array([Rjd[p,pi]-Rjd[q,qi],Rjd[p,qi],Rjd[q,pi]])
                gg=np.real(np.dot(g,g.T))
                ton=gg[0,0]-gg[1,1] 
                toff=gg[0,1]+gg[1,0]
                theta=0.5*np.arctan2(toff,ton+np.sqrt(ton**2+toff**2))
                # Givens update
                if abs(theta) > tftol:
                    encore=True
                    upds+=1
                    c=np.cos(theta) 
                    s=np.sin(theta)
                    G=np.matrix([[c,-s],[s,c]])
                    pair =np.array([p,q])
                    V[:,pair]=V[:,pair]*G
                    Rjd[pair,:]=G.T*Rjd[pair,:]
                    Rjd[:,np.concatenate([pi,qi])]=np.append(
                                                c*Rjd[:,pi]+s*Rjd[:,qi],
                                                -s*Rjd[:,pi]+c*Rjd[:,qi],axis=1)
        updates+=upds
    print 'Updated '+str(updates)+' times.'
    
    #compute estimated signal matrix
    Se=np.dot(V.T,Z)
    #compute estimated mixing matrix
    Ae=np.dot(np.linalg.pinv(W),V)
    
    return Se,Ae
