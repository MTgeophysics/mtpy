#module for MT calculations

import math
import numpy as np
import os
import fnmatch
import time
import datetime
import shutil
import LatLongUTMconversion as utm2ll
import scipy as sp
import scipy.signal as sps

#short spaces 3 spaces
tsp='   '
#long spaces 6 spaces
lsp='      '

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
#    for ii in range(len(bfield)):
#        bfield[ii]=(float(bfield[ii])/10.E7)*70000.*float(dlgain)*float(zadj)
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
#    for ii in range(len(efield)):
#        efield[ii]=(float(efield[ii])*float(dlgain))/(float(dlength)*float(egain))
    return efieldc


def rotateE(ex,ey,declination):
    """Rotate the electric field to geographic north, ex & ey are
    lists, declination is in degrees"""
    s=math.sin(-declination/180*math.pi)
    c=math.cos(-declination/180*math.pi)
    exrot=[]
    eyrot=[]
    for ii in range(len(ex)):
        exrot.append(round(ex[ii]*c+ey[ii]*s,8))
        eyrot.append(round(ex[ii]*c-ey[ii]*s,8))
    return exrot,eyrot

def rotateB(bx,by):
    """rotateB(bx,by) will Rotate the magnetic field such that Bx is pointing to
    magnetic north and by to geomagnetic east.  Assumes setup was orthogonal. 
    Returns list of strings that are 8 significant digits."""

    meanbx=sum(bx)/len(bx)
    meanby=sum(by)/len(by)
    z=math.atan(meanby/meanbx)
    rotang=z
        
    s=math.sin(rotang)
    c=math.cos(rotang)
    bxrot=[]
    byrot=[]
    for ii in range(len(bx)):
        bxrot.append(sigfigs(str(bx[ii]*c+by[ii]*s),8))
        byrot.append(sigfigs(str(bx[ii]*c-by[ii]*s),8))
    return bxrot,byrot
    
def padzeros(f,npad=None):
    """padzeros(f) will return a function that is padded with zeros to the next
    power of 2 for faster processing for fft or to length npad if given."""
    
    n=len(f)
    if npad==None:
        pow=np.log2(n)
        fpow=np.floor(pow)
        if pow!=fpow:
            npow=fpow+1
        else:
            npow=pow
        fpad=np.zeros(2**npow)
        fpad[0:n]=f[0:n]
    else:
        fpad=np.zeros(npad)
        fpad[0:n]=f[0:n]
        
    return fpad
    
def filter(f,fcutoff=10.,w=10.0,dt=.001):
    """Will apply a sinc filter of width w to the function f by multipling in
    the frequency domain. Returns filtered function"""

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
    """dctrend(f) will remove a dc trend from the function f."""
    
    n=len(f)
    dc=sum(np.array(f))/n
    fdc=np.array([f[ii]-dc for ii in range(n)])
    
    return fdc

def normalizeL2(f):
    """normalizeL2(f) will return the function f normalized by the L2 norm ->
    f/(sqrt(sum(abs(x_i)^2)))."""
    
    f=np.array(f)
    fnorm=f/np.sqrt(np.sum(np.abs(f)**2))
    
    return fnorm

def decimatef(farray,m):
    """Will decimate a function by the factor m. First an 8th order Cheybechev 
    type I filter with a cuttoff frequency of .8/m  is applied in both 
    directions to minimize any phase distortion and remove any aliasing. If m 
    is greater than 10, decimatef will be called multiple times."""
    
    
    n=len(farray)
    fdec=sps.resample(farray,n/m,window='hanning')
#    print 'Decimating using Chebyshev'
#    n=len(farray)
#    nout=np.ceil(n/m)
#    nfilt=8
#    rip=.05
#    
#    if m>10:
#        mlst=[]
#        mlst.append(10)
#        for mm in range(1,int(np.ceil(np.log10(m)))):
#            if mm<np.floor(np.log10(m)):
#                mlst.append(10)
#            else: 
#                mlst.append(np.floor(m/10**mm))
#        
#        #make a cheybeshev1 zero-phase filter with cuttoff frequency of .8/m
#        b,a=sps.iirfilter(8,.8/mlst[0],rp=.05,btype='low',ftype='cheby1',output='ba')
#        #check to make sure the filter will produce good results
#        nb=len(b)
#        na=len(a)
#        top=sum(np.exp(-1j*np.arange(nb)*np.pi*.8/mlst[0])*b)
#        bottom=sum(np.exp(-1j*np.arange(na)*np.pi*.8/mlst[0])*a)
#        H=abs(20*np.log10(abs(top/bottom)))
#        while b[0]<10**-5 or H>10**-6:
#            nfilt=nfilt-1
#            if nfilt==0:
#                ValueError("length of filter is 0 decimation value might be to large")
#            else:
#                b,a=sps.iirfilter(nfilt,.8/mlst[0],rp=rip,btype='low',ftype='cheby1',output='ba')
#                nb=len(b)
#            na=len(a)
#            top=sum(np.exp(-1j*np.arange(nb)*np.pi*.8/mlst[0])*b)
#            bottom=sum(np.exp(-1j*np.arange(na)*np.pi*.8/mlst[0])*a)
#            H=abs(20*np.log10(abs(top/bottom))+rip)
#            
#        
#        ffilt=sps.filtfilt(b,a,farray)
#        nout=np.ceil(n/mlst[0])
#        nbeg=n-mlst[0]*nout
#        fdec=np.array([ffilt[ii] for ii in np.arange(start=nbeg,stop=int(n),step=mlst[0])])
#        print 'Decimated once'
#        
#        for dd in range(1,len(mlst)):
#            #make a cheybeshev1 zero-phase filter with cuttoff frequency of .8/m
#            b,a=sps.iirfilter(8,.8/mlst[dd],rp=.05,btype='low',ftype='cheby1',output='ba')
#            #check to make sure the filter will produce good results
#            nb=len(b)
#            na=len(a)
#            top=sum(np.exp(-1j*np.arange(nb)*np.pi*.8/mlst[dd])*b)
#            bottom=sum(np.exp(-1j*np.arange(na)*np.pi*.8/mlst[dd])*a)
#            H=abs(20*np.log10(abs(top/bottom)))
#            while b[0]<10**-5 or H>10**-6:
#                nfilt=nfilt-1
#                if nfilt==0:
#                    ValueError("length of filter is 0 decimation value might be to large")
#                else:
#                    b,a=sps.iirfilter(nfilt,.8/mlst[dd],rp=rip,btype='low',ftype='cheby1',output='ba')
#                    nb=len(b)
#                na=len(a)
#                top=sum(np.exp(-1j*np.arange(nb)*np.pi*.8/mlst[dd])*b)
#                bottom=sum(np.exp(-1j*np.arange(na)*np.pi*.8/mlst[dd])*a)
#                H=abs(20*np.log10(abs(top/bottom))+rip)
#                
#            
#            ffilt=sps.filtfilt(b,a,fdec)
#            n=len(fdec)
#            nout=np.ceil(n/mlst[dd])
#            nbeg=n-mlst[dd]*nout
#            fdec=np.array([ffilt[ii] for ii in np.arange(start=nbeg,stop=int(n),step=mlst[dd])])
#            print 'Decimated '+str(dd+1)+' by '+str(mlst[dd])
#    else:
#        #make a cheybeshev1 zero-phase filter with cuttoff frequency of .8/m
#        b,a=sps.iirfilter(8,.8/m,rp=.05,btype='low',ftype='cheby1',output='ba')
#        #check to make sure the filter will produce good results
#        nb=len(b)
#        na=len(a)
#        top=sum(np.exp(-1j*np.arange(nb)*np.pi*.8/m)*b)
#        bottom=sum(np.exp(-1j*np.arange(na)*np.pi*.8/m)*a)
#        H=abs(20*np.log10(abs(top/bottom)))
#        while b[0]<10**-5 or H>10**-6:
#            nfilt=nfilt-1
#            if nfilt==0:
#                ValueError("length of filter is 0 decimation value might be to large")
#            else:
#                b,a=sps.iirfilter(nfilt,.8/m,rp=rip,btype='low',ftype='cheby1',output='ba')
#                nb=len(b)
#            na=len(a)
#            top=sum(np.exp(-1j*np.arange(nb)*np.pi*.8/m)*b)
#            bottom=sum(np.exp(-1j*np.arange(na)*np.pi*.8/m)*a)
#            H=abs(20*np.log10(abs(top/bottom))+rip)
#            
#        
#        ffilt=sps.filtfilt(b,a,farray)
#        nout=np.ceil(n/m)
#        nbeg=n-m*nout
#        fdec=np.array([ffilt[ii] for ii in np.arange(start=nbeg,stop=int(n),step=m)])
            
    return fdec
    
def openMTfile(filename,gain=2.5,egain=10,dlength=[100,100],magtype='lp',zadj=1):
    """Open an MT file, convert counts to units and return an 1D-array.  
    Make sure filename includes the entire path.  gain (verylow=2.5,low=1,
    high=.1), egain same as gain, dlength is dipole length in m of EX,EY. 
    magtype is the magnetic measurment type 'lp' for long period and 'bb' for 
    broadband, zadj is an adjustment for lp instruments on the z channel"""
        
    mtfile=np.loadtxt(filename)
    s=filename.find('.')
    if magtype=='lp':
        if filename[s+1:s+2]=='BY':
            convmtfile=convertlpB(mtfile,gain)
        elif filename[s+1:s+2]=='BX':
            convmtfile=convertlpB(mtfile,gain)
        elif filename[s+1:s+2]=='BZ':
            convmtfile=convertlpB(mtfile,gain,zadj=zadj)
    elif magtype=='bb':
        convmtfile=mtfile
    elif filename[s+1:s+3]=='EX':
        convmtfile=convertE(mtfile,gain,egain,dlength[0])
    elif filename[s+1:s+3]=='EY':
        convmtfile=convertE(mtfile,gain,egain,dlength[1])
    return convmtfile
    
    
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
    
    
    
def combineFewFiles(dirpath,station,starttime,endtime,cachelength,
                    complst=['BX','BY','BZ','EX','EY'],d=0,fdict=None):
    """
    combineFewFiles(dirpath,station,starttime,endtime,cachelength,
    complst=['BX','BY','BZ','EX','EY'],d=0)
    Will combine files in a directory path (dirpath) that have a given start and 
    end time in the form of (HHMMSS).  It looks for files at cachelength then 
    decimates the data by d. 
    Input:
        dirpath = directory path to folder where files to combine reside
        station = station name
        starttime = start time as 6 character string hhmmss
        endtime = end time as 6 character string hhmmss
        cachelength = cacherate of each file as 6 character string hhmmss
        complst = components to combine
        d = decimation factor, if no decimation enter as 0 or 1
        
    Output:
        cfilelst = list of combined files with full path for each component
                   files are saved as dirpath\CombHHtoHHd(d)\stationHHtoHH.comp
        fileslst = list of files that were combined including length of each 
                   file.
    """
    
    #set variables
    hbegin=starttime
    hend=endtime
    sbegin=hbegin[0:2]
    send=hend[0:2]
    cachesize=cachelength
    complst=complst
    station=station

    #compute total time
    tottime=int(hend)-int(hbegin)
    count=tottime/int(cachesize)
    
    #set pattern to search for to combine files
    patclst=[]
    for ii in range(count):
        pat=str(int(hbegin)+ii*int(cachesize))
        if len(pat)==1:
            patc='00000'+pat
        elif len(pat)==2:
            patc='0000'+pat
        elif len(pat)==3:
            patc='000'+pat
        elif len(pat)==4:
            patc='00'+pat
        elif len(pat)==5:
            patc='0'+pat
        elif len(pat)==6:
            patc=pat
        patclst.append(patc)
    
    
    #create combpath
    cfilenlst=[]
    if d==1 or d==0:
        dstr=''
    else:
        dstr='d'+str(d)
    #make a copy of the list cause can be a global variable that python can 
    #change, not good, should think about tuples perhaps
    complstt=list(complst)  
    combpath=os.path.join(dirpath,'Comb'+sbegin+'to'+send+dstr)
    returnlst=[os.path.join(combpath,station+sbegin+'to'+send+dstr+'.'+comp) 
                                                        for comp in complstt]
    if os.path.exists(combpath):
        for fn in os.listdir(combpath):
            for comp in complstt:
                if fnmatch.fnmatch(fn,'*.'+comp):
                    cfilenlst.append(os.path.join(combpath,fn))
                    complstt.remove(comp)
                    print 'File already combined: '+os.path.join(combpath,fn)
      
    #check to see if all components have been combined already
    if len(complstt)==0:
        return returnlst,['File already combined']
   
   #if not combine them
    else:
        #-----------------if no decimation copy files into one------------------
        if d==0 or d==1:
            if not os.path.exists(combpath):
                os.mkdir(combpath)
            fileslst=[]
            cfilenlst=[]
            for comp in complstt:
                fileextlst=[]
                cfilen=os.path.join(combpath,station+sbegin+'to'+
                                    send+'.'+comp)
                ##if the file already exists pass
                if os.path.isfile(cfilen):
                    cfilenlst.append(cfilen)
                    fileslst.append('Already combined')
                    print 'File already combined '+cfilen
                    
                #if there is a filter dictionary filter the data
                if fdict!=None:
                    print 'Applying Adaptive Notch Filter'
                    cfilenlst.append(cfilen)
                    combinemtfile=np.empty(0)
                    cfile=open(cfilen,'a')
                    for cc,ll in enumerate(range(count)):
                        ext='*'+patclst[ll]+'.'+comp
                        for filename in os.listdir(dirpath):
                            if fnmatch.fnmatch(filename,ext):
                                sfile=adaptiveNotchFilter(
                                        np.loadtxt(
                                        os.path.join(dirpath,filename)),
                                        **fdict)[0]
                                ls=len(sfile)
                                ##this is slow at the moment figure out a better
                                ##way to write these files
                                cfile.writelines(['%.9g' % ii+'\n' 
                                                    for ii in sfile])
#                                if cc==0:
#                                    combinemtfile=sfile
#                                else:
#                                    combinemtfile=np.hstack(sfile)
                                fileextlst.append(['File: '+ filename+' -> Length: '+str(ls)])
                                print 'Opened file: ',filename,' -> Length= ',ls
                    cfile.close()
                    print comp+' Finish Time: ',time.ctime()
#                    np.savetxt(cfilen,combinemtfile,fmt='%.9g')
                    print 'Combined file to: ',cfilen
                    fileslst.append(fileextlst)
                #if no decimation and filtering needed combine files quickly
                else:
                    cfilenlst.append(cfilen)
                    cfile=file(cfilen,'w')
                    fileextlst=[]
                    for ll in range(count):
                        ext='*'+patclst[ll]+'.'+comp
                        for filename in os.listdir(dirpath):
                            if fnmatch.fnmatch(filename,ext):
                                sfile=open(dirpath+os.sep+filename,'r')
                                sfn=sfile.readlines()
                                #make sure file position is 0
                                if sfile.tell()!=0: 
                                    sfile.seek(0)
                                shutil.copyfileobj(sfile, cfile,-1)
                                fileextlst.append(['File: '+filename+ \
                                                    ' -> Length: '+str((len(sfn)))])
                                print 'Opened file: ',filename,'Length= ',len(sfn)
                    print 'Combined file to: ',cfilen
                    fileslst.append(fileextlst)
                    cfile.close()
            return returnlst,fileslst
        #if decimation
        elif d!=0:
            if not os.path.exists(combpath):
                os.mkdir(combpath)
            fileslst=[]
            cfilenlst=[]
            for comp in complstt:
                fileextlst=[]
                cfilen=os.path.join(combpath,station+sbegin+'to'+
                                    send+dstr+'.'+comp)
                if os.path.isfile(cfilen):
                    cfilenlst.append(cfilen)
                    fileslst.append('Already combined')
                    print 'File already combined '+cfilen
                #if there is a filter dictionary filter the data
                if fdict!=None:
                    cfilenlst.append(cfilen)
                    countd=1
                    combinemtfile=np.empty(0)
                    for ll in range(count):
                        ext='*'+patclst[ll]+'.'+comp
                        for filename in os.listdir(dirpath):
                            if fnmatch.fnmatch(filename,ext):
                                sfile=adaptiveNotchFilter(decimatef(np.loadtxt(
                                                os.path.join(dirpath,filename))
                                                ,d),**fdict)[0]
                                ls=len(sfile)
                                if countd==1:
                                    combinemtfile=sfile
                                else:
                                    combinemtfile=np.hstack(sfile)
                                fileextlst.append(['File: '+ filename+' -> Length: '+str(ls)])
                                print 'Opened file: ',filename,' -> Length= ',ls
                                countd+=1
                    print comp+' Finish Time: ',time.ctime()
                    np.savetxt(cfilen,combinemtfile,fmt='%.9g')
                    print 'Combined file to: ',cfilen
                    fileslst.append(fileextlst)
                else:
                    cfilenlst.append(cfilen)
                    countd=1
                    combinemtfile=np.empty(0)
                    for ll in range(count):
                        ext='*'+patclst[ll]+'.'+comp
                        for filename in os.listdir(dirpath):
                            if fnmatch.fnmatch(filename,ext):
                                sfile=np.loadtxt(dirpath+os.sep+filename)
                                ls=len(sfile)
                                if countd==1:
                                    combinemtfile=decimatef(sfile,d)
                                else:
                                    combinemtfile=np.hstack((combinemtfile,
                                                             decimatef(sfile,
                                                                       d)))
                                fileextlst.append(['File: '+ filename+' -> Length: '+str(ls)])
                                print 'Opened file: ',filename,' -> Length= ',ls
                                countd+=1
                    print comp+' Finish Time: ',time.ctime()
                    np.savetxt(cfilen,combinemtfile,fmt='%.9g')
                    print 'Combined file to: ',cfilen
                    fileslst.append(fileextlst)
            return returnlst,fileslst
        
def combineFiles(dirpath,station,daylst,cacherate,
                 complst=['EX','EY','BX','BY','BZ'],dec=0,fdict=None):
    """
    combineFiles(dirpath,station,daydict,complst=['EX','EY','BX','BY','BZ'])
    will combine files from different days into one file.
    
    Inputs:
        dirpath = directory path where station files are
        station = name of station
        daylst = list of information to combine files with a dictionary for 
                each day with keys:
                    day = utc day (thee character string)
                    start = start time in ('hhmmss')
                    stop = end time in('hhmmss')
                    filt = 'Filtered' if time series was filtered
        cacherate = length of file ('hhmmss')
        dec = decimation factor (0 for none)
        complst = list of components to combine
        fdict = dictionary for adaptiveNotchFilter
    
    Outputs:
        cfilelst = list of combined files with full path as 
                dirpath/station/Combdaylst[0]todaylst[-1] 
    """
    
    complstd=list(complst)
    
    sday=daylst[0]['day']
    eday=daylst[-1]['day']
    
    if dec==0 or dec==1:    
        decstr=''
        dec=0
    else:
        decstr='d'+str(dec)
    try:
        daylst[0]['filt']
        savepath=os.path.join(dirpath,station,'Comb'+sday+'to'+eday+decstr+'Filt')
    except KeyError:
        savepath=os.path.join(dirpath,station,'Comb'+sday+'to'+eday+decstr)
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    
    #check to see if the files are already combined    
    cfilelst=[os.path.join(savepath,station+sday+'to'+eday+decstr+
             '.'+comp) for comp in complstd] 
    #remove any components that are already combined
    for cfile in cfilelst:
        if os.path.isfile(cfile)==True:
            complstd.remove(cfile[-2:])
   
    #if there are components that are not already combined combine them
    if len(complstd)>0:
        #make a dictionary for combining later        
        cfiledict={}
        for comp in complstd:
            cfiledict[comp]=[]
        #for each day combine files
        for dd,day in enumerate(daylst):
            #make a directory path            
            try:
                cdirpath=os.path.join(dirpath,station,day['day'],day['filt'])
            except KeyError:
                cdirpath=os.path.join(dirpath,station,day['day'])
            
            #combine files using function combineFewFiles
            clst,flst=combineFewFiles(cdirpath,station,day['start'],day['stop'],
                                      cacherate,complst=complstd,d=dec,
                                      fdict=fdict)

            #put files into a list for each component for easy combining
            for comp in complstd:
                ext='*.'+comp
                for fnc in clst:
                    if fnmatch.fnmatch(fnc,ext):
                        cfiledict[comp].append(fnc) 
        #combine day files into a days file
        fileslst=[]
        #ext='*'+patc+'.'+complst[jj]
        for comp in complstd:
            cfilen=os.path.join(savepath,station+sday+'to'+eday+decstr+'.'+comp)
            cfile=file(cfilen,'w')
            fileextlst=[]
            for fn in cfiledict[comp]:
                sfile=open(fn,'r')
                sfn=sfile.readlines()
                if sfile.tell()!=0: #make sure file position is 0
                    sfile.seek(0)
                shutil.copyfileobj(sfile, cfile,-1)
                fileextlst.append(['File: '+fn+' -> Length: '+
                                    str((len(sfn)))])
                print 'Opened file: ',fn,'Length= ',len(sfn)
                
            print 'Combined file to: ',cfilen
            fileslst.append(fileextlst)
        #return the filenames of the combined files and the files used to 
        #combine
        return cfilelst,fileslst
    #if the files are already combined return the filenames
    else:
        print 'Files already combined'
        return cfilelst,[]
    
def adaptiveNotchFilter(bx,df=100,notches=[50,100],notchradius=.5,freqrad=.9,rp=.1):
    """
    adaptiveNotchFilter(bx,df,notches=[50,100],notchradius=.3,freqrad=.9)
    will apply a notch filter to the array bx by finding the nearest peak around
    the supplied notch locations.  The filter is a zero-phase Chebyshev type 1 
    bandstop filter with minimal ripples.
    
    Inputs:
        bx = array to filter
        df = sampling frequency
        notches = list of notches to filter
        notchradius = radius of the notch in frequency domain
        freqrad = radius for searching for peak about notch
        rp = ripple of Chebyshev type 1 filter, lower numbers means less ripples

    Outputs:
        bx = filtered array
        filtlst = location of notches and power difference
    """
    bx=np.array(bx)
    
    if type(notches)!=list:
        notches=[notches]
    
    df=float(df)         #make sure df is a float
    dt=1./df             #sampling rate
    n=len(bx)            #length of array
    dfn=df/n             #frequency step
    dfnn=freqrad/dfn          #radius of frequency search
    fn=notchradius               #filter radius

    BX=np.fft.fft(bx)
    freq=np.fft.fftfreq(n,dt)
    
    filtlst=[]
    for notch in notches:
        fspot=int(round(notch/dfn))
        nspot=np.where(abs(BX)==max(abs(BX[fspot-dfnn:fspot+dfnn])))[0][0]
        dbstop=np.log10(abs(BX[nspot])-abs(BX).mean())
        if np.nan_to_num(dbstop)==0.0 or dbstop<3:
            filtlst.append('No need to filter \n')
            pass
        else:
            filtlst.append([freq[nspot],dbstop])
            ws=2*np.array([freq[nspot]-fn,freq[nspot]+fn])/df
            wp=2*np.array([freq[nspot]-2*fn,freq[nspot]+2*fn])/df
            ord,wn=sps.cheb1ord(wp,ws,1,dbstop)
            b,a=sps.cheby1(1,.5,wn,btype='bandstop')
            bx=sps.filtfilt(b,a,bx)
    
    return bx,filtlst

        

def removePeriodicNoise(filename,dt,noiseperiods,cacherate=10,save='n'):
    """
    removePeriodicNoise will take average a window of length noise period and 
    average the signal for as many windows that can fit within the data.  This
    averaged window is convolved with a series of delta functions at each window
    location to create a noise time series. This is then subtracted from the 
    data to get a 'noise free' time series.
    
    Inputs:
        filename = name of file to have periodic noise removed from
                    can be an array
        dt = time sample rate (s)
        noiseperiods = a list of estimated periods with a range of values to 
                        look around [[noiseperiod1,df1]...] where df1 is a 
                        fraction value find the peak about noiseperiod1 must be
                        less than 1.
        cacherate = length of file (min)
    
    Outputs:
        bxnf = filtered time series
        pn = periodic noise
        fitlst = list of peaks found
        
    """
    
    if type(noiseperiods)!=list:
        noiseperiods=[noiseperiods]
    
    dt=float(dt)    
    #sampling frequency    
    df=1./dt
    #length of array
    T=int(cacherate*60*df)
    #frequency step
    dfn=df/T
    #make a frequency array that describes BX
    pfreq=np.fft.fftfreq(int(T),dt) 
    
    filtlst=[]
    pnlst=[]
    for kk,nperiod in enumerate(noiseperiods):
        #if the nperiod is the first one load file or make an array of the input
        if kk==0:
            #load file
            if type(filename) is str:
                bx=np.loadtxt(filename)
                m=len(bx)
            else:
                bx=np.array(filename)
                m=len(bx)
        #else copy the already filtered array
        else:
            bx=bxnf.copy()
            m=len(bx)
        #get length of array
        T=len(bx)
        
        #get noise period in points along frequency axis
        nperiodnn=round((1./nperiod[0])/dfn)
        
        #get region to look around to find exact peak
        try:
            dfnn=nperiodnn*nperiod[1]
        except IndexError:
            dfnn=.2*nperiodnn
        
        
        #comput FFT of input to find peak value
        BX=np.fft.fft(bx)
        #if dfnn is not 0 then look for max with in region nperiod+-dfnn
        if dfnn!=0:
            nspot=np.where(abs(BX)==
                            max(abs(BX[nperiodnn-dfnn:nperiodnn+dfnn])))[0][0]
        else:
            nspot=nperiodnn
        #output the peak frequency found
        filtlst.append('Found peak at : '+str(pfreq[nspot])+' Hz \n')
        
        #make nperiod the peak period in data points
        nperiod=(1./pfreq[nspot])/dt

        #create list of time instances for windowing
        #nlst=np.arange(start=nperiod,stop=T-nperiod,step=nperiod,dtype='int')
        nlst=np.arange(start=0,stop=m,step=nperiod,dtype='int')
        
        #convolve a series of delta functions with average of periodic window
        dlst=np.zeros(T)               #delta function list
        dlst[0]=1
        winlst=np.zeros((len(nlst),int(nperiod)))
        for nn,ii in enumerate(nlst):
            if T-ii<nperiod:
                dlst[ii]=1
            else:
                winlst[nn]=bx[ii:ii+int(nperiod)]
                dlst[ii]=1
                
        #compute median window to remove any influence of outliers
        medwin=np.median(winlst,axis=0)
        
        #make a time series by convolving
        pn=np.convolve(medwin,dlst)[0:T]

        #remove noise from data
        bxnf=(bx-pn)
        
        pnlst.append(pn)
    if len(pnlst)>1:
        pn=np.sum(pnlst,axis=0)
    else:
        pn=np.array(pn)
    if save=='y':
        savepath=os.path.join(os.path.dirname(filename),'Filtered')
        if not os.path.exists(savepath):
            os.mkdir(savepath)
        #savepathCN=os.path.join(savepath,'CN')
        np.savetxt(os.path.join(savepath,filename),bxnf,fmt='%.7g')
        #np.savetxt(os.path.join(savepathCN,filename[0:5]+'CN'+filename[5:]),pn)
    else:
        return bxnf,pn,filtlst

def imp2resphase(z,freq,zvar=None,ffactor=1):
    """imp2resphase(z,zvar,freq) will convert impedances 
    z[4,3,len(freq)] to resistivities (ohm-m) and phases (deg) as well as the 
    errors of each.  Note the phase is calculated using arctan2 putting the 
    phase in the correct quadrant.  The xy component is placed into the positive
    quadrant by adding 180 deg. 
    INPUTS: 
        z=[[zxxreal,zxximag,zxxvar],...]], shape=(4,3,len(freq) OR 
        z=[[[zxxr+1j*zxxi,zxyr+1j*zxyi],[zyx,zyy]]], shape=(len(freq),2,2) 
        freq=frequency array 
        zvar=[[[zxxvar,zyyvar],[zyxvar,zyyvar]]], shape=(len(freq),2,2)  
        ffactor=fudge factor to make sure resistivities are the right order
    OUTPUTS: 
        [[resxx,resxxvar],[resxy,resxyvar]...],[[phasexx,phasexxvar]...]"""
    
    #if in first format (4,3,n) change to the second format with shape (n,2,2)
    s=z.shape
    if s[1]==3:
        if zvar==None:
            znew=[]
            zvar=[]
            for jj in range(s[2]):
                zxx=z[0,0,jj]+1j*z[0,1,jj]
                zxy=z[1,0,jj]+1j*z[1,1,jj]
                zyx=z[2,0,jj]+1j*z[2,1,jj]
                zyy=z[3,0,jj]+1j*z[3,1,jj]
                znew.append([[zxx,zxy],[zyx,zyy]])
                zvar.append([[z[0,2,jj],z[1,2,jj]],[z[2,2,jj],z[3,2,jj]]])
        z=np.array(znew)
        zvar=np.array(zvar)
    
    period=1/freq
    if period[0]>period[-1]:
        z=z[::-1,:,:]
        zvar=zvar[::-1,:,:]
        period=period[::-1]
        freq=freq[::-1]
        print 'Flipped arrays so frequency is decending.'
    
    if zvar==None:
        zvar=np.zeros_like(z.real)
    resxy=[]
    resxyerr=[]
    resyx=[]
    resyxerr=[]
    resxx=[]
    resxxerr=[]
    resyy=[]
    resyyerr=[]
    
    phasexy=[]
    phasexyerr=[]
    phaseyx=[]
    phaseyxerr=[]
    phasexx=[]
    phasexxerr=[]
    phaseyy=[]
    phaseyyerr=[]
    for jj in range(len(z)):
        wt=.2/(freq[jj])*float(ffactor)
        resxx.append(wt*abs(z[jj,0,0])**2)
        resxy.append(wt*abs(z[jj,0,1])**2)
        resyx.append(wt*abs(z[jj,1,0])**2)
        resyy.append(wt*abs(z[jj,1,1])**2)
        
        resxxerr.append(wt*(abs(z[jj,0,0])+zvar[jj,0,0])**2
                                                            -resxx[jj])
        resxyerr.append(wt*(abs(z[jj,0,1])+zvar[jj,0,1])**2
                                                            -resxy[jj])
        resyxerr.append(wt*(abs(z[jj,1,0])+zvar[jj,1,0])**2
                                                            -resyx[jj])
        resyyerr.append(wt*(abs(z[jj,1,1])+zvar[jj,1,1])**2
                                                            -resyy[jj])
        
        phasexx.append(np.arctan2(z[jj,0,0].imag,z[jj,0,0].real)\
                                                            *(180/np.pi))
        phasexy.append(np.arctan2(z[jj,0,1].imag,z[jj,0,1].real)\
                                                            *(180/np.pi))
        phaseyx.append(np.arctan2(z[jj,1,0].imag,z[jj,1,0].real)\
                                                            *(180/np.pi))
        phaseyy.append(np.arctan2(z[jj,1,1].imag,z[jj,1,1].real)\
                                                            *(180/np.pi))
        
        phasexxerr.append(np.arcsin(zvar[jj,0,0]/abs(z[jj,0,0]))\
                                                            *(180/np.pi))
        phasexyerr.append(np.arcsin(zvar[jj,0,1]/abs(z[jj,0,1]))\
                                                            *(180/np.pi))
        phaseyxerr.append(np.arcsin(zvar[jj,1,0]/abs(z[jj,1,0]))\
                                                            *(180/np.pi))
        phaseyyerr.append(np.arcsin(zvar[jj,1,1]/abs(z[jj,1,1]))\
                                                            *(180/np.pi))
    res=[[resxx,resxxerr],[resxy,resxyerr],[resyx,resyxerr],[resyy,resyyerr]]
    phase=[[phasexx,phasexxerr],[phasexy,phasexyerr],[phaseyx,phaseyxerr],
           [phaseyy,phaseyyerr]]
           
        
    return np.array(res),np.array(phase)

def sigfigs(numstr,digits=8,fmt='g'):
    """sigfigs(numstr,digits=8,fmt='g') will return a string with the proper 
    amount of significant digits for the input number, can be str or float."""
    numfmt='%.'+str(digits)+fmt
    if type(numstr) is float:
        numstr=str(numstr)
        
    if numstr.find('nan')>=0 or numstr.find('NaN')>=0 or numstr.find('+Inf')>=0:
        signum=numfmt % -6.66666666
    else:
        signum=numfmt % float(numstr)
    
    if len(signum)>digits:
        signum=signum[0:digits]
        
    return signum
    
def writeedi(z,zvar,freq,stationinfofile,station,edidir=None,rrstation=None,
             birrpdict=None):
    """writeedi2(z,zvar,freq,stationinfofile,station,edidir=None,rrstation=None,
    birrpdict=None) will write an .edi file for a station. 
    Returns full filepath."""
    
    #get station info from station info file
    statdict=getStationInfo(stationinfofile,station)
    if freq[0]<freq[-1]:
        z=z[:,:,::-1]
        print 'Flipped array so frequency is decending'
    nfreq=len(freq)

        
    #list of orientation components
    orilst=['HX','HY','EX','EY','RX','RY']
    
    if edidir==None:
        dirpath=os.getcwd()
    else:
        dirpath=edidir
    edipath=os.path.join(dirpath,station+'.edi')
    edifid=file(edipath,'w')
    
    #---------------------------------------------------------------------------
    #write header information
    edifid.write('>HEAD \n')
    edifid.write(tsp+'DATAID="'+station+'"'+'\n')
    
    #acquired by:
    try:
        edifid.write(tsp+'ACQBY="'+statdict['acqby']+'"'+'\n')
    except KeyError:
        edifid.write(tsp+'ACQBY="'+'Adelaide University'+'"'+'\n')
    #aqcuired date
    try:
        mdate=statdict['date'].split('/')
        mday=int(mdate[0])
        mmonth=int(mdate[1])
        if len(mdate[2])<3:
            myear=int('20'+mdate[2])
        else:
            myear=int(mdate[2])
        md=datetime.date(myear,mmonth,mday)
        edifid.write(tsp+'ACQDATE='+datetime.date.strftime(md,'%B %d, %Y')+'\n')
    except KeyError:
        edifid.write(tsp+'ACQDATE='+'\n')
    #date edi file written
    edifid.write(tsp+'FILEDATE='+datetime.date.strftime(datetime.date.today(),
                                                       '%B %d, %Y')+'\n')
    #survey location
    try: 
        edifid.write(tsp+'PROSPECT="'+statdict['location']+'"'+'\n')
    except KeyError:
        edifid.write(tsp+'PROSPECT=" "'+'\n')
    #station name
        edifid.write(tsp+'LOC="'+station+'"'+'\n')
    #latitude
    try:
        edifid.write(tsp+'LAT='+'%2.8g' % float(statdict['lat'])+'\n')
    except KeyError:
        edifid.write(tsp+'LAT= \n')
    #longitude
    try: 
        edifid.write(tsp+'LONG='+'%2.8g' % float(statdict['long'])+'\n')
    except KeyError:
        edifid.write(tsp+'LONG= \n')
    #elevation
    try:
        edifid.write(tsp+'ELEV='+statdict['elev']+'\n')
    except:
        edifid.write(tsp+'ELEV= \n')
    edifid.write('\n')
    #---------------------------------------------------------------------------
    #Write info block
    edifid.write('>INFO'+tsp+'MAX LINES=1000'+'\n')
    
    #survey parameters
    edifid.write(tsp+'Survey Parameters: \n')
    try:
        edifid.write(lsp+'Sampling Frequency (Hz): '+statdict['df']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Cache Rate (HHMMSS): '+statdict['cacherate']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Data Logger Gain: '+statdict['dlgain']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Interface Box Gain: '+statdict['egain']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Instrument Box no: '+statdict['box no']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Coil Numbers (BX,BY,BZ): '+statdict['coil no(bx,by)']
                                                                    +'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Data Logger: '+statdict['dlbox']
                                                                +'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Hard Drive no: '+statdict['harddrive']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Interface Box no: '+statdict['interfacebox']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Battery no: '+statdict['battery']
                        +' Starting Voltage: '+statdict['start volt']
                        +' End Voltage '+statdict['end volt']+'\n')
    except KeyError:
        pass
    
    try:
        edifid.write(lsp+'Other Notes: '+statdict['notes']+'\n')
    except KeyError:
        pass
    
    #BIRRP parameters
    edifid.write('   Transfer Functions Computed using BIRRP 5.1  \n')
    if birrpdict!=None:
        edifid.write(lsp+'Coil Calibration File='+birrpdict['bbfile']+'\n')
        edifid.write(lsp+'Interaction Level (ILEV)='+birrpdict['ilev']+'\n')
        edifid.write(lsp+'Number of outputs (NOUT)='+birrpdict['nout']+'\n')
        edifid.write(lsp+'Number of inputs (NINP)='+birrpdict['ninp']+'\n')
        edifid.write(lsp+'Slepian Filter order (TBW)='+birrpdict['nref']+'\n')
        edifid.write(lsp+'Max length of fft window (NFFT)='+birrpdict['nfft']+'\n')
        edifid.write(lsp+'Maximum number of fft sections (NSCTMAX)='
                                    +birrpdict['nsctmax']+'\n')
        edifid.write(lsp+'Small leverage point control (UIN)='+birrpdict['uin']
                                                                +'\n')
        edifid.write(lsp+'Large leverage point control (AINUIN)='
                                +birrpdict['ainuin']+'\n')
        edifid.write(lsp+'Electric coherence threshold (C2THRESHE)='
                                    +birrpdict['c2threshe']+'\n')
        edifid.write(lsp+'Z component (NZ)='+birrpdict['nz']+'\n')
        edifid.write(lsp+'Coherence threshold z channel (C2threse1)='
                                +birrpdict['c2threse1']+'\n')
        edifid.write(lsp+'Order of prewhitening filter (NAR)='+birrpdict['nar']
                                                                +'\n')
        edifid.write(lsp+'Electric channel rotation angles (THETAE)='
                                        +birrpdict['thetae']+'\n')
        edifid.write(lsp+'Magnetic channel rotation angles (THETAB)='
                                        +birrpdict['thetab']+'\n')
        edifid.write(lsp+'Final channel rotation angles (THETAF)='
                                        +birrpdict['thetaf']+'\n')
        
    #remote reference
    if rrstation!=None:
        rrdict=getStationInfo(stationinfofile,rrstation)
        edifid.write(lsp+'Remote Reference Station: '+rrstation+'\n')
        edifid.write(lsp+'Remote Reference Lat='
                                    +'%2.8g' % float(rrdict['lat'])+'\n')
        edifid.write(lsp+'Remote Reference Long='
                                    +'%2.8g' % float(rrdict['long'])+'\n')
        edifid.write(lsp+'Remote Reference Elev='+rrdict['elev']+'\n')
    
    edifid.write('\n')
    
    #---------------------------------------------------------------------------
    #write define measurement block
    edifid.write('>=DEFINEMEAS'+'\n'+'\n')
    edifid.write(tsp+'MAXCHAN=6'+'\n')
    edifid.write(tsp+'MAXRUN=999'+'\n')
    edifid.write(tsp+'MAXMEAS=99999'+'\n')
    edifid.write(tsp+'UNITS=M'+'\n')
    edifid.write(tsp+'REFTYPY=CART'+'\n')
    edifid.write(tsp+'REFLAT='+'%2.8g' % float(statdict['lat'])+'\n')
    edifid.write(tsp+'REFLONG='+'%2.8g' % float(statdict['long'])+'\n')
    edifid.write(tsp+'REFELEV='+statdict['elev']+'\n')
    edifid.write('\n'+'\n')
    
    edifid.write('>HMEAS ID=1001.001 CHTYPE=HX X=0 Y=0 AZM=0'+'\n')
    edifid.write('>HMEAS ID=1002.001 CHTYPE=HY X=0 Y=0 AZM=90'+'\n')
    edifid.write('>EMEAS ID=1003.001 CHTYPE=EX X=0 Y=0 X2='+statdict['ex']
                                                            +' Y2=0'+'\n')
    edifid.write('>EMEAS ID=1004.001 CHTYPE=EY X=0 Y=0 X2=0 Y2='+statdict['ey']
                                                            +'\n')
    edifid.write('>HMEAS ID=1005.001 CHTYPE=RX X=0 Y=0 AZM=0'+'\n')
    edifid.write('>HMEAS ID=1006.001 CHTYPE=RY X=0 Y=0 AZM=90'+'\n')
    edifid.write('\n')
    
    #---------------------------------------------------------------------------
    #write mtsect block
    edifid.write('>=MTSECT \n')
    edifid.write(tsp+'SECTID='+station+'\n')
    edifid.write(tsp+'NFREQ='+nfreq+'\n')
    edifid.write(tsp+orilst[0]+'=1001.001'+'\n')
    edifid.write(tsp+orilst[1]+'=1002.001'+'\n')
    edifid.write(tsp+orilst[2]+'=1003.001'+'\n')
    edifid.write(tsp+orilst[3]+'=1004.001'+'\n')
    edifid.write(tsp+orilst[4]+'=1005.001'+'\n')
    edifid.write(tsp+orilst[5]+'=1006.001'+'\n')
    edifid.write('\n')
    edifid.write('>!****FREQUENCIES****!'+'\n')
    if freq[0]<freq[-1]:
        order='INC'
    else:
        order='DEC'
    edifid.write('>FREQ'+tsp+'NFREQ='+nfreq+tsp+'ORDER='+order+tsp+'// '+
                                                                    nfreq+'\n')
    for kk in range(int(nfreq)):
        edifid.write(tsp+'%2.6f' % freq[kk])
        if np.remainder(float(kk)+1,5.)==0:
            edifid.write('\n')
    edifid.write('\n')
    edifid.write('>!****IMPEDANCES****!'+'\n')
    
    implst=[['ZXXR',0,0],['ZXXI',0,1],['ZXX.VAR',0,2],['ZXYR',1,0],['ZXYI',1,1],\
        ['ZXY.VAR',1,2],['ZYXR',2,0],['ZYXI',2,1], ['ZYX.VAR',2,2],\
        ['ZYYR',3,0],['ZYYI',3,1],['ZYY.VAR',3,2]]
    #write new impedances and variances
    for jj,imp in enumerate(implst):
        mm=imp[1]
        nn=imp[2]
        edifid.write('>'+imp[0]+' // '+nfreq+'\n')
        for kk in range(int(nfreq)):
            edifid.write(tsp+'%2.6e' % z[mm,nn,kk])
            if np.remainder(float(kk)+1,5.)==0:
                edifid.write('\n')
        edifid.write('\n')
    edifid.write('\n')
    
    #---------------------------------------------------------------------------
    #write tipper info
    
    edifid.write('>!****TIPPER****!'+'\n')
    tiplst=['TXR','TXI','TX.VAR','TY','TYI','TY.VAR']
    for tip in tiplst:
        edifid.write('>'+tip+' // '+nfreq+'\n')
        for kk in range(int(nfreq)):
            edifid.write(tsp+'%2.6e' % 0.0)
            if np.remainder(float(kk),5.)==np.remainder(float(nfreq),5.):
                edifid.write('\n')
        edifid.write('\n')
    edifid.write('\n')
    edifid.write('>END')
    edifid.close()
    
    return edipath
    
def rewriteedi(edifile,znew=None,zvarnew=None,freqnew=None,newfile='y',
               tipnew=None,tipvarnew=None,thetar=0,dr='n'):
    """
    rewriteedi(edifile) will rewrite an edifile say if it needs to be rotated 
    or distortion removed.

    Inputs:
        edifile = full path to edifile to be rewritten
        znew = impedance tensor if a new one has been created
        zvarnew = errors in impedance tensor if a new one has been created
        freqnew = new frequency list if one has been created
        newfile = 'y' for yes or 'n' for no if you want a new file to be
                 created.
        tipnew = new tipper array
        tipvarnew = new tipper error array
        thetar = rotation angle counter clockwise (N=0, E=-90)
        dr = 'n' if no distortion removal and 'y' for distortion removal
    
    Outputs:
        nedi = dirpath(edifile)+basename(edifile)+rw or dr if dr='y'
    
    """

    #get direcotry path make one if not there    
    dirpath=os.path.dirname(edifile)
    if newfile=='y':
        if dr=='y':
            drdirpath=os.path.join(dirpath,'DR')
        else:
            drdirpath=os.path.join(dirpath,'Rot{0:.0f}'.format(thetar))
        if not os.path.exists(drdirpath):
            os.mkdir(drdirpath)
            print 'Made directory: ',drdirpath
        else:
            pass
        #copy files into a common directory
        drpath=os.path.dirname(dirpath)
        count=1
        while count<10:
            drpath=os.path.dirname(drpath)
            for folder in os.listdir(drpath):
                if folder.find('EDIfiles')>=0:
                    edifolderyn='y'
                    drcopypath=os.path.join(drpath,'EDIfiles')
                    if os.path.exists(drcopypath):
                        if not os.path.exists(os.path.join(drcopypath,'DR')):
                            os.mkdir(os.path.join(drcopypath,'DR'))
                            print 'Made directory ',os.path.join(drcopypath,'DR')
                else:
                    edifolderyn='n'
                    drcopypath=os.path.dirname(drdirpath)
            count+=1
        
        #get new file name
        if dr=='y':
            newedifile=os.path.basename(edifile)[:-4]+'dr.edi'
        else:
            newedifile=os.path.basename(edifile)[:-4]+'rw.edi'
        #open a new edifile
        newedifid=open(os.path.join(drdirpath,newedifile),'w')
        
    else:
        newedifile=os.path.basename(edifile)[:-4]+'c.edi'
        newedifid=open(os.path.join(dirpath,newedifile),'w')
        edifolderyn='n'
        
    if znew==None:
        edidict=readedi(edifile)
        #rotate data if desired
        if thetar!=0:
            znew=np.zeros_like(edidict['z'])
            zvarnew=np.zeros_like(edidict['z'])
            #convert thetar into radians
            if abs(thetar)>2*np.pi:
                thetar=thetar*(np.pi/180)
            #make rotation matrix
            rotmatrix=np.array([[np.cos(thetar), np.sin(thetar)],
                         [-np.sin(thetar), np.cos(thetar)]])
            #rotate the data
            for rr in range(len(edidict['frequency'])):
                znew[rr]=np.dot(rotmatrix,np.dot(edidict['z'][rr],rotmatrix.T))
                zvarnew[rr]=np.dot(rotmatrix,np.dot(edidict['zvar'][rr],
                                        rotmatrix.T))
        else:
            znew=edidict['z']
            zvarnew=edidict['zvar']
    
    #open existing edi file
    edifid=open(edifile,'r')
    #read existing edi file and find where impedances start and number of freq
    edilines=edifid.readlines()
    spot=[]
    for ii,line in enumerate(edilines):
        if line.find('>FREQ')>=0:
            spot.append(ii)
            linelst=line.split('//')
            nfreq=int(linelst[-1].rstrip())
        if line.find('IMPEDANCES')>=0:
            spot.append(ii)
        if line.find('TIPPER')>=0:
            spot.append(ii)
    
    #if there is a new frequency list
    if freqnew!=None:
        #get number of frequencies    
        nfreq=len(freqnew)
        for ll,line in enumerate(edilines[0:spot[0]]):
            if line.find('NFREQ=')>=0:
                newedifid.write(tsp+'NFREQ='+str(int(nfreq))+'\n')
            else:
                newedifid.write(line)
        
        #get order of frequencies        
        if freqnew[0]<freqnew[-1]:
            order='INC'
        else:
            order='DEC'
        
        #write frequency header
        newedifid.write('>FREQ'+tsp+'NFREQ='+str(int(nfreq))+tsp+'ORDER='+order+
                        tsp+'// '+str(int(nfreq))+'\n')
        for kk in range(int(nfreq)):
            newedifid.write(tsp+'%2.8f' % freqnew[kk])
            if np.remainder(float(kk)+1,5.)==0:
                newedifid.write('\n')
        newedifid.write('\n')
        newedifid.write('>!****IMPEDANCES****!'+' Rotated {0:.0f} counterclockwise\n'.format(thetar*180/np.pi))
    #from the start of file to where impedances start write from existing edifile
    else:
        for line in edilines[0:spot[1]]:
            newedifid.write(line)
    newedifid.write('>!****IMPEDANCES****!'+' Rotated {0:.0f} counterclockwise\n'.format(thetar*180/np.pi))
    
    #create an implst to make code simpler
    implst=[['ZXXR',0,0],['ZXXI',0,0],['ZXX.VAR',0,0],['ZXYR',0,1],['ZXYI',0,1],\
        ['ZXY.VAR',0,1],['ZYXR',1,0],['ZYXI',1,0], ['ZYX.VAR',1,0],\
        ['ZYYR',1,1],['ZYYI',1,1],['ZYY.VAR',1,1]]
    #write new impedances and variances
    for jj,imp in enumerate(implst):
        mm=imp[1]
        nn=imp[2]
        newedifid.write('>'+imp[0]+' // '+str(int(nfreq))+'\n')
        if imp[0].find('.')==-1:
            if imp[0].find('R')>=0:
                for kk in range(int(nfreq)):
                    newedifid.write(tsp+'%+.7E' % znew.real[kk,mm,nn])
                    if kk>0:
                        if np.remainder(float(kk)+1,5.)==0:
                            newedifid.write('\n')
                    else:
                        pass
            elif imp[0].find('I')>=0:
                for kk in range(int(nfreq)):
                    newedifid.write(tsp+'%+.7E' % znew.imag[kk,mm,nn])
                    if kk>0:
                        if np.remainder(float(kk)+1,5.)==0:
                            newedifid.write('\n')
                    else:
                        pass
        else:
            for kk in range(int(nfreq)):
                newedifid.write(tsp+'%+.7E' % zvarnew[kk,mm,nn])
                if kk>0:
                    if np.remainder(float(kk)+1,5.)==0:
                        newedifid.write('\n')
                else:
                    pass
        newedifid.write('\n')
    
    if tipnew!=None:
        newedifid.write('>!****TIPPER****!'+'\n')
        tiplst=[['TXR',0],['TXI',0],['TX.VAR',0],['TYR',1],['TYI',1],
                ['TY.VAR',1]]
        if len(tipnew)==0:
            tipnew=np.zeros((2,3,float(nfreq)))
        else:
            tipnew=np.array(tipnew)
        for jj,tip in enumerate(tiplst):
            mm=tip[1]
            newedifid.write('>'+tip[0]+' // '+str(int(nfreq))+'\n')
            if tip[0].find('.')==-1:
                if tip[0].find('R')>=0:
                    for kk in range(int(nfreq)):
                        newedifid.write(tsp+'%+.7E' % tipnew.real[kk,mm])
                        if kk>0:
                            if np.remainder(float(kk)+1,5.)==0:
                                newedifid.write('\n')
                        else:
                            pass
                elif tip[0].find('I')>=0:
                    for kk in range(int(nfreq)):
                        newedifid.write(tsp+'%+.7E' % tipnew.imag[kk,mm])
                        if kk>0:
                            if np.remainder(float(kk)+1,5.)==0:
                                newedifid.write('\n')
                        else:
                            pass
            else:
                for kk in range(int(nfreq)):
                    newedifid.write(tsp+'%+.7E' % tipvarnew.real[kk,mm])
                    if kk>0:
                        if np.remainder(float(kk)+1,5.)==0:
                            newedifid.write('\n')
                    else:
                        pass
            newedifid.write('\n')
        newedifid.write('>END')
    else:
        #print nlinesleft
        for line in edilines[spot[2]:]:
            newedifid.write(line)
    
    edifid.close()
    newedifid.close()
    #copy file to a common folder
    if edifolderyn=='y':
        if newfile=='y':
            if dr=='y':
                shutil.copy(os.path.join(drdirpath,newedifile),
                            os.path.join(drcopypath,'DR',newedifile))
            else:
                shutil.copy(os.path.join(drdirpath,newedifile),
                            os.path.join(drcopypath,'RW',newedifile))
        else:
            shutil.copy(os.path.join(dirpath,newedifile),
                    os.path.join(drcopypath,newedifile))
    else:
        pass
    if newfile=='y':
        return os.path.join(drdirpath,newedifile)
    else:
        return os.path.join(dirpath,newedifile)
    
def getnum(numberlst):
    """get number from string list and put it into an array"""
    
    nlst=[]
    for jj in range(len(numberlst)):
        numlst=numberlst[jj].rstrip().split()
        for ii in range(len(numlst)):
            num=numlst[ii]
            if len(num)>=2:
                try:
                    nlst.append(float(num))
                except ValueError:
                    nlst.append(-6.666)
    return np.array(nlst)
    
def readedi(filename):
    """readedi(edifile) will read in an edi file written in a format given by
    format given by http://www.dias.ie/mtnet/docs/ediformat.txt. 
    Returns: lat,lon,frequency,Z[zreal+i*zimag],Zvar,tipper,tippervar (if 
    applicable)
    
    Input: 
        filename = full path to edifile
        
    Output:
        edidict = dictionary with keys:
            station = station name
            lat = latitude 
            lon = longitude
            frequency = frequency array
            z = impedance tensor as (nfreq,2,2)
            zvar = variance of impedance tensor as (nfreq,2,2)
            tipper = tipper as (nfreq,2,1)
            tippervar = tipper variance as (nfreq,2,1)
    
    """

    edifid=file(filename)
    edilines=edifid.readlines()
    existlst=[]
    tfind=0
    for ii in range(len(edilines)):
        strline=str(edilines[ii])
        if strline.find('HEAD')>0:
            existlst.append(['head',ii])
        if strline.find('INFO',0,6)>0:
            existlst.append(['info',ii])
        if strline.find('DEFINE')>0:
            existlst.append(['define',ii])
        if strline.find('MTSECT')>0:
            existlst.append(['mtsect',ii])
        if strline.find('FREQUENCIES')>0:
            existlst.append(['frequencies',ii])
        if strline.find('ROTATION')>0:
            existlst.append(['rotation',ii])
        if strline.find('IMPEDANCES')>0:
            existlst.append(['impedances',ii])
        if strline.find('TIPPER')>0:
            existlst.append(['tipper',ii])
        if strline.find('END')>0:
            existlst.append(['end',ii])
            
#    print existlst

    #get things into lists
    for ii in range(len(existlst)):
        if existlst[ii][0].find('head')>=0:
            header=edilines[existlst[ii][1]:existlst[ii+1][1]]
            for ll in range(len(header)):
                if header[ll].find('LAT')>=0:
                    latfind='y'
                    latstr=header[ll].strip()
                    latlst=latstr.split('=')[1]
                    #change from hh:mm:ss to decimal deg
                    if latlst.find(':')>=0:
                        latlst=latlst.split(':')
                        latstr=str(int(latlst[0]))+'.'+str(float(latlst[1])/60
                            +float(latlst[2])/3600)[2:]
                        lat=float(latstr)
                    else:
                        lat=float(latlst)
                if header[ll].find('LONG')>=0:
                    lonfind='y'
                    lonstr=header[ll].strip()
                    lonlst=lonstr.split('=')[1]
                    #change from hh:mm:ss to decimal deg
                    if lonlst.find(':')>=0:
                        lonlst=lonlst.split(':')
                        lonstr=str(int(lonlst[0]))+'.'+str(float(lonlst[1])/60
                            +float(lonlst[2])/3600)[2:]
                        lon=float(lonstr)
                    else:
                        lon=float(lonlst)
                if header[ll].find("LOC")>=0:
                    locstr=header[ll].strip()
                    locstr=locstr.split('=')[1]
                    station=locstr.replace('"','')
                elif header[ll].find("DATAID")>=0:
                    locstr=header[ll].strip()
                    locstr=locstr.split('=')[1]
                    station=locstr.replace('"','')
                else:
                    station=os.path.basename(filename)
                    
            #print ii, ' Got Header'
        elif existlst[ii][0].find('info')>=0:
            #print ii, ' Got info'
            info=edilines[existlst[ii][1]:existlst[ii+1][1]]
        elif existlst[ii][0].find('define')>=0:
            definem=edilines[existlst[ii][1]:existlst[ii+1][1]]
            try:
                latfind
                lonfind
            except NameError:
                for ll in range(len(definem)):
                    if definem[ll].find('LAT')>=0:
                        latstr=definem[ll].strip()
                        latlst=latstr.split('=')[1]
                        #change from hh:mm:ss to decimal deg
                        if latlst.find(':')>=0:
                            latlst=latlst.split(':')
                            latstr=str(int(latlst[0]))+'.'+str(float(latlst[1])/60
                                +float(latlst[2])/3600)[2:]
                            lat=float(latstr)
                        else:
                            lat=float(latlst)
                    if definem[ll].find('LONG')>=0:
                        lonstr=definem[ll].strip()
                        lonlst=lonstr.split('=')[1]
                        #change from hh:mm:ss to decimal deg
                        if lonlst.find(':')>=0:
                            lonlst=lonlst.split(':')
                            lonstr=str(int(lonlst[0]))+'.'+str(float(lonlst[1])/60
                                +float(lonlst[2])/3600)[2:]
                            lon=float(lonstr)
                        else:
                            lon=float(lonlst)
            
            #print ii, ' Got define measurement'
            
        elif existlst[ii][0].find('mtsect')>=0:
            #print ii, ' Got mtsect'
            mtsect=edilines[existlst[ii][1]:existlst[ii+1][1]]
        elif existlst[ii][0].find('frequencies')>=0:
            #print ii, ' Got frequencies'
            frequencylst=edilines[existlst[ii][1]+2:existlst[ii+1][1]]
            freq=getnum(frequencylst)
        elif existlst[ii][0].find('rotation')>=0:
            #print ii, ' Got rotations'
            rotlst=edilines[existlst[ii][1]+2:existlst[ii+1][1]]
            rot=getnum(rotlst)
        elif existlst[ii][0].find('impedances')>=0:
            
            impedancelst=edilines[existlst[ii][1]+1:existlst[ii+1][1]]
            imfindlst=[]
            for jj in range(len(impedancelst)):
                if impedancelst[jj].find('>')>=0:
                    imfindlst.append(jj)
            zxxr=getnum(impedancelst[imfindlst[0]+1:imfindlst[1]])
            zxxi=getnum(impedancelst[imfindlst[1]+1:imfindlst[2]])
            zxxv=getnum(impedancelst[imfindlst[2]+1:imfindlst[3]])
            zxyr=getnum(impedancelst[imfindlst[3]+1:imfindlst[4]])
            zxyi=getnum(impedancelst[imfindlst[4]+1:imfindlst[5]])
            zxyv=getnum(impedancelst[imfindlst[5]+1:imfindlst[6]])
            zyxr=getnum(impedancelst[imfindlst[6]+1:imfindlst[7]])
            zyxi=getnum(impedancelst[imfindlst[7]+1:imfindlst[8]])
            zyxv=getnum(impedancelst[imfindlst[8]+1:imfindlst[9]])
            zyyr=getnum(impedancelst[imfindlst[9]+1:imfindlst[10]])
            zyyi=getnum(impedancelst[imfindlst[10]+1:imfindlst[11]])
            zyyv=getnum(impedancelst[imfindlst[11]+1:imfindlst[11]+imfindlst[1]-imfindlst[0]])
            #print ii, ' Got impedances'
        elif existlst[ii][0].find('tipper')>=0:
            tflst=['trot','txr','txi','txvar','tx.var','tyr','tyi','tyvar',
                   'ty.var']
            tfind=1
            tipperlst=edilines[existlst[ii][1]+1:existlst[ii+1][1]]
            tipfindlst=[]
            for nn in range(len(tipperlst)):
                for tt in tflst:
                    if tipperlst[nn].find(tt.upper())>0:
                        tipfindlst.append(nn)
#                if tipperlst[nn].find('>')>=0:
#                    tipfindlst.append(nn)
#            print tipfindlst,len(tipfindlst)
            if len(tipfindlst)==6:
                txr=getnum(tipperlst[tipfindlst[0]+1:tipfindlst[1]])
                txi=getnum(tipperlst[tipfindlst[1]+1:tipfindlst[2]])
                txv=getnum(tipperlst[tipfindlst[2]+1:tipfindlst[3]])
                tyr=getnum(tipperlst[tipfindlst[3]+1:tipfindlst[4]])
                tyi=getnum(tipperlst[tipfindlst[4]+1:tipfindlst[5]])
                tyv=getnum(tipperlst[tipfindlst[5]+1:tipfindlst[5]+
                                                tipfindlst[1]-tipfindlst[0]])
            elif len(tipfindlst)==7:
                trot=getnum(tipperlst[tipfindlst[0]+1:tipfindlst[1]])
                txr=getnum(tipperlst[tipfindlst[1]+1:tipfindlst[2]])
                txi=getnum(tipperlst[tipfindlst[2]+1:tipfindlst[3]])
                txv=getnum(tipperlst[tipfindlst[3]+1:tipfindlst[4]])
                tyr=getnum(tipperlst[tipfindlst[4]+1:tipfindlst[5]])
                tyi=getnum(tipperlst[tipfindlst[5]+1:tipfindlst[6]])
                tyv=getnum(tipperlst[tipfindlst[6]+1:tipfindlst[6]+
                                                tipfindlst[1]-tipfindlst[0]])
            #print ii, ' Got tipper'
    #put things into a dictionary to return values
    edidict={}
    edidict['station']=station
    edidict['lat']=lat
    edidict['lon']=lon
    edidict['frequency']=np.array(freq)
    edidict['z']=[]
    edidict['zvar']=[]
    edidict['tipper']=[]
    edidict['tippervar']=[]
    for ff in range(len(freq)):
        edidict['z'].append(np.array([[zxxr[ff]+1j*zxxi[ff],zxyr[ff]+1j*zxyi[ff]],[zyxr[ff]+1j*zyxi[ff],zyyr[ff]+1j*zyyi[ff]]]))
        edidict['zvar'].append(np.array([[zxxv[ff],zxyv[ff]],[zyxv[ff],zyyv[ff]]]))
    if tfind==1:
        for ff in range(len(txr)):        
            edidict['tipper'].append(np.array([txr[ff]+1j*txi[ff],tyr[ff]+1j*tyi[ff]]))
            edidict['tippervar'].append(np.array([txv[ff],tyv[ff]]))
    else:
        edidict['tipper'].append(np.array([0.0+1j*0.0,0.0+1j*0.0]))
        edidict['tippervar'].append(np.array([0.0+1j*0.0,0.0+1j*0.0]))
    #make things into arrays for easy manipulation
    edidict['z']=np.array(edidict['z'])
    edidict['zvar']=np.array(edidict['zvar'])
    edidict['tipper']=np.array(edidict['tipper'])
    edidict['tippervar']=np.array(edidict['tippervar'])
    try:
        edidict['rotation']=np.array(rot)
    except UnboundLocalError:
        edidict['rotation']=0
    try:
        edidict['tiprotation']=np.array(trot)
    except UnboundLocalError:
        edidict['tiprotation']=0
    
    if edidict['frequency'][0]<edidict['frequency'][-1]:
        edidict['frequency']=edidict['frequency'][::-1]
        edidict['z']=edidict['z'][::-1,:,:]
        edidict['zvar']=edidict['zvar'][::-1,:,:]
        edidict['tipper']=edidict['tipper'][::-1,:]
        edidict['tippervar']=edidict['tippervar'][::-1,:]
        edidict['rotation']=edidict['rotation']
        edidict['tiprotation']=edidict['tiprotation']
        print 'Flipped arrays so frequency is decending'
    return edidict

def combineEdifiles(edifile1,edifile2,nread1,nread2):
    """
    combineEdifiles(edifile1,edifile2,read1,read2) will combine edifile1 with
    edifile2 according to read1 and read2. It will combine frequencies,
    impedance and tipper.
    
    Note nread1 is from the start for edifile1 and nread2 is from end for 
    edifile2.
    
    Inputs:
        edifile1 = full path to edifile with the largest frequencies
        edifile2 = full path to edifile with smallest frequencies, longest 
                    periods
        nread1 = integer of number of frequencies to read in from edifile1.  
                Index is from the start, ie [0:nread1]
        nread2 = integer of number of frequencies to read in from edifile2.
                Index is from end of frequencies, ie [-nread2:]
    Outputs:
        newedifile = full path to new edifile, which is saved to edifile1 with
                    a c before .edi
    """
    
    edidict1=readedi(edifile1)
    edidict2=readedi(edifile2)
    
    znew=np.append(edidict1['z'][:nread1],edidict2['z'][-nread2:],axis=0)
    zvarnew=np.append(edidict1['zvar'][:nread1],edidict2['zvar'][-nread2:],axis=0)
    freqnew=np.append(edidict1['frequency'][:nread1],
                      edidict2['frequency'][-nread2:],axis=0)
    tipnew=np.append(edidict1['tipper'][:nread1],edidict2['tipper'][-nread2:],axis=0)
    tipvarnew=np.append(edidict1['tippervar'][:nread1],
                        edidict2['tippervar'][-nread2:],axis=0)
    
    newedifile=rewriteedi(edifile1,znew,zvarnew,freqnew=freqnew,newfile='n',
                             tipnew=tipnew,tipvarnew=tipvarnew)
    return newedifile
def sil(iniline):
    """sil(iniline) will split a single line written in an .ini file
    for burpinterface and return the list of strings."""
    
    inistr=iniline.replace('\n','')
    linelst=inistr.split('=')
    
    return linelst[1]

def readStationInfo(stationinfofile):
    """readStationInfo(stationinfofile) will read in a .txt (tab 
    delimited) or .csv(comma delimited)file that has the following information: 
    station name, latitude, longitude, elevation, date collected, 
    components measured (a number: 4 for ex,ey,bx,by, 5 for 
    ex,ey,bx,by,bz, 6 for ex,ey,bx,by,bz,tp), magnetic type (bb or lp), dipole 
    lengths (m), data logger gain (very low=2.5,low=1,high=.1),interface box 
    gain (10 or 100), bz correction for longperiod data and the bx,by,bz 
    components as they were measured in the field. Use headers: station, lat,
    lon, elev, date, mcomps, magtype, ex, ey, dlgain, igain, lpbzcor, magori. 
    Returns a list of the dictionaries with found information."""
    
    #figure out what type of file it is and choose the appropriate delimeter
    filebase=os.path.basename(stationinfofile)
    if filebase.find('.txt')>=0:
        delimiter='\t'
    elif filebase.find('.csv')>=0:
        delimiter=','
        
    #get header information    
    infofid=file(stationinfofile,'r')
    infolines=infofid.readlines()
    hdrlinestr=infolines[0].rstrip()
    hdrlinestr=hdrlinestr.lower()
    
    if filebase.find('.csv')>=0:
        if hdrlinestr.find(',')==-1:
            delimiter=' '
            
    hdrlineslst=hdrlinestr.split(delimiter)
 
    info=[]
    for ii in range(1,len(infolines)):
        infodict={}
        for jj in range(len(hdrlineslst)):
            infoline=infolines[ii]
            infoline=infoline.split(delimiter)
            infostr=infoline[jj].strip()
            infostr=infostr.replace('"','')
            infodict[hdrlineslst[jj]]=infostr
        info.append(infodict)
    
    return info

def getStationInfo(stationinfofile,station,mapversion=23):
    """getStationInfo(stationinfofile,station) returns info for the nominated
    station from the stationinfofile as a dictionary with key words in the 
    hdrinfo."""
    
    info=readStationInfo(stationinfofile)
    for ii in range(len(info)):
        if info[ii]['station'].lower()==station.lower():
            #convert UTM to decimal latitude and longitude using 
            #LatLongUTMconversion.py
            #note: the system is wgs84, if you use a different system change the
            #number as per LatLongUTMconversion.py
            try:
                info[ii]['easting']
                info[ii]['northing']
                lat,lon=utm2ll.UTMtoLL(mapversion,float(info[ii]['northing']),
                           float(info[ii]['easting']),info[ii]['zone'])
                info[ii]['lat']=sigfigs(str(lat),12)
                info[ii]['long']=sigfigs(str(lon),12)
            except KeyError:
                #if there is no decimal place in lat or long probably northing
                #and easting so convert to lats and longs
                if info[ii]['lat'].find('.')==-1 or \
                    info[ii]['long'].find('.')==-1:
                    pass
                elif info[ii]['lat'].find(':')>=0 or info[ii]['long'].find(':')>=0:
                    #if in HH:MM:SS.SS format convert to decimal
                    #convert lat
                    latlst=info[ii]['lat'].split(':')
                    latstr=str(int(latlst[0]))+'.'+str(float(latlst[1])/60
                            +float(latlst[2])/3600)[2:]
                    info[ii]['lat']=latstr
                    #convert long
                    lonlst=info[ii]['long'].split(':')
                    lonstr=str(int(lonlst[0]))+'.'+str(float(lonlst[1])/60
                            +float(lonlst[2])/3600)[2:]
                    info[ii]['long']=lonstr
                        
                
                else:
                    pass
            #if there is no bz correction for long period add one for consistancy
            try:
                info[ii]['lpbzcor']
            except KeyError:
                info[ii]['lpbzcor']='0'
            return info[ii]
        else:
            pass
        
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

def removeStaticShift(edifile,stol=.2,dm=1000):
    """
    removeStaticShift(edifile,stol=.2,dm=1000) will remove static shift by 
    calculating the median of respones of near by stations, within dm.  If the 
    ratio of the station response to the median on either side of 1+-stol then 
    the impedance tensor for that electric component will be corrected for 
    static shift.
    
    Inputs:
        edifile = full path to edi file. Note nearby stations will be looked 
                  for in the dirname of edifile.  So have all edis in one folder
        stol = ratio tolerance.  If there is no static shift the ratio between
               the response and the median response should be 1, but noise and 
               other factors can be present so a tolerance arount 1 is assumed.
        dm = nearby station radius in meters.  If there is no station within 
             that radius then no static shift will be corrected for.
    
    Outputs:
        newedifile = full path to new edifile
    """
    
    #get directory name where all edi files should be
    edipath=os.path.dirname(edifile)
    
    #make a list of filename from edipath
    edifilelst=[os.path.join(edipath,dedifile) 
                for dedifile in os.listdir(edipath) if dedifile.find('.')!=-1]    
                
    #read in the edifile into a dictionary
    edidict=readedi(edifile)
    
    #compute resistivity of edifile
    res,phase=imp2resphase(edidict['z'],edidict['frequency'])
    
    #loop over files to find nearby stations
    reslst=[]
    statlst=[]
    zone,northing,easting=utm2ll.LLtoUTM(23,edidict['lat'],edidict['lon'])   
    for kk,kedi in enumerate(edifilelst):
        kedidict=readedi(kedi)
        zone,dn,de=utm2ll.LLtoUTM(23,kedidict['lat'],kedidict['lon'])        
        deltad=np.sqrt((dn-northing)**2+(de-easting)**2)
        if deltad<=dm:
            kres,kphase=imp2resphase(kedidict['z'],kedidict['frequency'])
            reslst.append(kres)
            statlst.append(kk)
            
    #make resistivity list an array for easy manipulating
    reslst=np.array(reslst)
    
    #calculate static shift of x-components
    staticx=np.mean(res[1,0,:]/np.median(reslst[:,1,0,:],axis=0))
    
    #see if it is within the tolerance level    
    if staticx<1-stol or staticx>1+stol:
        edidict['z'][:,0,:]=edidict['z'][:,0,:]/np.sqrt(staticx)
        print 'X-Shifted '+edidict['station']+' by '+str(1/np.sqrt(staticx))
        xyn='y'
    else:
        xyn='n'
    
    #calculate static shift of y-components
    staticy=np.mean(res[2,0,:]/np.median(reslst[:,2,0,:],axis=0))
    
    #see if it is within the tolerance level
    if staticy<1-stol or staticy>1+stol:
        edidict['z'][:,1,:]=edidict['z'][:,1,:]/np.sqrt(staticy)
        print 'Y-Shifted '+edidict['station']+' by '+str(1/np.sqrt(staticy)) 
        yyn='y'
    else:yyn='n'
    
    #if there was a correction write a new edi file
    if xyn=='y' or yyn=='y':
        nedi=rewriteedi(edifile,edidict['z'],edidict['zvar'])
        savepath=os.path.join(edipath,'SS')
        if not os.path.exists(savepath):
            os.mkdir(savepath)
        newedifile=os.path.join(savepath,os.path.basename(nedi)[0:4]+'s.edi')
        shutil.move(nedi,newedifile)
        
        print 'Wrote: ',newedifile
        return newedifile
        
    #if no correction was made return the same edifile
    else:
        print 'No Static Shift Correction for ',edidict['station']
        return edifile
        
def makeDayFoldersFromObservatoryData(obsdirpath,savepath,startday='000',d=10,
                                      ndays=10,df=1,complst=['BX','BY'],
                                      station='CNB',rowskip=6):
    """
    makeDayFoldersFromObservatoryData will take observatory data and put it
    into day folders for processing.
    
    Input:
        obsdirpath = path to observatory data.  If the observatory data comes
                     in one long file input the file, program know's its a file
                     by searching for the dot.
        savepath = path to save the day files, should be in the same directory
                   as collected data. ie datapath/observatory_station_name
        startday = UTC start day as a string of length 3
        d = decimation factor of collected data sampling frequency to get to 
            observatory sampling frequency.  So if you sampled at 1000 Hz and 
            the observatory data is at 1 sec, d=1000
        ndays = number of days, only used if observatory data is in one big file
        df = sampling frequency of observatory data in Hz
        complst = components to make into files, can add 'BZ' if you like
        station = station name of observatory data
        rowskip = number of rows to skip from observatory file, usually the 
                  header ends at row 6.
    
    Outputs:
        files will be saved to savepath/day/Comb00to24ddf/station00to24ddf.cnb
    """

    #make save directory if it doesn't exist
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    
    #start day
    sday=startday
    
    #number of points per day
    dl=86400*df
    
    #folder to save each day to
    folder='Comb00to24d'+str(d)
    
    if obsdirpath.find('.')>0:
        for comp in complst:
            bx=np.loadtxt(os.path.join(obsdirpath,station+comp.lower()+'.txt'))
            for dd in range(ndays):
                lbx=bx[dd*dl:(dd+1)*dl]
                day=str(sday+dd)
                if not os.path.exists(os.path.join(savepath,day)):
                    os.mkdir(os.path.join(savepath,day))
                if not os.path.exists(os.path.join(savepath,day,folder)):
                    os.mkdir(os.path.join(savepath,day,folder))
                np.savetxt(os.path.join(savepath,day,folder,
                                        station+'00to24d10.'+comp),
                                        lbx,fmt='%.2f')
    else:
        for ii,filename in enumerate(os.listdir(obsdirpath)):
            #get day string            
            dayi=int(sday)+ii
            if dayi<10:
                day='00'+str(dayi)
            if dayi<100 and dayi>=10:
                day='0'+str(dayi)
            else:
                day=str(dayi)            
            
            #make folder
            if not os.path.exists(os.path.join(savepath,day)):
                os.mkdir(os.path.join(savepath,day))

            #make comb folder
            if not os.path.exists(os.path.join(savepath,day,folder)):
                os.mkdir(os.path.join(savepath,day,folder))
            #load file 
            bx,by,bz=np.loadtxt(os.path.join(obsdirpath,filename),
                                delimiter=',',skiprows=rowskip,usecols=(1,2,3),
                                unpack=True)
            if len(bx)!=dl:
                raise ValueError('File length is not a full day, check time column in'+\
                                'original file '+os.pathjoin(obsdirpath,filename))
            #savefiles
            savefnpath=os.path.join(savepath,day,folder,
                                    station+'00to24d'+str(d)+'.')
            for ii,bn in enumerate([bx,by,bz]):
                np.savetxt(savefnpath+complst[ii],bn,fmt='%.7g')
                print 'Saved file: ',savefnpath+complst[ii]
                                