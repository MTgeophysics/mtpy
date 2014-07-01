# -*- coding: utf-8 -*-
"""
Module of tool that can be applied to BIRRP outputs

Created on %(date)s

@author: %(Jared Peacock)s
"""

import numpy as np
import os
import fnmatch
import datetime
from scipy import interpolate
import time
import mtpy.legacy.mttools as mt
import shutil
import subprocess
#import winsound as ws

#short spaces 3 spaces
tsp='   '
#long spaces 6 spaces
lsp='      '

def finishBeep():
    pass
    #ws.Beep(1000,80)
    #ws.Beep(300,150)
#    ws.Beep(150,40)
#    ws.Beep(200,40)
#    ws.Beep(300,80)
#    ws.Beep(1000,100)
    #for ii in np.arange(50,1200,100):
        #ws.Beep(ii,25)
#    ws.Beep(523,50)
#    ws.Beep(100,25)
#    ws.Beep(264,50)

def readProDict(processinginfofile,dirpath,filtered=None,rrfiltered=None,
                kwargs=None):
    """
    readProDict will read in the information from processing infofile and output
    as a list of dictionaries.  
    Inputs:
        processinginfofile=tab delimited text file with appropriate information.
        dirpath = directory path to where station folders are
        
    Outputs:
        plst = list of dictionaries with key words related to headers of txt 
                file
    
    """
    
    pfid=open(processinginfofile,'r')
    plines=pfid.readlines()
    pkeys=plines[0].rstrip()
    pkeys=pkeys.split('\t')
    plst=[]
    for pline in plines[1:]:
        pstr=pline.rstrip()
        pstr=pstr.split('\t')
        if len(pstr)>1:
            pdict={}
            for kk,pkey in enumerate(pkeys):
                if pstr[kk].find('"')>=0:
                    pstr[kk]=pstr[kk].replace('"','')
                pdict[pkey]=pstr[kk]
        pdict['dirpath']=dirpath
        if filtered!=None:
            pdict['filtered']='Filtered'
        if rrfiltered!=None:
            pdict['rrfiltered']='Filtered'
        plst.append(pdict)
        
    pfid.close()
    
    return plst

def scriptfilePrep(prepdict):
    """
    scriptfilePrep(prepdict) will calculate parameters, combine files and 
    convert if need be. Keys for prepdict are:
        dirpath = dirpath to station folders -> d:\Data
        cdirpath = full path to time series files if not entered will default to
                    dirpath\station\day
        cdirpathr = full path to remote reference time series if not entered
                    dirpath\rrstation\day
        station = station name
        day = UTC day string of length 3 -> 035
        rrstation = remote reference station name
        cacherate = length of data files string of hhmmss -> 1 hour data files
                    is 010000, 10 min is 001000
        start = start time for combining files in hhmmss
        stop = stop time for combining files in hhmmss
        rrstart = start time of remote reference in hhmmss
        rrstop = stop time of remote reference in hhmmss
        dec = decimation factor (0 for no decimation, needs to be an integer)
        df = sampling frequency in Hz
        magori = orientation of magnetic components as N,E,Z -> for measured BX 
                 to East input as [BY,BX,BZ] omit BZ if not measured
        elecori = orientation of electric components as N,E -> for EX measured
                  to N input as [EX,EY]
        rrmagori = orientation of remote reference magnetometers as N,E,Z omit
                   Z if not measured
        filtered = 'Filtered' if filtered, don't enter if not filtered
        rrfiltered = 'Filtered' if remote reference is filtered 
        fdict = dictionary for adaptive notch filter see 
                MTtools.adaptiveNotchFilter() for details
        
    returns procesing dictionary with keys:
        cfilelst = list of files combined full path
        rrcfilelst = list of remote reference files combined full path
        station = station
        deltat = sampling rate (+) for time (-) for frequency
        nfft = max window length
        nsctmax = max number of windows
        nskipr = number of points to skip for remote reference
    """
    #print the dictionary going into combining files
    
#    pkeys=prepdict.keys()
#    for key in pkeys.sort():
#        print key+': '+str(prepdict[key])
        
    station=prepdict['station']
    if len(prepdict['rrstation'].split(';'))>1:    
        rrstation=prepdict['rrstation'].split(';')
    else:
        rrstation=[prepdict['rrstation'],prepdict['rrstation'],
                  prepdict['rrstation']]

    try:
        fdictc=prepdict['fdict']
        if type(fdictc) is str:
            if len(fdictc)>5:
                
                fd1=fdictc[1:-1].split(',')
                fd={}
                for ff in fd1:
                    ff=ff.split(':')
                    if ff[1].find('[')>=0:
                        fv=ff[1][1:-1].split(';')
                        fd[ff[0][1:-1]]=[float(vv) for vv in fv]
                    else:
                        fd[ff[0][1:-1]]=float(ff[1])
                fdictc=fd
            elif fdictc=='-1':
                fdictc=None
    except KeyError:
        fdictc=None
    #cacherate
    cacherate=prepdict['cacherate']
    
    #list of components to combine
    combcomplst=prepdict['elecori'].split(',')+prepdict['magori'].split(',')
    complstp=prepdict['elecori'].split(',')+prepdict['magori'].split(',')
    combcomplstr=prepdict['rrmagori'].split(',')
    complstpr=prepdict['rrmagori'].split(',')
    dec=int(prepdict['dec'])
    if dec==0:
        dec=1
    
    if len(prepdict['day'].split(';'))>1 or len(prepdict['day'].split(','))>1:
        
        dayslst=prepdict['day'].split(';')
        dstartlst=prepdict['start'].split(';')
        dstoplst=prepdict['stop'].split(';')
        drrstartlst=prepdict['rrstart'].split(';')
        drrstoplst=prepdict['rrstop'].split(';')
        cfilelst=[]
        rrcfilelst=[]
        nread=[]
        nskipr=[]
        
        for ii,days in enumerate(dayslst):
            combcomplst=prepdict['elecori'].split(',')+prepdict['magori'].split(',')
            combcomplstr=prepdict['rrmagori'].split(',')
            if len(days.split(','))>1:
                daylst=days.split(',')
                startlst=dstartlst[ii].split(',')
                stoplst=dstoplst[ii].split(',')
                dlst=[]
                nreadi=0
                for dd in range(len(daylst)):
                    ddict={}
                    ddict['day']=daylst[dd]
                    ddict['start']=startlst[dd]
                    ddict['stop']=stoplst[dd]
                    try:
                        ddict['filt']=prepdict['filtered']
                    except KeyError:
                        pass
                    dlst.append(ddict)
                    nreadi+=float(ddict['stop'][0:2])-float(ddict['start'][0:2])
                nreadi=int(nreadi*3600*float(prepdict['df'])/float(dec))
                cfilelsti,fileslsti=mt.combineFiles(prepdict['dirpath'],
                                                    prepdict['station'],
                                                    dlst,cacherate,
                                                    complst=combcomplst,
                                                    dec=dec,fdict=fdictc)
                #remote reference
                if rrstation[ii]!=station:
                    rrstartlst=drrstartlst[ii].split(',')
                    rrstoplst=drrstoplst[ii].split(',')
                    rrdlst=[]
                    #nreadr=0
                    for dd in range(len(daylst)):
                        rrddict={}
                        rrddict['day']=daylst[dd]
                        rrddict['start']=rrstartlst[dd]
                        rrddict['stop']=rrstoplst[dd]
                        try:
                           rrddict['filt']=prepdict['rrfiltered']
                        except KeyError:
                            pass
                        rrdlst.append(rrddict)
                        
                    rrcfilelsti,rrfileslsti=mt.combineFiles(prepdict['dirpath'],
                                                            rrstation[ii],
                                                            rrdlst,
                                                            cacherate,
                                                            complst=combcomplstr,
                                                            dec=dec,
                                                            fdict=fdictc)
                    #get number of points to skip for remote reference
                    nskipri=int((float(dlst[0]['start'][0:2])-\
                            float(rrdlst[0]['start'][0:2]))*\
                           (float(prepdict['df'])/float(dec)))
                else:
                    rrcfilelsti=[]
                    for cfile in cfilelsti:
                        for rcomp in combcomplstr:
                            if cfile.find(rcomp)>=0:
                                rrcfilelsti.append(cfile)
                    nskipri=0 
            #if multiple files but not continuous.  
            else:
                day=days
                start=dstartlst[ii]
                stop=dstoplst[ii]
                try:
                    filt=prepdict['filtered']
                except KeyError:
                    filt=''
                        
                #get number of points to read
                nreadi=int((float(stop[0:2])-float(start[0:2]))*3600*\
                        float(prepdict['df'])/dec)
                
                #make a directory path
                if filt!='':
                    cdirpath=os.path.join(prepdict['dirpath'],station,day,filt)
                else:
                    cdirpath=os.path.join(prepdict['dirpath'],station,day)
                #combine files
                #print 'combining ',cdirpath,start,stop,cacherate,combcomplst,dec
                cfilelsti,fileslsti=mt.combineFewFiles(cdirpath,
                                                       station,
                                                       start,stop,cacherate,
                                                       complst=combcomplst,
                                                       d=dec,fdict=fdictc)

                #remote reference
                if rrstation[ii]!=station:
                    rrday=days
                    rrstart=drrstartlst[ii]
                    rrstop=drrstoplst[ii]
                    try:
                        filt=prepdict['rrfiltered']
                    except KeyError:
                        filt=''
                    if filt!='':
                        rrcdirpath=os.path.join(prepdict['dirpath'],
                                                rrstation[ii],rrday,filt)
                    else:
                        rrcdirpath=os.path.join(prepdict['dirpath'],
                                                rrstation[ii],rrday)  
        
                    rrcfilelsti,rrfileslsti=mt.combineFewFiles(rrcdirpath,
                                                               rrstation[ii],
                                                               rrstart,rrstop,
                                                               cacherate,
                                                               complst=
                                                               combcomplstr,
                                                               d=dec,
                                                               fdict=fdictc)
                    #get number of points to skip for remote reference
                    nskipri=(int((float(dstartlst[ii][0:2])-\
                                float(drrstartlst[ii][0:2]))*\
                                (float(prepdict['df'])/float(dec))))
                else:
                    rrcfilelsti=[]
                    for cfile in cfilelsti:
                        for rcomp in combcomplstr:
                            if cfile.find(rcomp)>=0:
                                rrcfilelsti.append(cfile)
                    nskipri=0  
            #append files to a list
            cfilelst.append(cfilelsti)
            rrcfilelst.append(rrcfilelsti)
            nread.append(nreadi)
            nskipr.append(nskipri)                                        
    #else normal input for one day    
    else:  
        #dirpath for timeseries
        try:
            cdirpath=prepdict['cdirpath']
        except KeyError:
            try:
                cdirpath=os.path.join(prepdict['dirpath'],prepdict['station'],
                                      prepdict['day'],prepdict['filtered'])
            except KeyError:
                cdirpath=os.path.join(prepdict['dirpath'],prepdict['station'],
                                      prepdict['day'])
        #dirpath for remote reference station 
        try:
            cdirpathr=prepdict['cdirpathr']
        except KeyError:
            try:
                cdirpathr=os.path.join(prepdict['dirpath'],rrstation[0],
                                       prepdict['day'],prepdict['rrfiltered'])
            except KeyError:
                cdirpathr=os.path.join(prepdict['dirpath'],rrstation[0],
                                       prepdict['day'])
        
        #start and stop time for both data and remote reference time series
        stime=prepdict['start']
        etime=prepdict['stop']
        
        stimer=prepdict['rrstart']
        etimer=prepdict['rrstop']
    
        #Combine time series files
        cfilelst,fileslst=mt.combineFewFiles(cdirpath,station,stime,etime,
                                             cacherate,complst=combcomplst,
                                             d=dec,fdict=fdictc)
        
        #combine remote reference files
        if rrstation[0]!=station:
            rrcfilelst,rrfileslst=mt.combineFewFiles(cdirpathr,rrstation[0],
                                                     stimer,etimer,cacherate,
                                                     complst=combcomplstr,
                                                     d=dec,fdict=fdictc)
            #get number of points to skip for remote reference
            if float(dec)==float(0):
                dec=1
            else:
                pass 
            
            nskipr=int((float(stime[0:2])-float(stimer[0:2]))*\
                (float(prepdict['df'])/float(dec)))                                                          
        else:
            rrcfilelst=[]
            for cfile in cfilelst:
                for rcomp in combcomplstr:
                    if cfile.find(rcomp)>=0:
                        rrcfilelst.append(cfile)
            nskipr=0
             
                        
        #get number of points to read, length of fft and maximum decimations
        if float(dec)==float(0):
            dec=1
        else:
            pass        
        nread=int(3600*(float(prepdict['df'])/float(dec))*\
             (float(etime[0:2])-float(stime[0:2])))
        

    #sampling frequency
    if float(prepdict['dec'])==0.0:
        prepdict['dec']=str(1)
    deltat='-'+str(float(prepdict['df'])/float(prepdict['dec'])) 
    
    try:
        nfft=prepdict['nfft']
    except KeyError:
        #get max length of time window
        npow=np.floor(np.log2(float(max(nread)))-16)
        if npow>6:
            nfftpow=17
        elif npow>=2 and npow<=6:
            nfftpow=16
        elif npow>=-2 and npow<2:
            nfftpow=15
        elif npow>=-6 and npow<-2:
            nfftpow=14
        elif npow>=-6 and npow<-2:
            nfftpow=13
        nfft=2**nfftpow
    try:
        nsctmax=prepdict['nsctmax']
    except KeyError:
        #get max length of time window
        npow=np.floor(np.log2(float(max(nread)))-16)
        #get max number of windows
        nsctmax=nfftpow-4
            
    
    #make sure the cfilelst is in order according to complst
    cfarray=np.array(cfilelst)
    try:
        nds,ndf=cfarray.shape
    except ValueError:
        nds=0
    
    if nds==0:
        cflst=[]
        for comp in complstp:
            for cfile in cfarray:
                if fnmatch.fnmatch(cfile,'*.'+comp):
                    cflst.append(cfile)
    else:
        cflst=[]
        for ii in range(nds):
            cfnlst=[]
            for comp in complstp:
                for cfile in cfarray[ii]:    
                    if fnmatch.fnmatch(cfile,'*.'+comp):
                        cfnlst.append(cfile)
            cflst.append(cfnlst)
    
    #make sure the rrcfilelst is in order according to complst
    rrcfarray=np.array(rrcfilelst)
    try:
        rrnds,rrndf=rrcfarray.shape
    except ValueError:
        rrnds=0
    
    if rrnds==0:
        rrcflst=[]
        for compr in complstpr:
            for rcfile in rrcfarray:
                if fnmatch.fnmatch(rcfile,'*.'+compr):
                    rrcflst.append(rcfile)
    else:
        rrcflst=[]
        for ii in range(rrnds):
            rrcfnlst=[]
            for compr in complstpr:
                for rcfile in rrcfarray[ii]:    
                    if fnmatch.fnmatch(rcfile,'*.'+compr):
                        rrcfnlst.append(rcfile)
            rrcflst.append(rrcfnlst)    
    
    processingdict={}
    processingdict['cfilelst']=cflst
    processingdict['rrcfilelst']=rrcflst
    processingdict['station']=station
    processingdict['deltat']=deltat
    processingdict['nfft']=nfft
    processingdict['nsctmax']=nsctmax
    processingdict['nread']=nread
    
    for pkey in prepdict.keys():
        processingdict[pkey]=prepdict[pkey]
    try:
        processingdict['nskipr']=prepdict['nskipr']
    except KeyError:
        try:
            if len(nskipr)!=len(processingdict['cfilelst']):
                nskipr=[0 for ii in range(len(processingdict['cfilelst']))]
        except TypeError:
            nskipr=nskipr
        processingdict['nskipr']=nskipr

    return processingdict
    
    
def writeScriptfile(processingdict):
    """
    writeScriptfile(processingdict) will write a script file for BIRRP using 
    info in processingdict which is a dictionary with keys:
        cfilelst = list of combined filenames (full path) ordered as:
            [electric N, electric E, mag N, mag E, mag Z] can omit Z if not 
            measured 
        rrcfilelst = list of combined remote reference filenames (full 
            path) ordered as: [mag N, mag E, mag Z] omit z if not measured 
        station = station name
        magtype = bb for broadband lp for longperiod
        any BIRRP parameters you desire, (if not input default values will be 
        assigned) namely:
            ilev = processing mode 0 for basic and 1 for advanced RR-2 stage
            nout = Number of Output time series (2 or 3-> for BZ)
            ninp = Number of input time series for E-field (1,2,3) 
            nref = Number of reference channels (2 for MT)
            nrr = bounded remote reference (0) or 2 stage bounded influence (1)
            tbw = Time bandwidth for Sepian sequence
            deltat = Sampling rate (+) for (s), (-) for (Hz)
            nfft = Length of FFT (should be even)
            nsctinc = section increment divisor (2 to divide by half)
            nsctmax = Number of windows used in FFT
            nf1 = 1st frequency to extract from FFT window (>=3)
            nfinc = frequency extraction increment 
            nfsect = number of frequencies to extract
            mfft = AR filter factor, window divisor (2 for half)
            uin = Quantile factor determination
            ainlin =  Residual rejection factor low end (usually 0)
            ainuin = Residual rejection factor high end (.95-.99)
            c2threshb = Coherence threshold for magnetics (0 if undesired)
            c2threshe = Coherence threshold for electrics (0 if undesired)
            nz= Threshold for Bz (0=separate from E, 1=E threshold, 2=E and B) 
                Input if 3 B components else None
            c2thresh1= Squared coherence for Bz, input if NZ=0, Nout=3
            perlo = longest period to apply coherence threshold over
            perhi = shortes period to apply coherence threshold over
            ofil = Output file root(usually three letters, can add full path)
            nlev = Output files (0=Z; 1=Z,qq; 2=Z,qq,w; 3=Z,qq,w,d)
            nprej = number of frequencies to reject
            prej = frequencies to reject (+) for period, (-) for frequency
            npcs = Number of independent data to be processed (1 for one 
                   segement)
            nar = Prewhitening Filter (3< >15) or 0 if not desired',
            imode = Output file mode (0=ascii; 1=binary; 2=headerless ascii; 
                    3=ascii in TS mode',
            jmode = nput file mode (0=user defined; 1=start time 
                    YYYY-MM-DD HH:MM:SS)',
            nread = Number of points to be read for each data set  
                    (if segments>1 -> npts1,npts2...)',
            nfil = Filter parameters (0=none; >0=input parameters; <0=filename
            nskip = Skip number of points in time series (0) if no skip, 
                    (if segements >1 -> input1,input2...)',
            nskipr = Number of points to skip over (0) if none,
                     (if segements >1 -> input1,input2...)',
            thetae = Rotation angles for electrics (relative to geomagnetic 
                     North)(N,E,rot)',
            thetab = Rotation angles for magnetics (relative to geomagnetic 
                     North)(N,E,rot)',
            thetar = Rotation angles for calculation (relative to geomagnetic 
                     North)(N,E,rot)'
                     
        => see BIRRP Manual for more details on the parameters
        => see A. D. Chave and D. J. Thomson [1989,2003,2004] for more
            information on Bounded influence and robust processing.
    
    ouputs:
        scriptfile (full path typically(\dirpath\station\day\CombSTtoETddec\BF)
        birrpdict = dictionary of all birrp parameters in script file
    
    """
     
    #===================================================================
    # Write a script file for BIRRP, Chave et al. [2004]            
    #===================================================================
    #print processingdict
    #compute how many timeseries and days there are 
    #ndf = # of files per day
    #nds = # of day
    cfarray=np.array(processingdict['cfilelst'])    
    
    try:
        nds,ndf=cfarray.shape
    except ValueError:
        ndf=cfarray.shape[0]
        nds=0
    
    if nds==0:
        bfpath=os.path.join(os.path.dirname(processingdict['cfilelst'][0]),
                            'BF')
        npcs=1
    elif nds==1:
        nds=0
        npcs=1
        processingdict['cfilelst']=processingdict['cfilelst'][0]
        processingdict['rrcfilelst']=processingdict['rrcfilelst'][0]  
        
        bfpath=os.path.join(os.path.dirname(processingdict['cfilelst'][0]),
                            'BF')

    else:
        bfpath=os.path.join(os.path.dirname(processingdict['cfilelst'][0][0]),
                            'BF')
        npcs=int(nds)
    #make a directory to put BIRRP Files (BF)
    if not os.path.exists(bfpath):
        os.mkdir(bfpath)
        print 'Made directory: ', bfpath
#    else:
#        nruns=len([folder for folder in os.listdir(bfpath) if bfpath.find('.'==-1)])
#        bfpath=os.path.join(bfpath,'Run{0}'.format(nruns+1))
    #output file stem, full path
    ofil=os.path.join(bfpath,processingdict['station']) 
    
    #are the frequencies ok always y in 0 interactive mode
    yesno='y'
    
    
    #see if user input some values if not put default ones in
    try:
        magtype=processingdict['magtype']
    except KeyError:
        magtype=input('input magnetometer type broadband (bb), long period(lp)')

    #mode to process: default is basic
    try:
        ilev=int(processingdict['ilev'])
    except KeyError:
        ilev=0
    
    #number of output channels
    try:
        nout=int(processingdict['nout'])
    except KeyError:
        if nds!=0:
            if ndf==5:
                nout=3
            else:
                nout=2
        elif ndf==5:
            nout=3
        else:
            nout=2
    
    #number of input channels default is 2
    try:
        ninp=int(processingdict['ninp'])
    except KeyError:
        ninp=2
        
    #time bandwidth window size            
    try:
        tbw=float(processingdict['tbw'])
    except KeyError:
        tbw=2
        
    #------------Options for Advanced mode-------------------------
    if ilev==1:
        #Advanced: number of remote reference channels
        try:
            nref=int(processingdict['nref'])
        except KeyError:
            nref=2
            
        #Advanced: remote reference type processing
        try:
            nrr=int(processingdict['nrr'])
        except KeyError:
            nrr=1            
            
        #Advanced: magnetic coherence threshold
        try:
            c2threshb=float(processingdict['c2threshb'])
        except KeyError:
            c2threshb=0
            
        #Advanced: window increment divisor
        try:
            nsctinc=int(processingdict['nsctinc'])
        except KeyError:
            nsctinc=2
            
        #first frequency to extract
        try:
            nf1=int(processingdict['nf1'])
        except KeyError:
            nf1=int(tbw)+2
        
        #frequency increment
        try:
            nfinc=float(processingdict['nfinc'])
        except KeyError:
            nfinc=int(tbw)
            
        #number of frequencies to extract
        try:
            nfsect=int(processingdict['nfsect'])
        except KeyError:
            nfsect=2
            
        #number AR filter is divided by
        try:
            mfft=float(processingdict['mfft'])
        except KeyError:
            mfft=2
        
        #Advanced: lower bound of leverage point rejection
        try:
            ainlin=float(processingdict['ainlin'])
        except KeyError:
            ainlin=.0001
        
        #Advanced: coherence threshold low period
        try:
            perlo=float(processingdict['perlo'])
        except KeyError:
            perlo=1000
            
        #Advanced: coherenct threshold high period
        try:
            perhi=float(processingdict['perhi'])
        except KeyError:
            perhi=.001
            
        #Advanced:  number of frequencies to reject
        try:
            nprej=int(processingdict['nprej'])
        except KeyError:
            nprej=0
        
        #Advanced
        try:
            prej=processingdict['prej'].split(',')
            if type(prej) is list:
                prej=[float(ff) for ff in prej]
            if nprej!=len(prej):
                nprej=len(prej)
        except KeyError:
            prej=[]
            
    #---------------------Options for Basic Mode--------------------------
    
    #time series sampling rate
    try:
        deltat=float(processingdict['deltat'])
    except KeyError:
        deltat=-100
    
    #max length of fft window    
    try:
        nfft=int(processingdict['nfft'])
    except KeyError:
        nfft=2**16

    #maximum number of sections
    try:
        nsctmax=int(processingdict['nsctmax'])
    except KeyError:
        nsctmax=12
    
    #quantile factor            
    try:
        uin=float(processingdict['uin'])
    except KeyError:
        uin=0
    
    #upper bound of leverage point rejection
    try:
        ainuin=float(processingdict['ainuin'])
    except KeyError:
        ainuin=.9999
 
    #electric channel coherence threshold
    try:
        c2threshe=float(processingdict['c2threshe'])
    except KeyError:
        c2threshe=0

    #Bz coherency threshold mode
    try:
        nz=int(processingdict['nz'])
    except KeyError:
        if nout==3:
            nz=0
        else:
            nz=None
    
    #Bz coherence threshold
    try:
        c2threshe1=float(processingdict['c2threshe1'])
    except KeyError:
        if nout==3:
            c2threshe1=0
        else:
            c2threshe1=None
    
    try:
        nlev=int(processingdict['nlev'])
    except KeyError:
        nlev=0
    
    try: 
        nread=processingdict['nread']
    except KeyError:
        nread=1440000
    
    try:
        nar=int(processingdict['nar'])
    except KeyError:
        nar=5
        
    try:
        imode=processingdict['imode']
    except KeyError:
        imode=0
        
    try:
        jmode=int(processingdict['jmode'])
    except KeyError:
        jmode=0    
    
    try:
        nfil=int(processingdict['nfil'])
    except KeyError:
        nfil=0
    
    try:
        nskip=processingdict['nskip']
    except KeyError:
        if nds!=0:
            nskip=[0 for ii in range(nds)]
        else:
            nskip=0
        
    try:
        nskipr=processingdict['nskipr']
        #print '1) ', nskipr
        if nds==0:
            nskipr=nskipr
        else:
            if type(nskipr) is list:
                pass
    
            else:        
                nskipr=[nskipr for ii in range(nds)]
    except KeyError:
        if nds==0:
            nskipr=0
        else:
            nskipr=[0 for ii in range(nds)]
            #print '5) ',nskipr
    
    
    try:
        thetae=processingdict['thetae']
    except KeyError:
        #rotation angles
        if magtype=='bb':
            thetae='0,90,180'
        else:
            thetae='0,90,0'
    
    try:
        thetab=processingdict['thetab']
    except KeyError:
        thetab='0,90,0'
    
    try:
        thetaf=processingdict['thetaf']
    except KeyError:
        thetaf='0,90,0'

    #===================================================================
    # Write values to a .script file
    #===================================================================
    #print '+++ ',nskipr
    #print ndf,nds
    #write to a file
    scriptfile=ofil+'.script'
    fid=file(scriptfile,'w')
    if ilev==0: 
        fid.write('{0:d} \n'.format(ilev))
        fid.write('{0:d} \n'.format(nout))
        fid.write('{0:d} \n'.format(ninp))
        fid.write('{0:.3f} \n'.format(tbw))
        fid.write('{0:.3f} \n'.format(deltat))
        fid.write('{0:d},{1:d} \n'.format(nfft,nsctmax))
        fid.write('y \n')
        fid.write('{0:.5f},{1:.5f} \n'.format(uin,ainuin))
        fid.write('{0:.3f} \n'.format(c2threshe))
        #parameters for bz component if ninp=3
        if nout==3:
            if c2threshe==0:
                fid.write('{0:d} \n'.format(0))
                fid.write('{0:.3f} \n'.format(c2threshe1))
            else:
                fid.write('{0:d} \n'.format(nz))
                fid.write('{0:.3f} \n'.format(c2threshe1))
        else:
            pass
        fid.write(ofil+'\n')
        fid.write('{0:d} \n'.format(nlev))
        
    elif ilev==1:
        print 'Writing Advanced mode'
        fid.write('{0:d} \n'.format(ilev))
        fid.write('{0:d} \n'.format(nout))
        fid.write('{0:d} \n'.format(ninp))
        fid.write('{0:d} \n'.format(nref))
        if nref>3:
            nrrlst=np.array([len(rrlst) 
                            for rrlst in processingdict['rrcfilelst']])
            nr3=len(np.where(nrrlst==3)[0])
            nr2=len(np.where(nrrlst==2)[0])
            fid.write('{0:d},{1:d} \n'.format(nr3,nr2))
        fid.write('{0:d} \n'.format(nrr))
        #if remote referencing
        if int(nrr)==0:
            fid.write('{0:.3f} \n'.format(tbw))
            fid.write('{0:.3f} \n'.format(deltat))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nfft,nsctinc,nsctmax))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nf1,nfinc,nfsect))
            fid.write('y \n')
            fid.write('{0:.2g} \n'.format(mfft))        
            fid.write('{0:.5g},{1:.5g},{2:.5g} \n'.format(uin,ainlin,ainuin))
            fid.write('{0:.3f} \n'.format(c2threshe))
            #parameters for bz component if ninp=3
            if nout==3:
                if c2threshe!=0:
                    fid.write('{0:d} \n'.format(nz))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                else:
                    fid.write('{0:d} \n'.format(0))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                if c2threshe1!=0.0 or c2threshe!=0.0:
                    fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            else:
                if c2threshe!=0.0:
                    fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            fid.write(ofil+'\n')
            fid.write('{0:d} \n'.format(nlev))
            fid.write('{0:d} \n'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej=[prej]
                fid.writelines(['{0:.5g} \n'.format(nn) for nn in prej])
        #if 2 stage processing
        elif int(nrr)==1:
            fid.write('{0:.5g} \n'.format(tbw))
            fid.write('{0:.5g} \n'.format(deltat))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nfft,nsctinc,nsctmax))        
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nf1,nfinc,nfsect))
            fid.write('y \n')
            fid.write('{0:.2g} \n'.format(mfft))        
            fid.write('{0:.5g},{1:.5g},{2:.5g} \n'.format(uin,ainlin,ainuin))
            fid.write('{0:.3f} \n'.format(c2threshb))        
            fid.write('{0:.3f} \n'.format(c2threshe))
            if nout==3:
                if c2threshb!=0 or c2threshe!=0:
                    fid.write('{0:d} \n'.format(nz))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                elif c2threshb==0 and c2threshe==0:
                    fid.write('{0:d} \n'.format(0))
                    fid.write('{0:.3f} \n'.format(0))
            if c2threshb!=0.0 or c2threshe!=0.0:
                fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            fid.write(ofil+'\n')
            fid.write('{0:d} \n'.format(nlev))
            fid.write('{0:d} \n'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej=[prej]
                fid.writelines(['{0:.5g} \n'.format(nn) for nn in prej])
        
    fid.write('{0:d} \n'.format(npcs))    
    fid.write('{0:d} \n'.format(nar))    
    fid.write('{0:d} \n'.format(imode))    
    fid.write('{0:d} \n'.format(jmode))    
    #write in filenames 
    if npcs!=1:
        fid.write(str(nread[0])+'\n')
        #write filenames
        for tfile in processingdict['cfilelst'][0]:
            fid.write(str(nfil)+'\n')
            fid.write(tfile+'\n')
            fid.write(str(nskip[0])+'\n')
        for rfile in processingdict['rrcfilelst'][0]:
            fid.write(str(nfil)+'\n')
            fid.write(rfile+'\n')
            fid.write(str(nskipr[0])+'\n')
        for nn in range(1,npcs):
            fid.write(str(nread[nn])+'\n')            
            #write filenames
            for tfile in processingdict['cfilelst'][nn]:
                fid.write(tfile+'\n')
                fid.write(str(nskip[0])+'\n')
            for rfile in processingdict['rrcfilelst'][nn]:
                fid.write(rfile+'\n')
                fid.write(str(nskipr[nn])+'\n')
        
    else:
        if type(nread) is list:
            fid.write(str(nread[0])+'\n')
        else:
            fid.write(str(nread)+'\n')
        #write filenames
        if nds==0:
            for tfile in processingdict['cfilelst']:
                fid.write(str(nfil)+'\n')
                fid.write(tfile+'\n')
                if type(nskip) is list:
                    fid.write(str(nskip[0])+'\n')
                else:
                    fid.write(str(nskip)+'\n')
            for rfile in processingdict['rrcfilelst']:
                fid.write(str(nfil)+'\n')
                fid.write(rfile+'\n')
                if type(nskipr) is list:
                    fid.write(str(nskipr[0])+'\n')
                else:
                    fid.write(str(nskipr)+'\n')
        else:
            for tfile in processingdict['cfilelst'][0]:
                fid.write(str(nfil)+'\n')
                fid.write(tfile+'\n')
                if type(nskip) is list:
                    fid.write(str(nskip[0])+'\n')
                else:
                    fid.write(str(nskip)+'\n')
            for rfile in processingdict['rrcfilelst'][0]:
                fid.write(str(nfil)+'\n')
                fid.write(rfile+'\n')
                if type(nskipr) is list:
                    fid.write(str(nskipr[0])+'\n')
                else:
                    fid.write(str(nskipr)+'\n')
    #write rotation angles
    fid.write(thetae.replace(',',' ')+'\n')
    fid.write(thetab.replace(',',' ')+'\n')
    fid.write(thetaf.replace(',',' '))    
    fid.close()
    
    birrpdict={}
    
    if ilev==0:
        birrpdict['ilev']=ilev
        birrpdict['nout']=nout
        birrpdict['ninp']=ninp
        birrpdict['tbw']=tbw
        birrpdict['nfft']=nfft
        birrpdict['nsctmax']=nsctmax
        birrpdict['uin']=uin
        birrpdict['ainuin']=ainuin
        birrpdict['c2threshe']=c2threshe
        birrpdict['nz']=nz
        birrpdict['c2threshe1']=c2threshe1
        birrpdict['ofil']=ofil
        birrpdict['nlev']=nlev
        birrpdict['nar']=nar
        birrpdict['imode']=imode
        birrpdict['jmode']=jmode
        birrpdict['nfil']=nfil
        birrpdict['nskip']=nskip
        birrpdict['nskipr']=nskipr
        birrpdict['thetae']=thetae
        birrpdict['thetab']=thetab
        birrpdict['thetaf']=thetaf
    elif ilev==1:
        birrpdict['ilev']=ilev
        birrpdict['nout']=nout
        birrpdict['ninp']=ninp
        birrpdict['nref']=nref
        birrpdict['nrr']=nrr
        birrpdict['tbw']=tbw
        birrpdict['nfft']=nfft
        birrpdict['nsctinc']=nsctinc
        birrpdict['nsctmax']=nsctmax
        birrpdict['nf1']=nf1
        birrpdict['nfinc']=nfinc
        birrpdict['nfsect']=nfsect
        birrpdict['uin']=uin
        birrpdict['ainlin']=ainlin
        birrpdict['ainuin']=ainuin
        if nrr==1:
            birrpdict['c2threshb']=c2threshb
            birrpdict['c2threshe']=c2threshe
            if c2threshe==0 and c2threshb==0:
                birrpdict['nz']=0
            else:
                birrpdict['nz']=0
                birrpdict['perlo']=perlo
                birrpdict['perhi']=perhi
        elif nrr==0:
            birrpdict['c2threshb']=0
            birrpdict['c2threshe']=c2threshe
        birrpdict['nprej']=nprej
        birrpdict['prej']=prej
        birrpdict['c2threshe1']=c2threshe1
        birrpdict['ofil']=ofil
        birrpdict['nlev']=nlev
        birrpdict['nar']=nar
        birrpdict['imode']=imode
        birrpdict['jmode']=jmode
        birrpdict['nfil']=nfil
        birrpdict['nskip']=nskip
        birrpdict['nskipr']=nskipr
        birrpdict['thetae']=thetae
        birrpdict['thetab']=thetab
        birrpdict['thetaf']=thetaf

    print 'Made .script file: '+ofil+'.script'
    
    return scriptfile,birrpdict
    
def callBIRRP(scriptfile,birrploc):
    """
    callBIRRP(scriptfile,birrploc) will call BIRRP from a command window and 
    run the script file provided.
    
    Inputs: 
        scriptfile = full path script file
        birrploc = path to BIRRP on your computer, can be the full path to a 
                   BIRRP executable C:\BIRRP\birrp5.exe
    
    Outputs:
        bfpath = path to BIRRP output files
    """
    #if the birrploc is just a path to the directory
    if os.path.isfile(birrploc)==False:

        birrpstr=os.path.join(birrploc,'birrp5<')
        
    else:
        if birrploc.find('.exe')>=0:
            birrpstr=birrploc[:-4]+'<'
        else:
            birrpstr=birrploc
    #----Run BIRRP---
    stimeburp=time.time()
    station=os.path.basename(scriptfile)[:-7]
    print 'Started BIRRP at: ',time.ctime() 
    
    subprocess.os.system(birrpstr+'"'+scriptfile+'"')
    
    #print run time
    etimeburp=time.time()
    dtburp=etimeburp-stimeburp
    
    print 'BIRRP Run time (m:s) for ' +station+': %.0f' % int(
            np.floor(dtburp/60))+':%2f' % int(
            np.round(np.remainder(dtburp,60)))
    
    bfpath=os.path.dirname(scriptfile)
    return bfpath
    
def convertBIRRPoutputs(bfpath,stationinfofile,station,rrstation=None,
                        bbfile=None,birrpdict=None,ffactor=1,edipath=None):
    """
    convertBIRRPoutputs(bfpath,stationinfodict) will take the outputs of BIRRP
    and manipulate them into .edi, .dat, .coh, .imp files as well as generate 
    plots of apparent resistivity and phase and coherence.
    
    Inputs:
        bfpath = path to BIRRP output files
        stationinfodict = dictionary containing pertanent information with keys:
            
        station = station to process
        rrstation = remote reference station name
        bbfile = broad band calibration file path
        birrpdict = dictionary of birrp parameters used
        ffactor = scaling factor for apparent resistivities and phase
        edipath = path to save edifile to
        
    Outputs:
        cohfile
        datfile
        edifile
        impfile
    """
    
    spath=os.path.dirname(bfpath)
    count=1
    day=None
    if os.path.basename(spath)[6]=='t':
        st=os.path.basename(spath)[4:6]
    else:
        st=''
    while day==None and count<10:
        try:
            dayf=float(os.path.basename(spath))
            day=os.path.basename(spath)
        except ValueError:
            spath=os.path.dirname(spath)
            day=None
            count+=1
    if day==None:
        day=''
    
    brpexcept='BIRRP Ran incorrectly check script file'
    if rrstation==None:
        rrstation=station
    #get station info
    statdict=mt.getStationInfo(stationinfofile,station)  

    #try to write a coherence file
    try:
        cohfile=writecoh(bfpath,tsp='   ')
    except IndexError:
        print brpexcept+' did not produce .coh file'
        cohfile=' '
    
    #try to write an edifile
    try:
        edifile=writeedi(bfpath,station,stationinfofile=stationinfofile,
                         rrstation=rrstation,
                         birrpdict=birrpdict,
                         bbfile=bbfile,tsp='   ',ffactor=ffactor)
    except IndexError:
        print brpexcept+' did not find all files to produce a complete .edi file'
        edifile=' '
                         
    magtype=statdict['magtype']
    dlgain=statdict['dlgain']
    egain=statdict['egain']
    dlen=[float(statdict['ex']),float(statdict['ey'])]
    
    #try to write a dat file
    try:
        datfile=writedat(os.path.join(bfpath,station+'.j'),
                         egain=egain,dlgain=dlgain,dlen=dlen,
                         magtype=magtype,bbfile=bbfile,tsp='   ',
                         ffactor=ffactor)
    except IndexError:
        print brpexcept+' did not produce a .dat file'
        datfile=' '

    #try to write an imp file    
    try:
        impfile=writeimp(bfpath,egain=egain,dlgain=dlgain,dlen=dlen,
                         magtype=magtype,bbfile=bbfile,tsp='   ',
                         ffactor=ffactor)
    except IndexError:
        print brpexcept+' did not produce a .imp file'
        impfile=' '
    
    stationdir=os.path.dirname(bfpath)
    count=0
    while os.path.basename(stationdir).lower()!=station.lower() and count<12:
        stationdir=os.path.dirname(stationdir)
        count+=1
    
    if len(edifile)>2:    
        #make a folder to put all .edi files into
        if edipath==None:
            edipath=os.path.join(os.path.dirname(stationdir),'EDIfiles')
        #make folder if it doesn't exist
        if not os.path.exists(edipath):
            os.mkdir(edipath)
            print 'Made path: '+edipath
        
        #copy .edi file to common folder
        shutil.copy(edifile,os.path.join(edipath,station+day+st+'.edi'))
    else:
        print 'No edi file'

    return cohfile,edifile,datfile,impfile
    
def plotBFfiles(edifile,cohfile=None,cfilelst=None,save='y',
                ffactor=1,show='n',sdict=None):
    """
    plotFiles(edifile,cohfile,cfilelst) will plot the apparent resisitivty
    and phase calculated from edifile as 2 seperate plots for the different 
    polarizations.  It will plot the coherency between components and if input
    will plot the time series.
    
    Input:
        edifile = full path to edifile
        cohfile = full path to coherence file
        cfilelst = list of paths to combined files
        saveyn= 'y' to save plots, 'n' to not save plots
        seeyn = 'y' if you want to see the plots, 'n' if you just want to save
    
    Outputs:
        
    """
    
    station=os.path.basename(edifile)[:-4]
    stationdir=os.path.dirname(edifile)
    sname=[]
    day=''
    try:    
        if os.path.basename(os.path.dirname(stationdir))[6]=='t':
            st=os.path.basename(os.path.dirname(stationdir))[4:6]
        else:
            st=''
    except:
        st=''
    while os.path.basename(stationdir).lower()!=station.lower():
#        if os.path.basename(stationdir).find('to')>=0:
#            combstr=os.path.basename(stationdir)
#            tspot=combstr.find('to')
#            start=combstr[tspot-2:tspot]+'0000'
#            stop=combstr[tspot+2:tspot+4]+'0000'
        try:
            float(os.path.basename(stationdir))
            day=os.path.basename(stationdir)
        except ValueError:
            day=''
        sname.append(os.path.basename(stationdir))
        stationdir=os.path.dirname(stationdir)
    
    if save=='y':
	import mtpy.imaging.mtplottools as mtplot
        #plot the coherence and apparent resistivity and phase
        if not os.path.exists(os.path.join(os.path.dirname(stationdir),'Plots')):
            os.mkdir(os.path.join(os.path.dirname(stationdir),'Plots'))
            print 'Made path: '+os.path.join(os.path.dirname(stationdir),'Plots')
        
        plotpath=os.path.join(os.path.dirname(stationdir),'Plots')
        if cohfile!=None:
            fig1=mtplot.plotcoh(cohfile,fignum=1,
                                savefigfilename=os.path.join(plotpath,
                                                             station+day+st+
                                                             'coh.pdf'))
            print 'Save plots: '+os.path.join(plotpath,station+'coh.pdf')

        #plot and save apparent resistivity and phase
        svstr=''
        if sdict!=None:
            for key in sdict.keys():
                svstr+=key+'_'+str(sdict[key])
        
        svfn=os.path.join(plotpath,station+day+st+svstr+'Res.pdf')
        
        fig2=mtplot.plotResPhase(edifile,fignum=2,plotnum=2,ffactor=ffactor,
                                 savefigfilename=svfn)
        
        print '\t '+os.path.join(plotpath,station+day+st+svstr+'Resxy.pdf')
        print '\t '+os.path.join(plotpath,station+day+st+svstr+'Resxx.pdf')
        if cfilelst!=None:
            mtplot.plotTS(cfilelst,fignum=4,
                          savefigname=os.path.join(plotpath,station+'TS.pdf'))
    
    if show=='y':
	import mtpy.imaging.mtplottools as mtplot
        if cohfile!=None:
            mtplot.plotcoh(cohfile,fignum=1)
        mtplot.plotResPhase(edifile,fignum=2,plotnum=2,ffactor=ffactor)
        if cfilelst!=None:
            mtplot.plotTS(cfilelst,fignum=4)
    
def writeLogfile(bfpath,station,cohfile,datfile,impfile,scriptfile):
    """
    writeLogfile will write a log file of how a station was processed
    """
    
    logfile=os.path.join(bfpath,station+'.log')
    logfid=file(logfile,'a')
    todaysdate=datetime.datetime.today()
    logfid.write('==========================================================='+'\n')
    logfid.write('Processing log for station: ' +os.path.join(bfpath,station) +'\n')
    logfid.write('Processed on: ' 
                +todaysdate.strftime("%b %d, %Y at %I:%M:%S %p local time")+'\n')
    logfid.write('-----BIRRP Parameters----- \n')
    birrpfid=file(scriptfile,'r')
    birrplines=birrpfid.readlines()
    for bb in range(len(birrplines)):
        line=birrplines[bb].rstrip()
        logfid.write(line+'\n')
    birrpfid.close()
    logfid.write('-----Impedance Calculated from BIRRP-----'+'\n')
    try:    
        impfid=file(impfile,'r')
        implines=impfid.readlines()
        for mm in range(len(implines)):
            logfid.write(implines[mm].rstrip()+'\n')
        impfid.close()
    except IOError:
        pass
    
    logfid.write('-----Coherence Calculated from BIRRP-----'+'\n')
    try:    
        cohfid=file(cohfile,'r')
        cohlines=cohfid.readlines()
        for nn in range(len(cohlines)):
            logfid.write(cohlines[nn].rstrip()+'\n')
        cohfid.close()
    except IOError:
        pass
    logfid.write('-----Resistivity and Phase Calculated from BIRRP-----'+'\n')
    try:    
        datfid=file(datfile,'r')
        datlines=datfid.readlines()
        for pp in range(len(datlines)):
            logfid.write(datlines[pp].rstrip()+'\n')
        datfid.close()
    except IOError:
        pass
    logfid.close()
    
    return logfile

def runBIRRPpp(dirpath,processingdict,stationinfofile,birrploc,ffactor=1,
               edipath=None):
    """
    runBIRRPpp will processes a station from start to finish on seperate 
    processors
    
    Inputs:
        dirpath = full path to station folders
        processinginfodict = dictionary of processing info
        
        stationinfofile = full path to station info file containing info about
                          survey parameters
        birrploc = location of BIRRP and bbconv.csv
        ffactor = factor to scale apparent resistivities and phases
        edipath = common folder to copy edifiles
        
    Outputs:
        dictionary with keys:
            cohfile
            datfile
            edifile
            impfile
            logfile
            scriptfile
        
    """
    
    
    #===============================================================================
    # Run scriptfilePrep
    #===============================================================================
    prepdict=scriptfilePrep(processingdict)


            
    try:
        prepdict['magtype']=processingdict['magtype']
    except KeyError:
        pass
        
    try:
        prepdict['nout']=processingdict['nout']
    except KeyError:
        pass
    
    try:
        prepdict['ninp']=processingdict['ninp']
    except KeyError:
        pass
    
    try:
        prepdict['npcs']=processingdict['npcs']
    except KeyError:
        pass
    
    try:
        prepdict['tbw']=processingdict['tbw']
    except KeyError:
        pass
    
    try:
        prepdict['uin']=processingdict['uin']
    except KeyError:
        pass
    
    try:
        prepdict['ainuin']=processingdict['ainuin']
    except KeyError:
        pass
    
    try:
        prepdict['c2threshe']=processingdict['c2threshe']
    except KeyError:
        pass
    
    try:
        prepdict['nz']=processingdict['nz']
    except KeyError:
        pass
    
    try:
        prepdict['c2threshe1']=processingdict['c2threshe1']
    except KeyError:
        pass
    
    try:
        prepdict['nlev']=processingdict['nlev']
    except KeyError:
        pass
    
    try:
        prepdict['nar']=processingdict['nar']
    except KeyError:
        pass
        
    try:
        prepdict['imode']=processingdict['imode']
    except KeyError:
        pass
        
    try:
        prepdict['jmode']=processingdict['jmode']
    except KeyError:
        pass    
    
    try:
        prepdict['nfil']=processingdict['nfil']
    except KeyError:
        pass
    
    try:
        prepdict['nskip']=processingdict['nskip']
    except KeyError:
        pass
    
    try:
        prepdict['nskipr']=processingdict['nskipr']
    except KeyError:
        pass
    
    try:
        prepdict['thetae']=processingdict['thetae']
    except KeyError:
        pass
    
    try:
        prepdict['thetab']=processingdict['thetab']
    except KeyError:
        pass
    
    try:
        prepdict['thetaf']=processingdict['thetaf']
    except KeyError:
        pass
                         
#    print prepdict
    #===========================================================================
    # Run writeScriptfile
    #===========================================================================
    scriptfile,birrpdict=writeScriptfile(prepdict)
    #===========================================================================
    # Run callBIRRP 
    #===========================================================================

    bfpath=callBIRRP(scriptfile,birrploc)


    #===========================================================================
    # Run convertBIRRPoutputs
    #===========================================================================
    
    station=processingdict['station']
    rrstation=processingdict['rrstation']
    if birrploc.find('.')>=0:
        bbfile=os.path.join(os.path.dirname(birrploc),'BBConv.txt')
    else:
        bbfile=os.path.join(birrploc,'BBConv.txt')
    birrpdict['bbfile']=bbfile
    #print birrpdict
    cohfile,edifile,datfile,impfile=convertBIRRPoutputs(bfpath,
                                                        stationinfofile,
                                                        station,
                                                        rrstation=rrstation,
                                                        bbfile=bbfile,
                                                        birrpdict=birrpdict,
                                                        ffactor=ffactor,
                                                        edipath=edipath)
    
    #===========================================================================
    # Run writeLogfile
    #===========================================================================
    
    logfile=writeLogfile(bfpath,station,cohfile,datfile,impfile,scriptfile)
    
    returndict={}
    returndict['cohfile']=cohfile
    returndict['datfile']=datfile
    returndict['edifile']=edifile
    returndict['impfile']=impfile
    returndict['logfile']=logfile
    returndict['scriptfile']=scriptfile
    returndict['cfilelst']=prepdict['cfilelst']
    returndict['rrcfilelst']=prepdict['rrcfilelst']
    
#    #need to delete the day folder due to naming convention
#    if type(returndict['rrcfilelst'][0]) is str or\
#        type(returndict['rrcfilelst'][0]) is np.string_:
#        rrfind=returndict['rrcfilelst'][0].find('Comb')
#        try:
#            float(returndict['rrcfilelst'][0][rrfind+6])
#            cbn=os.path.dirname(returndict['rrcfilelst'][0])
#            for cfile in os.listdir(cbn):
#                os.remove(os.path.join(cbn,cfile))
#            os.rmdir(cbn)
#            print 'Removed directory: ',cbn
#        except ValueError:
#            pass
#        except WindowsError:
#            pass
#    elif type(returndict['rrcfilelst'][0]) is list:
#        for cc in range(len(returndict['rrcfilelst'])):
#            rrfind=returndict['rrcfilelst'][cc][0].find('Comb')
#            try:
#                float(returndict['rrcfilelst'][cc][0][rrfind+6])
#                cbn=os.path.dirname(returndict['rrcfilelst'][cc][0])
#                for cfile in os.listdir(cbn):
#                    os.remove(os.path.join(cbn,cfile))
#                os.rmdir(cbn)
#                print 'Removed directory: ',cbn
#            except ValueError:
#                pass
#            except WindowsError:
#                pass
    
    
    #finishBeep()
    
    return returndict


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
    
    return signum

def sil(iniline):
    """sil(iniline) will split a single line written in an .ini file
    for burpinterface and return the list of strings."""
    
    inistr=iniline.replace('\n','')
    linelst=inistr.split('=')
    
    return linelst[1]

def read2c2(filename):
    """read2c2(filename) will read in a .2c2 file output from BIRRP and return a list 
    containing [period],[freq],[coh],[zcoh]. Note if any of the coherences are 
    negative a value of 0 will be given to them."""
    period=[]
    freq=[]
    coh1=[]
    zcoh1=[]
    fid=file(filename,'r')
    fidlines=fid.readlines()
    for ii in range(len(fidlines)):
        cohline=fidlines[ii]
        cohstr=cohline.split(' ')
        cohlst=[]
        for jj in range(len(cohstr)):
            if len(cohstr[jj])>3:
                if float(cohstr[jj])<0:
                    cohlst.append(sigfigs('0.00',digits=7,fmt='f'))
                else:
                    cohlst.append(sigfigs(cohstr[jj],digits=7,fmt='f'))
        period.append(cohlst[0])
        freq.append(cohlst[1])
        coh1.append(cohlst[2])
        zcoh1.append(cohlst[4])
    return period,freq,coh1,zcoh1

def writecoh(dirpath,tsp='   '):
    """writecoh(dirpath) will write a coherence file using the BIRRP outputs 
    residing in the dirpath folder.  The output file is tab delimited and if any 
    values are negative they are put to 0. Each value has 7 significant digits.
    Returns written to filename."""

    #locate file names
    cohfilenames=[]
    for filename in os.listdir(dirpath):
        if fnmatch.fnmatch(filename,'*.1r.2c2'):
            cohfilenames.append(filename)
        elif fnmatch.fnmatch(filename,'*.2r.2c2'):
            cohfilenames.append(filename)
        elif fnmatch.fnmatch(filename,'*.3r.2c2'):
            cohfilenames.append(filename)
        else:
            pass
    ofilloc=cohfilenames[0].find('.')
    ofil=cohfilenames[0][0:ofilloc]
    if len(cohfilenames)==3:
        period,freq,coh1,zcoh1=read2c2(os.path.join(dirpath,cohfilenames[0]))
        period,freq,coh2,zcoh2=read2c2(os.path.join(dirpath,cohfilenames[1]))
        period,freq,coh3,zcoh3=read2c2(os.path.join(dirpath,cohfilenames[2]))
        cohfid=file(dirpath+os.sep+ofil+'.coh','w')
        cohfid.write(ofil+'\n')
        cohfid.write('period'+tsp+'freq'+tsp+'coh1'+tsp+'zcoh1'+tsp+'coh2'+tsp
                                           +'zcoh2'+tsp+'coh3'+tsp+'zcoh3'+'\n')
        for ff in range(len(period)):
            try:
                c1=coh1[ff]
                zc1=zcoh1[ff]
            except IndexError:
                c1='0.0000000'
                zc1='0.0000000'
            try:
                c2=coh2[ff]
                zc2=zcoh2[ff]
            except IndexError:
                c2='0.0000000'
                zc2='0.0000000'
            try:
                c3=coh3[ff]
                zc3=zcoh3[ff]
            except IndexError:
                c3='0.0000000'
                zc3='0.0000000'
            cohfid.write(period[ff]+tsp+freq[ff]+tsp+c1+tsp+zc1+tsp
                        +c2+tsp+zc2+tsp+c3+tsp+zc3+'\n')
        cohfid.close() 
    elif len(cohfilenames)==2:
        period,freq,coh1,zcoh1=read2c2(os.path.join(dirpath,cohfilenames[0]))
        period,freq,coh2,zcoh2=read2c2(os.path.join(dirpath,cohfilenames[1]))
        cohfid=file(dirpath+os.sep+ofil+'.coh','w')
        cohfid.write(ofil+'\n')
        cohfid.write('period'+tsp+'freq'+tsp+'coh1'+tsp+'zcoh1'+tsp+'coh2'+tsp
                                        +'zcoh2'+tsp+'coh3'+tsp+'zcoh3'+'\n')
        for ff in range(len(period)):
            try:
                c1=coh1[ff]
                zc1=zcoh1[ff]
            except IndexError:
                c1='0.0000000'
                zc1='0.0000000'
            try:
                c2=coh2[ff]
                zc2=zcoh2[ff]
            except IndexError:
                c2='0.0000000'
                zc2='0.0000000'
            cohfid.write(period[ff]+tsp+freq[ff]+tsp+c1+tsp+zc1+tsp
                        +c2+tsp+zc2+tsp+'0.000000'+tsp+'0.000000'+'\n')
        cohfid.close() 
    else:
        print 'Somethings wrong.  Either no coherence files or there are \r'
        'too many in the dirpath, retard'
    
    return dirpath+os.sep+ofil+'.coh'

def readcoh(filename):
    """readcoh(filename) will read a coherence file writen by writecoh.  The 
    output is period,frequency,coh1,zcoh1,coh2,zcoh2,coh3,zcoh3"""
    
    fid=file(filename,'r')
    lines=fid.readlines()
    station=lines[0].replace('\n','')
    period=[]
    freq=[]
    coh1=[]
    zcoh1=[]
    coh2=[]
    zcoh2=[]
    coh3=[]
    zcoh3=[]
    for ii in range(len(lines)-2):
        cohstr=lines[ii+2].replace('\n','')
        cohlst=cohstr.split(tsp)
        period.append(float(cohlst[0]))
        freq.append(float(cohlst[1]))
        coh1.append(float(cohlst[2]))
        zcoh1.append(float(cohlst[3]))
        coh2.append(float(cohlst[4]))
        zcoh2.append(float(cohlst[5]))
        coh3.append(float(cohlst[6]))
        zcoh3.append(float(cohlst[7]))
    
    return station,np.array(period),np.array(freq),np.array(coh1),\
        np.array(zcoh1),np.array(coh2),np.array(zcoh2),np.array(coh3),\
        np.array(zcoh3)
def bbcalfunc(bbfile,nfreqlst):
    """bbcalfunc(bbfile,nfreqlst) will generate a function fitting the 
    calibration data in bbfile to the frequencies in nfreqlst using a linear (old:cubic)
    interpolation algorithm.  The output is the real and imaginary functions
    as a function of the nfreqlst (unit depends on calibration file).

    Calibration file must contain values in 3 columns: freq, real, imag

    """ 
    
    fid=file(bbfile,'r')
    fidlines=fid.readlines()
    #define the delimiter
    if bbfile.find('.txt')>=0:
        delimiter='\t'
    elif bbfile.find('.csv')>=0:
        delimiter=','
    
    freq=[]
    breal=[]
    bimag=[]
    for ii in range(1,len(fidlines)):
        linestr=fidlines[ii]
        linestr=linestr.rstrip()
        linelst=linestr.split(delimiter)
        if len(linelst)>2:
            #frequency linear:!!
            freq.append(float(linelst[0]))
            breal.append(float(linelst[1]))
            bimag.append(float(linelst[2]))
        else:
            pass
    
    freq=np.array(freq)
    #frequencies in log space:
    logfreqs = np.log(freq)
    #linear values:
    breal=np.array(breal)
    bimag=np.array(bimag)
    #frequency of interest in log space:
    freqs_of_interest = np.array(nfreqlst)
    brep = []
    bimp = []

    for idx,f in enumerate(freqs_of_interest):
        logfreq = np.log(f)
        #find the calibration value closest to the current freq, assume it's lower
        closest_lower = np.abs(logfreq-logfreqs).argmin()
        
        #if it coincides with the highest frequency/last entry:
        if closest_lower == len(logfreqs)-1:
            brep.append(breal[-1])
            bimp.append(bimag[-1])
        # or the lowest
        elif closest_lower == 0:
            brep.append(breal[0])
            bimp.append(bimag[0])
        else:
            #in case the closest frequency value is not lower but higher, 
            #take the freq value below as lower bound for the interval:        
            if logfreqs[closest_lower] > logfreq:
                closest_lower -= 1
            
            #define the interval:
            logfreq1 = logfreqs[closest_lower]
            logfreq2 = logfreqs[closest_lower+1]

            #take the interval values:
            realval1 = breal[closest_lower]
            realval2 = breal[closest_lower+1]
            imagval1 = bimag[closest_lower]
            imagval2 = bimag[closest_lower+1]


            loginterval = logfreq2 - logfreq1
            logfreq = np.log(freq)
            weight = (logfreq2-logfreq)/loginterval

            #for low frequencies take the log of the values to get into loglog space:
            if freq <= 5:                    
                logrealval1 = np.log(realval1)
                logrealval2 = np.log(realval2)
                logimagval1 = np.log(imagval1)
                logimagval2 = np.log(imagval2)
             
                interpval_real = np.exp(weight*logrealval1 +  (1-weight) * logrealval2)
                interpval_imag = np.exp(weight*logimagval1 +  (1-weight) * logimagval2)
                

            else:
                interpval_real = weight*realval1 +  (1-weight) * realval2
                interpval_imag = weight*imagval1 +  (1-weight) * imagval2


            brep.append(interpval_real)
            bimp.append(interpval_imag)
 
    brep = np.array(brep)
    bimp = np.array(bimp)

    return brep,bip
    
def readj(jfn,egain=1,dlgain=1,dlen=[50,50],magtype='bb',
          bbfile=r'c:\Peacock\PHD\BIRRP\bbconv.txt',ffactor=1):
    """
    readj will read in the .j file output by BIRRP, which is better than 
    reading in the .irj.rf files.
    """   
    print 'Reading ',jfn    
    
    jfid=file(jfn,'rb')
    jlines=jfid.readlines()        
    jfid.close()

    tz=-1
    for ii,jline in enumerate(jlines):
        if jline.find('ZXX')==0:
            jj=ii
        elif jline.find('TZX')==0:
            tz=ii
            print 'Found Tipper'
    
    #get impedance tensor components and period
    try:
        nf=int(jlines[jj+1].strip())
    except UnboundLocalError:
        print 'Did not produce a legit .j file'
    z=np.zeros((4,3,nf))
    period=np.zeros((4,nf))
    pnot=[]
    prange=range(4)
#    for ii in range(2,10,2):
    for kk in range(4):
        for ll,nn in enumerate(range(jj+2*(kk+1)+kk*nf,jj+2*(kk+1)+(kk+1)*nf)):
            jline=jlines[nn].strip().split()
            per=float(jline[0])
            zr=float(jline[1])
            zi=float(jline[2])
            zerr=float(jline[3])
            if per!=-999:
                period[kk,ll]=float(jline[0])
            else:
                try:
                    prange.remove(kk)
                except ValueError:
                    pass
                pnot.append((kk,ll))
            #get z_real
            if zr!=-999 or zr!=float('Inf'):
                z[kk,0,ll]=zr
            else:
                pass
            
            #get z_imaginary
            if zi!=-999 or zi!=float('Inf'):
                z[kk,1,ll]=zi
            else:
                pass
            
            #get z_error
            if zerr!=-999 or zerr!=float('Inf'):
                z[kk,2,ll]=zerr
            else:
                pass
    z=np.nan_to_num(z)
    z=np.array(z)*ffactor
    #convert electric channels to microV/m
    z[0:2,:,:]=z[0:2,:,:]*float(dlgain)/(float(dlen[0])*float(egain))
    z[2:,:,:]=z[2:,:,:]*float(dlgain)/(float(dlen[1])*float(egain))

    #get tipper components
    tip=np.zeros((2,3,nf))
    if tz!=-1:
        for kk in range(2):
            for ll,nn in enumerate(range(tz+2*(kk+1)+kk*nf,tz+2*(kk+1)+(kk+1)*nf)):           
                jline=jlines[nn].strip().split()
                tr=float(jline[1])
                ti=float(jline[2])
                terr=float(jline[3])
                #get Tipper_real
                if tr!=-999 or tr!=float('Inf'):
                    tip[kk,0,ll]=tr
                else:
                    pass
                #get Tipper_real
                if ti!=-999 or ti!=float('Inf'):
                    tip[kk,1,ll]=ti
                else:
                    pass
                #get Tipper_real
                if terr!=-999 or terr!=float('Inf'):
                    tip[kk,2,ll]=terr
                else:
                    pass
    else:
        pass
    
    #find the frequency with no -999 in it
    try:
        period=period[prange[0],:]
    except IndexError:
        prange=0
        plen=len(period[0])
        for ll,per in enumerate(period[1:]):
            if len(per)>plen:
                prange=ll
        period=period[prange,:]
    if len(pnot)>0:
        nclst=['Z_xx','Z_xy','Z_yx','Z_yy']
        for mm in pnot:
            print 'No response for {1} at T={0:.3f}'.format(period[mm[1]],nclst[mm[0]])
    #flip array so frequency is decreasing, period is increasing
    freq=np.array(1./period)
    if freq[0]<freq[-1]:
        freq=freq[::-1]
        period=period[::-1]
        z=z[:,:,::-1]
        if type(tip)!=list:
            tip=tip[:,:,::-1]
    
    #convert magnetics
    if magtype.lower()=='bb':
        zconv=bbconvz(z,freq,bbfile,dlgain)
    if magtype.lower()=='lp':
        zconv=lpconvz(z,dlgain)
    
    ofil=os.path.basename(jfn)[:-4]
        
    return ofil,period,freq,zconv,tip
                

def readrf(dirpath,egain=1,dlgain=1,dlen=[1,1],magtype='bb',
           bbfile=r'c:\Peacock\PHD\BIRRP\bbconv.txt',ffactor=1):
    """readrf(dirpath) will read in .irj.rf file output by BIRRP that reside in 
    the folder nominated by dirpath. Returns: period, freq, 
    z[ijreal,ijimag,ijvar] as array.  If any of the numbers are NaN or +Inf they
    are set to zero"""
    impfn=[]
    tipfn=[]
    for filename in os.listdir(dirpath):
        if fnmatch.fnmatch(filename, '*.1r1.rf'):
            impfn.append(filename)
        elif fnmatch.fnmatch(filename, '*.1r2.rf'):
            impfn.append(filename)
        elif fnmatch.fnmatch(filename, '*.2r1.rf'):
            impfn.append(filename)
        elif fnmatch.fnmatch(filename, '*.2r2.rf'):
            impfn.append(filename)
        elif fnmatch.fnmatch(filename, '*.3r1.rf'):
            tipfn.append(filename)
        elif fnmatch.fnmatch(filename, '*.3r2.rf'):
            tipfn.append(filename)
    if len(impfn)!=4:
        raise ValueError('Somethings a miss, did not find all .rf files') 
    if len(tipfn)==2:
        print 'Got Tipper files'
    
    
    #get ofil from filenames
    ofilloc=impfn[0].find('.')
    ofil=impfn[0][0:ofilloc]
    z=[]
    period=[]
    freq=[]
    tip=[]
    for ll in range(4):
        impfid=file(os.path.join(dirpath,impfn[ll]),'r')
        implines=impfid.readlines()
        period=[]
        freq=[]
        zijreal=[]
        zijimag=[]
        zijerr=[]
        for ii in range(len(implines)):
            line=implines[ii]
            line=line.rstrip()
            impstr=line.split(' ')
            implst=[]
            for kk in range(len(impstr)):
                if len(impstr[kk])>=3:
                    if impstr=='NaN':
                        implst.append(0.0)
                    elif impstr=='+Inf':
                        implst.append(0.0)
                    else:
                        implst.append(float(impstr[kk]))
            period.append(implst[0])
            freq.append(implst[1])
            zijreal.append(implst[2])
            zijimag.append(implst[3])
            zijerr.append(implst[4])
        z.append([zijreal,zijimag,zijerr])
    try:
        z=np.array(z)*ffactor
        #convert electric channels to microV/m
        z[0:2,:,:]=z[0:2,:,:]*float(dlgain)/(float(dlen[0])*float(egain))
        z[2:,:,:]=z[2:,:,:]*float(dlgain)/(float(dlen[1])*float(egain))
    except ValueError:
        raise ValueError('BIRRP has output uneven file lengths, try running '+
                         'again with slight change in parameters.')

    #get tipper
    if len(tipfn)==2: 
        for tfn in tipfn:
            tipfid=file(os.path.join(dirpath,tfn),'r')
            tiplines=tipfid.readlines()
            tipijreal=[]
            tipijimag=[]
            tipijerr=[]
            for ii in range(len(tiplines)):
                line=tiplines[ii]
                line=line.rstrip()
                tipstr=line.split(' ')
                tiplst=[]
                for kk in range(len(tipstr)):
                    if len(tipstr[kk])>=3:
                        if tipstr=='NaN':
                            tiplst.append(0.0)
                        elif impstr=='+Inf':
                            tiplst.append(0.0)
                        else:
                            tiplst.append(float(tipstr[kk]))
                tipijreal.append(tiplst[2])
                tipijimag.append(tiplst[3])
                tipijerr.append(tiplst[4])
            tip.append([tipijreal,tipijimag,tipijerr])
        tip=np.array(tip)
        
    #flip array so frequency is decreasing, period is increasing
    freq=np.array(freq)
    period=np.array(period)
    if freq[0]<freq[-1]:
        freq=freq[::-1]
        period=period[::-1]
        z=z[:,:,::-1]
        if type(tip)!=list:
            tip=tip[:,:,::-1]
    
    #convert magnetics
    if magtype.lower()=='bb':
        zconv=bbconvz(z,freq,bbfile,dlgain)
    if magtype.lower()=='lp':
        zconv=lpconvz(z,dlgain)
        
        
    return ofil,period,freq,zconv,tip
    
def bbconvz(z,freq,bbfile,dlgain):
    """bbconvz(dirpath,bbfile,dlgain)will convert the .rf output files of BIRRP 
    into correct units using broadband calibrations given by file bbcalfile that
    has format log(freq),breal,bimag. dlgain is the data logger gain (verylow=
    2.5,low=1,high=.1) The output is period,freq,z[[zijr,ziji,zije]].
    Note that it assumes the efields are already converted."""
    
    if type(dlgain) is str:
        dlgain=float(dlgain)
    
    if type(z)!=np.ndarray:    
        z=np.array(z)
    
    if type(freq)!=np.ndarray:
        freq=np.array(freq)
    
    #interpolate values from frequency
    brep,bip=bbcalfunc(bbfile,freq)
    
    zconv=[]
    for jj in range(len(z)):
        zr=[]
        zi=[]
        ze=[]
        for ii in range(len(freq)):
            ztest=float(z[jj,0,ii])+1j*float(z[jj,1,ii])
            bcal=(brep[ii]+1j*bip[ii])
            zcal=ztest*bcal*(1./dlgain)
            zr.append(zcal.real)
            zi.append(zcal.imag)
            ze.append(float(z[jj,2,ii])*abs(bcal)*(1./dlgain))
        zconv.append([zr,zi,ze])
    
    return np.array(zconv)
    
def lpconvz(z,dlgain=1):
    """
    Convert the magnetic field from counts to units of microV/nT.
    bfield is a list of numbers. dlain is amount of gain applied
    by data logger(verylow=2.5,low=1, high=.1)
    
    Inputs:
        bfield = 1D array of long period data
        dlgain = data logger gain (very low= 2.5,low = 1, high = .1)
    
    Outputs:
        bfieldc = scaled bfield 1D array
    """
    zconv=(z*10.E7)/(70000.*float(dlgain))
    
    return zconv

def writeimp(dirpath,egain=10,dlgain=1,dlen=[50,50],magtype='bb',
           bbfile=r'c:\Peacock\PHD\BIRRP\bbconv.txt',tsp='   ',ffactor=1):
    """
    writebbimp(dirpath,magtype,dlgain,bbfile=None) writes a tab delimited .imp 
    file of converted impedances from.rf outputs of BIRRP and calibrates using 
    file bbfile and dlgain which is the data logger gain 
    (verylow=2.5,low=1,high=.1). Returns written to filename.
    """
    
    ofil,period,freq,z,tip=readrf(dirpath,egain=egain,dlgain=dlgain,
                                  dlen=dlen,magtype=magtype,bbfile=bbfile)
        
    impfid=file(os.path.join(dirpath,ofil+'.imp'),'w')
    impfid.write(ofil+'\n')
    if type(tip)!=list:
        impfid.write('period'+tsp+'rxx'+tsp+'ixx'+tsp+'rxy'+tsp+'ixy'+tsp+'ryx'
                        +tsp+'iyx'+tsp+'ryy'+tsp+'iyy'+tsp+'txxr'+tsp+'txxi'+
                        tsp+'tyyr'+tsp+'tyyi'+'\n')
    else:
        impfid.write('period'+tsp+'rxx'+tsp+'ixx'+tsp+'rxy'+tsp+'ixy'+tsp+'ryx'
                        +tsp+'iyx'+tsp+'ryy'+tsp+'iyy'+'\n')
    for ii in range(len(period)):
        per=sigfigs(str(period[ii]),digits=7,fmt='e')
        rxx=sigfigs(str(z[0,0,ii]),digits=7,fmt='e')
        ixx=sigfigs(str(z[0,1,ii]),digits=7,fmt='e')
        rxy=sigfigs(str(z[1,0,ii]),digits=7,fmt='e')
        ixy=sigfigs(str(z[1,1,ii]),digits=7,fmt='e')
        ryx=sigfigs(str(z[2,0,ii]),digits=7,fmt='e')
        iyx=sigfigs(str(z[2,1,ii]),digits=7,fmt='e')
        ryy=sigfigs(str(z[3,0,ii]),digits=7,fmt='e')
        iyy=sigfigs(str(z[3,1,ii]),digits=7,fmt='e')
        if type(tip)!=list:
            txxr=sigfigs(str(tip[0,0,ii]),digits=7,fmt='e')
            txxi=sigfigs(str(tip[0,1,ii]),digits=7,fmt='e')
            tyyr=sigfigs(str(tip[0,0,ii]),digits=7,fmt='e')
            tyyi=sigfigs(str(tip[0,1,ii]),digits=7,fmt='e')
            impfid.write(per+tsp+rxx+tsp+ixx+tsp+rxy+tsp+ixy+tsp+ryx+tsp+iyx
                            +tsp+ryy+tsp+iyy+tsp+txxr+tsp+txxi+tsp+tyyr+tsp+
                            tyyi+'\n')
        else:
            impfid.write(per+tsp+rxx+tsp+ixx+tsp+rxy+tsp+ixy+tsp+ryx+tsp+iyx
                                                        +tsp+ryy+tsp+iyy+'\n')
    impfid.close()
    #print 'Wrote .imp file to:'+os.path.join(dirpath,ofil+'.imp') 
    return os.path.join(dirpath,ofil+'.imp')

def readimp(filename):
    """readimp(filename) will read in the impedances from a .imp file written 
    by writeimp. the output is ofil,period,z[zr,zi]"""
    
    impfid=file(filename,'r')
    implines=impfid.readlines()
    ofil=implines[0][0:3]
    period=[]
    z=[]
    tip=[]
    for ii in np.arange(start=2,stop=len(implines),step=1):
        line=implines[ii].rstrip()
        line=line.split(tsp)
        period.append(float(line[0]))
        zxxr=float(line[1])
        zxxi=float(line[2])
        zxyr=float(line[3])
        zxyi=float(line[4])
        zyxr=float(line[5])
        zyxi=float(line[6])
        zyyr=float(line[7])
        zyyi=float(line[8])
        z.append([[zxxr+1j*zxxi,zxyr+1j*zxyi],[zyxr+1j*zyxi,zyyr+1j*zyyi]])
        try:
            txxr=float(line[9])
            txxi=float(line[10])
            tyyr=float(line[11])
            tyyi=float(line[12])
            tip.append([[txxr+1j*txxi],[tyyr+1j*tyyi]])
        except IndexError:
            pass
  
    return ofil,np.array(period),np.array(z),np.array(tip) 
    
def writedat(dirpath,egain=10,dlgain=1,dlen=[50,50],magtype='bb',
           bbfile=r'c:\Peacock\PHD\BIRRP\bbconv.txt',tsp='   ',ffactor=1):
    """
    writebbdat(dirpath,magtype,df,dlgain,bbfile=None) will write a .dat 
    (resistivity and phase) file from .rf output 
    files from birrp after converting the broadband magnetic channels. dirpath
    is where the .rf files reside, bb file is where the conversion file
    resides and df is the sampling frequency. Returns written to filename.
    """
    
#    ofil,period,freq,z,tip=readrf(dirpath,egain=egain,dlgain=dlgain,
#                                  dlen=dlen,magtype=magtype,bbfile=bbfile)
    ofil,period,freq,z,tip=readj(dirpath,
                                  egain=egain,dlgain=dlgain,dlen=dlen,
                                  magtype=magtype,bbfile=bbfile)
                                  
    res,phase=mt.imp2resphase(z,freq,ffactor=ffactor)
    
    r=['R11','R12','R21','R22']
    datfid=file(os.path.join(os.path.dirname(dirpath),ofil+'.dat'),'w')
    datfid.write(ofil+'\n')
    for ii in range(len(res)):
        datfid.write(r[ii]+'\n')
        datfid.write('Period'+tsp+'App Res (Ohm.m)'+tsp+'App Res Err'
                                        +tsp+'Phase(deg)'+tsp+'Phase err'+'\n')
        for jj in range(len(period)):
            datfid.write(sigfigs(str(period[jj]),digits=7,fmt='e')
                +tsp+sigfigs(str(res[ii][0][jj]),digits=7,fmt='e')+tsp
                +sigfigs(str(res[ii][1][jj]),digits=7,fmt='e')+tsp
                +sigfigs(str(phase[ii][0][jj]),digits=7,fmt='e')
                +tsp+sigfigs(str(phase[ii][1][jj]),digits=7,fmt='e')+'\n')
                
#    res=[[resxx,resxxerr],[resxy,resxyerr],[resyx,resyxerr],[resyy,resyyerr]]
#    phase=[[phasexx,phasexxerr],[phasexy,phasexyerr],[phaseyx,phaseyxerr],
#           [phaseyy,phaseyyerr]]
    
    datfid.close()
    return os.path.join(dirpath,ofil+'.dat')

def readdat(filename):
    """readdat(filename) will read in a .dat file written by writedat and output 
    ofil, period,[resistivity,resistivityerr],[phase,phaseerr]"""

    datfid=file(filename,'r')
    datlines=datfid.readlines()
    ofil=datlines[0][0:3]
    rfind=[]
    for ii in np.arange(start=1,stop=len(datlines),step=1):
        if datlines[ii].find('R',0,2)>=0:
            rfind.append(ii)
    numper=rfind[1]-rfind[0]-2
    numcomp=int(len(datlines)/numper)
    if numcomp!=4:
        raise ValueError('Not right amoung of components in .dat file, check file')

    res=[]
    reserr=[]
    phase=[]
    phaseerr=[]
    for ii in range(numcomp):
        resij=[]
        resijerr=[]
        phaseij=[]
        phaseijerr=[]
        period=[]
        for jj in np.arange(start=rfind[ii]+2,stop=rfind[ii]+numper+2,step=1):
            line=datlines[jj].replace('\n','')
            line=line.split(tsp)
            period.append(float(line[0]))
            resij.append(float(line[1]))
            resijerr.append(float(line[2]))
            phaseij.append(float(line[3]))
            phaseijerr.append(float(line[4]))
        res.append(resij)
        reserr.append(resijerr)
        phase.append(phaseij)
        phaseerr.append(phaseijerr)
    return ofil,period,[res,reserr],[phase,phaseerr]
    
def writeedi(dirpath,station,stationinfofile=None,rrstation=None,
             birrpdict=None,bbfile=r'c:\Peacock\PHD\BIRRP\bbconv.txt',
             tsp='   ',ffactor=1):
    """
    writeedi2(dirpath,stationinfofile,station,bbfile=
    r'c:\Peacock\PHD\BIRRP\bbconv.txt') will write an .edi file for a station 
    processed by BIRRP given the station info file, station and bbfile if
    applicable. Returns the full path to the .edi file.
    """
    
    #long spaces 6 spaces
    lsp=2*tsp
    
    #get station info from station info file
    if stationinfofile==None or stationinfofile=='None':
        statdict={}
    else:
        statdict=mt.getStationInfo(stationinfofile,station)
    
    #convert outputs of BIRRP to impedances depending on magtype
    try: 
        magtype=statdict['magtype'].lower()
    except KeyError:
        magtype=input('Enter magtype bb for broadband, lp for longperiod. '+\
                    'be sure to put quotes around answer. \n')
        statdict['lat']=input('Enter latitude. \n')
        statdict['long']=input('Enter longitude. \n')
        statdict['ex']=str(input('Enter Ex length (m). \n'))
        statdict['ey']=str(input('Enter Ey length (m). \n'))
        statdict['elev']=str(input('Enter Elevation (m). \n'))
        statdict['egain']=str(input('Enter interface box gain (10 or 100). \n'))
    try:
        dlgain=statdict['dlgain']
    except KeyError:
        dlgain=input('Enter data logger gain (verylow=2.5,low=1,high=.1)'+\
                    'be sure to put quotes around answer. \n')
        
#    ofil,period,freq,z,tip=readrf(dirpath,
#                                  egain=float(statdict['egain']),
#                                  dlgain=float(dlgain),
#                                  dlen=[float(statdict['ex']),
#                                        float(statdict['ey'])],
#                                  magtype=magtype,
#                                  bbfile=bbfile,ffactor=ffactor)
    ofil,period,freq,z,tip=readj(os.path.join(dirpath,station+'.j'),
                                  egain=float(statdict['egain']),
                                  dlgain=float(dlgain),
                                  dlen=[float(statdict['ex']),
                                        float(statdict['ey'])],
                                  magtype=magtype,
                                  bbfile=bbfile,ffactor=ffactor)
    
    if freq[0]<freq[-1]:
        freq=freq[::-1]
        period=period[::-1]
        z=z[:,:,::-1]
        if type(tip)!=list:
            tip=tip[:,:,::-1]
        print 'Flipped array so frequency is decending'
                    
    nfreq=str(len(freq))
        
    #list of orientation components
    orilst=['HX','HY','EX','EY','RX','RY']
    
    #open an edifile to write to
    if ofil!=station:
        ofil=station
    edifid=file(os.path.join(dirpath,ofil+'.edi'),'w')
    
    #---------------------------------------------------------------------------
    #write header information
    edifid.write('>HEAD \n')
    edifid.write(tsp+'DATAID="'+ofil+'"'+'\n')
    
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
        elif len(mdate[2])>4:
            myear=int(mdate[2][0:4])
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
        edifid.write(tsp+'LOC="'+ofil+'"'+'\n')
    
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
        edifid.write(lsp+'Coil Calibration File='+bbfile+'\n')
        try:
            edifid.write(lsp+'Interaction Level (ILEV)='+str(birrpdict['ilev'])+
                        '\n')
            ilevyn=int(birrpdict['ilev'])
        except KeyError:
            pass
        try:
            edifid.write(lsp+'Number of outputs (NOUT)='+str(birrpdict['nout'])+
                       '\n')
        except KeyError:
            pass
        try:
            edifid.write(lsp+'Number of inputs (NINP)='+str(birrpdict['ninp'])+
                        '\n')
        except KeyError:
            pass
        if ilevyn==1:
            try:
                edifid.write(lsp+'Number of Remote Reference time series (NREF)='
                             +str(birrpdict['nref'])+'\n')
            except KeyError:
                pass
            
            try:
                edifid.write(lsp+'Remote reference(0) or bounded influence(1)'+
                            '(NRR)='+str(birrpdict['nrr'])+'\n')
            except KeyError:
                pass
        try:
            edifid.write(lsp+'Slepian Filter order (TBW)='+str(birrpdict['tbw'])+
                        '\n')
        except KeyError:
            pass
        try:
            edifid.write(lsp+'Max length of fft window (NFFT)='+
                        str(birrpdict['nfft'])+'\n')
        except KeyError:
            pass
        if ilevyn==1:
            try:
                edifid.write(lsp+'Section Increment divisor (NSCTINC)='
                             +str(birrpdict['nsctinc'])+'\n')
            except KeyError:
                pass
        
        try:
            edifid.write(lsp+'Maximum number of fft sections (NSCTMAX)='+
                        str(birrpdict['nsctmax'])+'\n')
        except KeyError:
            pass
        
        if ilevyn==1:
            try:
                edifid.write(lsp+'First frequency extracted (NF1)='
                             +str(birrpdict['nf1'])+'\n')
            except KeyError:
                pass
            try:
                edifid.write(lsp+'Frequency increment per window (NFINC)='
                             +str(birrpdict['nfinc'])+'\n')
            except KeyError:
                pass
            try:
                edifid.write(lsp+'Number of frequencies per window (NFSECT)='
                             +str(birrpdict['nf1'])+'\n')
            except KeyError:
                pass
        try:
            edifid.write(lsp+'Small leverage point control (UIN)='+
                        str(birrpdict['uin'])+'\n')
        except KeyError:
            pass
        if ilevyn==1:
            try:
                edifid.write(lsp+'Lower leverage point control (AINLIN)='
                             +str(birrpdict['ainlin'])+'\n')
            except KeyError:
                pass
        try:
            edifid.write(lsp+'Large leverage point control (AINUIN)='+
                        str(birrpdict['ainuin'])+'\n')
        except KeyError:
            pass
        
        if ilevyn==1:
            try:
                edifid.write(lsp+'Magnetic coherence threshold (C2THRESHEB)='
                             +str(birrpdict['c2threshb'])+'\n')
            except KeyError:
                pass
        try:
            edifid.write(lsp+'Electric coherence threshold (C2THRESHE)='+
                        str(birrpdict['c2threshe'])+'\n')
        except KeyError:
            pass
        try:
            edifid.write(lsp+'Z component (NZ)='+str(birrpdict['nz'])+'\n')
        except KeyError:
            pass
        
        try:  
            edifid.write(lsp+'Coherence threshold z channel (c2threshe1)='+
                        str(birrpdict['c2threshe1'])+'\n')
        except KeyError:
            pass
        if ilevyn==1:
            try:
                edifid.write(lsp+'Low and high periods for coherence threshold '
                            +'(PERLO,PERHI)='+str(birrpdict['perlo'])+','+
                            str(birrpdict['perhi'])+'\n')
            except KeyError:
                pass
            try:
                edifid.write(lsp+'Number of periods to reject (NPREJ)='
                             +str(birrpdict['nprej'])+'\n')
            except KeyError:
                pass
            try:
                edifid.write(lsp+'Periods to reject (PREJ)='
                             +str(birrpdict['prej'])+'\n')
            except KeyError:
                pass
        try:
            edifid.write(lsp+'Order of prewhitening filter (NAR)='+
                        str(birrpdict['nar'])+'\n')
        except KeyError:
            pass
        try:
            edifid.write(lsp+'Electric channel rotation angles (THETAE)='+
                        birrpdict['thetae']+'\n')
        except KeyError:
            pass
        try:
            edifid.write(lsp+'Magnetic channel rotation angles (THETAB)='+
                        birrpdict['thetab']+'\n')
        except KeyError:
            pass
        try:
            edifid.write(lsp+'Final channel rotation angles (THETAF)='+
                        birrpdict['thetaf']+'\n')
        except KeyError:
            pass
            
    #remote reference
    if rrstation!=None:
        rrstation=rrstation.replace(';',',')
        rrstationlst=rrstation.split(',')
    else:
        rrstationlst=[station]
    if len(rrstationlst)<=1:
        rrstation=rrstationlst[0]
        if rrstation!=None and rrstation!=station:
            if stationinfofile==None or stationinfofile=='None':
                pass
            else:
                rrdict=mt.getStationInfo(stationinfofile,rrstation)
                edifid.write(lsp+'Remote Reference Station: '+rrstation+'\n')
                edifid.write(lsp+'Remote Reference Lat='
                                            +'%2.8g' % float(rrdict['lat'])+'\n')
                edifid.write(lsp+'Remote Reference Long='
                                            +'%2.8g' % float(rrdict['long'])+'\n')
                edifid.write(lsp+'Remote Reference Elev='+rrdict['elev']+'\n')
        else:
            edifid.write(lsp+'Remote Reference Station: '+station+'\n')
            edifid.write(lsp+'Remote Reference Lat='
                                        +'%2.8g' % float(statdict['lat'])+'\n')
            edifid.write(lsp+'Remote Reference Long='
                                        +'%2.8g' % float(statdict['long'])+'\n')
            edifid.write(lsp+'Remote Reference Elev='+statdict['elev']+'\n')
    else:
        for rrs in rrstationlst:
            rfind=np.where(np.array(rrstationlst)==rrs)[0]
            if len(rfind)>1:
                for rf in range(len(rfind)):
                    try:
                        rrstationlst.__delitem__(rfind[rf])
                    except IndexError:
                        break
            if rrs!=station:
                if stationinfofile==None or stationinfofile=='None':
                    pass
                else:
                    rrdict=mt.getStationInfo(stationinfofile,rrs)
                    edifid.write(lsp+'Remote Reference Station: '+rrs+'\n')
                    edifid.write(lsp+'Remote Reference Lat='
                                                +'%2.8g' % float(rrdict['lat'])+'\n')
                    edifid.write(lsp+'Remote Reference Long='
                                                +'%2.8g' % float(rrdict['long'])+'\n')
                    edifid.write(lsp+'Remote Reference Elev='+rrdict['elev']+'\n')
            else:
                pass
    
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
    edifid.write(tsp+'SECTID='+ofil+'\n')
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
            znum='{0:+.6e}'.format(z[mm,nn,kk])
            if znum.find('INF')>=0:
                znum='{0:+.6e}'.format(-6.666)  
            edifid.write(tsp+znum)
            if np.remainder(float(kk)+1,5.)==0:
                edifid.write('\n')
        edifid.write('\n')
    edifid.write('\n')
    
    #---------------------------------------------------------------------------
    #write tipper info
    
    edifid.write('>!****TIPPER****!'+'\n')
    tiplst=[['TXR',0,0],['TXI',0,1],['TX.VAR',0,2],['TYR',1,0],['TYI',1,1],
            ['TY.VAR',1,2]]
    if len(tip)==0:
        tip=np.zeros((2,3,float(nfreq)))
        ntip=int(nfreq)
    else:
        tip=np.array(tip)
        ntip=tip.shape[2]
    for jj,tcomp in enumerate(tiplst):
        mm=tcomp[1]
        nn=tcomp[2]
        edifid.write('>'+tcomp[0]+' // '+str(ntip)+'\n')
        for kk in range(int(ntip)):
            tipnum='{0:+.6e}'.format(tip[mm,nn,kk])
            if tipnum.find('INF')>=0:
                znum='{0:+.6e}'.format(-6.666)   
            edifid.write(tsp+tipnum)
            if np.remainder(float(kk)+1,5.)==0:
                edifid.write('\n')
        edifid.write('\n')
    edifid.write('\n')
    edifid.write('>END')
    edifid.close()
    
    return os.path.join(dirpath,ofil+'.edi')
    
def readini(inifilename):
    """readini(inifilename) will read in an inifile and return a dictionary of 
    initial parameters for the burpinterface program."""
    
    inifid=file(inifilename,'r')
    inilines=inifid.readlines()
    inidict={}

    if len(inilines)==38:
        inidict['defdirpath']=sil(inilines[0])
        inidict['station']=sil(inilines[1])
        inidict['magtype']=sil(inilines[2])
        inidict['lpyn']=sil(inilines[3])
        inidict['eyn']=sil(inilines[4])
        inidict['mcomps']=sil(inilines[5])
        inidict['magdec']=sil(inilines[6])
        inidict['df']=sil(inilines[7])
        inidict['cacherate']=sil(inilines[8])
        inidict['dlength']=sil(inilines[9])
        inidict['dlgain']=sil(inilines[10])
        inidict['egain']=sil(inilines[11])
        inidict['lpbzcor']=sil(inilines[12])
        inidict['bbfile']=sil(inilines[13])
        inidict['magori']=sil(inilines[14])
        inidict['birrploc']=sil(inilines[15])
        inidict['ilev']=sil(inilines[16])
        inidict['nout']=sil(inilines[17])
        inidict['ninp']=sil(inilines[18])
        inidict['tbw']=sil(inilines[19])
        inidict['nfft']=sil(inilines[20])
        inidict['nsctmax']=sil(inilines[21])
        inidict['uin']=sil(inilines[22])
        inidict['ainuin']=sil(inilines[23])
        inidict['c2threshe']=sil(inilines[24])
        inidict['nz']=sil(inilines[25])
        inidict['c2threshe1']=sil(inilines[26])
        inidict['ofil']=sil(inilines[27])
        inidict['nlev']=sil(inilines[28])
        inidict['nar']=sil(inilines[29])
        inidict['imode']=sil(inilines[30])
        inidict['jmode']=sil(inilines[31])
        inidict['nfil']=sil(inilines[32])
        inidict['complstr']=sil(inilines[33])
        inidict['thetae']=sil(inilines[34])
        inidict['thetab']=sil(inilines[35])
        inidict['thetaf']=sil(inilines[36])
        inidict['stationinfofile']=sil(inilines[37])
    elif len(inilines)==37:
        inidict['defdirpath']=sil(inilines[0])
        inidict['station']=sil(inilines[1])
        inidict['magtype']=sil(inilines[2])
        inidict['lpyn']=sil(inilines[3])
        inidict['eyn']=sil(inilines[4])
        inidict['mcomps']=sil(inilines[5])
        inidict['magdec']=sil(inilines[6])
        inidict['df']=sil(inilines[7])
        inidict['cacherate']=sil(inilines[8])
        inidict['dlength']=sil(inilines[9])
        inidict['dlgain']=sil(inilines[10])
        inidict['egain']=sil(inilines[11])
        inidict['lpbzcor']=sil(inilines[12])
        inidict['bbfile']=sil(inilines[13])
        inidict['magori']=sil(inilines[14])
        inidict['birrploc']=sil(inilines[15])
        inidict['ilev']=sil(inilines[16])
        inidict['nout']=sil(inilines[17])
        inidict['ninp']=sil(inilines[18])
        inidict['tbw']=sil(inilines[19])
        inidict['nfft']=sil(inilines[20])
        inidict['nsctmax']=sil(inilines[21])
        inidict['uin']=sil(inilines[22])
        inidict['ainuin']=sil(inilines[23])
        inidict['c2threshe']=sil(inilines[24])
        inidict['nz']=sil(inilines[25])
        inidict['c2threshe1']=sil(inilines[26])
        inidict['ofil']=sil(inilines[27])
        inidict['nlev']=sil(inilines[28])
        inidict['nar']=sil(inilines[29])
        inidict['imode']=sil(inilines[30])
        inidict['jmode']=sil(inilines[31])
        inidict['nfil']=sil(inilines[32])
        inidict['complstr']=sil(inilines[33])
        inidict['thetae']=sil(inilines[34])
        inidict['thetab']=sil(inilines[35])
        inidict['thetaf']=sil(inilines[36])
        inidict['stationinfofile']=None
    else:
        raise ValueError('Length of .ini file is not correct (len= %.3g' % 
                len(inilines)+') check to make sure all values are entered')
          
    return inidict
    
def writeini(filepath,argsdict):
    """writeini(filepath,station,argsdict) will write an .ini file from the to
    the filepath as station.ini.  The argsdict must be len 40 or will not write,
    the variables should be: dict[defdirpath,station,magtype,lpyn,eyn,mcomps,magdec,
    df,cacherate, dlength,dlgain,egain,lpbzcor,bbcal,magori,birrploc,ilev,nout,
    ninp,tbw,nfft,nsctmax,uin,ainuin,c2threshe,nz,c2threshe1,ofil,
    nlev,nar,imode,jmode,nfil,complstr,thetae,thetab,thetaf].""" 

    if len(argsdict)!=38:
        raise ValueError('Length of argsdict is not correct (len= %.3g' % 
                len(argsdict)+') check to make sure all values are entered.')
    else:
        inifid=file(os.path.join(filepath,argsdict['station']+'.ini'),'w')
        inifid.write('Default directory path='+argsdict['defdirpath']+'\n')
        inifid.write('Station='+argsdict['station']+'\n')
        inifid.write('Magnetic type='+argsdict['magtype']+'\n')
        inifid.write('Long period magnetic fields converted to microV/nT='
                                                        +argsdict['lpyn']+'\n')
        inifid.write('Electric fields converted to microV/m='+
                                                        argsdict['eyn']+'\n')
        inifid.write('Measurement components in order for BIRRP='+
                                                        argsdict['mcomps']+'\n')
        inifid.write('Magnetic declination='+argsdict['magdec']+'\n')
        inifid.write('Sampling frequency (Hz)='+argsdict['df']+'\n')
        inifid.write('Cache rate (hhmmss)='+argsdict['cacherate']+'\n')
        inifid.write('Dipole lengths Ex,Ey (m)='+argsdict['dlength']+'\n')
        inifid.write('Data logger gain='+argsdict['dlgain']+'\n')
        inifid.write('Interface box gain='+argsdict['egain']+'\n')
        inifid.write('Long period Bz correction='+argsdict['lpbzcor']+'\n')
        inifid.write('Broadband calibration file='+argsdict['bbfile']+'\n')
        inifid.write('Magnectic measured orientation relative to Cartesian '\
                                        +'Bx,By,Bz='+argsdict['magori']+'\n')
        inifid.write('BIRRP location='+argsdict['birrploc']+'\n')
        inifid.write('Interaction Level='+argsdict['ilev']+'\n')
        inifid.write('Number of output time series='+argsdict['nout']+'\n')
        inifid.write('Number of input time series='+argsdict['ninp']+'\n')
        inifid.write('Number of reference timeseries='+argsdict['tbw']+'\n')
        inifid.write('Length of FFT window='+argsdict['nfft']+'\n')
        inifid.write('Number of windows used in FFT='+argsdict['nsctmax']+'\n')
        inifid.write('Residual rejection factor low end (usually 0)='
                                                        +argsdict['uin']+'\n')
        inifid.write('Residual rejection factor high end (.95-.99)='
                                                    +argsdict['ainuin']+'\n')
        inifid.write('Coherence threshold (0 if not desired)='
                                                    +argsdict['c2threshe']+'\n')
        inifid.write('Threshold for Bz='+argsdict['nz']+'\n')
        inifid.write('Squared coherence for Bz='+argsdict['c2threshe1']+'\n')
        inifid.write('Output file root='+argsdict['ofil']+'\n')
        inifid.write('Output files mode='+argsdict['nlev']+'\n')
        inifid.write('Prewhitening filter (3< >15) or 0 if not desired='
                                                        +argsdict['nar']+'\n')
        inifid.write('Output file mode='+argsdict['imode']+'\n')
        inifid.write('Input file mode='+argsdict['jmode']+'\n')
        inifid.write('Filter parameters='+argsdict['nfil']+'\n')
        inifid.write('Remote Reference component list (in order)='+argsdict['complstr']+'\n')
        inifid.write('Rotation angles for electrics (relative to geomagnetic '\
                                    +'North)(N,E,rot)='+argsdict['thetae']+'\n')
        inifid.write('Rotation angles for magnetics (relative to geomagnetic '\
                                    +'North)(N,E,rot)='+argsdict['thetab']+'\n')
        inifid.write('Rotation angles for calculation (relative to geomagnetic'\
                                    +'North)(N,E,rot)='+argsdict['thetaf']+'\n')
        inifid.write('Station Info File='+argsdict['stationinfofile']+'\n')
        inifid.close()
