# -*- coding: utf-8 -*-
"""
Created on Mon Apr 02 11:54:33 2012

@author: a1185872
"""

import os
import numpy as np
import Z
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
from matplotlib.patches import Ellipse,Rectangle,Arrow
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.colorbar import *
import matplotlib.gridspec as gridspec


#tolerance to find frequencies
ptol=.15

#error of data in percentage
zerr=.05
#errormap values which is multiplied by zerr to get a total error
zxxerrmap=10
zxyerrmap=1
zyxerrmap=1
zyyerrmap=10
zerrmap=[zxxerrmap,zxyerrmap,zyxerrmap,zyyerrmap]

#==============================================================================
# Colormaps for plots
#==============================================================================
#phase tensor map
ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmap=LinearSegmentedColormap('ptcmap',ptcmapdict,256)

#phase tensor map for difference (reverse)
ptcmapdictr={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,1.0,0.0),(1.0,1.0,0.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmapr=LinearSegmentedColormap('ptcmapr',ptcmapdictr,256)

#resistivity tensor map for calculating delta
ptcmapdict2={'red':((0.0,1.0,0.0),(1.0,1.0,0.0)),
            'green':((0.0,0.5,0.5),(1.0,0.5,0.5)),
            'blue':((0.0,0.5,0.5),(1.0,0.5,0.5))}
ptcmap2=LinearSegmentedColormap('ptcmap2',ptcmapdict2,256)

#resistivity tensor map for calcluating resistivity difference
rtcmapdict={'red':((0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,1.0,0.0)),
            'green':((0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,0.0,0.0)),
            'blue':((0.0,0.0,1.0),(0.5,1.0,1.0),(1.0,0.0,0.0))}
rtcmap=LinearSegmentedColormap('rtcmap',rtcmapdict,256)

#resistivity tensor map for calcluating apparent resistivity
rtcmapdictr={'red':((0.0,1.0,1.0),(0.5,1.0,1.0),(1.0,0.0,0.0)),
            'green':((0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,0.0,0.0)),
            'blue':((0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,1.0,1.0))}
rtcmapr=LinearSegmentedColormap('rtcmapr',rtcmapdictr,256)

#==============================================================================
#  define some helping functions
#==============================================================================
#make a class to pick periods
class ListPeriods:
    def __init__(self,fig):
        self.plst=[]
        self.fig=fig
        self.count=1
    def connect(self):
        self.cid=self.fig.canvas.mpl_connect('button_press_event',
                                                self.onclick)
    def onclick(self,event):
        print '{0} Period: {1:.5g}'.format(self.count,event.xdata)
        self.plst.append(event.xdata)
        self.count+=1

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid)

def readWLOutFile(outfn,ncol=5):
    """
    read .out file from winglink
    
    Inputs:
        outfn = full path to .out file from winglink
        
    Outputs:
        dx,dy,dz = cell nodes in x,y,z directions (note x is to the East here
                    and y is to the north.)
    """
    
    rfid=file(outfn,'r')
    rlines=rfid.readlines()
    nx,ny,nz,rhostart,values=rlines[0].strip().split()
    nx=int(nx)
    ny=int(ny)
    nz=int(nz)
    dx=np.zeros(nx)
    dy=np.zeros(ny)
    dz=np.zeros(nz)
    
    #indices to read from 
    rx=int(round(nx/float(ncol))+1)
    ry=rx+int(np.ceil(ny/float(ncol)))
    rz=ry+int(np.ceil(nz/float(ncol)))
    
#    print 'rx:{0}; ry:{1}; rz:{2}'.format(rx,ry,rz)
    
    #get x distances
    for ii,ll in enumerate(rlines[1:rx]):
        line=ll.strip().split()
        for jj,kk in enumerate(line):
            dx[ii*ncol+jj]=float(kk)
    #get y distances
    for ii,ll in enumerate(rlines[rx:ry]):
        line=ll.strip().split()
#        print line
        for jj,kk in enumerate(line):
            dy[ii*ncol+jj]=float(kk)
    
    #get z distances
    for ii,ll in enumerate(rlines[ry:rz]):
        line=ll.strip().split()
        for jj,kk in enumerate(line):
            dz[ii*ncol+jj]=float(kk)
      
    #make cells to the west and south negative      
    dx[0:nx/2]=-dx[0:nx/2]
    dy[0:ny/2]=-dy[0:ny/2]
    dz[0:nz/2]=-dz[0:nz/2]
            
    return dx,dy,dz
    
def readSitesFile(sitesfn):
    """
    read sites_ file output from winglink
    
    Input: 
        sitesfn = full path to the sites file output by winglink
        
    Output:
        slst = list of dictionaries for each station.  Keys include:
            station = station name
            dx = number of blocks from center of grid in East-West direction
            dy = number of blocks from center of grid in North-South direction
            dz = number of blocks from center of grid vertically
            number = block number in the grid
        sitelst = list of station names 
    """
    
    sfid=file(sitesfn,'r')
    slines=sfid.readlines()
    
    slst=[]
    sitelst=[]
    for ss in slines:
        sdict={}
        sline=ss.strip().split()
        sdict['station']=sline[0][0:-4]
        sdict['dx']=int(sline[1])-1
        sdict['dy']=int(sline[2])-1
        sdict['dz']=int(sline[3])-1
        sdict['something']=int(sline[4])
        sdict['number']=int(sline[5])
        slst.append(sdict)
        sitelst.append(sline[0][0:-4])
    return slst,sitelst
    
def getXY(sitesfn,outfn,ncol=5):
    """
    get x (e-w) and y (n-s) position of station and put in middle of cell
    
    Input:
        sitesfn = full path to sites file output from winglink
        outfn = full path to .out file output from winglink
        ncol = number of columns the data is in
        
    Outputs:
        xarr = array of relative distance for each station from center of the
                grid.  Note this is E-W direction
        yarr = array of relative distance for each station from center of the
                grid.  Note this is N-S direction
                
    """
    
    slst,sitelst=readSitesFile(sitesfn)
    
    dx,dy,dz=readWLOutFile(outfn,ncol=ncol)
    
    ns=len(slst)
    nxh=len(dx)/2
    nyh=len(dy)/2
    xarr=np.zeros(ns)
    yarr=np.zeros(ns)
    
    
    for ii,sdict in enumerate(slst):
        xx=sdict['dx']
        yy=sdict['dy']
        if xx<nxh:
            xarr[ii]=dx[xx:nxh].sum()-dx[xx]/2
        else:
            xarr[ii]=dx[nxh:xx].sum()+dx[xx]/2                    
        if yy<nyh:
            yarr[ii]=-1*(dy[yy:nyh].sum()-dy[yy]/2)
        else:
            yarr[ii]=-1*(dy[nyh:yy].sum()+dy[yy]/2)   

    return xarr,yarr  

def getPeriods(edipath,errthresh=10):
    """
    Plots periods for all stations in edipath and the plot is interactive, just
    clikc on the period you want to select and it will appear in the console,
    it will also be saved to lp.plst.  To sort this list type lp.plst.sort()
    
    The x's mark a conformation that the station contains that period.  So 
    when looking for the best periods to invert for look for a dense line of 
    x's
    
    Inputs:
        edipath = path to where all your edi files are.  Note that only the 
            impedance components are supported so if you have spectra data, 
            export them from wingling to have impedance information.
        errthresh = threshold on the error in impedance estimation, this just 
                    gives an indication on where bad stations and bad periods
                    are, anything above this level will be colored in red.
    
    Outputs:
        periodlst = list of periods for each station
        errorlst = error in the impedance determinant for each station at 
                   each period.
        lp = data type lp has attributes: 
            plst = period list of chosen periods, again to sort this list type
                    lp.plst.sort().  this will then be the input to make the 
                    data file later.
        
    """
    
    plt.rcParams['font.size']=10
    plt.rcParams['figure.subplot.left']=.13
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.1
    plt.rcParams['figure.subplot.top']=.95
    plt.rcParams['figure.subplot.wspace']=.25
    plt.rcParams['figure.subplot.hspace']=.05    
    
    periodlst=[]
    errorlst=[]
    
    fig1=plt.figure(5)
    ax=fig1.add_subplot(1,1,1)
    for edi in os.listdir(edipath):
        if edi.find('.edi')>0:
            z1=Z.Z(os.path.join(edipath,edi))
            periodlst.append(z1.period)
            zdet=np.array([np.sqrt(abs(np.linalg.det(zz))) for zz in z1.z])
            error=np.array([np.sqrt(abs(np.linalg.det(zz))) for zz in z1.zvar])
            perror=(error/zdet)*100            
            errorlst.append(perror)
            #make a plot to pick frequencies from showing period and percent 
            #error
            ax.scatter(z1.period,perror,marker='x',picker=5)
            pfind=np.where(perror>errthresh)[0]
            if len(pfind)>0: 
                print 'Error greater than {0:.3f} for '.format(errthresh)+z1.station
                for jj in pfind:
                    ax.scatter(z1.period[jj],perror[jj],marker='x',color='r')
                    ax.text(z1.period[jj],perror[jj]*1.05,z1.station,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size':8,'color':'red'})
                    print jj,z1.period[jj]
                    
    ax.set_xscale('log')
    ax.set_xlim(10**np.floor(np.log10(z1.period[0])),
                10**np.ceil(np.log10(z1.period[-1])))
    ax.set_ylim(0,3*errthresh)
    ax.set_yscale('log')
    ax.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})
    ax.set_ylabel('Percent Error',fontdict={'size':12,'weight':'bold'})
    ax.grid('on',which='both')
    
    lp=ListPeriods(fig1)
    lp.connect()
        
    return periodlst,errorlst,lp
    
            
    
    
def writeWSDataFile(sitesfn,outfn,periodlst,edipath,zerr=.05,ptol=.15,
                    zerrmap=[10,1,1,10],savepath=None,ncol=5,units='mv'):
    """
    writes a data file for WSINV3D from winglink outputs
    
    Inputs:
        sitesfn = sites filename (full path)
        outfn = winglink .out file
        periodlst = periods to extract from edifiles can get them from 
                    using the function getPeriods.
        edipath = path to edifiles
        savepath = directory or full path to save data file to, default path
                    is dirname sitesfn.  saves as: savepath/WSDataFile.dat
        zerr = percent error to give to impedance tensor components
        ptol = percent tolerance to locate frequencies in case edi files don't
                have the same frequencies.  Need to add interpolation
        zerrmap = multiple to multiply err of zxx,zxy,zyx,zyy by.
                  Note the total error is zerr*zerrmap[ii]
        ncol = number of columns in outfn
        
    
    Outputs:
        datafn = full path to data file, saved in dirname(sitesfn) or savepath
                where savepath can be a directory or full filename
    """
    
    #get units correctly
    if units=='mv':
        zconv=1./796.

    #create the output filename
    if savepath==None:
        ofile=os.path.join(os.path.dirname(sitesfn),'WSDataFile.dat')
    elif savepath.find('.')==-1:
        ofile=os.path.join(savepath,'WSDataFile.dat')
    else:
        ofile=savepath
        
    #read in stations from sites file
    sitelst,slst=readSitesFile(sitesfn)
    
    #get x and y locations on a relative grid
    xlst,ylst=getXY(sitesfn,outfn,ncol=ncol)
    
    #define some lengths
    ns=len(slst)
    nperiod=len(periodlst)
    
    #make an array to put data into for easy writing
    zarr=np.zeros((ns,nperiod,4),dtype='complex')
    
    #--------find frequencies---------------------------------------------------
    linelst=[]
    for ss,s1 in enumerate(slst):
        edifn=os.path.join(edipath,s1+'.edi')
        #find the edi file if not named the same as in the site file
        count=1
        while os.path.isfile(edifn)==False:
            edifn=os.path.join(edipath,s1[0:-count]+'.edi')
            count+=1
            if count>6:
                raise IOError('Could not find an .edi file for:'+s1)
    
        z1=Z.Z(os.path.join(edifn))
        if s1[0:len(z1.station)]==z1.station:
            sdict={}
            fspot={}
            for ff,f1 in enumerate(periodlst):
                for kk,f2 in enumerate(z1.period):
                    if f2>=(1-ptol)*f1 and f2<=(1+ptol)*f1:
                        zderr=np.array([abs(z1.zvar[kk,nn,mm])/
                                        abs(z1.z[kk,nn,mm])*100 
                                        for nn in range(2) for mm in range(2)])
                        fspot['{0:.6g}'.format(f1)]=(kk,f2,zderr[0],zderr[1],
                                                      zderr[2],zderr[3])
                        zarr[ss,ff,:]=z1.z[kk].reshape(4,)
                        
            print z1.station, len(fspot)
            sdict['fspot']=fspot
            sdict['station']=z1.station
            linelst.append(sdict)
    
    #-----Write data file-------------------------------------------------------
    
    ofid=file(ofile,'w')
    ofid.write('{0:d} {1:d} {2:d}\n'.format(ns,nperiod,8))
    
    #write N-S locations
    ofid.write('Station_Location: N-S \n')
    for ii in range(ns/8+1):
        for ll in range(8):
            try:
                ofid.write('{0:+.4e} '.format(ylst[ii*8+ll]))
            except IndexError:
                pass
        ofid.write('\n')
    
    #write E-W locations
    ofid.write('Station_Location: E-W \n')
    for ii in range(ns/8+1):
        for ll in range(8):
            try:
                ofid.write('{0:+.4e} '.format(xlst[ii*8+ll]))
            except IndexError:
                pass
        ofid.write('\n')
        
    #write impedance tensor components
    for ii,p1 in enumerate(periodlst):
        ofid.write('DATA_Period: {0:3.6f}\n'.format(p1))
        for ss in range(ns):
            zline=zarr[ss,ii,:]
            for jj in range(4):
                ofid.write('{0:+.4e} '.format(zline[jj].real*zconv))
                ofid.write('{0:+.4e} '.format(-zline[jj].imag*zconv))
            ofid.write('\n')
    
    #write error as a percentage of Z
    for ii,p1 in enumerate(periodlst):
        ofid.write('ERROR_Period: {0:3.6f}\n'.format(p1))
        for ss in range(ns):
            zline=zarr[ss,ii,:]
            for jj in range(4):
                ofid.write('{0:+.4e} '.format(zline[jj].real*zerr*zconv))
                ofid.write('{0:+.4e} '.format(zline[jj].imag*zerr*zconv))
            ofid.write('\n')
            
    #write error maps
    for ii,p1 in enumerate(periodlst):
        ofid.write('ERMAP_Period: {0:3.6f}\n'.format(p1))
        for ss in range(ns):
            zline=zarr[ss,ii,:]
            for jj in range(4):
                ofid.write('{0:.5e} '.format(zerrmap[jj]))
                ofid.write('{0:.5e} '.format(zerrmap[jj]))
            ofid.write('\n')
    ofid.close()
    print 'Wrote file to: '+ofile
    
    #write out places where errors are larger than error tolerance
    errfid=file(os.path.join(os.path.dirname(ofile),'DataErrorLocations.txt'),
                'w')
    errfid.write('Errors larger than error tolerance of: \n')
    errfid.write('Zxx={0} Zxy={1} Zyx={2} Zyy={3} \n'.format(zerrmap[0]*zerr,
                 zerrmap[1]*zerr,zerrmap[2]*zerr,zerrmap[3]*zerr))
    errfid.write('-'*20+'\n')
    errfid.write('station  T=period(s) Zij err=percentage \n')
    for pfdict in linelst:
        for kk,ff in enumerate(pfdict['fspot']):
            if pfdict['fspot'][ff][2]>zerr*100*zerrmap[0]:
                errfid.write(pfdict['station']+'  T='+ff+\
                        ' Zxx err={0:.3f} \n'.format(pfdict['fspot'][ff][2])) 
            if pfdict['fspot'][ff][3]>zerr*100*zerrmap[1]:
                errfid.write(pfdict['station']+'  T='+ff+\
                        ' Zxy err={0:.3f} \n'.format(pfdict['fspot'][ff][3])) 
            if pfdict['fspot'][ff][4]>zerr*100*zerrmap[2]:
                errfid.write(pfdict['station']+'  T='+ff+\
                        ' Zyx err={0:.3f} \n'.format(pfdict['fspot'][ff][4]))
            if pfdict['fspot'][ff][5]>zerr*100*zerrmap[3]:
                errfid.write(pfdict['station']+'  T='+ff+\
                        ' Zyy err={0:.3f} \n'.format(pfdict['fspot'][ff][5])) 
    errfid.close()
    print 'Wrote errors lager than tolerance to: '
    print os.path.join(os.path.dirname(ofile),'DataErrorLocations.txt')
                
    
    return ofile,linelst


def writeInit3DFile(outfn,rhostart=100,ncol=5,savepath=None):
    """
    Makes an init3d file for WSINV3D
    
    Inputs:
        outfn = full path to .out file from winglink
        rhostart = starting homogeneous half space in Ohm-m
        ncol = number of columns for data to be written in
        savepath = full path to save the init file
        
    Output:
        ifile = full path to init file
    """
    
    #create the output filename
    if savepath==None:
        ifile=os.path.join(os.path.dirname(outfn),'init3d')
    elif savepath.find('.')==-1:
        ifile=os.path.join(savepath,'init3d')
    else:
        ifile=savepath
        
    dx,dy,dz=readWLOutFile(outfn,ncol=ncol)
    
    nx=len(dx)
    ny=len(dy)
    nz=len(dz)
    
    ifid=file(ifile,'w')
    ifid.write('#Initial model \n')
    ifid.write('{0} {1} {2} 1 \n'.format(ny,nx,nz))
        
    #write y locations
    for ii in range(ny/8+1):
        for jj in range(8):
            try:
                ifid.write('{0:+.4e} '.format(abs(dy[ii*8+jj])))
            except IndexError:
                pass
        ifid.write('\n')
        
    #write x locations
    for ii in range(nx/8+1):
        for jj in range(8):
            try:
                ifid.write('{0:+.4e} '.format(abs(dx[ii*8+jj])))
            except IndexError:
                pass
        ifid.write('\n')
        
    #write z locations
    for ii in range(nz/8+1):
        for jj in range(8):
            try:
                ifid.write('{0:+.4e} '.format(abs(dz[ii*8+jj])))
            except IndexError:
                pass
        ifid.write('\n')
        
    ifid.write('{0} \n'.format(rhostart))
    
    ifid.close()
    
    print 'Wrote init file to: '+ifile
    
    return ifile
        
def writeStartupFile(datafn,initialfn=None,outputfn=None,savepath=None,
                    apriorfn=None,modells=[5,0.3,0.3,0.3],targetrms=1.0,
                    control=None,maxiter=10,errortol=None,staticfn=None,
                    lagrange=None):
    """
    makes a startup file for WSINV3D t.  Most of these parameters are not input
    
    Inputs:
        datafn = full path to the data file written for inversion
        initialfn = full path to init file
        outputfn = output stem to which the _model and _resp will be written
        savepath = full path to save the startup file to
        apriorfn = full path to apriori model
        modells = smoothing parameters 
        targetrms = target rms
        control = something
        maxiter = maximum number of iterations
        errotol = error tolerance for the computer?
        staticfn = full path to static shift file name
        lagrange = starting lagrange multiplier
        
    Outputs:
        sfile = full path to startup file
        
    """
    
    #create the output filename
    if savepath==None:
        sfile=os.path.join(os.path.dirname(datafn),'startup')
    elif savepath.find('.')==-1:
        sfile=os.path.join(savepath,'startup')
    else:
        sfile=savepath
    
    sfid=file(sfile,'w')
    
    sfid.write('DATA_FILE'+' '*11+'../'+os.path.basename(datafn)+'\n')
 
    if outputfn==None:
        sfid.write('OUTPUT_FILE'+' '*9+'Iter_ \n')
    else:
        sfid.write('OUTPUT_FILE'+' '*9+outputfn+' \n')
        
    if initialfn==None:
        sfid.write('INITIAL_MODEL_FILE'+' '*2+'../init3d \n')
    else:
        sfid.write('INITIAL_MODEL_FILE'+' '*2+initialfn+' \n')
        
    if apriorfn==None:
        sfid.write('PRIOR_MODEL_FILE'+' '*4+'default \n')
    else:
        sfid.write('PRIOR_MODEL_FILE'+' '*4+apriorfn+' \n')
        
    if control==None:
        sfid.write('CONTROL_MODEL_INDEX'+' '+'default \n')
    else:
        sfid.write('CONTROL_MODEL_INDEX'+' '+control+' \n')
        
    sfid.write('TARGET_RMS'+' '*10+'{0} \n'.format(targetrms))
    
    sfid.write('MAX_NO_ITERATION'+' '*4+'{0} \n'.format(maxiter))

    sfid.write('MODEL_LENGTH_SCALE'+' '*2+
                '{0} {1:.1f} {1:.1f} {1:.1f} \n'.format(modells[0],modells[1],
                                                        modells[2],modells[3]))
        
    if initialfn==None:
        sfid.write('LAGRANGE_INFO'+' '*7+'default \n')
    else:
         sfid.write('LAGRANGE_INFO'+' '*7+lagrange+' \n')
    

    if initialfn==None:
        sfid.write('ERROR_TOL_LEVEL'+' '*5+'default \n')
    else:
         sfid.write('ERROR_TOL_LEVEL'+' '*5+errortol+' \n')
         
    if initialfn==None:
        sfid.write('STATIC_FILE'+' '*9+'default \n')
    else:
         sfid.write('STATIC_FILE'+' '*9+staticfn+' \n')

    sfid.close()
    
    print 'Wrote startup file to: '+sfile
    
    return sfile
    
def readDataFile(datafn,sitesfn=None,units='mv'):
    """
    read in data file
    
    Inputs:
        datafn = full path to data file
        sitesfn = full path to sites file output by winglink
        units = 'mv' always
        
    Outputs:
       period = list of periods used for the inversion
       zarr = array of impedance values 
               (number of staitons x number of periods x 2 x 2)
       zerr = array of errors in impedance component
       nsarr = station locations relative distance from center of grid in N-S
       ewarr = station locations relative distance from center of grid in E-W
       sitelst = list of sites used in data         
    """
    
    if units=='mv':
        zconv=796.
    else:
        zconv=1
    
        
    dfid=file(datafn,'r')
    dlines=dfid.readlines()

    #get size number of stations, number of frequencies, number of Z components    
    ns,nf,nz=np.array(dlines[0].strip().split(),dtype='int')
    nsstart=2
    
    findlst=[]
    for ii,dline in enumerate(dlines[1:50],1):
        if dline.find('Station_Location: N-S')==0:
            findlst.append(ii)
        elif dline.find('Station_Location: E-W')==0:
            findlst.append(ii)
        elif dline.find('DATA_Period:')==0:
            findlst.append(ii)
            
    ncol=len(dlines[nsstart].strip().split())
#    print ncol
#    nsstop=nsstart+ns/ncol+1
#    ewstart=nsstop+1
#    ewstop=ewstart+ns/ncol+1
#    zstart=ewstop
#    print nsstop,ewstart,ewstop,zstart
    
    #get site names if entered a sites file
    if sitesfn!=None:
        slst,sitelst=readSitesFile(sitesfn)
    else:
        sitelst=np.arange(ns)

    #get N-S locations
    nsarr=np.zeros(ns)
    for ii,dline in enumerate(dlines[findlst[0]+1:findlst[1]],0):
        dline=dline.strip().split()
        for jj in range(ncol):
            try:
                nsarr[ii*ncol+jj]=float(dline[jj])
            except IndexError:
                pass
            except ValueError:
                break
            
    #get E-W locations
    ewarr=np.zeros(ns)
    for ii,dline in enumerate(dlines[findlst[1]+1:findlst[2]],0):
        dline=dline.strip().split()
        for jj in range(8):
            try:
                ewarr[ii*ncol+jj]=float(dline[jj])
            except IndexError:
                pass
            except ValueError:
                break
    #make some empty array to put stuff into
    period=np.zeros(nf)
    zarr=np.zeros((ns,nf,2,2),dtype=np.complex)
    zerr=np.zeros_like(zarr)
    zerrmap=np.zeros_like(zarr)

    #get data
    pcount=0
    zcount=0
    for ii,dl in enumerate(dlines[findlst[2]:findlst[2]+nf*(ns+1)]):
        if dl.find('DATA_Period')==0:
            period[pcount]=float(dl.strip().split()[1])
            kk=0
            pcount+=1
            if ii==0:
                pass
            else:
                zcount+=1
        else:
            zline=np.array(dl.strip().split(),dtype=np.float)*zconv
            zarr[kk,zcount,:,:]=np.array([[zline[0]-1j*zline[1],
                                                zline[2]-1j*zline[3]],
                                                [zline[4]-1j*zline[5],
                                                 zline[6]-1j*zline[7]]])
            kk+=1
    
    #if the data file is made from this program or is the input data file than
    #get the errors from that file
    if len(dlines)>2*nf*ns:
        print 'Getting Error'
        pecount=0
        zecount=0
        for ii,dl in enumerate(dlines[findlst[2]+nf*(ns+1):findlst[2]+2*nf*(ns+1)]):
            if dl.find('ERROR_Period')==0:
                kk=0
                pecount+=1
                if ii==0:
                    pass
                else:
                    zecount+=1
            else:
                zline=np.array(dl.strip().split(),dtype=np.float)*zconv
                zerr[kk,zecount,:,:]=np.array([[zline[0]-1j*zline[1],
                                                    zline[2]-1j*zline[3]],
                                                    [zline[4]-1j*zline[5],
                                                     zline[6]-1j*zline[7]]])
                kk+=1
                
    #get errormap values
    if len(dlines)>3*nf*ns:
        print 'Getting Error Map'
        pmcount=0
        zmcount=0
        for ii,dl in enumerate(dlines[findlst[2]+2*nf*(ns+1):findlst[2]+3*nf*(ns+1)]):
            if dl.find('ERMAP_Period')==0:
                kk=0
                pmcount+=1
                if ii==0:
                    pass
                else:
                    zmcount+=1
            else:
                #account for end of file empty lines
                if len(dl.split())>2:
                    zline=np.array(dl.strip().split(),dtype=np.float)
                    zerrmap[kk,zmcount,:,:]=np.array([[zline[0]-1j*zline[1],
                                                        zline[2]-1j*zline[3]],
                                                        [zline[4]-1j*zline[5],
                                                         zline[6]-1j*zline[7]]])
                    kk+=1
    
    #multiply errmap and error and convert from Ohm to mv/km nT
    zerr=zerr*zerrmap                                           

        
    return period,zarr,zerr,nsarr,ewarr,sitelst
    
def plotDataResPhase(datafn,respfn=None,sitesfn=None,plottype='1',plotnum=1,
                     dpi=150,units='mv',colormode='color'):
    """
    plot responses from the data file and if there is a response file
    
    Inputs:
        datafn = fullpath to data file
        respfn = full path to respsonse file, if not input, just the data is
                 plotted. Can be a list of response files from the same 
                 inversion
        plottype= '1' to plot each station in a different window
                  [stations] for list of stations to plot (stations are numbers)
        plotnum = 1 for just xy,yx
                  2 for all components
    """
    #plot in color mode or black and white
    if colormode=='color':
        #color for data
        cted=(0,0,1)
        ctmd=(1,0,0)
        mted='*'
        mtmd='*'
        
        #color for occam model
        ctem=(0,.3,1.0)
        ctmm=(1,.3,0)
        mtem='+'
        mtmm='+'
        
    elif colormode=='bw':
        #color for data
        cted=(0,0,0)
        ctmd=(0,0,0)
        mted='*'
        mtmd='v'
        
        #color for occam model
        ctem=(0.6,.6,.6)
        ctmm=(.6,.6,.6)
        mtem='+'
        mtmm='x'
    
    
    #load the data file     
    period,dz,dzerr,north,east,slst=readDataFile(datafn,sitesfn=sitesfn,
                                                 units=units)
    #get shape of impedance tensors
    ns,nf=dz.shape[0],dz.shape[1]

    #read in response files
    if respfn!=None:
        rzlst=[]
        rzerrlst=[]
        if type(respfn) is not list:
            respfn=[respfn]
        for rfile in respfn:
            period,rz,rzerr,north,east,slst=readDataFile(rfile,sitesfn=sitesfn,
                                                         units=units)
            rzlst.append(rz)
            rzerrlst.append(rzerr)
    else:
        rzlst=[]
    #get number of response files
    nr=len(rzlst)
    
    if type(plottype) is list:
        ns=len(plottype)
      
    plt.rcParams['font.size']=10
    plt.rcParams['figure.subplot.left']=.13
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.1
    plt.rcParams['figure.subplot.top']=.92
    plt.rcParams['figure.subplot.wspace']=.25
    plt.rcParams['figure.subplot.hspace']=.05
    
    fontdict={'size':12,'weight':'bold'}    
    gs=gridspec.GridSpec(2,2,height_ratios=[2,1.5],hspace=.1)    
    
    
    if plottype!='1':
        pstationlst=[]
        if type(plottype) is not list:
            plottype=[plottype]
        for ii,station in enumerate(slst):
            if type(station) is str:
                for pstation in plottype:
                    if station.find(str(pstation))>=0:
                        pstationlst.append(ii)
            else:
                for pstation in plottype:
                    if station==int(pstation):
                        pstationlst.append(ii)
    else:
        pstationlst=np.arange(ns)
    
    for jj in pstationlst:
        print 'Plotting: '+str(slst[jj])
        
        #check for masked points
        dz[jj][np.where(dz[jj]==7.95204E5-7.95204E5j)]=0.0+0.0j
        dzerr[jj][np.where(dz[jj]==7.95204E5-7.95204E5j)]=1.0+1.0j
        
        #convert to apparent resistivity and phase
        rp=Z.ResPhase(dz[jj],period,zvar=dzerr[jj])
        
        #find locations where points have been masked
        nzxx=np.where(rp.resxx!=0)[0]
        nzxy=np.where(rp.resxy!=0)[0]
        nzyx=np.where(rp.resyx!=0)[0]
        nzyy=np.where(rp.resyy!=0)[0]
        
        if respfn!=None:
            plotr=True
        else:
            plotr=False
        
        #make figure for xy,yx components
        if plotnum==1: 
            fig=plt.figure(jj,[10,12],dpi=dpi)
            gs.update(hspace=.1,wspace=.15,left=.1)
        elif plotnum==2:
            fig=plt.figure(jj,[12,12],dpi=dpi)
            gs.update(hspace=.1,wspace=.15,left=.07)
        
        #---------plot the apparent resistivity-----------------------------------
        if plotnum==1:
            ax=fig.add_subplot(gs[0,:])
            ax2=fig.add_subplot(gs[1,:],sharex=ax)
            ax.yaxis.set_label_coords(-.055, 0.5)
            ax2.yaxis.set_label_coords(-.055, 0.5)
        elif plotnum==2:
            ax=fig.add_subplot(gs[0,0])
            ax2=fig.add_subplot(gs[1,0],sharex=ax)
            ax.yaxis.set_label_coords(-.075, 0.5)
            ax2.yaxis.set_label_coords(-.075, 0.5)
        
        fig.suptitle(str(slst[jj]),fontdict={'size':15,'weight':'bold'})
        erxy=ax.errorbar(period[nzxy],rp.resxy[nzxy],marker=mted,ms=4,
                         mfc='None',mec=cted,mew=1,ls=':',
                         yerr=rp.resxyerr[nzxy],ecolor=cted,color=cted)
        eryx=ax.errorbar(period[nzyx],rp.resyx[nzyx],marker=mtmd,ms=4,
                         mfc='None',mec=ctmd,mew=1,ls=':',
                         yerr=rp.resyxerr[nzyx],ecolor=ctmd,color=ctmd)
        if plotr==True:
            for rr in range(nr):
                if colormode=='color':   
                    cxy=(0,.4+float(rr)/(3*nr),0)
                    cyx=(.7+float(rr)/(4*nr),.13,.63-float(rr)/(4*nr))
                elif colormode=='bw':
                    cxy=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))                    
                    cyx=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))
                
                rpr=Z.ResPhase(rzlst[rr][jj],period,zvar=rzerrlst[rr][jj])
                
#                rms=np.sqrt(np.sum([abs(np.linalg.det(rp.z[ll])-
#                                        np.linalg.det(rpr.z[ll]))**2 
#                            for ll in range(len(rp.period))])/len(rp.period))
                rms=np.sqrt(np.mean([abs(abs(np.linalg.det(rp.z[ll]))-
                                    abs(np.linalg.det(rpr.z[ll]))) 
                                    for ll in range(len(rp.period))]))
                print rms
                erxyr=ax.errorbar(period[nzxy],rpr.resxy[nzxy],marker=mtem,
                                  ms=8,mfc='None',mec=cxy,mew=1,ls='--',
                                  yerr=rpr.resxyerr[nzxy],
                                  ecolor=cxy,color=cxy)
                eryxr=ax.errorbar(period[nzyx],rpr.resyx[nzyx],marker=mtmm,
                                  ms=8,mfc='None',mec=cyx,mew=1,ls='--',
                                  yerr=rpr.resyxerr[nzyx],
                                  ecolor=cyx,color=cyx)
        #ax.set_xlabel('Period (s)',fontdict=fontdict)
        pylab.setp( ax.get_xticklabels(), visible=False)
        ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                   fontdict=fontdict)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
                 xmax=10**(np.ceil(np.log10(period[-1]))))
        ax.grid(True,alpha=.25)
        if plotr==True:
            ax.legend((erxy[0],eryx[0],erxyr[0],eryxr[0]),
                      ('Data $E_x/B_y$','Data $E_y/B_x$',
                      'Mod $E_x/B_y$','Mod $E_y/B_x$'),
                      loc=0, markerscale=1,borderaxespad=.01,labelspacing=.07,
                      handletextpad=.2,borderpad=.02)
        else:
            ax.legend((erxy[0],eryx[0]),('$E_x/B_y$','$E_y/B_x$'),loc=0,
                        markerscale=1,borderaxespad=.01,labelspacing=.07,
                        handletextpad=.2,borderpad=.02)
        
        #-----Plot the phase----------------------------------------------------
        
        ax2.errorbar(period[nzxy],rp.phasexy[nzxy],marker=mted,ms=4,mfc='None',
                     mec=cted,mew=1,ls=':',yerr=rp.phasexyerr[nzxy],ecolor=cted,
                     color=cted)
        ax2.errorbar(period[nzyx],np.array(rp.phaseyx[nzyx])+180,marker=mtmd,
                     ms=4,mfc='None',mec=ctmd,mew=1,ls=':',
                     yerr=rp.phaseyxerr[nzyx],
                     ecolor=ctmd,color=ctmd)
        if plotr==True:
            for rr in range(nr):
                if colormode=='color':   
                    cxy=(0,.4+float(rr)/(3*nr),0)
                    cyx=(.7+float(rr)/(4*nr),.13,.63-float(rr)/(4*nr))
                elif colormode=='bw':
                    cxy=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))                    
                    cyx=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))
                rpr=Z.ResPhase(rzlst[rr][jj],period,zvar=rzerrlst[rr][jj])
                ax2.errorbar(period[nzxy],rpr.phasexy[nzxy],marker=mtem,ms=8,
                             mfc='None',mec=cxy,mew=1,ls='--',
                             yerr=rp.phasexyerr[nzxy],
                             ecolor=cxy,color=cxy)
                ax2.errorbar(period[nzyx],np.array(rpr.phaseyx[nzyx])+180,
                             marker=mtmm,ms=8,mfc='None',mec=cyx,mew=1,ls='--',
                             yerr=rp.phaseyxerr[nzyx],ecolor=cyx,color=cyx)
        ax2.set_xlabel('Period (s)',fontdict)
        ax2.set_ylabel('Phase (deg)',fontdict)
        ax2.set_xscale('log')
        #ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
        #         xmax=10**(np.ceil(np.log10(period[-1]))))
        #check the phase to see if any point are outside of [0:90]    
        if min(rp.phasexy)<0 or min(rp.phaseyx+180)<0:
            pymin=min([min(rp.phasexy),min(rp.phaseyx+180)])
            if pymin>0:
                pymin=0
        else:
            pymin=0
        
        if max(rp.phasexy)>90 or max(rp.phaseyx+180)>90:
            pymax=min([max(rp.phasexy),max(rp.phaseyx+180)])
            if pymax<91:
                pymax=90
        else:
            pymax=90
        
        ax2.set_ylim(ymin=pymin,ymax=pymax)        
        ax2.yaxis.set_major_locator(MultipleLocator(30))
        ax2.yaxis.set_minor_locator(MultipleLocator(1))
        ax2.grid(True,alpha=.25)
        
        if plotnum==2:
            #---------plot the apparent resistivity-----------------------------------
            ax3=plt.subplot(gs[0,1])
            ax3.yaxis.set_label_coords(-.1, 0.5)
            erxx=ax3.errorbar(period[nzxx],rp.resxx[nzxx],marker=mted,ms=4,
                              mfc='None',mec=cted,mew=1,ls=':',
                              yerr=rp.resxxerr[nzxx],
                              ecolor=cted,color=cted)
            eryy=ax3.errorbar(period[nzyy],rp.resyy[nzyy],marker=mtmd,ms=4,
                              mfc='None',mec=ctmd,mew=1,ls=':',
                              yerr=rp.resyyerr[nzyy],
                              ecolor=ctmd,color=ctmd)
            if plotr==True:
                for rr in range(nr):
                    if colormode=='color':   
                        cxy=(0,.4+float(rr)/(3*nr),0)
                        cyx=(.7+float(rr)/(4*nr),.13,.63-float(rr)/(4*nr))
                    elif colormode=='bw':
                        cxy=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))                    
                        cyx=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))
                    rpr=Z.ResPhase(rzlst[rr][jj],period,zvar=rzerrlst[rr][jj])
                    erxxr=ax3.errorbar(period[nzxx],rpr.resxx[nzxx],
                                       marker=mtem,ms=8,mfc='None',mec=cxy,
                                       mew=1,ls='--',yerr=rpr.resxxerr[nzxx],
                                       ecolor=cxy,color=cxy)
                    eryyr=ax3.errorbar(period[nzyy],rpr.resyy[nzyy],
                                       marker=mtmm,ms=8,mfc='None',mec=cyx,
                                       mew=1,ls='--',yerr=rpr.resyyerr[nzyy],
                                       ecolor=cyx,color=cyx)

            ax3.set_yscale('log')
            ax3.set_xscale('log')
            pylab.setp( ax3.get_xticklabels(), visible=False)
            ax3.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
                     xmax=10**(np.ceil(np.log10(period[-1]))))
            ax3.grid(True,alpha=.25)
            if plotr==True:
                ax3.legend((erxx[0],eryy[0],erxxr[0],eryyr[0]),
                          ('Data $E_x/B_x$','Data $E_y/B_y$',
                          'Mod $E_x/B_x$','Mod $E_y/B_y$'),
                          loc=0, markerscale=1,borderaxespad=.01,
                          labelspacing=.07,handletextpad=.2,borderpad=.02)
            else:
                ax3.legend((erxx[0],eryy[0]),('$E_x/B_x$','$E_y/B_y$'),loc=0,
                            markerscale=1,borderaxespad=.01,labelspacing=.07,
                            handletextpad=.2,borderpad=.02)
            
            #-----Plot the phase----------------------------------------------------
            ax4=plt.subplot(gs[1,1],sharex=ax3)
            
            ax4.yaxis.set_label_coords(-.1, 0.5)
            ax4.errorbar(period[nzxx],rp.phasexx[nzxx],marker=mted,ms=4,
                         mfc='None',mec=cted,mew=1,ls=':',
                         yerr=rp.phasexxerr[nzxx],ecolor=cted,color=cted)
            ax4.errorbar(period[nzyy],np.array(rp.phaseyy[nzyy]),marker=mtmd,
                         ms=4,mfc='None',mec=ctmd,mew=1,ls=':',
                         yerr=rp.phaseyyerr[nzyy],
                         ecolor=ctmd,color=ctmd)
            if plotr==True:
                for rr in range(nr):
                    if colormode=='color':   
                        cxy=(0,.4+float(rr)/(3*nr),0)
                        cyx=(.7+float(rr)/(4*nr),.13,.63-float(rr)/(4*nr))
                    elif colormode=='bw':
                        cxy=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))                    
                        cyx=(1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))
                    rpr=Z.ResPhase(rzlst[rr][jj],period,zvar=rzerrlst[rr][jj])
                    ax4.errorbar(period[nzxx],rpr.phasexx[nzxx],marker=mtem,
                                 ms=8,mfc='None',mec=cxy,mew=1,ls='--',
                                 yerr=rp.phasexxerr[nzxx],
                                 ecolor=cxy,color=cxy)
                    ax4.errorbar(period[nzyy],np.array(rpr.phaseyy[nzyy]),
                                 marker=mtmm,ms=8,mfc='None',mec=cyx,mew=1,
                                 ls='--',yerr=rp.phaseyyerr[nzyy], 
                                 ecolor=cyx,color=cyx)
            ax4.set_xlabel('Period (s)',fontdict)
            #ax4.set_ylabel('Imepdance Phase (deg)',fontdict)
            ax4.set_xscale('log')
            #ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
            #         xmax=10**(np.ceil(np.log10(period[-1]))))
            ax4.set_ylim(ymin=-180,ymax=180)        
            ax4.yaxis.set_major_locator(MultipleLocator(30))
            ax4.yaxis.set_minor_locator(MultipleLocator(5))
            ax4.grid(True,alpha=.25)

def plotTensorMaps(datafn,respfn=None,sitesfn=None,periodlst=None,
                   esize=(1,1,5,5),ecolor='phimin',
                   colormm=[(0,90),(0,1),(0,4),(-2,2)],
                   xpad=.500,units='mv',dpi=150):
    """
    plot phase tensor maps for data and or response, each figure is of a
    different period.  If response is input a third column is added which is 
    the residual phase tensor showing where the model is not fitting the data 
    well.  The data is plotted in km in units of ohm-m.
    
    Inputs:
        datafn = full path to data file
        respfn = full path to response file, if none just plots data
        sitesfn = full path to sites file
        periodlst = indicies of periods you want to plot
        esize = size of ellipses as:
                0 = phase tensor ellipse
                1 = phase tensor residual
                2 = resistivity tensor ellipse
                3 = resistivity tensor residual
        ecolor = 'phimin' for coloring with phimin or 'beta' for beta coloring
        colormm = list of min and max coloring for plot, list as follows:
                0 = phase tensor min and max for ecolor in degrees
                1 = phase tensor residual min and max [0,1]
                2 = resistivity tensor coloring as resistivity on log scale
                3 = resistivity tensor residual coloring as resistivity on 
                    linear scale
        xpad = padding of map from stations at extremities (km)
        units = 'mv' to convert to Ohm-m 
        dpi = dots per inch of figure
    """
    
    period,zd,zderr,nsarr,ewarr,sitelst=readDataFile(datafn,sitesfn=sitesfn,
                                                      units=units)
    
    if respfn!=None:
        period,zr,zrerr,nsarr,ewarr,sitelst=readDataFile(respfn,sitesfn=sitesfn,
                                                         units=units)
    
    if periodlst==None:
        periodlst=range(len(period))
        
    #put locations into an logical coordinate system in km
    nsarr=-nsarr/1000
    ewarr=-ewarr/1000

    #get coloring min's and max's    
    if colormm!=None:
        ptmin,ptmax=(colormm[0][0]*np.pi/180,colormm[0][1]*np.pi/180)
        ptrmin,ptrmax=colormm[1]
        rtmin,rtmax=colormm[2]
        rtrmin,rtrmax=colormm[3]
    else:
        pass
    
    #get ellipse sizes
    ptsize=esize[0]
    ptrsize=esize[1]
    rtsize=esize[2]
    rtrsize=esize[3]
        
    plt.rcParams['font.size']=10
    plt.rcParams['figure.subplot.left']=.03
    plt.rcParams['figure.subplot.right']=.98
    plt.rcParams['figure.subplot.bottom']=.1
    plt.rcParams['figure.subplot.top']=.90
    plt.rcParams['figure.subplot.wspace']=.005
    plt.rcParams['figure.subplot.hspace']=.005
    
    ns=zd.shape[0]
    
    for ff,per in enumerate(periodlst):
        print 'Plotting Period: {0:.5g}'.format(period[per])
        fig=plt.figure(per+1,dpi=dpi)

        #get phase tensor
        pt=Z.PhaseTensor(zd[:,per])

        #get resistivity tensor
        rt=Z.ResistivityTensor(zd[:,per],np.repeat(1./period[per],ns))
        
        if respfn!=None:
            #get phase tensor and residual phase tensor
            ptr=Z.PhaseTensor(zr[:,per])
            ptd=Z.PhaseTensorResidual(zd[:,per],zr[:,per])
            
            #get resistivity tensor and residual
            rtr=Z.ResistivityTensor(zr[:,per],np.repeat(1./period[per],ns))
            rtd=Z.ResistivityTensorResidual(zd[:,per],zr[:,per],
                                            np.repeat(1./period[per],ns))
            
            if colormm==None:
                if ecolor=='phimin':
                    ptmin,ptmax=(ptr.phimin.min()/(np.pi/2),
                                 ptr.phimin.max()/(np.pi/2))
                elif ecolor=='beta':
                    ptmin,ptmax=(ptr.beta.min(),ptr.beta.max())
                    
                ptrmin,ptrmax=(ptd.ecolor.min(),ptd.ecolor.max())
                rtmin,rtmax=(np.log10(rtr.rhodet.min()),
                             np.log10(rtr.rhodet.max()))
                rtrmin,rtrmax=rtd.rhodet.min(),rtd.rhodet.max()
            #make subplots            
            ax1=fig.add_subplot(2,3,1,aspect='equal')
            ax2=fig.add_subplot(2,3,2,aspect='equal')
            ax3=fig.add_subplot(2,3,3,aspect='equal')
            ax4=fig.add_subplot(2,3,4,aspect='equal')
            ax5=fig.add_subplot(2,3,5,aspect='equal')
            ax6=fig.add_subplot(2,3,6,aspect='equal')
            
            for jj in range(ns):
                #-----------plot data phase tensors---------------
                eheightd=pt.phimin[jj]/ptr.phimax.max()*ptsize
                ewidthd=pt.phimax[jj]/ptr.phimax.max()*ptsize
            
                ellipd=Ellipse((ewarr[jj],nsarr[jj]),width=ewidthd,
                               height=eheightd,angle=pt.azimuth[jj])
                #color ellipse:
                if ecolor=='phimin':
                    cvar=(pt.phimin[jj]/(np.pi/2)-ptmin)/(ptmax-ptmin)
                    if abs(cvar)>1:
                        ellipd.set_facecolor((1,0,.1))
                    else:
                        ellipd.set_facecolor((1,1-abs(cvar),.1))
                if ecolor=='beta':
                    cvar=(abs(pt.beta[jj])-ptmin)/(ptmax-ptmin)
                    if abs(cvar)>1:
                        ellipd.set_facecolor((1,1,.1))
                    else:
                        ellipd.set_facecolor((1-cvars,1-cvars,1))
                
                ax1.add_artist(ellipd)
                
                #----------plot response phase tensors---------------------
                eheightr=ptr.phimin[jj]/ptr.phimax.max()*ptsize
                ewidthr=ptr.phimax[jj]/ptr.phimax.max()*ptsize
            
                ellipr=Ellipse((ewarr[jj],nsarr[jj]),width=ewidthr,
                               height=eheightr,angle=ptr.azimuth[jj])
                #color ellipse:
                if ecolor=='phimin':
                    cvar=(ptr.phimin[jj]/(np.pi/2)-ptmin)/(ptmax-ptmin)
                    if abs(cvar)>1:
                        ellipr.set_facecolor((1,0,.1))
                    else:
                        ellipr.set_facecolor((1,1-abs(cvar),.1))
                if ecolor=='beta':
                    cvar=(abs(ptr.beta[jj])-ptmin)/(ptmax-ptmin)
                    if abs(cvar)>1:
                        ellipr.set_facecolor((1,1,.1))
                    else:
                        ellipr.set_facecolor((1-cvars,1-cvars,1))
                ax2.add_artist(ellipr)
                
                #--------plot residual phase tensors-------------
                eheight=ptd.phimin[jj]/ptd.phimax.max()*ptrsize
                ewidth=ptd.phimax[jj]/ptd.phimax.max()*ptrsize
            
                ellip=Ellipse((ewarr[jj],nsarr[jj]),width=ewidth,
                               height=eheight,angle=ptd.azimuth[jj]-90)
                #color ellipse:
                cvar=(ptd.ecolor[jj]-ptrmin)/(ptrmax-ptrmin)
                if abs(cvar)>1:
                    ellip.set_facecolor((0,0,0))
                else:
                    ellip.set_facecolor((abs(cvar),.5,.5))
                
                ax3.add_artist(ellip)
                
                #-----------plot data resistivity tensors---------------
                rheightd=rt.rhomin[jj]/rtr.rhomax.max()*rtsize
                rwidthd=rt.rhomax[jj]/rtr.rhomax.max()*rtsize
            
                rellipd=Ellipse((ewarr[jj],nsarr[jj]),width=rwidthd,
                               height=rheightd,angle=rt.rhoazimuth[jj])
                #color ellipse:
                cvar=(np.log10(rt.rhodet[jj])-rtmin)/(rtmax-rtmin)
                if cvar>.5:
                    if cvar>1:
                        rellipd.set_facecolor((0,0,1))
                    else:
                        rellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                else:
                    if cvar<-1:
                        rellipd.set_facecolor((1,0,0))
                    else:
                        rellipd.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
                
                ax4.add_artist(rellipd)
                
                #----------plot response resistivity tensors---------------------
                rheightr=rtr.rhomin[jj]/rtr.rhomax.max()*rtsize
                rwidthr=rtr.rhomax[jj]/rtr.rhomax.max()*rtsize
            
                rellipr=Ellipse((ewarr[jj],nsarr[jj]),width=rwidthr,
                               height=rheightr,angle=rtr.rhoazimuth[jj])

                #color ellipse:
                cvar=(np.log10(rtr.rhodet[jj])-rtmin)/(rtmax-rtmin)
                if cvar>.5:
                    if cvar>1:
                        rellipr.set_facecolor((0,0,1))
                    else:
                        rellipr.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                else:
                    if cvar<-1:
                        rellipr.set_facecolor((1,0,0))
                    else:
                        rellipr.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
                
                ax5.add_artist(rellipr)
                
                #--------plot residual resistivity tensors-------------
                rheight=rtd.rhomin[jj]/rtd.rhomax.max()*rtrsize
                rwidth=rtd.rhomax[jj]/rtd.rhomax.max()*rtrsize
            
                rellip=Ellipse((ewarr[jj],nsarr[jj]),width=rwidth,
                               height=rheight,angle=rtd.azimuth[jj]-90)
                #color ellipse:
                cvar=(rtd.rhodet[jj]-rtrmin)/(rtrmax-rtrmin)
                if cvar<0:
                    if cvar<-1:
                        rellip.set_facecolor((0,0,1))
                    else:
                        rellip.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                else:
                    if cvar>1:
                        rellip.set_facecolor((1,0,0))
                    else:
                        rellip.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
                    
                ax6.add_artist(rellip)
                
            for aa,ax in enumerate([ax1,ax2,ax3,ax4,ax5,ax6]):
                ax.set_xlim(ewarr.min()-xpad,ewarr.max()+xpad)
                ax.set_ylim(nsarr.min()-xpad,nsarr.max()+xpad)
                ax.grid('on')
                if aa<3:
                    pylab.setp(ax.get_xticklabels(),visible=False)
                if aa==0 or aa==3:
                    pass
                else:
                    pylab.setp(ax.get_yticklabels(),visible=False)
                
                cbax=make_axes(ax,shrink=.9,pad=.05,orientation='vertical')
                if aa==0 or aa==1:
                    cbx=ColorbarBase(cbax[0],cmap=ptcmap,
                                     norm=Normalize(vmin=ptmin*180/np.pi,
                                                    vmax=ptmax*180/np.pi),
                                     orientation='vertical',format='%.2g')
                    
                    cbx.set_label('Phase (deg)',
                                  fontdict={'size':7,'weight':'bold'})
                if aa==2:
                    cbx=ColorbarBase(cbax[0],cmap=ptcmap2,
                                     norm=Normalize(vmin=ptrmin,
                                                    vmax=ptrmax),
                                     orientation='vertical',format='%.2g')
                    
                    cbx.set_label('$\Delta_{\Phi}$',
                                  fontdict={'size':7,'weight':'bold'})
                if aa==3 or aa==4:
                    cbx=ColorbarBase(cbax[0],cmap=rtcmapr,
                                     norm=Normalize(vmin=10**rtmin,
                                                    vmax=10**rtmax),
                                     orientation='vertical',format='%.2g')
                    
                    cbx.set_label('App. Res. ($\Omega \cdot$m)',
                                  fontdict={'size':7,'weight':'bold'})
                if aa==5:
                    cbx=ColorbarBase(cbax[0],cmap=rtcmap,
                                     norm=Normalize(vmin=rtrmin,
                                                    vmax=rtrmax),
                                     orientation='vertical',format='%.2g')
                    
                    cbx.set_label('$\Delta_{rho}$',
                                  fontdict={'size':7,'weight':'bold'})
                 
            plt.show()
        
        #----Plot Just the data------------------        
        else:
            if colormm==None:
                if ecolor=='phimin':
                    ptmin,ptmax=(pt.phimin.min()/(np.pi/2),
                                 pt.phimin.max()/(np.pi/2))
                elif ecolor=='beta':
                    ptmin,ptmax=(pt.beta.min(),pt.beta.max())
                    
                rtmin,rtmax=(np.log10(rt.rhodet.min()),
                             np.log10(rt.rhodet.max()))
            ax1=fig.add_subplot(1,2,1,aspect='equal')
            ax2=fig.add_subplot(1,2,2,aspect='equal')
            for jj in range(ns):
                #-----------plot data phase tensors---------------
                #check for nan in the array cause it messes with the max                
                pt.phimax=np.nan_to_num(pt.phimax)
                
                #scale the ellipse
                eheightd=pt.phimin[jj]/pt.phimax.max()*ptsize
                ewidthd=pt.phimax[jj]/pt.phimax.max()*ptsize
            
                #make the ellipse
                ellipd=Ellipse((ewarr[jj],nsarr[jj]),width=ewidthd,
                               height=eheightd,angle=pt.azimuth[jj])
                #color ellipse:
                if ecolor=='phimin':
                    cvar=(pt.phimin[jj]/(np.pi/2)-ptmin)/(ptmax-ptmin)
                    if abs(cvar)>1:
                        ellipd.set_facecolor((1,0,.1))
                    else:
                        ellipd.set_facecolor((1,1-abs(cvar),.1))
                if ecolor=='beta':
                    cvar=(abs(pt.beta[jj])-ptmin)/(ptmax-ptmin)
                    if abs(cvar)>1:
                        ellipd.set_facecolor((1,1,.1))
                    else:
                        ellipd.set_facecolor((1-cvars,1-cvars,1))
                
                ax1.add_artist(ellipd)
                
                #-----------plot data resistivity tensors---------------
                rt.rhomax=np.nan_to_num(rt.rhomax)
                rheightd=rt.rhomin[jj]/rt.rhomax.max()*rtsize
                rwidthd=rt.rhomax[jj]/rt.rhomax.max()*rtsize
            
                rellipd=Ellipse((ewarr[jj],nsarr[jj]),width=rwidthd,
                               height=rheightd,angle=rt.rhoazimuth[jj])
                #color ellipse:
                cvar=(np.log10(rt.rhodet[jj])-rtmin)/(rtmax-rtmin)
                if cvar>.5:
                    if cvar>1:
                        rellipd.set_facecolor((0,0,1))
                    else:
                        rellipd.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                else:
                    if cvar<-1:
                        rellipd.set_facecolor((1,0,0))
                    else:
                        rellipd.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
                
                ax2.add_artist(rellipd)
                
            for aa,ax in enumerate([ax1,ax2]):
                ax.set_xlim(ewarr.min()-xpad,ewarr.max()+xpad)
                ax.set_ylim(nsarr.min()-xpad,nsarr.max()+xpad)
                ax.grid('on')
                ax.set_xlabel('easting (km)',fontdict={'size':10,
                              'weight':'bold'})

                if aa==1:
                    pylab.setp(ax.get_yticklabels(),visible=False)
                else:
                    ax.set_ylabel('northing (km)',fontdict={'size':10,
                              'weight':'bold'})
#                cbax=make_axes(ax,shrink=.8,pad=.15,orientation='horizontal',
#                               anchor=(.5,1))
                #l,b,w,h
#                cbax=fig.add_axes([.1,.95,.35,.05])
                if aa==0:
                    cbax=fig.add_axes([.12,.97,.31,.02])
                    cbx=ColorbarBase(cbax,cmap=ptcmap,
                                     norm=Normalize(vmin=ptmin*180/np.pi,
                                                    vmax=ptmax*180/np.pi),
                                     orientation='horizontal',format='%.2g')
                    
                    cbx.set_label('Phase (deg)',
                                  fontdict={'size':7,'weight':'bold'})
                if aa==1:
                    cbax=fig.add_axes([.59,.97,.31,.02])
                    cbx=ColorbarBase(cbax,cmap=rtcmapr,
                                     norm=Normalize(vmin=10**rtmin,
                                                    vmax=10**rtmax),
                                     orientation='horizontal',format='%.2g')
                    
                    cbx.set_label('App. Res. ($\Omega \cdot$m)',
                                  fontdict={'size':7,'weight':'bold'})
                    cbx.set_ticks((10**rtmin,10**rtmax))
            plt.show()
            

def readModelFile(mfile,ncol=7):
    """
    read in a model file
    """            
    
    mfid=file(mfile,'r')
    mlines=mfid.readlines()

    #get info at the beggining of file
    info=mlines[0].strip().split()
    infodict=dict([(info[0][1:],info[1]),(info[2],info[3]),(info[4],info[5])])
    
    #get lengths of things
    nx,ny,nz,nn=np.array(mlines[1].strip().split(),dtype=np.int)
    
    #make empty arrays to put stuff into
    xarr=np.zeros(nx)
    yarr=np.zeros(ny)
    zarr=np.zeros(nz)
#    resarr=np.zeros((nx,ny,nz))
    resarr=np.zeros(nx*ny*nz)
    
    #make indices to get information from
    ls=2
    lx=nx/ncol+ls
    ly=lx+ny/ncol
    lz=ly+nz/ncol+1
    
    print lx,ly,lz
    
    #get x positions
    for ii,mline in enumerate(mlines[ls:lx+1]):
        xline=mline.strip().split()
        for ll in range(ncol):
            try:
                xarr[ll+ii*ncol]=float(xline[ll])/1000.
            except IndexError:
                break
    xarr[nx/2:]=-xarr[nx/2:]

    #get y positions
    for ii,mline in enumerate(mlines[lx:ly+1]):
        yline=mline.strip().split()
        for ll in range(ncol):
            try:
                yarr[ll+ii*ncol]=float(yline[ll])/1000.
            except IndexError:
                break
    yarr[ny/2:]=-yarr[ny/2:]
            
    #get z positions        
    for ii,mline in enumerate(mlines[ly+1:lz+1]):
        zline=mline.strip().split()
        for ll in range(ncol):
            try:
                zarr[ll+ii*ncol]=-float(zline[ll])/1000.
            except IndexError:
                break
        
#    get resistivity values
    mm=0
    for kk in range(nz):
        for jj in range(ny):
            for ii in range(nx):
                resarr[mm]=float(mlines[lz+1+mm].strip())
                mm+=1
    
#    resarr=resarr.reshape(nz,ny,nx)
    return xarr,yarr,zarr,resarr,infodict
        
                
    
        
    
                        
        
        
        
    
    
