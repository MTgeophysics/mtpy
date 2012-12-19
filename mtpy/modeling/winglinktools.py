# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:19:30 2011

@author: a1185872
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator

def readOutputFile(outputfile):
    """readOutputFile will read an output file from winglink and output data
    in the form of a dictionary.

    Input:
        outputfile = the full path and filename of outputfile

    Output:
        idict = dictionary with keys of station name
                each idict[station name] is a dictionary with keys
                corresponding to modeled and observed responses:
                    'obsresxy','obsphasexy','modresxy','modphasexy','obsresyx',
                    'obsphaseyx','modresyx','modphaseyx','obshzres',
                    'obshzphase','modhzres','modhzphase','period'
        rplst = list of dictionaries for each station with keywords:
            'station' = station name
            'offset' = relative offset,
            'resxy' = TE resistivity and error as row 0 and 1 respectively,
            'resyx'= TM resistivity and error as row 0 and 1 respectively,
            'phasexy'= TE phase and error as row 0 and 1 respectively,
            'phaseyx'= Tm phase and error as row 0 and 1 respectively,
            'realtip'= Real Tipper and error as row 0 and 1 respectively,
            'imagtip'= Imaginary Tipper and error as row 0 and 1 respectively
        plst = periodlst as the median of all stations.
        stationlst = list of stations in order from profile
        title = list of parameters for plotting as [title,profile,inversiontype]

    """

    ofid=open(outputfile,'r')
    lines=ofid.readlines()

    idict={}
    stationlst=[]

    #get title line
    titleline=lines[1].replace('"','')
    titleline=titleline.rstrip().split(',')
    title=titleline[1].split(':')[1]
    profile=titleline[0].split(':')[1]
    inversiontype=lines[2].rstrip()

    dkeys=['obsresyx','obsphaseyx','modresyx','modphaseyx','obsresxy',
           'obsphasexy','modresxy','modphasexy','obshzres','obshzphase',
           'modhzres','modhzphase','period']

    for line in lines[3:]:
        if line.find('Data for station')==0:
            station=line.rstrip().split(':')[1][1:]
            idict[station]={}
            stationlst.append(station)
            print 'Read in station: ',station
            for key in dkeys:
                idict[station][key]=[]
        elif line.find('RMS')==0:
            idict[station]['rms']=float(line.strip().split(' = ')[1])
        elif line.find('==')==0:
            pass
        else:
            linelst=line.split()
            if len(linelst)==len(dkeys):
                for kk,key in enumerate(dkeys):
                    try:
                        if key.find('phase')>=0:
                            idict[station][key].append(-1*float(linelst[kk]))
                        else:
                            idict[station][key].append(float(linelst[kk]))
                    except ValueError:
                        idict[station][key].append(0)
            else:
                pass

    #get data into a more useful format that takes into account any masking of
    #data points.

    #get the median of period lists for survey
    plst=np.median(np.array([idict[station]['period'] for station in stationlst]),
                   axis=0)
    #length of period
    nperiod=len(plst)

    #make a dictionary of period indicies
    pdict=dict([('%2.4g' % key,ii) for ii,key in enumerate(plst)])

    #make a dictionary of indicies for spots to put res_ij and phase_ij
    wldict={}
    for dkey in dkeys:
        if dkey[0:3].find('obs')==0:
            wldict[dkey]=(dkey[3:],0)
        elif dkey[0:3].find('mod')==0:
            wldict[dkey]=(dkey[3:],1)

    #make empty arrays to put things into
    asize=(2,nperiod)
    rplst=[{'station':station,
            'resxy':np.zeros(asize),
            'resyx':np.zeros(asize),
            'phasexy':np.zeros(asize),
            'phaseyx':np.zeros(asize),
            'hzres':np.zeros(asize),
            'hzphase':np.zeros(asize),
            } for ii,station in enumerate(stationlst)]

    #put information into the corresponding arrays
    for rpdict in rplst:
        station=rpdict['station']
        for kk in range(nperiod):
            ii=pdict['%2.4g' % idict[station]['period'][kk]]
            for dkey in dkeys[:-1]:
                rkey,jj=wldict[dkey]
                try:
                    rpdict[rkey][jj,ii]=idict[station][dkey][kk]
                except ValueError:
                    pass
                except IndexError:
                    rpdict[rkey][jj,ii]=1
    return idict,rplst,plst,stationlst,[title,profile,inversiontype]

def plotResponses(outputfile,maxcol=8,plottype='all',**kwargs):
    """
    plotResponse will plot the responses modeled from winglink against the
    observed data.

    Inputs:
        outputfile = full path and filename to output file
        maxcol = maximum number of columns for the plot
        plottype = 'all' to plot all on the same plot
                   '1' to plot each respones in a different figure
                   station to plot a single station or enter as a list of
                   stations to plot a few stations [station1,station2].  Does
                   not have to be verbatim but should have similar unique
                   characters input pb01 for pb01cs in outputfile
    Outputs:
        None
    """

    idict,rplst,plst,stationlst,titlelst=readOutputFile(outputfile)
    nstations=len(idict)

    #plot all responses onto one plot
    if plottype=='all':
        maxcol=8
        nrows=int(np.ceil(nstations/float(maxcol)))

        fig=plt.figure(1,[14,10])
        gs=gridspec.GridSpec(nrows,1,wspace=.15,left=.03)
        count=0
        for rr in range(nrows):
            g1=gridspec.GridSpecFromSubplotSpec(6,maxcol,subplot_spec=gs[rr],
                                                    hspace=.15,wspace=.05)
            count=rr*(maxcol)
            for cc in range(maxcol):
                try:
                    station=stationlst[count+cc]
                except IndexError:
                    break
                #plot resistivity
                axr=plt.Subplot(fig,g1[:4,cc])
                fig.add_subplot(axr)
                axr.loglog(idict[station]['period'],idict[station]['obsresxy'],
                           's',ms=2,color='b',mfc='b')
                axr.loglog(idict[station]['period'],idict[station]['modresxy'],
                           '*', ms=5,color='r',mfc='r')
                axr.loglog(idict[station]['period'],idict[station]['obsresyx'],
                           'o',ms=2,color='c',mfc='c')
                axr.loglog(idict[station]['period'],idict[station]['modresyx'],
                           'x',ms=5,color='m',mfc='m')
                axr.set_title(station+'; rms= %.2f' % idict[station]['rms'],
                              fontdict={'size':12,'weight':'bold'})
                axr.grid(True)
                axr.set_xticklabels(['' for ii in range(10)])
                if cc>0:
                    axr.set_yticklabels(['' for ii in range(6)])


                #plot phase
                axp=plt.Subplot(fig,g1[-2:,cc])
                fig.add_subplot(axp)
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['obsphasexy']),
                             's',ms=2,color='b',mfc='b')
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['modphasexy']),
                             '*',ms=5,color='r',mfc='r')
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['obsphaseyx']),
                             'o',ms=2,color='c',mfc='c')
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['modphaseyx']),
                             'x',ms=5,color='m',mfc='m')
                axp.set_ylim(0,90)
                axp.grid(True)
                axp.yaxis.set_major_locator(MultipleLocator(30))
                axp.yaxis.set_minor_locator(MultipleLocator(5))

                if cc==0 and rr==0:
                    axr.legend(['$Obs_{xy}$','$Mod_{xy}$','$Obs_{yx}$',
                                '$Mod_{yx}$'],
                                loc=2,markerscale=1,borderaxespad=.05,
                                labelspacing=.08,
                                handletextpad=.15,borderpad=.05)
                if cc==0:
                    axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                   fontdict={'size':12,'weight':'bold'})
                    axp.set_ylabel('Phase (deg)',
                                   fontdict={'size':12,'weight':'bold'})
                    axr.yaxis.set_label_coords(-.08,.5)
                    axp.yaxis.set_label_coords(-.08,.5)

                if cc>0:
                    axr.set_yticklabels(['' for ii in range(6)])
                    axp.set_yticklabels(['' for ii in range(6)])
                if rr==nrows-1:
                    axp.set_xlabel('Period (s)',
                                   fontdict={'size':12,'weight':'bold'})

    #plot each respones in a different figure
    elif plottype=='1':
        gs=gridspec.GridSpec(6,2,wspace=.05)
        for ii,station in enumerate(stationlst):
            fig=plt.figure(ii+1,[7,8])

            #plot resistivity
            axr=fig.add_subplot(gs[:4,:])

            axr.loglog(idict[station]['period'],idict[station]['obsresxy'],
                       's',ms=2,color='b',mfc='b')
            axr.loglog(idict[station]['period'],idict[station]['modresxy'],
                       '*', ms=5,color='r',mfc='r')
            axr.loglog(idict[station]['period'],idict[station]['obsresyx'],
                       'o',ms=2,color='c',mfc='c')
            axr.loglog(idict[station]['period'],idict[station]['modresyx'],
                       'x',ms=5,color='m',mfc='m')
            axr.set_title(station+'; rms= %.2f' % idict[station]['rms'],
                          fontdict={'size':12,'weight':'bold'})
            axr.grid(True)
            axr.set_xticklabels(['' for ii in range(10)])

            #plot phase
            axp=fig.add_subplot(gs[-2:,:])
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphasexy']),
                         's',ms=2,color='b',mfc='b')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphasexy']),
                         '*',ms=5,color='r',mfc='r')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphaseyx']),
                         'o',ms=2,color='c',mfc='c')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphaseyx']),
                         'x',ms=5,color='m',mfc='m')
            axp.set_ylim(0,90)
            axp.grid(True)
            axp.yaxis.set_major_locator(MultipleLocator(10))
            axp.yaxis.set_minor_locator(MultipleLocator(1))

            axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_ylabel('Phase (deg)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})
            axr.legend(['$Obs_{xy}$','$Mod_{xy}$','$Obs_{yx}$',
                                '$Mod_{yx}$'],
                                loc=2,markerscale=1,borderaxespad=.05,
                                labelspacing=.08,
                                handletextpad=.15,borderpad=.05)
            axr.yaxis.set_label_coords(-.05,.5)
            axp.yaxis.set_label_coords(-.05,.5)

    else:
        pstationlst=[]
        if type(plottype) is not list:
            plottype=[plottype]
        for station in stationlst:
            for pstation in plottype:
                if station.find(pstation)>=0:
                    print 'plotting ',station
                    pstationlst.append(station)

        gs=gridspec.GridSpec(6,2,wspace=.05,left=.1)
        for ii,station in enumerate(pstationlst):
            fig=plt.figure(ii+1,[7,7])

            #plot resistivity
            axr=fig.add_subplot(gs[:4,:])

            axr.loglog(idict[station]['period'],idict[station]['obsresxy'],
                       's',ms=2,color='b',mfc='b')
            axr.loglog(idict[station]['period'],idict[station]['modresxy'],
                       '*', ms=5,color='r',mfc='r')
            axr.loglog(idict[station]['period'],idict[station]['obsresyx'],
                       'o',ms=2,color='c',mfc='c')
            axr.loglog(idict[station]['period'],idict[station]['modresyx'],
                       'x',ms=5,color='m',mfc='m')
            axr.set_title(station+'; rms= %.2f' % idict[station]['rms'],
                          fontdict={'size':12,'weight':'bold'})
            axr.grid(True)
            axr.set_xticklabels(['' for ii in range(10)])

            #plot phase
            axp=fig.add_subplot(gs[-2:,:])
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphasexy']),
                         's',ms=2,color='b',mfc='b')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphasexy']),
                         '*',ms=5,color='r',mfc='r')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphaseyx']),
                         'o',ms=2,color='c',mfc='c')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphaseyx']),
                         'x',ms=5,color='m',mfc='m')
            axp.set_ylim(0,90)
            axp.grid(True)
            axp.yaxis.set_major_locator(MultipleLocator(10))
            axp.yaxis.set_minor_locator(MultipleLocator(1))

            axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_ylabel('Phase (deg)',
                           fontdict={'size':12,'weight':'bold'})
            axp.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})
            axr.legend(['$Obs_{xy}$','$Mod_{xy}$','$Obs_{yx}$',
                                '$Mod_{yx}$'],
                                loc=2,markerscale=1,borderaxespad=.05,
                                labelspacing=.08,
                                handletextpad=.15,borderpad=.05)
            axr.yaxis.set_label_coords(-.05,.5)
            axp.yaxis.set_label_coords(-.05,.5)

def readModelFile(modelfile,profiledirection='ew'):
    """
    readModelFile reads in the XYZ txt file output by Winglink.

    Inputs:
        modelfile = fullpath and filename to modelfile
        profiledirection = 'ew' for east-west predominantly, 'ns' for
                            predominantly north-south.  This gives column to
                            fix
    """

    mfid=open(modelfile,'r')
    lines=mfid.readlines()
    nlines=len(lines)

    X=np.zeros(nlines)
    Y=np.zeros(nlines)
    Z=np.zeros(nlines)
    rho=np.zeros(nlines)
    clst=[]
    #file starts from the bottom of the model grid in X Y Z Rho coordinates
    if profiledirection=='ew':
        for ii,line in enumerate(lines):
            linestr=line.split()
            X[ii]=float(linestr[0])
            Y[ii]=float(linestr[1])
            Z[ii]=float(linestr[2])
            rho[ii]=float(linestr[3])
            if ii>0:
                if X[ii]<X[ii-1]:
                    clst.append(ii)

    clst=np.array(clst)
    cspot=np.where(np.remainder(clst,clst[0])!=0)[0]


    return X,Y,Z,rho,clst


def readWLOutFile(outfn,ncol=5):
    """
    read .out file from winglink

    Inputs:
        outfn = full path to .out file from winglink

    Outputs:
        dx,dy,dz = cell widths in x,y,z directions (note x is to the East here
                    and y is to the north.)
    """

    wingLinkDataFH = file(outfn,'r')
    raw_data       = wingLinkDataFH.read().strip().split()

    nx = int(raw_data[0])
    ny = int(raw_data[1])
    nz = int(raw_data[2])


    dx=np.zeros(nx)
    dy=np.zeros(ny)
    dz=np.zeros(nz)

    for x_idx in range(nx):
      dx[x_idx] = raw_data[x_idx + 5]
    for y_idx in range(ny):
      dy[y_idx] = raw_data[y_idx + 5 + nx]
    for z_idx in range(nz):
      dz[z_idx] = raw_data[z_idx + 5 + nx + ny]

    #dx[0:nx/2]=-dx[0:nx/2]
    #dy[0:ny/2]=-dy[0:ny/2]


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


def readSitesFile2(sitesfn):
    """
    read sites_ file output from winglink

    Input:
        sitesfn = full path to the sites file output by winglink

    Output:
        sites_dict
        Dictionary with station names as keys. Each value is again a dictionary with keys:
        idx_east
        (index of the station-containing mesh block in east direction - starts at 1 for westernmost block)
        idx_south
        (index of the station-containing mesh block in south direction - starts at 1 for northernmost block)
        idx_z
        (index of the station-containing mesh block in down direction)

    """
    stations_dict = {}

    sfid=file(sitesfn,'r')
    slines=sfid.readlines()
    for ss in slines:
        current_dict={}
        sline   = ss.strip().split()
        station = sline[0][0:-4]
        current_dict['idx_east']  = int(sline[1])
        current_dict['idx_south'] = int(sline[2])
        current_dict['idx_z']     = int(sline[1])
        stations_dict[station] = current_dict

    return stations_dict


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

    #redundant, but necessary, since station order in sites file can be messed up
    xy_dict = {}

    for ii,sdict in enumerate(slst):
        xx=sdict['dx']
        yy=sdict['dy']
        station = sdict['station']

        if xx<nxh:
            xarr[ii]=dx[xx:nxh].sum()-dx[xx]/2
        else:
            xarr[ii]=dx[nxh:xx].sum()+dx[xx]/2
        if yy<nyh:
            yarr[ii]=-1*(dy[yy:nyh].sum()-dy[yy]/2)
        else:
            yarr[ii]=-1*(dy[nyh:yy].sum()+dy[yy]/2)


        xy_dict[station] = [xarr[ii], yarr[ii]  ]


    return xarr,yarr, xy_dict


def getmeshblockcoordinates(WL_outfile):
    """
    return a list of 3 lists, which again contain the X/Y/Z coordinate of a mesh block

    Orientation is X-North, Y-East, Z-Down.
    Horizontal origin is in the center of the mesh,
    Indexing starts at the lower left (SouthWest) corner

    Referring to a station, which has an entry in the WingLink 'Sites' file of (7,12,1), you get the coordinates as

       ( thislist[0][-12], thislist[1][6],, thislist[2][0] )

    [ explanation:
        1.Winglink refers to the UPPER left corner
        2.WingLink starts with East
        3.WingLink begins indexing at 1
    ]

    """

    east_blockwidths, north_blockwidths, z_blockwidths = readWLOutFile(WL_outfile)

    n_east  = len(east_blockwidths)
    n_north = len(north_blockwidths)
    n_down  = len(z_blockwidths)

    coord_list_xyz =[]

    total_width_ew = np.sum(east_blockwidths)
    center_ew      = total_width_ew/2.

    total_width_ns = np.sum(north_blockwidths)
    center_ns      = total_width_ns/2.

    total_depth    = np.sum(z_blockwidths)

    #depths
    lo_depths = []
    current_depth = z_blockwidths[0]/2.
    lo_depths.append(current_depth)

    for idx_z in range(n_down-1):
        current_depth += (z_blockwidths[idx_z]/2. + z_blockwidths[idx_z+1]/2.)
        lo_depths.append(current_depth)


    lo_norths = []
    current_north = north_blockwidths[0]/2.
    lo_norths.append(current_north)
    for idx_n in range(n_north-1):
        current_north += (north_blockwidths[idx_n]/2. + north_blockwidths[idx_n+1]/2.)
        lo_norths.append(current_north)

    lo_norths_centered = list(np.array(lo_norths)-center_ns)
    coord_list_xyz.append(lo_norths_centered)

    lo_easts = []
    current_east= east_blockwidths[0]/2.
    lo_easts.append(current_east)
    for idx_e in range(n_east-1):
        current_east+= (east_blockwidths[idx_e]/2. + east_blockwidths[idx_e+1]/2.)
        lo_easts.append(current_east)

    lo_easts_centered = list(np.array(lo_easts)-center_ew)
    coord_list_xyz.append(lo_easts_centered)

    coord_list_xyz.append(lo_depths)


    return coord_list_xyz

