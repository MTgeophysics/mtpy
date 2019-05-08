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

    ofid = open(outputfile, 'r')
    lines = ofid.readlines()

    idict = {}
    stationlst = []

    # get title line
    titleline = lines[1].replace('"', '')
    titleline = titleline.rstrip().split(',')
    title = titleline[1].split(':')[1]
    profile = titleline[0].split(':')[1]
    inversiontype = lines[2].rstrip()

    dkeys = ['obsresyx', 'obsphaseyx', 'modresyx', 'modphaseyx', 'obsresxy',
             'obsphasexy', 'modresxy', 'modphasexy', 'obshzres', 'obshzphase',
             'modhzres', 'modhzphase', 'period']

    for line in lines[3:]:
        if line.find('Data for station') == 0:
            station = line.rstrip().split(':')[1][1:]
            idict[station] = {}
            stationlst.append(station)
            print('Read in station: ', station)
            for key in dkeys:
                idict[station][key] = []
        elif line.find('RMS') == 0:
            idict[station]['rms'] = float(line.strip().split(' = ')[1])
        elif line.find('==') == 0:
            pass
        else:
            linelst = line.split()
            if len(linelst) == len(dkeys):
                for kk, key in enumerate(dkeys):
                    try:
                        if key.find('phase') >= 0:
                            idict[station][key].append(-1 * float(linelst[kk]))
                        else:
                            idict[station][key].append(float(linelst[kk]))
                    except ValueError:
                        idict[station][key].append(0)
            else:
                pass

    # get data into a more useful format that takes into account any masking of
    # data points.

    # get the median of period lists for survey
    plst = np.median(np.array([idict[station]['period'] for station in stationlst]),
                     axis=0)
    # length of period
    nperiod = len(plst)

    # make a dictionary of period indicies
    pdict = dict([('%2.4g' % key, ii) for ii, key in enumerate(plst)])

    # make a dictionary of indicies for spots to put res_ij and phase_ij
    wldict = {}
    for dkey in dkeys:
        if dkey[0:3].find('obs') == 0:
            wldict[dkey] = (dkey[3:], 0)
        elif dkey[0:3].find('mod') == 0:
            wldict[dkey] = (dkey[3:], 1)

    # make empty arrays to put things into
    asize = (2, nperiod)
    rplst = [{'station': station,
              'resxy': np.zeros(asize),
              'resyx': np.zeros(asize),
              'phasexy': np.zeros(asize),
              'phaseyx': np.zeros(asize),
              'hzres': np.zeros(asize),
              'hzphase': np.zeros(asize),
              } for ii, station in enumerate(stationlst)]

    # put information into the corresponding arrays
    for rpdict in rplst:
        station = rpdict['station']
        for kk in range(nperiod):
            ii = pdict['%2.4g' % idict[station]['period'][kk]]
            for dkey in dkeys[:-1]:
                rkey, jj = wldict[dkey]
                try:
                    rpdict[rkey][jj, ii] = idict[station][dkey][kk]
                except ValueError:
                    pass
                except IndexError:
                    rpdict[rkey][jj, ii] = 1
    return idict, rplst, plst, stationlst, [title, profile, inversiontype]


def plotResponses(outputfile, maxcol=8, plottype='all', **kwargs):
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

    idict, rplst, plst, stationlst, titlelst = readOutputFile(outputfile)
    nstations = len(idict)

    # plot all responses onto one plot
    if plottype == 'all':
        maxcol = 8
        nrows = int(np.ceil(nstations / float(maxcol)))

        fig = plt.figure(1, [14, 10])
        gs = gridspec.GridSpec(nrows, 1, wspace=.15, left=.03)
        count = 0
        for rr in range(nrows):
            g1 = gridspec.GridSpecFromSubplotSpec(6, maxcol, subplot_spec=gs[rr],
                                                  hspace=.15, wspace=.05)
            count = rr * (maxcol)
            for cc in range(maxcol):
                try:
                    station = stationlst[count + cc]
                except IndexError:
                    break
                # plot resistivity
                axr = plt.Subplot(fig, g1[:4, cc])
                fig.add_subplot(axr)
                axr.loglog(idict[station]['period'], idict[station]['obsresxy'],
                           's', ms=2, color='b', mfc='b')
                axr.loglog(idict[station]['period'], idict[station]['modresxy'],
                           '*', ms=5, color='r', mfc='r')
                axr.loglog(idict[station]['period'], idict[station]['obsresyx'],
                           'o', ms=2, color='c', mfc='c')
                axr.loglog(idict[station]['period'], idict[station]['modresyx'],
                           'x', ms=5, color='m', mfc='m')
                axr.set_title(station + '; rms= %.2f' % idict[station]['rms'],
                              fontdict={'size': 12, 'weight': 'bold'})
                axr.grid(True)
                axr.set_xticklabels(['' for ii in range(10)])
                if cc > 0:
                    axr.set_yticklabels(['' for ii in range(6)])

                # plot phase
                axp = plt.Subplot(fig, g1[-2:, cc])
                fig.add_subplot(axp)
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['obsphasexy']),
                             's', ms=2, color='b', mfc='b')
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['modphasexy']),
                             '*', ms=5, color='r', mfc='r')
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['obsphaseyx']),
                             'o', ms=2, color='c', mfc='c')
                axp.semilogx(idict[station]['period'],
                             np.array(idict[station]['modphaseyx']),
                             'x', ms=5, color='m', mfc='m')
                axp.set_ylim(0, 90)
                axp.grid(True)
                axp.yaxis.set_major_locator(MultipleLocator(30))
                axp.yaxis.set_minor_locator(MultipleLocator(5))

                if cc == 0 and rr == 0:
                    axr.legend(['$Obs_{xy}$', '$Mod_{xy}$', '$Obs_{yx}$',
                                '$Mod_{yx}$'],
                               loc=2, markerscale=1, borderaxespad=.05,
                               labelspacing=.08,
                               handletextpad=.15, borderpad=.05)
                if cc == 0:
                    axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                   fontdict={'size': 12, 'weight': 'bold'})
                    axp.set_ylabel('Phase (deg)',
                                   fontdict={'size': 12, 'weight': 'bold'})
                    axr.yaxis.set_label_coords(-.08, .5)
                    axp.yaxis.set_label_coords(-.08, .5)

                if cc > 0:
                    axr.set_yticklabels(['' for ii in range(6)])
                    axp.set_yticklabels(['' for ii in range(6)])
                if rr == nrows - 1:
                    axp.set_xlabel('Period (s)',
                                   fontdict={'size': 12, 'weight': 'bold'})

    # plot each respones in a different figure
    elif plottype == '1':
        gs = gridspec.GridSpec(6, 2, wspace=.05)
        for ii, station in enumerate(stationlst):
            fig = plt.figure(ii + 1, [7, 8])

            # plot resistivity
            axr = fig.add_subplot(gs[:4, :])

            axr.loglog(idict[station]['period'], idict[station]['obsresxy'],
                       's', ms=2, color='b', mfc='b')
            axr.loglog(idict[station]['period'], idict[station]['modresxy'],
                       '*', ms=5, color='r', mfc='r')
            axr.loglog(idict[station]['period'], idict[station]['obsresyx'],
                       'o', ms=2, color='c', mfc='c')
            axr.loglog(idict[station]['period'], idict[station]['modresyx'],
                       'x', ms=5, color='m', mfc='m')
            axr.set_title(station + '; rms= %.2f' % idict[station]['rms'],
                          fontdict={'size': 12, 'weight': 'bold'})
            axr.grid(True)
            axr.set_xticklabels(['' for ii in range(10)])

            # plot phase
            axp = fig.add_subplot(gs[-2:, :])
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphasexy']),
                         's', ms=2, color='b', mfc='b')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphasexy']),
                         '*', ms=5, color='r', mfc='r')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphaseyx']),
                         'o', ms=2, color='c', mfc='c')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphaseyx']),
                         'x', ms=5, color='m', mfc='m')
            axp.set_ylim(0, 90)
            axp.grid(True)
            axp.yaxis.set_major_locator(MultipleLocator(10))
            axp.yaxis.set_minor_locator(MultipleLocator(1))

            axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size': 12, 'weight': 'bold'})
            axp.set_ylabel('Phase (deg)',
                           fontdict={'size': 12, 'weight': 'bold'})
            axp.set_xlabel('Period (s)', fontdict={
                           'size': 12, 'weight': 'bold'})
            axr.legend(['$Obs_{xy}$', '$Mod_{xy}$', '$Obs_{yx}$',
                        '$Mod_{yx}$'],
                       loc=2, markerscale=1, borderaxespad=.05,
                       labelspacing=.08,
                       handletextpad=.15, borderpad=.05)
            axr.yaxis.set_label_coords(-.05, .5)
            axp.yaxis.set_label_coords(-.05, .5)

    else:
        pstationlst = []
        if type(plottype) is not list:
            plottype = [plottype]
        for station in stationlst:
            for pstation in plottype:
                if station.find(pstation) >= 0:
                    print('plotting ', station)
                    pstationlst.append(station)

        gs = gridspec.GridSpec(6, 2, wspace=.05, left=.1)
        for ii, station in enumerate(pstationlst):
            fig = plt.figure(ii + 1, [7, 7])

            # plot resistivity
            axr = fig.add_subplot(gs[:4, :])

            axr.loglog(idict[station]['period'], idict[station]['obsresxy'],
                       's', ms=2, color='b', mfc='b')
            axr.loglog(idict[station]['period'], idict[station]['modresxy'],
                       '*', ms=5, color='r', mfc='r')
            axr.loglog(idict[station]['period'], idict[station]['obsresyx'],
                       'o', ms=2, color='c', mfc='c')
            axr.loglog(idict[station]['period'], idict[station]['modresyx'],
                       'x', ms=5, color='m', mfc='m')
            axr.set_title(station + '; rms= %.2f' % idict[station]['rms'],
                          fontdict={'size': 12, 'weight': 'bold'})
            axr.grid(True)
            axr.set_xticklabels(['' for ii in range(10)])

            # plot phase
            axp = fig.add_subplot(gs[-2:, :])
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphasexy']),
                         's', ms=2, color='b', mfc='b')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphasexy']),
                         '*', ms=5, color='r', mfc='r')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['obsphaseyx']),
                         'o', ms=2, color='c', mfc='c')
            axp.semilogx(idict[station]['period'],
                         np.array(idict[station]['modphaseyx']),
                         'x', ms=5, color='m', mfc='m')
            axp.set_ylim(0, 90)
            axp.grid(True)
            axp.yaxis.set_major_locator(MultipleLocator(10))
            axp.yaxis.set_minor_locator(MultipleLocator(1))

            axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                           fontdict={'size': 12, 'weight': 'bold'})
            axp.set_ylabel('Phase (deg)',
                           fontdict={'size': 12, 'weight': 'bold'})
            axp.set_xlabel('Period (s)', fontdict={
                           'size': 12, 'weight': 'bold'})
            axr.legend(['$Obs_{xy}$', '$Mod_{xy}$', '$Obs_{yx}$',
                        '$Mod_{yx}$'],
                       loc=2, markerscale=1, borderaxespad=.05,
                       labelspacing=.08,
                       handletextpad=.15, borderpad=.05)
            axr.yaxis.set_label_coords(-.05, .5)
            axp.yaxis.set_label_coords(-.05, .5)


def readModelFile(modelfile, profiledirection='ew'):
    """
    readModelFile reads in the XYZ txt file output by Winglink.    

    Inputs:
        modelfile = fullpath and filename to modelfile
        profiledirection = 'ew' for east-west predominantly, 'ns' for 
                            predominantly north-south.  This gives column to 
                            fix
    """

    mfid = open(modelfile, 'r')
    lines = mfid.readlines()
    nlines = len(lines)

    X = np.zeros(nlines)
    Y = np.zeros(nlines)
    Z = np.zeros(nlines)
    rho = np.zeros(nlines)
    clst = []
    # file starts from the bottom of the model grid in X Y Z Rho coordinates
    if profiledirection == 'ew':
        for ii, line in enumerate(lines):
            linestr = line.split()
            X[ii] = float(linestr[0])
            Y[ii] = float(linestr[1])
            Z[ii] = float(linestr[2])
            rho[ii] = float(linestr[3])
            if ii > 0:
                if X[ii] < X[ii - 1]:
                    clst.append(ii)

    clst = np.array(clst)
    cspot = np.where(np.remainder(clst, clst[0]) != 0)[0]

    return X, Y, Z, rho, clst
