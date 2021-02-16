# -*- coding: utf-8 -*-
"""
Created on Mon Apr 02 11:54:33 2012

@author: a1185872
"""

import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Ellipse, Rectangle, Arrow
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import mtpy1.core.z as Z
import mtpy1.utils.latlongutmconversion as ll2utm

# tolerance to find frequencies
ptol = 0.15

# error of data in percentage
zerr = 0.05
# errormap values which is multiplied by zerr to get a total error
zxxerrmap = 10
zxyerrmap = 1
zyxerrmap = 1
zyyerrmap = 10
zerrmap = [zxxerrmap, zxyerrmap, zyxerrmap, zyyerrmap]

# ==============================================================================
# Colormaps for plots
# ==============================================================================
# phase tensor map
ptcmapdict = {
    "red": ((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
    "green": ((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)),
    "blue": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
}
ptcmap = LinearSegmentedColormap("ptcmap", ptcmapdict, 256)

# phase tensor map for difference (reverse)
ptcmapdictr = {
    "red": ((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)),
    "green": ((0.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
    "blue": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
}
ptcmapr = LinearSegmentedColormap("ptcmapr", ptcmapdictr, 256)

# resistivity tensor map for calculating delta
ptcmapdict2 = {
    "red": ((0.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
    "green": ((0.0, 0.5, 0.5), (1.0, 0.5, 0.5)),
    "blue": ((0.0, 0.5, 0.5), (1.0, 0.5, 0.5)),
}
ptcmap2 = LinearSegmentedColormap("ptcmap2", ptcmapdict2, 256)

# resistivity tensor map for calcluating resistivity difference
rtcmapdict = {
    "red": ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 0.0)),
    "green": ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
    "blue": ((0.0, 0.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
}
rtcmap = LinearSegmentedColormap("rtcmap", rtcmapdict, 256)

# resistivity tensor map for calcluating apparent resistivity
rtcmapdictr = {
    "red": ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
    "green": ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
    "blue": ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)),
}
rtcmapr = LinearSegmentedColormap("rtcmapr", rtcmapdictr, 256)

# ==============================================================================
#  define some helping functions
# ==============================================================================
# make a class to pick periods


class ListPeriods:
    def __init__(self, fig):
        self.plst = []
        self.fig = fig
        self.count = 1

    def connect(self):
        self.cid = self.fig.canvas.mpl_connect("button_press_event", self.onclick)

    def onclick(self, event):
        print "{0} Period: {1:.5g}".format(self.count, event.xdata)
        self.plst.append(event.xdata)
        self.count += 1

    def disconnect(self):
        self.fig.canvas.mpl_disconnect(self.cid)


def readWLOutFile(outfn, ncol=5):
    """
    read .out file from winglink

    Inputs:
        outfn = full path to .out file from winglink

    Outputs:
        dx,dy,dz = cell nodes in x,y,z directions (note x is to the East here
                    and y is to the north.)
    """

    wingLinkDataFH = file(outfn, "r")
    raw_data = wingLinkDataFH.read().strip().split()

    nx = int(raw_data[0])
    ny = int(raw_data[1])
    nz = int(raw_data[2])

    dx = np.zeros(nx)
    dy = np.zeros(ny)
    dz = np.zeros(nz)

    for x_idx in range(nx):
        dx[x_idx] = raw_data[x_idx + 5]
    for y_idx in range(ny):
        dy[y_idx] = raw_data[y_idx + 5 + nx]
    for z_idx in range(nz):
        dz[z_idx] = raw_data[z_idx + 5 + nx + ny]

    # dx[0:nx/2]=-dx[0:nx/2]
    # dy[0:ny/2]=-dy[0:ny/2]

    return dx, dy, dz


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

    sfid = file(sitesfn, "r")
    slines = sfid.readlines()

    slst = []
    sitelst = []
    for ss in slines:
        sdict = {}
        sline = ss.strip().split()
        sdict["station"] = sline[0][0:-4]
        sdict["dx"] = int(sline[1]) - 1
        sdict["dy"] = int(sline[2]) - 1
        sdict["dz"] = int(sline[3]) - 1
        sdict["something"] = int(sline[4])
        sdict["number"] = int(sline[5])
        slst.append(sdict)
        sitelst.append(sline[0][0:-4])
    return slst, sitelst


def getXY(sitesfn, outfn, ncol=5):
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

    slst, sitelst = readSitesFile(sitesfn)

    dx, dy, dz = readWLOutFile(outfn, ncol=ncol)

    ns = len(slst)
    nxh = len(dx) / 2
    nyh = len(dy) / 2
    xarr = np.zeros(ns)
    yarr = np.zeros(ns)

    for ii, sdict in enumerate(slst):
        xx = sdict["dx"]
        yy = sdict["dy"]
        if xx < nxh:
            xarr[ii] = dx[xx:nxh].sum() - dx[xx] / 2
        else:
            xarr[ii] = dx[nxh:xx].sum() + dx[xx] / 2
        if yy < nyh:
            yarr[ii] = -1 * (dy[yy:nyh].sum() - dy[yy] / 2)
        else:
            yarr[ii] = -1 * (dy[nyh:yy].sum() + dy[yy] / 2)

    return xarr, yarr


def getPeriods(edilst, errthresh=10):
    """
    Plots periods for all stations in edipath and the plot is interactive, just
    click on the period you want to select and it will appear in the console,
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

    plt.rcParams["font.size"] = 10
    plt.rcParams["figure.subplot.left"] = 0.13
    plt.rcParams["figure.subplot.right"] = 0.98
    plt.rcParams["figure.subplot.bottom"] = 0.1
    plt.rcParams["figure.subplot.top"] = 0.95
    plt.rcParams["figure.subplot.wspace"] = 0.25
    plt.rcParams["figure.subplot.hspace"] = 0.05

    periodlst = []
    errorlst = []

    fig1 = plt.figure(5)
    ax = fig1.add_subplot(1, 1, 1)
    for edi in edilst:
        if not os.path.isfile(edi):
            print "Could not find " + edi
        else:
            z1 = Z.Z(edi)
            periodlst.append(z1.period)
            zdet = np.array([np.sqrt(abs(np.linalg.det(zz))) for zz in z1.z])
            error = np.array([np.sqrt(abs(np.linalg.det(zz))) for zz in z1.zvar])
            perror = (error / zdet) * 100
            errorlst.append(perror)
            # make a plot to pick frequencies from showing period and percent
            # error
            ax.scatter(z1.period, perror, marker="x", picker=5)
            pfind = np.where(perror > errthresh)[0]
            if len(pfind) > 0:
                print "Error greater than {0:.3f} for ".format(errthresh) + z1.station
                for jj in pfind:
                    ax.scatter(z1.period[jj], perror[jj], marker="x", color="r")
                    ax.text(
                        z1.period[jj],
                        perror[jj] * 1.05,
                        z1.station,
                        horizontalalignment="center",
                        verticalalignment="baseline",
                        fontdict={"size": 8, "color": "red"},
                    )
                    print jj, z1.period[jj]

    ax.set_xscale("log")
    ax.set_xlim(
        10 ** np.floor(np.log10(z1.period[0])), 10 ** np.ceil(np.log10(z1.period[-1]))
    )
    ax.set_ylim(0, 3 * errthresh)
    ax.set_yscale("log")
    ax.set_xlabel("Period (s)", fontdict={"size": 12, "weight": "bold"})
    ax.set_ylabel("Percent Error", fontdict={"size": 12, "weight": "bold"})
    ax.grid("on", which="both")

    lp = ListPeriods(fig1)
    lp.connect()

    plt.show()

    return periodlst, errorlst, lp


def make3DGrid(
    edilst,
    xspacing=500,
    yspacing=500,
    z1layer=10,
    xpad=5,
    ypad=5,
    zpad=5,
    xpadroot=5,
    ypadroot=5,
    zpadroot=2,
    zpadpow=(5, 15),
    nz=30,
    plotyn="y",
    plotxlimits=None,
    plotylimits=None,
    plotzlimits=None,
):
    """
    makes a grid from the edifiles to go into wsinv3d.  The defaults usually
    work relatively well, but it might take some effort to get a desired grid.

    Inputs:
    --------
        **edilst** : list
                     list of full paths to the .edi files to be included in
                     the inversion.

        **xspacing** : float
                       spacing of cells in the east-west direction in meters.
                       *default* is 500 (m)

        **yspacing** : float
                       spacing of cells in the north-south direction in meters.
                       *default* is 500 (m)

        **z1layer** : float
                      the depth of the first layer in the model in meters.
                      This is usually about 1/10th of your shallowest skin
                      depth.
                      *default* is 10 (m)

        **xpad** : int
                   number of cells to pad on either side in the east-west
                   direction.  The width of these cells grows exponentially
                   to the edge.
                   *default* is 5

        **ypad** : int
                   number of cells to pad on either side in the north-south
                   direction.  The width of these cells grows exponentially
                   to the edge.
                   *default* is 5

        **zpad** : int
                   number of cells to pad on either side in the vertical
                   direction.  This is to pad beneath the depth of
                   investigation and grows faster exponentially than the zone
                   of study.  The purpose is to decrease the number of cells
                   in the model.
                   *default* is 5

        **xpadroot** : float
                       the root number that is multiplied to itself for
                       calculating the width of the padding cells in the
                       east-west direction.
                       *default* is 5

        **ypadroot** : float
                       the root number that is multiplied to itself for
                       calculating the width of the padding cells in the
                       north-south direction.
                       *default* is 5

        **zpadroot** : float
                       the root number that is multiplied to itself for
                       calculating the width of the padding cells in the
                       vertical direction.
                       *default* is 2

        **zpadpow** : tuple (min,max)
                      the power to which zpadroot is raised for the padding
                      cells in the vertical direction.  Input as a tuple with
                      minimum power and maximum power.
                      *default* is (5,15)

        **nz** : int
                 number of layers in the vertical direction.  Remember that
                 the inversion code automatically adds 7 air layers to the
                 model which need to be used when estimating the memory that
                 it is going to take to run the model.
                 *default* is 30

        **plotyn** : [ 'y' | 'n' ]
                     if plotyn=='y' then a plot showing map view (east:north)
                     and a cross sectional view (east:vertical) plane

                     * 'y' to plot the grid with station locations

                     * 'n' to suppress the plotting.

        **plotxlimits** : tuple (xmin,xmax)
                         plot min and max distances in meters for the east-west
                         direction.  If not input, the xlimits will be set to
                         the furthest stations east and west.
                         *default* is None

        **plotylimits** : tuple (ymin,ymax)
                         plot min and max distances in meters for the east-west
                         direction. If not input, the ylimits will be set to
                         the furthest stations north and south.
                         *default* is None

        **plotzlimits** : tuple (zmin,zmax)
                         plot min and max distances in meters for the east-west
                         direction.  If not input, the zlimits will be set to
                         the nz layer and 0.
                         *default* is None

    Returns:
    --------
        xgrid,ygrid,zgrid,locations,slst
        **xgrid** : np.array
                    array of the east-west cell locations

        **ygrid** : np.array
                    array of the north-south cell locations

        **zgrid** : np.array
                    array of the vertical cell locations

        **locations** : np.array (ns,2)
                        array of station locations placed in the center of
                        the cells.
                        * column 1 is for east-west locations
                        * column 2 is for the north-south location

        **slst** : list
                   list of dictionaries for each station with keys:
                       * *'station'* for the station name
                       * *'east'* for easting in model coordinates
                       * *'east_c'* for easting in model coordinates to place
                                    the station at the center of the cell
                       * *'north'* for northing in model coordinates
                       * *'north_c'* for northing in model coordinates to place
                                    the station at the center of the cell


    :Example: ::

        >>> import mtpy.modeling.ws3dtools as ws
        >>> import os
        >>> edipath=r"/home/edifiles"
        >>> edilst=[os.path.join(edipath,edi) for os.listdir(edipath)]
        >>> xg,yg,zg,loc,statlst=ws.make3DGrid(edilst,plotzlimits=(-2000,200))

    """
    ns = len(edilst)
    locations = np.zeros((ns, 2))
    slst = []
    for ii, edi in enumerate(edilst):
        zz = Z.Z(edi)
        zone, east, north = ll2utm.LLtoUTM(23, zz.lat, zz.lon)
        locations[ii, 0] = east
        locations[ii, 1] = north
        slst.append({"station": zz.station, "east": east, "north": north})

    # estimate the mean distance to  get into relative coordinates
    xmean = locations[:, 0].mean()
    ymean = locations[:, 1].mean()

    # remove the average distance to get coordinates in a relative space
    locations[:, 0] -= xmean
    locations[:, 1] -= ymean
    for sdict in slst:
        sdict["east"] -= xmean
        sdict["north"] -= ymean

    # translate the stations so they are relative to 0,0
    xcenter = (locations[:, 0].max() - np.abs(locations[:, 0].min())) / 2
    ycenter = (locations[:, 1].max() - np.abs(locations[:, 1].min())) / 2

    # remove the average distance to get coordinates in a relative space
    locations[:, 0] -= xcenter
    locations[:, 1] -= ycenter
    for sdict in slst:
        sdict["east"] -= xcenter
        sdict["north"] -= ycenter

    # pickout the furtherst south and west locations
    # and put that station as the bottom left corner of the main grid
    xleft = locations[:, 0].min() - xspacing / 2
    xright = locations[:, 0].max() + xspacing / 2
    ybottom = locations[:, 1].min() - yspacing / 2
    ytop = locations[:, 1].max() + yspacing / 2

    # ---make a grid around the stations from the parameters above---
    # make grid in east-west direction
    midxgrid = np.arange(start=xleft, stop=xright + xspacing, step=xspacing)
    xpadleft = (
        np.round(-xspacing * 5 ** np.arange(start=0.5, stop=3, step=3.0 / xpad)) + xleft
    )
    xpadright = (
        np.round(xspacing * 5 ** np.arange(start=0.5, stop=3, step=3.0 / xpad)) + xright
    )
    xgridr = np.append(np.append(xpadleft[::-1], midxgrid), xpadright)

    # make grid in north-south direction
    midygrid = np.arange(start=ybottom, stop=ytop + yspacing, step=yspacing)
    ypadbottom = (
        np.round(-yspacing * 5 ** np.arange(start=0.5, stop=3, step=3.0 / xpad))
        + ybottom
    )
    ypadtop = (
        np.round(yspacing * 5 ** np.arange(start=0.5, stop=3, step=3.0 / xpad)) + ytop
    )
    ygridr = np.append(np.append(ypadbottom[::-1], midygrid), ypadtop)

    # make depth grid
    zgrid1 = z1layer * 2 ** np.round(np.arange(0, zpadpow[0], zpadpow[0] / (nz - zpad)))
    zgrid2 = z1layer * 2 ** np.round(
        np.arange(zpadpow[0], zpadpow[1], (zpadpow[1] - zpadpow[0]) / (zpad))
    )

    zgrid = np.append(zgrid1, zgrid2)

    # --Need to make an array of the individual cell dimensions for the wsinv3d
    xnodes = xgridr.copy()
    nx = xgridr.shape[0]
    xnodes[: nx / 2] = np.array(
        [abs(xgridr[ii] - xgridr[ii + 1]) for ii in range(int(nx / 2))]
    )
    xnodes[nx / 2 :] = np.array(
        [abs(xgridr[ii] - xgridr[ii + 1]) for ii in range(int(nx / 2) - 1, nx - 1)]
    )

    ynodes = ygridr.copy()
    ny = ygridr.shape[0]
    ynodes[: ny / 2] = np.array(
        [abs(ygridr[ii] - ygridr[ii + 1]) for ii in range(int(ny / 2))]
    )
    ynodes[ny / 2 :] = np.array(
        [abs(ygridr[ii] - ygridr[ii + 1]) for ii in range(int(ny / 2) - 1, ny - 1)]
    )

    # --put the grids into coordinates relative to the center of the grid
    xgrid = xnodes.copy()
    xgrid[: int(nx / 2)] = -np.array(
        [xnodes[ii : int(nx / 2)].sum() for ii in range(int(nx / 2))]
    )
    xgrid[int(nx / 2) :] = (
        np.array([xnodes[int(nx / 2) : ii + 1].sum() for ii in range(int(nx / 2), nx)])
        - xnodes[int(nx / 2)]
    )

    ygrid = ynodes.copy()
    ygrid[: int(ny / 2)] = -np.array(
        [ynodes[ii : int(ny / 2)].sum() for ii in range(int(ny / 2))]
    )
    ygrid[int(ny / 2) :] = (
        np.array([ynodes[int(ny / 2) : ii + 1].sum() for ii in range(int(ny / 2), ny)])
        - ynodes[int(ny / 2)]
    )

    # make sure that the stations are in the center of the cell as requested by
    # the code.
    for sdict in slst:
        # look for the closest grid line
        xx = [
            nn
            for nn, xf in enumerate(xgrid)
            if xf > (sdict["east"] - xspacing) and xf < (sdict["east"] + xspacing)
        ]

        # shift the station to the center in the east-west direction
        if xgrid[xx[0]] < sdict["east"]:
            sdict["east_c"] = xgrid[xx[0]] + xspacing / 2
        elif xgrid[xx[0]] > sdict["east"]:
            sdict["east_c"] = xgrid[xx[0]] - xspacing / 2

        # look for closest grid line
        yy = [
            mm
            for mm, yf in enumerate(ygrid)
            if yf > (sdict["north"] - yspacing) and yf < (sdict["north"] + yspacing)
        ]

        # shift station to center of cell in north-south direction
        if ygrid[yy[0]] < sdict["north"]:
            sdict["north_c"] = ygrid[yy[0]] + yspacing / 2
        elif ygrid[yy[0]] > sdict["north"]:
            sdict["north_c"] = ygrid[yy[0]] - yspacing / 2

    # =Plot the data if desired=========================
    if plotyn == "y":
        fig = plt.figure(1, figsize=[10, 10], dpi=300)

        # ---plot map view
        ax1 = fig.add_subplot(1, 2, 1, aspect="equal")

        for sdict in slst:
            # make sure the station is in the center of the cell

            ax1.scatter(sdict["east_c"], sdict["north_c"], marker="v")

        for xp in xgrid:
            ax1.plot([xp, xp], [ygrid.min(), ygrid.max()], color="k")

        for yp in ygrid:
            ax1.plot([xgrid.min(), xgrid.max()], [yp, yp], color="k")

        if plotxlimits is None:
            ax1.set_xlim(
                locations[:, 0].min() - 10 * xspacing,
                locations[:, 0].max() + 10 * xspacing,
            )
        else:
            ax1.set_xlim(plotxlimits)

        if plotylimits is None:
            ax1.set_ylim(
                locations[:, 1].min() - 50 * yspacing,
                locations[:, 1].max() + 50 * yspacing,
            )
        else:
            ax1.set_ylim(plotylimits)

        ax1.set_ylabel("Northing (m)", fontdict={"size": 10, "weight": "bold"})
        ax1.set_xlabel("Easting (m)", fontdict={"size": 10, "weight": "bold"})

        # ----plot depth view
        ax2 = fig.add_subplot(1, 2, 2, aspect="auto")

        for xp in xgrid:
            ax2.plot([xp, xp], [-zgrid.sum(), 0], color="k")

        for sdict in slst:
            ax2.scatter(sdict["east_c"], 0, marker="v")

        for zz, zp in enumerate(zgrid):
            ax2.plot(
                [xgrid.min(), xgrid.max()],
                [-zgrid[0:zz].sum(), -zgrid[0:zz].sum()],
                color="k",
            )

        if plotzlimits is None:
            ax2.set_ylim(-zgrid1.max(), 200)
        else:
            ax2.set_ylim(plotzlimits)

        if plotxlimits is None:
            ax2.set_xlim(
                locations[:, 0].min() - xspacing, locations[:, 0].max() + xspacing
            )
        else:
            ax2.set_xlim(plotxlimits)

        ax2.set_ylabel("Depth (m)", fontdict={"size": 10, "weight": "bold"})
        ax2.set_xlabel("Easting (m)", fontdict={"size": 10, "weight": "bold"})

        plt.show()

    print "-" * 15
    print "   Number of stations = {0}".format(len(slst))
    print "   Dimensions: "
    print "      e-w = {0}".format(xgrid.shape[0])
    print "      n-s = {0}".format(ygrid.shape[0])
    print "       z  = {0}".format(zgrid.shape[0])
    print "   Extensions: "
    print "      e-w = {0:.1f} (m)".format(xgrid.__abs__().sum())
    print "      n-s = {0:.1f} (m)".format(ygrid.__abs__().sum())
    print "      0-z = {0:.1f} (m)".format(zgrid.__abs__().sum())
    print "-" * 15
    return ynodes, xnodes, zgrid, locations, slst


def writeWSDataFile(
    periodlst,
    edilst,
    sitesfn=None,
    outfn=None,
    sitelocations=None,
    zerr=0.05,
    ptol=0.15,
    zerrmap=[10, 1, 1, 10],
    savepath=None,
    ncol=5,
    units="mv",
):
    """
    writes a data file for WSINV3D from winglink outputs

    Inputs:
    --------
        **periodlst** :list
                        periods to extract from edifiles, can get them from
                        using the function getPeriods.

        **edilst** : list
                    list of full paths to .edi files to use for inversion

        **sitelocations**  : np.array (ns,2)
                            array of station locations where [:,0] corresponds
                            to the east-west location and [:,1] corresponds to
                            the north-south location.  This can be found from
                            Make3DGrid.  Locations are in meters in grid
                            coordinates.

        **sitesfn** : string
                     if you used Winglink to make the model then you need to
                     input the sites filename (full path)

        **outfn** : string
                    if you used Winglink to make the model need to input the
                    winglink .out file (full path)

        **savepath** : string
                       directory or full path to save data file to, default
                       path is dirname sitesfn.
                       saves as: savepath/WSDataFile.dat
                       *Need to input if you did not use Winglink*

        **zerr** : float
                  percent error to give to impedance tensor components in
                  decimal form --> 10% = 0.10
                  *default* is .05

        **ptol** : float
                   percent tolerance to locate frequencies in case edi files
                   don't have the same frequencies.  Need to add interpolation.
                   *default* is 0.15

        **zerrmap** :  tuple (zxx,zxy,zyx,zyy)
                       multiple to multiply err of zxx,zxy,zyx,zyy by.
                       Note the total error is zerr*zerrmap[ii]

        **ncol** : int
                   number of columns in outfn, sometimes it outputs different
                   number of columns.


    Returns:
    --------

        **datafn** : full path to data file, saved in dirname(sitesfn) or
                     savepath where savepath can be a directory or full
                     filename
    """

    ns = len(edilst)

    # get units correctly
    if units == "mv":
        zconv = 1.0 / 796.0

    # create the output filename
    if savepath is None:
        ofile = os.path.join(os.path.dirname(sitesfn), "WSDataFile.dat")
    elif savepath.find(".") == -1:
        ofile = os.path.join(savepath, "WSDataFile.dat")
    else:
        ofile = savepath

    # if there is a site file from someone who naively used winglink
    if sitesfn is not None:
        # read in stations from sites file
        sitelst, slst = readSitesFile(sitesfn)

        # get x and y locations on a relative grid
        xlst, ylst = getXY(sitesfn, outfn, ncol=ncol)

    # if the user made a grid in python or some other fashion
    if sitelocations is not None:
        if isinstance(sitelocations[0], dict):
            xlst = np.zeros(ns)
            ylst = np.zeros(ns)
            slst = []
            for dd, sd in enumerate(sitelocations):
                xlst[dd] = sd["east_c"]
                ylst[dd] = sd["north_c"]
                slst.append(sd["station"])
        else:
            xlst = sitelocations[:, 0]
            ylst = sitelocations[:, 1]

    # define some lengths
    nperiod = len(periodlst)

    # make an array to put data into for easy writing
    zarr = np.zeros((ns, nperiod, 4), dtype="complex")

    # --------find frequencies-------------------------------------------------
    linelst = []
    for ss, edi in enumerate(edilst):
        if not os.path.isfile(edi):
            raise IOError("Could not find " + edi)

        z1 = Z.Z(edi)
        sdict = {}
        fspot = {}
        for ff, f1 in enumerate(periodlst):
            for kk, f2 in enumerate(z1.period):
                if f2 >= (1 - ptol) * f1 and f2 <= (1 + ptol) * f1:
                    zderr = np.array(
                        [
                            abs(z1.zvar[kk, nn, mm]) / abs(z1.z[kk, nn, mm]) * 100
                            for nn in range(2)
                            for mm in range(2)
                        ]
                    )
                    fspot["{0:.6g}".format(f1)] = (
                        kk,
                        f2,
                        zderr[0],
                        zderr[1],
                        zderr[2],
                        zderr[3],
                    )
                    zarr[ss, ff, :] = z1.z[kk].reshape(4,)

        print z1.station, len(fspot)
        sdict["fspot"] = fspot
        sdict["station"] = z1.station
        linelst.append(sdict)

    # -----Write data file-----------------------------------------------------

    ofid = file(ofile, "w")
    ofid.write("{0:d} {1:d} {2:d}\n".format(ns, nperiod, 8))

    # write N-S locations
    ofid.write("Station_Location: N-S \n")
    for ii in range(ns / 8 + 1):
        for ll in range(8):
            try:
                ofid.write("{0:+.4e} ".format(ylst[ii * 8 + ll]))
            except IndexError:
                pass
        ofid.write("\n")

    # write E-W locations
    ofid.write("Station_Location: E-W \n")
    for ii in range(ns / 8 + 1):
        for ll in range(8):
            try:
                ofid.write("{0:+.4e} ".format(xlst[ii * 8 + ll]))
            except IndexError:
                pass
        ofid.write("\n")

    # write impedance tensor components
    for ii, p1 in enumerate(periodlst):
        ofid.write("DATA_Period: {0:3.6f}\n".format(p1))
        for ss in range(ns):
            zline = zarr[ss, ii, :]
            for jj in range(4):
                ofid.write("{0:+.4e} ".format(zline[jj].real * zconv))
                ofid.write("{0:+.4e} ".format(-zline[jj].imag * zconv))
            ofid.write("\n")

    # write error as a percentage of Z
    for ii, p1 in enumerate(periodlst):
        ofid.write("ERROR_Period: {0:3.6f}\n".format(p1))
        for ss in range(ns):
            zline = zarr[ss, ii, :]
            for jj in range(4):
                ofid.write("{0:+.4e} ".format(zline[jj].real * zerr * zconv))
                ofid.write("{0:+.4e} ".format(zline[jj].imag * zerr * zconv))
            ofid.write("\n")

    # write error maps
    for ii, p1 in enumerate(periodlst):
        ofid.write("ERMAP_Period: {0:3.6f}\n".format(p1))
        for ss in range(ns):
            zline = zarr[ss, ii, :]
            for jj in range(4):
                ofid.write("{0:.5e} ".format(zerrmap[jj]))
                ofid.write("{0:.5e} ".format(zerrmap[jj]))
            ofid.write("\n")
    ofid.close()
    print "Wrote file to: " + ofile

    # write out places where errors are larger than error tolerance
    errfid = file(os.path.join(os.path.dirname(ofile), "DataErrorLocations.txt"), "w")
    errfid.write("Errors larger than error tolerance of: \n")
    errfid.write(
        "Zxx={0} Zxy={1} Zyx={2} Zyy={3} \n".format(
            zerrmap[0] * zerr, zerrmap[1] * zerr, zerrmap[2] * zerr, zerrmap[3] * zerr
        )
    )
    errfid.write("-" * 20 + "\n")
    errfid.write("station  T=period(s) Zij err=percentage \n")
    for pfdict in linelst:
        for kk, ff in enumerate(pfdict["fspot"]):
            if pfdict["fspot"][ff][2] > zerr * 100 * zerrmap[0]:
                errfid.write(
                    pfdict["station"]
                    + "  T="
                    + ff
                    + " Zxx err={0:.3f} \n".format(pfdict["fspot"][ff][2])
                )
            if pfdict["fspot"][ff][3] > zerr * 100 * zerrmap[1]:
                errfid.write(
                    pfdict["station"]
                    + "  T="
                    + ff
                    + " Zxy err={0:.3f} \n".format(pfdict["fspot"][ff][3])
                )
            if pfdict["fspot"][ff][4] > zerr * 100 * zerrmap[2]:
                errfid.write(
                    pfdict["station"]
                    + "  T="
                    + ff
                    + " Zyx err={0:.3f} \n".format(pfdict["fspot"][ff][4])
                )
            if pfdict["fspot"][ff][5] > zerr * 100 * zerrmap[3]:
                errfid.write(
                    pfdict["station"]
                    + "  T="
                    + ff
                    + " Zyy err={0:.3f} \n".format(pfdict["fspot"][ff][5])
                )
    errfid.close()
    print "Wrote errors lager than tolerance to: "
    print os.path.join(os.path.dirname(ofile), "DataErrorLocations.txt")

    return ofile, linelst


def writeInit3DFile_wl(outfn, rhostart=100, ncol=5, savepath=None):
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

    # create the output filename
    if savepath is None:
        ifile = os.path.join(os.path.dirname(outfn), "init3d")
    elif savepath.find(".") == -1:
        ifile = os.path.join(savepath, "init3d")
    else:
        ifile = savepath

    dx, dy, dz = readWLOutFile(outfn, ncol=ncol)

    nx = len(dx)
    ny = len(dy)
    nz = len(dz)

    init_modelFH = open(ifile, "w")
    init_modelFH.write("#Initial model \n")
    init_modelFH.write("%i %i %i 1 \n" % (ny, nx, nz))

    # write y locations
    y_string = ""
    y_counter = 0
    for y_idx in range(ny):
        y_string += "%.3e  " % (dy[y_idx])
        y_counter += 1
        if y_counter == 8:
            y_string += "\n"
            y_counter = 0
    if ny % 8:
        y_string += "\n"
    init_modelFH.write(y_string)

    # write x locations
    x_string = ""
    x_counter = 0
    for x_idx in range(nx):
        x_string += "%.3e  " % (dx[x_idx])
        x_counter += 1
        if x_counter == 8:
            x_string += "\n"
            x_counter = 0
    if nx % 8:
        x_string += "\n"
    init_modelFH.write(x_string)

    # write z locations
    z_string = ""
    z_counter = 0
    for z_idx in range(nz):
        z_string += "%.3e  " % (dz[z_idx])
        z_counter += 1
        if z_counter == 8:
            z_string += "\n"
            z_counter = 0
    if nz % 8:
        z_string += "\n"
    init_modelFH.write(z_string)

    init_modelFH.write("%i \n" % int(rhostart))

    init_modelFH.close()

    print "Wrote init file to: " + ifile

    return ifile


def writeInit3DFile(
    xgrid,
    ygrid,
    zgrid,
    savepath,
    reslst=100,
    title="Initial File for WSINV3D",
    resmodel=None,
):
    """
    will write an initial file for wsinv3d.  At the moment can only make a
    layered model that can then be manipulated later.  Input for a layered
    model is in layers which is [(layer1,layer2,resistivity index for reslst)]

    Note that x is assumed to be S --> N, y is assumed to be W --> E and
    z is positive downwards.

    Also, the xgrid, ygrid and zgrid are assumed to be the relative distance
    between neighboring nodes.  This is needed because wsinv3d builds the
    model from the bottom NW corner assuming the cell width from the init file.

    Therefore the first line or index=0 is the southern most row of cells, so
    if you build a model by hand the the layer block will look upside down if
    you were to picture it in map view. Confusing, perhaps, but that is the
    way it is.

    Argumens:
    ----------

        **xgrid** : np.array(nx)
                    block dimensions (m) in the N-S direction. **Note** that
                    the code reads the grid assuming that index=0 is the
                    southern most point.

        **ygrid** : np.array(ny)
                    block dimensions (m) in the E-W direction.  **Note** that
                    the code reads in the grid assuming that index=0 is the
                    western most point.

        **zgrid** : np.array(nz)
                    block dimensions (m) in the vertical direction.  This is
                    positive downwards.

        **savepath** : string
                      Path to the director where the initial file will be saved
                      as savepath/init3d

        **reslst** : float or list
                    The start resistivity as a float or a list of resistivities
                    that coorespond to the starting resistivity model
                    **resmodel**.  This must be input if you input **resmodel**

        **title** : string
                    Title that goes into the first line of savepath/init3d

        **resmodel** : np.array((nx,ny,nz))
                        Starting resistivity model.  Each cell is allocated an
                        integer value that cooresponds to the index value of
                        **reslst**.  **Note** again that the modeling code
                        assumes that the first row it reads in is the southern
                        most row and the first column it reads in is the
                        western most column.  Similarly, the first plane it
                        reads in is the Earth's surface.

    Returns:
    --------

        **initfn** : full path to initial file



    """

    if not isinstance(reslst, list):
        reslst = [reslst]

    if os.path.isdir(savepath) == True:
        ifn = os.path.join(savepath, "init3d")

    else:
        ifn = os.path.join(savepath)

    ifid = file(ifn, "w")
    ifid.write("# " + title + "\n".upper())
    ifid.write(
        "{0} {1} {2} {3}\n".format(
            xgrid.shape[0], ygrid.shape[0], zgrid.shape[0], len(reslst)
        )
    )

    # write S --> N node block
    for ii, xx in enumerate(xgrid):
        ifid.write("{0:>12}".format("{:.1f}".format(abs(xx))))
        if ii != 0 and np.remainder(ii + 1, 5) == 0:
            ifid.write("\n")
        elif ii == xgrid.shape[0] - 1:
            ifid.write("\n")

    # write W --> E node block
    for jj, yy in enumerate(ygrid):
        ifid.write("{0:>12}".format("{:.1f}".format(abs(yy))))
        if jj != 0 and np.remainder(jj + 1, 5) == 0:
            ifid.write("\n")
        elif jj == ygrid.shape[0] - 1:
            ifid.write("\n")

    # write top --> bottom node block
    for kk, zz in enumerate(zgrid):
        ifid.write("{0:>12}".format("{:.1f}".format(abs(zz))))
        if kk != 0 and np.remainder(kk + 1, 5) == 0:
            ifid.write("\n")
        elif kk == zgrid.shape[0] - 1:
            ifid.write("\n")

    # write the resistivity list
    for ff in reslst:
        ifid.write("{0:.1f} ".format(ff))
    ifid.write("\n")

    #    else:
    if resmodel is None:
        ifid.close()
    else:
        # get similar layers
        l1 = 0
        layers = []
        for zz in range(zgrid.shape[0] - 1):
            if (resmodel[:, :, zz] == resmodel[:, :, zz + 1]).all() == False:
                layers.append((l1, zz))
                l1 = zz + 1
        # need to add on the bottom layers
        layers.append((l1, zgrid.shape[0] - 1))

        # write out the layers from resmodel
        for ll in layers:
            ifid.write("{0} {1}\n".format(ll[0] + 1, ll[1] + 1))
            for xx in range(xgrid.shape[0]):
                for yy in range(ygrid.shape[0]):
                    ifid.write("{0:.0f} ".format(resmodel[xx, yy, ll[0]]))
                ifid.write("\n")

    print "Wrote file to: " + ifn
    return ifn


def readInit3D(initfn):
    """
    read an initial file and return the pertinent information including grid
    positions in coordinates relative to the center point (0,0) and
    starting model.

    Arguments:
    ----------

        **initfn** : full path to initializing file.

    Returns:
    --------

        **xgrid** : np.array(nx)
                    array of nodes in S --> N direction

        **ygrid** : np.array(ny)
                    array of nodes in the W --> E direction

        **zgrid** : np.array(nz)
                    array of nodes in vertical direction positive downwards

        **resistivitivityModel** : dictionary
                    dictionary of the starting model with keys as layers

        **reslst** : list
                    list of resistivity values in the model

        **titlestr** : string
                       title string

    """

    ifid = file(initfn, "r")
    ilines = ifid.readlines()
    ifid.close()

    titlestr = ilines[0]

    # get size of dimensions, remembering that x is N-S, y is E-W, z is + down
    nsize = ilines[1].strip().split()
    nx = int(nsize[0])
    ny = int(nsize[1])
    nz = int(nsize[2])

    # initialize empy arrays to put things into
    xnodes = np.zeros(nx)
    ynodes = np.zeros(ny)
    znodes = np.zeros(nz)
    resmodel = np.zeros((nx, ny, nz))

    # get the grid line locations
    nn = 2
    xx = 0
    while xx < nx:
        iline = ilines[nn].strip().split()
        for xg in iline:
            xnodes[xx] = float(xg)
            xx += 1
        nn += 1

    yy = 0
    while yy < ny:
        iline = ilines[nn].strip().split()
        for yg in iline:
            ynodes[yy] = float(yg)
            yy += 1
        nn += 1

    zz = 0
    while zz < nz:
        iline = ilines[nn].strip().split()
        for zg in iline:
            znodes[zz] = float(zg)
            zz += 1
        nn += 1

    # put the grids into coordinates relative to the center of the grid
    xgrid = xnodes.copy()
    xgrid[: int(nx / 2)] = -np.array(
        [xnodes[ii : int(nx / 2)].sum() for ii in range(int(nx / 2))]
    )
    xgrid[int(nx / 2) :] = (
        np.array([xnodes[int(nx / 2) : ii + 1].sum() for ii in range(int(nx / 2), nx)])
        - xnodes[int(nx / 2)]
    )

    ygrid = ynodes.copy()
    ygrid[: int(ny / 2)] = -np.array(
        [ynodes[ii : int(ny / 2)].sum() for ii in range(int(ny / 2))]
    )
    ygrid[int(ny / 2) :] = (
        np.array([ynodes[int(ny / 2) : ii + 1].sum() for ii in range(int(ny / 2), ny)])
        - ynodes[int(ny / 2)]
    )

    zgrid = np.array([znodes[: ii + 1].sum() for ii in range(nz)])

    # get the resistivity values
    reslst = [float(rr) for rr in ilines[nn].strip().split()]
    nn += 1

    # get model
    iline = ilines[nn].strip().split()
    if len(iline) == 0 or len(iline) == 1:
        return xgrid, ygrid, zgrid, reslst, titlestr, resmodel
    else:
        while nn < len(ilines):

            iline = ilines[nn].strip().split()
            if len(iline) == 2:
                l1 = int(iline[0]) - 1
                l2 = int(iline[1])
                nn += 1
                xx = 0
            elif len(iline) == 0:
                break
            else:
                yy = 0
                while yy < ny:
                    resmodel[xx, yy, l1:l2] = int(iline[yy])
                    #                        if l1==20:
                    #                            print nn,xx,yy,l1,l2,iline[yy]
                    yy += 1
                xx += 1
                nn += 1

        return xgrid, ygrid, zgrid, reslst, titlestr, resmodel, xnodes, ynodes, znodes


def writeStartupFile(
    datafn,
    initialfn=None,
    outputfn=None,
    savepath=None,
    apriorfn=None,
    modells=[5, 0.3, 0.3, 0.3],
    targetrms=1.0,
    control=None,
    maxiter=10,
    errortol=None,
    staticfn=None,
    lagrange=None,
):
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

    # create the output filename
    if savepath is None:
        sfile = os.path.join(os.path.dirname(datafn), "startup")
    elif savepath.find(".") == -1:
        sfile = os.path.join(savepath, "startup")
    else:
        sfile = savepath

    sfid = file(sfile, "w")

    sfid.write("DATA_FILE" + " " * 11 + "../" + os.path.basename(datafn) + "\n")

    if outputfn is None:
        sfid.write("OUTPUT_FILE" + " " * 9 + "Iter_ \n")
    else:
        sfid.write("OUTPUT_FILE" + " " * 9 + outputfn + " \n")

    if initialfn is None:
        sfid.write("INITIAL_MODEL_FILE" + " " * 2 + "../init3d \n")
    else:
        sfid.write("INITIAL_MODEL_FILE" + " " * 2 + initialfn + " \n")

    if apriorfn is None:
        sfid.write("PRIOR_MODEL_FILE" + " " * 4 + "default \n")
    else:
        sfid.write("PRIOR_MODEL_FILE" + " " * 4 + apriorfn + " \n")

    if control is None:
        sfid.write("CONTROL_MODEL_INDEX" + " " + "default \n")
    else:
        sfid.write("CONTROL_MODEL_INDEX" + " " + control + " \n")

    sfid.write("TARGET_RMS" + " " * 10 + "{0} \n".format(targetrms))

    sfid.write("MAX_NO_ITERATION" + " " * 4 + "{0} \n".format(maxiter))

    sfid.write(
        "MODEL_LENGTH_SCALE"
        + " " * 2
        + "{0} {1:.1f} {1:.1f} {1:.1f} \n".format(
            modells[0], modells[1], modells[2], modells[3]
        )
    )

    if lagrange is None:
        sfid.write("LAGRANGE_INFO" + " " * 7 + "default \n")
    else:
        sfid.write("LAGRANGE_INFO" + " " * 7 + lagrange + " \n")

    if errortol is None:
        sfid.write("ERROR_TOL_LEVEL" + " " * 5 + "default \n")
    else:
        sfid.write("ERROR_TOL_LEVEL" + " " * 5 + errortol + " \n")

    if staticfn is None:
        sfid.write("STATIC_FILE" + " " * 9 + "default \n")
    else:
        sfid.write("STATIC_FILE" + " " * 9 + staticfn + " \n")

    sfid.close()

    print "Wrote startup file to: " + sfile

    return sfile


def readDataFile(datafn, sitesfn=None, units="mv"):
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

    if units == "mv":
        zconv = 796.0
    else:
        zconv = 1

    dfid = file(datafn, "r")
    dlines = dfid.readlines()

    # get size number of stations, number of frequencies, number of Z
    # components
    ns, nf, nz = np.array(dlines[0].strip().split(), dtype="int")
    nsstart = 2

    findlst = []
    for ii, dline in enumerate(dlines[1:50], 1):
        if dline.find("Station_Location: N-S") == 0:
            findlst.append(ii)
        elif dline.find("Station_Location: E-W") == 0:
            findlst.append(ii)
        elif dline.find("DATA_Period:") == 0:
            findlst.append(ii)

    ncol = len(dlines[nsstart].strip().split())
    #    print ncol
    #    nsstop=nsstart+ns/ncol+1
    #    ewstart=nsstop+1
    #    ewstop=ewstart+ns/ncol+1
    #    zstart=ewstop
    #    print nsstop,ewstart,ewstop,zstart

    # get site names if entered a sites file
    if sitesfn is not None:
        slst, sitelst = readSitesFile(sitesfn)
    else:
        sitelst = np.arange(ns)

    # get N-S locations
    nsarr = np.zeros(ns)
    for ii, dline in enumerate(dlines[findlst[0] + 1 : findlst[1]], 0):
        dline = dline.strip().split()
        for jj in range(ncol):
            try:
                nsarr[ii * ncol + jj] = float(dline[jj])
            except IndexError:
                pass
            except ValueError:
                break

    # get E-W locations
    ewarr = np.zeros(ns)
    for ii, dline in enumerate(dlines[findlst[1] + 1 : findlst[2]], 0):
        dline = dline.strip().split()
        for jj in range(8):
            try:
                ewarr[ii * ncol + jj] = float(dline[jj])
            except IndexError:
                pass
            except ValueError:
                break
    # make some empty array to put stuff into
    period = np.zeros(nf)
    zarr = np.zeros((ns, nf, 2, 2), dtype=np.complex)
    zerr = np.zeros_like(zarr)
    zerrmap = np.zeros_like(zarr)

    # get data
    pcount = 0
    zcount = 0
    for ii, dl in enumerate(dlines[findlst[2] : findlst[2] + nf * (ns + 1)]):
        if dl.find("DATA_Period") == 0:
            period[pcount] = float(dl.strip().split()[1])
            kk = 0
            pcount += 1
            if ii == 0:
                pass
            else:
                zcount += 1
        else:
            zline = np.array(dl.strip().split(), dtype=np.float) * zconv
            zarr[kk, zcount, :, :] = np.array(
                [
                    [zline[0] - 1j * zline[1], zline[2] - 1j * zline[3]],
                    [zline[4] - 1j * zline[5], zline[6] - 1j * zline[7]],
                ]
            )
            kk += 1

    # if the data file is made from this program or is the input data file than
    # get the errors from that file
    if len(dlines) > 2 * nf * ns:
        print "Getting Error"
        pecount = 0
        zecount = 0
        for ii, dl in enumerate(
            dlines[findlst[2] + nf * (ns + 1) : findlst[2] + 2 * nf * (ns + 1)]
        ):
            if dl.find("ERROR_Period") == 0:
                kk = 0
                pecount += 1
                if ii == 0:
                    pass
                else:
                    zecount += 1
            else:
                zline = np.array(dl.strip().split(), dtype=np.float) * zconv
                zerr[kk, zecount, :, :] = np.array(
                    [
                        [zline[0] - 1j * zline[1], zline[2] - 1j * zline[3]],
                        [zline[4] - 1j * zline[5], zline[6] - 1j * zline[7]],
                    ]
                )
                kk += 1

    # get errormap values
    if len(dlines) > 3 * nf * ns:
        print "Getting Error Map"
        pmcount = 0
        zmcount = 0
        for ii, dl in enumerate(
            dlines[findlst[2] + 2 * nf * (ns + 1) : findlst[2] + 3 * nf * (ns + 1)]
        ):
            if dl.find("ERMAP_Period") == 0:
                kk = 0
                pmcount += 1
                if ii == 0:
                    pass
                else:
                    zmcount += 1
            else:
                # account for end of file empty lines
                if len(dl.split()) > 2:
                    zline = np.array(dl.strip().split(), dtype=np.float)
                    zerrmap[kk, zmcount, :, :] = np.array(
                        [
                            [zline[0] - 1j * zline[1], zline[2] - 1j * zline[3]],
                            [zline[4] - 1j * zline[5], zline[6] - 1j * zline[7]],
                        ]
                    )
                    kk += 1

    # multiply errmap and error and convert from Ohm to mv/km nT
    zerr = zerr * zerrmap

    return period, zarr, zerr, nsarr, ewarr, sitelst


def plotDataResPhase(
    datafn,
    respfn=None,
    sitesfn=None,
    plottype="1",
    plotnum=1,
    dpi=150,
    units="mv",
    colormode="color",
):
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

    # plot in color mode or black and white
    if colormode == "color":
        # color for data
        cted = (0, 0, 1)
        ctmd = (1, 0, 0)
        mted = "*"
        mtmd = "*"

        # color for occam model
        ctem = (0, 0.3, 1.0)
        ctmm = (1, 0.3, 0)
        mtem = "+"
        mtmm = "+"

    elif colormode == "bw":
        # color for data
        cted = (0, 0, 0)
        ctmd = (0, 0, 0)
        mted = "*"
        mtmd = "v"

        # color for occam model
        ctem = (0.6, 0.6, 0.6)
        ctmm = (0.6, 0.6, 0.6)
        mtem = "+"
        mtmm = "x"

    # load the data file
    period, dz, dzerr, north, east, slst = readDataFile(
        datafn, sitesfn=sitesfn, units=units
    )
    # get shape of impedance tensors
    ns, nf = dz.shape[0], dz.shape[1]

    # read in response files
    if respfn is not None:
        rzlst = []
        rzerrlst = []
        if not isinstance(respfn, list):
            respfn = [respfn]
        for rfile in respfn:
            period, rz, rzerr, north, east, slst = readDataFile(
                rfile, sitesfn=sitesfn, units=units
            )
            rzlst.append(rz)
            rzerrlst.append(rzerr)
    else:
        rzlst = []
    # get number of response files
    nr = len(rzlst)

    if isinstance(plottype, list):
        ns = len(plottype)

    plt.rcParams["font.size"] = 10
    plt.rcParams["figure.subplot.left"] = 0.13
    plt.rcParams["figure.subplot.right"] = 0.98
    plt.rcParams["figure.subplot.bottom"] = 0.1
    plt.rcParams["figure.subplot.top"] = 0.92
    plt.rcParams["figure.subplot.wspace"] = 0.25
    plt.rcParams["figure.subplot.hspace"] = 0.05

    fontdict = {"size": 12, "weight": "bold"}
    gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1.5], hspace=0.1)

    if plottype != "1":
        pstationlst = []
        if not isinstance(plottype, list):
            plottype = [plottype]
        for ii, station in enumerate(slst):
            if isinstance(station, str):
                for pstation in plottype:
                    if station.find(str(pstation)) >= 0:
                        pstationlst.append(ii)
            else:
                for pstation in plottype:
                    if station == int(pstation):
                        pstationlst.append(ii)
    else:
        pstationlst = np.arange(ns)

    for jj in pstationlst:
        print "Plotting: " + str(slst[jj])

        # check for masked points
        dz[jj][np.where(dz[jj] == 7.95204e5 - 7.95204e5j)] = 0.0 + 0.0j
        dzerr[jj][np.where(dz[jj] == 7.95204e5 - 7.95204e5j)] = 1.0 + 1.0j

        # convert to apparent resistivity and phase
        rp = Z.ResPhase(dz[jj], period, zvar=dzerr[jj])

        # find locations where points have been masked
        nzxx = np.where(rp.resxx != 0)[0]
        nzxy = np.where(rp.resxy != 0)[0]
        nzyx = np.where(rp.resyx != 0)[0]
        nzyy = np.where(rp.resyy != 0)[0]

        if respfn is not None:
            plotr = True
        else:
            plotr = False

        # make figure for xy,yx components
        if plotnum == 1:
            fig = plt.figure(jj, [10, 12], dpi=dpi)
            gs.update(hspace=0.1, wspace=0.15, left=0.1)
        elif plotnum == 2:
            fig = plt.figure(jj, [12, 12], dpi=dpi)
            gs.update(hspace=0.1, wspace=0.15, left=0.07)

        # ---------plot the apparent resistivity--------------------------------
        if plotnum == 1:
            ax = fig.add_subplot(gs[0, :])
            ax2 = fig.add_subplot(gs[1, :], sharex=ax)
            ax.yaxis.set_label_coords(-0.055, 0.5)
            ax2.yaxis.set_label_coords(-0.055, 0.5)
        elif plotnum == 2:
            ax = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0], sharex=ax)
            ax.yaxis.set_label_coords(-0.075, 0.5)
            ax2.yaxis.set_label_coords(-0.075, 0.5)

        fig.suptitle(str(slst[jj]), fontdict={"size": 15, "weight": "bold"})
        erxy = ax.errorbar(
            period[nzxy],
            rp.resxy[nzxy],
            marker=mted,
            ms=4,
            mfc="None",
            mec=cted,
            mew=1,
            ls=":",
            yerr=rp.resxyerr[nzxy],
            ecolor=cted,
            color=cted,
        )
        eryx = ax.errorbar(
            period[nzyx],
            rp.resyx[nzyx],
            marker=mtmd,
            ms=4,
            mfc="None",
            mec=ctmd,
            mew=1,
            ls=":",
            yerr=rp.resyxerr[nzyx],
            ecolor=ctmd,
            color=ctmd,
        )
        if plotr == True:
            for rr in range(nr):
                if colormode == "color":
                    cxy = (0, 0.4 + float(rr) / (3 * nr), 0)
                    cyx = (
                        0.7 + float(rr) / (4 * nr),
                        0.13,
                        0.63 - float(rr) / (4 * nr),
                    )
                elif colormode == "bw":
                    cxy = (
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                    )
                    cyx = (
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                    )

                rpr = Z.ResPhase(rzlst[rr][jj], period, zvar=rzerrlst[rr][jj])

                #                rms=np.sqrt(np.sum([abs(np.linalg.det(rp.z[ll])-
                #                                        np.linalg.det(rpr.z[ll]))**2
                #                            for ll in range(len(rp.period))])/len(rp.period))
                rms = np.sqrt(
                    np.mean(
                        [
                            (
                                np.sqrt(abs(np.linalg.det(rp.z[ll])))
                                - np.sqrt(abs(np.linalg.det(rpr.z[ll])))
                            )
                            ** 2
                            for ll in range(len(rp.period))
                        ]
                    )
                )
                print "RMS = {:.2f}".format(rms)
                erxyr = ax.errorbar(
                    period[nzxy],
                    rpr.resxy[nzxy],
                    marker=mtem,
                    ms=8,
                    mfc="None",
                    mec=cxy,
                    mew=1,
                    ls="--",
                    yerr=rpr.resxyerr[nzxy],
                    ecolor=cxy,
                    color=cxy,
                )
                eryxr = ax.errorbar(
                    period[nzyx],
                    rpr.resyx[nzyx],
                    marker=mtmm,
                    ms=8,
                    mfc="None",
                    mec=cyx,
                    mew=1,
                    ls="--",
                    yerr=rpr.resyxerr[nzyx],
                    ecolor=cyx,
                    color=cyx,
                )
        # ax.set_xlabel('Period (s)',fontdict=fontdict)
        pylab.setp(ax.get_xticklabels(), visible=False)
        ax.set_ylabel("App. Res. ($\mathbf{\Omega \cdot m}$)", fontdict=fontdict)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(
            xmin=10 ** (np.floor(np.log10(period[0]))),
            xmax=10 ** (np.ceil(np.log10(period[-1]))),
        )
        ax.grid(True, alpha=0.25)
        if plotr == True:
            ax.legend(
                (erxy[0], eryx[0], erxyr[0], eryxr[0]),
                ("Data $E_x/B_y$", "Data $E_y/B_x$", "Mod $E_x/B_y$", "Mod $E_y/B_x$"),
                loc=0,
                markerscale=1,
                borderaxespad=0.01,
                labelspacing=0.07,
                handletextpad=0.2,
                borderpad=0.02,
            )
        else:
            ax.legend(
                (erxy[0], eryx[0]),
                ("$E_x/B_y$", "$E_y/B_x$"),
                loc=0,
                markerscale=1,
                borderaxespad=0.01,
                labelspacing=0.07,
                handletextpad=0.2,
                borderpad=0.02,
            )

        # -----Plot the phase---------------------------------------------------

        ax2.errorbar(
            period[nzxy],
            rp.phasexy[nzxy],
            marker=mted,
            ms=4,
            mfc="None",
            mec=cted,
            mew=1,
            ls=":",
            yerr=rp.phasexyerr[nzxy],
            ecolor=cted,
            color=cted,
        )
        ax2.errorbar(
            period[nzyx],
            np.array(rp.phaseyx[nzyx]) + 180,
            marker=mtmd,
            ms=4,
            mfc="None",
            mec=ctmd,
            mew=1,
            ls=":",
            yerr=rp.phaseyxerr[nzyx],
            ecolor=ctmd,
            color=ctmd,
        )
        if plotr == True:
            for rr in range(nr):
                if colormode == "color":
                    cxy = (0, 0.4 + float(rr) / (3 * nr), 0)
                    cyx = (
                        0.7 + float(rr) / (4 * nr),
                        0.13,
                        0.63 - float(rr) / (4 * nr),
                    )
                elif colormode == "bw":
                    cxy = (
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                    )
                    cyx = (
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                        1 - 1.25 / (rr + 2.0),
                    )
                rpr = Z.ResPhase(rzlst[rr][jj], period, zvar=rzerrlst[rr][jj])
                ax2.errorbar(
                    period[nzxy],
                    rpr.phasexy[nzxy],
                    marker=mtem,
                    ms=8,
                    mfc="None",
                    mec=cxy,
                    mew=1,
                    ls="--",
                    yerr=rp.phasexyerr[nzxy],
                    ecolor=cxy,
                    color=cxy,
                )
                ax2.errorbar(
                    period[nzyx],
                    np.array(rpr.phaseyx[nzyx]) + 180,
                    marker=mtmm,
                    ms=8,
                    mfc="None",
                    mec=cyx,
                    mew=1,
                    ls="--",
                    yerr=rp.phaseyxerr[nzyx],
                    ecolor=cyx,
                    color=cyx,
                )
        ax2.set_xlabel("Period (s)", fontdict)
        ax2.set_ylabel("Phase (deg)", fontdict)
        ax2.set_xscale("log")
        # ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
        #         xmax=10**(np.ceil(np.log10(period[-1]))))
        # check the phase to see if any point are outside of [0:90]
        if min(rp.phasexy) < 0 or min(rp.phaseyx + 180) < 0:
            pymin = min([min(rp.phasexy), min(rp.phaseyx + 180)])
            if pymin > 0:
                pymin = 0
        else:
            pymin = 0

        if max(rp.phasexy) > 90 or max(rp.phaseyx + 180) > 90:
            pymax = min([max(rp.phasexy), max(rp.phaseyx + 180)])
            if pymax < 91:
                pymax = 90
        else:
            pymax = 90

        ax2.set_ylim(ymin=pymin, ymax=pymax)
        ax2.yaxis.set_major_locator(MultipleLocator(30))
        ax2.yaxis.set_minor_locator(MultipleLocator(1))
        ax2.grid(True, alpha=0.25)

        if plotnum == 2:
            # ---------plot the apparent resistivity----------------------------
            ax3 = plt.subplot(gs[0, 1])
            ax3.yaxis.set_label_coords(-0.1, 0.5)
            erxx = ax3.errorbar(
                period[nzxx],
                rp.resxx[nzxx],
                marker=mted,
                ms=4,
                mfc="None",
                mec=cted,
                mew=1,
                ls=":",
                yerr=rp.resxxerr[nzxx],
                ecolor=cted,
                color=cted,
            )
            eryy = ax3.errorbar(
                period[nzyy],
                rp.resyy[nzyy],
                marker=mtmd,
                ms=4,
                mfc="None",
                mec=ctmd,
                mew=1,
                ls=":",
                yerr=rp.resyyerr[nzyy],
                ecolor=ctmd,
                color=ctmd,
            )
            if plotr == True:
                for rr in range(nr):
                    if colormode == "color":
                        cxy = (0, 0.4 + float(rr) / (3 * nr), 0)
                        cyx = (
                            0.7 + float(rr) / (4 * nr),
                            0.13,
                            0.63 - float(rr) / (4 * nr),
                        )
                    elif colormode == "bw":
                        cxy = (
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                        )
                        cyx = (
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                        )
                    rpr = Z.ResPhase(rzlst[rr][jj], period, zvar=rzerrlst[rr][jj])
                    erxxr = ax3.errorbar(
                        period[nzxx],
                        rpr.resxx[nzxx],
                        marker=mtem,
                        ms=8,
                        mfc="None",
                        mec=cxy,
                        mew=1,
                        ls="--",
                        yerr=rpr.resxxerr[nzxx],
                        ecolor=cxy,
                        color=cxy,
                    )
                    eryyr = ax3.errorbar(
                        period[nzyy],
                        rpr.resyy[nzyy],
                        marker=mtmm,
                        ms=8,
                        mfc="None",
                        mec=cyx,
                        mew=1,
                        ls="--",
                        yerr=rpr.resyyerr[nzyy],
                        ecolor=cyx,
                        color=cyx,
                    )

            ax3.set_yscale("log")
            ax3.set_xscale("log")
            pylab.setp(ax3.get_xticklabels(), visible=False)
            ax3.set_xlim(
                xmin=10 ** (np.floor(np.log10(period[0]))),
                xmax=10 ** (np.ceil(np.log10(period[-1]))),
            )
            ax3.grid(True, alpha=0.25)
            if plotr == True:
                ax3.legend(
                    (erxx[0], eryy[0], erxxr[0], eryyr[0]),
                    (
                        "Data $E_x/B_x$",
                        "Data $E_y/B_y$",
                        "Mod $E_x/B_x$",
                        "Mod $E_y/B_y$",
                    ),
                    loc=0,
                    markerscale=1,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.02,
                )
            else:
                ax3.legend(
                    (erxx[0], eryy[0]),
                    ("$E_x/B_x$", "$E_y/B_y$"),
                    loc=0,
                    markerscale=1,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.02,
                )

            # -----Plot the phase-----------------------------------------------
            ax4 = plt.subplot(gs[1, 1], sharex=ax3)

            ax4.yaxis.set_label_coords(-0.1, 0.5)
            ax4.errorbar(
                period[nzxx],
                rp.phasexx[nzxx],
                marker=mted,
                ms=4,
                mfc="None",
                mec=cted,
                mew=1,
                ls=":",
                yerr=rp.phasexxerr[nzxx],
                ecolor=cted,
                color=cted,
            )
            ax4.errorbar(
                period[nzyy],
                np.array(rp.phaseyy[nzyy]),
                marker=mtmd,
                ms=4,
                mfc="None",
                mec=ctmd,
                mew=1,
                ls=":",
                yerr=rp.phaseyyerr[nzyy],
                ecolor=ctmd,
                color=ctmd,
            )
            if plotr == True:
                for rr in range(nr):
                    if colormode == "color":
                        cxy = (0, 0.4 + float(rr) / (3 * nr), 0)
                        cyx = (
                            0.7 + float(rr) / (4 * nr),
                            0.13,
                            0.63 - float(rr) / (4 * nr),
                        )
                    elif colormode == "bw":
                        cxy = (
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                        )
                        cyx = (
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                        )
                    rpr = Z.ResPhase(rzlst[rr][jj], period, zvar=rzerrlst[rr][jj])
                    ax4.errorbar(
                        period[nzxx],
                        rpr.phasexx[nzxx],
                        marker=mtem,
                        ms=8,
                        mfc="None",
                        mec=cxy,
                        mew=1,
                        ls="--",
                        yerr=rp.phasexxerr[nzxx],
                        ecolor=cxy,
                        color=cxy,
                    )
                    ax4.errorbar(
                        period[nzyy],
                        np.array(rpr.phaseyy[nzyy]),
                        marker=mtmm,
                        ms=8,
                        mfc="None",
                        mec=cyx,
                        mew=1,
                        ls="--",
                        yerr=rp.phaseyyerr[nzyy],
                        ecolor=cyx,
                        color=cyx,
                    )
            ax4.set_xlabel("Period (s)", fontdict)
            # ax4.set_ylabel('Imepdance Phase (deg)',fontdict)
            ax4.set_xscale("log")
            # ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
            #         xmax=10**(np.ceil(np.log10(period[-1]))))
            ax4.set_ylim(ymin=-180, ymax=180)
            ax4.yaxis.set_major_locator(MultipleLocator(30))
            ax4.yaxis.set_minor_locator(MultipleLocator(5))
            ax4.grid(True, alpha=0.25)


def plotTensorMaps(
    datafn,
    respfn=None,
    sitesfn=None,
    periodlst=None,
    esize=(1, 1, 5, 5),
    ecolor="phimin",
    colormm=[(0, 90), (0, 1), (0, 4), (-2, 2)],
    xpad=0.500,
    units="mv",
    dpi=150,
):
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

    period, zd, zderr, nsarr, ewarr, sitelst = readDataFile(
        datafn, sitesfn=sitesfn, units=units
    )

    if respfn is not None:
        period, zr, zrerr, nsarr, ewarr, sitelst = readDataFile(
            respfn, sitesfn=sitesfn, units=units
        )

    if periodlst is None:
        periodlst = range(len(period))

    # put locations into an logical coordinate system in km
    nsarr = -nsarr / 1000
    ewarr = -ewarr / 1000

    # get coloring min's and max's
    if colormm is not None:
        ptmin, ptmax = (colormm[0][0] * np.pi / 180, colormm[0][1] * np.pi / 180)
        ptrmin, ptrmax = colormm[1]
        rtmin, rtmax = colormm[2]
        rtrmin, rtrmax = colormm[3]
    else:
        pass

    # get ellipse sizes
    ptsize = esize[0]
    ptrsize = esize[1]
    rtsize = esize[2]
    rtrsize = esize[3]

    plt.rcParams["font.size"] = 10
    plt.rcParams["figure.subplot.left"] = 0.03
    plt.rcParams["figure.subplot.right"] = 0.98
    plt.rcParams["figure.subplot.bottom"] = 0.1
    plt.rcParams["figure.subplot.top"] = 0.90
    plt.rcParams["figure.subplot.wspace"] = 0.005
    plt.rcParams["figure.subplot.hspace"] = 0.005

    ns = zd.shape[0]

    for ff, per in enumerate(periodlst):
        print "Plotting Period: {0:.5g}".format(period[per])
        fig = plt.figure(per + 1, dpi=dpi)

        # get phase tensor
        pt = Z.PhaseTensor(zd[:, per])

        # get resistivity tensor
        rt = Z.ResistivityTensor(zd[:, per], np.repeat(1.0 / period[per], ns))

        if respfn is not None:
            # get phase tensor and residual phase tensor
            ptr = Z.PhaseTensor(zr[:, per])
            ptd = Z.PhaseTensorResidual(zd[:, per], zr[:, per])

            # get resistivity tensor and residual
            rtr = Z.ResistivityTensor(zr[:, per], np.repeat(1.0 / period[per], ns))
            rtd = Z.ResistivityTensorResidual(
                zd[:, per], zr[:, per], np.repeat(1.0 / period[per], ns)
            )

            if colormm is None:
                if ecolor == "phimin":
                    ptmin, ptmax = (
                        ptr.phimin.min() / (np.pi / 2),
                        ptr.phimin.max() / (np.pi / 2),
                    )
                elif ecolor == "beta":
                    ptmin, ptmax = (ptr.beta.min(), ptr.beta.max())

                ptrmin, ptrmax = (ptd.ecolor.min(), ptd.ecolor.max())
                rtmin, rtmax = (np.log10(rtr.rhodet.min()), np.log10(rtr.rhodet.max()))
                rtrmin, rtrmax = rtd.rhodet.min(), rtd.rhodet.max()
            # make subplots
            ax1 = fig.add_subplot(2, 3, 1, aspect="equal")
            ax2 = fig.add_subplot(2, 3, 2, aspect="equal")
            ax3 = fig.add_subplot(2, 3, 3, aspect="equal")
            ax4 = fig.add_subplot(2, 3, 4, aspect="equal")
            ax5 = fig.add_subplot(2, 3, 5, aspect="equal")
            ax6 = fig.add_subplot(2, 3, 6, aspect="equal")

            for jj in range(ns):
                # -----------plot data phase tensors---------------
                eheightd = pt.phimin[jj] / ptr.phimax.max() * ptsize
                ewidthd = pt.phimax[jj] / ptr.phimax.max() * ptsize

                ellipd = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=ewidthd,
                    height=eheightd,
                    angle=pt.azimuth[jj],
                )
                # color ellipse:
                if ecolor == "phimin":
                    cvar = (pt.phimin[jj] / (np.pi / 2) - ptmin) / (ptmax - ptmin)
                    if abs(cvar) > 1:
                        ellipd.set_facecolor((1, 0, 0.1))
                    else:
                        ellipd.set_facecolor((1, 1 - abs(cvar), 0.1))
                if ecolor == "beta":
                    cvar = (abs(pt.beta[jj]) - ptmin) / (ptmax - ptmin)
                    if abs(cvar) > 1:
                        ellipd.set_facecolor((1, 1, 0.1))
                    else:
                        ellipd.set_facecolor((1 - cvars, 1 - cvars, 1))

                ax1.add_artist(ellipd)

                # ----------plot response phase tensors---------------------
                eheightr = ptr.phimin[jj] / ptr.phimax.max() * ptsize
                ewidthr = ptr.phimax[jj] / ptr.phimax.max() * ptsize

                ellipr = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=ewidthr,
                    height=eheightr,
                    angle=ptr.azimuth[jj],
                )
                # color ellipse:
                if ecolor == "phimin":
                    cvar = (ptr.phimin[jj] / (np.pi / 2) - ptmin) / (ptmax - ptmin)
                    if abs(cvar) > 1:
                        ellipr.set_facecolor((1, 0, 0.1))
                    else:
                        ellipr.set_facecolor((1, 1 - abs(cvar), 0.1))
                if ecolor == "beta":
                    cvar = (abs(ptr.beta[jj]) - ptmin) / (ptmax - ptmin)
                    if abs(cvar) > 1:
                        ellipr.set_facecolor((1, 1, 0.1))
                    else:
                        ellipr.set_facecolor((1 - cvars, 1 - cvars, 1))
                ax2.add_artist(ellipr)

                # --------plot residual phase tensors-------------
                eheight = ptd.phimin[jj] / ptd.phimax.max() * ptrsize
                ewidth = ptd.phimax[jj] / ptd.phimax.max() * ptrsize

                ellip = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=ewidth,
                    height=eheight,
                    angle=ptd.azimuth[jj] - 90,
                )
                # color ellipse:
                cvar = (ptd.ecolor[jj] - ptrmin) / (ptrmax - ptrmin)
                if abs(cvar) > 1:
                    ellip.set_facecolor((0, 0, 0))
                else:
                    ellip.set_facecolor((abs(cvar), 0.5, 0.5))

                ax3.add_artist(ellip)

                # -----------plot data resistivity tensors---------------
                rheightd = rt.rhomin[jj] / rtr.rhomax.max() * rtsize
                rwidthd = rt.rhomax[jj] / rtr.rhomax.max() * rtsize

                rellipd = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=rwidthd,
                    height=rheightd,
                    angle=rt.rhoazimuth[jj],
                )
                # color ellipse:
                cvar = (np.log10(rt.rhodet[jj]) - rtmin) / (rtmax - rtmin)
                if cvar > 0.5:
                    if cvar > 1:
                        rellipd.set_facecolor((0, 0, 1))
                    else:
                        rellipd.set_facecolor((1 - abs(cvar), 1 - abs(cvar), 1))
                else:
                    if cvar < -1:
                        rellipd.set_facecolor((1, 0, 0))
                    else:
                        rellipd.set_facecolor((1, 1 - abs(cvar), 1 - abs(cvar)))

                ax4.add_artist(rellipd)

                # ----------plot response resistivity tensors-------------------
                rheightr = rtr.rhomin[jj] / rtr.rhomax.max() * rtsize
                rwidthr = rtr.rhomax[jj] / rtr.rhomax.max() * rtsize

                rellipr = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=rwidthr,
                    height=rheightr,
                    angle=rtr.rhoazimuth[jj],
                )

                # color ellipse:
                cvar = (np.log10(rtr.rhodet[jj]) - rtmin) / (rtmax - rtmin)
                if cvar > 0.5:
                    if cvar > 1:
                        rellipr.set_facecolor((0, 0, 1))
                    else:
                        rellipr.set_facecolor((1 - abs(cvar), 1 - abs(cvar), 1))
                else:
                    if cvar < -1:
                        rellipr.set_facecolor((1, 0, 0))
                    else:
                        rellipr.set_facecolor((1, 1 - abs(cvar), 1 - abs(cvar)))

                ax5.add_artist(rellipr)

                # --------plot residual resistivity tensors-------------
                rheight = rtd.rhomin[jj] / rtd.rhomax.max() * rtrsize
                rwidth = rtd.rhomax[jj] / rtd.rhomax.max() * rtrsize

                rellip = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=rwidth,
                    height=rheight,
                    angle=rtd.azimuth[jj] - 90,
                )
                # color ellipse:
                cvar = (rtd.rhodet[jj] - rtrmin) / (rtrmax - rtrmin)
                if cvar < 0:
                    if cvar < -1:
                        rellip.set_facecolor((0, 0, 1))
                    else:
                        rellip.set_facecolor((1 - abs(cvar), 1 - abs(cvar), 1))
                else:
                    if cvar > 1:
                        rellip.set_facecolor((1, 0, 0))
                    else:
                        rellip.set_facecolor((1, 1 - abs(cvar), 1 - abs(cvar)))

                ax6.add_artist(rellip)

            for aa, ax in enumerate([ax1, ax2, ax3, ax4, ax5, ax6]):
                ax.set_xlim(ewarr.min() - xpad, ewarr.max() + xpad)
                ax.set_ylim(nsarr.min() - xpad, nsarr.max() + xpad)
                ax.grid("on")
                if aa < 3:
                    pylab.setp(ax.get_xticklabels(), visible=False)
                if aa == 0 or aa == 3:
                    pass
                else:
                    pylab.setp(ax.get_yticklabels(), visible=False)

                cbax = mcb.make_axes(ax, shrink=0.9, pad=0.05, orientation="vertical")
                if aa == 0 or aa == 1:
                    cbx = mcb.ColorbarBase(
                        cbax[0],
                        cmap=ptcmap,
                        norm=Normalize(
                            vmin=ptmin * 180 / np.pi, vmax=ptmax * 180 / np.pi
                        ),
                        orientation="vertical",
                        format="%.2g",
                    )

                    cbx.set_label("Phase (deg)", fontdict={"size": 7, "weight": "bold"})
                if aa == 2:
                    cbx = mcb.ColorbarBase(
                        cbax[0],
                        cmap=ptcmap2,
                        norm=Normalize(vmin=ptrmin, vmax=ptrmax),
                        orientation="vertical",
                        format="%.2g",
                    )

                    cbx.set_label(
                        "$\Delta_{\Phi}$", fontdict={"size": 7, "weight": "bold"}
                    )
                if aa == 3 or aa == 4:
                    cbx = mcb.ColorbarBase(
                        cbax[0],
                        cmap=rtcmapr,
                        norm=Normalize(vmin=10 ** rtmin, vmax=10 ** rtmax),
                        orientation="vertical",
                        format="%.2g",
                    )

                    cbx.set_label(
                        "App. Res. ($\Omega \cdot$m)",
                        fontdict={"size": 7, "weight": "bold"},
                    )
                if aa == 5:
                    cbx = mcb.ColorbarBase(
                        cbax[0],
                        cmap=rtcmap,
                        norm=Normalize(vmin=rtrmin, vmax=rtrmax),
                        orientation="vertical",
                        format="%.2g",
                    )

                    cbx.set_label(
                        "$\Delta_{rho}$", fontdict={"size": 7, "weight": "bold"}
                    )

            plt.show()

        # ----Plot Just the data------------------
        else:
            if colormm is None:
                if ecolor == "phimin":
                    ptmin, ptmax = (
                        pt.phimin.min() / (np.pi / 2),
                        pt.phimin.max() / (np.pi / 2),
                    )
                elif ecolor == "beta":
                    ptmin, ptmax = (pt.beta.min(), pt.beta.max())

                rtmin, rtmax = (np.log10(rt.rhodet.min()), np.log10(rt.rhodet.max()))
            ax1 = fig.add_subplot(1, 2, 1, aspect="equal")
            ax2 = fig.add_subplot(1, 2, 2, aspect="equal")
            for jj in range(ns):
                # -----------plot data phase tensors---------------
                # check for nan in the array cause it messes with the max
                pt.phimax = np.nan_to_num(pt.phimax)

                # scale the ellipse
                eheightd = pt.phimin[jj] / pt.phimax.max() * ptsize
                ewidthd = pt.phimax[jj] / pt.phimax.max() * ptsize

                # make the ellipse
                ellipd = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=ewidthd,
                    height=eheightd,
                    angle=pt.azimuth[jj],
                )
                # color ellipse:
                if ecolor == "phimin":
                    cvar = (pt.phimin[jj] / (np.pi / 2) - ptmin) / (ptmax - ptmin)
                    if abs(cvar) > 1:
                        ellipd.set_facecolor((1, 0, 0.1))
                    else:
                        ellipd.set_facecolor((1, 1 - abs(cvar), 0.1))
                if ecolor == "beta":
                    cvar = (abs(pt.beta[jj]) - ptmin) / (ptmax - ptmin)
                    if abs(cvar) > 1:
                        ellipd.set_facecolor((1, 1, 0.1))
                    else:
                        ellipd.set_facecolor((1 - cvars, 1 - cvars, 1))

                ax1.add_artist(ellipd)

                # -----------plot data resistivity tensors---------------
                rt.rhomax = np.nan_to_num(rt.rhomax)
                rheightd = rt.rhomin[jj] / rt.rhomax.max() * rtsize
                rwidthd = rt.rhomax[jj] / rt.rhomax.max() * rtsize

                rellipd = Ellipse(
                    (ewarr[jj], nsarr[jj]),
                    width=rwidthd,
                    height=rheightd,
                    angle=rt.rhoazimuth[jj],
                )
                # color ellipse:
                cvar = (np.log10(rt.rhodet[jj]) - rtmin) / (rtmax - rtmin)
                if cvar > 0.5:
                    if cvar > 1:
                        rellipd.set_facecolor((0, 0, 1))
                    else:
                        rellipd.set_facecolor((1 - abs(cvar), 1 - abs(cvar), 1))
                else:
                    if cvar < -1:
                        rellipd.set_facecolor((1, 0, 0))
                    else:
                        rellipd.set_facecolor((1, 1 - abs(cvar), 1 - abs(cvar)))

                ax2.add_artist(rellipd)

            for aa, ax in enumerate([ax1, ax2]):
                ax.set_xlim(ewarr.min() - xpad, ewarr.max() + xpad)
                ax.set_ylim(nsarr.min() - xpad, nsarr.max() + xpad)
                ax.grid("on")
                ax.set_xlabel("easting (km)", fontdict={"size": 10, "weight": "bold"})

                if aa == 1:
                    pylab.setp(ax.get_yticklabels(), visible=False)
                else:
                    ax.set_ylabel(
                        "northing (km)", fontdict={"size": 10, "weight": "bold"}
                    )
                #                cbax=mcb.make_axes(ax,shrink=.8,pad=.15,orientation='horizontal',
                #                               anchor=(.5,1))
                # l,b,w,h
                #                cbax=fig.add_axes([.1,.95,.35,.05])
                if aa == 0:
                    cbax = fig.add_axes([0.12, 0.97, 0.31, 0.02])
                    cbx = mcb.ColorbarBase(
                        cbax,
                        cmap=ptcmap,
                        norm=Normalize(
                            vmin=ptmin * 180 / np.pi, vmax=ptmax * 180 / np.pi
                        ),
                        orientation="horizontal",
                        format="%.2g",
                    )

                    cbx.set_label("Phase (deg)", fontdict={"size": 7, "weight": "bold"})
                if aa == 1:
                    cbax = fig.add_axes([0.59, 0.97, 0.31, 0.02])
                    cbx = mcb.ColorbarBase(
                        cbax,
                        cmap=rtcmapr,
                        norm=Normalize(vmin=10 ** rtmin, vmax=10 ** rtmax),
                        orientation="horizontal",
                        format="%.2g",
                    )

                    cbx.set_label(
                        "App. Res. ($\Omega \cdot$m)",
                        fontdict={"size": 7, "weight": "bold"},
                    )
                    cbx.set_ticks((10 ** rtmin, 10 ** rtmax))
            plt.show()


def readModelFile(mfile, ncol=7):
    """
    read in a model file as x-north, y-east, z-positive down
    """

    mfid = file(mfile, "r")
    mlines = mfid.readlines()

    # get info at the beggining of file
    info = mlines[0].strip().split()
    infodict = dict([(info[0][1:], info[1]), (info[2], info[3]), (info[4], info[5])])

    # get lengths of things
    nx, ny, nz, nn = np.array(mlines[1].strip().split(), dtype=np.int)

    # make empty arrays to put stuff into
    xarr = np.zeros(nx)
    yarr = np.zeros(ny)
    zarr = np.zeros(nz)
    resarr = np.zeros((nx, ny, nz))

    mm = 0
    nn = 2
    while mm < nx:
        xline = mlines[nn].strip().split()
        for xx in xline:
            xarr[mm] = float(xx)
            mm += 1
        nn += 1

    mm = 0
    while mm < ny:
        yline = mlines[nn].strip().split()
        for yy in yline:
            yarr[mm] = float(yy)
            mm += 1
        nn += 1

    mm = 0
    while mm < nz:
        zline = mlines[nn].strip().split()
        for zz in zline:
            zarr[mm] = float(zz)
            mm += 1
        nn += 1

    # put the grids into coordinates relative to the center of the grid
    nsarr = xarr.copy()
    nsarr[: int(nx / 2)] = -np.array(
        [xarr[ii : int(nx / 2)].sum() for ii in range(int(nx / 2))]
    )
    nsarr[int(nx / 2) :] = (
        np.array([xarr[int(nx / 2) : ii + 1].sum() for ii in range(int(nx / 2), nx)])
        - xarr[int(nx / 2)]
    )

    ewarr = yarr.copy()
    ewarr[: int(ny / 2)] = -np.array(
        [yarr[ii : int(ny / 2)].sum() for ii in range(int(ny / 2))]
    )
    ewarr[int(ny / 2) :] = (
        np.array([yarr[int(ny / 2) : ii + 1].sum() for ii in range(int(ny / 2), ny)])
        - yarr[int(ny / 2)]
    )

    zdepth = np.array([zarr[0 : ii + 1].sum() - zarr[0] for ii in range(nz)])

    mm = 0
    for kk in range(nz):
        for jj in range(ny):
            for ii in range(nx):
                resarr[(nx - 1) - ii, jj, kk] = float(mlines[nn + mm].strip())
                mm += 1

    return nsarr, ewarr, zdepth, resarr, infodict


def plotDepthSlice(
    datafn,
    modelfn,
    savepath=None,
    map_scale="km",
    ew_limits=None,
    ns_limits=None,
    depth_index=None,
    fig_dimensions=[4, 4],
    dpi=300,
    font_size=7,
    climits=(0, 4),
    cmap="jet_r",
    plot_grid="n",
    cb_dict={},
):
    """
    plot depth slices
    """

    # create a path to save figure to if it doesn't already exist
    if savepath is not None:
        if not os.path.exists(savepath):
            os.mkdir(savepath)

    # make map scale
    if map_scale == "km":
        dscale = 1000.0
    elif map_scale == "m":
        dscale = 1.0

    # read in data file to station locations
    period, zz, zzerr, ns, ew, slst = readDataFile(datafn)

    # scale the station locations to the desired units
    ns /= dscale
    ew /= dscale

    # read in model file
    x, y, z, resarr, idict = readModelFile(modelfn)

    # scale the model grid to desired units
    x /= dscale
    y /= dscale
    z /= dscale

    # create an list of depth slices to plot
    if depth_index is None:
        zrange = range(z.shape[0])
    elif isinstance(depth_index, int):
        zrange = [depth_index]
    elif isinstance(depth_index, list):
        zrange = depth_index

    # set the limits of the plot
    if ew_limits is None:
        xlimits = (np.floor(ew.min()), np.ceil(ew.max()))
    else:
        xlimits = ew_limits

    if ns_limits is None:
        ylimits = (np.floor(ns.min()), np.ceil(ns.max()))
    else:
        ylimits = ns_limits

    # make a mesh grid of north and east
    north1, east1 = np.meshgrid(x, y)

    fdict = {"size": font_size + 2, "weight": "bold"}

    cblabeldict = {
        -2: "$10^{-3}$",
        -1: "$10^{-1}$",
        0: "$10^{0}$",
        1: "$10^{1}$",
        2: "$10^{2}$",
        3: "$10^{3}$",
        4: "$10^{4}$",
        5: "$10^{5}$",
        6: "$10^{6}$",
        7: "$10^{7}$",
        8: "$10^{8}$",
    }

    plt.rcParams["font.size"] = font_size
    for ii in zrange:
        fig = plt.figure(ii, figsize=fig_dimensions, dpi=dpi)
        plt.clf()
        ax1 = fig.add_subplot(1, 1, 1, aspect="equal")
        ax1.pcolormesh(
            east1,
            north1,
            np.log10(np.rot90(resarr[:, :, ii], 3)),
            cmap=cmap,
            vmin=climits[0],
            vmax=climits[1],
        )

        # plot the stations
        for ee, nn in zip(ew, ns):
            ax1.text(
                ee,
                nn,
                "*",
                verticalalignment="center",
                horizontalalignment="center",
                fontdict={"size": 5, "weight": "bold"},
            )

        # set axis properties
        ax1.set_xlim(xlimits)
        ax1.set_ylim(ylimits)
        ax1.xaxis.set_minor_locator(MultipleLocator(100 * 1.0 / dscale))
        ax1.yaxis.set_minor_locator(MultipleLocator(100 * 1.0 / dscale))
        ax1.set_ylabel("Northing (" + map_scale + ")", fontdict=fdict)
        ax1.set_xlabel("Easting (" + map_scale + ")", fontdict=fdict)
        ax1.set_title(
            "Depth = {:.3f} ".format(z[ii]) + "(" + map_scale + ")", fontdict=fdict
        )

        # plot the grid if desired
        if plot_grid == "y":
            for xx in x:
                ax1.plot([y.min(), y.max()], [xx, xx], lw=0.1, color="k")
            for yy in y:
                ax1.plot([yy, yy], [x.min(), x.max()], lw=0.1, color="k")

        # plot the colorbar
        try:
            cb_dict["orientation"]
        except KeyError:
            cb_dict["orientation"] = "horizontal"

        if cb_dict["orientation"] == "horizontal":
            try:
                ax2 = fig.add_axes(cb_dict["position"])
            except KeyError:
                ax2 = fig.add_axes(
                    (
                        ax1.axes.figbox.bounds[3] - 0.225,
                        ax1.axes.figbox.bounds[1] + 0.05,
                        0.3,
                        0.025,
                    )
                )

        elif cb_dict["orientation"] == "vertical":
            try:
                ax2 = fig.add_axes(cb_dict["position"])
            except KeyError:
                ax2 = fig.add_axes(
                    (
                        ax1.axes.figbox.bounds[2] - 0.15,
                        ax1.axes.figbox.bounds[3] - 0.21,
                        0.025,
                        0.3,
                    )
                )

        cb = mcb.ColorbarBase(
            ax2,
            cmap=cmap,
            norm=Normalize(vmin=climits[0], vmax=climits[1]),
            orientation=cb_dict["orientation"],
        )

        if cb_dict["orientation"] == "horizontal":
            cb.ax.xaxis.set_label_position("top")
            cb.ax.xaxis.set_label_coords(0.5, 1.3)

        elif cb_dict["orientation"] == "vertical":
            cb.ax.yaxis.set_label_position("right")
            cb.ax.yaxis.set_label_coords(1.25, 0.5)
            cb.ax.yaxis.tick_left()
            cb.ax.tick_params(axis="y", direction="in")

        cb.set_label("Resistivity ($\Omega \cdot$m)", fontdict={"size": font_size})
        cb.set_ticks(np.arange(climits[0], climits[1] + 1))
        cb.set_ticklabels(
            [cblabeldict[cc] for cc in np.arange(climits[0], climits[1] + 1)]
        )

        if savepath is not None:

            fig.savefig(
                os.path.join(savepath, "Depth_{}_{:.4f}.png".format(ii, z[ii])), dpi=dpi
            )
            fig.clear()
            plt.close()

        else:
            pass


def computeMemoryUsage(nx, ny, nz, n_stations, n_zelements, n_period):
    """
    compute the memory usage of a model

    Arguments:
    ----------
        **nx** : int
                 number of cells in N-S direction

        **ny** : int
                 number of cells in E-W direction

        **nz** : int
                 number of cells in vertical direction including air layers (7)

        **n_stations** : int
                         number of stations

        **n_zelements** : int
                          number of impedence tensor elements either 4 or 8

        **n_period** : int
                       number of periods to invert for

    Returns:
    --------
        **mem_req** : float
                      approximate memory useage in GB
    """

    mem_req = 1.2 * (
        8 * (n_stations * n_period * n_zelements) ** 2
        + 8 * (nx * ny * nz * n_stations * n_period * n_zelements)
    )
    return mem_req * 1e-9
