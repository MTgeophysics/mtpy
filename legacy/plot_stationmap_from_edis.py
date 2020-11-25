#!/usr/bin/env python

"""
use mpl_toolkits.basemap to create a map to show MT station location

"""

import fnmatch
import os
import sys
import inspect
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap

# If use anaconda python, you can install basemap packages
#     conda install basemap
#     conda install -c conda-forge basemap-data-hires
#     OR
#     sudo pip install https://github.com/matplotlib/basemap/archive/master.zip


import mtpy.core.edi as EDI
import mtpy.utils.filehandling as MTfh

from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
_logger = MtPyLog.get_mtpy_logger(__name__)


def makemap(edilist, mapstretchfactor, symbolsize, labelsize, showlabel):

    this_fun_name = inspect.getframeinfo(inspect.currentframe())[2]
    _logger.debug("starting %s", this_fun_name)
    _logger.debug("mapstretchfactor, symbolsize, labelsize, showlabel")
    _logger.debug("%s %s %s %s", mapstretchfactor, symbolsize, labelsize, showlabel)

    lats = []
    lons = []
    names = []

    for i in edilist:
        e = EDI.Edi()
        e.read_edi_file(i)
        lats.append(e.lat)
        lons.append(e.lon)
        names.append(e.Header.dataid.lower())
        # names.append(e.head['dataid'].lower())

    coords = np.zeros((len(edilist), 2))
    coords[:, 0] = lats
    coords[:, 1] = lons

    latrange = max(lats) - min(lats)
    lonrange = max(lons) - min(lons)

    # center point for projection:
    c = [np.mean(lats), np.mean(lons)]

    # -----------------------
    # Matplotlib options
    mpl.rcParams["font.size"] = 10.0
    mpl.rcParams["axes.labelsize"] = 8.0
    mpl.rcParams["xtick.labelsize"] = 6.0
    mpl.rcParams["ytick.labelsize"] = 6.0

    plt.close("all")

    fig = plt.figure()  # figsize=(5.7,4.3))
    plt.subplots_adjust(
        left=0.2, right=0.9, top=0.90, bottom=0.1, wspace=0.15, hspace=0.05
    )
    ax = plt.subplot(111)

    stretch = float(mapstretchfactor)
    total_latmin = max((min(lats) - stretch * latrange), -90)
    total_lonmin = max((min(lons) - stretch * lonrange), -180)
    total_latmax = min((max(lats) + stretch * latrange), 90)
    total_lonmax = min((max(lons) + stretch * lonrange), 180)

    total_latrange = total_latmax - total_latmin
    total_lonrange = total_lonmax - total_lonmin

    # determine number of axes labels:
    maximumlabels = 5
    latnum = maximumlabels
    lonnum = maximumlabels
    lonlat_stretch = total_lonrange / total_latrange
    if int(lonlat_stretch) > 2:
        # significantly more long than lat
        # could be 0 factor = int(int(lonlat_stretch) / 2.)
        factor = lonlat_stretch / 2.0
        _logger.debug("factor=%s", factor)
        latnum = int(maximumlabels / factor) + 1
        lonnum = maximumlabels
    elif int(lonlat_stretch) < 0.5:
        # significantly more long than lat
        # factor = int(int(1. / lonlat_stretch) / 2.)
        factor = (1.0 / lonlat_stretch) / 2.0
        _logger.debug("factor=%s", factor)
        lonnum = int(maximumlabels / factor) + 1
        latnum = maximumlabels

    m = Basemap(
        projection="merc",
        lon_0=c[1],
        lat_0=c[0],
        lat_ts=c[0],
        llcrnrlat=total_latmin,
        urcrnrlat=total_latmax,
        llcrnrlon=total_lonmin,
        urcrnrlon=total_lonmax,
        rsphere=6371200.0,
        resolution="h",
    )  # ,ax=ax)

    lons.append(total_lonmin)
    lons.append(total_lonmax)
    lats.append(total_latmin)
    lats.append(total_latmax)

    xgrid, ygrid = m(*np.meshgrid(lons, lats))
    xdiff = xgrid.max() - xgrid.min()
    ydiff = ygrid.max() - ygrid.min()

    largest_extent = max(ydiff, xdiff)

    m.drawcoastlines(linewidth=0.25, ax=ax)
    m.drawcountries(linewidth=0.25, ax=ax)
    m.fillcontinents(color="coral", lake_color="aqua", ax=ax)
    m.drawmapboundary(fill_color="aqua", ax=ax)

    m.drawparallels(
        [
            round(i, 3)
            for i in np.linspace(total_latmin, total_latmax, latnum + 1, False)
        ][1:],
        labels=[1, 0, 0, 0],
        fmt="%.1f",
    )
    m.drawmeridians(
        [
            round(i, 3)
            for i in np.linspace(total_lonmin, total_lonmax, lonnum + 1, False)
        ][1:],
        labels=[0, 0, 0, 1],
        fmt="%.1f",
    )
    m.drawrivers()
    m.etopo()

    m.drawmapscale(
        total_lonmax - 0.15 * total_lonrange,
        total_latmax - 0.2 * total_latrange,
        c[1],
        c[0],
        2 * 10 ** (int(np.log10(largest_extent / 1000.0)) - 1),
        barstyle="simple",
        labelstyle="simple",
        fontsize=12,
        fontcolor="k",
        fillcolor1="r",
        fillcolor2="g",
        ax=None,
        format="%d",
    )

    for x, y, name in zip(xgrid[0, :], ygrid[:, 0], names):
        plt.plot(x, y, "v", ms=symbolsize, color="k", label=name)
        if showlabel is True:
            plt.text(
                x,
                y,
                name,
                fontsize=labelsize,
                ha="center",
                va="bottom",
                color="k",
                backgroundcolor="grey",
            )  # , labelfontsize=5)

    plt.title("locations of {0} MT stations".format(len(names)))

    f1 = "tmp/station_locations_map.png"
    f2 = "tmp/station_locations_map.svg"

    f1 = MTfh.make_unique_filename(f1)
    f2 = MTfh.make_unique_filename(f2)

    plt.savefig(f1, format="png", dpi=200)
    plt.savefig(f2, format="svg", transparent=True)

    plt.show()

    return os.path.abspath(f1), os.path.abspath(f2)


##########################################################################


def main():
    if len(sys.argv) < 2:
        _logger.debug(
            """\n\tusage: \n
            python plot_stationmap_from_edis.py <edi folder> [optional_arguments]

        list of optional arguments:

            - <symbolsize (int,24)>
            - <mapstretchfactor (float,1.)>
            - <labelsize (int,14)>
            - <showlabel (bool,True)>


        Output:

        2 files containing a geographical map (Mercator) with all station
        locations, which could be extracted from the EDI files in the folder
        provided:

        station_locations_map.png/svg

        (no overwriting of existing files with these names)


        """
        )
        return

    edifolder = sys.argv[1]
    edilist = []
    try:
        # if not op.isdir(edifolder):
        #     raise
        _logger.debug(edifolder)
        edifiles = os.listdir(edifolder)
        _logger.debug(edifiles)
        edilist = fnmatch.filter(os.listdir(edifolder), "*.[Ee][Dd][Ii]")
        edilist = [os.path.abspath(os.path.join(edifolder, i)) for i in edilist]
        if len(edilist) == 0:
            raise
    except:
        _logger.debug("No EDI files in folder %s", edifolder)
        return

    try:
        symbolsize = int(float(sys.argv[2]))
    except:
        _logger.debug("cannot read symbolsize value - using default")
        symbolsize = 24

    try:
        mapstretchfactor = float(sys.argv[3])
    except:
        _logger.debug("cannot read mapstretchfactor value - using default")
        mapstretchfactor = 1

    try:
        labelsize = int(float(sys.argv[4]))
    except:
        _logger.debug("cannot read labelsize value - using default")
        labelsize = 14

    showlabel = sys.argv[5]
    if showlabel.lower() == "true":
        showlabel = True
    else:
        showlabel = False

    f1, f2 = makemap(edilist, mapstretchfactor, symbolsize, labelsize, showlabel)

    _logger.info("wrote map files\n\t %s \n and\n\t %s\n", f1, f2)


##########################################################################
# Usage example:
# python mtpy/imaging/plot_stationmap_from_edis.py tests/data/edifiles/ 12 2.0 12 true
# python mtpy/imaging/plot_stationmap_from_edis.py  /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 12 2.0 12 true
# -------------------------------------------------------------------------
if __name__ == "__main__":
    main()
