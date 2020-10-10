"""
Description:
    Given a set of EDI files plot the Penetration Depth vs the station_location.
    Note that the values of periods within10% tolerance (ptol=0.1) are considered as equal.
    Setting a smaller value for ptol(=0.05) may result less MT sites data included.

Usage:
    python mtpy/imaging/penetration_depth3d.py /path2/edi_files_dir/  period_index

Author: fei.zhang@ga.gov.au
Date:   2017-01-23
"""

import glob
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import mtpy.core.mt as mt
from mtpy.imaging.penetration import get_index, load_edi_files, Depth3D
from mtpy.utils.mtpy_decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog
import logging

# mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'

# mpl.rcParams['figure.figsize'] = [20, 10]

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
_logger = MtPyLog.get_mtpy_logger(__name__)
#_logger.setLevel(logging.DEBUG)
_logger.setLevel(logging.INFO)



# logger =
# MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) #
# specific


# This is the major function to be maintained!!!
# use the Zcompotent=[det, zxy, zyx]
def plot_latlon_depth_profile(edi_dir, period, zcomponent='det', showfig=True,
                              savefig=True, savepath = None, fig_dpi=400,
                              fontsize=14, file_format='png',ptol=0.1):
    """
    MT penetration depth profile in lat-lon coordinates with pixelsize = 0.002
    :param savefig:
    :param showfig:
    :param edi_dir:
    :param period:
    :param zcomponent:
    :return:
    """

    # edi_dir = "/Softlab/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir="E:/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir=r"E:\Githubz\mtpy2\examples\data/edi2"

    # 1 get a list of edi files, which are suppose to be in a profile.
    # edifiles = glob.glob(os.path.join(edi_dir, '*.edi'))

    # logger.debug("edi files: %s", edifiles)

    edis = load_edi_files(edi_dir)

    image = Depth3D(edis=edis, period=period, rho=zcomponent, ptol=ptol)
    if isinstance(period, int):  # period is considered as an index
        image.plot(period_by_index=True, fontsize=fontsize)
    elif isinstance(period, float):  # period is considered as the actual value of period in second
        image.plot(fontsize=fontsize)
    else:
        raise Exception("Wrong type of the parameter period, %s" % period)

    if showfig is True:
        image.show()

    if savefig:
        if savepath is None:
            savepath = 'C:/tmp'
        savefn = 'P3Depth_Period%s.%s' % (image.get_period_fmt(),file_format)
        path2savefile = os.path.join(savepath, savefn)
        image.export_image(path2savefile, dpi=fig_dpi, bbox_inches='tight')

    # may want to remove the following 2 lines
    # plt.clf()
    # plt.close()
    return


@deprecated("this function is redundant as matplotlib.cm as the inverted version of color maps")
def reverse_colourmap(cmap, name='my_cmap_r'):
    """
    In: cmap, name
    Out: my_cmap_r

    Explanation: http://stackoverflow.com/questions/3279560/invert-colormap-in-matplotlib
    """
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1 - t[0], t[2], t[1]))
        reverse.append(sorted(data))

    linear_l = dict(list(zip(k, reverse)))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, linear_l)
    return my_cmap_r


@deprecated("please use get_index() in mtpy.imaging.penetration instead")
def get_index2(lat, lon, ref_lat, ref_lon, pixelsize):
    """ Mapping of lat lon to a grid
    :param lat:
    :param lon:
    :param ref_lon:
    :param ref_lat:
    :param pixelsize:
    :return:
    """

    index_x = (lon - ref_lon) / pixelsize
    index_y = (lat - ref_lat) / pixelsize

    return int(index_x), int(index_y)


#########################################################


def plot_bar3d_depth(edifiles, per_index, whichrho='det'):
    """
    plot 3D bar of penetration depths
    For a given freq/period index of a set of edifiles/dir,
    the station,periods, pendepth,(lat, lon) are extracted
    the geo-bounding box calculated, and the mapping from stations to grids
    is constructed and plotted.

    :param whichrho: z component either 'det', 'zxy' or 'zyx'
    :param edifiles: an edi_dir or list of edi_files
    :param per_index: period index number 0,1,2

    :return:
    """

    if os.path.isdir(edifiles):
        edi_dir = edifiles  # "E:/Githubz/mtpy2/tests/data/edifiles/"
        edifiles = glob.glob(os.path.join(edi_dir, '*.edi'))
        _logger.debug(edifiles)
    else:
        # Assume edifiles is [a list of files]
        pass

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    _logger.debug("The scaling parameter=%.6f" % scale_param)

    # per_index=0,1,2,....
    periods = []
    pen_depth = []
    stations = []
    latlons = []

    for afile in edifiles:
        mt_obj = mt.MT(afile)

        latlons.append((mt_obj.lat, mt_obj.lon))

        # the attribute Z
        zeta = mt_obj.Z

        if per_index >= len(zeta.freq):
            raise Exception(
                "Error: input period index must be less than the number of freqs in zeta.freq=%s", len(
                    zeta.freq))
        per = 1.0 / zeta.freq[per_index]
        periods.append(per)

        if whichrho == 'det':  # the 2X2 complex Z-matrix's determinant abs value
            # determinant value at the given period index
            det2 = np.abs(zeta.det[per_index])
            penetration_depth = -scale_param * np.sqrt(0.2 * per * det2 * per)
        elif whichrho == 'zxy':
            penetration_depth = - scale_param * \
                                np.sqrt(zeta.resistivity[per_index, 0, 1] * per)
        elif whichrho == 'zyx':
            penetration_depth = - scale_param * \
                                np.sqrt(zeta.resistivity[per_index, 1, 0] * per)

        pen_depth.append(penetration_depth)

        stations.append(mt_obj.station)

    # return (stations, periods, pen_depth, latlons)

    lats = [tup[0] for tup in latlons]
    lons = [tup[1] for tup in latlons]
    minlat = min(lats)
    maxlat = max(lats)
    minlon = min(lons)
    maxlon = max(lons)

    pixelsize = 0.002  # degree 0.001 = 100meters
    shift = 3
    ref_lat = minlat - shift * pixelsize
    ref_lon = minlon - shift * pixelsize

    xgrids = maxlon - minlon
    ygrids = maxlat - minlat

    # nx = xgrids / pixelsize
    # ny = ygrids / pixelsize

    # import matplotlib.pyplot as plt
    # import numpy as np

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')

    xpos = []  # a seq (1,2,3,4,5,6,7,8,9,10)
    ypos = []  # a seq [2,3,4,5,1,6,2,1,7,2]
    dz = []
    for iter, pair in enumerate(latlons):
        xpos.append(
            get_index(
                pair[0],
                pair[1],
                ref_lat,
                ref_lon,
                pixelsize)[0])
        ypos.append(
            get_index(
                pair[0],
                pair[1],
                ref_lat,
                ref_lon,
                pixelsize)[1])
        dz.append(np.abs(pen_depth[iter]))
        # dz.append(-np.abs(pen_depth[iter]))

    num_elements = len(xpos)
    zpos = np.zeros(num_elements)  # zpos = [0,0,0,0,0,0,0,0,0,0]
    dx = np.ones(num_elements)
    dy = np.ones(num_elements)
    # dz = [1,2,3,4,5,6,7,8,9,10]

    # print(xpos)
    # print(ypos)
    # print(zpos)
    #
    # print(dx)
    # print(dy)
    # print(dz)
    ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='r')

    # ax1

    plt.title(
        'Penetration Depth (Meter) Across Stations for period= %6.3f Seconds' %
        periods[0], fontsize=16)
    plt.xlabel('Longitude(deg-grid)', fontsize=16)
    plt.ylabel('Latitude(deg-grid)', fontsize=16)
    # plt.zlabel('Penetration Depth (m)')
    # bar_width = 0.4
    # plt.xticks(index + bar_width / 2, stations, rotation='horizontal', fontsize=16)
    # plt.legend()
    #
    # # plt.tight_layout()
    # plt.gca().xaxis.tick_top()
    # plt.show()

    plt.show()

# ======================
def get_penetration_depths_from_edi_file(edifile, rholist=['det']):
    """Compute the penetration depths of an edi file
    :param edifile: input edifile
    :param rholist: flag the method to compute penetration depth: det zxy zyx
    :return: a tuple:(station_lat, statoin_lon, periods_list, pendepth_list)
    """
    _logger.debug("processing the edi file %s", edifile)

    mt_obj = mt.MT(edifile)
    zeta = mt_obj.Z  # the attribute Z represent the impedance tensor 2X2 matrix
    freqs = zeta.freq  # frequencies

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    _logger.debug("the scale parameter should be 355.88127 =?= %s", scale_param)

    # The periods array
    periods = 1.0 / freqs

    if 'zxy' in rholist:
        # One of the 4-components: XY
        penetration_depth = scale_param * \
                            np.sqrt(zeta.resistivity[:, 0, 1] * periods)

    if 'zyx' in rholist:
        penetration_depth = scale_param * \
                            np.sqrt(zeta.resistivity[:, 1, 0] * periods)

    if 'det' in rholist:
        # determinant is |Zeta|**2
        det2 = np.abs(zeta.det[0])
        penetration_depth = scale_param * \
                            np.sqrt(0.2 * periods * det2 * periods)

    latlong_d = (mt_obj.lat, mt_obj.lon, periods, penetration_depth)
    return latlong_d

def create_penetration_depth_csv(edi_dir, outputcsv, zcomponent='det'):
    """ Loop over all edi files, and create a csv file with the columns:
    Header Lat, Lon, per0, per1,per2,.....

    TODO:
    calculate pen-depth for each period, and write into a file for each period, even if non-equal freq cross edi files.
    Moved this function into edi_collection.create_penetration_depth_csv()

    lat, lon, pendepth0, pendepth1, ...
    :param edi_dir: path_to_edifiles_dir
    :param zcomponent: det | zxy  | zyx
    :param outputcsv: path2output.csv file
    :return:
    """
    import csv

    if not os.path.isdir(edi_dir):
        _logger.error("input edi directory not exists", edi_dir)
        raise Exception("MTPy Exception: EDI Dir not exist")

    edi_files = glob.glob(os.path.join(edi_dir, "*.edi"))

    _logger.debug(edi_files)

    # the first period list as a reference for checking other stations period
    periods_list0 = None
    latlon_dep = []  # CSV to be returned
    for afile in edi_files:
        # for efile in edi_files[:2]:
        _logger.debug("processing %s", afile)
        lat, lon, periods, depths = get_penetration_depths_from_edi_file(afile)
        if periods_list0 is None:
            periods_list0 = periods  # initial value assignment
            #depth_string = ','.join(['%.2f' % num for num in depths])
            #latlon_dep.append((lat, lon, depth_string))
            latlon_dep.append(["Lat","Lon"] + list(periods)) #The first line header
            latlon_dep.append([lat, lon] + list(depths))


        # same length and same values.
        elif len(periods) == len(periods_list0) and (periods == periods_list0).all():
            # depth_string = ','.join(['%.2f' % num for num in depths])
            # latlon_dep.append((lat, lon, depth_string))
            latlon_dep.append([lat, lon] + list(depths))
        else:
            _logger.error(
                "MT Periods Not Equal !! %s %s VS %s %s",
                len(periods), periods,  len(periods_list0), periods_list0)
            # raise Exception ("MTPy Exception: Periods Not Equal")
            # pass this edi, let's continue

    # logger.debug(latlon_dep)

    if outputcsv is None:
        _logger.error("Output CSV file must be provided", outputcsv)

    _logger.info("Saving to csv file: %s", outputcsv)
    with open(outputcsv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(latlon_dep)

    return latlon_dep


def create_shapefile(edi_dir, outputfile=None, zcomponent='det'):
    """
    create a shapefile for station, penetration_depths
    :param edi_dir:
    :param outputfile:
    :param zcomponent:
    :return:
    """

    # TODO:

    return outputfile


def plot_many_periods(edidir, n_periods=5):
    from mtpy.core.edi_collection import EdiCollection

    edilist = glob.glob(os.path.join(edidir, '*.edi'))

    ediset = EdiCollection(edilist)
    for period_sec in ediset.all_unique_periods[:n_periods]:
        try:
            # This will enable the loop continue even though for some freq,
            #  cannot interpolate due to not enough data points
            plot_latlon_depth_profile(edidir, period_sec, zcomponent='det', showfig=False)
        except Exception as exwhy:
            print(str(exwhy))


# =============================================================================================
# Usage examples for small, med, large images
# python mtpy/imaging/penetration_depth3d.py tests/data/edifiles/ 2.857s
# python mtpy/imaging/penetration_depth3d.py /e/Datasets/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 0.58s
# python mtpy/imaging/penetration_depth3d.py /e/Datasets/MT_Datasets/GA_UA_edited_10s-10000s 16s [10s, 40s 341s]
#   OR  period index integer
# python mtpy/imaging/penetration_depth3d.py /e/Datasets/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 30
#   python mtpy/imaging/penetration_depth3d.py  examples/data/edi_files/
#   python mtpy/imaging/penetration_depth3d.py  examples/data/edi_files_2/   (non equal frequencies)
# =============================================================================================
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print(("Usage: python %s edi_dir period_sec " % sys.argv[0]))
        print("usage example: python mtpy/imaging/penetration_depth3d.py  examples/data/edi_files/ 10")
        print("usage example: python mtpy/imaging/penetration_depth3d.py  examples/data/edi_files/ 2.857s")
        sys.exit(1)
    elif len(sys.argv) == 2:
        edi_dir= sys.argv[1]

        bname =os.path.basename(os.path.normpath(edi_dir))
        print ("dir base name", bname)
        create_penetration_depth_csv(edi_dir, "/tmp/%s_MT_pen_depths.csv"%bname)

        # plot pendepth over multiple periods
        # plot_many_periods( edi_dir )

    elif len(sys.argv) > 2 and os.path.isdir(sys.argv[1]):
        edi_dir = sys.argv[1]

        per = sys.argv[2]
        print(("The input parameter for period was ", per))

        if per.endswith('s'):
            # try float case
            try:
                period_sec = float(per[:-1])
                print((" Using period value (second)", period_sec))
                plot_latlon_depth_profile(
                    edi_dir, period_sec, zcomponent='det')
            except Exception as ex:
                print(ex)
        else:
            try:
                period_index = int(per)
                print((" Using period index", period_index))
                plot_latlon_depth_profile(
                    edi_dir, period_index, zcomponent='det')
                sys.exit(0)  # done
            except Exception as why:
                raise Exception(
                    "Unable to plot the integer period index, because: %s" %
                    why)

                # plot_gridded_profile(edi_dir, period_index, zcomponent='det')   # 2D image
                # plot_latlon_depth_profile(edi_dir, period_index,zcomponent='det')
                # plot_bar3d_depth(edi_dir, period_index)
    else:
        print("Please provide an edi directory and period_index_OR_value")
