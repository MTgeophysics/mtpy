"""
Description:
    For a batch of MT_stations,  plot the Penetration Depth vs the station_location,
    for a given period (1/freq).

Usage:
    python examples/penetration_depth_profile3d.py /path2/edi_files_dir/  period_index

Author: fei.zhang@ga.gov.au
Date:   2017-01-23
"""

import sys
import os
import glob
import numpy as np
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
import matplotlib as mpl


mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'

mpl.rcParams['figure.figsize']=[20,10]

import mtpy.core.mt as mt

from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to config the logger
logger = MtPyLog().get_mtpy_logger(__name__)
#logger = MtPyLog(path2configfile='logging.yml').get_mtpy_logger(__name__) # specific

def plot_3d_profile(edi_dir, period_index, zcomponent='det'): #use the Zcompotent=[det, zxy, zyx]
    #edi_dir = "/Softlab/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir="E:/Githubz/mtpy2/tests/data/edifiles/"
    # edi_dir=r"E:\Githubz\mtpy2\examples\data/edi2"

    #1 get a list of edi files, which are suppose to be in a profile.
    edifiles = glob.glob(os.path.join(edi_dir, '*.edi'))

    logger.debug("edi files: %s", edifiles)

    (stations, periods, pendep, latlons) = get_penetration_depth(edifiles, period_index, whichrho=zcomponent)

    bbox=get_bounding_box(latlons)

    logger.debug("Bounding Box %s", bbox)

    xgrids = bbox[0][1] - bbox[0][0]
    ygrids = bbox[1][1] - bbox[1][0]

    logger.debug("xy grids %s %s", xgrids, ygrids)

    minlat = bbox[1][0]
    minlon = bbox[0][0]

    # Pixel size in Degree:  0.001=100meters, 0.01=1KM 1deg=100KM

    pixelsize = 0.002  # Degree 0.001=100meters, 0.01=1KM 1deg=100KM

    nx = int(np.ceil(xgrids / pixelsize))
    ny = int(np.ceil(ygrids / pixelsize))

    print(nx, ny)

    # make an image bigger than the (nx, ny)
    pad = 4
    nx2 = nx + pad
    ny2 = ny + pad

    # Z = 0.0* np.random.random((nx2,ny2))   # Test data
    # Z=  np.ones((nx2,ny2))
    zdep = np.zeros((ny2, nx2))
    # Z[10, 10]=12
    # Z[11, 20]=20
    # Z[13, 15]=30

    zdep[:, :] = np.nan

    logger.debug("zdep shape %s", zdep.shape)

    for iter, pair in enumerate(latlons):
        print pair
        (xi, yi) = get_index(pair[0], pair[1], minlat, minlon, pixelsize)
        zdep[zdep.shape[0] - yi, xi] = np.abs(pendep[iter])

    plt.imshow(zdep) #,  interpolation='none')
    #plt.imshow(zdep,  interpolation='spline36')
    plt.colorbar()

    plt.show()


# ================================== griddata interpolation

    print(zdep.shape[0], zdep.shape[1])
    grid_x, grid_y = np.mgrid[0:95:96j, 0:83:84j]  # how to param this mesh construction ?

    # print (grid_x, grid_y)
    points = np.zeros((len(latlons), 2))
    values = np.zeros(len(latlons))

    for iter, pair in enumerate(latlons):
        #     print pair

        (i, j) = get_index(pair[0], pair[1], minlat, minlon, pixelsize)
        points[iter, 0] = zdep.shape[0] - j
        points[iter, 1] = i
        values[iter] = np.abs(pendep[iter])

    #grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
    #grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
    grid_z = griddata(points, values, (grid_x, grid_y), method='cubic')

    # plt.imshow(grid_z)

    plt.subplot(111)

    plt.imshow(grid_z, origin='upper')
    plt.plot(points[:, 1], points[:, 0], 'kv', markersize=6) #stations

    plt.title('Cubic Interpolation of Penetration Depth at the Period=%s' % periods[0])  # Cubic
    plt.colorbar()
    plt.gcf().set_size_inches(6, 6)

    plt.show()


def get_bounding_box(latlons):
    """ get min max lat lon from the list of lat-lon-pairs points"""
    lats = [tup[0] for tup in latlons]
    lons = [tup[1] for tup in latlons]

    minlat = min(lats)
    maxlat = max(lats)

    print(minlat, maxlat)

    minlon = min(lons)
    maxlon = max(lons)

    print(minlon, maxlon)

    return ((minlon, maxlon), (minlat, maxlat))



def get_index(lat, lon, minlat, minlon, pixelsize, offset=1):
    index_x = (lon - minlon) / pixelsize
    index_y = (lat - minlat) / pixelsize

    ix = int(round(index_x))
    iy = int(round(index_y))

    print (ix, iy)

    return (ix + offset, iy + offset)


def get_penetration_depth(edi_file_list, per_index,  whichrho='det'): #whichrho=[zxy, zyx, zdeterminant]
    """ input period index and a list of edi files,
    return tuple of lists (stations, periods, penetrationdepth, lat-lons-pairs)
    """

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    # per_index=0,1,2,....
    periods = []
    pendep = []
    stations = []
    latlons = []

    for afile in edi_file_list:
        mt_obj = mt.MT(afile)

        latlons.append((mt_obj.lat, mt_obj.lon))

        # the attribute Z
        zeta = mt_obj.Z

        if per_index >= len(zeta.freq):
            logger.debug("number of frequecies= %s", len(zeta.freq))
            raise Exception("Index out_of_range Error: period index must be less than number of periods in zeta.freq")

        per = 1.0 / zeta.freq[per_index]
        periods.append(per)

        if whichrho == 'zxy':
            penetration_depth = - scale_param * np.sqrt(zeta.resistivity[per_index, 0, 1] * per)
        elif whichrho == 'zyx':
            penetration_depth = - scale_param * np.sqrt(zeta.resistivity[per_index, 1, 0] * per)
        elif whichrho == 'det':
            # determinant
            det2 = np.abs(zeta.det[0][per_index])  # determinant value at the given period index
            penetration_depth = -scale_param * np.sqrt(0.2 * per * det2 * per)
        else:
            logger.critical("unsupported method to compute penetration depth: %s", whichrho)
            sys.exit(100)

        pendep.append(penetration_depth)

        stations.append(mt_obj.station)

    # check_all_period_equal()
    return (stations, periods, pendep, latlons)



# =============================================================================================
# python examples/penetration_depth_3d_profile.py tests/data/edifiles/ 10
# python examples/penetration_depth_3d_profile.py examples/data/edi2/ 10
# =============================================================================================
if __name__=="__main__":

    if len(sys.argv)<2:
        print("Usage: %s edi_dir"%sys.argv[0])
        print ("python examples/penetration_depth_3d_profile.py tests/data/edifiles/ 10 ")
        sys.exit(1)
    elif os.path.isdir(sys.argv[1]):
        edi_dir = sys.argv[1]
        period_index= int(sys.argv[2])
        plot_3d_profile(edi_dir, period_index, zcomponent='det')

    else:
        print("Please provide an edi directory and period_index_list")

