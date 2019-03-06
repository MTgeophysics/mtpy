import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pytest


def demo():
    delta = 0.025
    x = y = np.arange(-3.0, 3.0, delta)
    X, Y = np.meshgrid(x, y)
    Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)

    Z = Z2 - Z1  # difference of Gaussians

    my_extent = (-5, 5, -5, 5)
    im = plt.imshow(Z, interpolation='bilinear', cmap=cm.RdYlGn,
                    origin='lower', extent=my_extent,
                    vmax=abs(Z).max(), vmin=-abs(Z).max())

    plt.show()


import sys
import os

from osgeo import gdal, osr

from gdalconst import *
import numpy as np


# import cartopy.crs as ccrs


def plot_geotiff(geofile='/e/Data/uncoverml/GA-cover2/PM_Gravity.tif', show=True):
    if not os.path.isfile(geofile):
        pytest.skip("file not found {}".format(geofile))
    # Register drivers
    gdal.AllRegister()

    # Open image
    ds = gdal.Open(geofile, GA_ReadOnly)

    if ds is None:
        raise Exception('Could not open image file %s' % (geofile))

    # get image size
    rows = ds.RasterYSize
    cols = ds.RasterXSize
    numbands = ds.RasterCount


    #     print 'rows= %s, cols= %s, number of bands = %s' %(str(rows), str(cols), str(numbands))
    #     print ("********************")

    # get projection and resolution info of the raster
    proj = ds.GetProjection()
    inproj = osr.SpatialReference()
    inproj.ImportFromWkt(proj)

    print(("Input Spatial Ref Projection : ", str(inproj)))

    # these need cartopy, pyepsg , still has error
    #   File "/usr/lib/python2.7/xml/etree/ElementTree.py", line 1517, in _raiseerror
    #  xml.etree.ElementTree.ParseError: not well-formed (invalid token): line 60, column 20
    # projcs = inproj.GetAuthorityCode('PROJCS')
    # projection = ccrs.epsg(projcs)
    # print(projection)


    transform = ds.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]


    #my_ext = (119.967, 121.525, -28.017, -26.955)

    my_ext = (transform[0], transform[0] + ds.RasterXSize * transform[1],
              transform[3] + ds.RasterYSize * transform[5], transform[3])
    print(my_ext)

    #     print ("Projection Info = %s"%(proj))
    #     print ("xOrigin = %s,  yOrigin = %s "%(xOrigin, yOrigin))
    #     print ("pixelWidth = %s,  pixelHeight = %s "%(pixelWidth, pixelHeight))

    # Read the data into numpy array
    numarray = []
    for i in range(1, numbands + 1):
        band = ds.GetRasterBand(i)  # the very first band is i=1
        data = band.ReadAsArray(0, 0, cols, rows)  # .astype('float32')
        numarray.append(data)

    # Once we're done, close properly the dataset
    ds = None

    f, ax = plt.subplots(1, 1, figsize=(10, 80))

    # ax.imshow(numarray[0]) # no georef info, just a gridded image origin is upper


    my_cmap = 'jet' # cm.RdYlGn

    # ValueError: Possible values are:
    # Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r,
    # CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r,
    # Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r,
    # Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r,
    # PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r,
    # RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r,
    # Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r,
    # YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r,
    # bone, bone_r, brg, brg_r, bwr, bwr_r, cool, cool_r, coolwarm, coolwarm_r,
    # copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r,
    # gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow,
    # gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot,
    # gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r,
    # jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r,
    # pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r,
    # spectral, spectral_r, spring, spring_r, summer, summer_r, terrain, terrain_r, viridis, viridis_r,
    # winter, winter_r

    ax.imshow(numarray[0], extent=my_ext, cmap=my_cmap)
    ax.set_title('%s\n' % ('Image ' + geofile))


    if show is True:
        plt.show()

    return ax  # use ax to do further plots-overlay


#############################################################
# python examples/sandpit/plot_geotiff_imshow.py 
# python examples/sandpit/plot_geotiff_imshow.py /e/Data/Dexport/CleanImages_LakeGorge_2015-10-02T100000.tiff 


if __name__ == "__main__":

    if len(sys.argv) > 1:
        tiffile = sys.argv[1]
        plot_geotiff(tiffile)
    else:  # show default geotif file
        plot_geotiff()
