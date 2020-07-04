"""
Displays geotiffs on matplotlib axes as basemap background images.

Revision History:
    brenainn.moushall@ga.gov.au 26-03-2020 15:06:38 AEDT:
        Add function for plotting a geotiff image on a prexisting,
        georefernced axes.
"""

import sys
import os

import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from osgeo import gdal, osr, ogr
from gdalconst import GA_ReadOnly


def plot_geotiff_on_axes(geotiff, axes, extents=None, epsg_code=None,
                         band_number=None, cmap='viridis'):
    """
    Plot a geotiff on a prexisting matplotlib axis that represents a
    georeferenced map. Doesn't return anything - the plotting is done
    on the axis in-place.

    Parameters
    ----------
    geotiff : str
        Full path to the geotiff file to display.
    axes : matplotlib.axes
        The axes to display the image on.
    extents : tuple of float, optional
        Bounding box within which to draw the image of format
        (left, bottom, right, top). *Must be in the same CRS as the
        axes*. The geotiff will be cropped using these extents and then
        displayed within the same extents. If not provided, extents
        will be retrieved from the axes x and y limits, i.e. the image
        will fill the entire axes.
    epsg_code : int, optional
        EPSG code of the axes map CRS. Must be provided if the axes and
        image CRS are different. This is used to reproject the given
        extents to the image CRS, so the extents can be used to crop the
        image. If not provided, it's assumed the CRS of the axes and
        image are the same.
    band_number : int, optional
        The band to display. If None, then all bands will be read.
        Leave as None to display RGB/A rasters.
    cmap : str or matplotlib.colors.Colormap, optional
        Used to color the image data. Defaults to 'viridis'.
    """
    def transform_point(src_srs, dst_srs, x, y):
        ctr = osr.CoordinateTransformation(src_srs, dst_srs)
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x, y)
        point.Transform(ctr)
        return point.GetX(), point.GetY()

    ds = gdal.Open(geotiff, GA_ReadOnly)

    # By default, will fill up as much of the grid axes as possible.
    if extents is None:
        xmin, xmax = axes.get_xlim()
        ymin, ymax = axes.get_ylim()
        extents = xmin, ymin, xmax, ymax
    l, b, r, t = extents

    if epsg_code is not None:
        # Reproject extents to those of geotiff
        ax_srs = osr.SpatialReference()
        ax_srs.ImportFromEPSG(epsg_code)
        img_srs = osr.SpatialReference()
        img_srs.ImportFromWkt(ds.GetProjection())
        l, b = transform_point(ax_srs, img_srs, l, b)
        r, t = transform_point(ax_srs, img_srs, r, t)

    # Work out geotiff extent.
    gt = ds.GetGeoTransform()
    x0, y0, x_pixel_size, y_pixel_size = gt[0], gt[3], gt[1], gt[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    gt_l = x0
    gt_b = y0 + (y_pixel_size * rows)
    gt_r = x0 + (x_pixel_size * cols)
    gt_t = y0

    # Crop the geotiff on edges that extend beyond the grid bounds.
    crop_l = min(l, gt_l)
    crop_b = min(b, gt_b)
    crop_r = min(r, gt_r)
    crop_t = min(t, gt_t)
    r1 = int((crop_t - y0) / y_pixel_size)
    c1 = int((crop_l - x0) / x_pixel_size)
    r2 = int((crop_b - y0) / y_pixel_size)
    c2 = int((crop_r - x0) / x_pixel_size)

    band_range = range(1, ds.RasterCount + 1) if band_number is None \
        else range(band_number, band_number + 1)
    data = []

    for i in band_range:
        band = ds.GetRasterBand(i)
        data.append(band.ReadAsArray(c1, r1, c2 - c1, r2 - r1))

    if len(data) == 1:
        data = data[0]
    else:
        data = np.stack(data, axis=2)

    if epsg_code is not None:
        crop_l, crop_b = transform_point(img_srs, ax_srs, crop_l, crop_b)
        crop_r, crop_t = transform_point(img_srs, ax_srs, crop_r, crop_t)

    axes.imshow(data, cmap=cmap, origin='upper', extent=(crop_l, crop_r, crop_b, crop_t))


def plot_geotiff(geofile='/e/Data/uncoverml/GA-cover2/PM_Gravity.tif', show=True):
    if not os.path.isfile(geofile):
        raise FileNotFoundError("Geotiff not found: {}".format(geofile))
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

    my_cmap = 'jet'  # cm.RdYlGn

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
