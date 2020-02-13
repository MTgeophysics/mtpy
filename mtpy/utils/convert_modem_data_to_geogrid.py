#! /usr/bin/env python
"""
Description:
    Convert input MODEM data  and resistivity rho files into a georeferenced raster/grid format,
    such as geotiff format, which can be visualized by GIS software.

CreationDate:   1/05/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     1/05/2019   FZ
    LastUpdate:     17/09/2019  FZ fix the geoimage coordinates, upside-down issues 
    LastUpdate:     dd/mm/yyyy
"""

import os
import sys
import argparse

from pyproj import Proj
import gdal
import osr
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mtpy.modeling.modem import Model, Data
from mtpy.utils import gis_tools
from mtpy.contrib.netcdf import nc
from mtpy.utils.mtpylog import MtPyLog
import mtpy.contrib.netcdf.modem_to_netCDF as modem2nc

_logger = MtPyLog.get_mtpy_logger(__name__)


def array2geotiff_writer(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array, epsg_code=4283):
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    # driver = gdal.GetDriverByName('AAIGrid')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg_code)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

    # output to ascii format
    format2 = 'AAIGrid'
    if newRasterfn.endswith(".tif"):
        _newRasterfn = newRasterfn[:-4]
    newRasterfn2 = "%s.asc" % _newRasterfn
    driver2 = gdal.GetDriverByName(format2)
    dst_ds_new = driver2.CreateCopy(newRasterfn2, outRaster)

    return newRasterfn

# TODO (BM): refactor this into test module
def test_array2geotiff(newRasterfn, epsg):
    """
    A  dummpy array of data to be written into a geotiff file. It looks like a image of "GDAL"
    :param newRasterfn:
    :param epsg: 4326, 4283
    :return:
    """
    # rasterOrigin = (-123.25745,45.43013)
    rasterOrigin = (149.298, -34.974)  # Longitude and Lattitude in Aussi continent
    pixelWidth = 0.01
    pixelHeight = -0.01  # this must be negative value, as a Geotiff image's origin is defined as the upper-left corner.

    # Define an image 2D-array: The black=0 pixels trace out GDAL; the bright=1 pixels are white background
    array = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
                      [1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                      [1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
                      [1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                      [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1],
                      [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])

    random = np.random.rand(array.shape[0], array.shape[1])

    array2 = 1000.0 * array + 10.0 * random
    print(array2)

    outfn = array2geotiff_writer(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array2,
                                 epsg_code=epsg)  # write to a raster file

    return outfn


def create_geogrid(data_file, model_file, out_dir,
                   xpad=6, ypad=6, zpad=10, grid_size=7500,
                   center_lat=None, center_lon=None, epsg_code=None,
                   depth_index=None):
    """Generate an output geotiff file and ASCII grid file.

    Args:
        data_file (str): Path to the ModEM .dat file. Used to get the
            grid center point.
        model_file (str): Path to the ModEM .rho file.
        xpad (int): TODO
        ypad (int): TODO
        zpad (int): TODO
        grid_size (int): Pixel resolution in meters.
        epsg_code (int): EPSG code of the model CRS. By default is
            inferred from the grid center point.
        depth_index: A list of integers, eg, [0,2,4] of the depth
            slice's index to be output. If None, all slices are
            selected.
        center_lat: Grid center latitude in degrees.
        center_lon: Grid center longitude in degrees.
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    model = Model()
    model.read_model_file(model_fn=model_file)

    data = Data()
    data.read_data_file(data_fn=data_file)
    center = data.center_point   # see data.py line1406
    center_lat = center.lat.item() if center_lat is None else None
    center_lon = center.lon.item() if center_lon is None else None

    _logger.info("Grid center (lat, lon) = ({}, {})".format(center_lat, center_lon))

    if epsg_code is None:
        zone_number, is_northern, utm_zone = gis_tools.get_utm_zone(center_lat, center_lon)
        epsg_code = gis_tools.get_epsg(center_lat, center_lon)
        _logger.info("Input data epsg code has been inferred as {}".format(epsg_code))

    # get the grid cells' centres (halfshift -cs/2?)
    mgce, mgcn, mgcz = [np.mean([arr[:-1], arr[1:]], axis=0)
                        for arr in [model.grid_east, model.grid_north, model.grid_z]]

    # Get X, y, Z paddings
    # BM: @FZ can we get a better explanation of these abbreviations?
    # padding off big-sized edge cells
    gce, gcn = mgce[xpad:-xpad], mgcn[ypad:-ypad]
    gcz = mgcz[:-zpad]

    _logger.info("E shape = {}, N shape = {}, Z shape = {}"
                 .format(gce.shape, gcn.shape, gcz.shape))

    _logger.info("Data center (east, north) = ({}, {})".format(center.east, center.north))

    #  May need to shift by half cellsize -cs/2
    # In [1]: -164848.1035642 -3750
    # Out[1]: -168598.1035642
    # In [2]: 5611364.73539792 - 3750
    # Out[2]: 5607614.73539792
    origin = (gce[0] + center.east - 0.5 * grid_size, gcn[-1] + center.north - 0.5 * grid_size)
    _logger.info("The Origin (UpperLeft Corner) =".format(origin))

    pixel_width = grid_size
    # This should be negative for geotiff spec, whose origin is at the
    # upper-left corner of image.
    pixel_height = -grid_size

    (target_gridx, target_gridy) = np.meshgrid(np.arange(gce[0], gce[-1], grid_size),
                                               np.arange(gcn[0], gcn[-1], grid_size))

    resgrid_nopad = model.res_model[xpad:-xpad, ypad:-ypad, 0:-zpad]

    depth_index = range(len(gcz)) if depth_index is None else depth_index

    _logger.info("Depth indicies = {}".format(depth_index))

    for di in depth_index:
        output_file = 'DepthSlice%1im.tif'.format(gcz[di])
        output_file = os.path.join(out_dir, output_file)
        # define interpolation function (interpolate in log10 measure-space)
        # See https://docs.scipy.org/doc/scipy-0.16.0/reference/interpolate.html
        interpfunc = RegularGridInterpolator((gce, gcn), np.log10(resgrid_nopad[:, :, di].T))
        # evaluate on the regular grid points, which to be output into geogrid formatted files
        newgridres = 10 ** interpfunc(np.vstack(
            [target_gridx.flatten(), target_gridy.flatten()]).T).reshape(target_gridx.shape)
        _logger.info("New interpolated resistivity grid shape at index {}: {} "
                     .format(di, newgridres.shape))

        # this original image may start from the lower left corner, if so must be flipped.
        # resis_data_flip = resis_data[::-1]  # flipped to ensure the image starts from the upper left corner

        array2geotiff_writer(output_file, origin, pixel_width, pixel_height, newgridres[::-1], epsg_code=epsg_code)

    return output_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('modem_data', help="ModEM data file")
    parser.add_argument('modem_model', help="ModEM model file")
    parser.add_argument('out_dir', help="output directory")
    parser.add_argument('--xpad', type=int, help='xpad value')
    parser.add_argument('--ypad', type=int, help='ypad value')
    parser.add_argument('--grid_size', type=int, help="pixel size in meters")
    parser.add_argument('--lat', type=float, help="grid center latitude in degrees")
    parser.add_argument('--lon', type=float, help="grid center longitude in degrees")
    parser.add_argument('--di', type=int, nargs='*', help="indicies for depth slices to convert, "
                        "eg., [0, 2, 5, 9]")
    parser.add_argument('--epsg', type=int, help="EPSG code for CRS of the model")

    args = parser.parse_args()

    create_geogrid(args.modem_data, args.modem_model, args.out_dir, xpad=args.xpad, ypad=args.ypad,
                   grid_size=args.grid_size, center_lat=args.lat, center_lon=args.lon,
                   depth_index=args.di, epsg_code=args.epsg)
