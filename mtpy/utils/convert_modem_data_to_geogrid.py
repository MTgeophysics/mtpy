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
import math

from pyproj import Proj
import gdal
import osr
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy import ndimage

from mtpy.modeling.modem import Model, Data
from mtpy.utils import gis_tools
from mtpy.utils.mtpylog import MtPyLog

_logger = MtPyLog.get_mtpy_logger(__name__)


def array2geotiff_writer(filename, origin, pixel_width, pixel_height, data,
                         angle=None, epsg_code=4283, center=None, rotate_origin=False):
    gt = [origin[0], pixel_width, 0, origin[1], 0, pixel_height]

    # perform rotation
    if angle:
        rot = math.radians(angle)

        if not rotate_origin:
            if center is None:
                raise ValueError("Cannot rotate about the center without center point")
            else:
                # BM: By default, rotation will be done about origin
                # (upper left). To rotate about center we have to
                # calculate a new origin by determining upper-left
                # point as though it were rotated around center.
                gt[0] = center.east + (origin[0] - center.east) * math.cos(rot) \
                    + (origin[1] - center.north) * math.sin(rot)
                gt[3] = center.north - (origin[0] - center.east) * math.sin(rot) \
                    + (origin[1] - center.north) * math.cos(rot)

        gt[1] = pixel_width * math.cos(rot)
        gt[2] = pixel_width * -math.sin(rot)
        gt[4] = pixel_height * math.sin(rot)
        gt[5] = pixel_height * math.cos(rot)
        filename = '{}_rotated_{}.tif'\
                   .format(os.path.splitext(filename)[0], angle)

    rows, cols = data.shape
    driver = gdal.GetDriverByName('GTiff')
    out_raster = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    out_raster.SetGeoTransform(gt)
    out_band = out_raster.GetRasterBand(1)
    out_band.SetNoDataValue(np.nan)
    out_band.WriteArray(data)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code)
    out_raster.SetProjection(srs.ExportToWkt())
    out_band.FlushCache()

    # output to ascii format
    ascii_filename = "{}.asc".format(os.path.splitext(filename)[0])
    driver2 = gdal.GetDriverByName('AAIGrid')
    driver2.CreateCopy(ascii_filename, out_raster)

    return filename, ascii_filename


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
                   depths=None, angle=None, rotate_origin=False,
                   list_depths=False):
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

    if list_depths:
        with np.printoptions(precision=0, suppress=True):
            print(gcz)
            # _logger.info(gcz) # Need to fix logging
        return

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

    def _nearest(array, value):
        """Get index for nearest element to value in an array.

        Args:
            array (np.ndarray): Array to get index for.
            value (float): Value to search for.

        Returns:
            int: Index of element closest to value.
        """
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array)
                or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
            return idx - 1
        else:
            return idx

    if depths:
        indicies = {_nearest(gcz, d) for d in depths}
    else:
        indicies = range(len(gcz))

    _logger.info("Depth indicies = {}".format(indicies))

    for di in indicies:
        output_file = 'DepthSlice{:.0f}m.tif'.format(gcz[di])
        output_file = os.path.join(out_dir, output_file)
        # define interpolation function (interpolate in log10 measure-space)
        # See https://docs.scipy.org/doc/scipy-0.16.0/reference/interpolate.html
        interpfunc = RegularGridInterpolator((gce, gcn), np.log10(resgrid_nopad[:, :, di].T))
        # evaluate on the regular grid points, which to be output into geogrid formatted files
        newgridres = 10 ** interpfunc(np.vstack(
            [target_gridx.flatten(), target_gridy.flatten()]).T).reshape(target_gridx.shape)
        _logger.info("New interpolated resistivity grid shape at index {}: {} "
                     .format(di, newgridres.shape))

        # This original image may start from the lower left corner, if
        # so must be flipped.
        # resis_data_flip = resis_data[::-1]

        array2geotiff_writer(output_file, origin, pixel_width, pixel_height, newgridres[::-1],
                             epsg_code=epsg_code, angle=angle, center=center,
                             rotate_origin=rotate_origin)

    return output_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('modem_data', help="ModEM data file")
    parser.add_argument('modem_model', help="ModEM model file")
    parser.add_argument('out_dir', help="output directory")
    parser.add_argument('--list-depths', action='store_true', default=False,
                        help='list depth of every slice in model then exit')
    parser.add_argument('--xpad', type=int, help='xpad value', default=6)
    parser.add_argument('--ypad', type=int, help='ypad value', default=6)
    parser.add_argument('--zpad', type=int, help='ypad value', default=10)
    parser.add_argument('--grid', type=int, help="pixel size in meters", default=7500)
    parser.add_argument('--epsg', type=int, help="EPSG code for CRS of the model")
    parser.add_argument('--lat', type=float, help="grid center latitude in degrees")
    parser.add_argument('--lon', type=float, help="grid center longitude in degrees")
    parser.add_argument('--depths', type=int, nargs='*', help="depths for slices to convert (in "
                        "meters) eg., '--di 3 22 105 782'")
    parser.add_argument('--angle', type=float, help="angle in degrees to rotate image by")
    parser.add_argument('--rotate-origin', action='store_true', default=False,
                        help='rotate around the original origin (upper left corner), '
                             'otherwise image will be rotated about center')
    args = parser.parse_args()

    create_geogrid(args.modem_data, args.modem_model, args.out_dir, xpad=args.xpad, ypad=args.ypad,
                   grid_size=args.grid, center_lat=args.lat, center_lon=args.lon,
                   depths=args.depths, epsg_code=args.epsg, angle=args.angle,
                   rotate_origin=args.rotate_origin, list_depths=args.list_depths)
