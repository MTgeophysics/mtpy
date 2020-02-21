#! /usr/bin/env python
"""
Description:
    Convert input MODEM data  and resistivity rho files into a georeferenced raster/grid format,
    such as geotiff format, which can be visualized by GIS software.

CreationDate:   1/05/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     01/05/2019  FZ
    LastUpdate:     17/09/2019  FZ fix the geoimage coordinates, upside-down issues
    LastUpdate:     14/02/2020  BM: cleanup, add rotation and indexing
                                    by depth
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

    # Apply rotation by tweaking geotransform. The data remains the
    # same but will appear roated in a viewer e.g. ArcGIS.
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
        filename = '{}_rot{}.tif'\
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


def list_depths(model_file, zpad=None):
    """
    Return a list of available depth slices in the model.

    Args:
        model_file (str): Path to ModEM .rho file.
        zpad (int, optional): Number of padding slices to remove from
            bottom of model. If None, model pad_z value is used.

    Returns:
        list of float: A list of available depth slices.
    """
    model = Model()
    model.read_model_file(model_fn=model_file)

    cz = np.mean([model.grid_z[:-1], model.grid_z[1:]], axis=0)
    zpad = model.pad_z if zpad is None else zpad
    return cz[:-zpad]


def _get_centres(arr):
    return np.mean([arr[:-1], arr[1:]], axis=0)


def _strip_padding(arr, pad, keep_start=False):
    if keep_start:
        return arr[:-pad]
    else:
        return arr[pad:-pad]


def create_geogrid(data_file, model_file, out_dir, x_pad=None, y_pad=None, z_pad=None,
                   x_res=None, y_res=None, center_lat=None, center_lon=None, epsg_code=None,
                   depths=None, angle=None, rotate_origin=False, log_scale=False):
    """Generate an output geotiff file and ASCII grid file.

    Args:
        data_file (str): Path to the ModEM .dat file. Used to get the
            grid center point.
        model_file (str): Path to the ModEM .rho file.
        out_dir (str): Path to directory for storing output data. Will
            be created if it does not exist.
        x_pad (int, optional): Number of east-west padding cells. This
            number of cells will be cropped from the east and west
            sides of the grid. If None, pad_east attribute of model
            will be used.
        y_pad (int, optional): Number of north-south padding cells. This
            number of cells will be cropped from the north and south
            sides of the grid. If None, pad_north attribute of model
            will be used.
        z_pad (int, optional): Number of depth padding cells. This
            number of cells (i.e. slices) will be cropped from the
            bottom of the grid. If None, pad_z attribute of model
            will be used.
        x_res (int, optional): East-west cell size in meters. If None,
            cell_size_east attribute of model will be used.
        y_res (int, optional): North-south cell size in meters. If
            None, cell_size_north of model will be used.
        epsg_code (int, optional): EPSG code of the model CRS. If None,
            is inferred from the grid center point.
        depths (list of int, optional): A list of integers,
            eg, [0, 100, 500], of the depth in metres of the slice to
            retrieve. Will find the closes slice to each depth
            specified. If None, all slices are selected.
        center_lat (float, optional): Grid center latitude in degrees.
            If None, the model's center point will be used.
        center_lon (float, optional): Grid center longitude in degrees.
            If None, the model's center point will be used.
        angle (float, optional): Angle in degrees to rotate image by.
            If None, no rotation is performed.
        rotate_origin (bool, optional): If True, image will be rotated
            around the origin (upper left point). If False, the image
            will be rotated around the center point.
        log_scale (bool, optional): If True, the data will be scaled using log10.
    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    model = Model()
    model.read_model_file(model_fn=model_file)

    data = Data()
    data.read_data_file(data_fn=data_file)
    center = data.center_point
    center_lat = center.lat.item() if center_lat is None else center_lat
    center_lon = center.lon.item() if center_lon is None else center_lon

    _logger.info("Grid center (lat, lon) = ({}, {})".format(center_lat, center_lon))

    if epsg_code is None:
        zone_number, is_northern, utm_zone = gis_tools.get_utm_zone(center_lat, center_lon)
        epsg_code = gis_tools.get_epsg(center_lat, center_lon)
        _logger.info("Input data epsg code has been inferred as {}".format(epsg_code))

    # Get the center point of the model grid cells to use as points
    #  in a resistivity grid.
    ce = _get_centres(model.grid_east)
    cn = _get_centres(model.grid_north)
    cz = _get_centres(model.grid_z)

    # Get X, Y, Z paddings
    x_pad = model.pad_east if x_pad is None else x_pad
    y_pad = model.pad_north if y_pad is None else y_pad
    z_pad = model.pad_z if z_pad is None else z_pad

    # Remove padding cells from the grid
    ce = _strip_padding(ce, x_pad)
    cn = _strip_padding(cn, y_pad)
    cz = _strip_padding(cz, z_pad, keep_start=True)

    # _logger.info("E shape = {}, N shape = {}, Z shape = {}"
    #             .format(ce.shape, cn.shape, cz.shape))
    print("E shape = {}, N shape = {}, Z shape = {}"
          .format(ce.shape, cn.shape, cz.shape))

    # _logger.info("Data center (east, north) = ({}, {})".format(center.east, center.north))
    print("Data center (east, north) = ({}, {})".format(center.east, center.north))

    x_res = model.cell_size_east if x_res is None else x_res
    y_res = model.cell_size_north if y_res is None else y_res

    # BM: The cells have been defined by their center point for making
    #  our grid and interpolating the resistivity model over it. For
    #  display purposes, GDAL expects the origin to be the upper-left
    #  corner of the image. So take the upper left-cell and shift it
    #  half a cell west and north so we get the upper-left corner of
    #  the grid as GDAL origin.
    origin = (ce[0] + center.east - x_res / 2, cn[-1] + center.north + y_res / 2)
    # _logger.info("The Origin (UpperLeft Corner) =".format(origin))
    print("The Origin (UpperLeft Corner) = {}".format(origin))

    target_gridx, target_gridy = np.meshgrid(np.arange(ce[0], ce[-1], x_res),
                                             np.arange(cn[0], cn[-1], y_res))

    resgrid_nopad = model.res_model[y_pad:-y_pad, x_pad:-x_pad, :-z_pad]

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
        indicies = {_nearest(cz, d) for d in depths}
    else:
        indicies = range(len(cz))

    # _logger.info("Depth indicies = {}".format(indicies))
    print("Depth indicies = {}".format(indicies))

    for di in indicies:
        # define interpolation function (interpolate in log10 measure-space)
        interpfunc = RegularGridInterpolator((ce, cn), np.log10(resgrid_nopad[:, :, di].T))
        # evaluate on the regular grid points, which to be output into geogrid formatted files
        newgridres = interpfunc(np.vstack(
            [target_gridx.flatten(), target_gridy.flatten()]).T).reshape(target_gridx.shape)
        if not log_scale:
            newgridres **= 10
            output_file = 'DepthSlice{:.0f}m.tif'.format(cz[di])
        else:
            output_file = 'DepthSlice{:.0f}m_log10.tif'.format(cz[di])
        output_file = os.path.join(out_dir, output_file)

        # _logger.info("New interpolated resistivity grid shape at index {}: {} "
        #              .format(di, newgridres.shape))

        array2geotiff_writer(output_file, origin, x_res, -y_res, newgridres[::-1],
                             epsg_code=epsg_code, angle=angle, center=center,
                             rotate_origin=rotate_origin)

    return output_file


if __name__ == '__main__':
    """
    Example usage:
        python convert_modem_data_to_geogrid.py model.dat model.rho \
                out_directory --grid 800 --depths 100 200 500 \
                --angle 50.0

    will output the slices closest to 100m, 200m and 500m at 800mx80m
    resolution and rotate them 50 degrees about the center.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('modem_data', help="ModEM data file")
    parser.add_argument('modem_model', help="ModEM model file")
    parser.add_argument('out_dir', help="output directory")
    parser.add_argument('--list-depths', action='store_true', default=False,
                        help='list depth of every slice in model then exit')
    parser.add_argument('--xpad', type=int, help='xpad value')
    parser.add_argument('--ypad', type=int, help='ypad value')
    parser.add_argument('--zpad', type=int, help='ypad value')
    parser.add_argument('--xres', type=int, help="cell width in meters")
    parser.add_argument('--yres', type=int, help="cell height in meters")
    parser.add_argument('--epsg', type=int, help="EPSG code for CRS of the model")
    parser.add_argument('--lat', type=float, help="grid center latitude in degrees")
    parser.add_argument('--lon', type=float, help="grid center longitude in degrees")
    parser.add_argument('--depths', type=int, nargs='*', help="depths for slices to convert (in "
                        "meters) eg., '--depths 3 22 105 782'")
    parser.add_argument('--angle', type=float, help="angle in degrees to rotate image by")
    parser.add_argument('--rotate-origin', action='store_true', default=False,
                        help='rotate around the original origin (upper left corner), '
                             'otherwise image will be rotated about center')
    parser.add_argument('--log-scale', action='store_true', default=False,
                        help='scale the data by taking the log10 of data')
    args = parser.parse_args()

    if args.list_depths:
        with np.printoptions(suppress=True):
            print(list_depths(args.modem_model, zpad=args.zpad))
    else:
        create_geogrid(args.modem_data, args.modem_model, args.out_dir, x_pad=args.xpad, y_pad=args.ypad,
                       z_pad=args.zpad, x_res=args.xres, y_res=args.yres, center_lat=args.lat,
                       center_lon=args.lon, depths=args.depths, epsg_code=args.epsg, angle=args.angle,
                       rotate_origin=args.rotate_origin, log_scale=args.log_scale)
