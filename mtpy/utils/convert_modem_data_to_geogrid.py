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

import os, sys
import argparse
from pyproj import Proj
import gdal, osr
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mtpy.modeling.modem import Model, Data
from mtpy.utils import gis_tools
from mtpy.contrib.netcdf import nc
from mtpy.utils.mtpylog import MtPyLog
import mtpy.contrib.netcdf.modem_to_netCDF as modem2nc

_logger = MtPyLog.get_mypy_logger(__name__)


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


def modem2geotiff_notused(data_file, model_file, output_file, source_proj=None):
    """
    First version code to generate an output geotiff file from a modems.rho model file and dat file
    superseded by the new fun create_geogrid()
    :param data_file: modem.dat
    :param model_file: modem.rho
    :param output_file: output.tif
    :param source_proj: None by default. The UTM zone infered from the input non-uniform grid parameters
    :return:
    """
    # Define Data and Model Paths
    data = Data()
    data.read_data_file(data_fn=data_file)

    # create a model object using the data object and read in model data
    model = Model(data_obj=data)
    model.read_model_file(model_fn=model_file)

    center = data.center_point
    if source_proj is None:
        zone_number, is_northern, utm_zone = gis_tools.get_utm_zone(center.lat.item(), center.lon.item())
        # source_proj = Proj('+proj=utm +zone=%d +%s +datum=%s' % (zone_number, 'north' if is_northern else 'south', 'WGS84'))

        epsg_code = gis_tools.get_epsg(center.lat.item(), center.lon.item())
        print("Input data epsg code is inferred as ", epsg_code)
    else:
        epsg_code = source_proj  # integer

    source_proj = Proj(init='epsg:' + str(epsg_code))

    resistivity_data = {
        'x': center.east.item() + (model.grid_east[1:] + model.grid_east[:-1]) / 2,
        'y': center.north.item() + (model.grid_north[1:] + model.grid_north[:-1]) / 2,
        'z': (model.grid_z[1:] + model.grid_z[:-1]) / 2,
        'resistivity': np.transpose(model.res_model, axes=(2, 0, 1))
    }

    # epsgcode= 4326 # 4326 output grid Coordinate systems: 4326 WGS84
    epsgcode = 4283  # 4283 https://spatialreference.org/ref/epsg/gda94/
    grid_proj = Proj(init='epsg:%s' % epsgcode)  # output grid Coordinate system
    # grid_proj = Proj(init='epsg:3112') # output grid Coordinate system 4326, 4283, 3112
    result = modem2nc.interpolate(resistivity_data, source_proj, grid_proj, center,
                                  modem2nc.median_spacing(model.grid_east), modem2nc.median_spacing(model.grid_north))

    print("result['latitude'] ==", result['latitude'])
    print("result['longitude'] ==", result['longitude'])
    print("result['depth'] ==", result['depth'])

    # origin=(result['longitude'][0],result['latitude'][0]) # which corner of the image?
    origin = (result['longitude'][0], result['latitude'][-1])
    pixel_width = result['longitude'][1] - result['longitude'][0]
    pixel_height = result['latitude'][0] - result['latitude'][
        1]  # This should be negative for geotiff with origin at the upper-left corner

    # write the depth_index
    depth_index = 1
    resis_data = result['resistivity'][depth_index, :,
                 :]  # this original image may start from the lower left corner, if so must be flipped.
    resis_data_flip = resis_data[::-1]  # flipped to ensure the image starts from the upper left corner 

    array2geotiff_writer(output_file, origin, pixel_width, pixel_height, resis_data_flip, epsg_code=epsgcode)

    return output_file


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


def get_user_input_params(user_option_file):
    """
    Get the user's input optional parameters from the input text file
    :param user_option_file: path to a file
    :return: a dictionary like  {"source_proj":28353,"depth_index":[0,1,2,10],}
    """

    # Parse the input file to get a disctionary.

    user_dict = {
        "source_proj": 28353,
        "depth_index": [0, 1, 2, 10],
    }
    return user_dict
#####################################################################################################################
# Quick test of this script
# cd /e/Githubz/mtpy
# export PYTHONPATH=/g/data/ha3/fxz547/Githubz/mtpy
# python mtpy/utils/convert_modem_data_to_geogrid.py examples/model_files/ModEM_2/Modular_MPI_NLCG_004.dat examples/model_files/ModEM_2/Modular_MPI_NLCG_004.rho
# python mtpy/utils/convert_modem_data_to_geogrid.py tmp/JinMing_GridData_sample/EFTF_MT_model/EF_NLCG_001.dat  tmp/JinMing_GridData_sample/EFTF_MT_model/EF_NLCG_001.rho
# python mtpy/utils/convert_modem_data_to_geogrid.py /c/Data/JinMing_GridData_sample/JM_model_002/EFTF_NLCG_002.dat /c/Data/JinMing_GridData_sample/JM_model_002/EFTF_NLCG_002.rho --user_option_file path2my.conf
#####################################################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('modem_data', help="ModEM data file")
    parser.add_argument('modem_model', help="ModEM model file")
    parser.add_argument('out_dir', help="output directory")
    parser.add_argument('--user_option_file', help="Uer Optional parameter file path")

    args = parser.parse_args()

    print("The cmdline input parameters are:", args)

    # Before calling the function create_geogrid(), a user should provide the right optional parameters.
    # Otherwise, default parameters will be used, which may not make sense
    if args.user_option_file is not None:
        user_option_dict = get_user_input_params(args.user_option_file)
        # eg {"source_proj":28353,"depth_index":[0,1,2,10],}
    else:
        user_option_dict= {}

    print("User Options:", user_option_dict)

    create_geogrid(args.modem_data, args.modem_model,  args.out_dir, user_options=user_option_dict)  #,depth_index=[0,1,2,10])

    # test_array2geotiff("test_geotiff_GDAL_img.tif", 4326)
