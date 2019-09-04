#! /usr/bin/env python
"""
Description:
    Convert input MODEM data  and resistivity rho files into a georeferenced raster/grid format,
    such as geotiff format, which can be visualized by GIS software.

CreationDate:   1/05/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     1/05/2019   FZ
    LastUpdate:     dd/mm/yyyy
"""

import os,sys
import argparse
from pyproj import Proj
import gdal, osr
import numpy as np
from mtpy.modeling.modem import Model, Data
from mtpy.utils import gis_tools
from mtpy.contrib.netcdf import nc
import mtpy.contrib.netcdf.modem_to_netCDF as modem2nc


def array2geotiff_writer(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()



def test_array2geotiff_writer(newRasterfn):
    #rasterOrigin = (-123.25745,45.43013)
    rasterOrigin = (149.298, -34.974)
    pixelWidth = 0.01
    pixelHeight = 0.01
    # Define an array: 0 = black 1=bright
    array = np.array([[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
                      [ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])


    random= np.random.rand(array.shape[0],array.shape[1])

    array2= 100*array + 10*random
    print (array2)
    #reversed_arr = array2[::-1]  # reverse the array to flip image upside down (water surface mirrored shadow)
    array2geotiff_writer(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array2)  # convert array to raster


def modem2geotiff(data_file, model_file, output_file, source_proj=None):
    """
    Generate an output geotiff file from a modems.dat file and related modems.rho model file
    :param data_file: modem.dat
    :param model_file: modem.rho
    :param output_file: output.tif
    :param source_proj: None by defult. The UTM zone infered from the input non-uniform grid parameters
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
        #source_proj = Proj('+proj=utm +zone=%d +%s +datum=%s' % (zone_number, 'north' if is_northern else 'south', 'WGS84'))

        epsg_code = gis_tools.get_epsg(center.lat.item(), center.lon.item())
        print("Input data epsg code is infered as ", epsg_code)
    else:
        epsg_code = source_proj  # integer

    source_proj = Proj(init='epsg:' + str(epsg_code))

    resistivity_data = {
        'x': center.east.item() + (model.grid_east[1:] + model.grid_east[:-1])/2,
        'y': center.north.item() + (model.grid_north[1:] + model.grid_north[:-1])/2,
        'z': (model.grid_z[1:] + model.grid_z[:-1])/2,
        'resistivity': np.transpose(model.res_model, axes=(2, 0, 1))
    }

    #grid_proj = Proj(init='epsg:4326') # output grid Coordinate systems: 4326 WGS84
    grid_proj = Proj(init='epsg:4283') # output grid Coordinate system , 4283 https://spatialreference.org/ref/epsg/gda94/
    # grid_proj = Proj(init='epsg:3112') # output grid Coordinate system 4326, 4283, 3112
    result = modem2nc.interpolate(resistivity_data, source_proj, grid_proj, center,
                         modem2nc.median_spacing(model.grid_east), modem2nc.median_spacing(model.grid_north))


    print("result['latitude'] ==", result['latitude'])
    print("result['longitude'] ==", result['longitude'])
    print("result['depth'] ==", result['depth'])

    origin=(result['longitude'][0],result['latitude'][0])
    pixel_width = result['longitude'][1] - result['longitude'][0]
    pixel_height = result['latitude'][1] - result['latitude'][0]

    # write the depth_index
    depth_index=1
    resis_data = result['resistivity'][depth_index,:,:]
    #resis_data_vflip = resis_data[::-1]  # flipped upside down 

    array2geotiff_writer(output_file,origin,pixel_width,pixel_height,resis_data)

    return output_file


#####################################################################################################################
# Section for quick test run of this script
# cd /e/Githubz/mtpy
# python mtpy/utils/convert_modem_data_to_geogrid.py examples/model_files/ModEM_2/Modular_MPI_NLCG_004.dat examples/model_files/ModEM_2/Modular_MPI_NLCG_004.rho
#####################################################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('modem_data', help="ModEM data file")
    parser.add_argument('modem_model', help="ModEM model file")
    parser.add_argument('--epsg', help="EPSG code for the modem source data CRS", type=int)
    parser.add_argument('--output-file', default="output.tif", help="Name of output file")
    args = parser.parse_args()

    modem2geotiff(args.modem_data, args.modem_model, args.output_file, args.epsg)


    # test_array2geotiff_writer("test_geotiff_GDAL_img.tif")
