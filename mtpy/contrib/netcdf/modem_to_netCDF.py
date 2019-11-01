#!/usr/bin/env python
"""
This code works on MTpy ModEM format output and generate a netCDF file.
The ModEM format generates an irregular grid in local units (east, north, depth).
Here we converted that into a regular grid in global CRS (in lon, lat, depth).
The regular spacing we used here is the statistical mode of the spacing our input data.
The interpolation is done in local grid before reprojecting.

@author Umma Zannat, GA, 2019
"""
import numpy as np
from scipy import stats
from scipy.interpolate import interp2d
from pyproj import Proj, transform

from mtpy.modeling.modem import Model, Data
from mtpy.utils import gis_tools
from mtpy.contrib.netcdf import nc

import argparse

def mid_point(arr):
    """
    Mid point of a one-dimensional array `arr`.
    Shape of `arr` should be just a tuple with one value in it.
    """
    [shape] = arr.shape
    mid = int(shape / 2)
    return arr[mid]


def uniform_interior_grid(arr, spacing, mid):
    """
    Make a regular grid in lat/lon space.
    """
    end_points = sorted((arr[0], arr[-1]))
    units = int((end_points[0] - mid) / spacing), int((end_points[1] - mid) / spacing) + 1

    candidates = [mid + i * spacing for i in range(units[0], units[1])]
    return np.array([x for x in candidates if end_points[0] < x and x < end_points[1]])


def median_spacing(arr):
    """
    The spacing that occurs the maximum time.
    """
    return stats.mode(arr[1:] - arr[:-1]).mode.item()

def lon_lat_grid_spacing(center, width, height, to_wgs84):
    """
    Returns center longitude and latitude, and spacing for longitude and latitude.
    """
    center_x, center_y = center.east.item(), center.north.item()
    center_lon, center_lat = to_wgs84(center_x, center_y)
    shifted_lon, _ = to_wgs84(center_x + width, center_y)
    _, shifted_lat = to_wgs84(center_x, center_y + height)

    return center_lon, center_lat, shifted_lon - center_lon, shifted_lat - center_lat


def interpolated_layer(x, y, layer):
    """
    Create the interpolation function for each layer.
    see https://docs.scipy.org/doc/scipy-0.16.0/reference/interpolate.html
    """
    return interp2d(x, y, layer)  # bounds_error=True


def converter(in_proj, out_proj):
    """
    Transfrom coordinates from one epsg to another.
    """
    def result(x, y):
        return transform(in_proj, out_proj, x, y)

    return result


def interpolate(resistivity_dict, source_proj, grid_proj, center, east_spacing, north_spacing):
    """
    Interpolate resistivity data to a regular grid.
    """

    to_grid = converter(source_proj, grid_proj)
    from_grid = converter(grid_proj, source_proj)

    center_lon, center_lat, width, height = lon_lat_grid_spacing(center, east_spacing, north_spacing, to_grid)

    lon_list = [to_grid(x, y)[0]
                for x in resistivity_dict['x']
                for y in resistivity_dict['y']]

    lat_list = [to_grid(x, y)[1]
                for x in resistivity_dict['x']
                for y in resistivity_dict['y']]

    interpolation_funcs = [interpolated_layer(resistivity_dict['x'],
                                              resistivity_dict['y'],
                                              resistivity_dict['resistivity'][z_index, :, :])
                           for z_index in range(resistivity_dict['z'].shape[0])]

    result = {
        'longitude': uniform_interior_grid(sorted(lon_list), width, center_lon),
        'latitude': uniform_interior_grid(sorted(lat_list), height, center_lat),
        'depth': resistivity_dict['z']}

    result['resistivity'] = np.zeros(tuple(result[key].shape[0]
                                           for key in ['depth', 'latitude', 'longitude']))

    def uniform_layer(interp_func, latitudes, longitudes):
        """
        Calculate the interpolated values for the layer.
        """
        lats, lons = latitudes.shape[0], longitudes.shape[0]

        ll_result = np.zeros((lats, lons))
        for j in range(lons):
            for i in range(lats):
                lon, lat = longitudes[j], latitudes[i]
                x, y = from_grid(lon, lat)

                ll_result[i, j] = interp_func(x, y)

        return ll_result

    for z_index in range(result['depth'].shape[0]):
        result['resistivity'][z_index, :, :] = uniform_layer(interpolation_funcs[z_index],
                                                             result['latitude'], result['longitude'])

    return result

def main(data_file, model_file, output_file, source_proj=None):
    """
    Generate an output netcdf file from data_file and model_file
    :param data_file: modem.dat
    :param model_file: modem.rho
    :param output_file: output.nc
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

    grid_proj = Proj(init='epsg:4326') # output grid Coordinate systems: 4326, 4283, 3112
    grid_proj = Proj(init='epsg:4283') # output grid Coordinate system 4326, 4283, 3112
    grid_proj = Proj(init='epsg:3112') # output grid Coordinate system 4326, 4283, 3112
    result = interpolate(resistivity_data, source_proj, grid_proj, center,
                         median_spacing(model.grid_east), median_spacing(model.grid_north))

    nc.write_resistivity_grid(output_file, grid_proj,
                              result['latitude'], result['longitude'], result['depth'],
                              result['resistivity'], z_label='depth')

#####################################################################################################################
# cd /e/Githubz/mtpy/mtpy/contrib
# python netcdf/modem_to_netCDF.py ../../examples/model_files/ModEM_2/Modular_MPI_NLCG_004.dat ../../examples/model_files/ModEM_2/Modular_MPI_NLCG_004.rho
#####################################################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('modem_data', help="ModEM data file")
    parser.add_argument('modem_model', help="ModEM model file")
    parser.add_argument('--epsg', help="EPSG code for source CRS", type=int)
    parser.add_argument('--output-file', default="output.nc", help="Name of output NetCDF file")
    args = parser.parse_args()

    main(args.modem_data, args.modem_model, args.output_file, args.epsg)
