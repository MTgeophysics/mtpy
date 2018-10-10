#!/usr/bin/env python

import os
from os.path import join, abspath

import numpy as np
from netCDF4 import Dataset
from osgeo import osr
from scipy.interpolate import interp2d
from pyproj import Proj, transform

from mtpy.modeling.modem import Model, Data


# Create a cell-centred NetCDF file
def create_dataset(filename):
    if os.path.exists(filename):
        print filename, 'already exists, removing'
        os.remove(filename)

    return Dataset(filename, 'w', format='NETCDF4')


# get the CRS
def get_spatial_ref(center):
    result = osr.SpatialReference()
    result.SetWellKnownGeogCS('WGS84')
    is_northern = center.lat.item() >= 0
    zone_number = int(1 + (center.lon.item() + 180.0) / 6.0)
    result.SetUTM(zone_number, is_northern)
    return result


def global_spatial_ref():
    result = osr.SpatialReference()
    result.ImportFromEPSG(4326)
    return result


def set_grid_mapping_attrs_projected(crs_var, spatial_ref):
    # currently unused
    crs_var.spatial_ref = spatial_ref.ExportToWkt()
    crs_var.crs_wkt = spatial_ref.ExportToWkt()
    crs_var.inverse_flattening = spatial_ref.GetInvFlattening()
    crs_var.long_name = spatial_ref.GetAttrValue('PROJCS')
    crs_var.semi_major_axis = spatial_ref.GetSemiMajor()
    crs_var.semi_minor_axis = spatial_ref.GetSemiMinor()
    crs_var.grid_mapping_name = spatial_ref.GetAttrValue('PROJECTION')
    crs_var.longitude_of_central_meridian = spatial_ref.GetProjParm(osr.SRS_PP_CENTRAL_MERIDIAN)
    crs_var.false_easting = spatial_ref.GetProjParm(osr.SRS_PP_FALSE_EASTING)
    crs_var.false_northing = spatial_ref.GetProjParm(osr.SRS_PP_FALSE_NORTHING)
    crs_var.latitude_of_projection_origin = spatial_ref.GetProjParm(osr.SRS_PP_LATITUDE_OF_ORIGIN)
    crs_var.scale_factor_at_central_meridian = spatial_ref.GetProjParm(osr.SRS_PP_SCALE_FACTOR)


def set_grid_mapping_attrs_geographic(crs_var, spatial_ref):
    crs_var.spatial_ref = spatial_ref.ExportToWkt()
    crs_var.crs_wkt = spatial_ref.ExportToWkt()
    crs_var.inverse_flattening = spatial_ref.GetInvFlattening()
    crs_var.long_name = spatial_ref.GetAttrValue('GEOGCS')
    crs_var.semi_major_axis = spatial_ref.GetSemiMajor()
    crs_var.grid_mapping_name = 'latitude_longitude'
    crs_var.units = spatial_ref.GetAttrValue('UNIT')


def write_to_file(filename, spatial_ref, resistivity_data):
    with create_dataset(filename) as dataset:
        dataset.description = 'Resistivity Model'

        # dimensions
        dataset.createDimension('latitude', resistivity_data['latitude'].shape[0])
        dataset.createDimension('longitude', resistivity_data['longitude'].shape[0])
        dataset.createDimension('depth', resistivity_data['depth'].shape[0])

        # variables
        x = dataset.createVariable('latitude', 'f4', ('latitude',))
        y = dataset.createVariable('longitude', 'f4', ('longitude',))
        z = dataset.createVariable('depth', 'f4', ('depth',))

        for var in [x, y]:
            var.units = 'degree'
        z.units = 'km'

        resistivity = dataset.createVariable('resistivity', 'f4', ('depth', 'latitude', 'longitude'))
        resistivity.grid_mapping = 'crs'
        resistivity.long_name = 'resistivity'
        resistivity.units = "ohm-m"

        # populate variables
        x[:] = resistivity_data['latitude']
        y[:] = resistivity_data['longitude']
        z[:] = resistivity_data['depth']

        resistivity[:, :, :] = resistivity_data['resistivity']

        # attach crs info
        crs_var = dataset.createVariable('crs', 'i4', ())
        set_grid_mapping_attrs_geographic(crs_var, spatial_ref)


def mid_point(arr):
    # shape should be just a tuple with one value in it
    [shape] = arr.shape
    mid = int(shape / 2)
    return arr[mid]


def uniform_interior_grid(arr, spacing, mid):
    end_points = sorted((arr[0], arr[-1]))
    units = int((end_points[0] - mid) / spacing), int((end_points[1] - mid) / spacing) + 1

    candidates = [mid + i * spacing for i in range(units[0], units[1])]
    return np.array([x for x in candidates if end_points[0] < x and x < end_points[1]])


def median_spacing(arr):
    return np.median(arr[1:] - arr[:-1])


def lon_lat_grid_spacing(center, width, height, to_wgs84):
    center_x, center_y = center.east.item(), center.north.item()
    center_lon, center_lat = to_wgs84(center_x, center_y)
    shifted_lon, _ = to_wgs84(center_x + width, center_y)
    _, shifted_lat = to_wgs84(center_x, center_y + height)

    return center_lon, center_lat, shifted_lon - center_lon, shifted_lat - center_lat


def interpolated_layer(x, y, layer):
    return interp2d(x, y, layer)  # bounds_error=True


def converter(in_spatial_ref, out_spatial_ref):
    in_proj = Proj(in_spatial_ref.ExportToProj4())
    out_proj = Proj(out_spatial_ref.ExportToProj4())

    def result(x, y):
        return transform(in_proj, out_proj, x, y)

    return result


def main():
    # Define Data and Model Paths
    MT_PATH = abspath(join(__file__, '../../../..'))

    data = Data()
    data.read_data_file(data_fn=join(MT_PATH, 'examples/model_files/ModEM/ModEM_Data.dat'))

    # create a model object using the data object and read in model data
    model = Model(data_obj=data)
    model.read_model_file(model_fn=join(MT_PATH, 'examples/model_files/ModEM/ModEM_Model_File.rho'))

    center = data.center_point

    resistivity_data = {
        'x': center.east.item() + (model.grid_east[1:] + model.grid_east[:-1])/2,
        'y': center.north.item() + (model.grid_north[1:] + model.grid_north[:-1])/2,
        'z': (model.grid_z[1:] + model.grid_z[:-1])/2,
        'resistivity': np.transpose(model.res_model, axes=(2, 0, 1))
    }

    spatial_ref = get_spatial_ref(center)
    global_ref = global_spatial_ref()

    to_wgs84 = converter(spatial_ref, global_ref)
    from_wgs84 = converter(global_ref, spatial_ref)

    center_lon, center_lat, width, height = lon_lat_grid_spacing(center, median_spacing(model.grid_east),
                                                                 median_spacing(model.grid_north), to_wgs84)

    lon_list = [to_wgs84(x, y)[0]
                for x in resistivity_data['x']
                for y in resistivity_data['y']]

    lat_list = [to_wgs84(x, y)[1]
                for x in resistivity_data['x']
                for y in resistivity_data['y']]

    interpolation_funcs = [interpolated_layer(resistivity_data['x'],
                                              resistivity_data['y'],
                                              resistivity_data['resistivity'][z_index, :, :])
                           for z_index in range(resistivity_data['z'].shape[0])]

    result = {
        'longitude': uniform_interior_grid(sorted(lon_list), width, center_lon),
        'latitude': uniform_interior_grid(sorted(lat_list), height, center_lat),
        'depth': resistivity_data['z']
    }

    result['resistivity'] = np.zeros(tuple(result[key].shape[0]
                                           for key in ['depth', 'latitude', 'longitude']))

    def uniform_layer(interp_func, latitudes, longitudes):
        lats, lons = latitudes.shape[0], longitudes.shape[0]

        result = np.zeros((lats, lons))
        for j in range(lons):
            for i in range(lats):
                lon, lat = longitudes[j], latitudes[i]
                x, y = from_wgs84(lon, lat)

                result[i, j] = interp_func(x, y)

        return result

    for z_index in range(result['depth'].shape[0]):
        print 'layer #', z_index + 1
        result['resistivity'][z_index, :, :] = uniform_layer(interpolation_funcs[z_index],
                                                             result['latitude'], result['longitude'])

    write_to_file('wgs84.nc', global_ref, result)


if __name__ == '__main__':
    main()
