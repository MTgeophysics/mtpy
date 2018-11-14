import os
from netCDF4 import Dataset

def create_dataset(filename, overwrite=True):
    if os.path.exists(filename):
        if overwrite:
            os.remove(filename)
        else:
            raise ValueError("file {} already exists".format(filename))

    return Dataset(filename, 'w', format='NETCDF4')

def set_grid_mapping_attrs_geographic(crs_var, spatial_ref):
    crs_var.spatial_ref = spatial_ref.ExportToWkt()
    crs_var.crs_wkt = spatial_ref.ExportToWkt()
    crs_var.inverse_flattening = spatial_ref.GetInvFlattening()
    crs_var.long_name = spatial_ref.GetAttrValue('GEOGCS')
    crs_var.semi_major_axis = spatial_ref.GetSemiMajor()
    crs_var.grid_mapping_name = 'latitude_longitude'
    crs_var.units = spatial_ref.GetAttrValue('UNIT')

def write_resistivity_grid(output_file, spatial_ref,
                           latitude, longitude, elevation, resistivity_data,
                           **kwargs):
    """ resistivity_data in (elevation, latitude, longitude) grid. """

    with create_dataset(output_file) as dataset:
        dataset.description = 'Resistivity Model'

        z_label = kwargs.get('z_label', 'elevation')

        # dimensions
        dataset.createDimension('latitude', latitude.shape[0])
        dataset.createDimension('longitude', longitude.shape[0])
        dataset.createDimension(z_label, elevation.shape[0])

        # variables
        x = dataset.createVariable('latitude', 'f4', ('latitude',))
        y = dataset.createVariable('longitude', 'f4', ('longitude',))
        z = dataset.createVariable(z_label, 'f4', (z_label,))

        for var in [x, y]:
            var.units = 'degree'
        z.units = 'm'

        resistivity = dataset.createVariable('resistivity', 'f4', (z_label, 'latitude', 'longitude'))
        resistivity.grid_mapping = 'crs'
        resistivity.long_name = 'resistivity'
        resistivity.units = "ohm-m"

        # populate variables
        x[:] = latitude
        y[:] = longitude
        z[:] = elevation

        resistivity[:, :, :] = resistivity_data

        # attach crs info
        crs_var = dataset.createVariable('crs', 'i4', ())
        set_grid_mapping_attrs_geographic(crs_var, spatial_ref)

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
