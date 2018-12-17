import os
from netCDF4 import Dataset

def create_dataset(filename, overwrite=True):
    if os.path.exists(filename):
        if overwrite:
            os.remove(filename)
        else:
            raise ValueError("file {} already exists".format(filename))

    return Dataset(filename, 'w', format='NETCDF4')

def set_grid_mapping_attrs_geographic(crs_var, epsg_code):
    # for now, we only support epsg 4326
    assert epsg_code == 4326
    wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
    crs_var.spatial_ref = wkt
    crs_var.crs_wkt = wkt
    crs_var.long_name = 'WGS 84'
    crs_var.grid_mapping_name = 'latitude_longitude'
    crs_var.units = 'degree'

def write_resistivity_grid(output_file, epsg_code,
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
        set_grid_mapping_attrs_geographic(crs_var, epsg_code)
