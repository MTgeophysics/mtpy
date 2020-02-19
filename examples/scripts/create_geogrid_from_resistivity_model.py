#! /usr/bin/env python
"""
Description:
    Example python script to create grid formats (geotiff and ascii)
    from ModEM model file.

    See: 'mtpy.utils.convert_modem_data_to_geogrid' for implementation
    and a command line interface.

References:

CreationDate:   6/11/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     6/11/2019   FZ
    LastUpdate:     17/02/2020  BM    Hook up revised code to script

"""
from mtpy.utils.convert_modem_data_to_geogrid import create_geogrid, list_depths

if __name__ == "__main__":
    # Path to model .dat file
    dat_file = "/path/to/model.dat"
    # Path to model .rho file
    rho_file = "/path/to/model.rho"
    # Path to output dir, will be created if it does not exist
    out_dir = "/path/to/out_dir"

    kwargs = {
        # Number of padding cells in respective direction. These will
        #  be removed from the grid. If None, the number of padding
        #  cells will be retrieved from the model itself.
        "x_pad": None,
        "y_pad": None,
        "z_pad": None,
        # Resolution in meters, e.g. x_res = 800 and y_res = 800 will
        #  produce a TIFF with pixel width and height of 800m x 800m.
        "x_res": None,
        "y_res": None,
        # Center lat/lon of model. Leave as 'None' to use .dat file
        #  center coordaintes.
        "center_lat": None,
        "center_lon": None,
        # EPSG code of survey area, e.g. 'EPSG:4326'. Leave as 'None'
        #  to derive from .dat file.
        "epsg_code": None,
        # A list of depths (in meters) to retrieve slices for, e.g.
        #  [100, 200, 500] will return the slices closest to 100m,
        #  200m and 500m. Leave as 'None' to retrieve all slices.
        "depths": [10, 500, 10000],
        # Angle in degrees to rotate TIFF by. 'None' will perform
        #  no rotation.
        "angle": None,
        # Whether or not rotate around the origin. If 'False', the
        #  pivot for rotation is the center point. If 'True, the pivot
        #  is the origin (upper left point).
        "rotate_origin": False,
        # Whether or not to scale the data logarithmically. If True,
        #  the log10 of the data will be taken. If False, the data will
        #  be untouched.
        'log_scale': False
    }

    create_geogrid(dat_file, rho_file, out_dir, **kwargs)

    # Call this function to list available depths:
    # list_depths(rho_file, kwargs[zpad])

