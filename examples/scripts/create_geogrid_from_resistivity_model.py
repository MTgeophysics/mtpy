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
from mtpy.utils.convert_modem_data_to_geogrid import create_geogrid

if __name__ == "__main__":
    # Path to model .dat file
    dat_file = "/path/to/dat_file"
    # Path to model .rho file
    rho_file = "/path/to/rho_file"
    # Path to output dir, will be created if it does not exist
    out_dir = "/path/to/out_dir"

    kwargs = {
        # TODO: @FZ can we get an explanation of these padding params?
        "xpad": 6,
        "ypad": 6,
        "zpad": 10,
        # Resolution in meters, e.g. '800' will produce a TIFF with
        #  pixel width and height of 800m x 800m.
        "grid_size": 800,
        # Center latitude of model. Leave as 'None' to use .dat file
        #  center latitiude.
        "center_lat": None,
        # Center longitude of model. Leave as 'None' to use .dat file
        #  center longitude.
        "center_lon": None,
        # EPSG code of survey area, e.g. 'EPSG:4326'. Leave as 'None'
        #  to derive from .dat file.
        "epsg_code": None,
        # A list of depths (in meters) to retrieve slices for, e.g.
        #  [100, 200, 500] will return the slices closest to 100m,
        #  200m and 500m. Leave as 'None' to retrieve all slices.
        "depths": [0, 100, 500],
        # Angle in degrees to rotate TIFF by. 'None' will perform
        #  no rotation.
        "angle": 45.0,
        # Whether or not rotate around the origin. If 'False', the
        #  pivot for rotation is the center point. If 'True, the pivot
        #  is the origin (upper left point).
        "rotate_origin": False
    }

    create_geogrid(dat_file, rho_file, out_dir, **kwargs)
