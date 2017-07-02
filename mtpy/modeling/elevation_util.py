"""
Description:
   Utility functions to handle elevation data

Author: fei.zhang@ga.gov.au

Date:
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spi

import mtpy.utils.latlon_utm_conversion as utm2ll


# ==============================================================================
# Add in elevation to the model
# ==============================================================================

def read_surface_ascii(ascii_fn):
    """
    read in surface which is ascii format ()
    unlike original function, returns list of lat, long and elevation (no projections)

    The ascii format is assumed to be:
    ncols        2743
    nrows        2019
    xllcorner    111.791666666667 (lon of lower left)
    yllcorner    -45.341666666667 (lat of lower left)
    cellsize     0.016666666667
    NODATA_value  -9999
    elevation data W --> E
    N
    |
    V
    S
    """
    dfid = file(ascii_fn, 'r')
    d_dict = {}
    skiprows = 0
    for ii in range(6):
        dline = dfid.readline()
        dline = dline.strip().split()
        key = dline[0].strip().lower()
        value = float(dline[1].strip())
        d_dict[key] = value
        # check if key is an integer
        try:
            int(key)
        except:
            skiprows += 1
    dfid.close()

    x0 = d_dict['xllcorner']
    y0 = d_dict['yllcorner']
    nx = int(d_dict['ncols'])
    ny = int(d_dict['nrows'])
    cs = d_dict['cellsize']

    elevation = np.loadtxt(ascii_fn, skiprows=skiprows)[::-1]

    # create lat and lon arrays from the dem fle
    lon = np.arange(x0, x0 + cs * (nx), cs)
    lat = np.arange(y0, y0 + cs * (ny), cs)
    lon = np.linspace(x0, x0 + cs * (nx - 1), nx)
    lat = np.linspace(y0, y0 + cs * (ny - 1), ny)

    return lon, lat, elevation
    # return lat, lon, elevation # FZ: switch lat-lon to match with MT's coordinate definition


# --> read in ascii dem file
def read_dem_ascii(ascii_fn, cell_size=500, model_center=(0, 0), rot_90=0, epsg=None):
    """
    read in dem which is ascii format

    The ascii format is assumed to be:
    ncols         3601
    nrows         3601
    xllcorner     -119.00013888889
    yllcorner     36.999861111111
    cellsize      0.00027777777777778
    NODATA_value  -9999
    elevation data W --> E
    N
    |
    V
    S
    """
    dfid = file(ascii_fn, 'r')
    d_dict = {}
    for ii in range(6):
        dline = dfid.readline()
        dline = dline.strip().split()
        key = dline[0].strip().lower()
        value = float(dline[1].strip())
        d_dict[key] = value

    x0 = d_dict['xllcorner']
    y0 = d_dict['yllcorner']
    nx = int(d_dict['ncols'])
    ny = int(d_dict['nrows'])
    cs = d_dict['cellsize']

    # read in the elevation data
    elevation = np.zeros((nx, ny))

    for ii in range(1, int(ny) + 2):
        dline = dfid.readline()
        if len(str(dline)) > 1:
            # needs to be backwards because first line is the furthest north
            # row.
            elevation[
            :, -ii] = np.array(dline.strip().split(' '), dtype='float')
        else:
            break

    # create lat and lon arrays from the dem fle
    lon = np.arange(x0, x0 + cs * (nx), cs)
    lat = np.arange(y0, y0 + cs * (ny), cs)

    # calculate the lower left and uper right corners of the grid in meters
    ll_en = utm2ll.LLtoUTM(23, lat[0], lon[0])
    ur_en = utm2ll.LLtoUTM(23, lat[-1], lon[-1])

    # estimate cell sizes for each dem measurement
    d_east = abs(ll_en[1] - ur_en[1]) / nx
    d_north = abs(ll_en[2] - ur_en[2]) / ny

    # calculate the number of new cells according to the given cell size
    # if the given cell size and cs are similar int could make the value 0,
    # hence the need to make it one if it is 0.
    num_cells = max([1, int(cell_size / np.mean([d_east, d_north]))])

    # make easting and northing arrays in meters corresponding to lat and lon
    east = np.arange(ll_en[1], ur_en[1], d_east)
    north = np.arange(ll_en[2], ur_en[2], d_north)

    # resample the data accordingly
    new_east = east[np.arange(0, east.shape[0], num_cells)]
    new_north = north[np.arange(0, north.shape[0], num_cells)]

    try:
        new_x, new_y = np.meshgrid(np.arange(0, east.shape[0], num_cells),
                                   np.arange(0, north.shape[0], num_cells),
                                   indexing='ij')
    except TypeError:
        new_x, new_y = [arr.T for arr in np.meshgrid(np.arange(0, east.shape[0], num_cells),
                                                     np.arange(0, north.shape[0], num_cells))]
    elevation = elevation[new_x, new_y]

    # estimate the shift of the DEM to relative model coordinates
    shift_east = new_east.mean() - model_center[0]
    shift_north = new_north.mean() - model_center[1]

    # shift the easting and northing arrays accordingly so the DEM and model
    # are collocated.
    new_east = (new_east - new_east.mean()) + shift_east
    new_north = (new_north - new_north.mean()) + shift_north

    # need to rotate cause I think I wrote the dem backwards
    if rot_90 == 1 or rot_90 == 3:
        elevation = np.rot90(elevation, rot_90)
        return new_north, new_east, elevation
    else:
        elevation = np.rot90(elevation, rot_90)

        return new_east, new_north, elevation


def interpolate_elevation(elev_east, elev_north, elevation, model_east,
                          model_north, pad=3):
    """
    interpolate the elevation onto the model grid.

    Arguments:
    ---------------

        *elev_east* : np.ndarray(num_east_nodes)
                      easting grid for elevation model

        *elev_north* : np.ndarray(num_north_nodes)
                      northing grid for elevation model

        *elevation* : np.ndarray(num_east_nodes, num_north_nodes)
                     elevation model assumes x is east, y is north
                     Units are meters

        *model_east* : np.ndarray(num_east_nodes_model)
                     relative easting grid of resistivity model

        *model_north* : np.ndarray(num_north_nodes_model)
                     relative northin grid of resistivity model

        *pad* : int
                number of cells to repeat elevation model by.  So for pad=3,
                then the interpolated elevation model onto the resistivity
                model grid will have the outer 3 cells will be repeats of
                the adjacent cell.  This is to extend the elevation model
                to the resistivity model cause most elevation models will
                not cover the entire area.

    Returns:
    --------------

        *interp_elev* : np.ndarray(num_north_nodes_model, num_east_nodes_model)
                        the elevation model interpolated onto the resistivity
                        model grid.

    """
    # need to line up the elevation with the model
    grid_east, grid_north = np.broadcast_arrays(elev_east[:, None],
                                                elev_north[None, :])
    # interpolate onto the model grid
    interp_elev = spi.griddata((grid_east.ravel(), grid_north.ravel()),
                               elevation.ravel(),
                               (model_east[:, None],
                                model_north[None, :]),
                               method='linear',
                               fill_value=elevation.mean())

    interp_elev[0:pad, pad:-pad] = interp_elev[pad, pad:-pad]
    interp_elev[-pad:, pad:-pad] = interp_elev[-pad - 1, pad:-pad]
    interp_elev[:, 0:pad] = interp_elev[:, pad].repeat(pad).reshape(
        interp_elev[:, 0:pad].shape)
    interp_elev[:, -pad:] = interp_elev[:, -pad - 1].repeat(pad).reshape(
        interp_elev[:, -pad:].shape)

    # transpose the modeled elevation to align with x=N, y=E
    interp_elev = interp_elev.T

    return interp_elev


def make_elevation_model(interp_elev, model_nodes_z, elevation_cell=30,
                         pad=3, res_air=1e12, fill_res=100, res_sea=0.3):
    """
    Take the elevation data of the interpolated elevation model and map that
    onto the resistivity model by adding elevation cells to the existing model.

    ..Note: that if there are large elevation gains, the elevation cell size
            might need to be increased.

    Arguments:
    -------------
        *interp_elev* : np.ndarray(num_nodes_north, num_nodes_east)
                        elevation model that has been interpolated onto the
                        resistivity model grid. Units are in meters.

        *model_nodes_z* : np.ndarray(num_z_nodes_of_model)
                          vertical nodes of the resistivity model without
                          topography.  Note these are the nodes given in
                          relative thickness, not the grid, which is total
                          depth.  Units are meters.

        *elevation_cell* : float
                           height of elevation cells to be added on.  These
                           are assumed to be the same at all elevations.
                           Units are in meters

        *pad* : int
                number of cells to look for maximum and minimum elevation.
                So if you only want elevations within the survey area,
                set pad equal to the number of padding cells of the
                resistivity model grid.

        *res_air* : float
                    resistivity of air.  Default is 1E12 Ohm-m

        *fill_res* : float
                     resistivity value of subsurface in Ohm-m.

    Returns:
    -------------
        *elevation_model* : np.ndarray(num_north_nodes, num_east_nodes,
                                       num_elev_nodes+num_z_nodes)
                         Model grid with elevation mapped onto it.
                         Where anything above the surface will be given the
                         value of res_air, everything else will be fill_res

        *new_nodes_z* : np.ndarray(num_z_nodes+num_elev_nodes)
                        a new array of vertical nodes, where any nodes smaller
                        than elevation_cell will be set to elevation_cell.
                        This can be input into a modem.Model object to
                        rewrite the model file.

    """

    # calculate the max elevation within survey area
    elev_max = interp_elev[pad:-pad, pad:-pad].max()

    # need to set sea level to 0 elevation
    elev_min = max([0, interp_elev[pad:-pad, pad:-pad].min()])

    # scale the interpolated elevations to fit within elev_max, elev_min
    interp_elev[np.where(interp_elev > elev_max)] = elev_max
    # interp_elev[np.where(interp_elev < elev_min)] = elev_min

    # calculate the number of elevation cells needed
    num_elev_cells = int((elev_max - elev_min) / elevation_cell)
    print 'Number of elevation cells: {0}'.format(num_elev_cells)

    # find sea level if it is there
    if elev_min < 0:
        sea_level_index = num_elev_cells - \
                          abs(int((elev_min) / elevation_cell)) - 1
    else:
        sea_level_index = num_elev_cells - 1

    print 'Sea level index is {0}'.format(sea_level_index)

    # make an array of just the elevation for the model
    # north is first index, east is second, vertical is third
    elevation_model = np.ones((interp_elev.shape[0],
                               interp_elev.shape[1],
                               num_elev_cells + model_nodes_z.shape[0]))

    elevation_model[:, :, :] = fill_res

    # fill in elevation model with air values.  Remeber Z is positive down, so
    # the top of the model is the highest point and index 0 is highest
    # elevation
    for nn in range(interp_elev.shape[0]):
        for ee in range(interp_elev.shape[1]):
            # need to test for ocean
            if interp_elev[nn, ee] < 0:
                # fill in from bottom to sea level, then rest with air
                elevation_model[nn, ee, 0:sea_level_index] = res_air
                dz = sea_level_index + \
                     abs(int((interp_elev[nn, ee]) / elevation_cell)) + 1
                elevation_model[nn, ee, sea_level_index:dz] = res_sea
            else:
                dz = int((elev_max - interp_elev[nn, ee]) / elevation_cell)
                elevation_model[nn, ee, 0:dz] = res_air

    # make new z nodes array
    new_nodes_z = np.append(np.repeat(elevation_cell, num_elev_cells),
                            model_nodes_z)

    new_nodes_z[np.where(new_nodes_z < elevation_cell)] = elevation_cell

    return elevation_model, new_nodes_z

