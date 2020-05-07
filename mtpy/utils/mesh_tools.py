# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 09:35:31 2017

@author: Alison Kirkby

functions to assist with mesh generation

"""

import numpy as np
import mtpy.utils.filehandling as mtfh
from mtpy.utils import gis_tools
import scipy.interpolate as spi



def grid_centre(grid_edges):
    """
    calculate the grid centres from an array that defines grid edges
    :param grid_edges: array containing grid edges
    :returns: grid_centre: centre points of grid
    """
    return np.mean([grid_edges[1:], grid_edges[:-1]], axis=0)
    

def rotate_mesh(grid_east,grid_north,origin,
                rotation_angle,return_centre = False):
    """
    rotate a mesh defined by grid_east and grid_north.
    
    :param grid_east: 1d array defining the edges of the mesh in the east-west direction
    :param grid_north: 1d array defining the edges of the mesh in the north-south direction
    :param origin: real-world position of the (0,0) point in grid_east, grid_north
    :param rotation_angle: angle in degrees to rotate the grid by
    :param return_centre: True/False option to return points on centre of grid instead of grid edges
    
    :return: grid_east, grid_north - 2d arrays describing the east and north coordinates

    """
    x0,y0 = origin
    
    # centre of grid in relative coordinates
    if return_centre:
        gce, gcn = [np.mean([arr[1:], arr[:-1]], axis=0)
                  for arr in [grid_east,grid_north]]
    else:
        gce, gcn = grid_east, grid_north
        
    # coordinates (2d array)
    coords = np.array([arr.flatten() for arr in np.meshgrid(gce,gcn)])


    
    if rotation_angle != 0:
        # create the rotation matrix
        cos_ang = np.cos(np.deg2rad(rotation_angle))
        sin_ang = np.sin(np.deg2rad(rotation_angle))
        rot_matrix = np.array([[cos_ang, sin_ang],
                               [-sin_ang, cos_ang]])

        # rotate the relative grid coordinates
        new_coords = np.array(np.dot(rot_matrix, coords))
    else:
        new_coords = coords
    

    # location of grid centres in real-world coordinates to interpolate elevation onto
    xg = (new_coords[0] + x0).reshape(len(gcn),len(gce))
    yg = (new_coords[1] + y0).reshape(len(gcn),len(gce))
    
    return xg, yg


def interpolate_elevation_to_grid(grid_east, grid_north, epsg=None, utm_zone=None,
                                  surfacefile=None, surface=None, method='linear',
                                  fast=True):
    """
    # Note: this documentation is outdated and seems to be copied from
    #  model.interpolate_elevation2. It needs to be updated. This
    #  funciton does not update a dictionary but returns an array of
    #  elevation data.

    project a surface to the model grid and add resulting elevation data
    to a dictionary called surface_dict. Assumes the surface is in lat/long
    coordinates (wgs84)
    The 'fast' method extracts a subset of the elevation data that falls within the
    mesh-bounds and interpolates them onto mesh nodes. This approach significantly
    speeds up (~ x5) the interpolation procedure.

    **returns**
    nothing returned, but surface data are added to surface_dict under
    the key given by surfacename.

    **inputs**
    choose to provide either surface_file (path to file) or surface (tuple).
    If both are provided then surface tuple takes priority.

    surface elevations are positive up, and relative to sea level.
    surface file format is:

    ncols         3601
    nrows         3601
    xllcorner     -119.00013888889 (longitude of lower left)
    yllcorner     36.999861111111  (latitude of lower left)
    cellsize      0.00027777777777778
    NODATA_value  -9999
    elevation data W --> E
    N
    |
    V
    S

    Alternatively, provide a tuple with:
    (lon,lat,elevation)
    where elevation is a 2D array (shape (ny,nx)) containing elevation
    points (order S -> N, W -> E)
    and lon, lat are either 1D arrays containing list of longitudes and
    latitudes (in the case of a regular grid) or 2D arrays with same shape
    as elevation array containing longitude and latitude of each point.

    other inputs:
    surfacename = name of surface for putting into dictionary
    surface_epsg = epsg number of input surface, default is 4326 for lat/lon(wgs84)
    method = interpolation method. Default is 'nearest', if model grid is
    dense compared to surface points then choose 'linear' or 'cubic'
    """
    # read the surface data in from ascii if surface not provided
    if surfacefile:
        lon, lat, elev = mtfh.read_surface_ascii(surfacefile)
    elif surface:
        lon, lat, elev = surface
    else:
        raise ValueError("'surfacefile' or 'surface' must be provided")

    # if lat/lon provided as a 1D list, convert to a 2d grid of points
    if len(lon.shape) == 1:
        # BM: There seems to be an issue using dense grids (X and Y
        #  become arrays of (N, N)), and get flattened to 1D array
        #  of N^2 in point projection below.
        #  Interpolation then can't be performed because
        #  there's a dimension mismatch between lon/lat and elev.
        #  This issue doesn't happen when using 'fast' method below, so
        #  when not using fast, get a sparse grid instead.
        if fast:
            lon, lat = np.meshgrid(lon, lat)
        else:
            lon, lat = np.meshgrid(lon, lat, sparse=True)
            lat = lat.T

    if len(grid_east.shape) == 1:
        grid_east, grid_north = np.meshgrid(grid_east, grid_north)

    if(fast):
        buffer = 1  # use a buffer of 1 degree around mesh-bounds
        mlatmin, mlonmin = gis_tools.project_point_utm2ll(grid_east.min(),
                                                          grid_north.min(),
                                                          epsg=epsg,
                                                          utm_zone=utm_zone)

        mlatmax, mlonmax = gis_tools.project_point_utm2ll(grid_east.max(), 
                                                          grid_north.max(),
                                                          epsg=epsg,
                                                          utm_zone=utm_zone)
        subsetIndices = (lon >= mlonmin - buffer) & \
                        (lon <= mlonmax + buffer) & \
                        (lat >= mlatmin - buffer) & \
                        (lat <= mlatmax + buffer)
        lon = lon[subsetIndices]
        lat = lat[subsetIndices]
        elev = elev[subsetIndices]

    # end if
    projected_points = gis_tools.project_point_ll2utm(lat, lon, epsg=epsg,
                                                      utm_zone=utm_zone)
    # elevation in model grid
    # first, get lat,lon points of surface grid
    points = np.vstack([arr.flatten() for arr in [projected_points.easting,
                                                  projected_points.northing]]).T
    # corresponding surface elevation points
    values = elev.flatten()
    # xi, the model grid points to interpolate to
    xi = np.vstack([arr.flatten() for arr in [grid_east, grid_north]]).T
    # elevation on the centre of the grid nodes
    elev_mg = spi.griddata(points, values, xi, 
                           method=method).reshape(grid_north.shape)

    return elev_mg


def get_nearest_index(array, value):
    """
    Return the index of the nearest value to the provided value in an array:
    
        inputs:
            array = array or list of values
            value = target value
            
    """
    array = np.array(array)
    
    abs_diff = np.abs(array - value)
    
    return np.where(abs_diff==np.amin(abs_diff))[0][0]
    


def make_log_increasing_array(z1_layer, target_depth, n_layers, 
                              increment_factor=0.9):
    """
    create depth array with log increasing cells, down to target depth,
    inputs are z1_layer thickness, target depth, number of layers (n_layers)
    """        
    
    # make initial guess for maximum cell thickness
    max_cell_thickness = target_depth
    # make initial guess for log_z
    log_z = np.logspace(np.log10(z1_layer), 
                        np.log10(max_cell_thickness),
                        num=n_layers)
    counter = 0
    
    while np.sum(log_z) > target_depth:
        max_cell_thickness *= increment_factor
        log_z = np.logspace(np.log10(z1_layer), 
                            np.log10(max_cell_thickness),
                            num=n_layers) 
        counter += 1
        if counter > 1e6:
            break        

    return log_z


def get_padding_cells(cell_width, max_distance, num_cells, stretch):
    """
    get padding cells, which are exponentially increasing to a given 
    distance.  Make sure that each cell is larger than the one previously.
    
    Arguments
    -------------
    
        **cell_width** : float
                         width of grid cell (m)
                         
        **max_distance** : float
                           maximum distance the grid will extend (m)
                           
        **num_cells** : int
                        number of padding cells
                        
        **stretch** : float
                      base geometric factor
                        
    Returns
    ----------------
    
        **padding** : np.ndarray
                      array of padding cells for one side
    
    """

    # compute scaling factor
    scaling = ((max_distance)/(cell_width*stretch))**(1./(num_cells-1)) 
    
    # make padding cell
    padding = np.zeros(num_cells)
    for ii in range(num_cells):
        # calculate the cell width for an exponential increase
        exp_pad = np.round((cell_width*stretch)*scaling**ii, -2)
        
        # calculate the cell width for a geometric increase by 1.2
        mult_pad = np.round((cell_width*stretch)*((1-stretch**(ii+1))/(1-stretch)), -2)
        
        # take the maximum width for padding
        padding[ii] = max([exp_pad, mult_pad])

    return padding


def get_padding_from_stretch(cell_width, pad_stretch, num_cells):
    """
    get padding cells using pad stretch factor
    
    """
    nodes = np.around(cell_width * (np.ones(num_cells)*pad_stretch)**np.arange(num_cells),-2)
    
    return np.array([nodes[:i].sum() for i in range(1,len(nodes)+1)])
    
    

def get_padding_cells2(cell_width, core_max, max_distance, num_cells):
    """
    get padding cells, which are exponentially increasing to a given 
    distance.  Make sure that each cell is larger than the one previously.
    """
    # check max distance is large enough to accommodate padding
    max_distance = max(cell_width*num_cells, max_distance)

    cells = np.around(np.logspace(np.log10(core_max),np.log10(max_distance),num_cells), -2)
    cells -= core_max
        
    return cells
    
    
def get_station_buffer(grid_east,grid_north,station_east,station_north,buf=10e3):
    """
    get cells within a specified distance (buf) of the stations
    returns a 2D boolean (True/False) array
    
    """
    first = True
    for xs,ys in np.vstack([station_east,station_north]).T:
        xgrid,ygrid = np.meshgrid(grid_east,grid_north)
        station_distance = ((xs - xgrid)**2 + (ys - ygrid)**2)**0.5
        if first:
            where = station_distance < buf
            first = False
        else:
            where = np.any([where,station_distance < buf],axis=0)
            
    return where
    
    
