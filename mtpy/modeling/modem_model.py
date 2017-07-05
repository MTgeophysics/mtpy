#!/usr/bin/env python
"""
Description:
    Define the Model class.
    This module is refactored from modem.py which had over 5000 lines too big to manage and edit


Author: fei.zhang@ga.gov.au

Date:   2017-06-05
"""

__author__ = 'fei.zhang@ga.gov.au'

# from __future__ import print_function

import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spi

import mtpy.modeling.ws3dinv as ws
import mtpy.utils.gocad as mtgocad
import mtpy.utils.latlon_utm_conversion as utm2ll
from mtpy.modeling.modem_data import Data
import mtpy.modeling.elevation_util as elev_util
from mtpy.utils.mtpylog import MtPyLog
try:
    from evtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')

    print ('Note: if you are using Windows you should build evtk first with'
           'either MinGW or cygwin using the command: \n'
           '    python setup.py build -compiler=mingw32  or \n'
           '    python setup.py build -compiler=cygwin')

logger = MtPyLog().get_mtpy_logger(__name__)


class ModelError(Exception):
    """Raise for ModEM Model class specific exception"""
    pass


# ==============================================================================

class Model(object):
    """
    make and read a FE mesh grid

    The mesh assumes the coordinate system where:
        x == North
        y == East
        z == + down

    All dimensions are in meters.

    ..note:: ModEM assumes all coordinates are relative to North and East, and
             does not accommodate mesh rotations, therefore, here the rotation
             is of the stations, which essentially does the same thing.  You
             will need to rotate you data to align with the 'new' coordinate
             system.

    ==================== ======================================================
    Attributes           Description    
    ==================== ======================================================
    cell_size_east       mesh block width in east direction
                         *default* is 500
    cell_size_north      mesh block width in north direction
                         *default* is 500
    edi_list             list of .edi files to invert for
    grid_east            overall distance of grid nodes in east direction 
    grid_north           overall distance of grid nodes in north direction 
    grid_z               overall distance of grid nodes in z direction 
    model_fn             full path to initial file name
    n_layers             total number of vertical layers in model
    nodes_east           relative distance between nodes in east direction 
    nodes_north          relative distance between nodes in north direction 
    nodes_z              relative distance between nodes in east direction 
    pad_east             number of cells for padding on E and W sides
                         *default* is 7
    pad_north            number of cells for padding on S and N sides
                         *default* is 7
    pad_root_east        padding cells E & W will be pad_root_east**(x)
    pad_root_north       padding cells N & S will be pad_root_north**(x) 
    pad_z                number of cells for padding at bottom
                         *default* is 4
    res_list             list of resistivity values for starting model
    res_model            starting resistivity model
    mesh_rotation_angle  Angle to rotate the grid to. Angle is measured
                         positve clockwise assuming North is 0 and east is 90.
                         *default* is None
    save_path            path to save file to  
    station_fn           full path to station file
    station_locations    location of stations
    title                title in initial file
    z1_layer             first layer thickness
    z_bottom             absolute bottom of the model *default* is 300,000 
    z_target_depth       Depth of deepest target, *default* is 50,000
    _utm_grid_size_east  size of a UTM grid in east direction. 
                         *default* is 640000 meters
    _utm_grid_size_north size of a UTM grid in north direction. 
                         *default* is 888960 meters

    ==================== ======================================================

    ..note:: If the survey steps across multiple UTM zones, then a 
                 distance will be added to the stations to place them in 
                 the correct location.  This distance is 
                 _utm_grid_size_north and _utm_grid_size_east.  You should 
                 these parameters to place the locations in the proper spot
                 as grid distances and overlaps change over the globe.

    ==================== ======================================================
    Methods              Description
    ==================== ======================================================
    make_mesh            makes a mesh from the given specifications
    plot_mesh            plots mesh to make sure everything is good
    write_initial_file   writes an initial model file that includes the mesh
    ==================== ======================================================
    """

    def __init__(self, **kwargs):  # edi_list=None,

        #        self.edi_list = edi_list
        self.Data = kwargs.pop('Data', None)

        # size of cells within station area in meters
        self.cell_size_east = kwargs.pop('cell_size_east', 500)
        self.cell_size_north = kwargs.pop('cell_size_north', 500)

        # padding cells on either side
        self.pad_east = kwargs.pop('pad_east', 7)
        self.pad_north = kwargs.pop('pad_north', 7)
        self.pad_z = kwargs.pop('pad_z', 4)

        # root of padding cells
        self.pad_stretch_h = kwargs.pop('pad_stretch_h', 1.2)
        self.pad_stretch_v = kwargs.pop('pad_stretch_v', 1.2)

        self.z1_layer = kwargs.pop('z1_layer', 10)
        self.z_target_depth = kwargs.pop('z_target_depth', 50000)
        self.z_bottom = kwargs.pop('z_bottom', 300000)

        # number of vertical layers
        self.n_layers = kwargs.pop('n_layers', 30)

        # number of air layers
        self.n_airlayers = kwargs.pop('n_airlayers', 0)
        # sea level in grid_z coordinates. Auto adjusts when topography read in?
        self.sea_level = 0.

        # strike angle to rotate grid to
        self.mesh_rotation_angle = kwargs.pop('mesh_rotation_angle', 0)

        # --> attributes to be calculated
        # station information
        if self.Data is not None:
            self.station_locations = self.Data.station_locations
        else:
            self.station_locations = None

        # grid nodes/cells sizes
        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        # grid lines locations arrays
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None

        # dictionary to contain any surfaces (e.g. topography)
        self.surfaces = {}

        # size of a utm grid - Not really used if we use epsg=3112 GA LLC
        self._utm_grid_size_north = 888960.0
        self._utm_grid_size_east = 640000.0
        self._utm_cross = False
        self._utm_ellipsoid = 23

        self.epsg = kwargs.pop('epsg', None)
        self.center_position_EN = kwargs.pop('center_position_EN', None)

        # if data object is provided, get epsg and center position from them
        if self.Data is not None:
            self.epsg = self.Data.epsg
            #for att in ('epsg', 'center_position_EN'):  # center_position_EN is not initially created in Data
                # attvalue = getattr(self.Data, att)
                # if attvalue is not None:
                #     setattr(self, att, attvalue)

        # resistivity model
        self.res_model = kwargs.pop('res_model', None)

        self.grid_center = None

        # initial file stuff
        self.model_fn = kwargs.pop('model_fn', None)
        self.save_path = kwargs.pop('save_path', None)
        self.model_fn_basename = kwargs.pop('model_fn_basename','ModEM_Model.ws')

        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
            self.model_fn_basename = os.path.basename(self.model_fn)

        self.title = 'Model File written by MTpy.modeling.modem'
        self.res_scale = kwargs.pop('res_scale', 'loge')


    def _reset_defaults_for_reading(self):
        """
        Reset all the defaults for input parameters prior to reading a model
        """
        # size of cells within station area in meters
        self.cell_size_east = None
        self.cell_size_north = None

        self.z1_layer = None
        self.z_target_depth = None
        self.z_bottom = None

        # number of vertical layers
        self.n_layers = None

        # number of air layers
        self.n_airlayers = None
        # sea level in grid_z coordinates. Auto adjusts when topography read in
        self.sea_level = 0.

    def make_mesh(self, update_data_center=True):
        """ 
        create finite element mesh according to user-input parameters.

        The mesh is built by first finding the center of the station area.  
        Then cells are added in the north and east direction with width
        cell_size_east and cell_size_north to cover all the station area.

        Padding cells are then added to extend the model to reduce
        edge effects.  The number of cells are pad_east and pad_north and the
        increase in size is by pad_stretch_h (pad_root_).

        The station locations are then computed as the center of the nearest cell
        as required by the code (what inversion code?)

        The vertical cells are built to increase in size exponentially with depth.
        The first cell depth is first_layer_thickness and should be
        about 1/10th the shortest skin depth.
        The layers then increase on a log scale to z_target_depth?
        Further deep, the model is padded with pad_z number of cells to extend the depth of the model.

        Note: Air-layers should NOT be added in this method, but when calling add_topography().
        air layers are added on top of the model constructed here.
        """

        # find the edges of the grid: bounding box of the survey area.
        nc_extra = 2
        west = self.station_locations['rel_east'].min() - self.cell_size_east * nc_extra
        east = self.station_locations['rel_east'].max() + self.cell_size_east * nc_extra
        south = self.station_locations['rel_north'].min() - self.cell_size_north * nc_extra
        north = self.station_locations['rel_north'].max() + self.cell_size_north * nc_extra

        # rounding appropriately.
        westr = np.round(west, -2)
        eastr = np.round(east, -2)
        southr = np.round(south, -2)
        northr = np.round(north, -2)
        #        # adjust center position (centre may be moved by rounding)
        #        self.Data.center_position_EN[0] += (westr + eastr - west - east)/2.
        #        self.Data.center_position_EN[1] += (southr + northr - south - north)/2.
        # -------make a grid around the stations from the parameters above-----
        # --> make grid in east-west direction
        # cells within station area
        east_gridr = np.arange(start=westr, stop=eastr + self.cell_size_east,
                               step=self.cell_size_east)
        print ("FZ: type of east_gridr = ", type(east_gridr))
        print ("FZ: east_gridr_1 = ", east_gridr)
        mean_egrid = np.mean(east_gridr)
        print("mean_egrid =", mean_egrid)
        if self.Data.rotation_angle == 0:
            self.Data.center_position_EN[0] -= mean_egrid
            self.station_locations['rel_east'] += mean_egrid
        east_gridr -= mean_egrid
        print ("FZ: east_gridr_2 = ", east_gridr)
        # padding cells in the east-west direction
        for ii in range(1, self.pad_east + 1):
            east_0 = float(east_gridr[-1])
            west_0 = float(east_gridr[0])
            # add_size = np.round(self.cell_size_east * self.pad_stretch_h * ii, -2) # -2 round to decimal left
            add_size = np.round(self.cell_size_east * self.pad_stretch_h ** ii, 2)
            pad_w = west_0 - add_size
            pad_e = east_0 + add_size
            east_gridr = np.insert(east_gridr, 0, pad_w)
            east_gridr = np.append(east_gridr, pad_e)

        # --> For some inversion code, need to make sure none of the stations lie on the nodes
        # this section would make the cell-sizes become unequal
        shift_station = 0.0  # originally = 0.02
        for s_east in sorted(self.station_locations['rel_east']):
            try:
                node_index = np.where(abs(s_east - east_gridr) <
                                      shift_station * self.cell_size_east)[0][0]
                if s_east - east_gridr[node_index] > 0:
                    east_gridr[node_index] -= shift_station * self.cell_size_east
                elif s_east - east_gridr[node_index] < 0:
                    east_gridr[node_index] += shift_station * self.cell_size_east
            except IndexError:
                continue

        # --> make grid in north-south direction
        # N-S cells with in station area
        north_gridr = np.arange(start=southr, stop=northr + self.cell_size_north,
                                step=self.cell_size_north)
        if self.Data.rotation_angle == 0:
            self.Data.center_position_EN[1] -= np.mean(north_gridr)
            self.station_locations['rel_north'] += np.mean(north_gridr)
        north_gridr -= np.mean(north_gridr)
        # padding cells in the east-west direction
        for ii in range(1, self.pad_north + 1):
            south_0 = float(north_gridr[0])
            north_0 = float(north_gridr[-1])
            # add_size = np.round(self.cell_size_north *self.pad_stretch_h * ii, -2)
            add_size = np.round(self.cell_size_north * self.pad_stretch_h ** ii, 2)
            pad_s = south_0 - add_size
            pad_n = north_0 + add_size
            north_gridr = np.insert(north_gridr, 0, pad_s)
            north_gridr = np.append(north_gridr, pad_n)

        # --> need to make sure none of the stations lie on the nodes
        for s_north in sorted(self.station_locations['rel_north']):
            try:
                node_index = np.where(abs(s_north - north_gridr) <
                                      shift_station * self.cell_size_north)[0][0]
                if s_north - north_gridr[node_index] > 0:
                    north_gridr[node_index] -= shift_station * self.cell_size_north
                elif s_north - north_gridr[node_index] < 0:
                    north_gridr[node_index] += shift_station * self.cell_size_north
            except IndexError:
                continue

        # keep the following section. may want to use later.
        # --> make depth gridz using logspace, target the depth to z_target_depth
        # log_z = np.logspace(np.log10(self.z1_layer),
        #                     np.log10(self.z_target_depth),
        #                     num=self.n_layers - self.pad_z - self.n_airlayers + 1)

        # deriv the z_cell size (vertical layers thickness)
        # log_z = log_z[1:] - log_z[:-1]  # the first layer thickness will not be equal to the intended z1_layer !!
        # #z_nodes = np.array([zz - zz % 10 ** np.floor(np.log10(zz)) for zz in log_z])  # why this to make round numbers?
        # z_nodes = log_z  # FZ: try not using the dubious code above.

        # FZ: use simple formula. relation between first z1_layer, stretch_v and target depth:
        p = self.pad_stretch_v
        nzf = np.log10((p - 1) * self.z_target_depth / self.z1_layer) / np.log10(p) - 1
        nz = int(nzf)
        if (nz > self.n_layers):
            self.n_layers = nz  # adjust z layers to prevent too big numbers.

        numz = self.n_layers - self.pad_z + 1  # - self.n_airlayers
        factorz = 1.2  # first few layers excluding the air_layers.
        exp_list = [self.z1_layer * (factorz ** nz) for nz in xrange(0, numz)]
        log_z = np.array(exp_list)
        z_nodes = log_z

        print("log_z logspace:", log_z)

        print("cell_sizes log_z = ", log_z)
        print("vs z_nodes = ", z_nodes)

        # index of top of padding
        itp = len(z_nodes) - 1

        print("index of top of padding itp=", itp)

        # padding cells in the end of the vertical direction
        for ii in range(1, self.pad_z + 1):
            z_0 = np.float(z_nodes[itp])
            # wrong: pad_d = np.round(z_0 * self.pad_stretch_v * ii, -2)
            pad_d = np.round(z_0 * self.pad_stretch_v ** ii, 2)
            z_nodes = np.append(z_nodes, pad_d)

        # JM said there should be no air layer in the mesh model ???
        # add air layers and define ground surface level.
        # initial layer thickness is same as z1_layer
        # add_air = self.n_airlayers
        add_air = 0  # set this =0 will not add any air layers below.
        z_nodes = np.hstack([[self.z1_layer] * add_air, z_nodes])

        # make an array of sum values as coordinates of the horizontal lines
        z_grid = np.array([z_nodes[:ii].sum()
                           for ii in range(z_nodes.shape[0] + 1)])

        # z_grid point at zero level
        # wrong: the following line does not make any sense if no air layer was added above.
        # incorrrect: self.sea_level = z_grid[self.n_airlayers]
        self.sea_level = z_grid[add_air]
        print("FZ:***1 sea_level = ", self.sea_level)

        # Need to make an array of the individual cell dimensions for modem
        east_nodes = east_gridr[1:] - east_gridr[:-1]
        north_nodes = north_gridr[1:] - north_gridr[:-1]

        # compute grid center
        center_east = -east_nodes.__abs__().sum() / 2
        center_north = -north_nodes.__abs__().sum() / 2
        center_z = 0
        self.grid_center = np.array([center_north, center_east, center_z])

        # make nodes/cells attributes
        self.nodes_east = east_nodes
        self.nodes_north = north_nodes
        self.nodes_z = z_nodes

        # make grid lines
        self.grid_east = east_gridr
        self.grid_north = north_gridr
        self.grid_z = z_grid

        # print ('self.nodes_z', self.nodes_z)  # FZ: cell sizes
        # print ('self.grid_z', self.grid_z)  # FZ: grid location
        # if desired, update the data center position (need to first project
        # east/north back to lat/lon) and rewrite to file
        if update_data_center:  # update the data file's centre position, reprojected back to degrees
            try:
                self.Data.center_position = self.Data.project_xy(self.Data.center_position_EN[0],
                                                                 self.Data.center_position_EN[1])
            except:
                pass

            print("make_mesh(): writing a data file, without topo, nor air layers.")
            self.Data.write_data_file(fill=False)

        # --> print out useful information
        print '-' * 15
        print '   Number of stations = {0}'.format(len(self.station_locations))
        print '   Dimensions: '
        print '      e-w = {0}'.format(east_gridr.shape[0])
        print '      n-s = {0}'.format(north_gridr.shape[0])
        print '       z  = {0} (including air layers: {1})'.format(z_grid.shape[0], self.n_airlayers)
        print '   Extensions: '
        print '      e-w = {0:.1f} (m)'.format(east_nodes.__abs__().sum())
        print '      n-s = {0:.1f} (m)'.format(north_nodes.__abs__().sum())
        print '      0-z = {0:.1f} (m)'.format(self.nodes_z.__abs__().sum())

        print '  Stations rotated by: {0:.1f} deg clockwise positive from N'.format(self.mesh_rotation_angle)
        print ''
        print ' ** Note ModEM does not accommodate mesh rotations, it assumes'
        print '    all coordinates are aligned to geographic N, E'
        print '    therefore rotating the stations will have a similar effect'
        print '    as rotating the mesh.'
        print '-' * 15

        if self._utm_cross is True:
            print '{0} {1} {2}'.format('-' * 25, 'NOTE', '-' * 25)
            print '   Survey crosses UTM zones, be sure that stations'
            print '   are properly located, if they are not, adjust parameters'
            print '   _utm_grid_size_east and _utm_grid_size_north.'
            print '   these are in meters and represent the utm grid size'
            print ' Example: '
            print ' >>> modem_model._utm_grid_size_east = 644000'
            print ' >>> modem_model.make_mesh()'
            print ''
            print '-' * 56

        return


    def add_topography(self, topographyfile=None, topographyarray=None, interp_method='nearest',
                       air_resistivity=1e17, sea_resistivity=0.3):
        """
        if air_layers is non-zero, will add topo: read in topograph file, make a surface model.
        Call project_stations_on_topography in the end, which will re-write the .dat file.

        If n_airlayers is zero, then canot add topo data, only bathymetry is needed.
        """
        # first, get surface data
        if topographyfile is not None:
            self.project_surface(surfacefile=topographyfile,
                                 surfacename='topography',
                                 method=interp_method)
        if topographyarray is not None:
            self.surface_dict['topography'] = topographyarray

        if self.n_airlayers > 0:
            # cell size is topomax/n_airlayers, rounded to nearest 1 s.f.
            cs = np.amax(self.surface_dict['topography'])/float(self.n_airlayers)
            #  cs = np.ceil(0.1*cs/10.**int(np.log10(cs)))*10.**(int(np.log10(cs))+1)
            cs = np.ceil(cs)

            # add air layers
            new_airlayers = np.linspace(
                0, self.n_airlayers, self.n_airlayers + 1) * cs
            add_z = new_airlayers[-1] - self.grid_z[self.n_airlayers]
            self.grid_z[self.n_airlayers + 1:] += add_z
            self.grid_z[:self.n_airlayers + 1] = new_airlayers

            # adjust the nodes
            self.nodes_z = self.grid_z[1:] - self.grid_z[:-1]

            # adjust sea level
            # wrong? self.sea_level = self.grid_z[self.n_airlayers]
            self.sea_level = self.grid_z[self.n_airlayers]
            print("FZ:***2 sea_level = ", self.sea_level)

            # assign topography
            self.assign_resistivity_from_surfacedata(
                'topography', air_resistivity, where='above')
        else:
            print "No air layers specified, so cannot add topography !!!"
            print "Only bathymetry shall be added below !!!"

        # assign sea water
        # first make a mask for all-land =1, which will be modified later according to air, water
        self.covariance_mask = np.ones_like(self.res_model)  # of grid size (xc, yc, zc)

        # assign model areas below sea level but above topography, as seawater
        # get grid centres
        gcz = np.mean([self.grid_z[:-1], self.grid_z[1:]], axis=0)

        # convert topography to local grid coordinates
        topo = self.sea_level - self.surface_dict['topography']
        # assign values
        for j in range(len(self.res_model)):
            for i in range(len(self.res_model[j])):
                # assign all sites above the topography to air
                ii1 = np.where(gcz <= topo[j, i])
                if len(ii1) > 0:
                    self.covariance_mask[j, i, ii1[0]] = 0.
                # assign sea water to covariance and model res arrays
                ii = np.where(
                    np.all([gcz > self.sea_level, gcz <= topo[j, i]], axis=0))
                if len(ii) > 0:
                    self.covariance_mask[j, i, ii[0]] = 9.
                    self.res_model[j, i, ii[0]] = sea_resistivity

        self.covariance_mask = self.covariance_mask[::-1]
        self.project_stations_on_topography()

        print ("FZ: write model file after air layers added ***** ")
        self.write_model_file(save_path=self.save_path)

        return

    def project_surface(self, surfacefile=None, surface=None, surfacename=None,
                        surface_epsg=4326, method='nearest'):
        """
        project a surface to the model grid and add resulting elevation data 
        to a dictionary called surface_dict. Assumes the surface is in lat/long
        coordinates (wgs84), if not, need to supply the epsg of the surface xy 
        points

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
        # initialise a dictionary to contain the surfaces
        if not hasattr(self, 'surface_dict'):
            self.surface_dict = {}

        # read the surface data in from ascii if surface not provided
        if surface is None:
            surface = elev_util.read_surface_ascii(surfacefile)

        x, y, elev = surface

        # if lat/lon provided as a 1D list, convert to a 2d grid of points
        if len(x.shape) == 1:
            x, y = np.meshgrid(x, y)

        epsg_from, epsg_to = surface_epsg, self.Data.epsg
        xs, ys = utm2ll.project(x, y, epsg_from, epsg_to)

        # get centre position of model grid in real world coordinates
        x0, y0 = [np.median(self.station_locations[dd] - self.station_locations['rel_' + dd]) for dd in
                  ['east', 'north']]

        # centre points of model grid in real world coordinates
        xg, yg = [np.mean([arr[1:], arr[:-1]], axis=0)
                  for arr in [self.grid_east + x0, self.grid_north + y0]]

        # elevation in model grid
        # first, get lat,lon points of surface grid
        points = np.vstack([arr.flatten() for arr in [xs, ys]]).T
        # corresponding surface elevation points
        values = elev.flatten()
        # xi, the model grid points to interpolate to
        xi = np.vstack([arr.flatten() for arr in np.meshgrid(xg, yg)]).T
        # elevation on the centre of the grid nodes
        elev_mg = spi.griddata(
            points, values, xi, method=method).reshape(len(yg), len(xg))

        print(" Elevation data over the meshgrid *** ",type(elev_mg), len(yg), len(xg), elev_mg.shape)

        np.savetxt('E:/tmp/elev_mg.txt', elev_mg, fmt='%10.5f')

        plt.imshow(elev_mg)
        plt.imshow(elev_mg[::-1])  # water shadow flip the image
        plt.show()

        # get a name for surface
        if surfacename is None:
            if surfacefile is not None:
                surfacename = os.path.basename(surfacefile)
            else:
                ii = 1
                surfacename = 'surface%01i' % ii
                while surfacename in self.surface_dict.keys():
                    ii += 1
                    surfacename = 'surface%01i' % ii

        # add surface to a dictionary of surface elevation data
        self.surface_dict[surfacename] = elev_mg

        return

    def assign_resistivity_from_surfacedata(self, surfacename, resistivity_value, where='above'):
        """
        assign resistivity value to all points above or below a surface
        requires the surface_dict attribute to exist and contain data for
        surface key (can get this information from ascii file using 
        project_surface)

        **inputs**
        surfacename = name of surface (must correspond to key in surface_dict)
        resistivity_value = value to assign
        where = 'above' or 'below' - assign resistivity above or below the 
                surface
        """

        gcz = np.mean([self.grid_z[:-1], self.grid_z[1:]], axis=0)

        # convert to positive down, relative to the top of the grid
        surfacedata = self.sea_level - self.surface_dict[surfacename]

        # define topography, so that we don't overwrite cells above topography
        # first check if topography exists
        if 'topography' in self.surface_dict.keys():
            # second, check topography isn't the surface we're trying to assign
            # resistivity for
            if surfacename == 'topography':
                topo = np.zeros_like(surfacedata)
            else:
                topo = self.sea_level - self.surface_dict['topography']
        # if no topography, assign zeros
        else:
            topo = self.sea_level + np.zeros_like(surfacedata)

        # assign resistivity value
        for j in range(len(self.res_model)):
            for i in range(len(self.res_model[j])):
                if where == 'above':
                    ii = np.where((gcz <= surfacedata[j, i]) & (
                        gcz > topo[j, i]))[0]
                else:
                    ii = np.where(gcz > surfacedata[j, i])[0]
                self.res_model[j, i, ii] = resistivity_value

    def project_stations_on_topography(self, air_resistivity=1e17):
        """
        This method is used in add_topography().
        It will Re-write the data file to change the elevation column.
        And update covariance mask according topo elevation model.
        :param air_resistivity:
        :return:
        """

        sx = self.station_locations['rel_east']
        sy = self.station_locations['rel_north']

        # find index of each station on grid
        for sname in self.station_locations['station']:
            ss = np.where(self.station_locations['station'] == sname)[0][0]
            # relative locations of stations
            sx, sy = self.station_locations['rel_east'][ss], \
                     self.station_locations['rel_north'][ss]
            # indices of stations on model grid
            sxi = np.where((sx <= self.grid_east[1:]) & (
                sx > self.grid_east[:-1]))[0][0]
            syi = np.where((sy <= self.grid_north[1:]) & (
                sy > self.grid_north[:-1]))[0][0]

            # first check if the site is in the sea
            if np.any(self.covariance_mask[::-1][syi, sxi] == 9):
                szi = np.amax(
                    np.where(self.covariance_mask[::-1][syi, sxi] == 9)[0])
            # second, check if there are any air cells
            elif np.any(self.res_model[syi, sxi] > 0.95 * air_resistivity):
                szi = np.amin(
                    np.where((self.res_model[syi, sxi] < 0.95 * air_resistivity))[0])
            # otherwise place station at the top of the model
            else:
                szi = 0

            # print("FZ:*** szi=", szi)
            # FZ: debug here to assign topography value for .dat file.

            # topoval = self.grid_z[szi]
            # self.station_locations['elev'][ss] = topoval # + 1.  # why +1 in elev ???
            # self.Data.data_array['elev'][ss] = topoval # + 1.

            # use topo elevation directly
            print(sname, ss, sxi, syi, szi)
            topoval = self.surface_dict['topography'][syi,sxi]
            print(sname,ss, sxi, syi, szi, topoval)

            self.station_locations['elev'][ss] = topoval # + 1.  # why +1 in elev ???
            self.Data.data_array['elev'][ss] = topoval # + 1.

        # This will shift stations' location to be relative to the defined mesh-grid centre
        self.Data.station_locations = self.station_locations

        print ("Re-write data file after adding topo")
        self.Data.write_data_file(fill=False)  # (Xi, Yi, Zi) of each station-i may be shifted

        # debug self.Data.write_data_file(save_path='/e/tmp', fill=False)
        print("FZ:*** what is self.grid_z=", self.grid_z.shape, self.grid_z)

        return

    def plot_mesh(self, east_limits=None, north_limits=None, z_limits=None,
                  **kwargs):
        """Plot the mesh to show model grid

        Arguments:
        ----------
            **east_limits** : tuple (xmin,xmax)
                             plot min and max distances in meters for the 
                             E-W direction.  If None, the east_limits
                             will be set to furthest stations east and west.
                             *default* is None

            **north_limits** : tuple (ymin,ymax)
                             plot min and max distances in meters for the 
                             N-S direction.  If None, the north_limits
                             will be set to furthest stations north and south.
                             *default* is None

            **z_limits** : tuple (zmin,zmax)
                            plot min and max distances in meters for the 
                            vertical direction.  If None, the z_limits is
                            set to the number of layers.  Z is positive down
                            *default* is None
        """

        fig_size = kwargs.pop('fig_size', [6, 6])
        fig_dpi = kwargs.pop('fig_dpi', 300)
        fig_num = kwargs.pop('fig_num', 1)

        station_marker = kwargs.pop('station_marker', 'v')
        marker_color = kwargs.pop('station_color', 'b')
        marker_size = kwargs.pop('marker_size', 2)

        line_color = kwargs.pop('line_color', 'k')
        line_width = kwargs.pop('line_width', .5)

        plt.rcParams['figure.subplot.hspace'] = .3
        plt.rcParams['figure.subplot.wspace'] = .3
        plt.rcParams['figure.subplot.left'] = .12
        plt.rcParams['font.size'] = 7

        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        plt.clf()

        # make a rotation matrix to rotate data
        # cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
        # sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))

        # turns out ModEM has not accomodated rotation of the grid, so for
        # now we will not rotate anything.
        cos_ang = 1
        sin_ang = 0

        # --->plot map view
        ax1 = fig.add_subplot(1, 2, 1, aspect='equal')

        # plot station locations
        plot_east = self.station_locations['rel_east']
        plot_north = self.station_locations['rel_north']

        ax1.scatter(plot_east,
                    plot_north,
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        east_line_xlist = []
        east_line_ylist = []
        north_min = self.grid_north.min()
        north_max = self.grid_north.max()
        for xx in self.grid_east:
            east_line_xlist.extend([xx * cos_ang + north_min * sin_ang,
                                    xx * cos_ang + north_max * sin_ang])
            east_line_xlist.append(None)
            east_line_ylist.extend([-xx * sin_ang + north_min * cos_ang,
                                    -xx * sin_ang + north_max * cos_ang])
            east_line_ylist.append(None)
        ax1.plot(east_line_xlist,
                 east_line_ylist,
                 lw=line_width,
                 color=line_color)

        north_line_xlist = []
        north_line_ylist = []
        east_max = self.grid_east.max()
        east_min = self.grid_east.min()
        for yy in self.grid_north:
            north_line_xlist.extend([east_min * cos_ang + yy * sin_ang,
                                     east_max * cos_ang + yy * sin_ang])
            north_line_xlist.append(None)
            north_line_ylist.extend([-east_min * sin_ang + yy * cos_ang,
                                     -east_max * sin_ang + yy * cos_ang])
            north_line_ylist.append(None)
        ax1.plot(north_line_xlist,
                 north_line_ylist,
                 lw=line_width,
                 color=line_color)

        # if east_limits == None:
        #     ax1.set_xlim(plot_east.min() - 50 * self.cell_size_east,
        #                  plot_east.max() + 50 * self.cell_size_east)
        # else:
        #     ax1.set_xlim(east_limits)
        #
        # if north_limits == None:
        #     ax1.set_ylim(plot_north.min() - 50 * self.cell_size_north,
        #                  plot_north.max() + 50 * self.cell_size_north)
        # else:
        #     ax1.set_ylim(north_limits)

        ax1.set_xlim(east_min, east_max)
        ax1.set_ylim(north_min, north_max)

        ax1.set_ylabel('Northing (m)', fontdict={'size': 9, 'weight': 'bold'})
        ax1.set_xlabel('Easting (m)', fontdict={'size': 9, 'weight': 'bold'})

        # ---------------------------------------
        # plot depth view
        ax2 = fig.add_subplot(1, 2, 2, aspect='auto', sharex=ax1)

        # plot the grid
        east_line_xlist = []
        east_line_ylist = []
        for xx in self.grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend([0,
                                    self.grid_z.max()])
            east_line_ylist.append(None)
        ax2.plot(east_line_xlist,
                 east_line_ylist,
                 lw=line_width,
                 color=line_color)

        z_line_xlist = []
        z_line_ylist = []
        for zz in self.grid_z:
            z_line_xlist.extend([self.grid_east.min(),
                                 self.grid_east.max()])
            z_line_xlist.append(None)
            z_line_ylist.extend([zz, zz])
            z_line_ylist.append(None)
        ax2.plot(z_line_xlist,
                 z_line_ylist,
                 lw=line_width,
                 color=line_color)

        # --> plot stations
        ax2.scatter(plot_east,
                    [0] * self.station_locations.shape[0],
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        if z_limits == None:
            ax2.set_ylim(self.z_target_depth, -200)
        else:
            ax2.set_ylim(z_limits)
        #
        # if east_limits == None:
        #     ax2.set_xlim(plot_east.min() - 50 * self.cell_size_east,
        #                  plot_east.max() + 50 * self.cell_size_east)
        # else:
        #     ax2.set_xlim(east_limits)

        ax2.set_ylabel('Depth (m)', fontdict={'size': 9, 'weight': 'bold'})
        ax2.set_xlabel('Easting (m)', fontdict={'size': 9, 'weight': 'bold'})

        plt.show()

        return

    def write_model_file(self, **kwargs):
        """
        will write an initial file for ModEM.  

        Note that x is assumed to be S --> N, y is assumed to be W --> E and
        z is positive downwards.  This means that index [0, 0, 0] is the 
        southwest corner of the first layer.  Therefore if you build a model
        by hand the layer block will look as it should in map view. 

        Also, the xgrid, ygrid and zgrid are assumed to be the relative 
        distance between neighboring nodes.  This is needed because wsinv3d 
        builds the  model from the bottom SW corner assuming the cell width
        from the init file.



        Key Word Arguments:
        ----------------------

            **nodes_north** : np.array(nx)
                        block dimensions (m) in the N-S direction. 
                        **Note** that the code reads the grid assuming that
                        index=0 is the southern most point.

            **nodes_east** : np.array(ny)
                        block dimensions (m) in the E-W direction.  
                        **Note** that the code reads in the grid assuming that
                        index=0 is the western most point.

            **nodes_z** : np.array(nz)
                        block dimensions (m) in the vertical direction.  
                        This is positive downwards.

            **save_path** : string
                          Path to where the initial file will be saved
                          to savepath/model_fn_basename

            **model_fn_basename** : string
                                    basename to save file to
                                    *default* is ModEM_Model.ws
                                    file is saved at savepath/model_fn_basename

            **title** : string
                        Title that goes into the first line 
                        *default* is Model File written by MTpy.modeling.modem 

            **res_model** : np.array((nx,ny,nz))
                        Prior resistivity model. 

                        .. note:: again that the modeling code 
                        assumes that the first row it reads in is the southern
                        most row and the first column it reads in is the 
                        western most column.  Similarly, the first plane it 
                        reads in is the Earth's surface.

            **res_scale** : [ 'loge' | 'log' | 'log10' | 'linear' ]
                            scale of resistivity.  In the ModEM code it 
                            converts everything to Loge, 
                            *default* is 'loge'

        """

        keys = ['nodes_east', 'nodes_north', 'nodes_z', 'title',
                'res_model', 'save_path', 'model_fn', 'model_fn_basename']

        for key in keys:
            try:
                setattr(self, key, kwargs[key])
            except KeyError:
                if self.__dict__[key] is None:
                    pass

        if self.save_path is not None:
            self.model_fn = os.path.join(self.save_path,
                                         self.model_fn_basename)
        if self.model_fn is None:
            if self.save_path is None:
                self.save_path = os.getcwd()
                self.model_fn = os.path.join(self.save_path,
                                             self.model_fn_basename)
            elif os.path.isdir(self.save_path) == True:
                self.model_fn = os.path.join(self.save_path,
                                             self.model_fn_basename)
            else:
                self.save_path = os.path.dirname(self.save_path)
                self.model_fn = self.save_path

        if self.res_model is None or type(self.res_model) is float or \
                        type(self.res_model) is int:
            res_model = np.zeros((self.nodes_north.shape[0],
                                  self.nodes_east.shape[0],
                                  self.nodes_z.shape[0]))

            if self.res_model is None:
                res_model[:, :, :] = 100.0
            else:
                res_model[:, :, :] = self.res_model

            self.res_model = res_model  # initial resistivity values

        if not hasattr(self, 'covariance_mask'):
            self.covariance_mask = np.ones_like(self.res_model)

        # --> write file
        ifid = file(self.model_fn, 'w')
        ifid.write('# {0}\n'.format(self.title.upper()))
        ifid.write('{0:>5}{1:>5}{2:>5}{3:>5} {4}\n'.format(self.nodes_north.shape[0],
                                                           self.nodes_east.shape[0],
                                                           self.nodes_z.shape[0],
                                                           0,
                                                           self.res_scale.upper()))

        # write S --> N node block
        for ii, nnode in enumerate(self.nodes_north):
            ifid.write('{0:>12.3f}'.format(abs(nnode)))

        ifid.write('\n')

        # write W --> E node block
        for jj, enode in enumerate(self.nodes_east):
            ifid.write('{0:>12.3f}'.format(abs(enode)))
        ifid.write('\n')

        # write top --> bottom node block
        for kk, zz in enumerate(self.nodes_z):
            ifid.write('{0:>12.3f}'.format(abs(zz)))
        ifid.write('\n')

        # write the resistivity in log e format
        if self.res_scale.lower() == 'loge':
            write_res_model = np.log(self.res_model[::-1, :, :])
        elif self.res_scale.lower() == 'log' or \
                        self.res_scale.lower() == 'log10':
            write_res_model = np.log10(self.res_model[::-1, :, :])
        elif self.res_scale.lower() == 'linear':
            write_res_model = self.res_model[::-1, :, :]

        # write out the layers from resmodel
        for zz in range(self.nodes_z.shape[0]):
            ifid.write('\n')
            for ee in range(self.nodes_east.shape[0]):
                for nn in range(self.nodes_north.shape[0]):
                    ifid.write('{0:>13.5E}'.format(
                        write_res_model[nn, ee, zz]))
                ifid.write('\n')

        if self.grid_center is None:
            # compute grid center
            center_east = -self.nodes_east.__abs__().sum() / 2
            center_north = -self.nodes_north.__abs__().sum() / 2
            center_z = 0
            self.grid_center = np.array([center_north, center_east, center_z])

        # Finally, write grid center coordinate and mesh rotation angle
        ifid.write('\n{0:>16.3f}{1:>16.3f}{2:>16.3f}\n'.format(self.grid_center[0],
                                                               self.grid_center[1], self.grid_center[2]))

        if self.mesh_rotation_angle is None:
            ifid.write('{0:>9.3f}\n'.format(0))
        else:
            ifid.write('{0:>9.3f}\n'.format(self.mesh_rotation_angle))
        ifid.close()

        print 'Wrote file to: {0}'.format(self.model_fn)

    def read_model_file(self, model_fn=None):
        """
        read an initial file and return the pertinent information including
        grid positions in coordinates relative to the center point (0,0) and 
        starting model.

        Note that the way the model file is output, it seems is that the 
        blocks are setup as 

        ModEM:                           WS:
        ----------                      ----- 
        0-----> N_north                 0-------->N_east
        |                               |
        |                               |
        V                               V
        N_east                          N_north


        Arguments:
        ----------

            **model_fn** : full path to initializing file.

        Outputs:
        --------

            **nodes_north** : np.array(nx)
                        array of nodes in S --> N direction

            **nodes_east** : np.array(ny) 
                        array of nodes in the W --> E direction

            **nodes_z** : np.array(nz)
                        array of nodes in vertical direction positive downwards

            **res_model** : dictionary
                        dictionary of the starting model with keys as layers

            **res_list** : list
                        list of resistivity values in the model

            **title** : string
                         title string

        """

        if model_fn is not None:
            self.model_fn = model_fn

        if self.model_fn is None:
            raise ModelError('model_fn is None, input a model file name')

        if not os.path.isfile(self.model_fn):
            raise ModelError(
                'Cannot find {0}, check path'.format(self.model_fn))

        self.save_path = os.path.dirname(self.model_fn)

        ifid = file(self.model_fn, 'r')
        ilines = ifid.readlines()
        ifid.close()

        self.title = ilines[0].strip()

        # get size of dimensions, remembering that x is N-S, y is E-W, z is +
        # down
        nsize = ilines[1].strip().split()
        n_north = int(nsize[0])
        n_east = int(nsize[1])
        n_z = int(nsize[2])
        log_yn = nsize[4]

        # get nodes
        self.nodes_north = np.array([np.float(nn)
                                     for nn in ilines[2].strip().split()])
        self.nodes_east = np.array([np.float(nn)
                                    for nn in ilines[3].strip().split()])
        self.nodes_z = np.array([np.float(nn)
                                 for nn in ilines[4].strip().split()])

        self.res_model = np.zeros((n_north, n_east, n_z))

        # get model
        count_z = 0
        line_index = 6
        count_e = 0
        while count_z < n_z:
            iline = ilines[line_index].strip().split()
            # blank lines spit the depth blocks, use those as a marker to
            # set the layer number and start a new block
            if len(iline) == 0:
                count_z += 1
                count_e = 0
                line_index += 1
            # 3D grid model files don't have a space at the end
            # additional condition to account for this.
            elif (len(iline) == 3) & (count_z == n_z - 1):
                count_z += 1
                count_e = 0
                line_index += 1
                # each line in the block is a line of N-->S values for an east
                # value
            else:
                north_line = np.array([float(nres) for nres in
                                       ilines[line_index].strip().split()])

                # Need to be sure that the resistivity array matches
                # with the grids, such that the first index is the
                # furthest south
                self.res_model[:, count_e, count_z] = north_line[::-1]

                count_e += 1
                line_index += 1

        # --> get grid center and rotation angle
        if len(ilines) > line_index:
            for iline in ilines[line_index:]:
                ilist = iline.strip().split()
                # grid center
                if len(ilist) == 3:
                    self.grid_center = np.array(ilist, dtype=np.float)
                # rotation angle
                elif len(ilist) == 1:
                    self.rotation_angle = np.float(ilist[0])
                else:
                    pass

        # --> make sure the resistivity units are in linear Ohm-m
        if log_yn.lower() == 'loge':
            self.res_model = np.e ** self.res_model
        elif log_yn.lower() == 'log' or log_yn.lower() == 'log10':
            self.res_model = 10 ** self.res_model

        # put the grids into coordinates relative to the center of the grid
        self.grid_north = np.array([self.nodes_north[0:ii].sum()
                                    for ii in range(n_north + 1)])
        self.grid_east = np.array([self.nodes_east[0:ii].sum()
                                   for ii in range(n_east + 1)])

        self.grid_z = np.array([self.nodes_z[:ii].sum()
                                for ii in range(n_z + 1)])
        # center the grids
        if self.grid_center is not None:
            self.grid_north += self.grid_center[0]
            self.grid_east += self.grid_center[1]
            self.grid_z += self.grid_center[2]

    def read_ws_model_file(self, ws_model_fn):
        """
        reads in a WS3INV3D model file
        """

        ws_model_obj = ws.WSModel(ws_model_fn)
        ws_model_obj.read_model_file()

        # set similar attributes
        for ws_key in ws_model_obj.__dict__.keys():
            for md_key in self.__dict__.keys():
                if ws_key == md_key:
                    setattr(self, ws_key, ws_model_obj.__dict__[ws_key])

        # compute grid center
        center_east = -self.nodes_east.__abs__().sum() / 2
        center_north = -self.nodes_norths.__abs__().sum() / 2
        center_z = 0
        self.grid_center = np.array([center_north, center_east, center_z])

    def write_vtk_file(self, vtk_save_path=None,
                       vtk_fn_basename='ModEM_model_res'):
        """
        write a vtk file to view in Paraview or other

        Arguments:
        -------------
            **vtk_save_path** : string
                                directory to save vtk file to.  
                                *default* is Model.save_path
            **vtk_fn_basename** : string
                                  filename basename of vtk file
                                  *default* is ModEM_model_res, evtk will add
                                  on the extension .vtr
        """

        if vtk_save_path is not None:
            vtk_fn = os.path.join(self.save_path, vtk_fn_basename)
        else:
            vtk_fn = os.path.join(vtk_save_path, vtk_fn_basename)

        # grids need to be n+1
        vtk_east = np.append(self.grid_east, 1.5 * self.grid_east[-1])
        vtk_north = np.append(self.grid_north, 1.5 * self.grid_north[-1])
        vtk_z = np.append(self.grid_z, 1.5 * self.grid_z[-1])
        gridToVTK(vtk_fn,
                  vtk_north,
                  vtk_east,
                  vtk_z,
                  pointData={'resistivity': self.res_model})

        print 'Wrote file to {0}'.format(vtk_fn)

    def write_gocad_sgrid_file(self, fn=None, origin=[0, 0, 0], clip=0, no_data_value=-99999):
        """
        write a model to gocad sgrid

        optional inputs:

        fn = filename to save to. File extension ('.sg') will be appended. 
             default is the model name with extension removed
        origin = real world [x,y,z] location of zero point in model grid
        clip = how much padding to clip off the edge of the model for export,
               provide one integer value or list of 3 integers for x,y,z directions
        no_data_value = no data value to put in sgrid

        """
        if not np.iterable(clip):
            clip = [clip, clip, clip]

            # determine save path
        savepath = None
        if fn is not None:
            savepath = os.path.dirname(fn)
            if len(savepath) == 0:
                savepath = None
        if savepath is None:
            savepath = os.path.dirname(self.model_fn)

        if fn is None:
            fn = os.path.join(os.path.dirname(self.model_fn),
                         os.path.basename(self.model_fn).split('.')[0])

        # number of cells in the ModEM model
        nyin, nxin, nzin = np.array(self.res_model.shape) + 1

        # get x, y and z positions
        gridedges = [self.grid_east[clip[0]:nxin - clip[0]] + origin[0],
                     self.grid_north[clip[1]:nyin - clip[1]] + origin[1],
                     -1. * self.grid_z[:nzin - clip[2]] - origin[2]]
        gridedges = np.meshgrid(*gridedges)

        # resistivity values, clipped to one smaller than grid edges
        resvals = self.res_model[clip[1]:nyin - clip[1] - 1,
                  clip[0]:nxin - clip[0] - 1, :nzin - clip[2] - 1]

        sgObj = mtgocad.Sgrid(resistivity=resvals, grid_xyz=gridedges,
                              fn=fn, workdir=savepath)
        sgObj.write_sgrid_file()

    def read_gocad_sgrid_file(self, sgrid_header_file, air_resistivity=1e39, sea_resistivity=0.3):
        """
        read a gocad sgrid file and put this info into a ModEM file.
        Note: can only deal with grids oriented N-S or E-W at this stage,
        with orthogonal coordinates

        """
        # read sgrid file
        sgObj = mtgocad.Sgrid()
        sgObj.read_sgrid_file(sgrid_header_file)
        self.sgObj = sgObj

        # check if we have a data object and if we do, is there a centre position
        # if not then assume it is the centre of the grid
        calculate_centre = True
        if self.Data is not None:
            if hasattr(self.Data, 'center_position_EN'):
                if self.Data.center_position_EN is not None:
                    centre = np.zeros(3)
                    centre[:2] = self.Data.center_position_EN
                    calculate_centre = False

                    # get resistivity model values
        self.res_model = sgObj.resistivity

        # get nodes and grid locations
        grideast, gridnorth, gridz = [
            np.unique(sgObj.grid_xyz[i]) for i in range(3)]
        gridz = np.abs(gridz)
        gridz.sort()
        if np.all(np.array([len(gridnorth), len(grideast), len(gridz)]) - 1 == np.array(self.res_model.shape)):
            self.grid_east, self.grid_north, self.grid_z = grideast, gridnorth, gridz
        else:
            print "Cannot read sgrid, can't deal with non-orthogonal grids or grids not aligned N-S or E-W"
            return

        # get nodes
        self.nodes_east = self.grid_east[1:] - self.grid_east[:-1]
        self.nodes_north = self.grid_north[1:] - self.grid_north[:-1]
        self.nodes_z = self.grid_z[1:] - self.grid_z[:-1]

        self.z1_layer = self.nodes_z[0]
        #        self.z_target_depth = None
        self.z_bottom = self.nodes_z[-1]

        # number of vertical layers
        self.n_layers = len(self.grid_z) - 1

        # number of air layers
        self.n_airlayers = sum(
            np.amax(self.res_model, axis=(0, 1)) > 0.9 * air_resistivity)

        # sea level in grid_z coordinates, calculate and adjust centre
        self.sea_level = self.grid_z[self.n_airlayers]

        print("FZ:***3 sea_level = ", self.sea_level)

        # get relative grid locations
        if calculate_centre:
            print "Calculating center position"
            centre = np.zeros(3)
            centre[0] = (self.grid_east.max() + self.grid_east.min()) / 2.
            centre[1] = (self.grid_north.max() + self.grid_north.min()) / 2.
        centre[2] = self.grid_z[self.n_airlayers]
        self.grid_east -= centre[0]
        self.grid_north -= centre[1]
        self.grid_z += centre[2]

    def write_xyres(self, location_type='EN', origin=[0, 0], model_epsg=None, depth_index='all',
                    savepath=None, outfile_basename='DepthSlice', log_res=False):
        """
        write files containing depth slice data (x, y, res for each depth)

        origin = x,y coordinate of zero point of ModEM_grid, or name of file
                 containing this info (full path or relative to model files)
        savepath = path to save to, default is the model object save path
        location_type = 'EN' or 'LL' xy points saved as eastings/northings or 
                        longitude/latitude, if 'LL' need to also provide model_epsg
        model_epsg = epsg number that was used to project the model
        outfile_basename = string for basename for saving the depth slices.
        log_res = True/False - option to save resistivity values as log10 
                               instead of linear

        """
        if savepath is None:
            savepath = self.save_path

        # make a directory to save the files
        savepath = os.path.join(savepath, outfile_basename)
        if not os.path.exists(savepath):
            os.mkdir(savepath)

        # try getting centre location info from file
        if type(origin) == str:
            try:
                origin = np.loadtxt(origin)
            except:
                print "Please provide origin as a list, array or tuple or as a valid filename containing this info"
                origin = [0, 0]

        # reshape the data
        x, y, z = [np.mean([arr[1:], arr[:-1]], axis=0) for arr in
                   [self.grid_east + origin[0], self.grid_north + origin[1], self.grid_z]]
        x, y = [arr.flatten() for arr in np.meshgrid(x, y)]

        # set format for saving data
        fmt = ['%.1f', '%.1f', '%.3e']

        # convert to lat/long if needed
        if location_type == 'LL':
            if np.any(origin) == 0:
                print "Warning, origin coordinates provided as zero, output lat/long are likely to be incorrect"
            x, y = utm2ll.project(x, y, model_epsg, 4326)
            # update format to accommodate lat/lon
            fmt[:2] = ['%.6f', '%.6f']

        # make depth indices into a list
        if depth_index == 'all':
            depthindices = range(len(z))
        elif np.iterable(depth_index):
            depthindices = np.array(depth_index).astype(int)
        else:
            depthindices = [depth_index]

        for k in depthindices:
            fname = os.path.join(savepath, outfile_basename + '_%1im.xyz' % z[k])
            vals = self.res_model[:, :, k].flatten()
            if log_res:
                vals = np.log10(vals)
                fmt[-1] = '%.3f'
            data = np.vstack([x, y, vals]).T
            np.savetxt(fname, data, fmt=fmt)


    def add_topography_to_model(self, dem_ascii_fn, model_fn, model_center=(0, 0),
                                rot_90=0, cell_size=500, elev_cell=30):
        """
        Add topography to an existing model from a dem in ascii format.

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

        Arguments:
        -------------
            *dem_ascii_fn* : string
                             full path to ascii dem file

            *model_fn* : string
                         full path to existing ModEM model file

            *model_center* : (east, north) in meters
                             Sometimes the center of the DEM and the center of the
                             model don't line up.  Use this parameter to line
                             everything up properly.

            *rot_90* : [ 0 | 1 | 2 | 3 ]
                       rotate the elevation model by rot_90*90 degrees.  Sometimes
                       the elevation model is flipped depending on your coordinate
                       system.

            *cell_size* : float (meters)
                          horizontal cell size of grid to interpolate elevation
                          onto.  This should be smaller or equal to the input
                          model cell size to be sure there is not spatial aliasing

            *elev_cell* : float (meters)
                          vertical size of each elevation cell.  This value should
                          be about 1/10th the smalles skin depth.

        Returns:
        ---------------
            *new_model_fn* : string
                             full path to model file that contains topography

        """
        # 1.) read in the dem and center it onto the resistivity model
        e_east, e_north, elevation = elev_util.read_dem_ascii(dem_ascii_fn, cell_size=cell_size,
                                                    model_center=model_center,
                                                    rot_90=3)
        plt.figure()
        plt.pcolormesh(e_east, e_north, elevation)
        m_obj = Model()
        m_obj.read_model_file(model_fn)
        # 2.) interpolate the elevation model onto the model grid
        m_elev = elev_util.interpolate_elevation(e_east, e_north, elevation,
                                       m_obj.grid_east, m_obj.grid_north, pad=3)
        # 3.) make a resistivity model that incoorporates topography
        mod_elev, elev_nodes_z = elev_util.make_elevation_model(m_elev, m_obj.nodes_z,
                                                      elevation_cell=elev_cell)
        plt.figure()
        #    plt.pcolormesh(m_obj.grid_east, m_obj.grid_north,m_elev)
        # 4.) write new model file
        m_obj.nodes_z = elev_nodes_z
        m_obj.res_model = mod_elev
        m_obj.write_model_file(model_fn_basename='{0}_topo.rho'.format(
            os.path.basename(m_obj.model_fn)[0:-4]))


    def change_data_elevation(self, data_fn, model_fn, new_data_fn=None, res_air=1e12):
        """
        At each station in the data file rewrite the elevation, so the station is
        on the surface, not floating in air.

        Arguments:
        ------------------
            *data_fn* : string
                        full path to a ModEM data file

            *model_fn* : string
                        full path to ModEM model file that has elevation
                        incoorporated.

            *new_data_fn* : string
                            full path to new data file name.  If None, then
                            new file name will add _elev.dat to input filename

            *res_air* : float
                        resistivity of air.  Default is 1E12 Ohm-m
        Returns:
        -------------
            *new_data_fn* : string
                            full path to new data file.
        """

        d_obj = Data()
        d_obj.read_data_file(data_fn)

        m_obj = Model()
        m_obj.read_model_file(model_fn)

        for key in d_obj.mt_dict.keys():
            mt_obj = d_obj.mt_dict[key]
            e_index = np.where(m_obj.grid_east > mt_obj.grid_east)[0][0]
            n_index = np.where(m_obj.grid_north > mt_obj.grid_north)[0][0]
            z_index = np.where(
                m_obj.res_model[n_index, e_index, :] < res_air * .9)[0][0]
            s_index = np.where(d_obj.data_array['station'] == key)[0][0]
            d_obj.data_array[s_index]['elev'] = m_obj.grid_z[z_index]

            mt_obj.grid_elev = m_obj.grid_z[z_index]

        if new_data_fn is None:
            new_dfn = '{0}{1}'.format(data_fn[:-4], '_elev.dat')
        else:
            new_dfn = new_data_fn

        d_obj.write_data_file(save_path=os.path.dirname(new_dfn),
                              fn_basename=os.path.basename(new_dfn),
                              compute_error=False,
                              fill=False)

        return new_dfn

#####################################################################################
