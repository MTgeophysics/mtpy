"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
from __future__ import print_function

import os
import sys

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats as stats, interpolate as spi

import mtpy
import mtpy.utils.calculator as mtcc
from mtpy.modeling import ws3dinv as ws
from mtpy.utils import mesh_tools as mtmesh, gis_tools as gis_tools, filehandling as mtfh
from mtpy.utils.mtpylog import MtPyLog
from .exception import ModelError
import mtpy.utils.gocad as mtgocad

try:
    from evtk.hl import gridToVTK
except ImportError:
    print('If you want to write a vtk file for 3d viewing, you need download '
          'and install evtk from https://bitbucket.org/pauloh/pyevtk')

__all__ = ['Model']


class Model(object):
    """
    make and read a FE mesh grid

    The mesh assumes the coordinate system where:
        x == North
        y == East
        z == + down

    All dimensions are in meters.

    The mesh is created by first making a regular grid around the station area,
    then padding cells are added that exponentially increase to the given
    extensions.  Depth cell increase on a log10 scale to the desired depth,
    then padding cells are added that increase exponentially.

    Arguments
    -------------
        **station_object** : mtpy.modeling.modem.Stations object
                            .. seealso:: mtpy.modeling.modem.Stations

    Examples
    -------------

    :Example 1 --> create mesh first then data file: ::

        >>> import mtpy.modeling.modem as modem
        >>> import os
        >>> # 1) make a list of all .edi files that will be inverted for
        >>> edi_path = r"/home/EDI_Files"
        >>> edi_list = [os.path.join(edi_path, edi)
                        for edi in os.listdir(edi_path)
        >>> ...         if edi.find('.edi') > 0]
        >>> # 2) Make a Stations object
        >>> stations_obj = modem.Stations()
        >>> stations_obj.get_station_locations_from_edi(edi_list)
        >>> # 3) make a grid from the stations themselves with 200m cell spacing
        >>> mmesh = modem.Model(station_obj)
        >>> # change cell sizes
        >>> mmesh.cell_size_east = 200,
        >>> mmesh.cell_size_north = 200
        >>> mmesh.ns_ext = 300000 # north-south extension
        >>> mmesh.ew_ext = 200000 # east-west extension of model
        >>> mmesh.make_mesh()
        >>> # check to see if the mesh is what you think it should be
        >>> msmesh.plot_mesh()
        >>> # all is good write the mesh file
        >>> msmesh.write_model_file(save_path=r"/home/modem/Inv1")
        >>> # create data file
        >>> md = modem.Data(edi_list, station_locations=mmesh.station_locations)
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")

    :Example 2 --> Rotate Mesh: ::

        >>> mmesh.mesh_rotation_angle = 60
        >>> mmesh.make_mesh()

    .. note:: ModEM assumes all coordinates are relative to North and East, and
             does not accommodate mesh rotations, therefore, here the rotation
             is of the stations, which essentially does the same thing.  You
             will need to rotate you data to align with the 'new' coordinate
             system.

    ==================== ======================================================
    Attributes           Description
    ==================== ======================================================
    _logger              python logging object that put messages in logging
                         format defined in logging configure file, see MtPyLog
                         more information
    cell_number_ew       optional for user to specify the total number of sells
                         on the east-west direction. *default* is None
    cell_number_ns       optional for user to specify the total number of sells
                         on the north-south direction. *default* is None
    cell_size_east       mesh block width in east direction
                         *default* is 500
    cell_size_north      mesh block width in north direction
                         *default* is 500
    grid_center          center of the mesh grid
    grid_east            overall distance of grid nodes in east direction
    grid_north           overall distance of grid nodes in north direction
    grid_z               overall distance of grid nodes in z direction
    model_fn             full path to initial file name
    model_fn_basename    default name for the model file name
    n_air_layers         number of air layers in the model. *default* is 0
    n_layers             total number of vertical layers in model
    nodes_east           relative distance between nodes in east direction
    nodes_north          relative distance between nodes in north direction
    nodes_z              relative distance between nodes in east direction
    pad_east             number of cells for padding on E and W sides
                         *default* is 7
    pad_north            number of cells for padding on S and N sides
                         *default* is 7
    pad_num              number of cells with cell_size with outside of
                         station area.  *default* is 3
    pad_method           method to use to create padding:
                         extent1, extent2 - calculate based on ew_ext and
                         ns_ext
                         stretch - calculate based on pad_stretch factors
    pad_stretch_h        multiplicative number for padding in horizontal
                         direction.
    pad_stretch_v        padding cells N & S will be pad_root_north**(x)
    pad_z                number of cells for padding at bottom
                         *default* is 4
    ew_ext               E-W extension of model in meters
    ns_ext               N-S extension of model in meters
    res_scale            scaling method of res, supports
                           'loge' - for log e format
                           'log' or 'log10' - for log with base 10
                           'linear' - linear scale
                         *default* is 'loge'
    res_list             list of resistivity values for starting model
    res_model            starting resistivity model
    res_initial_value    resistivity initial value for the resistivity model
                         *default* is 100
    mesh_rotation_angle  Angle to rotate the grid to. Angle is measured
                         positve clockwise assuming North is 0 and east is 90.
                         *default* is None
    save_path            path to save file to
    sea_level            sea level in grid_z coordinates. *default* is 0
    station_locations    location of stations
    title                title in initial file
    z1_layer             first layer thickness
    z_bottom             absolute bottom of the model *default* is 300,000
    z_target_depth       Depth of deepest target, *default* is 50,000
    ==================== ======================================================


    ==================== ======================================================
    Methods              Description
    ==================== ======================================================
    add_topography       if air_layers is non-zero, will add topo: read in
                         topograph file, make a surface model. Call
                         project_stations_on_topography in the end, which will
                         re-write the .dat file.
                         If n_airlayers is zero, then cannot add topo data,
                         only bathymetry is needed.
    add_topography_to_mesh    For a given mesh grid, use the topofile data to
                              define resistivity model. No new air layers will
                              be added. Just identify the max elev height as
                              the ref point
                              read in topograph file, make a surface model.
                              Call project_stations_on_topography in the end,
                              which will re-write the .dat file.
    add_topography_to_model   Add topography to an existing model from a dem
                              in ascii format.
    add_topography_to_model2    if air_layers is non-zero, will add topo: read
                                in topograph file, make a surface model.
                                Call project_stations_on_topography in the end,
                                which will re-write the .dat file.
                                If n_airlayers is zero, then cannot add topo
                                data, only bathymetry is needed.
    assign_resistivity_from_surfacedata     assign resistivity value to all
                                            points above or below a surface
                                            requires the surface_dict attribute
                                            to exist and contain data for
                                            surface key (can get this
                                            information from ascii file using
                                            project_surface)
    get_parameters       get important model parameters to write to a file for
                         documentation later.
    interpolate_elevation   interpolate the elevation onto the model grid.
    interpolate_elevation2  project a surface to the model grid and add
                            resulting elevation data to a dictionary called
                            surface_dict. Assumes the surface is in lat/long
                            coordinates (wgs84)
    make_mesh            makes a mesh from the given specifications
    make_mesh_from_center   The mesh is built by first finding the center of
                            the station area. Then cells are added in the north
                            and east direction with width cell_size_east and
                            cell_size_north to cover all the station area.
    make_z_mesh          Create a mesh grid for vertical Earth layers.
    make_z_mesh_exp      version of make_z_mesh method in order to create
                         exp-increasing cell sizes from the top layer
    make_z_mesh_new      new version of make_z_mesh. make_z_mesh and M
    plot_mesh            plots mesh to make sure everything is good
    plot_mesh_xy         add mesh grid lines in xy plan north-east map
    plot_mesh_xz         display the mesh in North-Depth aspect
    plot_topograph       display topography elevation data together with
                         station locations on a cell-index N-E map
    print_mesh_params    print out the mesh-related paramas
    print_model_file_summary    print out the summary of the model file
    project_stations_on_topography  This method is used in add_topography().
                                    It will Re-write the data file to change
                                    the elevation column. And update
                                    covariance mask according topo elevation
                                    model.
    project_surface      project a surface to the model grid and add resulting
                         elevation data to a dictionary called surface_dict.
                         Assumes the surface is in lat/long coordinates (wgs84),
                         if not, need to supply the epsg of the surface xy
                         points
    read_dem_ascii       read in dem which is ascii format
    read_model_file      read an initial file and return the pertinent
                         information including grid positions in coordinates
                         relative to the center point (0,0) and starting model.
    read_ws_model_file   reads in a WS3INV3D model file
    write_model_file     writes an initial model file that includes the mesh
    write_vtk_file       write a vtk file to view in Paraview or other
    ==================== ======================================================
    """

    def __init__(self, station_object=None, data_object=None, **kwargs):
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)

        self.station_locations = None
        self.data_obj = None
        if station_object is not None:
            self.station_locations = station_object
            self._logger.info("Use Station object as input, all functions that "
                              "uses data_objects are no longer available.")
        elif data_object is not None:
            self.data_obj = data_object
            self.station_locations = self.data_obj.station_locations
            self._logger.info("Use Data object as input.")
        # else:
        #     raise AttributeError("please provide either Station object or Data object as input")

        # size of cells within station area in meters
        self.cell_size_east = 500
        self.cell_size_north = 500

        # FZ: added this for user input number of cells in the horizontal mesh
        self.cell_number_ew = None
        self.cell_number_ns = None

        # padding cells on either side
        self.pad_east = 7
        self.pad_north = 7
        self.pad_z = 4

        self.pad_num = 3

        self.ew_ext = 100000
        self.ns_ext = 100000

        # root of padding cells
        self.pad_stretch_h = 1.2
        self.pad_stretch_v = 1.2

        # method to use to create padding
        self.pad_method = 'extent1'
        self.z_mesh_method = 'original' # method to make z mesh, 'original','original_refactor','exp' or 'new'
                                        # use: code embedded in make_mesh function, or make_z_mesh or 'make_z_mesh_exp' or 'make_z_mesh_new' respectively
                                        # temporary fix until I have a chance to test all 4

        self.z1_layer = 10
        self.z_target_depth = 50000
        self.z_bottom = 300000

        # number of vertical layers
        self.n_layers = 30

        # number of air layers
        self.n_air_layers = 0
        # sea level in grid_z coordinates. Auto adjusts when topography read in?
        self.sea_level = 0.

        # strike angle to rotate grid to
        self.mesh_rotation_angle = 0

        # --> attributes to be calculated
        # grid nodes
        self._nodes_east = None
        self._nodes_north = None
        self._nodes_z = None

        # grid locations
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None
        self.grid_center = None

        # resistivity model
        self.res_initial_value = 100.0
        self.res_model = None

        # initial file stuff
        self.model_fn = None
        self.save_path = os.getcwd()
        self.model_fn_basename = 'ModEM_Model_File.rho'
        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
            self.model_fn_basename = os.path.basename(self.model_fn)

        self.title = 'Model File written by MTpy.modeling.modem'
        self.res_scale = 'loge'

        for key in kwargs.keys():
            if hasattr(self, key):
                setattr(self, key, kwargs[key])
            else:
                self._logger.warn("Argument {}={} is not supportted thus not been set.".format(key, kwargs[key]))

    # --> make nodes and grid symbiotic so if you set one the other one
    #     gets set as well
    # Nodes East
    @property
    def nodes_east(self):
        if self.grid_east is not None:
            self._nodes_east = np.array([abs(self.grid_east[ii + 1] - self.grid_east[ii])
                                         for ii in range(self.grid_east.size - 1)])
        return self._nodes_east

    @nodes_east.setter
    def nodes_east(self, nodes):
        nodes = np.array(nodes)
        self._nodes_east = nodes
        self.grid_east = np.array([-nodes.sum() / 2 + nodes[0:ii].sum()
                                   for ii in range(nodes.size)] + [nodes.sum() / 2])

    # Nodes North
    @property
    def nodes_north(self):
        if self.grid_north is not None:
            self._nodes_north = np.array([abs(self.grid_north[ii + 1] - self.grid_north[ii])
                                          for ii in range(self.grid_north.size - 1)])
        return self._nodes_north

    @nodes_north.setter
    def nodes_north(self, nodes):
        nodes = np.array(nodes)
        self._nodes_north = nodes
        self.grid_north = np.array([-nodes.sum() / 2 + nodes[0:ii].sum()
                                    for ii in range(nodes.size)] + [nodes.sum() / 2])

    @property
    def nodes_z(self):
        if self.grid_z is not None:
            self._nodes_z = np.array([abs(self.grid_z[ii + 1] - self.grid_z[ii])
                                      for ii in range(self.grid_z.size - 1)])

            return self._nodes_z

    @nodes_z.setter
    def nodes_z(self, nodes):
        nodes = np.array(nodes)
        self._nodes_z = nodes
        self.grid_z = np.array([nodes[0:ii].sum() for ii in range(nodes.size)] + [nodes.sum()])

    def make_mesh(self):
        """
        create finite element mesh according to user-input parameters.

        The mesh is built by:
            1. Making a regular grid within the station area.
            2. Adding pad_num of cell_width cells outside of station area
            3. Adding padding cells to given extension and number of padding
               cells.
            4. Making vertical cells starting with z1_layer increasing
               logarithmically (base 10) to z_target_depth and num_layers.
            5. Add vertical padding cells to desired extension.
            6. Check to make sure none of the stations lie on a node.
               If they do then move the node by .02*cell_width

        """

        # --> find the edges of the grid
        # calculate the extra width of padding cells
        # multiply by 1.5 because this is only for 1 side
        pad_width_east = self.pad_num * 1.5 * self.cell_size_east
        pad_width_north = self.pad_num * 1.5 * self.cell_size_north

        # get the extremities
        west = self.station_locations.rel_east.min() - pad_width_east
        east = self.station_locations.rel_east.max() + pad_width_east
        south = self.station_locations.rel_north.min() - pad_width_north
        north = self.station_locations.rel_north.max() + pad_width_north

        # round the numbers so they are easier to read
        west = np.round(west, -2)
        east = np.round(east, -2)
        south = np.round(south, -2)
        north = np.round(north, -2)

        # -------make a grid around the stations from the parameters above------

        # adjust the edges so we have a whole number of cells
        add_ew = ((east - west) % self.cell_size_east) / 2.
        add_ns = ((north - south) % self.cell_size_north) / 2.

        # --> make the inner grid first
        inner_east = np.arange(west + add_ew - self.cell_size_east,
                               east - add_ew + 2 * self.cell_size_east,
                               self.cell_size_east)
        inner_north = np.arange(south + add_ns + self.cell_size_north,
                                north - add_ns + 2 * self.cell_size_north,
                                self.cell_size_north)

        # compute padding cells
        # first validate ew_ext and ns_ext to ensure it is large enough
        if 'extent' in self.pad_method:
            self._validate_extent(inner_east.min(),inner_east.max(),
                                  inner_north.min(),inner_north.max())
            
            
        if self.pad_method == 'extent1':
            padding_east = mtmesh.get_padding_cells(self.cell_size_east,
                                                    self.ew_ext / 2 - east,
                                                    self.pad_east,
                                                    self.pad_stretch_h)
            padding_north = mtmesh.get_padding_cells(self.cell_size_north,
                                                     self.ns_ext / 2 - north,
                                                     self.pad_north,
                                                     self.pad_stretch_h)
        elif self.pad_method == 'extent2':
            padding_east = mtmesh.get_padding_cells2(self.cell_size_east,
                                                     inner_east[-1],
                                                     self.ew_ext / 2.,
                                                     self.pad_east)
            padding_north = mtmesh.get_padding_cells2(self.cell_size_north,
                                                      inner_north[-1],
                                                      self.ns_ext / 2.,
                                                      self.pad_north)
        elif self.pad_method == 'stretch':
            padding_east = mtmesh.get_padding_from_stretch(self.cell_size_east,
                                                           self.pad_stretch_h,
                                                           self.pad_east)
            padding_north = mtmesh.get_padding_from_stretch(self.cell_size_north,
                                                            self.pad_stretch_h,
                                                            self.pad_north)
        else:
            raise NameError("Padding method \"{}\" is not supported".format(self.pad_method))

        # make the horizontal grid
        self.grid_east = np.append(np.append(-1 * padding_east[::-1] + inner_east.min(),
                                             inner_east),
                                   padding_east + inner_east.max())
        self.grid_north = np.append(np.append(-1 * padding_north[::-1] + inner_north.min(),
                                              inner_north),
                                    padding_north + inner_north.max())

        # --> need to make sure none of the stations lie on the nodes
        for s_east in sorted(self.station_locations.rel_east):
            try:
                node_index = np.where(abs(s_east - self.grid_east) <
                                      .02 * self.cell_size_east)[0][0]
                if s_east - self.grid_east[node_index] > 0:
                    self.grid_east[node_index] -= .02 * self.cell_size_east
                elif s_east - self.grid_east[node_index] < 0:
                    self.grid_east[node_index] += .02 * self.cell_size_east
            except IndexError:
                continue

        # --> need to make sure none of the stations lie on the nodes
        for s_north in sorted(self.station_locations.rel_north):
            try:
                node_index = np.where(abs(s_north - self.grid_north) <
                                      .02 * self.cell_size_north)[0][0]
                if s_north - self.grid_north[node_index] > 0:
                    self.grid_north[node_index] -= .02 * self.cell_size_north
                elif s_north - self.grid_north[node_index] < 0:
                    self.grid_north[node_index] += .02 * self.cell_size_north
            except IndexError:
                continue

        # --> make depth grid
        if self.z_mesh_method == 'original':
            log_z = np.logspace(np.log10(self.z1_layer),
                                np.log10(self.z_target_depth - np.logspace(np.log10(self.z1_layer),
                                                                           np.log10(self.z_target_depth),
                                                                           num=self.n_layers)[-2]),
                                num=self.n_layers - self.pad_z)
    
            z_nodes = np.array([np.round(zz, -int(np.floor(np.log10(zz)) - 1)) for zz in
                                log_z])

            # padding cells in the vertical
            z_padding = mtmesh.get_padding_cells(z_nodes[-1],
                                                 self.z_bottom - z_nodes.sum(),
                                                 self.pad_z,
                                                 self.pad_stretch_v)
            # make the blocks into nodes as oppose to total width
            z_padding = np.array([z_padding[ii + 1] - z_padding[ii]
                                  for ii in range(z_padding.size - 1)])
            
            self.nodes_z = np.append(z_nodes, z_padding)
        elif self.z_mesh_method == 'original_refactor':
            self.nodes_z,z_grid = self.make_z_mesh()
        elif self.z_mesh_method == 'exp':
            self.nodes_z,z_grid = self.make_z_mesh_exp()
        elif self.z_mesh_method == 'new':
            self.nodes_z,z_grid = self.make_z_mesh_new()

        else:
            raise NameError("Z mesh method \"{}\" is not supported".format(self.z_mesh_method))

        # compute grid center
        center_east = np.round(self.grid_east.min() - self.grid_east.mean(), -1)
        center_north = np.round(self.grid_north.min() - self.grid_north.mean(), -1)
        center_z = 0

        # this is the value to the lower left corner from the center.
        self.grid_center = np.array([center_north, center_east, center_z])

        # --> print out useful information
        self.print_mesh_params()

    def print_mesh_params(self, file=sys.stdout):
        # --> print out useful information
        print('-' * 15, file=file)
        print('\tNumber of stations = {0}'.format(len(self.station_locations.station)), file=file)
        print('\tDimensions: ', file=file)
        print('\t\te-w = {0}'.format(self.grid_east.size), file=file)
        print('\t\tn-s = {0}'.format(self.grid_north.size), file=file)
        print('\t\tz  = {0} (without 7 air layers)'.format(self.grid_z.size), file=file)
        print('\tExtensions: ', file=file)
        print('\t\te-w = {0:.1f} (m)'.format(self.nodes_east.__abs__().sum()), file=file)
        print('\t\tn-s = {0:.1f} (m)'.format(self.nodes_north.__abs__().sum()), file=file)
        print('\t\t0-z = {0:.1f} (m)'.format(self.nodes_z.__abs__().sum()), file=file)

        print('\tStations rotated by: {0:.1f} deg clockwise positive from N'.format(self.mesh_rotation_angle),
              file=file)
        print('', file=file)
        print(' ** Note ModEM does not accommodate mesh rotations, it assumes', file=file)
        print('    all coordinates are aligned to geographic N, E', file=file)
        print('    therefore rotating the stations will have a similar effect', file=file)
        print('    as rotating the mesh.', file=file)
        print('-' * 15, file=file)

    def make_mesh_from_center(self, update_data_center=True):
        """
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
        # first define some parameters. nc_extra_east and nc_extra_north is the number of cells outside the station
        # area (but with same cell size as inner cells - not padding). pad_east and pad_north is
        # number of padding cells, that increase with distance outward.
        nc_extra_east, pad_east = self.pad_east
        nc_extra_north, pad_north = self.pad_north

        if self.cell_number_ew is None:
            west = self.station_locations.rel_east.min() - self.cell_size_east * nc_extra_east
            east = self.station_locations.rel_east.max() + self.cell_size_east * nc_extra_east
        else:
            self._logger.debug("user specified cell number in east-west mesh %s" %
                               self.cell_number_ew)
            center_ew = 0.5 * (self.station_locations.rel_east.min() +
                               self.station_locations.rel_east.max())
            cell_num = int(self.cell_number_ew / 2)
            west = center_ew - self.cell_size_east * cell_num
            east = center_ew + self.cell_size_east * cell_num

        if self.cell_number_ns is None:
            south = self.station_locations.rel_north.min() - \
                    self.cell_size_north * nc_extra_north
            north = self.station_locations.rel_north.max() + \
                    self.cell_size_north * nc_extra_north
        else:
            self._logger.debug(
                "user specified cell number in north-south mesh %s" %
                self.cell_number_ns)
            center_ns = self.station_locations.rel_north.min() + \
                        self.station_locations.rel_north.max()
            center_ns = 0.5 * center_ns
            cell_num = int(self.cell_number_ns / 2)
            south = center_ns - self.cell_size_north * cell_num
            north = center_ns + self.cell_size_north * cell_num

        # rounding appropriately.
        west_r = np.round(west, -2)
        east_r = np.round(east, -2)
        south_r = np.round(south, -2)
        north_r = np.round(north, -2)
        #        # adjust center position (centre may be moved by rounding)
        #        self.Data.center_position_EN[0] += (westr + eastr - west - east)/2.
        #        self.Data.center_position_EN[1] += (southr + northr - south - north)/2.
        # -------make a grid around the stations from the parameters above-----
        # --> make grid in east-west direction
        # cells within station area
        east_grid_r = np.arange(start=west_r, stop=east_r + self.cell_size_east,
                                step=self.cell_size_east)

        self._logger.debug("FZ: east_gridr = %s" % east_grid_r)
        mean_egrid = np.mean(east_grid_r)
        self._logger.info("mean_egrid = %s" % mean_egrid)
        if self.data_obj.rotation_angle == 0:
            self.data_obj.center_position_utm[1] -= mean_egrid
            self.station_locations.rel_east += mean_egrid
        east_grid_r -= mean_egrid
        self._logger.debug("FZ: east_gridr_2 shifted centre = %s" % east_grid_r)

        # padding cells in the east-west direction
        for ii in range(1, pad_east + 1):
            east_0 = float(east_grid_r[-1])
            west_0 = float(east_grid_r[0])
            # add_size = mtcc.roundsf(self.cell_size_east * self.pad_stretch_h * ii, -2) # -2 round to decimal left
            # round to the nearest 2 significant figures
            add_size = mtcc.roundsf(self.cell_size_east * self.pad_stretch_h ** ii, 2)
            pad_w = west_0 - add_size
            pad_e = east_0 + add_size
            east_grid_r = np.insert(east_grid_r, 0, pad_w)
            east_grid_r = np.append(east_grid_r, pad_e)

        # --> For some inversion code, need to make sure none of the stations lie on the nodes
        # this section would make the cell-sizes become unequal
        shift_station = 0.0  # originally = 0.02
        for s_east in sorted(self.station_locations.rel_east):
            try:
                node_index = np.where(abs(s_east - east_grid_r) <
                                      shift_station * self.cell_size_east)[0][0]
                if s_east - east_grid_r[node_index] > 0:
                    east_grid_r[node_index] -= shift_station * self.cell_size_east
                elif s_east - east_grid_r[node_index] < 0:
                    east_grid_r[node_index] += shift_station * self.cell_size_east
            except IndexError:
                continue

        # --> make grid in north-south direction
        # N-S cells with in station area
        north_grid_r = np.arange(start=south_r, stop=north_r + self.cell_size_north,
                                 step=self.cell_size_north)
        if self.data_obj.rotation_angle == 0:
            self.data_obj.center_position_EN[1] -= np.mean(north_grid_r)
            self.station_locations.rel_north += np.mean(north_grid_r)
        north_grid_r -= np.mean(north_grid_r)
        # padding cells in the east-west direction
        for ii in range(1, pad_north + 1):
            south_0 = float(north_grid_r[0])
            north_0 = float(north_grid_r[-1])
            #            add_size = mtcc.roundsf(self.cell_size_north *self.pad_stretch_h * ii, -2)
            add_size = mtcc.roundsf(self.cell_size_north * self.pad_stretch_h ** ii, 2)
            pad_s = south_0 - add_size
            pad_n = north_0 + add_size
            north_grid_r = np.insert(north_grid_r, 0, pad_s)
            north_grid_r = np.append(north_grid_r, pad_n)

        # --> need to make sure none of the stations lie on the nodes
        for s_north in sorted(self.station_locations.rel_north):
            try:
                node_index = np.where(abs(s_north - north_grid_r) <
                                      shift_station * self.cell_size_north)[0][0]
                if s_north - north_grid_r[node_index] > 0:
                    north_grid_r[node_index] -= shift_station * self.cell_size_north
                elif s_north - north_grid_r[node_index] < 0:
                    north_grid_r[node_index] += shift_station * self.cell_size_north
            except IndexError:
                continue

        (z_nodes, z_grid) = self.make_z_mesh_new()

        # Need to make an array of the individual cell dimensions for modem
        east_nodes = east_grid_r[1:] - east_grid_r[:-1]
        north_nodes = north_grid_r[1:] - north_grid_r[:-1]

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
        self.grid_east = east_grid_r
        self.grid_north = north_grid_r
        self.grid_z = z_grid

        # print ('self.nodes_z', self.nodes_z)  # FZ: cell sizes
        # print ('self.grid_z', self.grid_z)  # FZ: grid location
        # if desired, update the data center position (need to first project
        # east/north back to lat/lon) and rewrite to file
        if update_data_center:  # update the data file's centre position, reprojected back to degrees
            try:
                self.data_obj.center_position = self.data_obj.project_xy(self.data_obj.center_position_EN[0],
                                                                         self.data_obj.center_position_EN[1])
            except:
                pass

            self._logger.info("writing a data file, without topo, nor air layers.")
            self.data_obj.write_data_file(fill=False)

        self.print_mesh_params()
        

    def make_z_mesh(self):
        """
        Create a mesh grid for vertical Earth layers.
        Refactored from the original make_mesh function, in order to modularize the logics.
        :return: (z_nodes, z_grid)
        """

        # keep the following section. may want to use later.
        # --> make depth gridz using logspace, target the depth to z_target_depth
        # log_z = np.logspace(np.log10(self.z1_layer),
        #                     np.log10(self.z_target_depth),
        #                     num=self.n_layers - self.pad_z - self.n_airlayers + 1)

        # derive the z_cell size (vertical layers thickness)
        # log_z = log_z[1:] - log_z[:-1]  # the first layer thickness will not be equal to the intended z1_layer !!
        # #z_nodes = np.array([zz - zz % 10 ** np.floor(np.log10(zz)) for zz in log_z])  # why this to make round numbers?
        # z_nodes = log_z
        # FZ: try not using these dubious code above.

        # FZ: use simple formula. relation between first z1_layer, stretch_v and target depth:
        p = self.pad_stretch_v
        nzf = np.log10((p - 1) * self.z_target_depth / self.z1_layer) / np.log10(p) - 1
        nz = int(nzf)
        if nz > self.n_layers:
            self.n_layers = nz  # adjust z layers to prevent too big numbers.

        num_z = self.n_layers - self.pad_z + 1  # - self.n_airlayers
        # numz = self.n_layers - self.pad_z + 1   - self.n_airlayers
        factor_z = 1.2  # first few layers excluding the air_layers.
        exp_list = [self.z1_layer * (factor_z ** nz) for nz in xrange(0, num_z)]
        log_z = np.array(exp_list)
        z_nodes = log_z

        self._logger.debug("cell_sizes log_z = %s" % log_z)
        self._logger.debug("and z_nodes = %s" % z_nodes)

        # index of top of padding
        itp = len(z_nodes) - 1

        self._logger.debug("index of top of padding itp= %s" % itp)

        # padding cells in the end of the vertical direction
        for ii in range(1, self.pad_z + 1):
            z_0 = np.float(z_nodes[itp])
            # wrong: pad_d = np.round(z_0 * self.pad_stretch_v * ii, -2)
            pad_d = np.round(z_0 * self.pad_stretch_v ** ii, 2)
            z_nodes = np.append(z_nodes, pad_d)

        # JM said there should be no air layer in this mesh building stage ???
        # add air layers and define ground surface level.
        # initial layer thickness is same as z1_layer
        # add_air = self.n_airlayers
        add_air = 0  # FZ: will add No air layers below if add_air=0
        z_nodes = np.hstack([[self.z1_layer] * add_air, z_nodes])

        # make an array of sum values as coordinates of the horizontal lines
        z_grid = np.array([z_nodes[:ii].sum()
                           for ii in range(z_nodes.shape[0] + 1)])

        # z_grid point at zero level
        # wrong: the following line does not make any sense if no air layer was added above.
        # incorrrect: self.sea_level = z_grid[self.n_airlayers]
        self.sea_level = z_grid[add_air]
        self._logger.debug("FZ:***1 sea_level = %s" % self.sea_level)

        return z_nodes, z_grid

    def make_z_mesh_exp(self, zfactor=1.2):
        """
        version of make_z_mesh method in order to create exp-increasing cell sizes from the top layer
        Create a mesh grid for vertical Earth layers.
        :return: (z_nodes, z_grid)
        """

        # FZ: use simple formula. relation between first z1_layer, stretch_v and target depth:

        # zfactor = first few vertical layers multiple factor

        # calculate the n_layers from target depth and first layer's thickness.
        if self.z_target_depth is not None:
            nf = np.log10((zfactor - 1) * self.z_target_depth / self.z1_layer) / np.log10(zfactor) - 1
            nz = int(nf)
            self._logger.debug("%s layers are needed to reach the target depth %s" %
                               (nz, self.z_target_depth))

        print("self.n_layers %s and calculated z_layers=%s" % (self.n_layers, nz))

        num_z = self.n_layers - self.pad_z  # speficied total_layers - padding_layers
        # numz should be close to nz.

        exp_list = [self.z1_layer * (zfactor ** n) for n in xrange(0, num_z)]
        log_z = np.array(exp_list)
        z_nodes = log_z  # a list of increasing cell sizes

        self._logger.debug("cell_sizes z_nodes = %s" % z_nodes)

        # index of top of padding
        itp = len(z_nodes) - 1

        self._logger.debug("index of top of padding itp= %s" % itp)

        # padding cells in the end of the vertical direction
        for ii in range(1, self.pad_z + 1):
            z_0 = np.float(z_nodes[itp])
            pad_d = np.round(z_0 * self.pad_stretch_v ** ii, 2)
            z_nodes = np.append(z_nodes, pad_d)

        # make an array of sum values as coordinates of the horizontal lines
        z_grid = np.array([z_nodes[:ii].sum() for ii in range(z_nodes.shape[0] + 1)])

        self._logger.debug("z grid lines = %s" % z_grid)
        self._logger.debug("Max vertical depth of the mesh= %s" % z_grid[-1])
        self._logger.debug("shape of vertical layers cells and grid lines %s, %s" %
                           (z_nodes.shape, z_grid.shape))

        return z_nodes, z_grid

    def make_z_mesh_new(self):
        """
        new version of make_z_mesh. make_z_mesh and M
        """

        # --> make depth grid
        # if n_airlayers < 0; set to 0
        nair = max(0, self.n_air_layers)

        log_z = mtcc.make_log_increasing_array(self.z1_layer, self.z_target_depth,
                                               self.n_layers - self.pad_z - nair)

        z_nodes = np.around(log_z[log_z < 100], decimals=-int(np.floor(np.log10(self.z1_layer))))
        z_nodes = np.append(z_nodes, np.around(log_z[log_z >= 100], decimals=-2))

        # index of top of padding
        itp = len(z_nodes) - 1

        # padding cells in the vertical direction
        for ii in range(1, self.pad_z + 1):
            z_0 = np.float(z_nodes[itp])
            pad_d = np.round(z_0 * self.pad_stretch_v ** ii, -2)
            z_nodes = np.append(z_nodes, pad_d)

        # add air layers and define ground surface level.
        # initial layer thickness is same as z1_layer
        z_nodes = np.hstack([[self.z1_layer] * self.n_air_layers, z_nodes])

        # make an array of absolute values
        z_grid = np.array([z_nodes[:ii].sum() for ii in range(z_nodes.shape[0] + 1)])

        return z_nodes, z_grid

    def add_topography_to_mesh(self, topography_file=None, topography_array=None,
                               interp_method='nearest',
                               air_resistivity=1e17, sea_resistivity=0.3):
        """
        For a given mesh grid, use the topofile data to define resistivity model.
        No new air layers will be added. Just identify the max elev height as the ref point

        read in topograph file, make a surface model.
        Call project_stations_on_topography in the end, which will re-write the .dat file.
        """

        sea_level = 0.0  # self.sea_level

        # first, get surface data
        if topography_file is not None:
            self.project_surface(surface_file=topography_file,
                                 surface_name='topography',
                                 method=interp_method)
        if topography_array is not None:
            self.surface_dict['topography'] = topography_array

        # update the z-centre as the max topo_elev (top air layer)
        self.grid_center[2] = -np.amax(self.surface_dict['topography'])

        # shift the grid_z lines to the new reference
        self.grid_z = self.grid_z + self.grid_center[2]

        self._logger.debug("New vertical grid lines = %s" % self.grid_z)

        self._logger.info("begin to self.assign_resistivity_from_surfacedata(...)")
        self.assign_resistivity_from_surfacedata('topography', air_resistivity, where='above')

        self._logger.info("begin to assign sea water resistivity")
        # first make a mask for all-land =1, which will be modified later according to air, water
        self.covariance_mask = np.ones_like(self.res_model)  # of grid size (xc, yc, zc)

        # assign model areas below sea level but above topography, as seawater
        # get the grid node centre points
        gcz = np.mean([self.grid_z[:-1], self.grid_z[1:]], axis=0)

        # convert topography to local grid coordinates
        topo = sea_level - self.surface_dict['topography']
        # assign values
        for j in range(len(self.res_model)):
            for i in range(len(self.res_model[j])):
                # assign all sites above the topography to air
                ii1 = np.where(gcz <= topo[j, i])
                if len(ii1) > 0:
                    self.covariance_mask[j, i, ii1[0]] = 0.
                # assign sea water to covariance and model res arrays
                ii = np.where(
                    np.all([gcz > sea_level, gcz <= topo[j, i]], axis=0))
                if len(ii) > 0:
                    self.covariance_mask[j, i, ii[0]] = 9.
                    self.res_model[j, i, ii[0]] = sea_resistivity

        self.covariance_mask = self.covariance_mask[::-1]

        self.station_grid_index = self.project_stations_on_topography()

        self._logger.debug("NEW res_model and cov_mask shapes: %s, %s" %
                           (self.res_model.shape,
                            self.covariance_mask.shape))

        return

    def add_topography(self, topography_file=None, topography_array=None, interp_method='nearest',
                       air_resistivity=1e17, sea_resistivity=0.3):
        """
        if air_layers is non-zero, will add topo: read in topograph file, make a surface model.
        Call project_stations_on_topography in the end, which will re-write the .dat file.

        If n_airlayers is zero, then cannot add topo data, only bathymetry is needed.
        """
        # first, get surface data
        if topography_file is not None:
            self.project_surface(surface_file=topography_file,
                                 surface_name='topography',
                                 method=interp_method)
        if topography_array is not None:
            self.surface_dict['topography'] = topography_array

        if self.n_air_layers is None or self.n_air_layers == 0:
            print("No air layers specified, so will not add air/topography !!!")
            print("Only bathymetry will be added below according to the topofile: sea-water low resistivity!!!")

        elif self.n_air_layers > 0:  # FZ: new logic, add equal blocksize air layers on top of the simple flat-earth grid
            # build air layers based on the inner core area
            padE = int(sum(self.pad_east))
            padN = int(sum(self.pad_north))
            topo_core = self.surface_dict['topography'][padN:-padN, padE:-padE]

            #            # compute the air cell size to be added = (topomax-topomin)/n_airlayers, rounded up to nearest whole number
            #            # use only the inner core area
            #            cs = (topo_core.max() - topo_core.min()) / float(self.n_airlayers)
            #            cs = np.ceil(cs)
            #            # define the bottom elevation of the bottom air layer, rounded to nearest whole number
            #            bottom_airlayer = int(round(topo_core.min()))
            #            # new air layers
            #            new_airlayers = -np.linspace(self.n_airlayers, 1, self.n_airlayers)*cs - bottom_airlayer


            # log increasing airlayers, in reversed order
            new_air_nodes = mtcc.make_log_increasing_array(self.z1_layer,
                                                           topo_core.max() - topo_core.min(),
                                                           self.n_air_layers,
                                                           increment_factor=0.999)[::-1]
            # sum to get grid cell locations
            new_airlayers = np.array([new_air_nodes[:ii].sum() for ii in range(len(new_air_nodes) + 1)])
            # round to nearest whole number and reverse the order
            new_airlayers = np.around(new_airlayers - topo_core.max())

            print("new_airlayers", new_airlayers)

            print("self.grid_z[0:2]", self.grid_z[0:2])

            # add new air layers, cut_off some tailing layers to preserve array size.
            self.grid_z = np.concatenate(
                [new_airlayers,
                 self.grid_z[self.n_air_layers + 1:] - self.grid_z[self.n_air_layers] + new_airlayers[-1]],
                axis=0)
            print(" NEW self.grid_z shape and values = ", self.grid_z.shape, self.grid_z)

            # adjust the nodes, which is simply the diff of adjacent grid lines
            self.nodes_z = self.grid_z[1:] - self.grid_z[:-1]

            # adjust sea level
            # wrong? self.sea_level = self.grid_z[self.n_airlayers]
            # self.sea_level = self.grid_z[self.n_airlayers]
            self.sea_level = 0  # keep as zero
            self._logger.debug("FZ:***2 sea_level = %s" % self.sea_level)

            # print (stop_here_for_debug)

        #        elif self.n_airlayers < 0: # if number of air layers < 0, auto calculate number of air layers required in air
        #
        #            # compute the air cell size to be added = topomax/n_airlayers, rounded to nearest 1 s.f.
        #            cs = np.amax(self.surface_dict['topography']) / float(self.n_airlayers)
        #            #  cs = np.ceil(0.1*cs/10.**int(np.log10(cs)))*10.**(int(np.log10(cs))+1)
        #            cs = np.ceil(cs)
        #
        #            # add air layers
        #            new_airlayers = np.linspace(
        #                0, self.n_airlayers, self.n_airlayers + 1) * cs
        #            add_z = new_airlayers[-1] - self.grid_z[self.n_airlayers]
        #            self.grid_z[self.n_airlayers + 1:] += add_z
        #            self.grid_z[:self.n_airlayers + 1] = new_airlayers
        #
        #            # adjust the nodes, which is simply the diff of adjacent grid lines
        #            self.nodes_z = self.grid_z[1:] - self.grid_z[:-1]
        #
        #            # adjust sea level
        #            # wrong? self.sea_level = self.grid_z[self.n_airlayers]
        #            self.sea_level = self.grid_z[self.n_airlayers]
        #            logger.debug("FZ:***2 sea_level = %s", self.sea_level)
        #
        #            # assign topography
        #            # self.assign_resistivity_from_surfacedata('topography', air_resistivity, where='above')
        #        else:
        #            pass

        # update the z-centre as the top air layer
        self.grid_center[2] = self.grid_z[0]

        self._logger.info("begin to self.assign_resistivity_from_surfacedata(...)")
        self.assign_resistivity_from_surfacedata('topography', air_resistivity, where='above')

        self._logger.info("begin to assign sea water resistivity")
        # first make a mask for all-land =1, which will be modified later according to air, water
        self.covariance_mask = np.ones_like(self.res_model)  # of grid size (xc, yc, zc)

        # assign model areas below sea level but above topography, as seawater
        # get grid node centres
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

        self.station_grid_index = self.project_stations_on_topography()

        self._logger.debug("NEW res_model and cov_mask shapes: %s, %s" %
                           (self.res_model.shape,
                            self.covariance_mask.shape))

        return

    def project_surface(self, surface_file=None, surface=None, surface_name=None,
                        surface_epsg=4326, method='nearest'):  # todo GA Version
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
            surface = mtfh.read_surface_ascii(surface_file)

        x, y, elev = surface

        # if lat/lon provided as a 1D list, convert to a 2d grid of points
        if len(x.shape) == 1:
            x, y = np.meshgrid(x, y)

        epsg_from, epsg_to = surface_epsg, self.data_obj.model_epsg
        xs, ys = mtpy.utils.gis_tools.epsg_project(x, y, epsg_from, epsg_to)

        # get centre position of model grid in real world coordinates
        x0, y0 = [np.median(self.station_locations.station_locations[dd] - self.station_locations.station_locations['rel_' + dd]) for dd in
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

        print(" Elevation data type and shape  *** ", type(elev_mg), elev_mg.shape, len(yg), len(xg))
        # <type 'numpy.ndarray'>  (65, 92), 65 92: it's 2D image with cell index as pixels
        # np.savetxt('E:/tmp/elev_mg.txt', elev_mg, fmt='%10.5f')

        # get a name for surface
        if surface_name is None:
            if surface_file is not None:
                surface_name = os.path.basename(surface_file)
            else:
                ii = 1
                surface_name = 'surface%01i' % ii
                while surface_name in self.surface_dict.keys():
                    ii += 1
                    surface_name = 'surface%01i' % ii

        # add surface to a dictionary of surface elevation data
        self.surface_dict[surface_name] = elev_mg

        return

    def assign_resistivity_from_surfacedata(self, surface_name, resistivity_value,
                                            where='above'):
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

        # FZ: should ref-define the self.res_model if its shape has changed after topo air layer are added

        gcz = np.mean([self.grid_z[:-1], self.grid_z[1:]], axis=0)

        self._logger.debug("gcz is the cells centre coordinates: %s, %s" %
                           (len(gcz), gcz))
        # convert to positive down, relative to the top of the grid
        surfacedata = self.sea_level - self.surface_dict[surface_name]
        # surfacedata = self.surface_dict[surfacename] - self.sea_level

        # define topography, so that we don't overwrite cells above topography
        # first check if topography exists
        if 'topography' in self.surface_dict.keys():
            # second, check topography isn't the surface we're trying to assign
            # resistivity for
            if surface_name == 'topography':
                # if it is, we need to define the upper limit as the top of the model
                top = np.zeros_like(surfacedata) + np.amin(self.grid_z)
            else:
                # if not, upper limit of resistivity assignment is the topography
                top = self.sea_level - self.surface_dict['topography']
        # if no topography, assign zeros
        else:
            top = self.sea_level + np.zeros_like(surfacedata)

        # assign resistivity value
        for j in range(len(self.res_model)):
            for i in range(len(self.res_model[j])):
                if where == 'above':
                    # needs to be above the surface but below the top (as defined before)
                    ii = np.where((gcz <= surfacedata[j, i]) & (gcz > top[j, i]))[0]

                else:  # for below the surface
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

        sx = self.station_locations.rel_east
        sy = self.station_locations.rel_north

        # find index of each station on grid
        station_index_x = []
        station_index_y = []
        for sname in self.station_locations.station:
            ss = np.where(self.station_locations.station == sname)[0][0]
            # relative locations of stations
            sx, sy = self.station_locations.rel_east[ss], \
                     self.station_locations.rel_north[ss]
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

            topoval = self.grid_z[szi]

            station_index_x.append(sxi)
            station_index_y.append(syi)

            #            # use topo elevation directly in modem.dat file
            #            !!! can't use topo elevation directly from topography file as the
            #                elevation needs to sit on the model mesh!
            #            topoval = self.surface_dict['topography'][syi, sxi]
            self._logger.debug("sname,ss, sxi, syi, szi, topoval: %s,%s,%s,%s,%s,%s"
                               % (sname, ss, sxi, syi, szi, topoval))

            # update elevation in station locations and data array, +1 m as
            # data elevation needs to be below the topography (as advised by Naser)
            self.station_locations.elev[ss] = topoval + 1.
            self.data_obj.data_array['elev'][ss] = topoval + 1.

        # This will shift stations' location to be relative to the defined mesh-grid centre
        self.data_obj.station_locations = self.station_locations

        self._logger.debug("Re-write data file after adding topo")
        self.data_obj.write_data_file(fill=False)  # (Xi, Yi, Zi) of each station-i may be shifted

        # debug self.Data.write_data_file(save_path='/e/tmp', fill=False)

        return station_index_x, station_index_y

    def plot_mesh(self, east_limits=None, north_limits=None, z_limits=None,
                  **kwargs):
        """
        Plot the mesh to show model grid

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
        # now we will not rotate anything (angle=0.0)
        cos_ang = 1
        sin_ang = 0

        # --->plot map view
        ax1 = fig.add_subplot(1, 2, 1, aspect='equal')

        # plot station locations
        plot_east = self.station_locations.rel_east
        plot_north = self.station_locations.rel_north

        # plot stations
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

        if east_limits is None:
            ax1.set_xlim(plot_east.min() - 10 * self.cell_size_east,
                         plot_east.max() + 10 * self.cell_size_east)
        else:
            ax1.set_xlim(east_limits)

        if north_limits is None:
            ax1.set_ylim(plot_north.min() - 10 * self.cell_size_north,
                         plot_north.max() + 10 * self.cell_size_east)
        else:
            ax1.set_ylim(north_limits)

        ax1.set_ylabel('Northing (m)', fontdict={'size': 9, 'weight': 'bold'})
        ax1.set_xlabel('Easting (m)', fontdict={'size': 9, 'weight': 'bold'})

        # ---------------------------------------
        # plot depth view along the east direction
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
                    [0] * self.station_locations.station.size,
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        if z_limits is None:
            ax2.set_ylim(self.z_target_depth, -200)
        else:
            ax2.set_ylim(z_limits)

        if east_limits is None:
            ax1.set_xlim(plot_east.min() - 10 * self.cell_size_east,
                         plot_east.max() + 10 * self.cell_size_east)
        else:
            ax1.set_xlim(east_limits)

        ax2.set_ylabel('Depth (m)', fontdict={'size': 9, 'weight': 'bold'})
        ax2.set_xlabel('Easting (m)', fontdict={'size': 9, 'weight': 'bold'})

        plt.show()

        return

    def plot_mesh_xy(self):
        """
        # add mesh grid lines in xy plan north-east map
        :return:
        """
        plt.figure(dpi=200)

        cos_ang = 1
        sin_ang = 0

        line_color = 'b'  # 'k'
        line_width = 0.5

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

        plt.plot(east_line_xlist, east_line_ylist, lw=line_width, color=line_color)

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

        plt.plot(north_line_xlist, north_line_ylist, lw=line_width, color=line_color)

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

        plt.xlim(east_min, east_max)
        plt.ylim(north_min, north_max)

        plt.ylabel('Northing (m)', fontdict={'size': 9, 'weight': 'bold'})
        plt.xlabel('Easting (m)', fontdict={'size': 9, 'weight': 'bold'})
        plt.title("Mesh grid in north-east dimension")

        plt.show()

        return

    def plot_mesh_xz(self):
        """
        display the mesh in North-Depth aspect
        :return:
        """
        station_marker = 'v'
        marker_color = 'b'
        marker_size = 2

        line_color = 'b'
        line_width = 0.5

        # fig = plt.figure(2, dpi=200)
        fig = plt.figure(dpi=200)
        plt.clf()
        ax2 = plt.gca()
        # ---------------------------------------
        # plot depth view along the north direction
        # ax2 = fig.add_subplot(1, 2, 2, aspect='auto', sharex=ax1)

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
        # ax2.scatter(plot_east, [0] * self.station_locations.shape[0],
        #            marker=station_marker, c=marker_color,s=marker_size)

        ax2.set_ylim(self.z_target_depth, -2000)

        #
        # if east_limits == None:
        #     ax2.set_xlim(plot_east.min() - 50 * self.cell_size_east,
        #                  plot_east.max() + 50 * self.cell_size_east)
        # else:
        #     ax2.set_xlim(east_limits)

        ax2.set_ylabel('Depth (m)', fontdict={'size': 9, 'weight': 'bold'})
        ax2.set_xlabel('Northing (m)', fontdict={'size': 9, 'weight': 'bold'})

        plt.show()

    def plot_topograph(self):
        """
        display topography elevation data together with station locations on a cell-index N-E map
        :return:
        """
        # fig_size = kwargs.pop('fig_size', [6, 6])
        # fig_dpi = kwargs.pop('fig_dpi', 300)
        # fig_num = kwargs.pop('fig_num', 1)
        #
        # station_marker = kwargs.pop('station_marker', 'v')
        # marker_color = kwargs.pop('station_color', 'b')
        # marker_size = kwargs.pop('marker_size', 2)
        #
        # line_color = kwargs.pop('line_color', 'k')
        # line_width = kwargs.pop('line_width', .5)
        #
        # plt.rcParams['figure.subplot.hspace'] = .3
        # plt.rcParams['figure.subplot.wspace'] = .3
        # plt.rcParams['figure.subplot.left'] = .12
        # plt.rcParams['font.size'] = 7

        # fig = plt.figure(3, dpi=200)
        fig = plt.figure(dpi=200)
        plt.clf()
        ax = plt.gca()

        # topography data image
        # plt.imshow(elev_mg) # this upside down
        # plt.imshow(elev_mg[::-1])  # this will be correct - water shadow flip of the image
        imgplot = plt.imshow(self.surface_dict['topography'],
                             origin='lower')  # the orgin is in the lower left corner SW.
        divider = make_axes_locatable(ax)
        # pad = separation from figure to colorbar
        cax = divider.append_axes("right", size="3%", pad=0.2)
        mycb = plt.colorbar(imgplot, cax=cax, use_gridspec=True)  # cmap=my_cmap_r, does not work!!
        mycb.outline.set_linewidth(2)
        mycb.set_label(label='Elevation (metre)', size=12)
        # make a rotation matrix to rotate data
        # cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
        # sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))

        # turns out ModEM has not accomodated rotation of the grid, so for
        # now we will not rotate anything.
        # cos_ang = 1
        # sin_ang = 0

        # --->plot map view
        # ax1 = fig.add_subplot(1, 2, 1, aspect='equal')

        # plot station locations in grid

        sgindex_x = self.station_grid_index[0]
        sgindex_y = self.station_grid_index[1]

        self._logger.debug("station grid index x: %s" % sgindex_x)
        self._logger.debug("station grid index y: %s" % sgindex_y)

        ax.scatter(sgindex_x, sgindex_y, marker='v', c='b', s=2)

        ax.set_xlabel('Easting Cell Index', fontdict={'size': 9, 'weight': 'bold'})
        ax.set_ylabel('Northing Cell Index', fontdict={'size': 9, 'weight': 'bold'})
        ax.set_title("Elevation and Stations in N-E Map (Cells)")

        plt.show()

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

            **res_starting_value** : float
                                     starting model resistivity value,
                                     assumes a half space in Ohm-m
                                     *default* is 100 Ohm-m

            **res_scale** : [ 'loge' | 'log' | 'log10' | 'linear' ]
                            scale of resistivity.  In the ModEM code it
                            converts everything to Loge,
                            *default* is 'loge'

        """
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

        if self.save_path is not None:
            self.model_fn = os.path.join(self.save_path,
                                         self.model_fn_basename)

        if self.model_fn is None:
            if self.save_path is None:
                self.save_path = os.getcwd()
                self.model_fn = os.path.join(self.save_path,
                                             self.model_fn_basename)
            elif os.path.isdir(self.save_path):
                self.model_fn = os.path.join(self.save_path,
                                             self.model_fn_basename)
            else:
                self.save_path = os.path.dirname(self.save_path)
                self.model_fn = self.save_path

        # get resistivity model
        if self.res_model is None:
            self.res_model = np.zeros((self.nodes_north.size,
                                       self.nodes_east.size,
                                       self.nodes_z.size))
            self.res_model[:, :, :] = self.res_initial_value

        elif type(self.res_model) in [float, int]:
            self.res_initial_value = self.res_model
            self.res_model = np.zeros((self.nodes_north.size,
                                       self.nodes_east.size,
                                       self.nodes_z.size))
            self.res_model[:, :, :] = self.res_initial_value

        # --> write file
        ifid = file(self.model_fn, 'w')
        ifid.write('# {0}\n'.format(self.title.upper()))
        ifid.write('{0:>5}{1:>5}{2:>5}{3:>5} {4}\n'.format(self.nodes_north.size,
                                                           self.nodes_east.size,
                                                           self.nodes_z.size,
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
        else:
            raise ModelError("resistivity scale \"{}\" is not supported.".format(self.res_scale))

        # write out the layers from resmodel
        for zz in range(self.nodes_z.size):
            ifid.write('\n')
            for ee in range(self.nodes_east.size):
                for nn in range(self.nodes_north.size):
                    ifid.write('{0:>13.5E}'.format(write_res_model[nn, ee, zz]))
                ifid.write('\n')

        if self.grid_center is None:
            # compute grid center
            center_east = -self.nodes_east.__abs__().sum() / 2
            center_north = -self.nodes_north.__abs__().sum() / 2
            center_z = 0
            self.grid_center = np.array([center_north, center_east, center_z])

        ifid.write('\n{0:>16.3f}{1:>16.3f}{2:>16.3f}\n'.format(self.grid_center[0],
                                                               self.grid_center[1], self.grid_center[2]))

        if self.mesh_rotation_angle is None:
            ifid.write('{0:>9.3f}\n'.format(0))
        else:
            ifid.write('{0:>9.3f}\n'.format(self.mesh_rotation_angle))
        ifid.close()

        self._logger.info('Wrote file to: {0}'.format(self.model_fn))

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

        if os.path.isfile(self.model_fn) is None:
            raise ModelError('Cannot find {0}, check path'.format(self.model_fn))

        self.save_path = os.path.dirname(self.model_fn)

        ifid = file(self.model_fn, 'r')
        ilines = ifid.readlines()
        ifid.close()

        self.title = ilines[0].strip()

        # get size of dimensions, remembering that x is N-S, y is E-W, z is + down
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
            # each line in the block is a line of N-->S values for an east value
            else:
                north_line = np.array([float(nres) for nres in iline])

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
                    self.mesh_rotation_angle = np.float(ilist[0])
                else:
                    pass

        # --> make sure the resistivity units are in linear Ohm-m
        if log_yn.lower() == 'loge':
            self.res_model = np.e ** self.res_model
        elif log_yn.lower() == 'log' or log_yn.lower() == 'log10':
            self.res_model = 10 ** self.res_model

        # center the grids
        if self.grid_center is None:
            self.grid_center = np.array([-self.nodes_north.sum() / 2,
                                         -self.nodes_east.sum() / 2,
                                         0.0])

        # need to shift the grid if the center is not symmetric
        shift_north = self.grid_center[0] + self.nodes_north.sum() / 2
        shift_east = self.grid_center[1] + self.nodes_east.sum() / 2

        # shift the grid.  if shift is + then that means the center is
        self.grid_north += shift_north
        self.grid_east += shift_east

        # get cell size
        self.cell_size_east = stats.mode(self.nodes_east)[0][0]
        self.cell_size_north = stats.mode(self.nodes_north)[0][0]

        # get number of padding cells
        self.pad_east = np.where(self.nodes_east[0:int(self.nodes_east.size / 2)]
                                 != self.cell_size_east)[0][-1]
        self.north_pad = np.where(self.nodes_north[0:int(self.nodes_north.size / 2)]
                                  != self.cell_size_north)[0][-1]

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
        center_north = -self.nodes_north.__abs__().sum() / 2
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

        if vtk_save_path is None:
            vtk_fn = os.path.join(self.save_path, vtk_fn_basename)
        else:
            vtk_fn = os.path.join(vtk_save_path, vtk_fn_basename)

        # use cellData, this makes the grid properly as grid is n+1
        gridToVTK(vtk_fn,
                  self.grid_north / 1000.,
                  self.grid_east / 1000.,
                  self.grid_z / 1000.,
                  cellData={'resistivity': self.res_model})

        self._logger.info('Wrote model file to {}'.format(vtk_fn))
        self.print_model_file_summary()

    def print_model_file_summary(self, file=sys.stdout):
        print('=' * 26, file=file)
        print('  model dimensions = {0}'.format(self.res_model.shape), file=file)
        print('     * north         {0}'.format(self.nodes_north.size), file=file)
        print('     * east          {0}'.format(self.nodes_east.size), file=file)
        print('     * depth         {0}'.format(self.nodes_z.size), file=file)
        print('=' * 26, file=file)

    def get_parameters(self):
        """
        get important model parameters to write to a file for documentation
        later.


        """

        parameter_list = ['cell_size_east',
                          'cell_size_north',
                          'ew_ext',
                          'ns_ext',
                          'pad_east',
                          'pad_north',
                          'pad_z',
                          'pad_num',
                          'z1_layer',
                          'z_target_depth',
                          'z_bottom',
                          'mesh_rotation_angle',
                          'res_initial_value',
                          'save_path']

        parameter_dict = {}
        for parameter in parameter_list:
            key = 'model.{0}'.format(parameter)
            parameter_dict[key] = getattr(self, parameter)

        parameter_dict['model.size'] = self.res_model.shape

        return parameter_dict


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
        if self.data_obj is not None:
            if hasattr(self.data_obj, 'center_position_EN'):
                if self.data_obj.center_position_EN is not None:
                    centre = np.zeros(3)
                    centre[:2] = self.data_obj.center_position_EN
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
            print("Cannot read sgrid, can't deal with non-orthogonal grids or grids not aligned N-S or E-W")
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
            print("Calculating center position")
            centre = np.zeros(3)
            centre[0] = (self.grid_east.max() + self.grid_east.min()) / 2.
            centre[1] = (self.grid_north.max() + self.grid_north.min()) / 2.
        centre[2] = self.grid_z[self.n_airlayers]
        self.grid_east -= centre[0]
        self.grid_north -= centre[1]
        self.grid_z += centre[2]



    # --> read in ascii dem file
    @staticmethod
    def read_dem_ascii(ascii_fn, cell_size=500, model_center=(0, 0),
                       rot_90=0, dem_rotation_angle=0):
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
        elevation = None
        with file(ascii_fn, 'r') as dfid:
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
                    # needs to be backwards because first line is the furthest north row.
                    elevation[:, -ii] = np.array(dline.strip().split(' '), dtype='float')
                else:
                    break

        # create lat and lon arrays from the dem fle
        lon = np.arange(x0, x0 + cs * nx, cs)
        lat = np.arange(y0, y0 + cs * ny, cs)

        # calculate the lower left and uper right corners of the grid in meters
        ll_en = gis_tools.project_point_ll2utm(lat[0], lon[0])
        ur_en = gis_tools.project_point_ll2utm(lat[-1], lon[-1])

        # estimate cell sizes for each dem measurement
        d_east = abs(ll_en[0] - ur_en[0]) / nx
        d_north = abs(ll_en[1] - ur_en[1]) / ny

        # calculate the number of new cells according to the given cell size
        # if the given cell size and cs are similar int could make the value 0,
        # hence the need to make it one if it is 0.
        num_cells = max([1, int(cell_size / np.mean([d_east, d_north]))])

        # make easting and northing arrays in meters corresponding to lat and lon
        east = np.arange(ll_en[0], ur_en[0], d_east)
        north = np.arange(ll_en[1], ur_en[1], d_north)

        # resample the data accordingly
        new_east = east[np.arange(0, east.size, num_cells)]
        new_north = north[np.arange(0, north.size, num_cells)]
        try:
            new_x, new_y = np.meshgrid(np.arange(0, east.shape[0], num_cells),
                                       np.arange(0, north.shape[0], num_cells),
                                       indexing='ij')
        except TypeError:
            new_x, new_y = [arr.T for arr in np.meshgrid(np.arange(0, east.shape[0], num_cells),
                                                         np.arange(0, north.shape[0], num_cells))]

        elevation = elevation[new_x, new_y]
        # make any null values set to minimum elevation, could be dangerous
        elevation[np.where(elevation == -9999.0)] = elevation[np.where(elevation != -9999.0)].min()

        # estimate the shift of the DEM to relative model coordinates
        mid_east = np.where(new_east >= model_center[0])[0][0]
        mid_north = np.where(new_north >= model_center[1])[0][0]

        new_east -= new_east[mid_east]
        new_north -= new_north[mid_north]

        # need to rotate cause I think I wrote the dem backwards
        if rot_90 == 1 or rot_90 == 3:
            elevation = np.rot90(elevation, rot_90)

        else:
            elevation = np.rot90(elevation, rot_90)

        if dem_rotation_angle != 0.0:
            cos_ang = np.cos(np.deg2rad(dem_rotation_angle))
            sin_ang = np.sin(np.deg2rad(dem_rotation_angle))
            rot_matrix = np.matrix(np.array([[cos_ang, sin_ang],
                                             [-sin_ang, cos_ang]]))

            new_coords = np.dot(rot_matrix, np.array([new_east, new_north]))
            new_east = new_coords[0]
            new_north = new_coords[1]

        return new_east, new_north, elevation

    @staticmethod
    def interpolate_elevation(elev_east, elev_north, elevation,
                              model_east, model_north, pad=3,
                              elevation_max=None):
        """
        interpolate the elevation onto the model grid.

        Arguments:
        ---------------

            **elev_east** : np.ndarray(num_east_nodes)
                          easting grid for elevation model

            **elev_north** : np.ndarray(num_north_nodes)
                          northing grid for elevation model

            **elevation** : np.ndarray(num_east_nodes, num_north_nodes)
                         elevation model assumes x is east, y is north
                         Units are meters

            **model_east** : np.ndarray(num_east_nodes_model)
                         relative easting grid of resistivity model

            **model_north** : np.ndarray(num_north_nodes_model)
                         relative northin grid of resistivity model

            **pad** : int
                    number of cells to repeat elevation model by.  So for pad=3,
                    then the interpolated elevation model onto the resistivity
                    model grid will have the outer 3 cells will be repeats of
                    the adjacent cell.  This is to extend the elevation model
                    to the resistivity model cause most elevation models will
                    not cover the entire area.

            **elevation_max** : float
                                maximum value for elevation
                                *default* is None, which will use
                                elevation.max()

        Returns:
        --------------

            **interp_elev** : np.ndarray(num_north_nodes_model, num_east_nodes_model)
                            the elevation model interpolated onto the resistivity
                            model grid.

        """
        # set a maximum on the elevation, used to get rid of singular high
        # points in the model
        if type(elevation_max) in [float, int]:
            max_find = np.where(elevation > float(elevation_max))
            elevation[max_find] = elevation_max

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

    def make_elevation_model(self, interp_elev, model_nodes_z,
                             elevation_cell=30, pad=3, res_air=1e12,
                             fill_res=100, res_sea=0.3):
        """
        Take the elevation data of the interpolated elevation model and map that
        onto the resistivity model by adding elevation cells to the existing model.

        ..Note: that if there are large elevation gains, the elevation cell size
                might need to be increased.

        Arguments:
        -------------
            **interp_elev** : np.ndarray(num_nodes_north, num_nodes_east)
                            elevation model that has been interpolated onto the
                            resistivity model grid. Units are in meters.

            **model_nodes_z** : np.ndarray(num_z_nodes_of_model)
                              vertical nodes of the resistivity model without
                              topography.  Note these are the nodes given in
                              relative thickness, not the grid, which is total
                              depth.  Units are meters.

            **elevation_cell** : float
                               height of elevation cells to be added on.  These
                               are assumed to be the same at all elevations.
                               Units are in meters

            **pad** : int
                    number of cells to look for maximum and minimum elevation.
                    So if you only want elevations within the survey area,
                    set pad equal to the number of padding cells of the
                    resistivity model grid.

            **res_air** : float
                        resistivity of air.  Default is 1E12 Ohm-m

            **fill_res** : float
                         resistivity value of subsurface in Ohm-m.

        Returns:
        -------------
            **elevation_model** : np.ndarray(num_north_nodes, num_east_nodes,
                                           num_elev_nodes+num_z_nodes)
                             Model grid with elevation mapped onto it.
                             Where anything above the surface will be given the
                             value of res_air, everything else will be fill_res

            **new_nodes_z** : np.ndarray(num_z_nodes+num_elev_nodes)
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
        self._logger.info('Number of elevation cells: {0}'.format(num_elev_cells))

        # find sea level if it is there
        if elev_min < 0:
            sea_level_index = num_elev_cells - abs(int(elev_min / elevation_cell)) - 1
        else:
            sea_level_index = num_elev_cells - 1

        self._logger.info('Sea level index is {0}'.format(sea_level_index))

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
                    dz = sea_level_index + abs(int((interp_elev[nn, ee]) / elevation_cell)) + 1
                    elevation_model[nn, ee, sea_level_index:dz] = res_sea
                else:
                    dz = int((elev_max - interp_elev[nn, ee]) / elevation_cell)
                    elevation_model[nn, ee, 0:dz] = res_air

        # make new z nodes array
        new_nodes_z = np.append(np.repeat(elevation_cell, num_elev_cells),
                                model_nodes_z)

        new_nodes_z[np.where(new_nodes_z < elevation_cell)] = elevation_cell

        return elevation_model, new_nodes_z

    def add_topography_to_model(self, dem_ascii_fn, write_file=True,
                                model_center=(0, 0), rot_90=0,
                                dem_rotation_angle=0, cell_size=500,
                                elev_cell=30, pad=1, elev_max=None):
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

        Arguments
        -------------
            **dem_ascii_fn** : string
                             full path to ascii dem file

            **model_fn** : string
                         full path to existing ModEM model file

            **model_center** : (east, north) in meters
                             Sometimes the center of the DEM and the center of the
                             model don't line up.  Use this parameter to line
                             everything up properly.

            **rot_90** : [ 0 | 1 | 2 | 3 ]
                       rotate the elevation model by rot_90*90 degrees.  Sometimes
                       the elevation model is flipped depending on your coordinate
                       system.

            **dem_rotation_angle: float (degrees from North)
                                  rotation angle to rotate station locations

            **cell_size** : float (meters)
                          horizontal cell size of grid to interpolate elevation
                          onto.  This should be smaller or equal to the input
                          model cell size to be sure there is not spatial aliasing

            **elev_cell** : float (meters)
                          vertical size of each elevation cell.  This value should
                          be about 1/10th the smalles skin depth.

        Returns
        ---------------
            **new_model_fn** : string
                             full path to model file that contains topography

        """
        # 1.) read in the dem and center it onto the resistivity model
        e_east, e_north, elevation = self.read_dem_ascii(dem_ascii_fn,
                                                         cell_size=cell_size,
                                                         model_center=model_center,
                                                         rot_90=rot_90,
                                                         dem_rotation_angle=dem_rotation_angle)

        # 2.) interpolate the elevation model onto the model grid
        m_elev = self.interpolate_elevation(e_east, e_north, elevation,
                                            self.grid_east, self.grid_north,
                                            pad=pad, elevation_max=elev_max)

        m_elev[np.where(m_elev == -9999.0)] = m_elev[np.where(m_elev != -9999.0)].min()
        # 3.) make a resistivity model that incoorporates topography
        mod_elev, elev_nodes_z = self.make_elevation_model(m_elev,
                                                           self.nodes_z,
                                                           elevation_cell=elev_cell)

        # 4.) write new model file
        self.nodes_z = elev_nodes_z
        self.res_model = mod_elev

        if write_file:
            self.save_path = os.path.dirname(self.model_fn)
            self.write_model_file(model_fn_basename='{0}_topo.rho'.format(
                os.path.basename(self.model_fn)[0:-4]))

            return self.model_fn

    def interpolate_elevation2(self, surfacefile=None, surface=None, surfacename=None,
                               method='nearest'):
        """
        project a surface to the model grid and add resulting elevation data
        to a dictionary called surface_dict. Assumes the surface is in lat/long
        coordinates (wgs84)

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
            surface = mtfh.read_surface_ascii(surfacefile)

        x, y, elev = surface

        # if lat/lon provided as a 1D list, convert to a 2d grid of points
        if len(x.shape) == 1:
            x, y = np.meshgrid(x, y)

        xs, ys, utm_zone = gis_tools.project_points_ll2utm(y, x,
                                                           epsg=self.station_locations.model_epsg,
                                                           utm_zone=self.station_locations.model_utm_zone
                                                           )

        # get centre position of model grid in real world coordinates
        x0, y0 = [np.median(
            self.station_locations.station_locations[dd] - self.station_locations.station_locations['rel_' + dd]) for dd
            in
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

        print(" Elevation data type and shape  *** ", type(elev_mg), elev_mg.shape, len(yg), len(xg))
        # <type 'numpy.ndarray'>  (65, 92), 65 92: it's 2D image with cell index as pixels
        # np.savetxt('E:/tmp/elev_mg.txt', elev_mg, fmt='%10.5f')

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

    def add_topography_to_model2(self, topographyfile=None, topographyarray=None,
                                 interp_method='nearest', air_resistivity=1e12):
        """
        if air_layers is non-zero, will add topo: read in topograph file, make a surface model.
        Call project_stations_on_topography in the end, which will re-write the .dat file.

        If n_airlayers is zero, then cannot add topo data, only bathymetry is needed.
        """
        # first, get surface data
        if topographyfile is not None:
            self.interpolate_elevation2(surfacefile=topographyfile,
                                        surfacename='topography',
                                        method=interp_method)
        if topographyarray is not None:
            self.surface_dict['topography'] = topographyarray

        if self.n_air_layers is None or self.n_air_layers == 0:
            self._logger.warn("No air layers specified, so will not add air/topography !!!")
            self._logger.warn("Only bathymetry will be added below according to the topofile: sea-water low resistivity!!!")

        elif self.n_air_layers > 0:  # FZ: new logic, add equal blocksize air layers on top of the simple flat-earth grid
            # build air layers based on the inner core area
            padE = self.pad_east
            padN = self.pad_north
            #            topo_core = self.surface_dict['topography'][padN:-padN,padE:-padE]
            gcx, gcy = [np.mean([arr[:-1], arr[1:]], axis=0) for arr in self.grid_east, self.grid_north]
            core_cells = mtmesh.get_station_buffer(gcx,
                                                   gcy,
                                                   self.station_locations.station_locations['rel_east'],
                                                   self.station_locations.station_locations['rel_north'],
                                                   buf=5 * (self.cell_size_east * 2 + self.cell_size_north ** 2) ** 0.5)
            topo_core = self.surface_dict['topography'][core_cells]

            # log increasing airlayers, in reversed order
            new_air_nodes = mtmesh.make_log_increasing_array(self.z1_layer,
                                                             topo_core.max() - topo_core.min(),
                                                             self.n_air_layers + 1,
                                                             increment_factor=0.999)[::-1]
            # sum to get grid cell locations
            new_airlayers = np.array([new_air_nodes[:ii].sum() for ii in range(len(new_air_nodes) + 1)])
            # round to nearest whole number and reverse the order
            new_airlayers = np.around(new_airlayers - topo_core.max())

            self._logger.debug("new_airlayers {}".format(new_airlayers))

            self._logger.debug("self.grid_z[0:2] {}".format(self.grid_z[0:2]))

            # add new air layers, cut_off some tailing layers to preserve array size.
            #            self.grid_z = np.concatenate([new_airlayers, self.grid_z[self.n_airlayers+1:] - self.grid_z[self.n_airlayers] + new_airlayers[-1]], axis=0)
            self.grid_z = np.concatenate([new_airlayers[:-1], self.grid_z + new_airlayers[-1]], axis=0)

        # print(" NEW self.grid_z shape and values = ", self.grid_z.shape, self.grid_z)
        #            print self.grid_z

        # update the z-centre as the top air layer
        self.grid_center[2] = self.grid_z[0]

        # update the resistivity model
        new_res_model = np.ones((self.nodes_north.size,
                                 self.nodes_east.size,
                                 self.nodes_z.size)) * self.res_initial_value
        new_res_model[:, :, self.n_air_layers + 1:] = self.res_model
        self.res_model = new_res_model

        #        logger.info("begin to self.assign_resistivity_from_surfacedata(...)")
        self.assign_resistivity_from_surfacedata('topography', air_resistivity, where='above')

        ##        logger.info("begin to assign sea water resistivity")
        #        # first make a mask for all-land =1, which will be modified later according to air, water
        #        self.covariance_mask = np.ones_like(self.res_model)  # of grid size (xc, yc, zc)
        #
        #        # assign model areas below sea level but above topography, as seawater
        #        # get grid node centres
        #        gcz = np.mean([self.grid_z[:-1], self.grid_z[1:]], axis=0)
        #
        #        # convert topography to local grid coordinates
        #        topo = -self.surface_dict['topography']
        #        # assign values
        #        for j in range(len(self.res_model)):
        #            for i in range(len(self.res_model[j])):
        #                # assign all sites above the topography to air
        #                ii1 = np.where(gcz <= topo[j, i])[0]
        #                if len(ii1) > 0:
        #                    self.covariance_mask[j, i, ii1] = 0.
        #                # assign sea water to covariance and model res arrays
        #                ii = np.where(
        #                    np.all([gcz > 0., gcz <= topo[j, i]], axis=0))[0]
        #                if len(ii) > 0:
        #                    self.covariance_mask[j, i, ii] = 9.
        #                    self.res_model[j, i, ii] = sea_resistivity
        #                    print "assigning sea", j, i, ii
        #
        #        self.covariance_mask = self.covariance_mask[::-1]

        #        self.station_grid_index = self.project_stations_on_topography()

        #        logger.debug("NEW res_model and cov_mask shapes: %s, %s", self.res_model.shape, self.covariance_mask.shape)

        return

    def _validate_extent(self,east,west,south,north,extent_ratio = 2.):
        """
        validate the provided ew_ext and ns_ext to make sure the model fits
        within these extents and allows enough space for padding according to 
        the extent ratio provided. If not, then update ew_ext and ns_ext parameters
        
        """
        inner_ew_ext = west - east
        inner_ns_ext = north - south
        
        if self.ew_ext < extent_ratio * inner_ew_ext:
            self._logger.warn("Provided or default ew_ext not sufficient to fit stations + padding, updating extent")
            self.ew_ext = np.ceil(extent_ratio * inner_ew_ext)

        if self.ns_ext < extent_ratio * inner_ns_ext:
            self._logger.warn("Provided or default ns_ext not sufficient to fit stations + padding, updating extent")
            self.ns_ext = np.ceil(extent_ratio * inner_ns_ext)
