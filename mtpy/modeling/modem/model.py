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
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats as stats, interpolate as spi

import mtpy.utils.calculator as mtcc
from mtpy.modeling import ws3dinv as ws
from mtpy.utils import mesh_tools as mtmesh, gis_tools as gis_tools, filehandling as mtfh
from mtpy.utils.mtpylog import MtPyLog
from .exception import ModelError
import mtpy.utils.gocad as mtgocad

try:
    from pyevtk.hl import gridToVTK
except ImportError:
    print('If you want to write a vtk file for 3d viewing, you need to '
          'install pyevtk')

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

    def __init__(self, stations_object=None, data_object=None, **kwargs):
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)

        self.station_locations = None
        self.data_obj = None

        if stations_object is not None:
            self.station_locations = stations_object# station location has to be moved
            # self.stations_obj = station_object.station_locations # station location has to be moved
            # self.data_obj = station_object # data_obj has to be updted
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

        self.z1_layer = 10
        self.z_layer_rounding = None
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
        self.grid_z = kwargs.pop('grid_z',None)
        if self.grid_z is not None:
            self.n_layers = len(self.grid_z)
            self.z_mesh_method = 'custom'
        else:
            self.z_mesh_method = 'new'
        if 'z_mesh_method' in list(kwargs.keys()):
            self.z_mesh_method = kwargs['z_mesh_method']
 
        # method to use to create padding
        self.pad_method = 'extent1'
       
        self.grid_center = None

        # resistivity model
        self.res_initial_value = 100.0
        self.res_model = None

        # initial file stuff
        self.model_fn = None
        self.save_path = None
        self.model_fn_basename = 'ModEM_Model_File.rho'
        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
            self.model_fn_basename = os.path.basename(self.model_fn)

        self.title = 'Model File written by MTpy.modeling.modem'
        self.res_scale = 'loge'

        for key in list(kwargs.keys()):
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
        self.grid_east = np.array([nodes[0:ii].sum()#-nodes.sum() / 2 + 
                                   for ii in range(nodes.size+1)])# + [shift])#[nodes.sum() / 2]

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
        self.grid_north = np.array([nodes[0:ii].sum()#-nodes.sum() / 2 + 
                                    for ii in range(nodes.size+1)])# + [shift])#[nodes.sum() / 2]

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

    # need some arrays for plotting that are the same length as the
    # resistivity model
    @property
    def plot_east(self):
        plot_east = np.array([self.nodes_east[0:ii].sum() 
                             for ii in range(self.nodes_east.size)])
        return plot_east-plot_east[-1]/2.
    
    @property
    def plot_north(self):
        plot_north = np.array([self.nodes_north[0:ii].sum() 
                          for ii in range(self.nodes_north.size)])
        return plot_north-plot_north[-1]/2.
    
    @property
    def plot_z(self):
        return np.array([self.nodes_z[0:ii].sum() 
                         for ii in range(self.nodes_z.size)])
    
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

        if self.z_mesh_method == 'custom':
            if self.grid_z is None:
                self.z_mesh_method = 'new'
                self._logger.warn('No grid_z provided, creating new z mesh using default method')
        
        if self.z_mesh_method == 'custom':
                self.nodes_z, z_grid = self.grid_z[1:]-self.grid_z[:-1], self.grid_z
        elif self.z_mesh_method == 'new':
            self.nodes_z, z_grid = self.make_z_mesh_new()
        else:
            raise NameError("Z mesh method \"{}\" is not supported".format(self.z_mesh_method))

        # compute grid center
        center_east = np.round(self.grid_east.min() - self.grid_east.mean(), -1)
        center_north = np.round(self.grid_north.min() - self.grid_north.mean(), -1)
        center_z = 0

        # this is the value to the lower left corner from the center.
        self.grid_center = np.array([center_north, center_east, center_z])
        
        # make the resistivity array
        self.res_model = np.zeros((self.nodes_north.size,
                                  self.nodes_east.size,
                                  self.nodes_z.size))
        self.res_model[:, :, :] = self.res_initial_value

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

    def make_z_mesh_new(self, n_layers=None):
        """
        new version of make_z_mesh. make_z_mesh and M
        """
        n_layers = self.n_layers if n_layers is None else n_layers

        # --> make depth grid
        # if n_airlayers < 0; set to 0
        log_z = mtcc.make_log_increasing_array(self.z1_layer,
                                               self.z_target_depth,
                                               n_layers - self.pad_z)

        if self.z_layer_rounding is not None:
            z_nodes = np.around(log_z, decimals=self.z_layer_rounding)
        else:
            # round any values less than 100 to the same s.f. as z1_layer
            z_nodes = np.around(log_z[log_z < 100],
                                decimals=-int(np.floor(np.log10(self.z1_layer))))
            # round any values greater than or equal to 100 to the nearest 100
            z_nodes = np.append(z_nodes, np.around(log_z[log_z >= 100],
                                                   decimals=-2))

        # index of top of padding
        #itp = len(z_nodes) - 1

        # padding cells in the vertical direction
        z_0 = np.float(z_nodes[-1])
        for ii in range(1, self.pad_z + 1):
            pad_d = np.round(z_0 * self.pad_stretch_v ** ii, -2)
            z_nodes = np.append(z_nodes, pad_d)
        # add air layers and define ground surface level.
        # initial layer thickness is same as z1_layer
        # z_nodes = np.hstack([[z1_layer] * n_air, z_nodes])

        # make an array of absolute values
        z_grid = np.array([z_nodes[:ii].sum() for ii in range(z_nodes.shape[0] + 1)])

        return z_nodes, z_grid
  
    def add_layers_to_mesh(self, n_add_layers=None, layer_thickness=None,
                           where='top'):
        """
        Function to add constant thickness layers to the top or bottom of mesh.
        Note: It is assumed these layers are added before the topography. If 
        you want to add topography layers, use function add_topography_to_model2

        :param n_add_layers: integer, number of layers to add
        :param layer_thickness: real value or list/array. Thickness of layers,
                                defaults to z1 layer. Can provide a single value
                                or a list/array containing multiple layer
                                thicknesses.
        :param where: where to add, top or bottom
   
        
        """
        # create array containing layers to add
        if layer_thickness is None:
            layer_thickness = self.z1_layer
        if np.iterable(layer_thickness):
            add_layers = np.insert(np.cumsum(layer_thickness),0,0)[:-1]
            layer_thickness = layer_thickness[-1]
            
            if n_add_layers != len(add_layers):
                self._logger.warn("Updating number of layers to reflect the length of the layer thickness array")
            n_add_layers = len(add_layers)
        else:
            add_layers = np.arange(0,n_add_layers*layer_thickness,layer_thickness)
            
        # create a new z grid
        self.grid_z = np.hstack([add_layers,self.grid_z + add_layers[-1] + layer_thickness])
        
        # update the number of layers
        self.n_layers = len(self.grid_z) - 1
        
        # add the extra layer to the res model
        self.res_model = np.vstack([self.res_model[:,:,:n_add_layers].T,self.res_model.T]).T

    def assign_resistivity_from_surfacedata(self, top_surface, bottom_surface, resistivity_value):
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

        # assign resistivity value
        for j in range(len(self.res_model)):
            for i in range(len(self.res_model[j])):
                ii = np.where((gcz > top_surface[j, i]) & (gcz <= bottom_surface[j, i]))[0]
                self.res_model[j, i, ii] = resistivity_value

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
        plot_names = kwargs.pop('plot_names', False)

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
        if plot_names:
            for s_arr in self.station_locations.station_locations:
                ax1.text(s_arr['rel_east'], s_arr['rel_north']+.05,
                         s_arr['station'], ha='center', va='baseline',
                         clip_on=True)
                
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
            east_line_ylist.extend([self.grid_z.min(),
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
                    self.station_locations.rel_elev,
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        if z_limits is None:
            ax2.set_ylim(self.z_target_depth, self.grid_z.min() - 200)
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

    def plot_topography(self):
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
        fig.clf()
        ax = fig.add_subplot(1, 1, 1, aspect='equal') 

        x, y = np.meshgrid(self.grid_east, self.grid_north)
        # topography data image
        # plt.imshow(elev_mg) # this upside down
        # plt.imshow(elev_mg[::-1])  # this will be correct - water shadow flip of the image

        imgplot = ax.pcolormesh(x, y, self.surface_dict['topography'])
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

#        sgindex_x = self.station_grid_index[0]
#        sgindex_y = self.station_grid_index[1]
#
#        self._logger.debug("station grid index x: %s" % sgindex_x)
#        self._logger.debug("station grid index y: %s" % sgindex_y)

        ax.scatter(self.station_locations.rel_east,
                   self.station_locations.rel_north,
                   marker='v', c='k', s=2)

        ax.set_xlabel('Easting (m)', fontdict={'size': 9, 'weight': 'bold'})
        ax.set_ylabel('Northing (m)', fontdict={'size': 9, 'weight': 'bold'})
        ax.set_title("Elevation and Stations Map")

        ax.scatter(self.station_locations.rel_east,
                   self.station_locations.rel_north, 
                   marker='v', c='b', s=2)
        ax.set_xlim((np.floor(self.station_locations.rel_east.min()) - 1000, 
                     np.ceil(self.station_locations.rel_east.max()) + 1000))
        ax.set_ylim((np.floor(self.station_locations.rel_north.min()) - 1000, 
                     np.ceil(self.station_locations.rel_north.max()) + 1000))


        plt.show()
        
    def plot_sealevel_resistivity(self):
        """
        create a quick pcolor plot of the resistivity at sea level with 
        stations, to check if we have stations in the sea
        
        """
        if self.res_model is None:
            print("Can't plot model, please read or create model file first")
            return        
        
        # index of sea level (zero level) in resistivity grid
        sli = mtcc.nearest_index(self.sea_level,self.grid_z)
        
        # make a figure
        plt.figure(figsize=(10,10))
        # plot the resistivity model (at sea level)
        plt.pcolormesh(self.grid_east,self.grid_north,self.res_model[:,:,sli],
                       vmin=1,vmax=1e4,norm=colors.LogNorm(),ec='0.5',lw=0.01,
                       cmap='bwr_r')
        
        # plot stations
        if self.station_locations is None:
            print("Can't plot stations, please read or create data file first")
        else:
            plt.plot(self.station_locations.rel_east,self.station_locations.rel_north,'.',color='k')
            
        # tidy up plot and make colorbar
        plt.gca().set_aspect(1)
        cbar=plt.colorbar(shrink=0.5)
        cbar.set_label('Resistivity, $\Omega$m')
        plt.xlabel('Grid East, relative (m)')
        plt.ylabel('Grid North, relative (m)')
        plt.title("Resistivity at sea level")
        
        
        
        
        
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
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

        self.save_path,self.model_fn, self.model_fn_basename = \
        mtfh.validate_save_file(savepath=self.save_path,
                                savefile=self.model_fn,
                                basename=self.model_fn_basename)

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
        with open(self.model_fn, 'w') as ifid:
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

            # not needed ifid.close()

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

        with open(self.model_fn, 'r') as ifid:
            ilines = ifid.readlines()

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
                print(iline)
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
        # use the grid centre from the model file
        shift_north = self.grid_center[0]# + self.nodes_north.sum() / 2
        shift_east = self.grid_center[1]# + self.nodes_east.sum() / 2
        shift_z = self.grid_center[2]

        # shift the grid.  if shift is + then that means the center is
        self.grid_north += shift_north
        self.grid_east += shift_east
        self.grid_z += shift_z

        # get cell size
        self.cell_size_east = stats.mode(self.nodes_east)[0][0]
        self.cell_size_north = stats.mode(self.nodes_north)[0][0]

        # get number of padding cells
        self.pad_east = np.where(self.nodes_east[0:int(self.nodes_east.size / 2)]
                                 != self.cell_size_east)[0].size
        self.pad_north = np.where(self.nodes_north[0:int(self.nodes_north.size / 2)]
                                  != self.cell_size_north)[0].size

    def read_ws_model_file(self, ws_model_fn):
        """
        reads in a WS3INV3D model file
        """

        ws_model_obj = ws.WSModel(ws_model_fn)
        ws_model_obj.read_model_file()

        # set similar attributes
        for ws_key in list(ws_model_obj.__dict__.keys()):
            for md_key in list(self.__dict__.keys()):
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
        if fn is not None:
            # if fn is a full path, convert to a file name
            fndir = os.path.basename(fn)
            if os.path.isdir(fndir):
                sg_basename = os.path.basename(fn)
            else:
                sg_basename = fn
        else:
            # create a basename if fn is None
            sg_basename = os.path.basename(self.model_fn).split('.')[0]
                
        self.save_path, fn, sg_basename = \
        mtfh.validate_save_file(savepath=self.save_path,
                                savefile=fn,
                                basename=sg_basename)

        if fn is None:
            fn = os.path.join(os.path.dirname(self.model_fn),
                              os.path.basename(self.model_fn).split('.')[0])

        # number of cells in the ModEM model
        nyin, nxin, nzin = np.array(self.res_model.shape) + 1
        
        
        gx,gy = mtmesh.rotate_mesh(self.grid_east[clip[0]:nxin - clip[0]],
                                   self.grid_north[clip[1]:nyin - clip[1]],
                                   origin[:2],self.mesh_rotation_angle)
        
        gz = -1.*self.grid_z[:nzin - clip[2]] - origin[2]
        
        gxm, gzm = np.meshgrid(gx, gz)
        gym, gzm = np.meshgrid(gy, gz)
        
        
        gxm = gxm.reshape(len(gz),len(gy),len(gx[0])).transpose(1,2,0)
        gym = gym.reshape(len(gz),len(gy),len(gx[0])).transpose(1,2,0)
        gzm = gzm.reshape(len(gz),len(gy),len(gx[0])).transpose(1,2,0)
        
        gridedges = (gxm,gym,gzm)

#        # get x, y and z positions
#        gridedges = [self.grid_east[clip[0]:nxin - clip[0]] + origin[0],
#                     self.grid_north[clip[1]:nyin - clip[1]] + origin[1],
#                     -1. * self.grid_z[:nzin - clip[2]] - origin[2]]
#        gridedges = np.meshgrid(*gridedges)

        # resistivity values, clipped to one smaller than grid edges
        resvals = self.res_model[clip[1]:nyin - clip[1] - 1,
                  clip[0]:nxin - clip[0] - 1, :nzin - clip[2] - 1]

        sgObj = mtgocad.Sgrid(resistivity=resvals, grid_xyz=gridedges,
                              fn=sg_basename, workdir=self.save_path)
        sgObj.write_sgrid_file()


    def read_gocad_sgrid_file(self, sgrid_header_file, air_resistivity=1e39, sea_resistivity=0.3,
                              sgrid_positive_up = True):
        """
        read a gocad sgrid file and put this info into a ModEM file.
        Note: can only deal with grids oriented N-S or E-W at this stage,
        with orthogonal coordinates

        """
        # read sgrid file
        sgObj = mtgocad.Sgrid()
        sgObj.read_sgrid_file(sgrid_header_file)
        self.sgObj = sgObj



        # get resistivity model values
        self.res_model = sgObj.resistivity

        # get nodes and grid locations
        grideast, gridnorth, gridz = [
            np.unique(sgObj.grid_xyz[i]) for i in range(3)]
        # check if sgrid is positive up and convert to positive down if it is
        # (ModEM grid is positive down)
        if sgrid_positive_up:
            gridz = -gridz
            
        gridz.sort()
        
        if np.all(np.array([len(gridnorth), len(grideast), len(gridz)]) - 1 == np.array(self.res_model.shape)):
            self.grid_east, self.grid_north, self.grid_z = grideast, gridnorth, gridz
        else:
            print("Cannot read sgrid, can't deal with non-orthogonal grids or grids not aligned N-S or E-W")
            return
        
        # check if we have a data object and if we do, is there a centre position
        # if not then assume it is the centre of the grid
        calculate_centre = True
        if self.data_obj is not None:
            if hasattr(self.data_obj, 'center_point'):
                if self.data_obj.center_point is not None:
                    centre = np.zeros(3)
                    centre[0] = self.data_obj.center_point['east'] 
                    centre[1] = self.data_obj.center_point['north']
                    calculate_centre = False
        # get relative grid locations
        if calculate_centre:
            print("Calculating center position")
            centre = np.zeros(3)
            centre[0] = (self.grid_east.max() + self.grid_east.min()) / 2.
            centre[1] = (self.grid_north.max() + self.grid_north.min()) / 2.
        centre[2] = self.grid_z[0]

        self.grid_east -= centre[0]
        self.grid_north -= centre[1]
        
        self.grid_center = np.array([self.grid_north[0],self.grid_east[0],self.grid_z[0]])       
        
        # get nodes
        # don't need to get nodes - as they are a property that auto-updates
#        self.nodes_east = self.grid_east[1:] - self.grid_east[:-1]
#        self.nodes_north = self.grid_north[1:] - self.grid_north[:-1]
#        self.nodes_z = self.grid_z[1:] - self.grid_z[:-1]

        
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


    def interpolate_elevation2(self, surfacefile=None, surface=None, get_surfacename=False,
                               method='nearest', fast=True):
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
        surface_epsg = epsg number of input surface, default is 4326 for lat/lon(wgs84)
        method = interpolation method. Default is 'nearest', if model grid is
        dense compared to surface points then choose 'linear' or 'cubic'

        """
        # initialise a dictionary to contain the surfaces
        if not hasattr(self, 'surface_dict'):
            self.surface_dict = {}

        # get centre position of model grid in real world coordinates
        x0, y0 = self.station_locations.center_point.east[0], self.station_locations.center_point.north[0]


        if self.mesh_rotation_angle is None:
            self.mesh_rotation_angle = 0
        
        xg,yg = mtmesh.rotate_mesh(self.grid_east,self.grid_north,
                                   [x0,y0],
                                   self.mesh_rotation_angle,
                                   return_centre = True)
        
        if surfacefile:
            elev_mg = mtmesh.interpolate_elevation_to_grid(
                    xg,yg, surfacefile=surfacefile, epsg=self.station_locations.model_epsg,
                    utm_zone=self.station_locations.model_utm_zone, method=method, fast=fast)
        elif surface:
            # Always use fast=False when reading from EDI data because
            #  we're already providing a subset of the grid.
            elev_mg = mtmesh.interpolate_elevation_to_grid(
                    xg, yg, surface=surface, epsg=self.station_locations.model_epsg,
                    utm_zone=self.station_locations.model_utm_zone,
                    method=method, fast=False)
        else:
            raise ValueError("'surfacefile' or 'surface' must be provided")

        print("Elevation data type and shape  *** ", type(elev_mg), elev_mg.shape, len(yg), len(xg))
        # <type 'numpy.ndarray'>  (65, 92), 65 92: it's 2D image with cell index as pixels
        # np.savetxt('E:/tmp/elev_mg.txt', elev_mg, fmt='%10.5f')

        # get a name for surface
        if get_surfacename:
            if surfacefile is not None:
                surfacename = os.path.basename(surfacefile)
            else:
                ii = 1
                surfacename = 'surface%01i' % ii
                while surfacename in list(self.surface_dict.keys()):
                    ii += 1
                    surfacename = 'surface%01i' % ii
            return elev_mg, surfacename
        else:
            return elev_mg


    def add_topography_from_data(self, data_object, interp_method='nearest', 
                                 air_resistivity=1e12, topography_buffer=None,
                                 airlayer_type='log_up'):
        """
        Wrapper around add_topography_to_model2 that allows creating
        a surface model from EDI data. The Data grid and station 
        elevations will be used to make a 'surface' tuple that will
        be passed to add_topography_to_model2 so a surface model
        can be interpolated from it.

        The surface tuple is of format (lon, lat, elev) containing
        station locations.

        Args:
            data_object (mtpy.modeling.ModEM.data.Data): A ModEm data
                object that has been filled with data from EDI files.
            interp_method (str, optional): Same as 
                add_topography_to_model2.
            air_resistivity (float, optional): Same as 
                add_topography_to_model2.
            topography_buffer (float): Same as 
                add_topography_to_model2.
            airlayer_type (str, optional): Same as 
                add_topography_to_model2.
        """
        lon = self.station_locations.lon
        lat = self.station_locations.lat
        elev = self.station_locations.elev
        surface = lon, lat, elev
        self.add_topography_to_model2(surface=surface, 
                                      interp_method=interp_method,
                                      air_resistivity=air_resistivity,
                                      topography_buffer=topography_buffer,
                                      airlayer_type=airlayer_type)


    def add_topography_to_model2(self, topographyfile=None, surface=None, 
                                 topographyarray=None, interp_method='nearest',
                                 air_resistivity=1e12, topography_buffer=None,
                                 airlayer_type = 'log_up', max_elev=None):
        """
        if air_layers is non-zero, will add topo: read in topograph file,
        make a surface model.
        
        Call project_stations_on_topography in the end, which will re-write 
        the .dat file.

        If n_airlayers is zero, then cannot add topo data, only bathymetry is needed.
        
        :param topographyfile: file containing topography (arcgis ascii grid)
        :param topographyarray: alternative to topographyfile - array of 
                                elevation values on model grid
        :param interp_method: interpolation method for topography,
                              'nearest', 'linear', or 'cubic'
        :param air_resistivity: resistivity value to assign to air
        :param topography_buffer: buffer around stations to calculate minimum
                                  and maximum topography value to use for 
                                  meshing
        :param airlayer_type: how to set air layer thickness - options are
                             'constant' for constant air layer thickness,
                             or 'log', for logarithmically increasing air 
                             layer thickness upward
        """
        # first, get surface data
        if topographyfile:
            self.surface_dict['topography'] = self.interpolate_elevation2(
                surfacefile=topographyfile, method=interp_method)
        elif surface:
            self.surface_dict['topography'] = self.interpolate_elevation2(
                surface=surface, method=interp_method)
        elif topographyarray:
            self.surface_dict['topography'] = topographyarray
        else:
            raise ValueError("'topographyfile', 'surface' or " +
                             "topographyarray must be provided")

        if self.n_air_layers is None or self.n_air_layers == 0:
            self._logger.warn("No air layers specified, so will not add air/topography !!!")
            self._logger.warn("Only bathymetry will be added below according to the topofile: sea-water low resistivity!!!")

        elif self.n_air_layers > 0:  # FZ: new logic, add equal blocksize air layers on top of the simple flat-earth grid
            # get grid centre
            gcx, gcy = [np.mean([arr[:-1], arr[1:]], axis=0) for arr in (self.grid_east, self.grid_north)]
            # get core cells
            if topography_buffer is None:
                topography_buffer = 5 * (self.cell_size_east ** 2 + self.cell_size_north ** 2) ** 0.5
            core_cells = mtmesh.get_station_buffer(gcx,
                                                   gcy,
                                                   self.station_locations.station_locations['rel_east'],
                                                   self.station_locations.station_locations['rel_north'],
                                                   buf=topography_buffer)
            topo_core = self.surface_dict['topography'][core_cells]
            topo_core_min = max(topo_core.min(),0)
            
            if airlayer_type == 'log_up':
                # log increasing airlayers, in reversed order
                new_air_nodes = mtmesh.make_log_increasing_array(self.z1_layer,
                                                                 topo_core.max() - topo_core_min,
                                                                 self.n_air_layers,
                                                                 increment_factor=0.999)[::-1]
            elif airlayer_type == 'log_down':
                # increase the number of layers
                self.n_layers += self.n_air_layers
                # make a new mesh
                n_layers = self.n_layers + self.n_air_layers
                self.nodes_z, z_grid = self.make_z_mesh_new(n_layers)

                # adjust level to topography min
                if max_elev is not None:
                    self.grid_z -= max_elev
                    ztops = np.where(self.surface_dict['topography'] > max_elev)
                    self.surface_dict['topography'][ztops] = max_elev
                else:
                    self.grid_z -= topo_core.max() 
                    
                                
            elif airlayer_type == 'constant':
                if max_elev is not None:
                    air_cell_thickness = np.ceil((max_elev - topo_core_min)/self.n_air_layers)
                else:
                    air_cell_thickness = np.ceil((topo_core.max() - topo_core_min)/self.n_air_layers)
                new_air_nodes = np.array([air_cell_thickness]*self.n_air_layers)

            if 'down' not in airlayer_type:
                 # sum to get grid cell locations
                new_airlayers = np.array([new_air_nodes[:ii].sum() 
                                          for ii in range(len(new_air_nodes) + 1)])
                # maximum topography cell on the grid
                topo_max_grid = topo_core_min + new_airlayers[-1]
                # round to nearest whole number and convert subtract the max elevation (so that sea level is at topo_core_min)
                new_airlayers = np.around(new_airlayers - topo_max_grid)
                # add new air layers, cut_off some tailing layers to preserve array size.
                self.grid_z = np.concatenate([new_airlayers[:-1], 
                                              self.grid_z + new_airlayers[-1]],
                                             axis=0)
                
            self._logger.debug("self.grid_z[0:2] {}".format(self.grid_z[0:2]))

        # update the z-centre as the top air layer
        self.grid_center[2] = self.grid_z[0]

        # update the resistivity model
        new_res_model = np.ones((self.nodes_north.size,
                                 self.nodes_east.size,
                                 self.nodes_z.size)) * self.res_initial_value

        if 'down' not in airlayer_type:
            new_res_model[:, :, self.n_air_layers:] = self.res_model
        
        self.res_model = new_res_model

        # assign topography
        top = np.zeros_like(self.surface_dict['topography']) + self.grid_z[0]
        bottom = -self.surface_dict['topography']
        self.assign_resistivity_from_surfacedata(top,bottom, 
                                                 air_resistivity)
        # assign bathymetry
        self.assign_resistivity_from_surfacedata(np.zeros_like(top),
                                                 bottom,
                                                 0.3)

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

    def _get_xyzres(self,location_type,origin,model_epsg,model_utm_zone,clip):
        # try getting centre location info from file
        if type(origin) == str:
            try:
                origin = np.loadtxt(origin)
            except:
                print("Please provide origin as a list, array or tuple or as a valid filename containing this info")
                origin = [0,0]
        
        # reshape the data and get grid centres
        x,y,z = [np.mean([arr[1:], arr[:-1]],axis=0) for arr in \
                [self.grid_east + origin[0], 
                 self.grid_north + origin[1], self.grid_z]]
        xsize, ysize = x.shape[0],y.shape[0]
        x, y, z = np.meshgrid(x[clip[0]:xsize-clip[0]],y[clip[1]:ysize-clip[1]],z)
        
        # set format for saving data
        fmt = ['%.1f','%.1f','%.3e']
        

        # convert to lat/long if needed
        if location_type == 'LL':
            if np.any(origin) == 0:
                print("Warning, origin coordinates provided as zero, output lat/long are likely to be incorrect")
            # project using epsg_project as preference as it is faster, but if pyproj not installed, use gdal
            try:
                import pyproj
                xp,yp = gis_tools.epsg_project(x,y,model_epsg,4326)
            except ImportError:
                xp,yp = np.zeros_like(x),np.zeros_like(y)
                for i in range(len(x)):
                    yp[i],xp[i] = gis_tools.project_point_utm2ll(x[i],y[i],model_utm_zone,epsg=model_epsg)
            # update format to accommodate lat/lon
            fmt[:2] = ['%.6f','%.6f']
        else:
            xp, yp = x, y
            
            
        resvals = self.res_model[clip[1]:ysize-clip[1],clip[0]:xsize-clip[0]]
            
        return xp, yp, z, resvals, fmt
    
    
    def write_xyzres(self,savefile=None,location_type='EN',origin=[0,0],model_epsg=None,log_res=False,model_utm_zone=None,clip=[0,0]):
        """
        save a model file as a space delimited x y z res file
    
        """
        xp, yp, z, resvals, fmt = self._get_xyzres(location_type,origin,model_epsg,model_utm_zone,clip)
        fmt.insert(2, '%.1f')
        xp, yp, z, resvals = xp.flatten(), yp.flatten(), z.flatten(), resvals.flatten()
            
        np.savetxt(savefile,np.vstack([xp,yp,z,resvals]).T,fmt=fmt)
            

    def write_xyres(self,savepath=None,location_type='EN',origin=[0,0],model_epsg=None,depth_index='all',
                    outfile_basename='DepthSlice',log_res=False,model_utm_zone=None,clip=[0,0]):
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
        clip = number of cells to clip on each of the east/west and north/south edges
        
        """
        if savepath is None:
            savepath = self.save_path
        
        # make a directory to save the files
        savepath = os.path.join(savepath,outfile_basename)
        if not os.path.exists(savepath):
            os.mkdir(savepath)
        

        xp, yp, z, resvals, fmt = self._get_xyzres(location_type,origin,model_epsg,model_utm_zone,clip)
        xp = xp[:,:,0].flatten()
        yp = yp[:,:,0].flatten()
        

        # make depth indices into a list
        if depth_index == 'all':
            depthindices = list(range(z.shape[2]))
        elif np.iterable(depth_index):
            depthindices = np.array(depth_index).astype(int)
        else:
            depthindices = [depth_index]
        

        
        for k in depthindices:
            fname = os.path.join(savepath,outfile_basename+'_%1im.xyz'%z[k])
            
            # get relevant depth slice
            vals = resvals[:,:,k].flatten()

            if log_res:
                vals = np.log10(vals)
                fmt[-1] = '%.3f'
            data = np.vstack([xp,yp,vals]).T

            np.savetxt(fname,data,fmt=fmt)

