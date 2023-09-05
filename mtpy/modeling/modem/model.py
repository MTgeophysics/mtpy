"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch
# revised by JP 2021 updating functionality and updating docs

"""

# =============================================================================
# Imports
# =============================================================================
from pathlib import Path
import numpy as np
from scipy import stats as stats
from loguru import logger

import mtpy.modeling.mesh_tools as mtmesh
import mtpy.modeling.gocad as mtgocad
import mtpy.utils.calculator as mtcc
import mtpy.utils.filehandling as mtfh

from .exception import ModelError
from mtpy.utils.gis_tools import project_point
from mtpy.modeling.plots.plot_mesh import PlotMesh
from mtpy.core.mt_location import MTLocation

from pyevtk.hl import gridToVTK

# =============================================================================


class Model:
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


    """

    def __init__(self, station_locations=None, center_point=None, **kwargs):
        self._logger = logger

        self.station_locations = None
        self.center_point = MTLocation()

        if station_locations is not None:
            self.station_locations = station_locations

        if center_point is not None:
            self.center_point = center_point
            self.model_epsg = self.center_point.utm_epsg

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
        self.z_layer_rounding = 0
        self.z_target_depth = 50000
        self.z_bottom = 300000

        # number of vertical layers
        self.n_layers = 30

        # number of air layers
        self.n_air_layers = 0
        # sea level in grid_z coordinates. Auto adjusts when topography read in?
        self.sea_level = 0.0

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
        self.grid_z = kwargs.pop("grid_z", None)
        if self.grid_z is not None:
            self.n_layers = len(self.grid_z)
            self.z_mesh_method = "custom"
        else:
            self.z_mesh_method = "new"
        if "z_mesh_method" in list(kwargs.keys()):
            self.z_mesh_method = kwargs["z_mesh_method"]

        # method to use to create padding
        self.pad_method = "extent1"

        self.grid_center = None
        self.surface_dict = {}

        # resistivity model
        self.res_initial_value = 100.0
        self.res_model = None

        # initial file stuff
        self.save_path = Path().cwd()
        self.model_fn_basename = "ModEM_Model_File.rho"

        self.title = "Model File written by MTpy.modeling.modem"
        self.res_scale = "loge"

        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self._logger.warning(
                    f"Argument {key}={value} is not supportted thus not been set."
                )

    def __str__(self):
        lines = ["ModEM Model Object:", "-" * 20]
        # --> print out useful information
        try:
            lines.append(
                f"\tNumber of stations = {len(self.station_locations.station)}"
            )
        except AttributeError:
            lines.append("\tNumber of stations = unknown")

        lines.append("\tMesh Parameter: ")
        lines.append(f"\t\tcell_size_east:    {self.cell_size_east}")
        lines.append(f"\t\tcell_size_north:   {self.cell_size_north}")
        lines.append(f"\t\tpad_east:          {self.pad_east}")
        lines.append(f"\t\tpad_north:         {self.pad_north}")
        lines.append(f"\t\tpad_num:           {self.pad_num}")
        lines.append(f"\t\tz1_layer:          {self.z1_layer}")
        lines.append(f"\t\tz_target_depth:    {self.z_target_depth}")
        lines.append(f"\t\tn_layers:          {self.n_layers}")
        lines.append(f"\t\tres_initial_value: {self.res_initial_value}")
        lines.append("\tDimensions: ")
        lines.append(f"\t\te-w: {self.grid_east.size}")
        lines.append(f"\t\tn-s: {self.grid_north.size}")
        lines.append(f"\t\tz:   {self.grid_z.size} (without 7 air layers)")
        lines.append("\tExtensions: ")
        lines.append(f"\t\te-w:  {self.nodes_east.__abs__().sum():.1f} (m)")
        lines.append(f"\t\tn-s:  {self.nodes_north.__abs__().sum():.1f} (m)")
        lines.append(f"\t\t0-z:  {self.nodes_z.__abs__().sum():.1f} (m)")
        if self.mesh_rotation_angle != 0:
            lines.append(
                f"\tStations rotated by: {self.mesh_rotation_angle:.1f} deg clockwise positive from N"
            )

            lines.append(
                " ** Note ModEM does not accommodate mesh rotations, it assumes"
            )
            lines.append("    all coordinates are aligned to geographic N, E")
            lines.append(
                "    therefore rotating the stations will have a similar effect"
            )
            lines.append("    as rotating the mesh.")
        lines.append("-" * 20)
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def save_path(self):
        return self._save_path

    @save_path.setter
    def save_path(self, save_path):
        if save_path is None:
            self._save_path = Path().cwd()
        else:
            self._save_path = Path(save_path)

        if not self._save_path.exists():
            self._save_path.mkdir()

    @property
    def model_fn(self):
        return self.save_path.joinpath(self.model_fn_basename)

    @model_fn.setter
    def model_fn(self, filename):
        if filename is not None:
            filename = Path(filename)
            self.save_path = filename.parent
            self.model_fn_basename = filename.name

    @property
    def model_epsg(self):
        return self.center_point.utm_epsg

    @model_epsg.setter
    def model_epsg(self, value):
        self.center_point.utm_epsg = value

    # --> make nodes and grid symbiotic so if you set one the other one
    #     gets set as well
    # Nodes East
    @property
    def nodes_east(self):
        if self.grid_east is not None:
            self._nodes_east = np.array(
                [
                    abs(self.grid_east[ii + 1] - self.grid_east[ii])
                    for ii in range(self.grid_east.size - 1)
                ]
            )
        return self._nodes_east

    @nodes_east.setter
    def nodes_east(self, nodes):
        nodes = np.array(nodes)
        self._nodes_east = nodes
        self.grid_east = np.array(
            [
                nodes[0:ii].sum() for ii in range(nodes.size + 1)
            ]  # -nodes.sum() / 2 +
        )  # + [shift])#[nodes.sum() / 2]

    # Nodes North
    @property
    def nodes_north(self):
        if self.grid_north is not None:
            self._nodes_north = np.array(
                [
                    abs(self.grid_north[ii + 1] - self.grid_north[ii])
                    for ii in range(self.grid_north.size - 1)
                ]
            )
        return self._nodes_north

    @nodes_north.setter
    def nodes_north(self, nodes):
        nodes = np.array(nodes)
        self._nodes_north = nodes
        self.grid_north = np.array(
            [
                nodes[0:ii].sum() for ii in range(nodes.size + 1)
            ]  # -nodes.sum() / 2 +
        )  # + [shift])#[nodes.sum() / 2]

    @property
    def nodes_z(self):
        if self.grid_z is not None:
            self._nodes_z = np.array(
                [
                    abs(self.grid_z[ii + 1] - self.grid_z[ii])
                    for ii in range(self.grid_z.size - 1)
                ]
            )

            return self._nodes_z

    @nodes_z.setter
    def nodes_z(self, nodes):
        nodes = np.array(nodes)
        self._nodes_z = nodes
        self.grid_z = np.array(
            [nodes[0:ii].sum() for ii in range(nodes.size)] + [nodes.sum()]
        )

    # need some arrays for plotting that are the same length as the
    # resistivity model
    @property
    def plot_east(self):
        plot_east = np.array(
            [self.nodes_east[0:ii].sum() for ii in range(self.nodes_east.size)]
        )
        return plot_east - plot_east[-1] / 2.0

    @property
    def plot_north(self):
        plot_north = np.array(
            [
                self.nodes_north[0:ii].sum()
                for ii in range(self.nodes_north.size)
            ]
        )
        return plot_north - plot_north[-1] / 2.0

    @property
    def plot_z(self):
        return np.array(
            [self.nodes_z[0:ii].sum() for ii in range(self.nodes_z.size)]
        )

    def make_mesh(self, verbose=True):
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
        west = self.station_locations.model_east.min() - pad_width_east
        east = self.station_locations.model_east.max() + pad_width_east
        south = self.station_locations.model_north.min() - pad_width_north
        north = self.station_locations.model_north.max() + pad_width_north

        # round the numbers so they are easier to read
        west = np.round(west, -2)
        east = np.round(east, -2)
        south = np.round(south, -2)
        north = np.round(north, -2)

        # -------make a grid around the stations from the parameters above------

        # adjust the edges so we have a whole number of cells
        add_ew = ((east - west) % self.cell_size_east) / 2.0
        add_ns = ((north - south) % self.cell_size_north) / 2.0

        # --> make the inner grid first
        inner_east = np.arange(
            west + add_ew - self.cell_size_east,
            east - add_ew + 2 * self.cell_size_east,
            self.cell_size_east,
        )
        inner_north = np.arange(
            south + add_ns + self.cell_size_north,
            north - add_ns + 2 * self.cell_size_north,
            self.cell_size_north,
        )

        # compute padding cells
        # first validate ew_ext and ns_ext to ensure it is large enough
        if "extent" in self.pad_method:
            self._validate_extent(
                inner_east.min(),
                inner_east.max(),
                inner_north.min(),
                inner_north.max(),
            )

        if self.pad_method == "extent1":
            padding_east = mtmesh.get_padding_cells(
                self.cell_size_east,
                self.ew_ext / 2 - east,
                self.pad_east,
                self.pad_stretch_h,
            )
            padding_north = mtmesh.get_padding_cells(
                self.cell_size_north,
                self.ns_ext / 2 - north,
                self.pad_north,
                self.pad_stretch_h,
            )
        elif self.pad_method == "extent2":
            padding_east = mtmesh.get_padding_cells2(
                self.cell_size_east,
                inner_east[-1],
                self.ew_ext / 2.0,
                self.pad_east,
            )
            padding_north = mtmesh.get_padding_cells2(
                self.cell_size_north,
                inner_north[-1],
                self.ns_ext / 2.0,
                self.pad_north,
            )
        elif self.pad_method == "stretch":
            padding_east = mtmesh.get_padding_from_stretch(
                self.cell_size_east, self.pad_stretch_h, self.pad_east
            )
            padding_north = mtmesh.get_padding_from_stretch(
                self.cell_size_north, self.pad_stretch_h, self.pad_north
            )
        else:
            raise NameError(
                'Padding method "{}" is not supported'.format(self.pad_method)
            )

        # make the horizontal grid
        self.grid_east = np.append(
            np.append(-1 * padding_east[::-1] + inner_east.min(), inner_east),
            padding_east + inner_east.max(),
        )
        self.grid_north = np.append(
            np.append(
                -1 * padding_north[::-1] + inner_north.min(), inner_north
            ),
            padding_north + inner_north.max(),
        )

        # --> need to make sure none of the stations lie on the nodes
        for s_east in sorted(self.station_locations.model_east):
            try:
                node_index = np.where(
                    abs(s_east - self.grid_east) < 0.02 * self.cell_size_east
                )[0][0]
                if s_east - self.grid_east[node_index] > 0:
                    self.grid_east[node_index] -= 0.02 * self.cell_size_east
                elif s_east - self.grid_east[node_index] < 0:
                    self.grid_east[node_index] += 0.02 * self.cell_size_east
            except IndexError:
                continue

        # --> need to make sure none of the stations lie on the nodes
        for s_north in sorted(self.station_locations.model_north):
            try:
                node_index = np.where(
                    abs(s_north - self.grid_north)
                    < 0.02 * self.cell_size_north
                )[0][0]
                if s_north - self.grid_north[node_index] > 0:
                    self.grid_north[node_index] -= 0.02 * self.cell_size_north
                elif s_north - self.grid_north[node_index] < 0:
                    self.grid_north[node_index] += 0.02 * self.cell_size_north
            except IndexError:
                continue

        if self.z_mesh_method == "custom":
            if self.grid_z is None:
                self.z_mesh_method = "new"
                self._logger.warn(
                    "No grid_z provided, creating new z mesh using default method"
                )

        if self.z_mesh_method == "custom":
            self.nodes_z, z_grid = (
                self.grid_z[1:] - self.grid_z[:-1],
                self.grid_z,
            )
        elif self.z_mesh_method == "new":
            self.nodes_z, z_grid = self.make_z_mesh()
        else:
            raise NameError(
                'Z mesh method "{}" is not supported'.format(
                    self.z_mesh_method
                )
            )

        # compute grid center
        center_east = np.round(
            self.grid_east.min() - self.grid_east.mean(), -1
        )
        center_north = np.round(
            self.grid_north.min() - self.grid_north.mean(), -1
        )
        center_z = 0

        # this is the value to the lower left corner from the center.
        self.grid_center = np.array([center_north, center_east, center_z])

        # make the resistivity array
        self.res_model = np.zeros(
            (self.nodes_north.size, self.nodes_east.size, self.nodes_z.size)
        )
        self.res_model[:, :, :] = self.res_initial_value

        # --> print out useful information
        if verbose:
            print(self.__str__())

    def make_z_mesh(self, n_layers=None):
        """
        new version of make_z_mesh. make_z_mesh and M
        """
        n_layers = self.n_layers if n_layers is None else n_layers

        # --> make depth grid
        # if n_airlayers < 0; set to 0
        log_z = mtcc.make_log_increasing_array(
            self.z1_layer, self.z_target_depth, n_layers - self.pad_z
        )

        if self.z_layer_rounding is not None:
            z_nodes = np.around(log_z, decimals=self.z_layer_rounding)
        else:
            # round any values less than 100 to the same s.f. as z1_layer
            z_nodes = np.around(
                log_z[log_z < 100],
                decimals=-int(np.floor(np.log10(self.z1_layer))),
            )
            # round any values greater than or equal to 100 to the nearest 100
            z_nodes = np.append(
                z_nodes, np.around(log_z[log_z >= 100], decimals=-2)
            )

        # index of top of padding
        # itp = len(z_nodes) - 1

        # padding cells in the vertical direction
        z_0 = float(z_nodes[-1])
        for ii in range(1, self.pad_z + 1):
            pad_d = np.round(z_0 * self.pad_stretch_v**ii, -2)
            z_nodes = np.append(z_nodes, pad_d)
        # add air layers and define ground surface level.
        # initial layer thickness is same as z1_layer
        # z_nodes = np.hstack([[z1_layer] * n_air, z_nodes])

        # make an array of absolute values
        z_grid = np.array(
            [z_nodes[:ii].sum() for ii in range(z_nodes.shape[0] + 1)]
        )

        return z_nodes, z_grid

    def add_layers_to_mesh(
        self, n_add_layers=None, layer_thickness=None, where="top"
    ):
        """
        Function to add constant thickness layers to the top or bottom of mesh.
        Note: It is assumed these layers are added before the topography. If
        you want to add topography layers, use function add_topography_to_model

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
            add_layers = np.insert(np.cumsum(layer_thickness), 0, 0)[:-1]
            layer_thickness = layer_thickness[-1]

            if n_add_layers != len(add_layers):
                self._logger.warn(
                    "Updating number of layers to reflect the length of the layer thickness array"
                )
            n_add_layers = len(add_layers)
        else:
            add_layers = np.arange(
                0, n_add_layers * layer_thickness, layer_thickness
            )

        # create a new z grid
        self.grid_z = np.hstack(
            [add_layers, self.grid_z + add_layers[-1] + layer_thickness]
        )

        # update the number of layers
        self.n_layers = len(self.grid_z) - 1

        # add the extra layer to the res model
        self.res_model = np.vstack(
            [self.res_model[:, :, :n_add_layers].T, self.res_model.T]
        ).T

    def assign_resistivity_from_surface_data(
        self, top_surface, bottom_surface, resistivity_value
    ):
        """
        assign resistivity value to all points above or below a surface
        requires the surface_dict attribute to exist and contain data for
        surface key (can get this information from ascii file using
        project_surface)

        **inputs**
        surface_name = name of surface (must correspond to key in surface_dict)
        resistivity_value = value to assign
        where = 'above' or 'below' - assign resistivity above or below the
                surface
        """

        # FZ: should ref-define the self.res_model if its shape has changed after topo air layer are added

        gcz = np.mean([self.grid_z[:-1], self.grid_z[1:]], axis=0)

        self._logger.debug(
            "gcz is the cells centre coordinates: %s, %s" % (len(gcz), gcz)
        )

        # assign resistivity value
        for j in range(len(self.res_model)):
            for i in range(len(self.res_model[j])):
                ii = np.where(
                    (gcz > top_surface[j, i]) & (gcz <= bottom_surface[j, i])
                )[0]
                self.res_model[j, i, ii] = resistivity_value

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
                          to save_path/model_fn_basename

            **model_fn_basename** : string
                                    basename to save file to
                                    *default* is ModEM_Model.ws
                                    file is saved at save_path/model_fn_basename

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

        # get resistivity model
        if self.res_model is None:
            self.res_model = np.zeros(
                (
                    self.nodes_north.size,
                    self.nodes_east.size,
                    self.nodes_z.size,
                )
            )
            self.res_model[:, :, :] = self.res_initial_value

        elif type(self.res_model) in [float, int]:
            self.res_initial_value = self.res_model
            self.res_model = np.zeros(
                (
                    self.nodes_north.size,
                    self.nodes_east.size,
                    self.nodes_z.size,
                )
            )
            self.res_model[:, :, :] = self.res_initial_value

        # --> write file
        with open(self.model_fn, "w") as ifid:
            ifid.write("# {0}\n".format(self.title.upper()))
            ifid.write(
                "{0:>5}{1:>5}{2:>5}{3:>5} {4}\n".format(
                    self.nodes_north.size,
                    self.nodes_east.size,
                    self.nodes_z.size,
                    0,
                    self.res_scale.upper(),
                )
            )

            # write S --> N node block
            for ii, nnode in enumerate(self.nodes_north):
                ifid.write("{0:>12.3f}".format(abs(nnode)))

            ifid.write("\n")

            # write W --> E node block
            for jj, enode in enumerate(self.nodes_east):
                ifid.write("{0:>12.3f}".format(abs(enode)))
            ifid.write("\n")

            # write top --> bottom node block
            for kk, zz in enumerate(self.nodes_z):
                ifid.write("{0:>12.3f}".format(abs(zz)))
            ifid.write("\n")

            # write the resistivity in log e format
            if self.res_scale.lower() == "loge":
                write_res_model = np.log(self.res_model[::-1, :, :])
            elif (
                self.res_scale.lower() == "log"
                or self.res_scale.lower() == "log10"
            ):
                write_res_model = np.log10(self.res_model[::-1, :, :])
            elif self.res_scale.lower() == "linear":
                write_res_model = self.res_model[::-1, :, :]
            else:
                raise ModelError(
                    'resistivity scale "{}" is not supported.'.format(
                        self.res_scale
                    )
                )

            # write out the layers from resmodel
            for zz in range(self.nodes_z.size):
                ifid.write("\n")
                for ee in range(self.nodes_east.size):
                    for nn in range(self.nodes_north.size):
                        ifid.write(
                            "{0:>13.5E}".format(write_res_model[nn, ee, zz])
                        )
                    ifid.write("\n")

            if self.grid_center is None:
                # compute grid center
                center_east = -self.nodes_east.__abs__().sum() / 2
                center_north = -self.nodes_north.__abs__().sum() / 2
                center_z = 0
                self.grid_center = np.array(
                    [center_north, center_east, center_z]
                )

            ifid.write(
                "\n{0:>16.3f}{1:>16.3f}{2:>16.3f}\n".format(
                    self.grid_center[0],
                    self.grid_center[1],
                    self.grid_center[2],
                )
            )

            if self.mesh_rotation_angle is None:
                ifid.write("{0:>9.3f}\n".format(0))
            else:
                ifid.write("{0:>9.3f}\n".format(self.mesh_rotation_angle))

            # not needed ifid.close()

        self._logger.info("Wrote file to: {0}".format(self.model_fn))

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
            raise ModelError("model_fn is None, input a model file name")

        if not self.model_fn.exists():
            raise ModelError(f"Cannot find {self.model_fn}, check path")

        with open(self.model_fn, "r") as ifid:
            ilines = ifid.readlines()

        self.title = ilines[0].strip()

        # get size of dimensions, remembering that x is N-S, y is E-W, z is + down
        nsize = ilines[1].strip().split()
        n_north = int(nsize[0])
        n_east = int(nsize[1])
        n_z = int(nsize[2])
        log_yn = nsize[4]

        # get nodes
        self.nodes_north = np.array(
            [float(nn) for nn in ilines[2].strip().split()]
        )
        self.nodes_east = np.array(
            [float(nn) for nn in ilines[3].strip().split()]
        )
        self.nodes_z = np.array(
            [float(nn) for nn in ilines[4].strip().split()]
        )

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
                    self.grid_center = np.array(ilist, dtype=float)
                # rotation angle
                elif len(ilist) == 1:
                    self.mesh_rotation_angle = float(ilist[0])
                else:
                    pass

        # --> make sure the resistivity units are in linear Ohm-m
        if log_yn.lower() == "loge":
            self.res_model = np.e**self.res_model
        elif log_yn.lower() == "log" or log_yn.lower() == "log10":
            self.res_model = 10**self.res_model

        # center the grids
        if self.grid_center is None:
            self.grid_center = np.array(
                [-self.nodes_north.sum() / 2, -self.nodes_east.sum() / 2, 0.0]
            )

        # need to shift the grid if the center is not symmetric
        # use the grid centre from the model file
        shift_north = self.grid_center[0]  # + self.nodes_north.sum() / 2
        shift_east = self.grid_center[1]  # + self.nodes_east.sum() / 2
        shift_z = self.grid_center[2]

        # shift the grid.  if shift is + then that means the center is
        self.grid_north += shift_north
        self.grid_east += shift_east
        self.grid_z += shift_z

        # get cell size
        self.cell_size_east = stats.mode(self.nodes_east)[0][0]
        self.cell_size_north = stats.mode(self.nodes_north)[0][0]

        # get number of padding cells
        self.pad_east = np.where(
            self.nodes_east[0 : int(self.nodes_east.size / 2)]
            != self.cell_size_east
        )[0].size
        self.pad_north = np.where(
            self.nodes_north[0 : int(self.nodes_north.size / 2)]
            != self.cell_size_north
        )[0].size

    def plot_mesh(self, **kwargs):
        """
        Plot model mesh

        :param plot_topography: DESCRIPTION, defaults to False
        :type plot_topography: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if "topography" in self.surface_dict.keys():
            kwargs["plot_topography"] = True
        return PlotMesh(self, **kwargs)

    @property
    def model_parameters(self):
        """
        get important model parameters to write to a file for documentation
        later.


        """

        parameter_list = [
            "cell_size_east",
            "cell_size_north",
            "ew_ext",
            "ns_ext",
            "pad_east",
            "pad_north",
            "pad_z",
            "pad_num",
            "z1_layer",
            "z_target_depth",
            "z_bottom",
            "mesh_rotation_angle",
            "res_initial_value",
            "save_path",
        ]

        parameter_dict = {}
        for parameter in parameter_list:
            key = "model.{0}".format(parameter)
            parameter_dict[key] = getattr(self, parameter)

        parameter_dict["model.size"] = self.res_model.shape

        return parameter_dict

    def write_gocad_sgrid_file(
        self, fn=None, origin=[0, 0, 0], clip=0, no_data_value=-99999
    ):
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
            fn = Path(fn)
            # if fn is a full path, convert to a file name
            fndir = fn.parent
            if fndir.is_dir():
                sg_basename = fn.name
            else:
                sg_basename = fn
        else:
            # create a basename if fn is None
            sg_basename = self.model_fn.stem

        self.save_path, fn, sg_basename = mtfh.validate_save_file(
            save_path=self.save_path, savefile=fn, basename=sg_basename
        )

        # number of cells in the ModEM model
        nyin, nxin, nzin = np.array(self.res_model.shape) + 1

        gx, gy = mtmesh.rotate_mesh(
            self.grid_east[clip[0] : nxin - clip[0]],
            self.grid_north[clip[1] : nyin - clip[1]],
            origin[:2],
            self.mesh_rotation_angle,
        )

        gz = -1.0 * self.grid_z[: nzin - clip[2]] - origin[2]

        gxm, gzm = np.meshgrid(gx, gz)
        gym, gzm = np.meshgrid(gy, gz)

        gxm = gxm.reshape(len(gz), len(gy), len(gx[0])).transpose(1, 2, 0)
        gym = gym.reshape(len(gz), len(gy), len(gx[0])).transpose(1, 2, 0)
        gzm = gzm.reshape(len(gz), len(gy), len(gx[0])).transpose(1, 2, 0)

        gridedges = (gxm, gym, gzm)

        # resistivity values, clipped to one smaller than grid edges
        resvals = self.res_model[
            clip[1] : nyin - clip[1] - 1,
            clip[0] : nxin - clip[0] - 1,
            : nzin - clip[2] - 1,
        ]

        sg_obj = mtgocad.Sgrid(
            resistivity=resvals,
            grid_xyz=gridedges,
            fn=sg_basename,
            workdir=self.save_path,
        )
        sg_obj.write_sgrid_file()

    def read_gocad_sgrid_file(
        self,
        sgrid_header_file,
        air_resistivity=1e39,
        sea_resistivity=0.3,
        sgrid_positive_up=True,
    ):
        """
        read a gocad sgrid file and put this info into a ModEM file.
        Note: can only deal with grids oriented N-S or E-W at this stage,
        with orthogonal coordinates

        """
        # read sgrid file
        sg_obj = mtgocad.Sgrid()
        sg_obj.read_sgrid_file(sgrid_header_file)
        self.sg_obj = sg_obj

        # get resistivity model values
        self.res_model = sg_obj.resistivity

        # get nodes and grid locations
        grideast, gridnorth, gridz = [
            np.unique(sg_obj.grid_xyz[i]) for i in range(3)
        ]
        # check if sgrid is positive up and convert to positive down if it is
        # (ModEM grid is positive down)
        if sgrid_positive_up:
            gridz = -gridz

        gridz.sort()

        if np.all(
            np.array([len(gridnorth), len(grideast), len(gridz)]) - 1
            == np.array(self.res_model.shape)
        ):
            self.grid_east, self.grid_north, self.grid_z = (
                grideast,
                gridnorth,
                gridz,
            )
        else:
            print(
                "Cannot read sgrid, can't deal with non-orthogonal grids or grids not aligned N-S or E-W"
            )
            return

        # check if we have a data object and if we do, is there a centre position
        # if not then assume it is the centre of the grid
        calculate_centre = True
        if self.data_obj is not None:
            if hasattr(self.data_obj, "center_point"):
                if self.data_obj.center_point is not None:
                    centre = np.zeros(3)
                    centre[0] = self.data_obj.center_point["east"]
                    centre[1] = self.data_obj.center_point["north"]
                    calculate_centre = False
        # get relative grid locations
        if calculate_centre:
            print("Calculating center position")
            centre = np.zeros(3)
            centre[0] = (self.grid_east.max() + self.grid_east.min()) / 2.0
            centre[1] = (self.grid_north.max() + self.grid_north.min()) / 2.0
        centre[2] = self.grid_z[0]

        self.grid_east -= centre[0]
        self.grid_north -= centre[1]

        self.grid_center = np.array(
            [self.grid_north[0], self.grid_east[0], self.grid_z[0]]
        )

        self.z1_layer = self.nodes_z[0]
        #        self.z_target_depth = None
        self.z_bottom = self.nodes_z[-1]

        # number of vertical layers
        self.n_layers = len(self.grid_z) - 1

        # number of air layers
        self.n_airlayers = sum(
            np.amax(self.res_model, axis=(0, 1)) > 0.9 * air_resistivity
        )

        # sea level in grid_z coordinates, calculate and adjust centre
        self.sea_level = self.grid_z[self.n_airlayers]

    def interpolate_elevation(
        self,
        surface_file=None,
        surface=None,
        get_surface_name=False,
        method="nearest",
        fast=True,
        shift_north=0,
        shift_east=0,
    ):
        """
        project a surface to the model grid and add resulting elevation data
        to a dictionary called surface_dict. Assumes the surface is in lat/long
        coordinates (wgs84)

        **returns**
        nothing returned, but surface data are added to surface_dict under
        the key given by surface_name.

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
        if not hasattr(self, "surface_dict"):
            self.surface_dict = {}

        # get centre position of model grid in real world coordinates
        x0, y0 = (
            self.center_point.east + shift_east,
            self.center_point.north + shift_north,
        )

        if self.mesh_rotation_angle is None:
            self.mesh_rotation_angle = 0

        xg, yg = mtmesh.rotate_mesh(
            self.grid_east,
            self.grid_north,
            [x0, y0],
            self.mesh_rotation_angle,
            return_centre=True,
        )
        if surface_file:
            elev_mg = mtmesh.interpolate_elevation_to_grid(
                xg,
                yg,
                surface_file=surface_file,
                utm_epsg=self.model_epsg,
                datum_epsg=self.center_point.datum_epsg,
                method=method,
                fast=fast,
            )
        elif surface:
            # Always use fast=False when reading from EDI data because
            #  we're already providing a subset of the grid.
            elev_mg = mtmesh.interpolate_elevation_to_grid(
                xg,
                yg,
                surface=surface,
                utm_epsg=self.model_epsg,
                datum_epsg=self.center_point.datum_epsg,
                method=method,
                fast=False,
            )
        else:
            raise ValueError("'surface_file' or 'surface' must be provided")

        # get a name for surface
        if get_surface_name:
            if surface_file is not None:
                surface_file = Path(surface_file)
                surface_name = surface_file.name
            else:
                ii = 1
                surface_name = "surface%01i" % ii
                while surface_name in list(self.surface_dict.keys()):
                    ii += 1
                    surface_name = "surface%01i" % ii
            return elev_mg, surface_name
        else:
            return elev_mg

    def add_topography_from_data(
        self,
        interp_method="nearest",
        air_resistivity=1e12,
        topography_buffer=None,
        airlayer_type="log_up",
    ):
        """
        Wrapper around add_topography_to_model that allows creating
        a surface model from EDI data. The Data grid and station
        elevations will be used to make a 'surface' tuple that will
        be passed to add_topography_to_model so a surface model
        can be interpolated from it.

        The surface tuple is of format (lon, lat, elev) containing
        station locations.

        Args:
            data_object (mtpy.modeling.ModEM.data.Data): A ModEm data
                object that has been filled with data from EDI files.
            interp_method (str, optional): Same as
                add_topography_to_model.
            air_resistivity (float, optional): Same as
                add_topography_to_model.
            topography_buffer (float): Same as
                add_topography_to_model.
            airlayer_type (str, optional): Same as
                add_topography_to_model.
        """
        lon = self.station_locations.longitude.to_numpy()
        lat = self.station_locations.latitude.to_numpy()
        elev = self.station_locations.elevation.to_numpy()
        surface = lon, lat, elev
        self.add_topography_to_model(
            surface=surface,
            interp_method=interp_method,
            air_resistivity=air_resistivity,
            topography_buffer=topography_buffer,
            airlayer_type=airlayer_type,
        )

    def add_topography_to_model(
        self,
        topography_file=None,
        surface=None,
        topography_array=None,
        interp_method="nearest",
        air_resistivity=1e12,
        topography_buffer=None,
        airlayer_type="log_up",
        max_elev=None,
        shift_east=0,
        shift_north=0,
    ):
        """
        if air_layers is non-zero, will add topo: read in topograph file,
        make a surface model.

        Call project_stations_on_topography in the end, which will re-write
        the .dat file.

        If n_airlayers is zero, then cannot add topo data, only bathymetry is needed.

        :param topography_file: file containing topography (arcgis ascii grid)
        :param topography_array: alternative to topography_file - array of
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
        if topography_file:
            self.surface_dict["topography"] = self.interpolate_elevation(
                surface_file=topography_file,
                method=interp_method,
                shift_east=shift_east,
                shift_north=shift_north,
            )
        elif surface:
            self.surface_dict["topography"] = self.interpolate_elevation(
                surface=surface,
                method=interp_method,
                shift_east=shift_east,
                shift_north=shift_north,
            )
        elif topography_array:
            self.surface_dict["topography"] = topography_array
        else:
            raise ValueError(
                "'topography_file', 'surface' or "
                + "topography_array must be provided"
            )

        if self.n_air_layers is None or self.n_air_layers == 0:
            self._logger.warn(
                "No air layers specified, so will not add air/topography !!!"
            )
            self._logger.warn(
                "Only bathymetry will be added below according to the topofile: sea-water low resistivity!!!"
            )

        elif (
            self.n_air_layers > 0
        ):  # FZ: new logic, add equal blocksize air layers on top of the simple flat-earth grid
            # get grid centre
            gcx, gcy = [
                np.mean([arr[:-1], arr[1:]], axis=0)
                for arr in (self.grid_east, self.grid_north)
            ]
            # get core cells
            if topography_buffer is None:
                topography_buffer = (
                    5
                    * (self.cell_size_east**2 + self.cell_size_north**2)
                    ** 0.5
                )
            core_cells = mtmesh.get_station_buffer(
                gcx,
                gcy,
                self.station_locations["model_east"],
                self.station_locations["model_north"],
                buf=topography_buffer,
            )
            topo_core = self.surface_dict["topography"][core_cells]
            topo_core_min = max(topo_core.min(), 0)

            if airlayer_type == "log_up":
                # log increasing airlayers, in reversed order
                new_air_nodes = mtmesh.make_log_increasing_array(
                    self.z1_layer,
                    topo_core.max() - topo_core_min,
                    self.n_air_layers,
                    increment_factor=0.999,
                )[::-1]
            elif airlayer_type == "log_down":
                # make a new mesh
                n_layers = self.n_layers + self.n_air_layers
                self.nodes_z, z_grid = self.make_z_mesh(n_layers)

                # adjust level to topography min
                if max_elev is not None:
                    self.grid_z -= max_elev
                    ztops = np.where(
                        self.surface_dict["topography"] > max_elev
                    )
                    self.surface_dict["topography"][ztops] = max_elev
                else:
                    self.grid_z -= topo_core.max()

            elif airlayer_type == "constant":
                if max_elev is not None:
                    air_cell_thickness = np.ceil(
                        (max_elev - topo_core_min) / self.n_air_layers
                    )
                else:
                    air_cell_thickness = np.ceil(
                        (topo_core.max() - topo_core_min) / self.n_air_layers
                    )
                new_air_nodes = np.array(
                    [air_cell_thickness] * self.n_air_layers
                )

            if "down" not in airlayer_type:
                # sum to get grid cell locations
                new_airlayers = np.array(
                    [
                        new_air_nodes[:ii].sum()
                        for ii in range(len(new_air_nodes) + 1)
                    ]
                )
                # maximum topography cell on the grid
                topo_max_grid = topo_core_min + new_airlayers[-1]
                # round to nearest whole number and convert subtract the max elevation (so that sea level is at topo_core_min)
                new_airlayers = np.around(new_airlayers - topo_max_grid)
                # add new air layers, cut_off some tailing layers to preserve array size.
                self.grid_z = np.concatenate(
                    [new_airlayers[:-1], self.grid_z + new_airlayers[-1]],
                    axis=0,
                )

            self._logger.debug("self.grid_z[0:2] {}".format(self.grid_z[0:2]))

        # update the z-centre as the top air layer
        self.grid_center[2] = self.grid_z[0]

        # update the resistivity model
        new_res_model = (
            np.ones(
                (
                    self.nodes_north.size,
                    self.nodes_east.size,
                    self.nodes_z.size,
                )
            )
            * self.res_initial_value
        )

        if "down" not in airlayer_type:
            new_res_model[:, :, self.n_air_layers :] = self.res_model

        self.res_model = new_res_model

        # assign topography
        top = np.zeros_like(self.surface_dict["topography"]) + self.grid_z[0]
        bottom = -self.surface_dict["topography"]
        self.assign_resistivity_from_surface_data(top, bottom, air_resistivity)
        # assign bathymetry
        self.assign_resistivity_from_surface_data(
            np.zeros_like(top), bottom, 0.3
        )

        return

    def _validate_extent(self, east, west, south, north, extent_ratio=2.0):
        """
        validate the provided ew_ext and ns_ext to make sure the model fits
        within these extents and allows enough space for padding according to
        the extent ratio provided. If not, then update ew_ext and ns_ext parameters

        """
        inner_ew_ext = west - east
        inner_ns_ext = north - south

        if self.ew_ext < extent_ratio * inner_ew_ext:
            self._logger.warn(
                "Provided or default ew_ext not sufficient to fit stations + padding, updating extent"
            )
            self.ew_ext = np.ceil(extent_ratio * inner_ew_ext)

        if self.ns_ext < extent_ratio * inner_ns_ext:
            self._logger.warn(
                "Provided or default ns_ext not sufficient to fit stations + padding, updating extent"
            )
            self.ns_ext = np.ceil(extent_ratio * inner_ns_ext)

    def _get_xyzres(self, location_type, origin, model_epsg, clip):
        # try getting centre location info from file
        if type(origin) == str:
            try:
                origin = np.loadtxt(origin)
            except:
                print(
                    "Please provide origin as a list, array or tuple or as a valid filename containing this info"
                )
                origin = [0, 0]

        # reshape the data and get grid centres
        x, y, z = [
            np.mean([arr[1:], arr[:-1]], axis=0)
            for arr in [
                self.grid_east + origin[0],
                self.grid_north + origin[1],
                self.grid_z,
            ]
        ]
        xsize, ysize = x.shape[0], y.shape[0]
        x, y, z = np.meshgrid(
            x[clip[0] : xsize - clip[0]], y[clip[1] : ysize - clip[1]], z
        )

        # set format for saving data
        fmt = ["%.1f", "%.1f", "%.3e"]

        # convert to lat/long if needed
        if location_type == "LL":
            if np.any(origin) == 0:
                print(
                    "Warning, origin coordinates provided as zero, output lat/long are likely to be incorrect"
                )
            # project using epsg_project as preference as it is faster, but if pyproj not installed, use gdal

            xp, yp = project_point(x, y, model_epsg, 4326)

            # update format to accommodate lat/lon
            fmt[:2] = ["%.6f", "%.6f"]
        else:
            xp, yp = x, y

        resvals = self.res_model[
            clip[1] : ysize - clip[1], clip[0] : xsize - clip[0]
        ]

        return xp, yp, z, resvals, fmt

    def write_xyzres(
        self,
        savefile=None,
        location_type="EN",
        origin=[0, 0],
        model_epsg=None,
        log_res=False,
        model_utm_zone=None,
        clip=[0, 0],
    ):
        """
        save a model file as a space delimited x y z res file

        """
        xp, yp, z, resvals, fmt = self._get_xyzres(
            location_type, origin, model_epsg, clip
        )
        fmt.insert(2, "%.1f")
        xp, yp, z, resvals = (
            xp.flatten(),
            yp.flatten(),
            z.flatten(),
            resvals.flatten(),
        )

        np.savetxt(savefile, np.vstack([xp, yp, z, resvals]).T, fmt=fmt)

    def write_xyres(
        self,
        save_path=None,
        location_type="EN",
        origin=[0, 0],
        model_epsg=None,
        depth_index="all",
        outfile_basename="DepthSlice",
        log_res=False,
        clip=[0, 0],
    ):
        """
        write files containing depth slice data (x, y, res for each depth)

        origin = x,y coordinate of zero point of ModEM_grid, or name of file
                 containing this info (full path or relative to model files)
        save_path = path to save to, default is the model object save path
        location_type = 'EN' or 'LL' xy points saved as eastings/northings or
                        longitude/latitude, if 'LL' need to also provide model_epsg
        model_epsg = epsg number that was used to project the model
        outfile_basename = string for basename for saving the depth slices.
        log_res = True/False - option to save resistivity values as log10
                               instead of linear
        clip = number of cells to clip on each of the east/west and north/south edges

        """
        if save_path is None:
            save_path = Path(self.save_path)
        else:
            save_path = Path(save_path)
        # make a directory to save the files
        save_path = save_path.joinpath(outfile_basename)
        if not save_path.exists():
            save_path.mkdir()

        xp, yp, z, resvals, fmt = self._get_xyzres(
            location_type, origin, model_epsg, clip
        )
        xp = xp[:, :, 0].flatten()
        yp = yp[:, :, 0].flatten()

        # make depth indices into a list
        if depth_index == "all":
            depthindices = list(range(z.shape[2]))
        elif np.iterable(depth_index):
            depthindices = np.array(depth_index).astype(int)
        else:
            depthindices = [depth_index]

        for k in depthindices:
            fname = save_path.joinpath(
                outfile_basename + "_%1im.xyz" % self.grid_z[k]
            )

            # get relevant depth slice
            vals = resvals[:, :, k].flatten()

            if log_res:
                vals = np.log10(vals)
                fmt[-1] = "%.3f"
            data = np.vstack([xp, yp, vals]).T

            np.savetxt(fname, data, fmt=fmt)

    def write_vtk_file(
        self,
        vtk_save_path=None,
        vtk_fn_basename="ModEM_model_res",
        shift_east=0,
        shift_north=0,
        shift_z=0,
        units="km",
        coordinate_system="nez+",
        label="resistivity",
    ):
        """
        Write a VTK file to plot in 3D rendering programs like Paraview

        :param vtk_save_path: directory to save vtk file to, defaults to None
        :type vtk_save_path: string or Path, optional
        :param vtk_fn_basename: filename basename of vtk file, note that .vtr
        extension is automatically added, defaults to "ModEM_stations"
        :type vtk_fn_basename: string, optional
        :type geographic: boolean, optional
        :param shift_east: shift in east directions in meters, defaults to 0
        :type shift_east: float, optional
        :param shift_north: shift in north direction in meters, defaults to 0
        :type shift_north: float, optional
        :param shift_z: shift in elevation + down in meters, defaults to 0
        :type shift_z: float, optional
        :param units: Units of the spatial grid [ km | m | ft ], defaults to "km"
        :type units: string, optional
        :type : string
        :param coordinate_system: coordinate system for the station, either the
        normal MT right-hand coordinate system with z+ down or the sinister
        z- down [ nez+ | enz- ], defaults to nez+
        :return: full path to VTK file
        :rtype: Path

        Write VTK file
        >>> model.write_vtk_file(vtk_fn_basename="modem_model")

        Write VTK file in geographic coordinates with z+ up
        >>> model.write_vtk_station_file(vtk_fn_basename="modem_model",
        >>> ...                          coordinate_system='enz-')
        """

        if isinstance(units, str):
            if units.lower() == "km":
                scale = 1.0 / 1000.00
            elif units.lower() == "m":
                scale = 1.0
            elif units.lower() == "ft":
                scale = 3.2808
        elif isinstance(units, (int, float)):
            scale = units

        if vtk_save_path is None:
            vtk_fn = self.save_path.joinpath(vtk_fn_basename)
        else:
            vtk_fn = Path(vtk_save_path).joinpath(vtk_fn_basename)

        # use cellData, this makes the grid properly as grid is n+1
        if coordinate_system == "nez+":
            vtk_x = (self.grid_north + shift_north) * scale
            vtk_y = (self.grid_east + shift_east) * scale
            vtk_z = (self.grid_z + shift_z) * scale
            cell_data = {label: self.res_model}

        elif coordinate_system == "enz-":
            vtk_y = (self.grid_north + shift_north) * scale
            vtk_x = (self.grid_east + shift_east) * scale
            vtk_z = -1 * (self.grid_z + shift_z) * scale
            cell_data = {label: np.rot90(self.res_model)}

        gridToVTK(vtk_fn.as_posix(), vtk_x, vtk_y, vtk_z, cellData=cell_data)

        self._logger.info("Wrote model file to {}".format(vtk_fn))

    def write_geosoft_xyz(
        self,
        save_fn,
        c_east=0,
        c_north=0,
        c_z=0,
        pad_north=0,
        pad_east=0,
        pad_z=0,
    ):
        """
        Write an XYZ file readable by Geosoft

        All input units are in meters.

        :param save_fn: full path to save file to
        :type save_fn: string or Path
        :param c_east: center point in the east direction, defaults to 0
        :type c_east: float, optional
        :param c_north: center point in the north direction, defaults to 0
        :type c_north: float, optional
        :param c_z: center point elevation, defaults to 0
        :type c_z: float, optional
        :param pad_north: number of cells to cut from the north-south edges, defaults to 0
        :type pad_north: int, optional
        :param pad_east: number of cells to cut from the east-west edges, defaults to 0
        :type pad_east: int, optional
        :param pad_z: number of cells to cut from the bottom, defaults to 0
        :type pad_z: int, optional


        """
        lines = [
            r"/ ------------------------------------------------------------------------------",
            r"/ XYZ  IMPORT [01/25/2021]",
            r"/ VOXEL   [.\electrical_resistivity.geosoft_voxel]",
            r"/ ------------------------------------------------------------------------------",
            r"/ X,Y,Z,Data",
        ]

        # --> write model xyz file
        for kk, zz in enumerate(self.grid_z[0:-pad_z]):
            for jj, yy in enumerate(self.grid_east[pad_east:-pad_east]):
                for ii, xx in enumerate(self.grid_north[pad_north:-pad_north]):
                    lines.append(
                        f"{yy + c_east:.3f} {xx + c_north:.3f} {-(zz + c_z):.3f} {self.res_model[ii, jj, kk]:.3f}"
                    )

        with open(save_fn, "w") as fid:
            fid.write("\n".join(lines))

    def write_out_file(
        self, save_fn, geographic_east, geographic_north, geographic_elevation
    ):
        """
        will write an .out file for LeapFrog.

        Note that y is assumed to be S --> N, e is assumed to be W --> E and
        z is positive upwards.  This means that index [0, 0, 0] is the
        southwest corner of the first layer.

        :param save_fn: full path to save file to
        :type save_fn: string or Path
        :param geographic_east: geographic center in easting (meters)
        :type geographic_east: float
        :param geographic_north: geographic center in northing (meters)
        :type geographic_north: float
        :param geographic_elevation: elevation of geographic center (meters)
        :type geographic_elevation: float
        :return: DESCRIPTION
        :rtype: TYPE

        """

        # get resistivity model
        if self.res_model is None:
            self.res_model = np.zeros(
                (
                    self.nodes_north.size,
                    self.nodes_east.size,
                    self.nodes_z.size,
                )
            )
            self.res_model[:, :, :] = self.res_initial_value

        elif type(self.res_model) in [float, int]:
            self.res_initial_value = self.res_model
            self.res_model = np.zeros(
                (
                    self.nodes_north.size,
                    self.nodes_east.size,
                    self.nodes_z.size,
                )
            )
            self.res_model[:, :, :] = self.res_initial_value

        shift_east = (
            geographic_east
            - (
                self.nodes_east[0]
                - self.nodes_east[1] / 2
                - self.grid_center[1] / 2
            )
        ) / 1000.0
        shift_north = (
            geographic_north
            + (
                self.nodes_north[0]
                - self.nodes_north[1] / 2
                - self.grid_center[0] / 2
            )
        ) / 1000.0

        shift_elevation = geographic_elevation / 1000.0

        # --> write file
        with open(save_fn, "w") as ifid:
            ifid.write("\n")
            ifid.write(
                "{0:>5}{1:>5}{2:>5}{3:>5} {4}\n".format(
                    self.nodes_east.size,
                    self.nodes_north.size,
                    self.nodes_z.size,
                    0,
                    "VAL",
                )
            )

            # write S --> N node block
            for ii, nnode in enumerate(self.nodes_east):
                ifid.write("{0:>12.3f}".format(abs(nnode)))

            ifid.write("\n")

            # write W --> E node block
            for jj, enode in enumerate(self.nodes_north):
                ifid.write("{0:>12.3f}".format(abs(enode)))
            ifid.write("\n")

            # write top --> bottom node block
            for kk, zz in enumerate(self.nodes_z):
                ifid.write("{0:>12.3f}".format(abs(zz)))
            ifid.write("\n")

            # write the resistivity in log e format
            write_res_model = self.res_model[::-1, :, :]

            # write out the layers from resmodel
            count = 1
            for zz in range(self.nodes_z.size):
                ifid.write(f"{count}\n")
                for nn in range(self.nodes_north.size):
                    for ee in range(self.nodes_east.size):
                        ifid.write(
                            "{0:>13.5E}".format(write_res_model[nn, ee, zz])
                        )
                    ifid.write("\n")
                count += 1

            # write footer
            ifid.write("\n")
            ifid.write("WINGLINK\n")
            ifid.write("  Project      (site name)\n")
            ifid.write("           1           1    (i j block numbers)\n")
            ifid.write(
                f"   {shift_east:.3f}       {shift_north:.3f}       (real world coordinates)\n"
            )
            ifid.write("  0.0000000E+00    (rotation)\n")
            ifid.write(f"   {shift_elevation:.3f}       (top elevation)\n")
            ifid.write("\n")

        self._logger.info("Wrote file to: {0}".format(save_fn))

    def write_ubc_files(self, basename, c_east=0, c_north=0, c_z=0):
        """
        Write a UBC .msh and .mod file

        :param save_fn: DESCRIPTION
        :type save_fn: TYPE
        :param c_east: DESCRIPTION, defaults to 0
        :type c_east: TYPE, optional
        :param c_north: DESCRIPTION, defaults to 0
        :type c_north: TYPE, optional
        :param c_z: DESCRIPTION, defaults to 0
        :type c_z: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE


        .. note:: not complete yet.
        """

        # write mesh first
        lines = [
            f"{self.nodes_east.size} {self.nodes_north.size} {self.nodes_z.size}"
        ]
        lines.append(
            str(self.nodes_east.tolist())
            .replace("[", "")
            .replace("]", "")
            .replace(",", "")
        )
        lines.append(
            str(self.nodes_north.tolist())
            .replace("[", "")
            .replace("]", "")
            .replace(",", "")
        )
        lines.append(
            str(self.nodes_z.tolist())
            .replace("[", "")
            .replace("]", "")
            .replace(",", "")
        )

        with open(self.save_path.joinpath(basename + ".msh"), "w") as fid:
            fid.write("\n".join(lines))
