# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:02:57 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import scipy.interpolate as spi

# =============================================================================
class Mesh:
    """
    deals only with the finite element mesh.  Builds a finite element mesh
    based on given parameters defined below.  The mesh reads in the station
    locations, finds the center and makes the relative location of the
    furthest left hand station 0.  The mesh increases in depth logarithmically
    as required by the physics of MT.  Also, the model extends horizontally
    and vertically with padding cells in order to fullfill the assumption of
    the forward operator that at the edges the structure is 1D.  Stations are
    place on the horizontal nodes as required by Wannamaker's forward
    operator.

    Mesh has the ability to create a mesh that incorporates topography given
    a elevation profile.  It adds more cells to the mesh with thickness
    z1_layer.  It then sets the values of the triangular elements according to
    the elevation value at that location.  If the elevation covers less than
    50% of the triangular cell, then the cell value is set to that of air

    .. note:: Mesh is inhereted by Regularization, so the mesh can also be
              be built from there, same as the example below.

    Arguments:
    -----------

    ======================= ===================================================
    Key Words/Attributes    Description
    ======================= ===================================================
    air_key                 letter associated with the value of air
                            *default* is 0
    air_value               value given to an air cell, *default* is 1E13
    cell_width              width of cells with in station area in meters
                            *default* is 100
    elevation_profile       elevation profile along the profile line.
                            given as np.ndarray(nx, 2), where the elements
                            are x_location, elevation.  If elevation profile
                            is given add_elevation is called automatically.
                            *default* is None
    mesh_fn                 full path to mesh file.
    mesh_values             letter values of each triangular mesh element
                            if the cell is free value is ?
    n_layers                number of vertical layers in mesh
                            *default* is 90
    num_x_pad_cells         number of horizontal padding cells outside the
                            the station area that will increase in size
                            by x_pad_multiplier. *default* is 7
    num_x_pad_small_cells   number of horizonal padding cells just outside
                            the station area with width cell_width.  This is
                            to extend the station area if needed.
                            *default* is 2
    num_z_pad_cells         number of vertical padding cells below
                            z_target_depth down to z_bottom. *default* is 5
    rel_station_locations   relative station locations within the mesh.  The
                            locations are relative to the center of the station
                            area.  *default* is None, filled later
    save_path               full path to save mesh file to.
                            *default* is current working directory.
    station_locations       location of stations in meters, can be on a
                            relative grid or in UTM.
    x_grid                  location of horizontal grid nodes in meters
    x_nodes                 relative spacing between grid nodes
    x_pad_multiplier        horizontal padding cells will increase by this
                            multiple out to the edge of the grid.
                            *default* is 1.5
    z1_layer                thickness of the first layer in the model.
                            Should be at least 1/4 of the first skin depth
                            *default* is 10
    z_bottom                bottom depth of the model (m).  Needs to be large
                            enough to be 1D at the edge.
                            *default* is 200000.0
    z_grid                  location of vertical nodes in meters
    z_nodes                 relative distance between vertical nodes in meters
    z_target_depth          depth to deepest target of interest.  Below this
                            depth cells will be padded to z_bottom
    ======================= ===================================================

    ======================= ===================================================
    Methods                 Description
    ======================= ===================================================
    add_elevation           adds elevation to the mesh given elevation
                            profile.
    build_mesh              builds the mesh given the attributes of Mesh.  If
                            elevation_profile is not None, add_elevation is
                            called inside build_mesh
    plot_mesh               plots the built mesh with station location.
    read_mesh_file          reads in an existing mesh file and populates the
                            appropriate attributes.
    write_mesh_file         writes a mesh file to save_path
    ======================= ===================================================


    :Example: ::

        >>> import mtpy.modeling.occam2d as occcam2d
        >>> edipath = r"/home/mt/edi_files"
        >>> slist = ['mt{0:03}'.format(ss) for ss in range(20)]
        >>> ocd = occam2d.Data(edi_path=edipath, station_list=slist)
        >>> ocd.save_path = r"/home/occam/Line1/Inv1"
        >>> ocd.write_data_file()
        >>> ocm = occam2d.Mesh(ocd.station_locations)
        >>> # add in elevation
        >>> ocm.elevation_profile = ocd.elevation_profile
        >>> # change number of layers
        >>> ocm.n_layers = 110
        >>> # change cell width in station area
        >>> ocm.cell_width = 200
        >>> ocm.build_mesh()
        >>> ocm.plot_mesh()
        >>> ocm.save_path = ocd.save_path
        >>> ocm.write_mesh_file()

    """

    def __init__(self, station_locations=None, **kwargs):

        self.station_locations = station_locations
        self.rel_station_locations = None
        self.n_layers = 90
        self.cell_width = 100
        self.num_x_pad_cells = 7
        self.num_z_pad_cells = 5
        self.x_pad_multiplier = 1.5
        self.z1_layer = 10.0
        self.z_bottom = 200000.0
        self.z_target_depth = 50000.0
        self.num_x_pad_small_cells = 2
        self.save_path = Path()
        self.mesh_fn = None
        self.elevation_profile = None

        self.x_nodes = None
        self.z_nodes = None
        self.x_grid = None
        self.z_grid = None
        self.mesh_values = None
        self.air_value = 1e13
        self.air_key = "0"

        for key, value in kwargs.items():
            setattr(self, key, value)

    def build_mesh(self):
        """
        Build the finite element mesh given the parameters defined by the
        attributes of Mesh.  Computes relative station locations by finding
        the center of the station area and setting the middle to 0.  Mesh
        blocks are built by calculating the distance between stations and
        putting evenly spaced blocks between the stations being close to
        cell_width.  This places a horizontal node at the station location.
        If the spacing between stations is smaller than
        cell_width, a horizontal node is placed between the stations to be
        sure the model has room to change between the station.

        If elevation_profile is given, add_elevation is called to add
        topography into the mesh.

        Populates attributes:
            * mesh_values
            * rel_station_locations
            * x_grid
            * x_nodes
            * z_grid
            * z_nodes

        :Example: ::
            >>> import mtpy.modeling.occam2d as occcam2d
            >>> edipath = r"/home/mt/edi_files"
            >>> slist = ['mt{0:03}'.format(ss) for ss in range(20)]
            >>> ocd = occam2d.Data(edi_path=edipath, station_list=slist)
            >>> ocd.save_path = r"/home/occam/Line1/Inv1"
            >>> ocd.write_data_file()
            >>> ocm = occam2d.Mesh(ocd.station_locations)
            >>> # add in elevation
            >>> ocm.elevation_profile = ocd.elevation_profile
            >>> # change number of layers
            >>> ocm.n_layers = 110
            >>> # change cell width in station area
            >>> ocm.cell_width = 200
            >>> ocm.build_mesh()
        """

        if self.station_locations is None:
            raise ValueError(
                "Need to input station locations to define "
                "a finite element mesh"
            )

        # be sure the station locations are sorted from left to right
        self.station_locations.sort()

        self.rel_station_locations = np.copy(self.station_locations)

        # center the stations around 0 like the mesh will be
        self.rel_station_locations -= self.rel_station_locations.mean()

        # 1) make horizontal nodes at station locations and fill in the cells
        #   around that area with cell width. This will put the station
        #   in the center of the regularization block as prescribed for occam
        # the first cell of the station area will be outside of the furthest
        # right hand station to reduce the effect of a large neighboring cell.
        self.x_grid = np.array(
            [
                self.rel_station_locations[0]
                - self.cell_width * self.x_pad_multiplier
            ]
        )

        for ii, offset in enumerate(self.rel_station_locations[:-1]):
            dx = self.rel_station_locations[ii + 1] - offset
            num_cells = int(np.ceil(dx / self.cell_width))
            # if the spacing between stations is smaller than mesh set cell
            # size to mid point between stations

            if num_cells == 0:
                cell_width = dx / 2.0
                num_cells = 1
            # calculate cell spacing so that they are equal between neighboring
            # stations
            else:
                cell_width = dx / num_cells

            if self.x_grid[-1] != offset:
                self.x_grid = np.append(self.x_grid, offset)
            for dd in range(num_cells):
                new_cell = offset + (dd + 1) * cell_width
                # make sure cells aren't too close together
                try:
                    if (
                        abs(self.rel_station_locations[ii + 1] - new_cell)
                        >= cell_width * 0.9
                    ):
                        self.x_grid = np.append(self.x_grid, new_cell)
                    else:
                        pass
                except IndexError:
                    pass

        self.x_grid = np.append(self.x_grid, self.rel_station_locations[-1])
        # add a cell on the right hand side of the station area to reduce
        # effect of a large cell next to it
        self.x_grid = np.append(
            self.x_grid,
            self.rel_station_locations[-1]
            + self.cell_width * self.x_pad_multiplier,
        )

        # --> pad the mesh with exponentially increasing horizontal cells
        #    such that the edge of the mesh can be estimated with a 1D model

        x_left = float(abs(self.x_grid[0] - self.x_grid[1]))
        x_right = float(abs(self.x_grid[-1] - self.x_grid[-2]))
        x_pad_cell = np.max([x_left, x_right])

        for ii in range(self.num_x_pad_cells):
            left_cell = self.x_grid[0]
            right_cell = self.x_grid[-1]
            pad_cell = x_pad_cell * self.x_pad_multiplier ** (ii + 1)
            self.x_grid = np.insert(self.x_grid, 0, left_cell - pad_cell)
            self.x_grid = np.append(self.x_grid, right_cell + pad_cell)

        # --> compute relative positions for the grid
        self.x_nodes = self.x_grid.copy()
        for ii, xx in enumerate(self.x_grid[:-1]):
            self.x_nodes[ii] = abs(self.x_grid[ii + 1] - xx)
        self.x_nodes = self.x_nodes[:-1]

        # 2) make vertical nodes so that they increase with depth
        # --> make depth grid
        log_z = np.logspace(
            np.log10(self.z1_layer),
            np.log10(
                self.z_target_depth
                - np.logspace(
                    np.log10(self.z1_layer),
                    np.log10(self.z_target_depth),
                    num=self.n_layers,
                )[-2]
            ),
            num=self.n_layers - self.num_z_pad_cells,
        )

        # round the layers to be whole numbers
        ztarget = np.array(
            [zz - zz % 10 ** np.floor(np.log10(zz)) for zz in log_z]
        )

        # --> create padding cells past target depth
        log_zpad = np.logspace(
            np.log10(self.z_target_depth),
            np.log10(
                self.z_bottom
                - np.logspace(
                    np.log10(self.z_target_depth),
                    np.log10(self.z_bottom),
                    num=self.num_z_pad_cells,
                )[-2]
            ),
            num=self.num_z_pad_cells,
        )
        # round the layers to be whole numbers
        zpadding = np.array(
            [zz - zz % 10 ** np.floor(np.log10(zz)) for zz in log_zpad]
        )
        zpadding.sort()

        # create the vertical nodes
        self.z_nodes = np.append(ztarget, zpadding)

        # calculate actual distances of depth layers
        self.z_grid = np.array(
            [
                self.z_nodes[: ii + 1].sum()
                for ii in range(self.z_nodes.shape[0])
            ]
        )

        self.mesh_values = np.zeros(
            (self.x_nodes.shape[0], self.z_nodes.shape[0], 4), dtype=str
        )
        self.mesh_values[:, :, :] = "?"

        # get elevation if elevation_profile is given
        if self.elevation_profile is not None:
            self.add_elevation(self.elevation_profile)

        print("=" * 55)
        print("{0:^55}".format("mesh parameters".upper()))
        print("=" * 55)
        print(
            "  number of horizontal nodes = {0}".format(self.x_nodes.shape[0])
        )
        print(
            "  number of vertical nodes   = {0}".format(self.z_nodes.shape[0])
        )
        print(
            "  Total Horizontal Distance  = {0:2f}".format(self.x_nodes.sum())
        )
        print(
            "  Total Vertical Distance    = {0:2f}".format(self.z_nodes.sum())
        )
        print("=" * 55)

    def add_elevation(self, elevation_profile=None):
        """
        the elevation model needs to be in relative coordinates and be a
        numpy.ndarray(2, num_elevation_points) where the first column is
        the horizontal location and the second column is the elevation at
        that location.

        If you have a elevation model use Profile to project the elevation
        information onto the profile line

        To build the elevation I'm going to add the elevation to the top
        of the model which will add cells to the mesh. there might be a better
        way to do this, but this is the first attempt. So I'm going to assume
        that the first layer of the mesh without elevation is the minimum
        elevation and blocks will be added to max elevation at an increment
        according to z1_layer

        .. note:: the elevation model should be symmetrical ie, starting
                  at the first station and ending on the last station, so for
                  now any elevation outside the station area will be ignored
                  and set to the elevation of the station at the extremities.
                  This is not ideal but works for now.

        Arguments:
        -----------
            **elevation_profile** : np.ndarray(2, num_elev_points)
                                    - 1st row is for profile location
                                    - 2nd row is for elevation values

        Computes:
        ---------
            **mesh_values** : mesh values, setting anything above topography
                              to the key for air, which for Occam is '0'
        """
        if elevation_profile is not None:
            self.elevation_profile = elevation_profile

        if self.elevation_profile is None:
            raise ValueError(
                "Need to input an elevation profile to "
                "add elevation into the mesh."
            )

        elev_diff = abs(
            elevation_profile[1].max() - elevation_profile[1].min()
        )
        num_elev_layers = int(elev_diff / self.z1_layer)

        # add vertical nodes and values to mesh_values
        self.z_nodes = np.append(
            [self.z1_layer] * num_elev_layers, self.z_nodes
        )
        self.z_grid = np.array(
            [
                self.z_nodes[: ii + 1].sum()
                for ii in range(self.z_nodes.shape[0])
            ]
        )
        # this assumes that mesh_values have not been changed yet and are all ?
        self.mesh_values = np.zeros(
            (self.x_grid.shape[0], self.z_grid.shape[0], 4), dtype=str
        )
        self.mesh_values[:, :, :] = "?"

        # --> need to interpolate the elevation values onto the mesh nodes
        # first center the locations about 0, this needs to be the same
        # center as the station locations.
        offset = elevation_profile[0] - elevation_profile[0].mean()
        elev = elevation_profile[1] - elevation_profile[1].min()

        func_elev = spi.interp1d(offset, elev, kind="linear")

        # need to figure out which cells and triangular cells need to be air
        xpad = self.num_x_pad_cells + 1
        for ii, xg in enumerate(self.x_grid[xpad:-xpad], xpad):
            # get the index value for z_grid by calculating the elevation
            # difference relative to the top of the model
            dz = elev.max() - func_elev(xg)
            # index of ground in the model for that x location
            zz = int(np.ceil(dz / self.z1_layer))
            if zz == 0:
                pass
            else:
                # --> need to figure out the triangular elements
                # top triangle
                zlayer = elev.max() - self.z_grid[zz]
                try:
                    xtop = xg + (self.x_grid[ii + 1] - xg) / 2
                    ytop = (
                        zlayer
                        + 3 * (self.z_grid[zz] - self.z_grid[zz - 1]) / 4
                    )
                    elev_top = func_elev(xtop)
                    # print xg, xtop, ytop, elev_top, zz
                    if elev_top > ytop:
                        self.mesh_values[ii, 0:zz, 0] = self.air_key
                    else:
                        self.mesh_values[ii, 0 : zz - 1, 0] = self.air_key
                except ValueError:
                    pass

                # left triangle
                try:
                    xleft = xg + (self.x_grid[ii + 1] - xg) / 4.0
                    yleft = (
                        zlayer + (self.z_grid[zz] - self.z_grid[zz - 1]) / 2.0
                    )
                    elev_left = func_elev(xleft)
                    # print xg, xleft, yleft, elev_left, zz
                    if elev_left > yleft:
                        self.mesh_values[ii, 0:zz, 1] = self.air_key
                except ValueError:
                    pass

                # bottom triangle
                try:
                    xbottom = xg + (self.x_grid[ii + 1] - xg) / 2
                    ybottom = (
                        zlayer + (self.z_grid[zz] - self.z_grid[zz - 1]) / 4
                    )
                    elev_bottom = func_elev(xbottom)
                    # print xg, xbottom, ybottom, elev_bottom, zz
                    if elev_bottom > ybottom:
                        self.mesh_values[ii, 0:zz, 2] = self.air_key
                except ValueError:
                    pass

                # right triangle
                try:
                    xright = xg + 3 * (self.x_grid[ii + 1] - xg) / 4
                    yright = (
                        zlayer + (self.z_grid[zz] - self.z_grid[zz - 1]) / 2
                    )
                    elev_right = func_elev(xright)
                    if elev_right > yright * 0.95:
                        self.mesh_values[ii, 0:zz, 3] = self.air_key
                except ValueError:
                    pass
        # --> need to fill out the padding cells so they have the same elevation
        #    as the extremity stations.
        for ii in range(xpad):
            self.mesh_values[ii, :, :] = self.mesh_values[xpad + 1, :, :]
        for ii in range(xpad + 1):
            self.mesh_values[-(ii + 1), :, :] = self.mesh_values[
                -xpad - 2, :, :
            ]

        print("{0:^55}".format("--- Added Elevation to Mesh --"))

    def plot_mesh(self, **kwargs):
        """
        Plot built mesh with station locations.

        =================== ===================================================
        Key Words           Description
        =================== ===================================================
        depth_scale         [ 'km' | 'm' ] scale of mesh plot.
                            *default* is 'km'
        fig_dpi             dots-per-inch resolution of the figure
                            *default* is 300
        fig_num             number of the figure instance
                            *default* is 'Mesh'
        fig_size            size of figure in inches (width, height)
                            *default* is [5, 5]
        fs                  size of font of axis tick labels, axis labels are
                            fs+2. *default* is 6
        ls                  [ '-' | '.' | ':' ] line style of mesh lines
                            *default* is '-'
        marker              marker of stations
                            *default* is r"$\blacktriangledown$"
        ms                  size of marker in points. *default* is 5
        plot_triangles      [ 'y' | 'n' ] to plot mesh triangles.
                            *default* is 'n'
        =================== ===================================================

        """
        fig_num = kwargs.pop("fig_num", "Mesh")
        fig_size = kwargs.pop("fig_size", [5, 5])
        fig_dpi = kwargs.pop("fig_dpi", 300)
        marker = kwargs.pop("marker", r"$\blacktriangledown$")
        ms = kwargs.pop("ms", 5)
        mc = kwargs.pop("mc", "k")
        lw = kwargs.pop("ls", 0.35)
        fs = kwargs.pop("fs", 6)
        plot_triangles = kwargs.pop("plot_triangles", "n")

        depth_scale = kwargs.pop("depth_scale", "km")

        # set the scale of the plot
        if depth_scale == "km":
            df = 1000.0
        elif depth_scale == "m":
            df = 1.0
        else:
            df = 1000.0

        plt.rcParams["figure.subplot.left"] = 0.12
        plt.rcParams["figure.subplot.right"] = 0.98
        plt.rcParams["font.size"] = fs

        if self.x_grid is None:
            self.build_mesh()

        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        ax = fig.add_subplot(1, 1, 1, aspect="equal")

        # plot the station marker
        # plots a V for the station cause when you use scatter the spacing
        # is variable if you change the limits of the y axis, this way it
        # always plots at the surface.
        for offset in self.rel_station_locations:
            ax.text(
                (offset) / df,
                0,
                marker,
                horizontalalignment="center",
                verticalalignment="baseline",
                fontdict={"size": ms, "color": mc},
            )

        # --> make list of column lines
        row_line_xlist = []
        row_line_ylist = []
        for xx in self.x_grid / df:
            row_line_xlist.extend([xx, xx])
            row_line_xlist.append(None)
            row_line_ylist.extend([0, self.z_grid[-1] / df])
            row_line_ylist.append(None)

        # plot column lines (variables are a little bit of a misnomer)
        ax.plot(row_line_xlist, row_line_ylist, color="k", lw=lw)

        # --> make list of row lines
        col_line_xlist = [self.x_grid[0] / df, self.x_grid[-1] / df]
        col_line_ylist = [0, 0]
        for yy in self.z_grid / df:
            col_line_xlist.extend([self.x_grid[0] / df, self.x_grid[-1] / df])
            col_line_xlist.append(None)
            col_line_ylist.extend([yy, yy])
            col_line_ylist.append(None)

        # plot row lines (variables are a little bit of a misnomer)
        ax.plot(col_line_xlist, col_line_ylist, color="k", lw=lw)

        if plot_triangles == "y":
            row_line_xlist = []
            row_line_ylist = []
            for xx in self.x_grid / df:
                row_line_xlist.extend([xx, xx])
                row_line_xlist.append(None)
                row_line_ylist.extend([0, self.z_grid[-1] / df])
                row_line_ylist.append(None)

            # plot columns
            ax.plot(row_line_xlist, row_line_ylist, color="k", lw=lw)

            col_line_xlist = []
            col_line_ylist = []
            for yy in self.z_grid / df:
                col_line_xlist.extend(
                    [self.x_grid[0] / df, self.x_grid[-1] / df]
                )
                col_line_xlist.append(None)
                col_line_ylist.extend([yy, yy])
                col_line_ylist.append(None)

            # plot rows
            ax.plot(col_line_xlist, col_line_ylist, color="k", lw=lw)

            diag_line_xlist = []
            diag_line_ylist = []
            for xi, xx in enumerate(self.x_grid[:-1] / df):
                for yi, yy in enumerate(self.z_grid[:-1] / df):
                    diag_line_xlist.extend([xx, self.x_grid[xi + 1] / df])
                    diag_line_xlist.append(None)
                    diag_line_xlist.extend([xx, self.x_grid[xi + 1] / df])
                    diag_line_xlist.append(None)

                    diag_line_ylist.extend([yy, self.z_grid[yi + 1] / df])
                    diag_line_ylist.append(None)
                    diag_line_ylist.extend([self.z_grid[yi + 1] / df, yy])
                    diag_line_ylist.append(None)

            # plot diagonal lines.
            ax.plot(diag_line_xlist, diag_line_ylist, color="k", lw=lw)

        # --> set axes properties
        ax.set_ylim(self.z_target_depth / df, -2000 / df)
        xpad = self.num_x_pad_cells - 1
        ax.set_xlim(self.x_grid[xpad] / df, -self.x_grid[xpad] / df)
        ax.set_xlabel(
            "Easting ({0})".format(depth_scale),
            fontdict={"size": fs + 2, "weight": "bold"},
        )
        ax.set_ylabel(
            "Depth ({0})".format(depth_scale),
            fontdict={"size": fs + 2, "weight": "bold"},
        )
        plt.show()

    def write_mesh_file(self, save_path=None, basename="Occam2DMesh"):
        """
        Write a finite element mesh file.

        Calls build_mesh if it already has not been called.

        Arguments:
        -----------
            **save_path** : string
                            directory path or full path to save file

            **basename** : string
                           basename of mesh file. *default* is 'Occam2DMesh'
        Returns:
        ----------
            **mesh_fn** : string
                          full path to mesh file

        :example: ::

            >>> import mtpy.modeling.occam2d as occam2d
            >>> edi_path = r"/home/mt/edi_files"
            >>> profile = occam2d.Profile(edi_path)
            >>> profile.plot_profile()
            >>> mesh = occam2d.Mesh(profile.station_locations)
            >>> mesh.build_mesh()
            >>> mesh.write_mesh_file(save_path=r"/home/occam2d/Inv1")
        """

        if save_path is not None:
            self.save_path = Path(save_path)

        if self.save_path is None:
            self.save_path = Path()

        self.mesh_fn = self.save_path.joinpath(basename)

        if self.x_nodes is None:
            self.build_mesh()

        mesh_lines = []
        nx = self.x_nodes.shape[0]
        nz = self.z_nodes.shape[0]
        mesh_lines.append("MESH FILE Created by mtpy.modeling.occam2d\n")
        mesh_lines.append(
            "   {0}  {1}  {2}  {0}  {0}  {3}\n".format(0, nx + 1, nz + 1, 2)
        )

        # --> write horizontal nodes
        node_str = ""
        for ii, xnode in enumerate(self.x_nodes):
            node_str += "{0:>9.1f} ".format(xnode)
            if np.remainder(ii + 1, 8) == 0:
                node_str += "\n"
                mesh_lines.append(node_str)
                node_str = ""

        node_str += "\n"
        mesh_lines.append(node_str)

        # --> write vertical nodes
        node_str = ""
        for ii, znode in enumerate(self.z_nodes):
            node_str += "{0:>9.1f} ".format(znode)
            if np.remainder(ii + 1, 8) == 0:
                node_str += "\n"
                mesh_lines.append(node_str)
                node_str = ""
        node_str += "\n"
        mesh_lines.append(node_str)

        # --> need a 0 after the nodes
        mesh_lines.append("    0\n")
        # --> write triangular mesh block values as ?
        for zz in range(self.z_nodes.shape[0]):
            for tt in range(4):
                mesh_lines.append("".join(self.mesh_values[:, zz, tt]) + "\n")

        with open(self.mesh_fn, "w") as fid:
            mfid.writelines(mesh_lines)

        print("Wrote Mesh file to {0}".format(self.mesh_fn))

    def read_mesh_file(self, mesh_fn):
        """
        reads an occam2d 2D mesh file

        Arguments:
        ----------
            **mesh_fn** : string
                          full path to mesh file

        Populates:
        -----------
            **x_grid** : array of horizontal locations of nodes (m)

            **x_nodes**: array of horizontal node relative distances
                        (column locations (m))

            **z_grid** : array of vertical node locations (m)

            **z_nodes** : array of vertical nodes
                          (row locations(m))

            **mesh_values** : np.array of free parameters

        To do:
        ------
            incorporate fixed values

        :Example: ::

            >>> import mtpy.modeling.occam2d as occam2d
            >>> mg = occam2d.Mesh()
            >>> mg.mesh_fn = r"/home/mt/occam/line1/Occam2Dmesh"
            >>> mg.read_mesh_file()
        """
        self.mesh_fn = mesh_fn

        with open(self.mesh_fn, "r") as mfid:
            mlines = mfid.readlines()

        nh = int(mlines[1].strip().split()[1])
        nv = int(mlines[1].strip().split()[2])

        self.x_nodes = np.zeros(nh)
        self.z_nodes = np.zeros(nv)
        self.mesh_values = np.zeros((nh, nv, 4), dtype=str)

        # get horizontal nodes
        h_index = 0
        v_index = 0
        m_index = 0
        line_count = 2

        # --> fill horizontal nodes
        for mline in mlines[line_count:]:
            mline = mline.strip().split()
            for m_value in mline:
                # print m_value, h_index
                self.x_nodes[h_index] = float(m_value)
                h_index += 1
                if h_index == nh - 1:
                    break
            line_count += 1
            if h_index == nh - 1:
                break

                # --> fill vertical nodes
        for mline in mlines[line_count:]:
            mline = mline.strip().split()
            for m_value in mline:
                self.z_nodes[v_index] = float(m_value)
                v_index += 1
                if v_index == nv - 1:
                    break
            line_count += 1
            if v_index == nv - 1:
                break

            line_count += 1

        # --> fill model values
        for ll, mline in enumerate(mlines[line_count + 1 :], line_count):
            mline = mline.strip()
            if m_index == nv or mline.lower().find("exception") > 0:
                break
            else:
                mlist = list(mline)
                if len(mlist) != nh - 1:
                    print("--- Line {0} in {1}".format(ll, self.mesh_fn))
                    print("Check mesh file too many columns")
                    print("Should be {0}, has {1}".format(nh, len(mlist)))
                    mlist = mlist[0:nh]
                for kk in range(4):
                    for jj, mvalue in enumerate(list(mlist)):
                        self.mesh_values[jj, m_index, kk] = mline[jj]
                m_index += 1

        # sometimes it seems that the number of nodes is not the same as the
        # header would suggest so need to remove the zeros
        self.x_nodes = self.x_nodes[np.nonzero(self.x_nodes)]
        if self.x_nodes.shape[0] != nh - 1:
            new_nh = self.x_nodes.shape[0]
            print(
                "The header number {0} should read {1}".format(nh - 1, new_nh)
            )
            self.mesh_values.resize(new_nh, nv, 4)
        else:
            new_nh = nh

        self.z_nodes = self.z_nodes[np.nonzero(self.z_nodes)]
        if self.z_nodes.shape[0] != nv - 1:
            new_nv = self.z_nodes.shape[0]
            print(f"The header number {nv-1} should read {new_nv}")
            self.mesh_values.resize(new_nh, nv, 4)

        # make x_grid and z_grid
        self.x_grid = self.x_nodes.copy()
        self.x_grid = np.append(self.x_grid, self.x_grid[-1])
        self.x_grid = np.array(
            [self.x_grid[:ii].sum() for ii in range(self.x_grid.shape[0])]
        )
        self.x_grid -= self.x_grid.mean()
        self.z_grid = np.array(
            [self.z_nodes[:ii].sum() for ii in range(self.z_nodes.shape[0])]
        )
