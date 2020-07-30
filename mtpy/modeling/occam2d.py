# -*- coding: utf-8 -*-
"""
Spin-off from 'occamtools'
(Created August 2011, re-written August 2013)

Tools for Occam2D

authors: JP/LK


Classes:
    - Data
    - Model
    - Setup
    - Run
    - Plot
    - Mask


Functions:
    - getdatetime
    - makestartfiles
    - writemeshfile
    - writemodelfile
    - writestartupfile
    - read_datafile
    - get_model_setup
    - blocks_elements_setup


"""
# ==============================================================================
import numpy as np
import scipy as sp
from scipy.stats import mode
import os
import os.path as op
import time
import matplotlib.colorbar as mcb
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.interpolate as spi
import mtpy.core.mt as mt
import mtpy.modeling.winglink as MTwl
import mtpy.analysis.geometry as MTgy
from mtpy.imaging.mtplottools import plot_errorbar
import mtpy.utils.calculator as mtcc
import mtpy.utils.mesh_tools as mtmesh
from mtpy.utils import gis_tools


# ==============================================================================

class Mesh():
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
        self.n_layers = kwargs.pop('n_layers', 90)
        self.cell_width = kwargs.pop('cell_width', 100)
        self.num_x_pad_cells = kwargs.pop('num_x_pad_cells', 7)
        self.num_z_pad_cells = kwargs.pop('num_z_pad_cells', 5)
        self.x_pad_multiplier = kwargs.pop('x_pad_multiplier', 1.5)
        self.z1_layer = kwargs.pop('z1_layer', 10.0)
        self.z_bottom = kwargs.pop('z_bottom', 200000.0)
        self.z_target_depth = kwargs.pop('z_target_depth', 50000.0)
        self.num_x_pad_small_cells = kwargs.pop('num_x_pad_small_cells', 2)
        self.save_path = kwargs.pop('save_path', None)
        self.mesh_fn = kwargs.pop('mesh_fn', None)
        self.elevation_profile = kwargs.pop('elevation_profile', None)

        self.x_nodes = None
        self.z_nodes = None
        self.x_grid = None
        self.z_grid = None
        self.mesh_values = None
        self.air_value = 1e13
        self.air_key = '0'

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
            raise OccamInputError('Need to input station locations to define '
                                  'a finite element mesh')

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
        self.x_grid = np.array([self.rel_station_locations[0] - self.cell_width * \
                                self.x_pad_multiplier])
        for ii, offset in enumerate(self.rel_station_locations[:-1]):
            dx = self.rel_station_locations[ii + 1] - offset
            num_cells = int(np.floor(dx / self.cell_width))
            # if the spacing between stations is smaller than mesh set cell
            # size to mid point between stations
            if num_cells == 0:
                cell_width = dx / 2.
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
                    if abs(self.rel_station_locations[ii + 1] - new_cell) >= cell_width * .9:
                        self.x_grid = np.append(self.x_grid, new_cell)
                    else:
                        pass
                except IndexError:
                    pass

        self.x_grid = np.append(self.x_grid, self.rel_station_locations[-1])
        # add a cell on the right hand side of the station area to reduce 
        # effect of a large cell next to it       
        self.x_grid = np.append(self.x_grid,
                                self.rel_station_locations[-1] + self.cell_width * \
                                self.x_pad_multiplier)
        
        # add an extra cell if there is an uneven number of cells
        if len(self.x_grid) % 2 == 0:
            self.x_grid = np.append(self.x_grid,
                                self.x_grid[-1] + self.cell_width * \
                                self.x_pad_multiplier)                                
        
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
        #        log_z = np.logspace(np.log10(self.z1_layer),
        #                            np.log10(self.z_target_depth-\
        #                                     np.logspace(np.log10(self.z1_layer),
        #                            np.log10(self.z_target_depth),
        #                            num=self.n_layers)[-2]),
        #                            num=self.n_layers-self.num_z_pad_cells)

        log_z = mtmesh.make_log_increasing_array(self.z1_layer,
                                                 self.z_target_depth,
                                                 self.n_layers - self.num_z_pad_cells,
                                                 increment_factor=0.99)

        # round the layers to be whole numbers
        ztarget = np.array([mtcc.roundsf(zz, 1) for zz in
                            log_z])  # zz-zz%10**np.floor(np.log10(zz))

        # --> create padding cells past target depth
        #        log_zpad = np.logspace(np.log10(self.z_target_depth),
        #                            np.log10(self.z_bottom-\
        #                                    np.logspace(np.log10(self.z_target_depth),
        #                            np.log10(self.z_bottom),
        #                            num=self.num_z_pad_cells)[-2]),
        #                            num=self.num_z_pad_cells)
        log_zpad = mtmesh.make_log_increasing_array(log_z[-1],
                                                    self.z_bottom,
                                                    self.num_z_pad_cells,
                                                    increment_factor=0.99)
        # round the layers to 1 sf
        zpadding = np.array([mtcc.roundsf(zz, 1) for zz in
                             log_zpad])  # zz-zz%10**np.floor(np.log10(zz))

        # create the vertical nodes
        self.z_nodes = np.append(ztarget, zpadding)

        # calculate actual distances of depth layers
        self.z_grid = np.array([self.z_nodes[:ii + 1].sum()
                                for ii in range(self.z_nodes.shape[0])])

        self.mesh_values = np.zeros((self.x_nodes.shape[0],
                                     self.z_nodes.shape[0], 4), dtype=str)
        self.mesh_values[:, :, :] = '?'

        # get elevation if elevation_profile is given
        if self.elevation_profile is not None:
            self.add_elevation(self.elevation_profile)

        print('=' * 55)
        print('{0:^55}'.format('mesh parameters'.upper()))
        print('=' * 55)
        print('  number of horizontal nodes = {0}'.format(self.x_nodes.shape[0]))
        print('  number of vertical nodes   = {0}'.format(self.z_nodes.shape[0]))
        print('  Total Horizontal Distance  = {0:2f}'.format(self.x_nodes.sum()))
        print('  Total Vertical Distance    = {0:2f}'.format(self.z_nodes.sum()))
        print('=' * 55)

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
            raise OccamInputError('Need to input an elevation profile to '
                                  'add elevation into the mesh.')

        elev_diff = abs(elevation_profile[1].max() - elevation_profile[1].min())
        num_elev_layers = int(elev_diff / self.z1_layer)

        # add vertical nodes and values to mesh_values
        self.z_nodes = np.append([self.z1_layer] * num_elev_layers, self.z_nodes)
        self.z_grid = np.array([self.z_nodes[:ii + 1].sum()
                                for ii in range(self.z_nodes.shape[0])])
        # this assumes that mesh_values have not been changed yet and are all ?
        self.mesh_values = np.zeros((self.x_grid.shape[0],
                                     self.z_grid.shape[0], 4), dtype=str)
        self.mesh_values[:, :, :] = '?'

        # --> need to interpolate the elevation values onto the mesh nodes
        # first center the locations about 0, this needs to be the same
        # center as the station locations.
        offset = elevation_profile[0] - elevation_profile[0].mean()
        elev = elevation_profile[1] - elevation_profile[1].min()

        func_elev = spi.interp1d(offset, elev, kind='linear')

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
                    ytop = zlayer + 3 * (self.z_grid[zz] - self.z_grid[zz - 1]) / 4
                    elev_top = func_elev(xtop)
                    # print xg, xtop, ytop, elev_top, zz
                    if elev_top > ytop:
                        self.mesh_values[ii, 0:zz, 0] = self.air_key
                    else:
                        self.mesh_values[ii, 0:zz - 1, 0] = self.air_key
                except ValueError:
                    pass

                # left triangle
                try:
                    xleft = xg + (self.x_grid[ii + 1] - xg) / 4.
                    yleft = zlayer + (self.z_grid[zz] - self.z_grid[zz - 1]) / 2.
                    elev_left = func_elev(xleft)
                    # print xg, xleft, yleft, elev_left, zz
                    if elev_left > yleft:
                        self.mesh_values[ii, 0:zz, 1] = self.air_key
                except ValueError:
                    pass

                # bottom triangle
                try:
                    xbottom = xg + (self.x_grid[ii + 1] - xg) / 2
                    ybottom = zlayer + (self.z_grid[zz] - self.z_grid[zz - 1]) / 4
                    elev_bottom = func_elev(xbottom)
                    # print xg, xbottom, ybottom, elev_bottom, zz
                    if elev_bottom > ybottom:
                        self.mesh_values[ii, 0:zz, 2] = self.air_key
                except ValueError:
                    pass

                # right triangle
                try:
                    xright = xg + 3 * (self.x_grid[ii + 1] - xg) / 4
                    yright = zlayer + (self.z_grid[zz] - self.z_grid[zz - 1]) / 2
                    elev_right = func_elev(xright)
                    if elev_right > yright * .95:
                        self.mesh_values[ii, 0:zz, 3] = self.air_key
                except ValueError:
                    pass
        # --> need to fill out the padding cells so they have the same elevation
        #    as the extremity stations.
        for ii in range(xpad):
            self.mesh_values[ii, :, :] = self.mesh_values[xpad + 1, :, :]
        for ii in range(xpad + 1):
            self.mesh_values[-(ii + 1), :, :] = self.mesh_values[-xpad - 2, :, :]

        print('{0:^55}'.format('--- Added Elevation to Mesh --'))

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
        fig_num = kwargs.pop('fig_num', 'Mesh')
        fig_size = kwargs.pop('fig_size', [5, 5])
        fig_dpi = kwargs.pop('fig_dpi', 300)
        marker = kwargs.pop('marker', r"$\blacktriangledown$")
        ms = kwargs.pop('ms', 5)
        mc = kwargs.pop('mc', 'k')
        lw = kwargs.pop('ls', .35)
        fs = kwargs.pop('fs', 6)
        plot_triangles = kwargs.pop('plot_triangles', 'n')

        depth_scale = kwargs.pop('depth_scale', 'km')

        # set the scale of the plot
        if depth_scale == 'km':
            df = 1000.
        elif depth_scale == 'm':
            df = 1.
        else:
            df = 1000.

        plt.rcParams['figure.subplot.left'] = .12
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['font.size'] = fs

        if self.x_grid is None:
            self.build_mesh()

        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        ax = fig.add_subplot(1, 1, 1, aspect='equal')

        # plot the station marker
        # plots a V for the station cause when you use scatter the spacing
        # is variable if you change the limits of the y axis, this way it
        # always plots at the surface.
        for offset in self.rel_station_locations:
            ax.text((offset) / df,
                    0,
                    marker,
                    horizontalalignment='center',
                    verticalalignment='baseline',
                    fontdict={'size': ms, 'color': mc})

        # --> make list of column lines
        row_line_xlist = []
        row_line_ylist = []
        for xx in self.x_grid / df:
            row_line_xlist.extend([xx, xx])
            row_line_xlist.append(None)
            row_line_ylist.extend([0, self.z_grid[-1] / df])
            row_line_ylist.append(None)

        # plot column lines (variables are a little bit of a misnomer)
        ax.plot(row_line_xlist,
                row_line_ylist,
                color='k',
                lw=lw)

        # --> make list of row lines
        col_line_xlist = [self.x_grid[0] / df, self.x_grid[-1] / df]
        col_line_ylist = [0, 0]
        for yy in self.z_grid / df:
            col_line_xlist.extend([self.x_grid[0] / df,
                                   self.x_grid[-1] / df])
            col_line_xlist.append(None)
            col_line_ylist.extend([yy, yy])
            col_line_ylist.append(None)

        # plot row lines (variables are a little bit of a misnomer)
        ax.plot(col_line_xlist,
                col_line_ylist,
                color='k',
                lw=lw)

        if plot_triangles == 'y':
            row_line_xlist = []
            row_line_ylist = []
            for xx in self.x_grid / df:
                row_line_xlist.extend([xx, xx])
                row_line_xlist.append(None)
                row_line_ylist.extend([0, self.z_grid[-1] / df])
                row_line_ylist.append(None)

            # plot columns
            ax.plot(row_line_xlist,
                    row_line_ylist,
                    color='k',
                    lw=lw)

            col_line_xlist = []
            col_line_ylist = []
            for yy in self.z_grid / df:
                col_line_xlist.extend([self.x_grid[0] / df,
                                       self.x_grid[-1] / df])
                col_line_xlist.append(None)
                col_line_ylist.extend([yy, yy])
                col_line_ylist.append(None)

            # plot rows
            ax.plot(col_line_xlist,
                    col_line_ylist,
                    color='k',
                    lw=lw)

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
            ax.plot(diag_line_xlist,
                    diag_line_ylist,
                    color='k',
                    lw=lw)

        # --> set axes properties
        ax.set_ylim(self.z_target_depth / df, -2000 / df)
        xpad = self.num_x_pad_cells - 1
        ax.set_xlim(self.x_grid[xpad] / df, -self.x_grid[xpad] / df)
        ax.set_xlabel('Easting ({0})'.format(depth_scale),
                      fontdict={'size': fs + 2, 'weight': 'bold'})
        ax.set_ylabel('Depth ({0})'.format(depth_scale),
                      fontdict={'size': fs + 2, 'weight': 'bold'})
        plt.show()

    def write_mesh_file(self, save_path=None, basename='Occam2DMesh'):
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
            self.save_path = save_path

        if self.save_path is None:
            self.save_path = os.getcwd()

        self.mesh_fn = os.path.join(self.save_path, basename)

        if self.x_nodes is None:
            self.build_mesh()

        mesh_lines = []
        nx = self.x_nodes.shape[0]
        nz = self.z_nodes.shape[0]
        mesh_lines.append('MESH FILE Created by mtpy.modeling.occam2d\n')
        mesh_lines.append("   {0}  {1}  {2}  {0}  {0}  {3}\n".format(0, nx + 1,
                                                                     nz + 1, 2))

        # --> write horizontal nodes
        node_str = ''
        for ii, xnode in enumerate(self.x_nodes):
            node_str += '{0:>9.1f} '.format(xnode)
            if np.remainder(ii + 1, 8) == 0:
                node_str += '\n'
                mesh_lines.append(node_str)
                node_str = ''

        node_str += '\n'
        mesh_lines.append(node_str)

        # --> write vertical nodes
        node_str = ''
        for ii, znode in enumerate(self.z_nodes):
            node_str += '{0:>9.1f} '.format(znode)
            if np.remainder(ii + 1, 8) == 0:
                node_str += '\n'
                mesh_lines.append(node_str)
                node_str = ''
        node_str += '\n'
        mesh_lines.append(node_str)

        # --> need a 0 after the nodes
        mesh_lines.append('    0\n')
        # --> write triangular mesh block values as ?
        for zz in range(self.z_nodes.shape[0]):
            for tt in range(4):
                mesh_lines.append(''.join(self.mesh_values[:, zz, tt]) + '\n')

        with open(self.mesh_fn, 'w') as mfid:
            mfid.writelines(mesh_lines)
        # mfid.close()

        print('Wrote Mesh file to {0}'.format(self.mesh_fn))

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

        with open(self.mesh_fn, 'r') as mfid:
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
#            print mline
            for m_value in mline:
                self.x_nodes[h_index] = float(m_value)
                h_index += 1
            line_count += 1
            # needs to be >= nh - 2 as python count starts from 0 and number of
            # horizontal columns is 1 less than listed at the top of the file
            if h_index >= nh - 2:
                break
            


        # --> fill vertical nodes
        for mline in mlines[line_count:]:
            mline = mline.strip().split()
            for m_value in mline:
                self.z_nodes[v_index] = float(m_value)
                v_index += 1
            line_count += 1
            if v_index >= nv - 2:
                break

        # --> fill model values
        for ll, mline in enumerate(mlines[line_count + 1:], line_count):
            mline = mline.strip()
            if m_index == nv or mline.lower().find('exception') > 0:
                break
            else:
                mlist = list(mline)
                if len(mlist) != nh - 1:
                    print('--- Line {0} in {1}'.format(ll, self.mesh_fn))
                    print('Check mesh file too many columns')
                    print('Should be {0}, has {1}'.format(nh, len(mlist)))
                    mlist = mlist[0:nh]
                for kk in range(4):
                    for jj, mvalue in enumerate(list(mlist)):
                        self.mesh_values[jj, m_index, kk] = mline[jj]
                m_index += 1

        # sometimes it seems that the number of nodes is not the same as the
        # header would suggest so need to remove the zeros
        self.x_nodes = self.x_nodes[np.nonzero(self.x_nodes)]
        if self.x_nodes.shape[0] != nh:
            new_nh = self.x_nodes.shape[0]
            print('The header number {0} should read {1}'.format(nh, new_nh))
            self.mesh_values = np.resize(self.mesh_values, (new_nh, nv, 4))
        else:
            new_nh = nh

        self.z_nodes = self.z_nodes[np.nonzero(self.z_nodes)]
        if self.z_nodes.shape[0] != nv:
            new_nv = self.z_nodes.shape[0]
            print('The header number {0} should read {1}'.format(nv, new_nv))
            self.mesh_values = np.resize(self.mesh_values, (new_nh, nv, 4))

        # make x_grid and z_grid
        self.x_grid = self.x_nodes.copy()
        self.x_grid = np.append(self.x_grid, self.x_grid[-1])
        self.x_grid = np.array([self.x_grid[:ii].sum()
                                for ii in range(self.x_grid.shape[0])])
        self.x_grid -= self.x_grid.mean()
        self.z_grid = np.array([self.z_nodes[:ii].sum()
                                for ii in range(self.z_nodes.shape[0])])


class Profile():
    """
    Takes data from .edi files to create a profile line for 2D modeling.
    Can project the stations onto a profile that is perpendicular to strike
    or a given profile direction.
    
    If _rotate_to_strike is True, the impedance tensor and tipper are rotated
    to align with the geoelectric strike angle.
    
    If _rotate_to_strike is True and geoelectric_strike is not given, 
    then it is calculated using the phase tensor.  First, 2D sections are
    estimated from the impedance tensort hen the strike is esitmated from the
    phase tensor azimuth + skew.  This angle is then used to project the 
    stations perpendicular to the strike angle.
    
    If you want to project onto an angle not perpendicular to strike, give
    profile_angle and set _rotate_to_strike to False.  This will project
    the impedance tensor and tipper to be perpendicular with the 
    profile_angle.

    Arguments:
    -----------
    
        **edi_path** : string
                       full path to edi files
                       
        **station_list** : list of stations to create profile for if None is
                           given all .edi files in edi_path will be used.
                           .. note:: that algorithm assumes .edi files are 
                                     named by station and it only looks for 
                                     the station within the .edi file name
                                     it does not match exactly, so if you have
                                     .edi files with similar names there
                                     might be some problems.
                                     
        **geoelectric_strike** : float
                                 geoelectric strike direction in degrees 
                                 assuming 0 is North and East is 90
                                 
        **profile_angle** : float
                            angle to project the stations onto a profile line
                            .. note:: the geoelectric strike angle and 
                                      profile angle should be orthogonal for
                                      best results from 2D modeling.
                                      
    
    ======================= ===================================================
    **Attributes**          Description    
    ======================= ===================================================
    edi_list                list of mtpy.core.mt.MT instances for each .edi
                            file found in edi_path 
    elevation_model         numpy.ndarray(3, num_elevation_points) elevation
                            values for the profile line (east, north, elev)
    geoelectric_strike      geoelectric strike direction assuming N == 0
    profile_angle           angle of profile line assuming N == 0
    profile_line            (slope, N-intercept) of profile line
    _profile_generated      [ True | False ] True if profile has already been 
                            generated
    edi_path                path to find .edi files
    station_list            list of stations to extract from edi_path
    num_edi                 number of edi files to create a profile for
    _rotate_to_strike       [ True | False] True to project the stations onto
                            a line that is perpendicular to geoelectric strike
                            also Z and Tipper are rotated to strike direction.
    ======================= ===================================================
 
    .. note:: change _rotate_to_strike to False if you want to project the 
              stations onto a given profile direction.  This will rotate
              Z and Tipper to be orthogonal to this direction
   
    ======================= ===================================================
    Methods                 Description
    ======================= ===================================================
    generate_profile        generates a profile for the given stations
    plot_profile            plots the profile line along with original station
                            locations to compare.  
    ======================= ===================================================
    
    :Example: ::
        
        >>> import mtpy.modeling.occam2d as occam
        >>> edi_path = r"/home/mt/edi_files"
        >>> station_list = ['mt{0:03}'.format(ss) for ss in range(0, 15)]
        >>> prof_line = occam.Profile(edi_path, station_list=station_list)
        >>> prof_line.plot_profile()
        >>> #if you want to project to a given strike
        >>> prof_line.geoelectric_strike = 36.7
        >>> prof_line.generate_profile()
        >>> prof_line.plot_profile()

        
    """

    def __init__(self, edi_path=None, **kwargs):

        self.edi_path = edi_path
        self.station_list = kwargs.pop('station_list', None)
        self.geoelectric_strike = kwargs.pop('geoelectric_strike', None)
        self.profile_angle = kwargs.pop('profile_angle', None)
        self.edi_list = []
        self._rotate_to_strike = True
        self.num_edi = 0
        self._profile_generated = False
        self.profile_line = None
        self.station_locations = None
        self.elevation_model = kwargs.pop('elevation_model', None)
        self.elevation_profile = None
        self.estimate_elevation = True
        self.optimize_line = kwargs.pop('optimize_line', False)

    @staticmethod
    def _optimized_path(eastings, northings, start=None):
        """
        TODO: should be refactored to somewhere else, is also used in
        modem.plot_slices.

        This function was adopted verbatim from the following link:

        https://stackoverflow.com/questions/45829155/sort-points-in-order-to-have-a-continuous-curve-using-python

        Order coordinates such that a continuous line can be formed.

        Parameters
        ----------
        eastings : np.ndarray
            1D array of float, easting coordinates for profile line
        northings : np.ndarray
            1D array of float, northing coordinates for profile line

        Returns
        -------
        list
            Index list for sorting E/N coords optimally (and by 
            extension, a list of EDI objects).
        """
        def distance(P1, P2):
            return ((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2) ** 0.5

        coords = np.column_stack((eastings, northings))
        coords_list = coords.tolist()
        # Assume start of profile is minimum easting
        if start is None:
            start = min(coords_list, key=lambda x: x[0])
        pass_by = coords_list
        path = [start]
        pass_by.remove(start)
        while pass_by:
            nearest = min(pass_by, key=lambda x: distance(path[-1], x))
            path.append(nearest)
            pass_by.remove(nearest)
        inds = [np.where(coords == p)[0][0] for p in path]
        return inds

    def _get_edi_list(self):
        """
        Get a list of edi files that coorespond to the station list.
        The ordering of stations is important as it affects the 
        ordering of the easting and northing points that are used
        in calculating the slope and intercept of the profile line.

        Three options for doing so:
            1. Provide self.station_list - the MT objects will be ordered
            by the station names provided.
            
            2. Set self.optimize_line to True. This will sort the stations
            by optimizing for the most continuous linear regression
            through the station locations. Assumes minimum easting
            is start of profile.

            3. Sort the station names alphabetically.
        
        each element of the list is a mtpy.core.mt.MT object
        """

        if self.station_list is not None:
            for station in self.station_list:
                for edi in os.listdir(self.edi_path):
                    if edi.find(station) == 0 and edi[-3:] == 'edi':
                        self.edi_list.append(mt.MT(os.path.join(self.edi_path,
                                                                edi)))
                        break
        elif self.optimize_line:
            edis = [mt.MT(os.path.join(self.edi_path, edi)) for
                    edi in os.listdir(self.edi_path) if edi.endswith('.edi')]
            eastings = np.array([edi.east for edi in edis])
            northings = np.array([edi.north for edi in edis])
            optimal_inds = Profile._optimized_path(eastings, northings)
            self.edi_list = np.asarray(edis)[optimal_inds].tolist()
            self.station_list = [os.path.splitext(os.path.basename(edi.fn))[0]
                                 for edi in self.edi_list]
        else:
            self.edi_list = [mt.MT(os.path.join(self.edi_path, edi)) for
                             edi in os.listdir(self.edi_path) if edi.endswith('.edi')]
            self.station_list = [os.path.splitext(os.path.basename(edi.fn))[0]
                                 for edi in self.edi_list]
        for edi in self.edi_list:
            if type(edi.Tipper.rotation_angle) is list:
                edi.Tipper.rotation_angle = np.array(edi.Tipper.rotation_angle)
        self.num_edi = len(self.edi_list)

    def generate_profile(self):
        """
        Generate linear profile by regression of station locations.

        If profile_angle is not None, then station are projected onto that 
        line.  Else, the a geoelectric strike is calculated from the data
        and the stations are projected onto an angle perpendicular to the
        estimated strike direction.  If _rotate_to_strike is True, the 
        impedance tensor and Tipper data are rotated to align with strike.
        Else, data is not rotated to strike.
        
        To project stations onto a given line, set profile_angle and 
        _rotate_to_strike to False.  This will project the stations onto 
        profile_angle and rotate the impedance tensor and tipper to be 
        perpendicular to the profile_angle.  
        """

        self._get_edi_list()

        strike_angles = np.zeros(self.num_edi)

        easts = np.zeros(self.num_edi)
        norths = np.zeros(self.num_edi)
        utm_zones = np.zeros(self.num_edi)
        
        if self.model_epsg is None:
            latlist = np.array([mtObj.lat for mtObj in self.edi_list])
            lonlist = np.array([mtObj.lon for mtObj in self.edi_list])
            lonc,latc = mtcc.centre_point(lonlist,latlist)
            self.model_epsg = gis_tools.get_epsg(latc,lonc)

        for ii, edi in enumerate(self.edi_list):
            # find strike angles for each station if a strike angle is not given
            if self.geoelectric_strike is None:
                try:
                    # check dimensionality to be sure strike is estimate for 2D
                    dim = MTgy.dimensionality(z_object=edi.Z)
                    # get strike for only those periods
                    gstrike = MTgy.strike_angle(edi.Z.z[np.where(dim == 2)])[:, 0]
                    strike_angles[ii] = np.median(gstrike)
                except:
                    pass

            if self.model_epsg is not None:
                edi.east,edi.north,edi.utm_zone = \
                gis_tools.project_point_ll2utm(edi.lat,
                                               edi.lon,
                                               epsg=self.model_epsg)
            easts[ii] = edi.east
            norths[ii] = edi.north
            utm_zones[ii] = int(edi.utm_zone[:-1])

        if len(self.edi_list) == 0:
            raise IOError('Could not find and .edi file in {0}'.format(self.edi_path))

        if self.geoelectric_strike is None:
            try:
                # might try mode here instead of mean
                self.geoelectric_strike = np.median(np.nonzero(strike_angles))
            except:
                # empty list or so....
                # can happen, if everyhing is just 1D
                self.geoelectric_strike = 0.

        # need to check the zones of the stations
        main_utmzone = mode(utm_zones)[0][0]

        if self.model_epsg is None:
            for ii, zone in enumerate(utm_zones):
                if zone == main_utmzone:
                    continue
                else:
                    print(('station {0} is out of main utm zone'.format(self.edi_list[ii].station) + \
                           ' will not be included in profile'))

        # check regression for 2 profile orientations:
        # horizontal (N=N(E)) or vertical(E=E(N))
        # use the one with the lower standard deviation
        profile1 = sp.stats.linregress(easts, norths)
        profile2 = sp.stats.linregress(norths, easts)
        profile_line = profile1[:2]
        # if the profile is rather E=E(N), the parameters have to converted
        # into N=N(E) form:
        if profile2[4] < profile1[4]:
            profile_line = (1. / profile2[0], -profile2[1] / profile2[0])
        self.profile_line = profile_line
        # profile_line = sp.polyfit(lo_easts, lo_norths, 1)
        if self.profile_angle is None:
            self.profile_angle = (90 - (np.arctan(profile_line[0]) * 180 / np.pi)) % 180

        # rotate Z according to strike angle,

        # if strike was explicitely given, use that value!

        # otherwise:
        # have 90 degree ambiguity in strike determination
        # choose strike which offers larger angle with profile
        # if profile azimuth is in [0,90].

        if self._rotate_to_strike is False:
            if 0 <= self.profile_angle < 90:
                if np.abs(self.profile_angle - self.geoelectric_strike) < 45:
                    self.geoelectric_strike += 90
            elif 90 <= self.profile_angle < 135:
                if self.profile_angle - self.geoelectric_strike < 45:
                    self.geoelectric_strike -= 90
            else:
                if self.profile_angle - self.geoelectric_strike >= 135:
                    self.geoelectric_strike += 90

        self.geoelectric_strike = self.geoelectric_strike % 180

        # rotate components of Z and Tipper to align with geoelectric strike
        # which should now be perpendicular to the profile strike
        if self._rotate_to_strike == True:
            self.profile_angle = self.geoelectric_strike + 90
            p1 = np.tan(np.deg2rad(90 - self.profile_angle))
            # need to project the y-intercept to the new angle
            p2 = (self.profile_line[0] - p1) * easts[0] + self.profile_line[1]
            self.profile_line = (p1, p2)
        
            for edi in self.edi_list:
                edi.Z.rotate(self.geoelectric_strike - edi.Z.rotation_angle)
                # rotate tipper to profile azimuth, not strike.
                try:
                    edi.Tipper.rotate((self.profile_angle - 90) % 180 -
                                      edi.Tipper.rotation_angle.mean())
                except AttributeError:
                    edi.Tipper.rotate((self.profile_angle - 90) % 180 -
                                      edi.Tipper.rotation_angle)

            print('=' * 72)
            print(('Rotated Z and Tipper to align with '
                   '{0:+.2f} degrees E of N'.format(self.geoelectric_strike)))
            print(('Profile angle is '
                   '{0:+.2f} degrees E of N'.format(self.profile_angle)))
            print('=' * 72)
        else:
            for edi in self.edi_list:
                edi.Z.rotate((self.profile_angle - 90) % 180 - edi.Z.rotation_angle)
                # rotate tipper to profile azimuth, not strike.
                try:
                    edi.Tipper.rotate((self.profile_angle - 90) % 180 -
                                      edi.Tipper.rotation_angle.mean())
                except AttributeError:
                    edi.Tipper.rotate((self.profile_angle - 90) % 180 -
                                      edi.Tipper.rotation_angle)

            print('=' * 72)
            print(('Rotated Z and Tipper to be perpendicular  with '
                   '{0:+.2f} profile angle'.format((self.profile_angle - 90) % 180)))
            print(('Profile angle is '
                   '{0:+.2f} degrees E of N'.format(self.profile_angle)))
            print('=' * 72)

        # --> project stations onto profile line
        projected_stations = np.zeros((self.num_edi, 2))
        self.station_locations = np.zeros(self.num_edi)

        # create profile vector
        profile_vector = np.array([1, self.profile_line[0]])
        # be sure the amplitude is 1 for a unit vector
        profile_vector /= np.linalg.norm(profile_vector)
        for ii, edi in enumerate(self.edi_list):
            station_vector = np.array([easts[ii], norths[ii] - self.profile_line[1]])
            position = np.dot(profile_vector, station_vector) * profile_vector
            self.station_locations[ii] = np.linalg.norm(position)
            edi.offset = np.linalg.norm(position)
            edi.projected_east = position[0]
            edi.projected_north = position[1] + self.profile_line[1]
            projected_stations[ii] = [position[0], position[1] + self.profile_line[1]]
        # set the first station to 0
        for edi in self.edi_list:
            edi.offset -= self.station_locations.min()
        self.station_locations -= self.station_locations.min()

        # Sort from West to East:
        index_sort = np.argsort(self.station_locations)
        if self.profile_angle == 0:
            # Exception: sort from North to South
            index_sort = np.argsort(norths)

        # sorting along the profile
        self.edi_list = [self.edi_list[ii] for ii in index_sort]
        self.station_locations = np.array([self.station_locations[ii]
                                           for ii in index_sort])
        if self.estimate_elevation == True:
            self.project_elevation()

        self._profile_generated = True

    def project_elevation(self, elevation_model=None):
        """
        projects elevation data into the profile
        
        Arguments:
        -------------
            **elevation_model** : np.ndarray(3, num_elevation_points)
                                  (east, north, elevation)
                                  for now needs to be in utm coordinates
                                  if None then elevation is taken from edi_list
                                  
        Returns:
        ----------
            **elevation_profile** : 
        """
        self.elevation_model = elevation_model

        # --> get an elevation model for the mesh
        if self.elevation_model == None:
            self.elevation_profile = np.zeros((2, len(self.edi_list)))
            self.elevation_profile[0, :] = np.array([ss
                                                     for ss in self.station_locations])
            self.elevation_profile[1, :] = np.array([edi.elev
                                                     for edi in self.edi_list])

        # --> project known elevations onto the profile line
        else:
            self.elevation_profile = np.zeros((2, self.elevation_model.shape[1]))
            # create profile vector
            profile_vector = np.array([1, self.profile_line[0]])
            # be sure the amplitude is 1 for a unit vector
            profile_vector /= np.linalg.norm(profile_vector)
            for ii in range(self.elevation_model.shape[1]):
                east = self.elevation_model[0, ii]
                north = self.elevation_model[1, ii]
                elev = self.elevation_model[2, ii]
                elev_vector = np.array([east, north - self.profile_line[1]])
                position = np.dot(profile_vector, elev_vector) * profile_vector
                self.elevation_profile[0, ii] = np.linalg.norm(position)
                self.elevation_profile[1, ii] = elev

    def plot_profile(self, **kwargs):
        """
        Plot the projected profile line along with original station locations
        to make sure the line projected is correct.
        
        ===================== =================================================
        Key Words             Description          
        ===================== =================================================
        fig_dpi               dots-per-inch resolution of figure
                              *default* is 300
        fig_num               number if figure instance
                              *default* is 'Projected Profile'
        fig_size              size of figure in inches (width, height)
                              *default* is [5, 5]
        fs                    [ float ] font size in points of axes tick labels
                              axes labels are fs+2
                              *default* is 6
        lc                    [ string | (r, g, b) ]color of profile line 
                              (see matplotlib.line for options)
                              *default* is 'b' -- blue
        lw                    float, width of profile line in points
                              *default* is 1
        marker                [ string ] marker for stations 
                              (see matplotlib.pyplot.plot) for options  
        mc                    [ string | (r, g, b) ] color of projected 
                              stations.  *default* is 'k' -- black
        ms                    [ float ] size of station marker
                              *default* is 5
        station_id            [min, max] index values for station labels
                              *default* is None
        ===================== =================================================
        
        :Example: ::
            >>> edipath = r"/home/mt/edi_files"
            >>> pr = occam2d.Profile(edi_path=edipath)
            >>> pr.generate_profile()
            >>> # set station labels to only be from 1st to 4th index 
            >>> # of station name
            >>> pr.plot_profile(station_id=[0,4])
        
        """

        fig_num = kwargs.pop('fig_num', 'Projected Profile')
        fig_size = kwargs.pop('fig_size', [5, 5])
        fig_dpi = kwargs.pop('fig_dpi', 300)
        marker = kwargs.pop('marker', 'v')
        ms = kwargs.pop('ms', 5)
        mc = kwargs.pop('mc', 'k')
        lc = kwargs.pop('lc', 'b')
        lw = kwargs.pop('ls', 1)
        fs = kwargs.pop('fs', 6)
        station_id = kwargs.pop('station_id', None)

        plt.rcParams['figure.subplot.left'] = .12
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['font.size'] = fs

        if self._profile_generated is False:
            self.generate_profile()

        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        ax = fig.add_subplot(1, 1, 1, aspect='equal')

        for edi in self.edi_list:
            m1, = ax.plot(edi.projected_east, edi.projected_north,
                          marker=marker, ms=ms, mfc=mc, mec=mc, color=lc)

            m2, = ax.plot(edi.east, edi.north, marker=marker,
                          ms=.5 * ms, mfc=(.6, .6, .6), mec=(.6, .6, .6),
                          color=lc)

            if station_id is None:
                ax.text(edi.projected_east, edi.projected_north * 1.00025,
                        edi.station,
                        ha='center', va='baseline',
                        fontdict={'size': fs, 'weight': 'bold'})
            else:
                ax.text(edi.projected_east, edi.projected_north * 1.00025,
                        edi.station[station_id[0]:station_id[1]],
                        ha='center', va='baseline',
                        fontdict={'size': fs, 'weight': 'bold'})

        peasts = np.array([edi.projected_east for edi in self.edi_list])
        pnorths = np.array([edi.projected_north for edi in self.edi_list])
        easts = np.array([edi.east for edi in self.edi_list])
        norths = np.array([edi.north for edi in self.edi_list])

        ploty = sp.polyval(self.profile_line, easts)
        ax.plot(easts, ploty, lw=lw, color=lc)
        ax.set_title('Original/Projected Stations')
        ax.set_ylim((min([norths.min(), pnorths.min()]) * .999,
                     max([norths.max(), pnorths.max()]) * 1.001))
        ax.set_xlim((min([easts.min(), peasts.min()]) * .98,
                     max([easts.max(), peasts.max()]) * 1.02))
        ax.set_xlabel('Easting (m)',
                      fontdict={'size': fs + 2, 'weight': 'bold'})
        ax.set_ylabel('Northing (m)',
                      fontdict={'size': fs + 2, 'weight': 'bold'})
        ax.grid(alpha=.5)
        ax.legend([m1, m2], ['Projected', 'Original'], loc='upper left',
                  prop={'size': fs})
        plt.show()
        plt.savefig('/tmp/profile_angle0.png')


class Regularization(Mesh):
    """
    Creates a regularization grid based on Mesh.  Note that Mesh is inherited
    by Regularization, therefore the intended use is to build a mesh with 
    the Regularization class.
    
    The regularization grid is what Occam calculates the inverse model on.
    Setup is tricky and can be painful, as you can see it is not quite fully
    functional yet, as it cannot incorporate topography yet.  It seems like 
    you'd like to have the regularization setup so that your target depth is 
    covered well, in that the regularization blocks to this depth are 
    sufficiently small to resolve resistivity structure at that depth.  
    Finally, you want the regularization to go to a half space at the bottom, 
    basically one giant block.
    
    Arguments:
    -----------
        **station_locations** : np.ndarray(n_stations)
                                array of station locations along a profile
                                line in meters.
                                
    ======================= ===================================================
    Key Words/Attributes    Description    
    ======================= ===================================================
    air_key                 letter associated with the value of air 
                            *default* is 0
    air_value               value given to an air cell, *default* is 1E13
    binding_offset          offset from the right side of the furthest left
                            hand model block in meters.  The regularization
                            grid is setup such that this should be 0.
    cell_width              width of cells with in station area in meters
                            *default* is 100
    description             description of the model for the model file.
                            *default* is 'simple inversion'
    elevation_profile       elevation profile along the profile line.
                            given as np.ndarray(nx, 2), where the elements
                            are x_location, elevation.  If elevation profile
                            is given add_elevation is called automatically.
                            *default* is None
    mesh_fn                 full path to mesh file.
    mesh_values             letter values of each triangular mesh element
                            if the cell is free value is ?
    model_columns
    model_name
    model_rows
    
    min_block_width         [ float ] minimum model block width in meters, 
                            *default* is 2*cell_width
    n_layers                number of vertical layers in mesh
                            *default* is 90 
    num_free_param          [ int ] number of free parameters in the model.
                            this is a tricky number to estimate apparently. 
    num_layers              [ int ] number of regularization layers.
    num_x_pad_cells         number of horizontal padding cells outside the
                            the station area that will increase in size
                            by x_pad_multiplier. *default* is 7
    num_x_pad_small_cells   number of horizonal padding cells just outside 
                            the station area with width cell_width.  This is 
                            to extend the station area if needed.  
                            *default* is 2 
    num_z_pad_cells         number of vertical padding cells below 
                            z_target_depth down to z_bottom. *default* is 5
    prejudice_fn            full path to prejudice file 
                            *default* is 'none'
    reg_basename            basename of regularization file (model file)
                            *default* is 'Occam2DModel'
    reg_fn                  full path to regularization file (model file)
                            *default* is save_path/reg_basename
    rel_station_locations   relative station locations within the mesh.  The
                            locations are relative to the center of the station
                            area.  *default* is None, filled later
    save_path               full path to save mesh and model file to. 
                            *default* is current working directory.
    statics_fn              full path to static shift file
                            Static shifts in occam may not work.
                            *default* is 'none'
    station_locations       location of stations in meters, can be on a 
                            relative grid or in UTM.
    trigger                 [ float ] multiplier to merge model blocks at 
                            depth.  A higher number increases the number of
                            model blocks at depth.  *default* is .65
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
        
    .. note:: regularization does not work with topography yet.  Having 
              problems calculating the number of free parameters.
    
    ========================= =================================================
    Methods                   Description
    ========================= =================================================
    add_elevation             adds elevation to the mesh given elevation
                              profile.  
    build_mesh                builds the mesh given the attributes of Mesh.  If
                              elevation_profile is not None, add_elevation is
                              called inside build_mesh
    build_regularization      builds the regularization grid from the build mesh
                              be sure to plot the grids before starting the
                              inversion to make sure coverage is appropriate.
    get_num_free_param        estimate the number of free parameters.  
                              **This is a work in progress**
    plot_mesh                 plots the built mesh with station location.  
    read_mesh_file            reads in an existing mesh file and populates the 
                              appropriate attributes.
    read_regularization_file  read in existing regularization file, populates
                              apporopriate attributes  
    write_mesh_file           writes a mesh file to save_path
    write_regularization_file writes a regularization file
    ======================= ===================================================
    
    :Example: ::

        >>> edipath = r"/home/mt/edi_files"
        >>> profile = occam2d.Profile(edi_path=edi_path)
        >>> profile.generate_profile()
        >>> reg = occam2d.Regularization(profile.station_locations)
        >>> reg.build_mesh()
        >>> reg.build_regularization()
        >>> reg.save_path = r"/home/occam2d/Line1/Inv1"
        >>> reg.write_regularization_file()
    
    """

    def __init__(self, station_locations=None, **kwargs):
        # Be sure to initialize Mesh        
        Mesh.__init__(self, station_locations, **kwargs)

        self.min_block_width = kwargs.pop('min_block_width',
                                          2 * np.median(self.cell_width))
        self.trigger = kwargs.pop('trigger', .75)
        self.model_columns = None
        self.model_rows = None
        self.binding_offset = None
        self.reg_fn = None
        self.reg_basename = 'Occam2DModel'
        self.model_name = 'model made by mtpy.modeling.occam2d'
        self.description = 'simple Inversion'
        self.num_param = None
        self.num_free_param = None
        self.statics_fn = kwargs.pop('statics_fn', 'none')
        self.prejudice_fn = kwargs.pop('prejudice_fn', 'none')
        self.num_layers = kwargs.pop('num_layers', None)

        # --> build mesh
        if self.station_locations is not None:
            self.build_mesh()
            self.build_regularization()

    def build_regularization(self):
        """
        Builds larger boxes around existing mesh blocks for the regularization.
        As the model deepens the regularization boxes get larger.  
        
        The regularization boxes are merged mesh cells as prescribed by the
        Occam method.
    
        """
        # list of the mesh columns to combine
        self.model_columns = []
        # list of mesh rows to combine
        self.model_rows = []

        # At the top of the mesh model blocks will be 2 combined mesh blocks
        # Note that the padding cells are combined into one model block
        var = (self.x_nodes.shape[0] - 2 * self.num_x_pad_cells)
        print ("******* var=", var)

        station_col = [2] *int((self.x_nodes.shape[0] - 2 * self.num_x_pad_cells) / 2)
        model_cols = [self.num_x_pad_cells] + station_col + [self.num_x_pad_cells]
        station_widths = [self.x_nodes[ii] + self.x_nodes[ii + 1] for ii in
                          range(self.num_x_pad_cells,
                                self.x_nodes.shape[0] - self.num_x_pad_cells, 2)]

        pad_width = self.x_nodes[0:self.num_x_pad_cells].sum()
        model_widths = [pad_width] + station_widths + [pad_width]
        num_cols = len(model_cols)

        #        model_thickness = np.append(self.z_nodes[0:self.z_nodes.shape[0]-
        #                                                        self.num_z_pad_cells],
        #                                    self.z_nodes[-self.num_z_pad_cells:].sum())
        model_thickness = np.hstack([[self.z_nodes[:2].sum()],
                                     self.z_nodes[2:-self.num_z_pad_cells],
                                     [self.z_nodes[-self.num_z_pad_cells:].sum()]])

        self.num_param = 0
        # --> now need to calulate model blocks to the bottom of the model
        columns = list(model_cols)
        widths = list(model_widths)
        for zz, thickness in enumerate(model_thickness):
            # index for model column blocks from first_row, start at 1 because
            # 0 is for padding cells            
            block_index = 1
            num_rows = 1
            if zz == 0:
                num_rows += 1
            if zz == len(model_thickness) - 1:
                num_rows = self.num_z_pad_cells
            while block_index + 1 < num_cols - 1:
                # check to see if horizontally merged mesh cells are not larger
                # than the thickness times trigger
                if thickness < self.trigger * (widths[block_index] + \
                                                       widths[block_index + 1]):
                    block_index += 1
                    continue
                # merge 2 neighboring cells to avoid vertical exaggerations
                else:
                    widths[block_index] += widths[block_index + 1]
                    columns[block_index] += columns[block_index + 1]
                    # remove one of the merged cells
                    widths.pop(block_index + 1)
                    columns.pop(block_index + 1)

                    num_cols -= 1
            self.num_param += num_cols

            self.model_columns.append(list(columns))
            self.model_rows.append([num_rows, num_cols])

        # calculate the distance from the right side of the furthest left
        # model block to the furthest left station which is half the distance
        # from the center of the mesh grid.
        self.binding_offset = self.x_grid[self.num_x_pad_cells + 1] + \
                              self.station_locations.mean()

        self.get_num_free_params()

        print('=' * 55)
        print('{0:^55}'.format('regularization parameters'.upper()))
        print('=' * 55)
        print('   binding offset       = {0:.1f}'.format(self.binding_offset))
        print('   number layers        = {0}'.format(len(self.model_columns)))
        print('   number of parameters = {0}'.format(self.num_param))
        print('   number of free param = {0}'.format(self.num_free_param))
        print('=' * 55)

    def get_num_free_params(self):
        """
        estimate the number of free parameters in model mesh.
        
        I'm assuming that if there are any fixed parameters in the block, then
        that model block is assumed to be fixed. Not sure if this is right
        cause there is no documentation.
        
        **DOES NOT WORK YET**
        """

        self.num_free_param = 0

        row_count = 0
        # loop over columns and rows of regularization grid
        for col, row in zip(self.model_columns, self.model_rows):
            rr = row[0]
            col_count = 0
            for ii, cc in enumerate(col):
                # make a model block from the index values of the regularization
                # grid
                model_block = self.mesh_values[row_count:row_count + rr,
                              col_count:col_count + cc, :]

                # find all the free triangular blocks within that model block
                find_free = np.where(model_block == '?')
                try:
                    # test to see if the number of free parameters is equal
                    # to the number of triangular elements with in the model
                    # block, if there is the model block is assumed to be free.
                    if find_free[0].size == model_block.size:
                        self.num_free_param += 1
                except IndexError:
                    pass
                col_count += cc
            row_count += rr

    def write_regularization_file(self, reg_fn=None, reg_basename=None,
                                  statics_fn='none', prejudice_fn='none',
                                  save_path=None):
        """
        Write a regularization file for input into occam.
        
        Calls build_regularization if build_regularization has not already
        been called.
        
        if reg_fn is None, then file is written to save_path/reg_basename
        
        Arguments:
        ----------
            **reg_fn** : string
                         full path to regularization file. *default* is None
                         and file will be written to save_path/reg_basename
                         
            **reg_basename** : string
                               basename of regularization file
            
            **statics_fn** : string
                             full path to static shift file
                             .. note:: static shift does not always work in
                                       occam2d.exe
            **prejudice_fn** : string
                               full path to prejudice file
            
            **save_path** : string
                            path to save regularization file.
                            *default* is current working directory
                                
        """
        if save_path is not None:
            self.save_path = save_path
        if reg_basename is not None:
            self.reg_basename = reg_basename
        if reg_fn is None:
            if self.save_path is None:
                self.save_path = os.getcwd()
            self.reg_fn = os.path.join(self.save_path, self.reg_basename)

        self.statics_fn = statics_fn
        self.prejudice_fn = prejudice_fn

        if self.model_columns is None:
            if self.binding_offset is None:
                self.build_mesh()
            self.build_regularization()

        reg_lines = []

        # --> write out header information
        reg_lines.append('{0:<18}{1}\n'.format('Format:',
                                               'occam2mtmod_1.0'.upper()))
        reg_lines.append('{0:<18}{1}\n'.format('Model Name:',
                                               self.model_name.upper()))
        reg_lines.append('{0:<18}{1}\n'.format('Description:',
                                               self.description.upper()))
        if os.path.dirname(self.mesh_fn) == self.save_path:
            reg_lines.append('{0:<18}{1}\n'.format('Mesh File:',
                                                   os.path.basename(self.mesh_fn)))
        else:
            reg_lines.append('{0:<18}{1}\n'.format('Mesh File:', self.mesh_fn))
        reg_lines.append('{0:<18}{1}\n'.format('Mesh Type:',
                                               'pw2d'.upper()))
        if os.path.dirname(self.statics_fn) == self.save_path:
            reg_lines.append('{0:<18}{1}\n'.format('Statics File:',
                                                   os.path.basename(self.statics_fn)))
        else:
            reg_lines.append('{0:<18}{1}\n'.format('Statics File:',
                                                   self.statics_fn))
        if os.path.dirname(self.prejudice_fn) == self.save_path:
            reg_lines.append('{0:<18}{1}\n'.format('Prejudice File:',
                                                   os.path.basename(self.prejudice_fn)))
        else:
            reg_lines.append('{0:<18}{1}\n'.format('Prejudice File:',
                                                   self.prejudice_fn))
        reg_lines.append('{0:<20}{1: .1f}\n'.format('Binding Offset:',
                                                    self.binding_offset))
        reg_lines.append('{0:<20}{1}\n'.format('Num Layers:',
                                               len(self.model_columns)))

        # --> write rows and columns of regularization grid
        for row, col in zip(self.model_rows, self.model_columns):
            reg_lines.append(''.join([' {0:>5}'.format(rr) for rr in row]) + '\n')
            reg_lines.append(''.join(['{0:>5}'.format(cc) for cc in col]) + '\n')

        reg_lines.append('{0:<18}{1}\n'.format('NO. EXCEPTIONS:', '0'))
        with open(self.reg_fn, 'w') as rfid:
            rfid.writelines(reg_lines)
            # rfid.close()

        print('Wrote Regularization file to {0}'.format(self.reg_fn))

    def read_regularization_file(self, reg_fn):
        """
        Read in a regularization file and populate attributes:
            * binding_offset            
            * mesh_fn    
            * model_columns
            * model_rows
            * prejudice_fn
            * statics_fn
        
        """
        self.reg_fn = reg_fn
        self.save_path = os.path.dirname(reg_fn)

        rfid = open(self.reg_fn, 'r')

        self.model_rows = []
        self.model_columns = []
        ncols = []

        rlines = rfid.readlines()

        for ii, iline in enumerate(rlines):
            # read header information
            if iline.find(':') > 0:
                iline = iline.strip().split(':')
                key = iline[0].strip().lower()
                key = key.replace(' ', '_').replace('file', 'fn')
                value = iline[1].strip()
                try:
                    setattr(self, key, float(value))
                except ValueError:
                    setattr(self, key, value)

                # append the last line
                if key.find('exception') > 0:
                    self.model_columns.append(ncols)

            # get mesh values
            else:
                iline = iline.strip().split()
                iline = [int(jj) for jj in iline]
                if len(iline) == 2:
                    if len(ncols) > 0:
                        self.model_columns.append(ncols)
                    self.model_rows.append(iline)
                    ncols = []
                elif len(iline) > 2:
                    ncols = ncols + iline

        # set mesh file name
        if not os.path.isfile(self.mesh_fn):
            self.mesh_fn = os.path.join(self.save_path, self.mesh_fn)

        # set statics file name
        if not os.path.isfile(self.mesh_fn):
            self.statics_fn = os.path.join(self.save_path, self.statics_fn)

        # set prejudice file name
        if not os.path.isfile(self.mesh_fn):
            self.prejudice_fn = os.path.join(self.save_path, self.prejudice_fn)


class Startup(object):
    """
    Reads and writes the startup file for Occam2D.
    
    .. note:: Be sure to look at the Occam 2D documentation for description
              of all parameters
    
    ========================= =================================================
    Key Words/Attributes      Description
    ========================= =================================================
    data_fn                   full path to data file
    date_time                 date and time the startup file was written
    debug_level               [ 0 | 1 | 2 ] see occam documentation
                              *default* is 1
    description               brief description of inversion run
                              *default* is 'startup created by mtpy'  
    diagonal_penalties        penalties on diagonal terms
                              *default* is 0
    format                    Occam file format
                              *default* is 'OCCAMITER_FLEX'
    iteration                 current iteration number
                              *default* is 0
    iterations_to_run         maximum number of iterations to run
                              *default* is 20
    lagrange_value            starting lagrange value
                              *default* is 5
    misfit_reached            [ 0 | 1 ] 0 if misfit has been reached, 1 if it
                              has.  *default* is 0
    misfit_value              current misfit value.  *default* is 1000
    model_fn                  full path to model file
    model_limits              limits on model resistivity values
                              *default* is None
    model_value_steps         limits on the step size of model values
                              *default* is None
    model_values              np.ndarray(num_free_params) of model values
    param_count               number of free parameters in model
    resistivity_start         starting resistivity value.  If model_values is
                              not given, then all values with in model_values
                              array will be set to resistivity_start
    roughness_type            [ 0 | 1 | 2 ] type of roughness
                              *default* is 1
    roughness_value           current roughness value.  
                              *default* is 1E10
    save_path                 directory path to save startup file to
                              *default* is current working directory  
    startup_basename          basename of startup file name. 
                              *default* is Occam2DStartup
    startup_fn                full path to startup file.
                              *default* is save_path/startup_basename  
    stepsize_count            max number of iterations per step
                              *default* is 8
    target_misfit             target misfit value.
                              *default* is 1.
    ========================= =================================================
    
    :Example: ::
    
        >>> startup = occam2d.Startup()
        >>> startup.data_fn = ocd.data_fn
        >>> startup.model_fn = profile.reg_fn
        >>> startup.param_count = profile.num_free_params
        >>> startup.save_path = r"/home/occam2d/Line1/Inv1"
    """

    def __init__(self, **kwargs):
        self.save_path = kwargs.pop('save_path', None)
        self.startup_basename = kwargs.pop('startup_basename', 'Occam2DStartup')
        self.startup_fn = kwargs.pop('startup_fn', None)
        self.model_fn = kwargs.pop('model_fn', None)
        self.data_fn = kwargs.pop('data_fn', None)
        self.format = kwargs.pop('format', 'OCCAMITER_FLEX')
        self.date_time = kwargs.pop('date_time', time.ctime())
        self.description = kwargs.pop('description', 'startup created by mtpy')
        self.iterations_to_run = kwargs.pop('iterations_to_run', 20)
        self.roughness_type = kwargs.pop('roughness_type', 1)
        self.target_misfit = kwargs.pop('target_misfit', 1.0)
        self.diagonal_penalties = kwargs.pop('diagonal_penalties', 0)
        self.stepsize_count = kwargs.pop('stepsize_count', 8)
        self.model_limits = kwargs.pop('model_limits', None)
        self.model_value_steps = kwargs.pop('model_value_steps', None)
        self.debug_level = kwargs.pop('debug_level', 1)
        self.iteration = kwargs.pop('iteration', 0)
        self.lagrange_value = kwargs.pop('lagrange_value', 5.0)
        self.roughness_value = kwargs.pop('roughness_value', 1e10)
        self.misfit_value = kwargs.pop('misfit_value', 1000)
        self.misfit_reached = kwargs.pop('misfit_reached', 0)
        self.param_count = kwargs.pop('param_count', None)
        self.resistivity_start = kwargs.pop('resistivity_start', 2)
        self.model_values = kwargs.pop('model_values', None)

    def write_startup_file(self, startup_fn=None, save_path=None,
                           startup_basename=None):
        """
        Write a startup file based on the parameters of startup class.  
        Default file name is save_path/startup_basename
        
        Arguments:
        -----------
            **startup_fn** : string
                             full path to startup file. *default* is None
            
            **save_path** : string
                            directory to save startup file. *default* is None
                            
            **startup_basename** : string
                                   basename of starup file. *default* is None
        
        """
        if save_path is not None:
            self.save_path = save_path

        if self.save_path is None:
            self.save_path = os.path.dirname(self.data_fn)
        if startup_basename is not None:
            self.startup_basename = startup_basename

        if startup_fn is None:
            self.startup_fn = os.path.join(self.save_path,
                                           self.startup_basename)

        # --> check to make sure all the important input are given
        if self.data_fn is None:
            raise OccamInputError('Need to input data file name')

        if self.model_fn is None:
            raise OccamInputError('Need to input model/regularization file name')

        if self.param_count is None:
            raise OccamInputError('Need to input number of model parameters')

        slines = []
        slines.append('{0:<20}{1}\n'.format('Format:', self.format))
        slines.append('{0:<20}{1}\n'.format('Description:', self.description))
        if os.path.dirname(self.model_fn) == self.save_path:
            slines.append('{0:<20}{1}\n'.format('Model File:',
                                                os.path.basename(self.model_fn)))
        else:
            slines.append('{0:<20}{1}\n'.format('Model File:', self.model_fn))
        if os.path.dirname(self.data_fn) == self.save_path:
            slines.append('{0:<20}{1}\n'.format('Data File:',
                                                os.path.basename(self.data_fn)))
        else:
            slines.append('{0:<20}{1}\n'.format('Data File:', self.data_fn))
        slines.append('{0:<20}{1}\n'.format('Date/Time:', self.date_time))
        slines.append('{0:<20}{1}\n'.format('Iterations to run:',
                                            self.iterations_to_run))
        slines.append('{0:<20}{1}\n'.format('Target Misfit:',
                                            self.target_misfit))
        slines.append('{0:<20}{1}\n'.format('Roughness Type:',
                                            self.roughness_type))
        slines.append('{0:<20}{1}\n'.format('Diagonal Penalties:',
                                            self.diagonal_penalties))
        slines.append('{0:<20}{1}\n'.format('Stepsize Cut Count:',
                                            self.stepsize_count))
        if self.model_limits is None:
            slines.append('{0:<20}{1}\n'.format('!Model Limits:', 'none'))
        else:
            slines.append('{0:<20}{1},{2}\n'.format('Model Limits:',
                                                    self.model_limits[0],
                                                    self.model_limits[1]))
        if self.model_value_steps is None:
            slines.append('{0:<20}{1}\n'.format('!Model Value Steps:', 'none'))
        else:
            slines.append('{0:<20}{1}\n'.format('Model Value Steps:',
                                                self.model_value_steps))
        slines.append('{0:<20}{1}\n'.format('Debug Level:', self.debug_level))
        slines.append('{0:<20}{1}\n'.format('Iteration:', self.iteration))
        slines.append('{0:<20}{1}\n'.format('Lagrange Value:',
                                            self.lagrange_value))
        slines.append('{0:<20}{1}\n'.format('Roughness Value:',
                                            self.roughness_value))
        slines.append('{0:<20}{1}\n'.format('Misfit Value:', self.misfit_value))
        slines.append('{0:<20}{1}\n'.format('Misfit Reached:',
                                            self.misfit_reached))
        slines.append('{0:<20}{1}\n'.format('Param Count:', self.param_count))

        # make an array of starting values if not are given
        if self.model_values is None:
            self.model_values = np.zeros(self.param_count)
            self.model_values[:] = self.resistivity_start

        if self.model_values.shape[0] != self.param_count:
            raise OccamInputError('length of model vaues array is not equal '
                                  'to param count {0} != {1}'.format(
                self.model_values.shape[0], self.param_count))

        # write out starting resistivity values
        sline = []
        for ii, mv in enumerate(self.model_values):
            sline.append('{0:^10.4f}'.format(mv))
            if np.remainder(ii + 1, 4) == 0:
                sline.append('\n')
                slines.append(''.join(list(sline)))
                sline = []
        slines.append(''.join(list(sline + ['\n'])))
        # --> write file
        with open(self.startup_fn, 'w') as sfid:
            sfid.writelines(slines)
            # sfid.close()

        print('Wrote Occam2D startup file to {0}'.format(self.startup_fn))


# ------------------------------------------------------------------------------
class Data(Profile):
    """
    Reads and writes data files and more.  
    
    Inherets Profile, so the intended use is to use Data to project stations 
    onto a profile, then write the data file.  
    
    ===================== =====================================================
    Model Modes           Description                     
    ===================== =====================================================
    1 or log_all          Log resistivity of TE and TM plus Tipper
    2 or log_te_tip       Log resistivity of TE plus Tipper
    3 or log_tm_tip       Log resistivity of TM plus Tipper
    4 or log_te_tm        Log resistivity of TE and TM
    5 or log_te           Log resistivity of TE
    6 or log_tm           Log resistivity of TM
    7 or all              TE, TM and Tipper
    8 or te_tip           TE plus Tipper
    9 or tm_tip           TM plus Tipper
    10 or te_tm           TE and TM mode
    11 or te              TE mode
    12 or tm              TM mode
    13 or tip             Only Tipper
    ===================== =====================================================
    
    
    **data** : is a list of dictioinaries containing the data for each station.
               keys include:
                   * 'station' -- name of station
                   * 'offset' -- profile line offset
                   * 'te_res' -- TE resisitivity in linear scale
                   * 'tm_res' -- TM resistivity in linear scale
                   * 'te_phase' -- TE phase in degrees
                   * 'tm_phase' --  TM phase in degrees in first quadrant
                   * 're_tip' -- real part of tipper along profile
                   * 'im_tip' -- imaginary part of tipper along profile
                   
               each key is a np.ndarray(2, num_freq)
               index 0 is for data
               index 1 is for error
    
    ===================== =====================================================
    Key Words/Attributes  Desctription
    ===================== =====================================================
    _data_header          header line in data file
    _data_string          full data string
    _profile_generated    [ True | False ] True if profile has already been
                          generated.
    _rotate_to_strike     [ True | False ] True to rotate data to strike
                          angle.  *default* is True
    data                  list of dictionaries of data for each station.
                          see above
    data_fn               full path to data file
    data_list             list of lines to write to data file
    edi_list              list of mtpy.core.mt instances for each .edi file
                          read
    edi_path              directory path where .edi files are
    edi_type              [ 'z' | 'spectra' ] for .edi format  
    elevation_model       model elevation np.ndarray(east, north, elevation) 
                          in meters
    elevation_profile     elevation along profile np.ndarray (x, elev) (m)
    fn_basename           data file basename *default* is OccamDataFile.dat
    freq                  list of frequencies to use for the inversion
    freq_max              max frequency to use in inversion. *default* is None
    freq_min              min frequency to use in inversion. *default* is None
    freq_num              number of frequencies to use in inversion
    geoelectric_strike    geoelectric strike angle assuming N = 0, E = 90
    masked_data           similar to data, but any masked points are now 0
    mode_dict             dictionary of model modes to chose from
    model_mode            model mode to use for inversion, see above
    num_edi               number of stations to invert for
    occam_dict            dictionary of occam parameters to use internally
    occam_format          occam format of data file.  
                          *default* is OCCAM2MTDATA_1.0
    phase_te_err          percent error in phase for TE mode. *default* is 5
    phase_tm_err          percent error in phase for TM mode. *default* is 5 
    profile_angle         angle of profile line realtive to N = 0, E = 90
    profile_line          m, b coefficients for mx+b definition of profile line
    res_te_err            percent error in resistivity for TE mode. 
                          *default* is 10
    res_tm_err            percent error in resistivity for TM mode.
                          *default* is 10
    save_path             directory to save files to
    station_list          list of station for inversion
    station_locations     station locations along profile line
    tipper_err            percent error in tipper. *default* is 5
    title                 title in data file.  
    ===================== =====================================================
    
    =========================== ===============================================
    Methods                     Description
    =========================== ===============================================
    _fill_data                  fills the data array that is described above    
    _get_data_list              gets the lines to write to data file
    _get_frequencies            gets frequency list to invert for 
    get_profile_origin          get profile origin in UTM coordinates
    mask_points                 masks points in data picked from 
                                plot_mask_points  
    plot_mask_points            plots data responses to interactively mask
                                data points.
    plot_resonse                plots data/model responses, returns 
                                PlotResponse data type.
    read_data_file              read in existing data file and fill appropriate
                                attributes.
    write_data_file             write a data file according to Data attributes
    =========================== ===============================================
    
    :Example Write Data File: ::
        >>> import mtpy.modeling.occam2d as occam2d
        >>> edipath = r"/home/mt/edi_files"
        >>> slst = ['mt{0:03}'.format(ss) for ss in range(1, 20)]
        >>> ocd = occam2d.Data(edi_path=edipath, station_list=slst)
        >>> # model just the tm mode and tipper
        >>> ocd.model_mode = 3
        >>> ocd.save_path = r"/home/occam/Line1/Inv1"
        >>> ocd.write_data_file()
        >>> # mask points
        >>> ocd.plot_mask_points()
        >>> ocd.mask_points()

    """

    def __init__(self, edi_path=None, **kwargs):
        Profile.__init__(self, edi_path, **kwargs)

        self.data_fn = kwargs.pop('data_fn', None)
        self.fn_basename = kwargs.pop('fn_basename', 'OccamDataFile.dat')
        self.save_path = kwargs.pop('save_path', None)
        self.freq = kwargs.pop('freq', None)
        self.model_mode = kwargs.pop('model_mode', '1')
        self.data = kwargs.pop('data', None)
        self.data_list = None
        self.model_epsg = kwargs.pop('model_epsg',None)

        self.res_te_err = kwargs.pop('res_te_err', 10)
        self.res_tm_err = kwargs.pop('res_tm_err', 10)
        self.phase_te_err = kwargs.pop('phase_te_err', 5)
        self.phase_tm_err = kwargs.pop('phase_tm_err', 5)
        self.tipper_err = kwargs.pop('tipper_err', 10)
        self.error_type = 'floor' # 'floor' or 'value'

        self.freq_min = kwargs.pop('freq_min', None)
        self.freq_max = kwargs.pop('freq_max', None)
        self.freq_num = kwargs.pop('freq_num', None)
        self.freq_tol = kwargs.pop('freq_tol', None)

        self.occam_format = 'OCCAM2MTDATA_1.0'
        self.title = 'MTpy-OccamDatafile'
        self.edi_type = 'z'
        self.masked_data = None

        self.occam_dict = {'1': 'log_te_res',
                           '2': 'te_phase',
                           '3': 're_tip',
                           '4': 'im_tip',
                           '5': 'log_tm_res',
                           '6': 'tm_phase',
                           '9': 'te_res',
                           '10': 'tm_res'}

        self.mode_dict = {'log_all': [1, 2, 3, 4, 5, 6],
                          'log_te_tip': [1, 2, 3, 4],
                          'log_tm_tip': [5, 6, 3, 4],
                          'log_te_tm': [1, 2, 5, 6],
                          'log_te': [1, 2],
                          'log_tm': [5, 6],
                          'all': [9, 2, 3, 4, 10, 6],
                          'te_tip': [9, 2, 3, 4],
                          'tm_tip': [10, 6, 3, 4],
                          'te_tm': [9, 2, 10, 6],
                          'te': [9, 2],
                          'tm': [10, 6],
                          'tip': [3, 4],
                          '1': [1, 2, 3, 4, 5, 6],
                          '2': [1, 2, 3, 4],
                          '3': [5, 6, 3, 4],
                          '4': [1, 2, 5, 6],
                          '5': [1, 2],
                          '6': [5, 6],
                          '7': [9, 2, 3, 4, 10, 6],
                          '8': [9, 2, 3, 4],
                          '9': [10, 6, 3, 4],
                          '10': [9, 2, 10, 6],
                          '11': [9, 2],
                          '12': [10, 6],
                          '13': [3, 4]}

        self._data_string = '{0:^6}{1:^6}{2:^6} {3: >8} {4: >8}\n'
        self._data_header = '{0:<6}{1:<6}{2:<6} {3:<8} {4:<8}\n'.format(
            'SITE', 'FREQ', 'TYPE', 'DATUM', 'ERROR')

    def read_data_file(self, data_fn=None):
        """
        Read in an existing data file and populate appropriate attributes
            * data
            * data_list
            * freq
            * station_list
            * station_locations
            
        Arguments:
        -----------
            **data_fn** : string
                          full path to data file
                          *default* is None and set to save_path/fn_basename
                
        :Example: ::
            
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Data()
            >>> ocd.read_data_file(r"/home/Occam2D/Line1/Inv1/Data.dat")
            
        """

        if data_fn is not None:
            self.data_fn = data_fn

        if os.path.isfile(self.data_fn) == False:
            raise OccamInputError('Could not find {0}'.format(self.data_fn))
        if self.data_fn is None:
            raise OccamInputError('data_fn is None, input filename')

        self.save_path = op.dirname(self.data_fn)

        print('Reading from {0}'.format(self.data_fn))

        dfid = open(self.data_fn, 'r')

        dlines = dfid.readlines()

        # get format of input data
        self.occam_format = dlines[0].strip().split(':')[1].strip()

        # get title
        title_str = dlines[1].strip().split(':')[1].strip()

        title_list = title_str.split(',')
        self.title = title_list[0]

        # get strike angle and profile angle
        if len(title_list) > 1:
            for t_str in title_list[1:]:
                t_list = t_str.split('=')
                if len(t_list) > 1:
                    key = t_list[0].strip().lower().replace(' ', '_')
                    if key == 'profile':
                        key = 'profile_angle'
                    elif key == 'strike':
                        key = 'geoelectric_strike'
                    value = t_list[1].split('deg')[0].strip()
                    print('    {0} = {1}'.format(key, value))
                    try:
                        setattr(self, key, float(value))
                    except ValueError:
                        setattr(self, key, value)

        # get number of sites
        nsites = int(dlines[2].strip().split(':')[1].strip())
        print('    {0} = {1}'.format('number of sites', nsites))

        # get station names
        self.station_list = np.array([dlines[ii].strip()
                                      for ii in range(3, nsites + 3)])

        # get offsets in meters
        self.station_locations = np.array([float(dlines[ii].strip())
                                           for ii in range(4 + nsites, 4 + 2 * nsites)])

        # get number of frequencies
        nfreq = int(dlines[4 + 2 * nsites].strip().split(':')[1].strip())
        print('    {0} = {1}'.format('number of frequencies', nfreq))

        # get frequencies
        self.freq = np.array([float(dlines[ii].strip())
                              for ii in range(5 + 2 * nsites, 5 + 2 * nsites + nfreq)])

        # get periods
        self.period = 1. / self.freq

        # -----------get data-------------------
        # set zero array size the first row will be the data and second the error
        asize = (2, self.freq.shape[0])

        # make a list of dictionaries for each station.
        self.data = [{'station': station,
                      'offset': offset,
                      'te_phase': np.zeros(asize),
                      'tm_phase': np.zeros(asize),
                      're_tip': np.zeros(asize),
                      'im_tip': np.zeros(asize),
                      'te_res': np.zeros(asize),
                      'tm_res': np.zeros(asize)}
                     for station, offset in zip(self.station_list,
                                                self.station_locations)]

        self.data_list = dlines[7 + 2 * nsites + nfreq:]
        for line in self.data_list:
            try:
                station, freq, comp, odata, oerr = line.split()
                # station index -1 cause python starts at 0
                ss = int(station) - 1

                # frequency index -1 cause python starts at 0
                ff = int(freq) - 1
                # data key
                key = self.occam_dict[comp]

                # put into array
                if int(comp) == 1 or int(comp) == 5:
                    self.data[ss][key[4:]][0, ff] = 10 ** float(odata)
                    # error
                    self.data[ss][key[4:]][1, ff] = float(oerr) * np.log(10)
                else:
                    self.data[ss][key][0, ff] = float(odata)
                    # error
                    self.data[ss][key][1, ff] = float(oerr)
            except ValueError:
                print('Could not read line {0}'.format(line))

    def _get_frequencies(self):
        """
        from the list of edi's get a frequency list to invert for.
        
        Uses Attributes:
        ------------
            **freq_min** : float (Hz)
                           minimum frequency to invert for.
                           *default* is None and will use the data to find min
            
            **freq_max** : float (Hz)
                           maximum frequency to invert for
                           *default* is None and will use the data to find max
                           
            **freq_num** : int
                           number of frequencies to invert for
                           *default* is None and will use the data to find num
        """

        if self.freq is not None:
            return

        # get all frequencies from all edi files
        lo_all_freqs = []
        for edi in self.edi_list:
            lo_all_freqs.extend(list(edi.Z.freq))

        # sort all frequencies so that they are in descending order,
        # use set to remove repeats and make an array
        all_freqs = np.array(sorted(list(set(lo_all_freqs)), reverse=True))

        # --> get min and max values if none are given
        if (self.freq_min is None) or (self.freq_min < all_freqs.min()) or \
                (self.freq_min > all_freqs.max()):
            self.freq_min = all_freqs.min()

        if (self.freq_max is None) or (self.freq_max > all_freqs.max()) or \
                (self.freq_max < all_freqs.max()):
            self.freq_max = all_freqs.max()

        # --> get all frequencies within the given range
        self.freq = all_freqs[np.where((all_freqs >= self.freq_min) &
                                       (all_freqs <= self.freq_max))]

        if len(self.freq) == 0:
            raise OccamInputError('No frequencies in user-defined interval '
                                  '[{0}, {1}]'.format(self.freq_min, self.freq_max))

        # check, if frequency list is longer than given max value
        if self.freq_num is not None:
            if int(self.freq_num) < self.freq.shape[0]:
                print(('Number of frequencies exceeds freq_num '
                       '{0} > {1} '.format(self.freq.shape[0], self.freq_num) +
                       'Trimming frequencies to {0}'.format(self.freq_num)))

                excess = self.freq.shape[0] / float(self.freq_num)
                if excess < 2:
                    offset = 0
                else:
                    stepsize = (self.freq.shape[0] - 1) / self.freq_num
                    offset = stepsize / 2.
                indices = np.array(np.around(np.linspace(offset,
                                                         self.freq.shape[0] - 1 - offset,
                                                         self.freq_num), 0), dtype='int')
                if indices[0] > (self.freq.shape[0] - 1 - indices[-1]):
                    indices -= 1
                self.freq = self.freq[indices]

    def _fill_data(self):
        """
        Read all Edi files. 
        Create a profile
        rotate impedance and tipper
        Extract frequencies. 

        Collect all information sorted according to occam specifications.

        Data of Z given in muV/m/nT = km/s
        Error is assumed to be 1 stddev.
        """

        # create a profile line, this sorts the stations by offset and rotates
        # data.
        self.generate_profile()
        self.plot_profile()

        # --> get frequencies to invert for
        self._get_frequencies()

        # set zero array size the first row will be the data and second the error
        asize = (2, self.freq.shape[0])

        # make a list of dictionaries for each station.
        self.data = [{'station': station,
                      'offset': offset,
                      'te_phase': np.zeros(asize),
                      'tm_phase': np.zeros(asize),
                      're_tip': np.zeros(asize),
                      'im_tip': np.zeros(asize),
                      'te_res': np.zeros(asize),
                      'tm_res': np.zeros(asize)}
                     for station, offset in zip(self.station_list,
                                                self.station_locations)]

        # loop over mt object in edi_list and use a counter starting at 1
        # because that is what occam starts at.
        for s_index, edi in enumerate(self.edi_list):

            if self.freq_tol is None:
                station_freq = edi.Z.freq
                interp_freq = self.freq[np.where((self.freq >= station_freq.min()) &
                                                 (self.freq <= station_freq.max()))]
                # interpolate data onto given frequency list
                z_interp, t_interp = edi.interpolate(interp_freq)
                #                z_interp._compute_res_phase()
                z_interp.compute_resistivity_phase()

                rho = z_interp.resistivity
                phi = z_interp.phase
                rho_err = z_interp.resistivity_err
                if t_interp is not None:
                    tipper = t_interp.tipper
                    tipper_err = t_interp.tipper_err
                else:
                    tipper = None
                    tipper_err = None
            else:
                station_freq = edi.Z.freq
                rho = edi.Z.resistivity
                phi = edi.Z.phase
                tipper = edi.Tipper.tipper
                tipper_err = edi.Tipper.tipper_err

            self.data[s_index]['station'] = edi.station
            self.data[s_index]['offset'] = edi.offset

            for freq_num, frequency in enumerate(self.freq):
                if self.freq_tol is not None:
                    try:
                        f_index = np.where((station_freq >= frequency * (1 - self.freq_tol)) &
                                           (station_freq <= frequency * (1 + self.freq_tol)))[0][0]

                    except IndexError:
                        f_index = None
                else:
                    # skip, if the listed frequency is not available for the station
                    if (frequency in interp_freq):
                        # find the respective frequency index for the station
                        f_index = np.abs(interp_freq - frequency).argmin()
                    else:
                        f_index = None

                if f_index == None:
                    continue

                # --> get te resistivity
                self.data[s_index]['te_res'][0, freq_num] = rho[f_index, 0, 1]
                # compute error
                if rho[f_index, 0, 1] != 0.0:
                    # --> get error from data
                    if self.res_te_err is None:
                        self.data[s_index]['te_res'][1, freq_num] = \
                            np.abs(rho_err[f_index, 0, 1])
                    # --> set generic error floor
                    elif self.error_type == 'floor':
                        self.data[s_index]['te_res'][1, freq_num] = \
                            max(rho[f_index, 0, 1]*self.res_te_err/100.,
                                rho_err[f_index, 0, 1])
                    else:
                        self.data[s_index]['te_res'][1, freq_num] = \
                        rho[f_index, 0, 1]*self.res_te_err/100.

                # --> get tm resistivity
                self.data[s_index]['tm_res'][0, freq_num] = rho[f_index, 1, 0]
                # compute error
                if rho[f_index, 1, 0] != 0.0:
                    # --> get error from data
                    if self.res_tm_err is None:
                        self.data[s_index]['tm_res'][1, freq_num] = \
                            np.abs(rho_err[f_index, 1, 0])
                    # --> set generic error floor
                    elif self.error_type == 'floor':
                        self.data[s_index]['tm_res'][1, freq_num] = \
                            max(rho[f_index, 1, 0]*self.res_tm_err / 100.,
                                rho_err[f_index, 1, 0])
                    else:
                        self.data[s_index]['tm_res'][1, freq_num] = \
                        rho[f_index, 1, 0]*self.res_tm_err / 100.

                # --> get te phase
                phase_te = phi[f_index, 0, 1]
                # be sure the phase is in the first quadrant
                if phase_te > 180:
                    phase_te -= 180
                # remove any remaining phase values that are out of the first
                # quadrant
                if ((phase_te > 90) or (phase_te < 0)):
                    phase_te = 0.

                self.data[s_index]['te_phase'][0, freq_num] = phase_te
                # compute error
                # if phi[f_index, 0, 1] != 0.0:
                # --> get error from data
                phase_te_errorval = \
                np.degrees(np.arcsin(.5 *
                                     self.data[s_index]['te_res'][1, freq_num]/\
                                     self.data[s_index]['te_res'][0, freq_num]))
                if self.phase_te_err is None:
                    self.data[s_index]['te_phase'][1, freq_num] =\
                    phase_te_errorval
                # --> set generic error floor
                elif self.error_type == 'floor':
                    self.data[s_index]['te_phase'][1, freq_num] = \
                        max((self.phase_te_err / 100.) * 57. / 2.,
                            phase_te_errorval)
                else:
                    self.data[s_index]['te_phase'][1, freq_num] = \
                    (self.phase_te_err / 100.) * 57. / 2.
                # --> get tm phase and be sure its in the first quadrant
                phase_tm = phi[f_index, 1, 0] % 180            

                # remove any remaining phase values that are out of the first
                # quadrant
                if ((phase_tm > 90) or (phase_tm < 0)):
                    phase_tm = 0.
                

                self.data[s_index]['tm_phase'][0, freq_num] = phase_tm
                # compute error
                # if phi[f_index, 1, 0] != 0.0:
                # --> get error from data
                phase_tm_errorval = \
                np.degrees(np.arcsin(.5 *
                                     self.data[s_index]['tm_res'][1, freq_num]/\
                                     self.data[s_index]['tm_res'][0, freq_num]))
                if self.phase_tm_err is None:
                    self.data[s_index]['tm_phase'][1, freq_num] = \
                        phase_tm_errorval
                # --> set generic error floor
                elif self.error_type == 'floor':
                    self.data[s_index]['tm_phase'][1, freq_num] = \
                        max((self.phase_tm_err / 100.) * 57. / 2.,
                            phase_tm_errorval)
                else:
                    self.data[s_index]['tm_phase'][1, freq_num] = \
                    (self.phase_tm_err / 100.) * 57. / 2.
                # --> get Tipper
                if tipper is not None:
                    self.data[s_index]['re_tip'][0, freq_num] = \
                        tipper[f_index, 0, 1].real
                    self.data[s_index]['im_tip'][0, freq_num] = \
                        tipper[f_index, 0, 1].imag

                    # get error
                    if self.tipper_err is not None:
                        self.data[s_index]['re_tip'][1, freq_num] = \
                            self.tipper_err / 100.
                        self.data[s_index]['im_tip'][1, freq_num] = \
                            self.tipper_err / 100.
                    else:
                        self.data[s_index]['re_tip'][1, freq_num] = \
                            tipper[f_index, 0, 1].real / tipper_err[f_index, 0, 1]
                        self.data[s_index]['im_tip'][1, freq_num] = \
                            tipper[f_index, 0, 1].imag / tipper_err[f_index, 0, 1]

    def _get_data_list(self):
        """
        Get all the data needed to write a data file.
        
        """
        self.data_list = []
        for ss, sdict in enumerate(self.data, 1):
            for ff in range(self.freq.shape[0]):
                for mmode in self.mode_dict[self.model_mode]:
                    # log(te_res)
                    if mmode == 1:
                        if sdict['te_res'][0, ff] != 0.0:
                            dvalue = np.log10(sdict['te_res'][0, ff])
                            derror = sdict['te_res'][1, ff]/\
                                     (sdict['te_res'][0,ff] * np.log(10))
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

                    # te_res
                    if mmode == 9:
                        if sdict['te_res'][0, ff] != 0.0:
                            dvalue = sdict['te_res'][0, ff]
                            derror = sdict['te_res'][1, ff]
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

                    # te_phase
                    if mmode == 2:
                        if sdict['te_phase'][0, ff] != 0.0:
                            dvalue = sdict['te_phase'][0, ff]
                            derror = sdict['te_phase'][1, ff]
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

                            # log(tm_res)
                    if mmode == 5:
                        if sdict['tm_res'][0, ff] != 0.0:
                            dvalue = np.log10(sdict['tm_res'][0, ff])
                            derror = sdict['tm_res'][1, ff] /\
                                     (sdict['tm_res'][0, ff]*np.log(10))
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

                    # tm_res
                    if mmode == 10:
                        if sdict['tm_res'][0, ff] != 0.0:
                            dvalue = sdict['tm_res'][0, ff]
                            derror = sdict['tm_res'][1, ff]
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

                    # tm_phase
                    if mmode == 6:
                        if sdict['tm_phase'][0, ff] != 0.0:
                            dvalue = sdict['tm_phase'][0, ff]
                            derror = sdict['tm_phase'][1, ff]
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

                    # Re_tip
                    if mmode == 3:
                        if sdict['re_tip'][0, ff] != 0.0:
                            dvalue = sdict['re_tip'][0, ff]
                            derror = sdict['re_tip'][1, ff]
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

                    # Im_tip
                    if mmode == 4:
                        if sdict['im_tip'][0, ff] != 0.0:
                            dvalue = sdict['im_tip'][0, ff]
                            derror = sdict['im_tip'][1, ff]
                            dstr = '{0:.4f}'.format(dvalue)
                            derrstr = '{0:.4f}'.format(derror)
                            line = self._data_string.format(ss, ff + 1, mmode,
                                                            dstr, derrstr)
                            self.data_list.append(line)

    def write_data_file(self, data_fn=None):
        """
        Write a data file.
        
        Arguments:
        -----------
            **data_fn** : string
                          full path to data file. 
                          *default* is save_path/fn_basename
                          
        If there data is None, then _fill_data is called to create a profile, 
        rotate data and get all the necessary data.  This way you can use 
        write_data_file directly without going through the steps of projecting
        the stations, etc.
        
        :Example: ::
            >>> edipath = r"/home/mt/edi_files"
            >>> slst = ['mt{0:03}'.format(ss) for ss in range(1, 20)]
            >>> ocd = occam2d.Data(edi_path=edipath, station_list=slst)
            >>> ocd.save_path = r"/home/occam/line1/inv1"
            >>> ocd.write_data_file()
        
        """

        if self.data is None:
            self._fill_data()

        # get the appropriate data to write to file
        self._get_data_list()

        if data_fn is not None:
            self.data_fn = data_fn
        else:
            if self.save_path is None:
                self.save_path = os.getcwd()
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

            self.data_fn = os.path.join(self.save_path, self.fn_basename)

        data_lines = []

        # --> header line
        data_lines.append('{0:<18}{1}\n'.format('FORMAT:', self.occam_format))

        # --> title line
        if self.profile_angle is None:
            self.profile_angle = 0
        if self.geoelectric_strike is None:
            self.geoelectric_strike = 0.0
        t_str = '{0}, Profile={1:.1f} deg, Strike={2:.1f} deg'.format(
            self.title, self.profile_angle, self.geoelectric_strike)
        data_lines.append('{0:<18}{1}\n'.format('TITLE:', t_str))

        # --> sites
        data_lines.append('{0:<18}{1}\n'.format('SITES:', len(self.data)))
        for sdict in self.data:
            data_lines.append('   {0}\n'.format(sdict['station']))

        # --> offsets
        data_lines.append('{0:<18}\n'.format('OFFSETS (M):'))
        for sdict in self.data:
            data_lines.append('   {0:.1f}\n'.format(sdict['offset']))
        # --> frequencies
        data_lines.append('{0:<18}{1}\n'.format('FREQUENCIES:',
                                                self.freq.shape[0]))
        for ff in self.freq:
            data_lines.append('   {0:<10.6e}\n'.format(ff))

        # --> data
        data_lines.append('{0:<18}{1}\n'.format('DATA BLOCKS:',
                                                len(self.data_list)))
        data_lines.append(self._data_header)
        data_lines += self.data_list

        dfid = open(self.data_fn, 'w')
        dfid.writelines(data_lines)
        dfid.close()

        print('Wrote Occam2D data file to {0}'.format(self.data_fn))

    def get_profile_origin(self):
        """
        get the origin of the profile in real world coordinates
        
        Author: Alison Kirkby (2013)
        
        NEED TO ADAPT THIS TO THE CURRENT SETUP.
        """

        x, y = self.easts, self.norths
        x1, y1 = x[0], y[0]
        [m, c1] = self.profile
        x0 = (y1 + (1.0 / m) * x1 - c1) / (m + (1.0 / m))
        y0 = m * x0 + c1
        self.profile_origin = [x0, y0]

    def plot_response(self, **kwargs):
        """
        plot data and model responses as apparent resistivity, phase and
        tipper.  See PlotResponse for key words.
        
        Returns:
        ---------
            **pr_obj** : PlotResponse object 
            
        :Example: ::
            >>> pr_obj = ocd.plot_response()
        
        """

        pr_obj = PlotResponse(self.data_fn, **kwargs)

        return pr_obj

    def plot_mask_points(self, data_fn=None, marker='h', res_err_inc=.25,
                         phase_err_inc=.05, **kwargs):
        """
        An interactive plotting tool to mask points an add errorbars
        
        Arguments:
        ----------
            **res_err_inc** : float
                            amount to increase the error bars. Input as a 
                            decimal percentage.  0.3 for 30 percent
                            *Default* is 0.2 (20 percent)
                            
            **phase_err_inc** : float
                              amount to increase the error bars. Input as a 
                              decimal percentage.  0.3 for 30 percent
                              *Default* is 0.05 (5 percent)
                            
            **marker** : string
                         marker that the masked points will be
                         *Default* is 'h' for hexagon
                           
        
        :Example: ::

            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Data()
            >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>> ocd.plot_mask_points()                   
            
        """

        if data_fn is not None:
            self.data_fn = data_fn

        pr_obj = self.plot_response(**kwargs)

        # make points an attribute of self which is a data type OccamPointPicker
        self.masked_data = OccamPointPicker(pr_obj.ax_list,
                                            pr_obj.line_list,
                                            pr_obj.err_list,
                                            phase_err_inc=phase_err_inc,
                                            res_err_inc=res_err_inc)
        plt.show()

    def mask_points(self, maskpoints_obj):
        """
        mask points and rewrite the data file
        
        NEED TO REDO THIS TO FIT THE CURRENT SETUP
        """

        mp_obj = maskpoints_obj
        m_data = list(self.data)
        # rewrite the data file
        # make a reverse dictionary for locating the masked points in the data
        # file
        rploc = dict([('{0}'.format(mp_obj.fndict[key]), int(key) - 1)
                      for key in list(mp_obj.fndict.keys())])

        # make a period dictionary to locate points changed
        frpdict = dict([('{0:.5g}'.format(fr), ff)
                        for ff, fr in enumerate(1. / self.freq)])

        # loop over the data list
        for dd, dat in enumerate(mp_obj.data):
            derror = self.points.error[dd]
            # loop over the 4 main entrie
            for ss, skey in enumerate(['resxy', 'resyx', 'phasexy', 'phaseyx']):
                # rewrite any coinciding points
                for frpkey in list(frpdict.keys()):
                    try:
                        ff = frpdict[frpkey]
                        floc = self.points.fdict[dd][ss][frpkey]

                        # CHANGE APPARENT RESISTIVITY
                        if ss == 0 or ss == 1:
                            # change the apparent resistivity value
                            if m_data[rploc[str(dd)]][skey][0][ff] != \
                                    np.log10(dat[ss][floc]):
                                if dat[ss][floc] == 0:
                                    m_data[rploc[str(dd)]][skey][0][ff] = 0.0
                                else:
                                    m_data[rploc[str(dd)]][skey][0][ff] = \
                                        np.log10(dat[ss][floc])

                            # change the apparent resistivity error value
                            if dat[ss][floc] == 0.0:
                                rerr = 0.0
                            else:
                                rerr = derror[ss][floc] / dat[ss][floc] / np.log(10)
                            if m_data[rploc[str(dd)]][skey][1][ff] != rerr:
                                m_data[rploc[str(dd)]][skey][1][ff] = rerr

                        # DHANGE PHASE
                        elif ss == 2 or ss == 3:
                            # change the phase value
                            if m_data[rploc[str(dd)]][skey][0][ff] != \
                                    dat[ss][floc]:
                                if dat[ss][floc] == 0:
                                    m_data[rploc[str(dd)]][skey][0][ff] = 0.0
                                else:
                                    m_data[rploc[str(dd)]][skey][0][ff] = \
                                        dat[ss][floc]

                            # change the apparent resistivity error value
                            if dat[ss][floc] == 0.0:
                                rerr = 0.0
                            else:
                                rerr = derror[ss][floc]
                            if m_data[rploc[str(dd)]][skey][1][ff] != rerr:
                                m_data[rploc[str(dd)]][skey][1][ff] = rerr
                    except KeyError:
                        pass


class Response(object):
    """
    Reads .resp file output by Occam.  Similar structure to Data.data.
    
    If resp_fn is given in the initialization of Response, read_response_file
    is called.
    
    Arguments:
    ------------
        **resp_fn** : string
                      full path to .resp file
                      
    Attributes:
    -------------
        **resp** : is a list of dictioinaries containing the data for each
                   station.  keys include:
                   
                   * 'te_res' -- TE resisitivity in linear scale
                   * 'tm_res' -- TM resistivity in linear scale
                   * 'te_phase' -- TE phase in degrees
                   * 'tm_phase' --  TM phase in degrees in first quadrant
                   * 're_tip' -- real part of tipper along profile
                   * 'im_tip' -- imaginary part of tipper along profile
                   
               each key is a np.ndarray(2, num_freq)
               index 0 is for model response
               index 1 is for normalized misfit
               
    :Example: ::
        >>> resp_obj = occam2d.Response(r"/home/occam/line1/inv1/test_01.resp")
        
    
    
    """

    def __init__(self, resp_fn=None, **kwargs):
        self.resp_fn = resp_fn

        self.resp = None
        self.occam_dict = {'1': 'log_te_res',
                           '2': 'te_phase',
                           '3': 're_tip',
                           '4': 'im_tip',
                           '5': 'log_tm_res',
                           '6': 'tm_phase',
                           '9': 'te_res',
                           '10': 'tm_res'}

        if resp_fn is not None:
            self.read_response_file()

    def read_response_file(self, resp_fn=None):
        """
        read in response file and put into a list of dictionaries similar 
        to Data
        """

        if resp_fn is not None:
            self.resp_fn = resp_fn

        if self.resp_fn is None:
            raise OccamInputError('resp_fn is None, please input response file')

        if not os.path.isfile(self.resp_fn):
            raise OccamInputError('Could not find {0}'.format(self.resp_fn))
        try:
            r_arr = np.loadtxt(self.resp_fn, dtype=[('station', np.int),
                                                    ('freq', np.int),
                                                    ('comp', np.int),
                                                    ('z', np.int),
                                                    ('data', np.float),
                                                    ('resp', np.float),
                                                    ('err', np.float)])
        except ValueError as e:
            try:
                # for ValueError: invalid literal for long() with base 10: ...
                # which is found on some Linux environments
                # tying to recover
                r_arr = np.loadtxt(self.resp_fn, dtype=[('station', np.float),
                                                        ('freq', np.float),
                                                        ('comp', np.float),
                                                        ('z', np.float),
                                                        ('data', np.float),
                                                        ('resp', np.float),
                                                        ('err', np.float)])
                r_arr = r_arr.astype([('station', np.int),
                                      ('freq', np.int),
                                      ('comp', np.int),
                                      ('z', np.int),
                                      ('data', np.float),
                                      ('resp', np.float),
                                      ('err', np.float)])
            except:
                raise OccamInputError("Filed to read file {}. numpy error message: {}".format(self.resp_fn, e.message))

        num_stat = r_arr['station'].max()
        num_freq = r_arr['freq'].max()

        # set zero array size the first row will be the data and second the error
        asize = (2, num_freq)

        # make a list of dictionaries for each station.
        self.resp = [{'te_phase': np.zeros(asize),
                      'tm_phase': np.zeros(asize),
                      're_tip': np.zeros(asize),
                      'im_tip': np.zeros(asize),
                      'te_res': np.zeros(asize),
                      'tm_res': np.zeros(asize)}
                     for ss in range(num_stat)]

        for line in r_arr:
            # station index -1 cause python starts at 0
            ss = line['station'] - 1

            # frequency index -1 cause python starts at 0
            ff = line['freq'] - 1
            # data key
            key = self.occam_dict[str(line['comp'])]
            # put into array
            if line['comp'] == 1 or line['comp'] == 5:
                self.resp[ss][key[4:]][0, ff] = 10 ** line['resp']
                # error
                self.resp[ss][key[4:]][1, ff] = line['err'] * np.log(10)
            else:
                self.resp[ss][key][0, ff] = line['resp']
                # error
                self.resp[ss][key][1, ff] = line['err']


class Model(Startup):
    """
    Read .iter file output by Occam2d.  Builds the resistivity model from 
    mesh and regularization files found from the .iter file.  The resistivity
    model is an array(x_nodes, z_nodes) set on a regular grid, and the values 
    of the model response are filled in according to the regularization grid.
    This allows for faster plotting.  
    
    Inherets Startup because they are basically the same object.
    
    Argument:
    ----------
        **iter_fn** : string
                      full path to .iter file to read. *default* is None.
                      
        **model_fn** : string
                       full path to regularization file. *default* is None
                       and found directly from the .iter file.  Only input
                       if the regularization is different from the file that
                       is in the .iter file.
                      
        **mesh_fn** : string
                      full path to mesh file. *default* is None
                      Found directly from the model_fn file.  Only input
                      if the mesh is different from the file that
                      is in the model file.
                      
    ===================== =====================================================
    Key Words/Attributes  Description    
    ===================== =====================================================
    data_fn               full path to data file
    iter_fn               full path to .iter file
    mesh_fn               full path to mesh file
    mesh_x                np.ndarray(x_nodes, z_nodes) mesh grid for plotting
    mesh_z                np.ndarray(x_nodes, z_nodes) mesh grid for plotting
    model_values          model values from startup file
    plot_x                nodes of mesh in horizontal direction
    plot_z                nodes of mesh in vertical direction
    res_model             np.ndarray(x_nodes, z_nodes) resistivity model 
                          values in linear scale
    ===================== =====================================================
    
    
    ===================== =====================================================
    Methods               Description     
    ===================== =====================================================
    build_model           get the resistivity model from the .iter file
                          in a regular grid according to the mesh file
                          with resistivity values according to the model file
    read_iter_file        read .iter file and fill appropriate attributes
    write_iter_file       write an .iter file incase you want to set it as the
                          starting model or a priori model
    ===================== =====================================================
         
    :Example: ::
        >>> model = occam2D.Model(r"/home/occam/line1/inv1/test_01.iter")
        >>> model.build_model()
                 
    """

    def __init__(self, iter_fn=None, model_fn=None, mesh_fn=None, **kwargs):
        Startup.__init__(self, **kwargs)
        self.iter_fn = iter_fn
        self.model_fn = model_fn
        self.mesh_fn = mesh_fn
        self.data_fn = kwargs.pop('data_fn', None)
        self.model_values = kwargs.pop('model_values', None)
        self.res_model = None
        self.plot_x = None
        self.plot_z = None
        self.mesh_x = None
        self.mesh_z = None

    def read_iter_file(self, iter_fn=None):
        """
        Read an iteration file.
        
        Arguments:
        ----------
            **iter_fn** : string
                        full path to iteration file if iterpath=None.  If 
                        iterpath is input then iterfn is just the name
                        of the file without the full path.

        Returns:
        --------
        
        :Example: ::
            
            >>> import mtpy.modeling.occam2d as occam2d
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam2d.Model(itfn)
            >>> ocm.read_iter_file()
            
        """

        if iter_fn is not None:
            self.iter_fn == iter_fn

        if self.iter_fn is None:
            raise OccamInputError('iter_fn is None, input iteration file')

        # check to see if the file exists
        if os.path.exists(self.iter_fn) == False:
            raise OccamInputError('Can not find {0}'.format(self.iter_fn))

        self.save_path = os.path.dirname(self.iter_fn)

        # open file, read lines, close file
        ifid = open(self.iter_fn, 'r')
        ilines = ifid.readlines()
        ifid.close()

        ii = 0
        # put header info into dictionary with similar keys
        while ilines[ii].lower().find('param') != 0:
            iline = ilines[ii].strip().split(':')
            key = iline[0].strip().lower()
            if key.find('!') != 0:
                key = key.replace(' ', '_').replace('file', 'fn').replace('/', '_')
                value = iline[1].strip()
                try:
                    setattr(self, key, float(value))
                except ValueError:
                    setattr(self, key, value)
            ii += 1

        # get number of parameters
        iline = ilines[ii].strip().split(':')
        key = iline[0].strip().lower().replace(' ', '_')
        value = int(iline[1].strip())
        setattr(self, key, value)

        self.model_values = np.zeros(self.param_count)
        kk = int(ii + 1)

        jj = 0
        mv_index = 0
        while jj < len(ilines) - kk:
            iline = np.array(ilines[jj + kk].strip().split(), dtype='float')
            self.model_values[mv_index:mv_index + iline.shape[0]] = iline
            jj += 1
            mv_index += iline.shape[0]

        # make sure data file is full path
        if os.path.isfile(self.data_fn) == False:
            self.data_fn = os.path.join(self.save_path, self.data_fn)

        # make sure model file is full path
        if os.path.isfile(self.model_fn) == False:
            self.model_fn = os.path.join(self.save_path, self.model_fn)

    def write_iter_file(self, iter_fn=None):
        """
        write an iteration file if you need to for some reason, same as 
        startup file
        """
        if iter_fn is not None:
            self.iter_fn = iter_fn

        self.write_startup_file(iter_fn)

    def build_model(self):
        """
        build the model from the mesh, regularization grid and model file
        
        """

        # first read in the iteration file
        self.read_iter_file()

        # read in the regulariztion file
        r1 = Regularization()
        r1.read_regularization_file(self.model_fn)
        r1.model_rows = np.array(r1.model_rows)
        self.model_rows = r1.model_rows
        self.model_columns = r1.model_columns

        # read in mesh file
        r1.read_mesh_file(r1.mesh_fn)

        # get the binding offset which is the right side of the furthest left
        # block, this helps locate the model in relative space
        bndgoff = r1.binding_offset

        # make sure that the number of rows and number of columns are the same
        assert len(r1.model_rows) == len(r1.model_columns)

        # initiate the resistivity model to the shape of the FE mesh
        self.res_model = np.zeros((r1.z_nodes.shape[0], r1.x_nodes.shape[0]))

        # read in the model and set the regularization block values to map onto
        # the FE mesh so that the model can be plotted as an image or regular
        # mesh.
        mm = 0
        for ii in range(len(r1.model_rows)):
            # get the number of layers to combine
            # this index will be the first index in the vertical direction
            ny1 = r1.model_rows[:ii, 0].sum()
            # the second index  in the vertical direction
            ny2 = ny1 + r1.model_rows[ii][0]
            # make the list of amalgamated columns an array for ease
            lc = np.array(r1.model_columns[ii])
            # loop over the number of amalgamated blocks
            for jj in range(len(r1.model_columns[ii])):
                # get first in index in the horizontal direction
                nx1 = lc[:jj].sum()
                # get second index in horizontal direction
                nx2 = nx1 + lc[jj]
                # put the apporpriate resistivity value into all the amalgamated
                # model blocks of the regularization grid into the forward model
                # grid
                self.res_model[ny1:ny2, nx1:nx2] = self.model_values[mm]
                mm += 1

        # make some arrays for plotting the model
        self.plot_x = np.array([r1.x_nodes[:ii + 1].sum()
                                for ii in range(len(r1.x_nodes))])
        self.plot_z = np.array([r1.z_nodes[:ii + 1].sum()
                                for ii in range(len(r1.z_nodes))])

        # center the grid onto the station coordinates
        x0 = bndgoff - self.plot_x[r1.model_columns[0][0]]
        self.plot_x += x0

        # flip the arrays around for plotting purposes
        # plotx = plotx[::-1] and make the first layer start at zero
        self.plot_z = self.plot_z[::-1] - self.plot_z[0]

        # make a mesh grid to plot in the model coordinates
        self.mesh_x, self.mesh_z = np.meshgrid(self.plot_x, self.plot_z)

        # flip the resmodel upside down so that the top is the stations
        self.res_model = np.flipud(self.res_model)


# ==============================================================================
# plot the MT and model responses            
# ==============================================================================
class PlotResponse():
    """
    Helper class to deal with plotting the MT response and occam2d model.
    
    Arguments:
    -------------
        **data_fn** : string
                      full path to data file
                      
        **resp_fn** : string or list
                      full path(s) to response file(s)   
                    
                     
    ==================== ======================================================
    Attributes/key words            description
    ==================== ======================================================
    ax_list              list of matplotlib.axes instances for use with
                         OccamPointPicker    
    color_mode           [ 'color' | 'bw' ] plot figures in color or 
                         black and white ('bw')
    cted                 color of Data TE marker and line
    ctem                 color of Model TE marker and line
    ctewl                color of Winglink Model TE marker and line
    ctmd                 color of Data TM marker and line
    ctmm                 color of Model TM marker and line
    ctmwl                color of Winglink Model TM marker and line
    e_capsize            size of error bar caps in points
    e_capthick           line thickness of error bar caps in points
    err_list             list of line properties of error bars for use with
                         OccamPointPicker
    fig_dpi              figure resolution in dots-per-inch 
    fig_list             list of dictionaries with key words
                         station --> station name
                         fig --> matplotlib.figure instance
                         axrte --> matplotlib.axes instance for TE app.res
                         axrtm --> matplotlib.axes instance for TM app.res
                         axpte --> matplotlib.axes instance for TE phase
                         axptm --> matplotlib.axes instance for TM phase
             
    fig_num              starting number of figure
    fig_size             size of figure in inches (width, height)
    font_size            size of axes ticklabel font in points
    line_list            list of matplotlib.Line instances for use with 
                         OccamPointPicker
    lw                   line width of lines in points
    ms                   marker size in points
    mted                 marker for Data TE mode
    mtem                 marker for Model TE mode
    mtewl                marker for Winglink Model TE
    mtmd                 marker for Data TM mode
    mtmm                 marker for Model TM mode
    mtmwl                marker for Winglink TM mode
    period               np.ndarray of periods to plot 
    phase_limits         limits on phase plots in degrees (min, max)
    plot_num             [ 1 | 2 ] 
                         1 to plot both modes in a single plot
                         2 to plot modes in separate plots (default)
    plot_tipper          [ 'y' | 'n' ] plot tipper data if desired
    plot_type            [ '1' | station_list]
                         '1' --> to plot all stations in different figures
                         station_list --> to plot a few stations, give names
                         of stations ex. ['mt01', 'mt07']
    plot_yn              [ 'y' | 'n']
                         'y' --> to plot on instantiation
                         'n' --> to not plot on instantiation
    res_limits           limits on resistivity plot in log scale (min, max)
    rp_list               list of dictionaries from read2Ddata
    station_list          station_list list of stations in rp_list
    subplot_bottom       subplot spacing from bottom (relative coordinates) 
    subplot_hspace       vertical spacing between subplots
    subplot_left         subplot spacing from left  
    subplot_right        subplot spacing from right
    subplot_top          subplot spacing from top
    subplot_wspace       horizontal spacing between subplots
    wl_fn                Winglink file name (full path)
    ==================== ======================================================
    
    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots the apparent resistiviy and phase of data and
                        model if given.  called on instantiation if plot_yn
                        is 'y'.
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figures        save all the matplotlib.figure instances in fig_list
    =================== =======================================================


    :Example: ::
        >>> data_fn = r"/home/occam/line1/inv1/OccamDataFile.dat"
        >>> resp_list = [r"/home/occam/line1/inv1/test_{0:02}".format(ii) 
                         for ii in range(2, 8, 2)]
        >>> pr_obj = occam2d.PlotResponse(data_fn, resp_list, plot_tipper='y')
        
    """

    def __init__(self, data_fn, resp_fn=None, **kwargs):

        self.data_fn = data_fn
        self.resp_fn = resp_fn
        if self.resp_fn is not None:
            if type(self.resp_fn) != list:
                self.resp_fn = [self.resp_fn]

        self.wl_fn = kwargs.pop('wl_fn', None)

        self.color_mode = kwargs.pop('color_mode', 'color')

        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        self.ax_list = []
        self.line_list = []
        self.err_list = []

        # color mode
        if self.color_mode == 'color':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0, .6, .3))
            self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')

            # color for Winglink model
            self.ctewl = kwargs.pop('ctewl', (0, .6, .8))
            self.ctmwl = kwargs.pop('ctmwl', (.8, .7, 0))
            self.mtewl = kwargs.pop('mtewl', 'x')
            self.mtmwl = kwargs.pop('mtmwl', 'x')

            # color of tipper
            self.ctipr = kwargs.pop('ctipr', self.cted)
            self.ctipi = kwargs.pop('ctipi', self.ctmd)

        # black and white mode
        elif self.color_mode == 'bw':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', '*')
            self.mtmd = kwargs.pop('mtmd', 'v')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')

            # color for Winglink model
            self.ctewl = kwargs.pop('ctewl', (0.3, 0.3, 0.3))
            self.ctmwl = kwargs.pop('ctmwl', (0.3, 0.3, 0.3))
            self.mtewl = kwargs.pop('mtewl', '|')
            self.mtmwl = kwargs.pop('mtmwl', '_')

            self.ctipr = kwargs.pop('ctipr', self.cted)
            self.ctipi = kwargs.pop('ctipi', self.ctmd)

        self.phase_limits = kwargs.pop('phase_limits', (-5, 95))
        self.res_limits = kwargs.pop('res_limits', None)
        self.tip_limits = kwargs.pop('tip_limits', (-.5, .5))

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = .1
        self.subplot_hspace = .15
        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_num = kwargs.pop('plot_num', 2)
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')
        self.plot_model_error = kwargs.pop('plot_model_err', 'y')
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        if self.plot_num == 1:
            self.ylabel_coord = kwargs.pop('ylabel_coords', (-.055, .5))
        elif self.plot_num == 2:
            self.ylabel_coord = kwargs.pop('ylabel_coords', (-.12, .5))

        self.fig_list = []

        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plot the data and model response, if given, in individual plots.
         
        """

        data_obj = Data()
        data_obj.read_data_file(self.data_fn)

        rp_list = data_obj.data
        nr = len(rp_list)

        # create station list
        self.station_list = [rp['station'] for rp in rp_list]

        # boolean for adding winglink output to the plots 0 for no, 1 for yes
        addwl = 0
        # read in winglink data file
        if self.wl_fn != None:
            addwl = 1
            self.subplot_hspace + .1
            wld, wlrp_list, wlplist, wlslist, wltlist = MTwl.readOutputFile(
                self.wl_fn)
            sdict = dict([(ostation, wlistation) for wlistation in wlslist
                          for ostation in self.station_list
                          if wlistation.find(ostation) >= 0])

        # set a local parameter period for less typing
        period = data_obj.period

        # ---------------plot each respones in a different figure---------------
        if self.plot_type == '1':
            pstation_list = list(range(len(self.station_list)))

        else:
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]

            pstation_list = []
            for ii, station in enumerate(self.station_list):
                for pstation in self.plot_type:
                    if station.find(pstation) >= 0:
                        pstation_list.append(ii)

        # set the grid of subplots
        if self.plot_tipper == 'y':
            gs = gridspec.GridSpec(3, 2,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=[2, 1.5, 1])
        else:
            gs = gridspec.GridSpec(2, 2,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=[2, 1.5])

        # --> set default font size
        plt.rcParams['font.size'] = self.font_size

        # loop over each station to plot
        for ii, jj in enumerate(pstation_list):
            fig = plt.figure(self.station_list[jj],
                             self.fig_size, dpi=self.fig_dpi)
            plt.clf()

            # --> set subplot instances
            # ---plot both TE and TM in same subplot---
            if self.plot_num == 1:
                axrte = fig.add_subplot(gs[0, :])
                axrtm = axrte
                axpte = fig.add_subplot(gs[1, :], sharex=axrte)
                axptm = axpte
                if self.plot_tipper == 'y':
                    axtipre = fig.add_subplot(gs[2, :], sharex=axrte)
                    axtipim = axtipre

            # ---plot TE and TM in separate subplots---
            elif self.plot_num == 2:
                axrte = fig.add_subplot(gs[0, 0])
                axrtm = fig.add_subplot(gs[0, 1])
                axpte = fig.add_subplot(gs[1, 0], sharex=axrte)
                axptm = fig.add_subplot(gs[1, 1], sharex=axrtm)
                if self.plot_tipper == 'y':
                    axtipre = fig.add_subplot(gs[2, 0], sharex=axrte)
                    axtipim = fig.add_subplot(gs[2, 1], sharex=axrtm)

            # plot the data, it should be the same for all response files
            # empty lists for legend marker and label
            rlistte = []
            llistte = []
            rlisttm = []
            llisttm = []
            # ------------Plot Resistivity----------------------------------
            # cut out missing data points first
            # --> data
            rxy = np.where(rp_list[jj]['te_res'][0] != 0)[0]
            ryx = np.where(rp_list[jj]['tm_res'][0] != 0)[0]

            # --> TE mode Data
            if len(rxy) > 0:
                rte_err = rp_list[jj]['te_res'][1, rxy] * \
                          rp_list[jj]['te_res'][0, rxy]
                rte = plot_errorbar(axrte,
                                    period[rxy],
                                    rp_list[jj]['te_res'][0, rxy],
                                    ls=':',
                                    marker=self.mted,
                                    ms=self.ms,
                                    color=self.cted,
                                    y_error=rte_err,
                                    lw=self.lw,
                                    e_capsize=self.e_capsize,
                                    e_capthick=self.e_capthick)

                rlistte.append(rte[0])
                llistte.append('$Obs_{TE}$')
            else:
                rte = [None, [None, None, None], [None, None, None]]

                # --> TM mode data
            if len(ryx) > 0:
                rtm_err = rp_list[jj]['tm_res'][1, ryx] * \
                          rp_list[jj]['tm_res'][0, ryx]
                rtm = plot_errorbar(axrtm,
                                    period[ryx],
                                    rp_list[jj]['tm_res'][0, ryx],
                                    ls=':',
                                    marker=self.mtmd,
                                    ms=self.ms,
                                    color=self.ctmd,
                                    y_error=rtm_err,
                                    lw=self.lw,
                                    e_capsize=self.e_capsize,
                                    e_capthick=self.e_capthick)

                rlisttm.append(rtm[0])
                llisttm.append('$Obs_{TM}$')
            else:
                rtm = [None, [None, None, None], [None, None, None]]
            # --------------------plot phase--------------------------------
            # cut out missing data points first
            # --> data
            pxy = np.where(rp_list[jj]['te_phase'][0] != 0)[0]
            pyx = np.where(rp_list[jj]['tm_phase'][0] != 0)[0]

            # --> TE mode data
            if len(pxy) > 0:
                pte = plot_errorbar(axpte,
                                    period[pxy],
                                    rp_list[jj]['te_phase'][0, pxy],
                                    ls=':',
                                    marker=self.mted,
                                    ms=self.ms,
                                    color=self.cted,
                                    y_error=rp_list[jj]['te_phase'][1, pxy],
                                    lw=self.lw,
                                    e_capsize=self.e_capsize,
                                    e_capthick=self.e_capthick)

            else:
                pte = [None, [None, None, None], [None, None, None]]

            # --> TM mode data
            if len(pyx) > 0:
                ptm = plot_errorbar(axptm,
                                    period[pyx],
                                    rp_list[jj]['tm_phase'][0, pyx],
                                    ls=':',
                                    marker=self.mtmd,
                                    ms=self.ms,
                                    color=self.ctmd,
                                    y_error=rp_list[jj]['tm_phase'][1, pyx],
                                    lw=self.lw,
                                    e_capsize=self.e_capsize,
                                    e_capthick=self.e_capthick)

            else:
                ptm = [None, [None, None, None], [None, None, None]]

            # append axis properties to lists that can be used by
            # OccamPointPicker
            self.ax_list.append([axrte, axrtm, axpte, axptm])
            self.line_list.append([rte[0], rtm[0], pte[0], ptm[0]])
            self.err_list.append([[rte[1][0], rte[1][1], rte[2][0]],
                                  [rtm[1][0], rtm[1][1], rtm[2][0]],
                                  [pte[1][0], pte[1][1], pte[2][0]],
                                  [ptm[1][0], ptm[1][1], ptm[2][0]]])

            # ---------------------plot tipper----------------------------------
            if self.plot_tipper == 'y':
                t_list = []
                t_label = []
                txy = np.where(rp_list[jj]['re_tip'][0] != 0)[0]
                tyx = np.where(rp_list[jj]['im_tip'][0] != 0)[0]
                # --> real tipper  data
                if len(txy) > 0:
                    per_list_p = []
                    tpr_list_p = []
                    per_list_n = []
                    tpr_list_n = []
                    for per, tpr in zip(period[txy],
                                        rp_list[jj]['re_tip'][0, txy]):
                        if tpr >= 0:
                            per_list_p.append(per)
                            tpr_list_p.append(tpr)
                        else:
                            per_list_n.append(per)
                            tpr_list_n.append(tpr)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipre.stem(per_list_p,
                                                              tpr_list_p,
                                                              markerfmt='^',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctipr)
                        plt.setp(m_line, 'markeredgecolor', self.ctipr)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctipr)
                        plt.setp(b_line, 'linewidth', .01)
                        t_list.append(m_line)
                        t_label.append('Real')
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipre.stem(per_list_n,
                                                              tpr_list_n,
                                                              markerfmt='v',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctipr)
                        plt.setp(m_line, 'markeredgecolor', self.ctipr)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctipr)
                        plt.setp(b_line, 'linewidth', .01)
                        if len(t_list) == 0:
                            t_list.append(m_line)
                            t_label.append('Real')

                else:
                    pass
                if len(tyx) > 0:
                    per_list_p = []
                    tpi_list_p = []
                    per_list_n = []
                    tpi_list_n = []
                    for per, tpi in zip(period[tyx],
                                        rp_list[jj]['im_tip'][0, tyx]):
                        if tpi >= 0:
                            per_list_p.append(per)
                            tpi_list_p.append(tpi)
                        else:
                            per_list_n.append(per)
                            tpi_list_n.append(tpi)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipim.stem(per_list_p,
                                                              tpi_list_p,
                                                              markerfmt='^',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctipi)
                        plt.setp(m_line, 'markeredgecolor', self.ctipi)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctipi)
                        plt.setp(b_line, 'linewidth', .01)
                        t_list.append(m_line)
                        t_label.append('Imag')
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipim.stem(per_list_n,
                                                              tpi_list_n,
                                                              markerfmt='v',
                                                              basefmt='k')
                        plt.setp(m_line, 'markerfacecolor', self.ctipi)
                        plt.setp(m_line, 'markeredgecolor', self.ctipi)
                        plt.setp(m_line, 'markersize', self.ms)
                        plt.setp(s_line, 'linewidth', self.lw)
                        plt.setp(s_line, 'color', self.ctipi)
                        plt.setp(b_line, 'linewidth', .01)
                        if len(t_list) <= 1:
                            t_list.append(m_line)
                            t_label.append('Imag')

                else:
                    pass

            # ------------------- plot model response --------------------------
            if self.resp_fn is not None:
                num_resp = len(self.resp_fn)
                for rr, rfn in enumerate(self.resp_fn):
                    resp_obj = Response()
                    resp_obj.read_response_file(rfn)

                    rp = resp_obj.resp
                    # create colors for different responses
                    if self.color_mode == 'color':
                        cxy = (0,
                               .4 + float(rr) / (3 * num_resp),
                               0)
                        cyx = (.7 + float(rr) / (4 * num_resp),
                               .13,
                               .63 - float(rr) / (4 * num_resp))
                    elif self.color_mode == 'bw':
                        cxy = (1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.))
                        cyx = (1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.))

                    # calculate rms's
                    rmslistte = np.hstack((rp[jj]['te_res'][1],
                                           rp[jj]['te_phase'][1]))
                    rmslisttm = np.hstack((rp[jj]['tm_res'][1],
                                           rp[jj]['tm_phase'][1]))
                    rmste = np.sqrt(np.sum([rms ** 2 for rms in rmslistte]) /
                                    len(rmslistte))
                    rmstm = np.sqrt(np.sum([rms ** 2 for rms in rmslisttm]) /
                                    len(rmslisttm))

                    # ------------Plot Resistivity------------------------------
                    # cut out missing data points first
                    # --> response
                    mrxy = np.where(rp[jj]['te_res'][0] != 0)[0]
                    mryx = np.where(rp[jj]['tm_res'][0] != 0)[0]

                    # --> TE mode Model Response
                    if len(mrxy) > 0:
                        r3 = plot_errorbar(axrte,
                                           period[mrxy],
                                           rp[jj]['te_res'][0, mrxy],
                                           ls='--',
                                           marker=self.mtem,
                                           ms=self.ms,
                                           color=cxy,
                                           y_error=None,
                                           lw=self.lw,
                                           e_capsize=self.e_capsize,
                                           e_capthick=self.e_capthick)

                        rlistte.append(r3[0])
                        llistte.append('$Mod_{TE}$ ' + '{0:.2f}'.format(rmste))
                    else:
                        pass

                    # --> TM mode model response
                    if len(mryx) > 0:
                        r4 = plot_errorbar(axrtm,
                                           period[mryx],
                                           rp[jj]['tm_res'][0, mryx],
                                           ls='--',
                                           marker=self.mtmm,
                                           ms=self.ms,
                                           color=cyx,
                                           y_error=None,
                                           lw=self.lw,
                                           e_capsize=self.e_capsize,
                                           e_capthick=self.e_capthick)
                        rlisttm.append(r4[0])
                        llisttm.append('$Mod_{TM}$ ' + '{0:.2f}'.format(rmstm))
                    else:
                        pass

                    # --------------------plot phase--------------------------------
                    # cut out missing data points first
                    # --> reponse
                    mpxy = np.where(rp[jj]['te_phase'][0] != 0)[0]
                    mpyx = np.where(rp[jj]['tm_phase'][0] != 0)[0]

                    # --> TE mode response
                    if len(mpxy) > 0:
                        p3 = plot_errorbar(axpte,
                                           period[mpxy],
                                           rp[jj]['te_phase'][0, mpxy],
                                           ls='--',
                                           ms=self.ms,
                                           color=cxy,
                                           y_error=None,
                                           lw=self.lw,
                                           e_capsize=self.e_capsize,
                                           e_capthick=self.e_capthick)

                    else:
                        pass

                    # --> TM mode response
                    if len(mpyx) > 0:
                        p4 = plot_errorbar(axptm,
                                           period[mpyx],
                                           rp[jj]['tm_phase'][0, mpyx],
                                           ls='--',
                                           marker=self.mtmm,
                                           ms=self.ms,
                                           color=cyx,
                                           y_error=None,
                                           lw=self.lw,
                                           e_capsize=self.e_capsize,
                                           e_capthick=self.e_capthick)
                    else:
                        pass

                    # ---------------------plot tipper--------------------------
                    if self.plot_tipper == 'y':
                        txy = np.where(rp[jj]['re_tip'][0] != 0)[0]
                        tyx = np.where(rp[jj]['im_tip'][0] != 0)[0]
                        # --> real tipper  data
                        if len(txy) > 0:
                            per_list_p = []
                            tpr_list_p = []
                            per_list_n = []
                            tpr_list_n = []
                            for per, tpr in zip(period[txy],
                                                rp[jj]['re_tip'][0, txy]):
                                if tpr >= 0:
                                    per_list_p.append(per)
                                    tpr_list_p.append(tpr)
                                else:
                                    per_list_n.append(per)
                                    tpr_list_n.append(tpr)
                            if len(per_list_p) > 0:
                                m_line, s_line, b_line = axtipre.stem(per_list_p,
                                                                      tpr_list_p,
                                                                      markerfmt='^',
                                                                      basefmt='k')
                                plt.setp(m_line, 'markerfacecolor', cxy)
                                plt.setp(m_line, 'markeredgecolor', cxy)
                                plt.setp(m_line, 'markersize', self.ms)
                                plt.setp(s_line, 'linewidth', self.lw)
                                plt.setp(s_line, 'color', cxy)
                                plt.setp(b_line, 'linewidth', .01)
                            if len(per_list_n) > 0:
                                m_line, s_line, b_line = axtipre.stem(per_list_n,
                                                                      tpr_list_n,
                                                                      markerfmt='v',
                                                                      basefmt='k')
                                plt.setp(m_line, 'markerfacecolor', cxy)
                                plt.setp(m_line, 'markeredgecolor', cxy)
                                plt.setp(m_line, 'markersize', self.ms)
                                plt.setp(s_line, 'linewidth', self.lw)
                                plt.setp(s_line, 'color', cxy)
                                plt.setp(b_line, 'linewidth', .01)

                        else:
                            pass
                        if len(tyx) > 0:
                            per_list_p = []
                            tpi_list_p = []
                            per_list_n = []
                            tpi_list_n = []
                            for per, tpi in zip(period[tyx],
                                                rp[jj]['im_tip'][0, tyx]):
                                if tpi >= 0:
                                    per_list_p.append(per)
                                    tpi_list_p.append(tpi)
                                else:
                                    per_list_n.append(per)
                                    tpi_list_n.append(tpi)
                            if len(per_list_p) > 0:
                                m_line, s_line, b_line = axtipim.stem(per_list_p,
                                                                      tpi_list_p,
                                                                      markerfmt='^',
                                                                      basefmt='k')
                                plt.setp(m_line, 'markerfacecolor', cyx)
                                plt.setp(m_line, 'markeredgecolor', cyx)
                                plt.setp(m_line, 'markersize', self.ms)
                                plt.setp(s_line, 'linewidth', self.lw)
                                plt.setp(s_line, 'color', cyx)
                                plt.setp(b_line, 'linewidth', .01)
                            if len(per_list_n) > 0:
                                m_line, s_line, b_line = axtipim.stem(per_list_n,
                                                                      tpi_list_n,
                                                                      markerfmt='v',
                                                                      basefmt='k')
                                plt.setp(m_line, 'markerfacecolor', cyx)
                                plt.setp(m_line, 'markeredgecolor', cyx)
                                plt.setp(m_line, 'markersize', self.ms)
                                plt.setp(s_line, 'linewidth', self.lw)
                                plt.setp(s_line, 'color', cyx)
                                plt.setp(b_line, 'linewidth', .01)

                        else:
                            pass
            # --------------add in winglink responses------------------------
            if addwl == 1:
                try:
                    wlrms = wld[sdict[self.station_list[jj]]]['rms']
                    axrte.set_title(self.station_list[jj] +
                                    '\n rms_occ_TE={0:.2f}'.format(rmste) +
                                    'rms_occ_TM={0:.2f}'.format(rmstm) +
                                    'rms_wl={0:.2f}'.format(wlrms),
                                    fontdict={'size': self.font_size,
                                              'weight': 'bold'})
                    for ww, wlistation in enumerate(wlslist):
                        if wlistation.find(self.station_list[jj]) == 0:
                            print('{0} was Found {0} in winglink file'.format(
                                self.station_list[jj], wlistation))
                            wlrpdict = wlrp_list[ww]

                    zrxy = [np.where(wlrpdict['te_res'][0] != 0)[0]]
                    zryx = [np.where(wlrpdict['tm_res'][0] != 0)[0]]

                    # plot winglink resistivity
                    r5 = axrte.loglog(wlplist[zrxy],
                                      wlrpdict['te_res'][1][zrxy],
                                      ls='-.',
                                      marker=self.mtewl,
                                      ms=self.ms,
                                      color=self.ctewl,
                                      mfc=self.ctewl,
                                      lw=self.lw)
                    r6 = axrtm.loglog(wlplist[zryx],
                                      wlrpdict['tm_res'][1][zryx],
                                      ls='-.',
                                      marker=self.mtmwl,
                                      ms=self.ms,
                                      color=self.ctmwl,
                                      mfc=self.ctmwl,
                                      lw=self.lw)

                    # plot winglink phase
                    axpte.semilogx(wlplist[zrxy],
                                   wlrpdict['te_phase'][1][zrxy],
                                   ls='-.',
                                   marker=self.mtewl,
                                   ms=self.ms,
                                   color=self.ctewl,
                                   mfc=self.ctewl,
                                   lw=self.lw)

                    axptm.semilogx(wlplist[zryx],
                                   wlrpdict['tm_phase'][1][zryx],
                                   ls='-.',
                                   marker=self.mtmwl,
                                   ms=self.ms,
                                   color=self.ctmwl,
                                   mfc=self.ctmwl,
                                   lw=self.lw)

                    rlistte.append(r5[0])
                    rlisttm.append(r6[0])
                    llistte.append('$WLMod_{TE}$ ' + '{0:.2f}'.format(wlrms))
                    llisttm.append('$WLMod_{TM}$ ' + '{0:.2f}'.format(wlrms))
                except (IndexError, KeyError):
                    print('Station not present')
            else:
                if self.plot_num == 1:
                    axrte.set_title(self.station_list[jj],
                                    fontdict={'size': self.font_size + 2,
                                              'weight': 'bold'})
                elif self.plot_num == 2:
                    fig.suptitle(self.station_list[jj],
                                 fontdict={'size': self.font_size + 2,
                                           'weight': 'bold'})

            # set the axis properties
            ax_list = [axrte, axrtm]
            for aa, axr in enumerate(ax_list):
                # set both axes to logarithmic scale
                axr.set_xscale('log', nonposx='clip')

                try:
                    axr.set_yscale('log', nonposy='clip')
                except ValueError:
                    pass

                # put on a grid
                axr.grid(True, alpha=.3, which='both', lw=.5 * self.lw)
                axr.yaxis.set_label_coords(self.ylabel_coord[0],
                                           self.ylabel_coord[1])

                # set resistivity limits if desired
                if self.res_limits != None:
                    axr.set_ylim(10 ** self.res_limits[0],
                                 10 ** self.res_limits[1])

                # set the tick labels to invisible
                plt.setp(axr.xaxis.get_ticklabels(), visible=False)
                if aa == 0:
                    axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

                # set legend based on the plot type
                if self.plot_num == 1:
                    if aa == 0:
                        axr.legend(rlistte + rlisttm, llistte + llisttm,
                                   loc=2, markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size': self.font_size + 1})
                elif self.plot_num == 2:
                    if aa == 0:
                        axr.legend(rlistte,
                                   llistte,
                                   loc=2, markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size': self.font_size + 1})

                    if aa == 1:
                        axr.legend(rlisttm,
                                   llisttm,
                                   loc=2, markerscale=1,
                                   borderaxespad=.05,
                                   labelspacing=.08,
                                   handletextpad=.15,
                                   borderpad=.05,
                                   prop={'size': self.font_size + 1})

            # set Properties for the phase axes
            for aa, axp in enumerate([axpte, axptm]):
                # set the x-axis to log scale
                axp.set_xscale('log', nonposx='clip')

                # set the phase limits
                axp.set_ylim(self.phase_limits)

                # put a grid on the subplot
                axp.grid(True, alpha=.3, which='both', lw=.5 * self.lw)

                # set the tick locations
                axp.yaxis.set_major_locator(MultipleLocator(10))
                axp.yaxis.set_minor_locator(MultipleLocator(2))

                # set the x axis label
                if self.plot_tipper == 'y':
                    plt.setp(axp.get_xticklabels(), visible=False)
                else:
                    axp.set_xlabel('Period (s)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

                # put the y label on the far left plot
                axp.yaxis.set_label_coords(self.ylabel_coord[0],
                                           self.ylabel_coord[1])
                if aa == 0:
                    axp.set_ylabel('Phase (deg)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

            # set axes properties of tipper axis
            if self.plot_tipper == 'y':
                for aa, axt in enumerate([axtipre, axtipim]):
                    axt.set_xscale('log', nonposx='clip')

                    # set tipper limits
                    axt.set_ylim(self.tip_limits)

                    # put a grid on the subplot
                    axt.grid(True, alpha=.3, which='both', lw=.5 * self.lw)

                    # set the tick locations
                    axt.yaxis.set_major_locator(MultipleLocator(.2))
                    axt.yaxis.set_minor_locator(MultipleLocator(.1))

                    # set the x axis label
                    axt.set_xlabel('Period (s)',
                                   fontdict={'size': self.font_size + 2,
                                             'weight': 'bold'})

                    axt.set_xlim(10 ** np.floor(np.log10(data_obj.period.min())),
                                 10 ** np.ceil(np.log10(data_obj.period.max())))

                    # put the y label on the far left plot
                    axt.yaxis.set_label_coords(self.ylabel_coord[0],
                                               self.ylabel_coord[1])
                    if aa == 0:
                        axt.set_ylabel('Tipper',
                                       fontdict={'size': self.font_size + 2,
                                                 'weight': 'bold'})
                        if self.plot_num == 2:
                            axt.text(axt.get_xlim()[0] * 1.25,
                                     self.tip_limits[1] * .9,
                                     'Real', horizontalalignment='left',
                                     verticalalignment='top',
                                     bbox={'facecolor': 'white'},
                                     fontdict={'size': self.font_size + 1})
                        else:
                            axt.legend(t_list, t_label,
                                       loc=2, markerscale=1,
                                       borderaxespad=.05,
                                       labelspacing=.08,
                                       handletextpad=.15,
                                       borderpad=.05,
                                       prop={'size': self.font_size + 1})
                    if aa == 1:
                        if self.plot_num == 2:
                            axt.text(axt.get_xlim()[0] * 1.25,
                                     self.tip_limits[1] * .9,
                                     'Imag', horizontalalignment='left',
                                     verticalalignment='top',
                                     bbox={'facecolor': 'white'},
                                     fontdict={'size': self.font_size + 1})

            # make sure the axis and figure are accessible to the user
            self.fig_list.append({'station': self.station_list[jj],
                                  'fig': fig, 'axrte': axrte, 'axrtm': axrtm,
                                  'axpte': axpte, 'axptm': axptm})

        # set the plot to be full screen well at least try
        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plot2DResponses()
            >>> #change color of te markers to a gray-blue
            >>> p1.cted = (.5, .5, .7)
            >>> p1.redraw_plot()
        """

        plt.close('all')
        self.plot()

    def save_figures(self, save_path, fig_fmt='pdf', fig_dpi=None,
                     close_fig='y'):
        """
        save all the figure that are in self.fig_list
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plot2DResponses()
            >>> p1.save_figures(r"/home/occam2d/Figures", fig_fmt='jpg')
        """

        if not os.path.exists(save_path):
            os.mkdir(save_path)

        for fdict in self.fig_list:
            svfn = '{0}_resp.{1}'.format(fdict['station'], fig_fmt)
            fdict['fig'].savefig(os.path.join(save_path, svfn),
                                 dpi=self.fig_dpi)
            if close_fig == 'y':
                plt.close(fdict['fig'])

            print("saved figure to {0}".format(os.path.join(save_path, svfn)))


# ==============================================================================
# plot model
# ==============================================================================
class PlotModel(Model):
    """
    plot the 2D model found by Occam2D.  The model is displayed as a meshgrid
    instead of model bricks.  This speeds things up considerably.  
    
    Inherets the Model class to take advantage of the attributes and methods
    already coded.
    
    Arguments:
    -----------
        **iter_fn** : string
                      full path to iteration file.  From here all the 
                      necessary files can be found assuming they are in the 
                      same directory.  If they are not then need to input
                      manually.
    
    
    ======================= ===============================================
    keywords                description
    ======================= ===============================================
    block_font_size         font size of block number is blocknum == 'on'
    blocknum                [ 'on' | 'off' ] to plot regulariztion block 
                            numbers.
    cb_pad                  padding between axes edge and color bar 
    cb_shrink               percentage to shrink the color bar
    climits                 limits of the color scale for resistivity
                            in log scale (min, max)
    cmap                    name of color map for resistivity values
    femesh                  plot the finite element mesh
    femesh_triangles        plot the finite element mesh with each block
                            divided into four triangles
    fig_aspect              aspect ratio between width and height of 
                            resistivity image. 1 for equal axes
    fig_dpi                 resolution of figure in dots-per-inch
    fig_num                 number of figure instance
    fig_size                size of figure in inches (width, height)
    font_size               size of axes tick labels, axes labels is +2
    grid                    [ 'both' | 'major' |'minor' | None ] string 
                            to tell the program to make a grid on the 
                            specified axes.
    meshnum                 [ 'on' | 'off' ] 'on' will plot finite element
                            mesh numbers
    meshnum_font_size       font size of mesh numbers if meshnum == 'on'
    ms                      size of station marker 
    plot_yn                 [ 'y' | 'n']
                            'y' --> to plot on instantiation
                            'n' --> to not plot on instantiation
    regmesh                 [ 'on' | 'off' ] plot the regularization mesh
                            plots as blue lines
    station_color           color of station marker
    station_font_color      color station label
    station_font_pad        padding between station label and marker
    station_font_rotation   angle of station label in degrees 0 is 
                            horizontal
    station_font_size       font size of station label
    station_font_weight     font weight of station label
    station_id              index to take station label from station name
    station_marker          station marker.  if inputing a LaTex marker
                            be sure to input as r"LaTexMarker" otherwise
                            might not plot properly
    subplot_bottom          subplot spacing from bottom  
    subplot_left            subplot spacing from left  
    subplot_right           subplot spacing from right
    subplot_top             subplot spacing from top
    title                   title of plot.  If None then the name of the
                            iteration file and containing folder will be
                            the title with RMS and Roughness.
    xlimits                 limits of plot in x-direction in (km) 
    xminorticks             increment of minor ticks in x direction
    xpad                    padding in x-direction in km
    ylimits                 depth limits of plot positive down (km)
    yminorticks             increment of minor ticks in y-direction
    ypad                    padding in negative y-direction (km)
    yscale                  [ 'km' | 'm' ] scale of plot, if 'm' everything
                            will be scaled accordingly.
    ======================= ===============================================
    
    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots resistivity model.  
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired 
                        location and format
    =================== ======================================================
    
    :Example: 
    ---------------
        >>> import mtpy.modeling.occam2d as occam2d
        >>> model_plot = occam2d.PlotModel(r"/home/occam/Inv1/mt_01.iter")
        >>> # change the color limits
        >>> model_plot.climits = (1, 4)
        >>> model_plot.redraw_plot()
        >>> #change len of station name
        >>> model_plot.station_id = [2, 5]
        >>> model_plot.redraw_plot()
        
    
    """

    def __init__(self, iter_fn=None, data_fn=None, **kwargs):
        Model.__init__(self, iter_fn, **kwargs)

        self.yscale = kwargs.pop('yscale', 'km')

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')

        self.xpad = kwargs.pop('xpad', 1.0)
        self.ypad = kwargs.pop('ypad', 1.0)

        self.ms = kwargs.pop('ms', 10)

        self.station_locations = None
        self.station_list = None
        self.station_id = kwargs.pop('station_id', None)
        self.station_font_size = kwargs.pop('station_font_size', 8)
        self.station_font_pad = kwargs.pop('station_font_pad', 1.0)
        self.station_font_weight = kwargs.pop('station_font_weight', 'bold')
        self.station_font_rotation = kwargs.pop('station_font_rotation', 60)
        self.station_font_color = kwargs.pop('station_font_color', 'k')
        self.station_marker = kwargs.pop('station_marker',
                                         r"$\blacktriangledown$")
        self.station_color = kwargs.pop('station_color', 'k')

        self.ylimits = kwargs.pop('ylimits', None)
        self.xlimits = kwargs.pop('xlimits', None)

        self.xminorticks = kwargs.pop('xminorticks', 5)
        self.yminorticks = kwargs.pop('yminorticks', 1)

        self.climits = kwargs.pop('climits', (0, 4))
        self.cmap = kwargs.pop('cmap', 'jet_r')
        if type(self.cmap) == str:
            self.cmap = cm.get_cmap(self.cmap)
        self.font_size = kwargs.pop('font_size', 8)

        self.femesh = kwargs.pop('femesh', 'off')
        self.femesh_triangles = kwargs.pop('femesh_triangles', 'off')
        self.femesh_lw = kwargs.pop('femesh_lw', .4)
        self.femesh_color = kwargs.pop('femesh_color', 'k')
        self.meshnum = kwargs.pop('meshnum', 'off')
        self.meshnum_font_size = kwargs.pop('meshnum_font_size', 3)

        self.regmesh = kwargs.pop('regmesh', 'off')
        self.regmesh_lw = kwargs.pop('regmesh_lw', .4)
        self.regmesh_color = kwargs.pop('regmesh_color', 'b')
        self.blocknum = kwargs.pop('blocknum', 'off')
        self.block_font_size = kwargs.pop('block_font_size', 3)
        self.grid = kwargs.pop('grid', None)

        self.cb_shrink = kwargs.pop('cb_shrink', .8)
        self.cb_pad = kwargs.pop('cb_pad', .01)

        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plotModel will plot the model output by occam2d in the iteration file.
        
        
        :Example: ::
            
            >>> import mtpy.modeling.occam2d as occam2d
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> model_plot = occam2d.PlotModel(itfn)
            >>> model_plot.ms = 20
            >>> model_plot.ylimits = (0,.350)
            >>> model_plot.yscale = 'm'
            >>> model_plot.spad = .10
            >>> model_plot.ypad = .125
            >>> model_plot.xpad = .025
            >>> model_plot.climits = (0,2.5)
            >>> model_plot.aspect = 'equal'
            >>> model_plot.redraw_plot()
            
        """
        # --> read in iteration file and build the model
        self.read_iter_file()
        self.build_model()

        # --> get station locations and names from data file
        d_object = Data()
        d_object.read_data_file(self.data_fn)
        setattr(self, 'station_locations', d_object.station_locations.copy())
        setattr(self, 'station_list', d_object.station_list.copy())

        # set the scale of the plot
        if self.yscale == 'km':
            df = 1000.
            pf = 1.0
        elif self.yscale == 'm':
            df = 1.
            pf = 1000.
        else:
            df = 1000.
            pf = 1.0

        # set some figure properties to use the maiximum space
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top

        # station font dictionary
        fdict = {'size': self.station_font_size,
                 'weight': self.station_font_weight,
                 'rotation': self.station_font_rotation,
                 'color': self.station_font_color}

        # plot the model as a mesh
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        # add a subplot to the figure with the specified aspect ratio
        ax = self.fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)

        # plot the model as a pcolormesh so the extents are constrained to
        # the model coordinates
        ax.pcolormesh(self.mesh_x / df,
                      self.mesh_z / df,
                      self.res_model,
                      cmap=self.cmap,
                      vmin=self.climits[0],
                      vmax=self.climits[1])

        # make a colorbar for the resistivity
        cbx = mcb.make_axes(ax, shrink=self.cb_shrink, pad=self.cb_pad)
        cb = mcb.ColorbarBase(cbx[0],
                              cmap=self.cmap,
                              norm=Normalize(vmin=self.climits[0],
                                             vmax=self.climits[1]))

        cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size': self.font_size + 1, 'weight': 'bold'})
        cb.set_ticks(np.arange(int(self.climits[0]), int(self.climits[1]) + 1))
        cb.set_ticklabels(['10$^{0}$'.format('{' + str(nn) + '}') for nn in
                           np.arange(int(self.climits[0]),
                                     int(self.climits[1]) + 1)])

        # set the offsets of the stations and plot the stations
        # need to figure out a way to set the marker at the surface in all
        # views.
        for offset, name in zip(self.station_locations, self.station_list):
            # plot the station marker
            # plots a V for the station cause when you use scatter the spacing
            # is variable if you change the limits of the y axis, this way it
            # always plots at the surface.
            ax.text(offset / df,
                    self.plot_z.min(),
                    self.station_marker,
                    horizontalalignment='center',
                    verticalalignment='baseline',
                    fontdict={'size': self.ms, 'color': self.station_color})

            # put station id onto station marker
            # if there is a station id index
            if self.station_id != None:
                ax.text(offset / df,
                        -self.station_font_pad * pf,
                        name[self.station_id[0]:self.station_id[1]],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=fdict)
            # otherwise put on the full station name found form data file
            else:
                ax.text(offset / df,
                        -self.station_font_pad * pf,
                        name,
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=fdict)

        # set the initial limits of the plot to be square about the profile line
        if self.ylimits == None:
            ax.set_ylim(abs(self.station_locations.max() -
                            self.station_locations.min()) / df,
                        -self.ypad * pf)
        else:
            ax.set_ylim(self.ylimits[1] * pf,
                        (self.ylimits[0] - self.ypad) * pf)
        if self.xlimits == None:
            ax.set_xlim(self.station_locations.min() / df - (self.xpad * pf),
                        self.station_locations.max() / df + (self.xpad * pf))
        else:
            ax.set_xlim(self.xlimits[0] * pf, self.xlimits[1] * pf)

        # set the axis properties
        ax.xaxis.set_minor_locator(MultipleLocator(self.xminorticks * pf))
        ax.yaxis.set_minor_locator(MultipleLocator(self.yminorticks * pf))

        # set axes labels
        ax.set_xlabel('Horizontal Distance ({0})'.format(self.yscale),
                      fontdict={'size': self.font_size + 2, 'weight': 'bold'})
        ax.set_ylabel('Depth ({0})'.format(self.yscale),
                      fontdict={'size': self.font_size + 2, 'weight': 'bold'})

        # put a grid on if one is desired
        if self.grid is not None:
            ax.grid(alpha=.3, which=self.grid, lw=.35)

        # set title as rms and roughness
        if type(self.title) is str:
            if self.title == 'on':
                titlestr = os.path.join(os.path.basename(
                    os.path.dirname(self.iter_fn)),
                    os.path.basename(self.iter_fn))
                ax.set_title('{0}: RMS={1:.2f}, Roughness={2:.0f}'.format(
                    titlestr, self.misfit_value, self.roughness_value),
                    fontdict={'size': self.font_size + 1,
                              'weight': 'bold'})
            else:
                ax.set_title('{0}; RMS={1:.2f}, Roughness={2:.0f}'.format(
                    self.title, self.misfit_value,
                    self.roughness_value),
                    fontdict={'size': self.font_size + 1,
                              'weight': 'bold'})
        else:
            print('RMS {0:.2f}, Roughness={1:.0f}'.format(self.misfit_value,
                                                          self.roughness_value))

        # plot forward model mesh
        # making an extended list seperated by None's speeds up the plotting
        # by as much as 99 percent, handy
        if self.femesh == 'on':
            row_line_xlist = []
            row_line_ylist = []
            for xx in self.plot_x / df:
                row_line_xlist.extend([xx, xx])
                row_line_xlist.append(None)
                row_line_ylist.extend([0, self.plot_zy[0] / df])
                row_line_ylist.append(None)

            # plot column lines (variables are a little bit of a misnomer)
            ax.plot(row_line_xlist,
                    row_line_ylist,
                    color='k',
                    lw=.5)

            col_line_xlist = []
            col_line_ylist = []
            for yy in self.plot_z / df:
                col_line_xlist.extend([self.plot_x[0] / df,
                                       self.plot_x[-1] / df])
                col_line_xlist.append(None)
                col_line_ylist.extend([yy, yy])
                col_line_ylist.append(None)

            # plot row lines (variables are a little bit of a misnomer)
            ax.plot(col_line_xlist,
                    col_line_ylist,
                    color='k',
                    lw=.5)

        if self.femesh_triangles == 'on':
            row_line_xlist = []
            row_line_ylist = []
            for xx in self.plot_x / df:
                row_line_xlist.extend([xx, xx])
                row_line_xlist.append(None)
                row_line_ylist.extend([0, self.plot_z[0] / df])
                row_line_ylist.append(None)

            # plot columns
            ax.plot(row_line_xlist,
                    row_line_ylist,
                    color='k',
                    lw=.5)

            col_line_xlist = []
            col_line_ylist = []
            for yy in self.plot_z / df:
                col_line_xlist.extend([self.plot_x[0] / df,
                                       self.plot_x[-1] / df])
                col_line_xlist.append(None)
                col_line_ylist.extend([yy, yy])
                col_line_ylist.append(None)

            # plot rows
            ax.plot(col_line_xlist,
                    col_line_ylist,
                    color='k',
                    lw=.5)

            diag_line_xlist = []
            diag_line_ylist = []
            for xi, xx in enumerate(self.plot_x[:-1] / df):
                for yi, yy in enumerate(self.plot_z[:-1] / df):
                    diag_line_xlist.extend([xx, self.plot_x[xi + 1] / df])
                    diag_line_xlist.append(None)
                    diag_line_xlist.extend([xx, self.plot_x[xi + 1] / df])
                    diag_line_xlist.append(None)

                    diag_line_ylist.extend([yy, self.plot_z[yi + 1] / df])
                    diag_line_ylist.append(None)
                    diag_line_ylist.extend([self.plot_z[yi + 1] / df, yy])
                    diag_line_ylist.append(None)

            # plot diagonal lines.
            ax.plot(diag_line_xlist,
                    diag_line_ylist,
                    color='k',
                    lw=.5)

        # plot the regularization mesh
        if self.regmesh == 'on':
            line_list = []
            for ii in range(len(self.model_rows)):
                # get the number of layers to combine
                # this index will be the first index in the vertical direction
                ny1 = self.model_rows[:ii, 0].sum()

                # the second index  in the vertical direction
                ny2 = ny1 + self.model_rows[ii][0]

                # make the list of amalgamated columns an array for ease
                lc = np.array(self.model_cols[ii])
                yline = ax.plot([self.plot_x[0] / df, self.plot_x[-1] / df],
                                [self.plot_z[-ny1] / df,
                                 self.plot_z[-ny1] / df],
                                color='b',
                                lw=.5)

                line_list.append(yline)

                # loop over the number of amalgamated blocks
                for jj in range(len(self.model_cols[ii])):
                    # get first in index in the horizontal direction
                    nx1 = lc[:jj].sum()

                    # get second index in horizontal direction
                    nx2 = nx1 + lc[jj]
                    try:
                        if ny1 == 0:
                            ny1 = 1
                        xline = ax.plot([self.plot_x[nx1] / df,
                                         self.plot_x[nx1] / df],
                                        [self.plot_z[-ny1] / df,
                                         self.plot_z[-ny2] / df],
                                        color='b',
                                        lw=.5)
                        line_list.append(xline)
                    except IndexError:
                        pass

        ##plot the mesh block numbers
        if self.meshnum == 'on':
            kk = 1
            for yy in self.plot_z[::-1] / df:
                for xx in self.plot_x / df:
                    ax.text(xx, yy, '{0}'.format(kk),
                            fontdict={'size': self.meshnum_font_size})
                    kk += 1

        ##plot regularization block numbers
        if self.blocknum == 'on':
            kk = 1
            for ii in range(len(self.model_rows)):
                # get the number of layers to combine
                # this index will be the first index in the vertical direction
                ny1 = self.model_rows[:ii, 0].sum()

                # the second index  in the vertical direction
                ny2 = ny1 + self.model_rows[ii][0]
                # make the list of amalgamated columns an array for ease
                lc = np.array(self.model_cols[ii])
                # loop over the number of amalgamated blocks
                for jj in range(len(self.model_cols[ii])):
                    # get first in index in the horizontal direction
                    nx1 = lc[:jj].sum()
                    # get second index in horizontal direction
                    nx2 = nx1 + lc[jj]
                    try:
                        if ny1 == 0:
                            ny1 = 1
                        # get center points of the blocks
                        yy = self.plot_z[-ny1] - (self.plot_z[-ny1] -
                                                  self.plot_z[-ny2]) / 2
                        xx = self.plot_x[nx1] - \
                             (self.plot_x[nx1] - self.plot_x[nx2]) / 2
                        # put the number
                        ax.text(xx / df, yy / df, '{0}'.format(kk),
                                fontdict={'size': self.block_font_size},
                                horizontalalignment='center',
                                verticalalignment='center')
                        kk += 1
                    except IndexError:
                        pass

        plt.show()

        # make attributes that can be manipulated
        self.ax = ax
        self.cbax = cb

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
                    fig_dpi=None, close_fig='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> model_plot.save_figure(r"/home/occam/figures", 
                                       file_format='jpg')
            
        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if not os.path.isdir(save_fn):
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, 'OccamModel.' +
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots the resistivity model found by Occam2D.")


# ==============================================================================
# plot L2 curve of iteration vs rms
# ==============================================================================
class PlotL2():
    """
    Plot L2 curve of iteration vs rms and rms vs roughness.
    
    Need to only input an .iter file, will read all similar .iter files
    to get the rms, iteration number and roughness of all similar .iter files.
    
    Arguments:
    ----------
        **iter_fn** : string
                      full path to an iteration file output by Occam2D.
                      
    ======================= ===================================================
    Keywords/attributes     Description
    ======================= ===================================================
    ax1                     matplotlib.axes instance for rms vs iteration
    ax2                     matplotlib.axes instance for roughness vs rms
    fig                     matplotlib.figure instance
    fig_dpi                 resolution of figure in dots-per-inch
    fig_num                 number of figure instance
    fig_size                size of figure in inches (width, height)
    font_size               size of axes tick labels, axes labels is +2
    plot_yn                 [ 'y' | 'n']
                            'y' --> to plot on instantiation
                            'n' --> to not plot on instantiation
    rms_arr                 structure np.array as described above
    rms_color               color of rms marker and line
    rms_lw                  line width of rms line
    rms_marker              marker for rms values
    rms_marker_size         size of marker for rms values
    rms_mean_color          color of mean line
    rms_median_color        color of median line
    rough_color             color of roughness line and marker
    rough_font_size         font size for iteration number inside roughness 
                            marker
    rough_lw                line width for roughness line 
    rough_marker            marker for roughness
    rough_marker_size       size of marker for roughness
    subplot_bottom          subplot spacing from bottom  
    subplot_left            subplot spacing from left  
    subplot_right           subplot spacing from right
    subplot_top             subplot spacing from top
    ======================= ===================================================
   
    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots L2 curve.  
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired 
                        location and format
    =================== ======================================================
     
    """

    def __init__(self, iter_fn, **kwargs):
        self.iter_path = os.path.dirname(iter_fn)
        self.iter_basename = os.path.basename(iter_fn)[:-7]
        self.iter_fn_list = []
        self.rms_arr = None
        self.rough_arr = None

        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .91
        self.subplot_bottom = .1

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.font_size = kwargs.pop('font_size', 8)

        self.rms_lw = kwargs.pop('rms_lw', 1)
        self.rms_marker = kwargs.pop('rms_marker', 'd')
        self.rms_color = kwargs.pop('rms_color', 'k')
        self.rms_marker_size = kwargs.pop('rms_marker_size', 5)
        self.rms_median_color = kwargs.pop('rms_median_color', 'red')
        self.rms_mean_color = kwargs.pop('rms_mean_color', 'orange')

        self.rough_lw = kwargs.pop('rough_lw', .75)
        self.rough_marker = kwargs.pop('rough_marker', 'o')
        self.rough_color = kwargs.pop('rough_color', 'b')
        self.rough_marker_size = kwargs.pop('rough_marker_size', 7)
        self.rough_font_size = kwargs.pop('rough_font_size', 6)

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def _get_iterfn_list(self):
        """
        get all iteration files for a given inversion
        
        """

        self.iter_fn_list = [os.path.join(self.iter_path, fn)
                             for fn in os.listdir(self.iter_path)
                             if fn.find(self.iter_basename) == 0 and
                             fn.find('.iter') > 0]

    def _get_values(self):
        """
        get rms and roughness values from iteration files
        """
        self._get_iterfn_list()
        self.rms_arr = np.zeros((len(self.iter_fn_list), 2))
        self.rough_arr = np.zeros((len(self.iter_fn_list), 2))

        for ii, itfn in enumerate(self.iter_fn_list):
            m_object = Model(itfn)
            m_object.read_iter_file()
            m_index = int(m_object.iteration)
            self.rms_arr[ii, 1] = float(m_object.misfit_value)
            self.rms_arr[ii, 0] = m_index
            self.rough_arr[ii, 1] = float(m_object.roughness_value)
            self.rough_arr[ii, 0] = m_index

            # sort by iteration number
            #        self.rms_arr = np.sort(self.rms_arr, axis=1)
            #        self.rough_arr = np.sort(self.rough_arr, axis=1)

    def plot(self):
        """
        plot L2 curve
        """

        self._get_values()
        nr = self.rms_arr.shape[0]
        med_rms = np.median(self.rms_arr[1:, 1])
        mean_rms = np.mean(self.rms_arr[1:, 1])

        # set the dimesions of the figure
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top

        # make figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        # make a subplot for RMS vs Iteration
        self.ax1 = self.fig.add_subplot(1, 1, 1)

        # plot the rms vs iteration
        l1, = self.ax1.plot(self.rms_arr[:, 0],
                            self.rms_arr[:, 1],
                            '-k',
                            lw=1,
                            marker='d',
                            ms=5)

        # plot the median of the RMS
        m1, = self.ax1.plot(self.rms_arr[:, 0],
                            np.repeat(med_rms, nr),
                            ls='--',
                            color=self.rms_median_color,
                            lw=self.rms_lw * .75)

        # plot the mean of the RMS
        m2, = self.ax1.plot(self.rms_arr[:, 0],
                            np.repeat(mean_rms, nr),
                            ls='--',
                            color=self.rms_mean_color,
                            lw=self.rms_lw * .75)

        # make subplot for RMS vs Roughness Plot
        self.ax2 = self.ax1.twiny()

        self.ax2.set_xlim(self.rough_arr[1:, 1].min(),
                          self.rough_arr[1:, 1].max())

        self.ax1.set_ylim(np.floor(self.rms_arr[1:, 1].min()),
                          self.rms_arr[1:, 1].max())

        # plot the rms vs roughness
        l2, = self.ax2.plot(self.rough_arr[:, 1],
                            self.rms_arr[:, 1],
                            ls='--',
                            color=self.rough_color,
                            lw=self.rough_lw,
                            marker=self.rough_marker,
                            ms=self.rough_marker_size,
                            mfc='white')

        # plot the iteration number inside the roughness marker
        for rms, ii, rough in zip(self.rms_arr[:, 1], self.rms_arr[:, 0],
                                  self.rough_arr[:, 1]):
            # need this because if the roughness is larger than this number
            # matplotlib puts the text out of bounds and a draw_text_image
            # error is raised and file cannot be saved, also the other
            # numbers are not put in.
            if rough > 1e8:
                pass
            else:
                self.ax2.text(rough,
                              rms,
                              '{0:.0f}'.format(ii),
                              horizontalalignment='center',
                              verticalalignment='center',
                              fontdict={'size': self.rough_font_size,
                                        'weight': 'bold',
                                        'color': self.rough_color})

        # make a legend
        self.ax1.legend([l1, l2, m1, m2],
                        ['RMS', 'Roughness',
                         'Median_RMS={0:.2f}'.format(med_rms),
                         'Mean_RMS={0:.2f}'.format(mean_rms)],
                        ncol=1,
                        loc='upper right',
                        columnspacing=.25,
                        markerscale=.75,
                        handletextpad=.15)

        # set the axis properties for RMS vs iteration
        self.ax1.yaxis.set_minor_locator(MultipleLocator(.1))
        self.ax1.xaxis.set_minor_locator(MultipleLocator(1))
        self.ax1.xaxis.set_major_locator(MultipleLocator(1))
        self.ax1.set_ylabel('RMS',
                            fontdict={'size': self.font_size + 2,
                                      'weight': 'bold'})
        self.ax1.set_xlabel('Iteration',
                            fontdict={'size': self.font_size + 2,
                                      'weight': 'bold'})
        self.ax1.grid(alpha=.25, which='both', lw=self.rough_lw)
        self.ax2.set_xlabel('Roughness',
                            fontdict={'size': self.font_size + 2,
                                      'weight': 'bold',
                                      'color': self.rough_color})

        for t2 in self.ax2.get_xticklabels():
            t2.set_color(self.rough_color)

        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
                    fig_dpi=None, close_fig='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, '_L2.' +
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots RMS vs Iteration computed by Occam2D")


# ==============================================================================
# plot pseudo section of data and model response                
# ==============================================================================
class PlotPseudoSection(object):
    """
    plot a pseudo section of the data and response if given
    
        
    Arguments:
    -------------
        **data_fn** : string
                      full path to data file.
        
        **resp_fn** : string
                      full path to response file.
    
    ==================== ======================================================
    key words            description
    ==================== ======================================================
    axmpte               matplotlib.axes instance for TE model phase
    axmptm               matplotlib.axes instance for TM model phase
    axmrte               matplotlib.axes instance for TE model app. res 
    axmrtm               matplotlib.axes instance for TM model app. res 
    axpte                matplotlib.axes instance for TE data phase 
    axptm                matplotlib.axes instance for TM data phase
    axrte                matplotlib.axes instance for TE data app. res.
    axrtm                matplotlib.axes instance for TM data app. res.
    cb_pad               padding between colorbar and axes
    cb_shrink            percentage to shrink the colorbar to
    fig                  matplotlib.figure instance
    fig_dpi              resolution of figure in dots per inch
    fig_num              number of figure instance
    fig_size             size of figure in inches (width, height)
    font_size            size of font in points
    label_list            list to label plots
    ml                   factor to label stations if 2 every other station
                         is labeled on the x-axis
    period               np.array of periods to plot
    phase_cmap           color map name of phase
    phase_limits_te      limits for te phase in degrees (min, max)
    phase_limits_tm      limits for tm phase in degrees (min, max)            
    plot_resp            [ 'y' | 'n' ] to plot response
    plot_tipper          [ 'y' | 'n' ] to plot tipper
    plot_yn              [ 'y' | 'n' ] 'y' to plot on instantiation
    res_cmap             color map name for resistivity
    res_limits_te        limits for te resistivity in log scale (min, max)
    res_limits_tm        limits for tm resistivity in log scale (min, max)
    rp_list               list of dictionaries as made from read2Dresp
    station_id           index to get station name (min, max)
    station_list          station list got from rp_list
    subplot_bottom       subplot spacing from bottom (relative coordinates) 
    subplot_hspace       vertical spacing between subplots
    subplot_left         subplot spacing from left  
    subplot_right        subplot spacing from right
    subplot_top          subplot spacing from top
    subplot_wspace       horizontal spacing between subplots
    ==================== ======================================================
    
    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots a pseudo-section of apparent resistiviy and phase
                        of data and model if given.  called on instantiation 
                        if plot_yn is 'y'.
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired 
                        location and format
    =================== =======================================================
                    
   :Example: ::
        
        >>> import mtpy.modeling.occam2d as occam2d
        >>> r_fn = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
        >>> d_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
        >>> ps_plot = occam2d.PlotPseudoSection(d_fn, r_fn) 
    
    """

    def __init__(self, data_fn, resp_fn=None, **kwargs):

        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.plot_resp = kwargs.pop('plot_resp', 'y')
        if self.resp_fn is None:
            self.plot_resp = 'n'

        self.label_list = [r'$\rho_{TE-Data}$', r'$\rho_{TE-Model}$',
                           r'$\rho_{TM-Data}$', r'$\rho_{TM-Model}$',
                           '$\phi_{TE-Data}$', '$\phi_{TE-Model}$',
                           '$\phi_{TM-Data}$', '$\phi_{TM-Model}$',
                           '$\Re e\{T_{Data}\}$', '$\Re e\{T_{Model}\}$',
                           '$\Im m\{T_{Data}\}$', '$\Im m\{T_{Model}\}$']

        self.phase_limits_te = kwargs.pop('phase_limits_te', (-5, 95))
        self.phase_limits_tm = kwargs.pop('phase_limits_tm', (-5, 95))
        self.res_limits_te = kwargs.pop('res_limits_te', (0, 3))
        self.res_limits_tm = kwargs.pop('res_limits_tm', (0, 3))
        self.tip_limits_re = kwargs.pop('tip_limits_re', (-1, 1))
        self.tip_limits_im = kwargs.pop('tip_limits_im', (-1, 1))

        self.phase_cmap = kwargs.pop('phase_cmap', 'jet')
        self.res_cmap = kwargs.pop('res_cmap', 'jet_r')
        self.tip_cmap = kwargs.pop('res_cmap', 'Spectral')

        self.ml = kwargs.pop('ml', 2)
        self.station_id = kwargs.pop('station_id', [0, 4])

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = .025
        self.subplot_hspace = .0
        self.subplot_right = .95
        self.subplot_left = .085
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_num = kwargs.pop('plot_num', 2)
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        self.cb_shrink = .7
        self.cb_pad = .015

        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None
        self.axmrte = None
        self.axmrtm = None
        self.axmpte = None
        self.axmptm = None
        self.axtpr = None
        self.axtpi = None
        self.axmtpr = None
        self.axmtpi = None

        self.te_res_arr = None
        self.tm_res_arr = None
        self.te_phase_arr = None
        self.tm_phase_arr = None
        self.tip_real_arr = None
        self.tip_imag_arr = None

        self.fig = None

        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plot pseudo section of data and response if given
        
        """
        if self.plot_resp == 'y':
            nr = 2
        else:
            nr = 1

        data_obj = Data()
        data_obj.read_data_file(self.data_fn)

        if self.resp_fn is not None:
            resp_obj = Response()
            resp_obj.read_response_file(self.resp_fn)

        ns = len(data_obj.station_list)
        nf = len(data_obj.period)
        ylimits = (data_obj.period.max(), data_obj.period.min())

        # make a grid for pcolormesh so you can have a log scale
        # get things into arrays for plotting
        offset_list = np.zeros(ns + 1)
        te_res_arr = np.ones((nf, ns, nr))
        tm_res_arr = np.ones((nf, ns, nr))
        te_phase_arr = np.zeros((nf, ns, nr))
        tm_phase_arr = np.zeros((nf, ns, nr))
        tip_real_arr = np.zeros((nf, ns, nr))
        tip_imag_arr = np.zeros((nf, ns, nr))

        for ii, d_dict in enumerate(data_obj.data):
            offset_list[ii] = d_dict['offset']
            te_res_arr[:, ii, 0] = d_dict['te_res'][0]
            tm_res_arr[:, ii, 0] = d_dict['tm_res'][0]
            te_phase_arr[:, ii, 0] = d_dict['te_phase'][0]
            tm_phase_arr[:, ii, 0] = d_dict['tm_phase'][0]
            tip_real_arr[:, ii, 0] = d_dict['re_tip'][0]
            tip_imag_arr[:, ii, 0] = d_dict['im_tip'][0]

        # read in response data
        if self.plot_resp == 'y':
            for ii, r_dict in enumerate(resp_obj.resp):
                te_res_arr[:, ii, 1] = r_dict['te_res'][0]
                tm_res_arr[:, ii, 1] = r_dict['tm_res'][0]
                te_phase_arr[:, ii, 1] = r_dict['te_phase'][0]
                tm_phase_arr[:, ii, 1] = r_dict['tm_phase'][0]
                tip_real_arr[:, ii, 1] = r_dict['re_tip'][0]
                tip_imag_arr[:, ii, 1] = r_dict['im_tip'][0]

        # need to make any zeros 1 for taking log10
        te_res_arr[np.where(te_res_arr == 0)] = 1.0
        tm_res_arr[np.where(tm_res_arr == 0)] = 1.0

        self.te_res_arr = te_res_arr
        self.tm_res_arr = tm_res_arr
        self.te_phase_arr = te_phase_arr
        self.tm_phase_arr = tm_phase_arr
        self.tip_real_arr = tip_real_arr
        self.tip_imag_arr = tip_imag_arr

        # need to extend the last grid cell because meshgrid expects n+1 cells
        offset_list[-1] = offset_list[-2] * 1.15
        # make a meshgrid for plotting
        # flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offset_list, data_obj.period[::-1])

        # make list for station labels
        sindex_1 = self.station_id[0]
        sindex_2 = self.station_id[1]
        slabel = [data_obj.station_list[ss][sindex_1:sindex_2]
                  for ss in range(0, ns, self.ml)]

        xloc = offset_list[0] + abs(offset_list[0] - offset_list[1]) / 5
        yloc = 1.10 * data_obj.period[1]

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.left'] = self.subplot_left

        log_labels_te = ['10$^{0}$'.format('{' + str(nn) + '}')
                         for nn in np.arange(int(self.res_limits_te[0]),
                                             int(self.res_limits_te[1]) + 1)]
        log_labels_tm = ['10$^{0}$'.format('{' + str(nn) + '}')
                         for nn in np.arange(int(self.res_limits_tm[0]),
                                             int(self.res_limits_tm[1]) + 1)]

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        if self.plot_resp == 'y':
            if self.plot_tipper == 'y':
                gs1 = gridspec.GridSpec(1, 3,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        wspace=self.subplot_wspace)
                gs4 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                       hspace=self.subplot_hspace,
                                                       wspace=0,
                                                       subplot_spec=gs1[2])
            else:
                gs1 = gridspec.GridSpec(1, 2,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        wspace=self.subplot_wspace)
            gs2 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                   hspace=self.subplot_hspace,
                                                   wspace=0,
                                                   subplot_spec=gs1[0])
            gs3 = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                   hspace=self.subplot_hspace,
                                                   wspace=0,
                                                   subplot_spec=gs1[1])

            # plot TE resistivity data
            self.axrte = plt.Subplot(self.fig, gs2[0, 0])
            self.fig.add_subplot(self.axrte)
            self.axrte.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(np.log10(te_res_arr[:, :, 0])),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_te[0],
                                  vmax=self.res_limits_te[1])

            # plot TE resistivity model
            self.axmrte = plt.Subplot(self.fig, gs2[0, 1])
            self.fig.add_subplot(self.axmrte)
            self.axmrte.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(np.log10(te_res_arr[:, :, 1])),
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits_te[0],
                                   vmax=self.res_limits_te[1])

            # plot TM resistivity data
            self.axrtm = plt.Subplot(self.fig, gs3[0, 0])
            self.fig.add_subplot(self.axrtm)
            self.axrtm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(np.log10(tm_res_arr[:, :, 0])),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_tm[0],
                                  vmax=self.res_limits_tm[1])

            # plot TM resistivity model
            self.axmrtm = plt.Subplot(self.fig, gs3[0, 1])
            self.fig.add_subplot(self.axmrtm)
            self.axmrtm.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(np.log10(tm_res_arr[:, :, 1])),
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits_tm[0],
                                   vmax=self.res_limits_tm[1])

            # plot TE phase data
            self.axpte = plt.Subplot(self.fig, gs2[1, 0])
            self.fig.add_subplot(self.axpte)
            self.axpte.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(te_phase_arr[:, :, 0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_te[0],
                                  vmax=self.phase_limits_te[1])

            # plot TE phase model
            self.axmpte = plt.Subplot(self.fig, gs2[1, 1])
            self.fig.add_subplot(self.axmpte)
            self.axmpte.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(te_phase_arr[:, :, 1]),
                                   cmap=self.phase_cmap,
                                   vmin=self.phase_limits_te[0],
                                   vmax=self.phase_limits_te[1])

            # plot TM phase data
            self.axptm = plt.Subplot(self.fig, gs3[1, 0])
            self.fig.add_subplot(self.axptm)
            self.axptm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(tm_phase_arr[:, :, 0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_tm[0],
                                  vmax=self.phase_limits_tm[1])

            # plot TM phase model
            self.axmptm = plt.Subplot(self.fig, gs3[1, 1])
            self.fig.add_subplot(self.axmptm)
            self.axmptm.pcolormesh(dgrid,
                                   fgrid,
                                   np.flipud(tm_phase_arr[:, :, 1]),
                                   cmap=self.phase_cmap,
                                   vmin=self.phase_limits_tm[0],
                                   vmax=self.phase_limits_tm[1])

            ax_list = [self.axrte, self.axmrte, self.axrtm, self.axmrtm,
                       self.axpte, self.axmpte, self.axptm, self.axmptm]

            if self.plot_tipper == 'y':
                # plot real tipper  data
                self.axtpr = plt.Subplot(self.fig, gs4[0, 0])
                self.fig.add_subplot(self.axtpr)
                self.axtpr.pcolormesh(dgrid,
                                      fgrid,
                                      np.flipud(tip_real_arr[:, :, 0]),
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                # plot real tipper  model
                self.axmtpr = plt.Subplot(self.fig, gs4[0, 1])
                self.fig.add_subplot(self.axmtpr)
                self.axmtpr.pcolormesh(dgrid,
                                       fgrid,
                                       np.flipud(tip_real_arr[:, :, 1]),
                                       cmap=self.tip_cmap,
                                       vmin=self.tip_limits_re[0],
                                       vmax=self.tip_limits_re[1])

                # plot imag tipper  data
                self.axtpi = plt.Subplot(self.fig, gs4[1, 0])
                self.fig.add_subplot(self.axtpi)
                self.axtpi.pcolormesh(dgrid,
                                      fgrid,
                                      np.flipud(tip_imag_arr[:, :, 0]),
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                # plot imag tipper  model
                self.axmtpi = plt.Subplot(self.fig, gs4[1, 1])
                self.fig.add_subplot(self.axmtpi)
                self.axmtpi.pcolormesh(dgrid,
                                       fgrid,
                                       np.flipud(tip_imag_arr[:, :, 1]),
                                       cmap=self.tip_cmap,
                                       vmin=self.tip_limits_re[0],
                                       vmax=self.tip_limits_re[1])

                ax_list.append(self.axtpr)
                ax_list.append(self.axmtpr)
                ax_list.append(self.axtpi)
                ax_list.append(self.axmtpi)

            # make everthing look tidy
            for xx, ax in enumerate(ax_list):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_list, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.set_xlim(offset_list.min(), offset_list.max())
                if np.remainder(xx, 2.0) == 1:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cbx = mcb.make_axes(ax,
                                        shrink=self.cb_shrink,
                                        pad=self.cb_pad)
                if xx == 2 or xx == 6 or xx == 8 or xx == 10:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                if xx < 4:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    if xx == 1:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                              norm=Normalize(vmin=self.res_limits_te[0],
                                                             vmax=self.res_limits_te[1]))
                        cb.set_ticks(np.arange(int(self.res_limits_te[0]),
                                               int(self.res_limits_te[1]) + 1))
                        cb.set_ticklabels(log_labels_te)
                    if xx == 3:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                              norm=Normalize(vmin=self.res_limits_tm[0],
                                                             vmax=self.res_limits_tm[1]))
                        cb.set_label('App. Res. ($\Omega \cdot$m)',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                        cb.set_label('Resistivity ($\Omega \cdot$m)',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                        cb.set_ticks(np.arange(int(self.res_limits_tm[0]),
                                               int(self.res_limits_tm[1]) + 1))
                        cb.set_ticklabels(log_labels_tm)
                else:
                    # color bar TE phase
                    if xx == 5:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                              norm=Normalize(vmin=self.phase_limits_te[0],
                                                             vmax=self.phase_limits_te[1]))
                    # color bar TM phase
                    if xx == 7:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                              norm=Normalize(vmin=self.phase_limits_tm[0],
                                                             vmax=self.phase_limits_tm[1]))
                        cb.set_label('Phase (deg)',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                    # color bar tipper Imag
                    if xx == 9:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                              norm=Normalize(vmin=self.tip_limits_re[0],
                                                             vmax=self.tip_limits_re[1]))
                        cb.set_label('Re{T}',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})
                    if xx == 11:
                        cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                              norm=Normalize(vmin=self.tip_limits_im[0],
                                                             vmax=self.tip_limits_im[1]))
                        cb.set_label('Im{T}',
                                     fontdict={'size': self.font_size + 1,
                                               'weight': 'bold'})

                ax.text(xloc, yloc, self.label_list[xx],
                        fontdict={'size': self.font_size + 1},
                        bbox={'facecolor': 'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 4:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size': self.font_size + 2,
                                            'weight': 'bold'})
                if xx > 3:
                    ax.set_xlabel('Station', fontdict={'size': self.font_size + 2,
                                                       'weight': 'bold'})

            plt.show()

        else:
            if self.plot_tipper == 'y':
                gs1 = gridspec.GridSpec(2, 3,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        hspace=self.subplot_hspace,
                                        wspace=self.subplot_wspace)

            else:
                gs1 = gridspec.GridSpec(2, 2,
                                        left=self.subplot_left,
                                        right=self.subplot_right,
                                        hspace=self.subplot_hspace,
                                        wspace=self.subplot_wspace)

            # plot TE resistivity data
            self.axrte = self.fig.add_subplot(gs1[0, 0])
            self.axrte.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(np.log10(te_res_arr[:, :, 0])),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_te[0],
                                  vmax=self.res_limits_te[1])

            # plot TM resistivity data
            self.axrtm = self.fig.add_subplot(gs1[0, 1])
            self.axrtm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(np.log10(tm_res_arr[:, :, 0])),
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits_tm[0],
                                  vmax=self.res_limits_tm[1])

            # plot TE phase data
            self.axpte = self.fig.add_subplot(gs1[1, 0])
            self.axpte.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(te_phase_arr[:, :, 0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_te[0],
                                  vmax=self.phase_limits_te[1])

            # plot TM phase data
            self.axptm = self.fig.add_subplot(gs1[1, 1])
            self.axptm.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(tm_phase_arr[:, :, 0]),
                                  cmap=self.phase_cmap,
                                  vmin=self.phase_limits_tm[0],
                                  vmax=self.phase_limits_tm[1])
            ax_list = [self.axrte, self.axrtm, self.axpte, self.axptm]
            if self.plot_tipper == 'y':
                # plot real tipper  data
                self.axtpr = plt.Subplot(self.fig, gs1[0, 2])
                self.fig.add_subplot(self.axtpr)
                self.axtpr.pcolormesh(dgrid,
                                      fgrid,
                                      np.flipud(tip_real_arr[:, :, 0]),
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                # plot real tipper  data
                self.axtpi = plt.Subplot(self.fig, gs1[1, 2])
                self.fig.add_subplot(self.axtpi)
                self.axtpi.pcolormesh(dgrid,
                                      fgrid,
                                      np.flipud(tip_imag_arr[:, :, 0]),
                                      cmap=self.tip_cmap,
                                      vmin=self.tip_limits_re[0],
                                      vmax=self.tip_limits_re[1])
                ax_list.append(self.axtpr)
                ax_list.append(self.axtpi)

            # make everything look tidy
            for xx, ax in enumerate(ax_list):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_list, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.grid(True, alpha=.25)
                ax.set_xlim(offset_list.min(), offset_list.max())
                cbx = mcb.make_axes(ax,
                                    shrink=self.cb_shrink,
                                    pad=self.cb_pad)
                if xx == 0:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits_te[0],
                                                         vmax=self.res_limits_te[1]))
                    cb.set_ticks(np.arange(self.res_limits_te[0],
                                           self.res_limits_te[1] + 1))
                    cb.set_ticklabels(log_labels_te)
                elif xx == 1:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits_tm[0],
                                                         vmax=self.res_limits_tm[1]))
                    cb.set_label('App. Res. ($\Omega \cdot$m)',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})
                    cb.set_ticks(np.arange(self.res_limits_tm[0],
                                           self.res_limits_tm[1] + 1))
                    cb.set_ticklabels(log_labels_tm)
                elif xx == 2:
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                          norm=Normalize(vmin=self.phase_limits_te[0],
                                                         vmax=self.phase_limits_te[1]))
                    cb.set_ticks(np.arange(self.phase_limits_te[0],
                                           self.phase_limits_te[1] + 1, 15))
                elif xx == 3:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                          norm=Normalize(vmin=self.phase_limits_tm[0],
                                                         vmax=self.phase_limits_tm[1]))
                    cb.set_label('Phase (deg)',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})
                    cb.set_ticks(np.arange(self.phase_limits_te[0],
                                           self.phase_limits_te[1] + 1, 15))

                # real tipper
                elif xx == 4:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                          norm=Normalize(vmin=self.tip_limits_re[0],
                                                         vmax=self.tip_limits_re[1]))
                    cb.set_label('Re{T}',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})
                # imag tipper
                elif xx == 5:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                          norm=Normalize(vmin=self.tip_limits_im[0],
                                                         vmax=self.tip_limits_im[1]))
                    cb.set_label('Im{T}',
                                 fontdict={'size': self.font_size + 1,
                                           'weight': 'bold'})

                ax.text(xloc, yloc, self.label_list[2 * xx],
                        fontdict={'size': self.font_size + 1},
                        bbox={'facecolor': 'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 2:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size': self.font_size + 2,
                                            'weight': 'bold'})
                if xx > 1:
                    ax.set_xlabel('Station', fontdict={'size': self.font_size + 2,
                                                       'weight': 'bold'})

            plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotPseudoSection()
            >>> #change color of te markers to a gray-blue
            >>> p1.res_cmap = 'seismic_r'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
                    fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, 'OccamPseudoSection.' +
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots a pseudo section of TE and TM modes for data and "
                "response if given.")


# ==============================================================================
# plot misfits as a pseudo-section
# ==============================================================================
class PlotMisfitPseudoSection(object):
    """
    plot a pseudo section of the data and response if given
    
        
    Arguments:
    -------------
        **rp_list** : list of dictionaries for each station with keywords:
                
                * *station* : string
                             station name
                
                * *offset* : float
                             relative offset
                
                * *resxy* : np.array(nf,4)
                            TE resistivity and error as row 0 and 1 respectively
                
                * *resyx* : np.array(fn,4)
                            TM resistivity and error as row 0 and 1 respectively
                
                * *phasexy* : np.array(nf,4)
                              TE phase and error as row 0 and 1 respectively
                
                * *phaseyx* : np.array(nf,4)
                              Tm phase and error as row 0 and 1 respectively
                
                * *realtip* : np.array(nf,4)
                              Real Tipper and error as row 0 and 1 respectively
                
                * *imagtip* : np.array(nf,4)
                              Imaginary Tipper and error as row 0 and 1 
                              respectively
                
                Note: that the resistivity will be in log10 space.  Also, there
                are 2 extra rows in the data arrays, this is to put the 
                response from the inversion.  
        
        **period** : np.array of periods to plot that correspond to the index
                     values of each rp_list entry ie. resxy.
    
    ==================== ==================================================
    key words            description
    ==================== ==================================================
    axmpte               matplotlib.axes instance for TE model phase
    axmptm               matplotlib.axes instance for TM model phase
    axmrte               matplotlib.axes instance for TE model app. res 
    axmrtm               matplotlib.axes instance for TM model app. res 
    axpte                matplotlib.axes instance for TE data phase 
    axptm                matplotlib.axes instance for TM data phase
    axrte                matplotlib.axes instance for TE data app. res.
    axrtm                matplotlib.axes instance for TM data app. res.
    cb_pad               padding between colorbar and axes
    cb_shrink            percentage to shrink the colorbar to
    fig                  matplotlib.figure instance
    fig_dpi              resolution of figure in dots per inch
    fig_num              number of figure instance
    fig_size             size of figure in inches (width, height)
    font_size            size of font in points
    label_list            list to label plots
    ml                   factor to label stations if 2 every other station
                         is labeled on the x-axis
    period               np.array of periods to plot
    phase_cmap           color map name of phase
    phase_limits_te      limits for te phase in degrees (min, max)
    phase_limits_tm      limits for tm phase in degrees (min, max)            
    plot_resp            [ 'y' | 'n' ] to plot response
    plot_yn              [ 'y' | 'n' ] 'y' to plot on instantiation

    res_cmap             color map name for resistivity
    res_limits_te        limits for te resistivity in log scale (min, max)
    res_limits_tm        limits for tm resistivity in log scale (min, max)
    rp_list               list of dictionaries as made from read2Dresp
    station_id           index to get station name (min, max)
    station_list          station list got from rp_list
    subplot_bottom       subplot spacing from bottom (relative coordinates) 
    subplot_hspace       vertical spacing between subplots
    subplot_left         subplot spacing from left  
    subplot_right        subplot spacing from right
    subplot_top          subplot spacing from top
    subplot_wspace       horizontal spacing between subplots
    ==================== ==================================================
    
    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots a pseudo-section of apparent resistiviy and phase
                        of data and model if given.  called on instantiation 
                        if plot_yn is 'y'.
    redraw_plot         call redraw_plot to redraw the figures, 
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired 
                        location and format
    =================== =======================================================
                    
   :Example: ::
        
        >>> import mtpy.modeling.occam2d as occam2d
        >>> ocd = occam2d.Occam2DData()
        >>> rfile = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
        >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
        >>> ps1 = ocd.plot2PseudoSection(resp_fn=rfile) 
    
    """

    def __init__(self, data_fn, resp_fn, **kwargs):

        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.label_list = [r'$\rho_{TE}$', r'$\rho_{TM}$',
                           '$\phi_{TE}$', '$\phi_{TM}$',
                           '$\Re e\{T\}$', '$\Im m\{T\}$']

        self.phase_limits_te = kwargs.pop('phase_limits_te', (-10, 10))
        self.phase_limits_tm = kwargs.pop('phase_limits_tm', (-10, 10))
        self.res_limits_te = kwargs.pop('res_limits_te', (-2, 2))
        self.res_limits_tm = kwargs.pop('res_limits_tm', (-2, 2))
        self.tip_limits_re = kwargs.pop('tip_limits_re', (-.2, .2))
        self.tip_limits_im = kwargs.pop('tip_limits_im', (-.2, .2))

        self.phase_cmap = kwargs.pop('phase_cmap', 'BrBG')
        self.res_cmap = kwargs.pop('res_cmap', 'BrBG_r')
        self.tip_cmap = kwargs.pop('tip_cmap', 'PuOr')
        self.plot_tipper = kwargs.pop('plot_tipper', 'n')

        self.ml = kwargs.pop('ml', 2)
        self.station_id = kwargs.pop('station_id', [0, 4])

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = .0025
        self.subplot_hspace = .0
        self.subplot_right = .95
        self.subplot_left = .085
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        self.cb_shrink = .7
        self.cb_pad = .015

        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None
        self.axtpr = None
        self.axtpi = None

        self.misfit_te_res = None
        self.misfit_te_phase = None
        self.misfit_tm_res = None
        self.misfit_tm_phase = None
        self.misfit_tip_real = None
        self.misfit_tip_imag = None

        self.fig = None
        self._data_obj = None

        if self.plot_yn == 'y':
            self.plot()

    def get_misfit(self):
        """
        compute misfit of MT response found from the model and the data.
        
        Need to normalize correctly
        """
        data_obj = Data()
        data_obj.read_data_file(self.data_fn)
        self._data_obj = data_obj

        resp_obj = Response()
        resp_obj.read_response_file(self.resp_fn)

        n_stations = len(data_obj.station_list)
        n_periods = len(data_obj.freq)

        self.misfit_te_res = np.zeros((n_periods, n_stations))
        self.misfit_te_phase = np.zeros((n_periods, n_stations))
        self.misfit_tm_res = np.zeros((n_periods, n_stations))
        self.misfit_tm_phase = np.zeros((n_periods, n_stations))
        self.misfit_tip_real = np.zeros((n_periods, n_stations))
        self.misfit_tip_imag = np.zeros((n_periods, n_stations))

        for rr, r_dict in zip(list(range(n_stations)), resp_obj.resp):
            self.misfit_te_res[:, rr] = r_dict['te_res'][1]
            self.misfit_tm_res[:, rr] = r_dict['tm_res'][1]
            self.misfit_te_phase[:, rr] = r_dict['te_phase'][1]
            self.misfit_tm_phase[:, rr] = r_dict['tm_phase'][1]
            self.misfit_tip_real[:, rr] = r_dict['re_tip'][1]
            self.misfit_tip_imag[:, rr] = r_dict['im_tip'][1]

        self.misfit_te_res = np.nan_to_num(self.misfit_te_res)
        self.misfit_te_phase = np.nan_to_num(self.misfit_te_phase)
        self.misfit_tm_res = np.nan_to_num(self.misfit_tm_res)
        self.misfit_tm_phase = np.nan_to_num(self.misfit_tm_phase)
        self.misfit_tip_real = np.nan_to_num(self.misfit_tip_real)
        self.misfit_tip_imag = np.nan_to_num(self.misfit_tip_imag)

    def plot(self):
        """
        plot pseudo section of data and response if given
        
        """

        self.get_misfit()

        ylimits = (self._data_obj.period.max(), self._data_obj.period.min())

        offset_list = np.append(self._data_obj.station_locations,
                                self._data_obj.station_locations[-1] * 1.15)

        # make a meshgrid for plotting
        # flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offset_list, self._data_obj.period[::-1])

        # make list for station labels
        ns = len(self._data_obj.station_list)
        sindex_1 = self.station_id[0]
        sindex_2 = self.station_id[1]
        slabel = [self._data_obj.station_list[ss][sindex_1:sindex_2]
                  for ss in range(0, ns, self.ml)]

        xloc = offset_list[0] + abs(offset_list[0] - offset_list[1]) / 5
        yloc = 1.10 * self._data_obj.period[1]

        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.hspace'] = self.subplot_hspace
        plt.rcParams['figure.subplot.wspace'] = self.subplot_wspace

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        if self.plot_tipper != 'y':
            self.axrte = self.fig.add_subplot(2, 2, 1)
            self.axrtm = self.fig.add_subplot(2, 2, 2, sharex=self.axrte)
            self.axpte = self.fig.add_subplot(2, 2, 3, sharex=self.axrte)
            self.axptm = self.fig.add_subplot(2, 2, 4, sharex=self.axrte)

        else:
            self.axrte = self.fig.add_subplot(2, 3, 1)
            self.axrtm = self.fig.add_subplot(2, 3, 2, sharex=self.axrte)
            self.axpte = self.fig.add_subplot(2, 3, 4, sharex=self.axrte)
            self.axptm = self.fig.add_subplot(2, 3, 5, sharex=self.axrte)
            self.axtpr = self.fig.add_subplot(2, 3, 3, sharex=self.axrte)
            self.axtpi = self.fig.add_subplot(2, 3, 6, sharex=self.axrte)

        # --> TE Resistivity
        self.axrte.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_te_res),
                              cmap=self.res_cmap,
                              vmin=self.res_limits_te[0],
                              vmax=self.res_limits_te[1])
        # --> TM Resistivity
        self.axrtm.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_tm_res),
                              cmap=self.res_cmap,
                              vmin=self.res_limits_tm[0],
                              vmax=self.res_limits_tm[1])
        # --> TE Phase
        self.axpte.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_te_phase),
                              cmap=self.phase_cmap,
                              vmin=self.phase_limits_te[0],
                              vmax=self.phase_limits_te[1])
        # --> TM Phase
        self.axptm.pcolormesh(dgrid,
                              fgrid,
                              np.flipud(self.misfit_tm_phase),
                              cmap=self.phase_cmap,
                              vmin=self.phase_limits_tm[0],
                              vmax=self.phase_limits_tm[1])

        ax_list = [self.axrte, self.axrtm, self.axpte, self.axptm]

        if self.plot_tipper == 'y':
            self.axtpr.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(self.misfit_tip_real),
                                  cmap=self.tip_cmap,
                                  vmin=self.tip_limits_re[0],
                                  vmax=self.tip_limits_re[1])
            self.axtpi.pcolormesh(dgrid,
                                  fgrid,
                                  np.flipud(self.misfit_tip_imag),
                                  cmap=self.tip_cmap,
                                  vmin=self.tip_limits_im[0],
                                  vmax=self.tip_limits_im[1])

            ax_list.append(self.axtpr)
            ax_list.append(self.axtpi)
            # make everthing look tidy
        for xx, ax in enumerate(ax_list):
            ax.semilogy()
            ax.set_ylim(ylimits)
            ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
            ax.xaxis.set_ticks(offset_list, minor=True)
            ax.xaxis.set_ticklabels(slabel)
            ax.set_xlim(offset_list.min(), offset_list.max())
            cbx = mcb.make_axes(ax,
                                shrink=self.cb_shrink,
                                pad=self.cb_pad)

            # te res
            if xx == 0:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                      norm=Normalize(vmin=self.res_limits_te[0],
                                                     vmax=self.res_limits_te[1]))
            # tm res
            elif xx == 1:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.res_cmap,
                                      norm=Normalize(vmin=self.res_limits_tm[0],
                                                     vmax=self.res_limits_tm[1]))
                cb.set_label('Log$_{10}$ App. Res. ($\Omega \cdot$m)',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})
            # te phase
            elif xx == 2:
                cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                      norm=Normalize(vmin=self.phase_limits_te[0],
                                                     vmax=self.phase_limits_te[1]))
            # tm phase
            elif xx == 3:
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.phase_cmap,
                                      norm=Normalize(vmin=self.phase_limits_tm[0],
                                                     vmax=self.phase_limits_tm[1]))
                cb.set_label('Phase (deg)',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})

            # real tipper
            elif xx == 4:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                      norm=Normalize(vmin=self.tip_limits_re[0],
                                                     vmax=self.tip_limits_re[1]))
                cb.set_label('Re{Tip}',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})
            # imag tipper
            elif xx == 5:
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(cbx[0], cmap=self.tip_cmap,
                                      norm=Normalize(vmin=self.tip_limits_im[0],
                                                     vmax=self.tip_limits_im[1]))
                cb.set_label('Im{Tip}',
                             fontdict={'size': self.font_size + 1,
                                       'weight': 'bold'})

            # make label for plot
            ax.text(xloc, yloc, self.label_list[xx],
                    fontdict={'size': self.font_size + 2},
                    bbox={'facecolor': 'white'},
                    horizontalalignment='left',
                    verticalalignment='top')

            if xx == 0 or xx == 2:
                ax.set_ylabel('Period (s)',
                              fontdict={'size': self.font_size + 2,
                                        'weight': 'bold'})
            if xx > 1:
                ax.set_xlabel('Station', fontdict={'size': self.font_size + 2,
                                                   'weight': 'bold'})

        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotPseudoSection()
            >>> #change color of te markers to a gray-blue
            >>> p1.res_cmap = 'seismic_r'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
                    fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, 'OccamMisfitPseudoSection.' +
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots a pseudo section of TE and TM modes for data and "
                "response if given.")


class OccamPointPicker(object):
    """
    This class helps the user interactively pick points to mask and add 
    error bars. 
    
    Useage:
    -------
    To mask just a single point right click over the point and a gray point 
    will appear indicating it has been masked
    
    To mask both the apparent resistivity and phase left click over the point.
    Gray points will appear over both the apparent resistivity and phase.  
    Sometimes the points don't exactly matchup, haven't quite worked that bug
    out yet, but not to worry it picks out the correct points
    
    To add error bars to a point click the middle or scroll bar button.  This
    only adds error bars to the point and does not reduce them so start out
    with reasonable errorbars.  You can change the increment that the error
    bars are increased with res_err_inc and phase_err_inc.
    
    .. note:: There is a bug when only plotting TE or TM that you cannot mask
              points in the phase.  I'm not sure where it comes from, but
              works with all modes.  So my suggestion is to make a data file 
              with all modes, mask data points and then rewrite that data file
              if you want to use just one of the modes.  That's the work 
              around for the moment.
    
    
    
    Arguments:
    ----------
        **ax_list** : list of the resistivity and phase axis that have been 
                    plotted as [axr_te,axr_tm,axp_te,axp_tm]
        
        **line_list** : list of lines used to plot the responses, not the error 
                      bars as [res_te,res_tm,phase_te,phase_tm]
        
        **err_list** : list of the errorcaps and errorbar lines as 
                   [[cap1,cap2,bar],...]
                 
        **res_err_inc** : increment to increase the errorbars for resistivity.
                        put .20 for 20 percent change. *Default* is .05
        
        **phase_err_inc** : increment to increase the errorbars for the phase
                          put .10 for 10 percent change. *Defualt* is .02 
                    
        **marker** : marker type for masked points.  See matplotlib.pyplot.plot
                    for options of markers.  *Default* is h for hexagon.
                   
    Attributes:
    -----------
    
        **ax_list** : axes list used to plot the data
        
        **line_list** : line list used to plot the data
        
        **err_list** : error list used to plot the data
        
        **data** : list of data points that were not masked for each plot.
        
        **fdict** : dictionary of frequency arrays for each plot and data set.
        
        **fndict** : dictionary of figure numbers to corresponed with data.
        
        **cid_list** : list of event ids.
        
        **res_err_inc** : increment to increase resistivity error bars
        
        **phase_inc** : increment to increase phase error bars
        
        **marker** : marker of masked points
        
        **fig_num** : figure numbers
        
        **data_list** : list of lines to write into the occam2d data file.
        
    :Example: ::
        
        >>> ocd = occam2d.Occam2DData()
        >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
        >>> ocd.plotMaskPoints()
    """

    def __init__(self, ax_list, line_list, err_list,
                 res_err_inc=.05, phase_err_inc=.02, marker='h'):

        # give the class some attributes
        self.ax_list = ax_list
        self.line_list = line_list
        self.err_list = err_list
        self.data = []
        self.error = []
        self.fdict = []
        self.fndict = {}
        # see if just one figure is plotted or multiple figures are plotted
        self.ax = ax_list[0][0]
        self.line = line_list[0][0]
        self.cidlist = []
        self.ax_num = None
        self.res_index = None
        self.phase_index = None
        self.fig_num = 0
        for nn in range(len(ax_list)):
            self.data.append([])
            self.error.append([])
            self.fdict.append([])

            # get data from lines and make a dictionary of frequency points for
            # easy indexing
            # line_find = False
            for ii, line in enumerate(line_list[nn]):
                try:
                    self.data[nn].append(line.get_data()[1])
                    self.fdict[nn].append(dict([('{0:.5g}'.format(kk), ff)
                                                for ff, kk in
                                                enumerate(line.get_data()[0])]))
                    self.fndict['{0}'.format(line.figure.number)] = nn

                    # set some events
                    # if line_find == False:
                    cid1 = line.figure.canvas.mpl_connect('pick_event', self)
                    cid2 = line.figure.canvas.mpl_connect('axes_enter_event',
                                                          self.inAxes)
                    cid3 = line.figure.canvas.mpl_connect('key_press_event',
                                                          self.on_close)
                    cid4 = line.figure.canvas.mpl_connect('figure_enter_event',
                                                          self.inFigure)
                    self.cidlist.append([cid1, cid2, cid3, cid4])

                    # set the figure number
                    self.fig_num = self.line.figure.number

                    # line_find = True

                except AttributeError:
                    self.data[nn].append([])
                    self.fdict[nn].append([])

            # read in the error in a useful way so that it can be translated to
            # the data file.  Make the error into an array
            for ee, err in enumerate(err_list[nn]):
                try:
                    errpath = err[2].get_paths()
                    errarr = np.zeros(len(list(self.fdict[nn][ee].keys())))
                    for ff, epath in enumerate(errpath):
                        errv = epath.vertices
                        errarr[ff] = abs(errv[0, 1] - self.data[nn][ee][ff])
                    self.error[nn].append(errarr)
                except AttributeError:
                    self.error[nn].append([])

        # set the error bar increment values
        self.res_err_inc = res_err_inc
        self.phase_err_inc = phase_err_inc

        # set the marker
        self.marker = marker

        # make a list of occam2d lines to write later
        self.data_list = []

    def __call__(self, event):
        """
        When the function is called the mouse events will be recorder for 
        picking points to mask or change error bars.  The axes is redrawn with
        a gray marker to indicate a masked point and/or increased size in 
        errorbars.
        
        Arguments:
        ----------
            **event** : type mouse_click_event
                
        Useage:
        -------
        
            **Left mouse button** will mask both resistivity and phase point
        
            **Right mouse button** will mask just the point selected
        
            **Middle mouse button** will increase the error bars
        
            **q** will close the figure.
        """
        self.event = event
        # make a new point that is an PickEvent type
        npoint = event.artist
        # if the right button is clicked mask the point
        if event.mouseevent.button == 3:
            # get the point that was clicked on
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]

            # set the x index from the frequency dictionary
            ll = self.fdict[self.fig_num][self.ax_num]['{0:.5g}'.format(xd[0])]

            # change the data to be a zero
            self.data[self.fig_num][self.ax_num][ll] = 0

            # reset the point to be a gray x
            self.ax.plot(xd, yd,
                         ls='None',
                         color=(.7, .7, .7),
                         marker=self.marker,
                         ms=4)

            # redraw the canvas
            self.ax.figure.canvas.draw()

        # if the left button is clicked change both resistivity and phase points
        elif event.mouseevent.button == 1:
            # get the point that was clicked on
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]

            # set the x index from the frequency dictionary
            ll = self.fdict[self.fig_num][self.ax_num]['{0:.5g}'.format(xd[0])]

            # set the data point to zero
            print(self.data[self.fig_num][self.ax_num][ll])
            self.data[self.fig_num][self.ax_num][ll] = 0

            # reset the point to be a gray x
            self.ax.plot(xd, yd,
                         ls='None',
                         color=(.7, .7, .7),
                         marker=self.marker,
                         ms=4)

            self.ax.figure.canvas.draw()

            # check to make sure there is a corresponding res/phase point
            try:
                kk = (self.ax_num + 2) % 4
                print(kk)
                # get the corresponding y-value
                yd2 = self.data[self.fig_num][kk][ll]

                # set that data point to 0 as well
                self.data[self.fig_num][kk][ll] = 0

                # make that data point a gray x
                self.ax_list[self.fig_num][kk].plot(xd, yd2,
                                                    ls='None',
                                                    color=(.7, .7, .7),
                                                    marker=self.marker,
                                                    ms=4)
                # redraw the canvas
                self.ax.figure.canvas.draw()
            except KeyError:
                print('Axis does not contain res/phase point')

        # if click the scroll button or middle button change increase the
        # errorbars by the given amount
        elif event.mouseevent.button == 2:
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]
            jj = self.ax_num

            # get x index
            ll = self.fdict[self.fig_num][jj]['{0:.5g}'.format(xd[0])]

            # make error bar array
            eb = self.err_list[self.fig_num][jj][2].get_paths()[ll].vertices

            # make ecap array
            ecapl = self.err_list[self.fig_num][jj][0].get_data()[1][ll]
            ecapu = self.err_list[self.fig_num][jj][1].get_data()[1][ll]

            # change apparent resistivity error
            if jj == 0 or jj == 1:
                nebu = eb[0, 1] - self.res_err_inc * eb[0, 1]
                nebl = eb[1, 1] + self.res_err_inc * eb[1, 1]
                ecapl = ecapl - self.res_err_inc * ecapl
                ecapu = ecapu + self.res_err_inc * ecapu

            # change phase error
            elif jj == 2 or jj == 3:
                nebu = eb[0, 1] - eb[0, 1] * self.phase_err_inc
                nebl = eb[1, 1] + eb[1, 1] * self.phase_err_inc
                ecapl = ecapl - ecapl * self.phase_err_inc
                ecapu = ecapu + ecapu * self.phase_err_inc

            # put the new error into the error array
            self.error[self.fig_num][jj][ll] = abs(nebu - \
                                                   self.data[self.fig_num][jj][ll])

            # set the new error bar values
            eb[0, 1] = nebu
            eb[1, 1] = nebl

            # reset the error bars and caps
            ncapl = self.err_list[self.fig_num][jj][0].get_data()
            ncapu = self.err_list[self.fig_num][jj][1].get_data()
            ncapl[1][ll] = ecapl
            ncapu[1][ll] = ecapu

            # set the values
            self.err_list[self.fig_num][jj][0].set_data(ncapl)
            self.err_list[self.fig_num][jj][1].set_data(ncapu)
            self.err_list[self.fig_num][jj][2].get_paths()[ll].vertices = eb

            # redraw the canvas
            self.ax.figure.canvas.draw()

    # get the axis number that the mouse is in and change to that axis
    def inAxes(self, event):
        """
        gets the axes that the mouse is currently in.
        
        Arguments:
        ---------
            **event**: is a type axes_enter_event
                
        Returns:
        --------
        
            **OccamPointPicker.jj** : index of resistivity axes for ax_list
            
            **OccamPointPicker.kk** : index of phase axes for ax_list
        
        """

        self.event2 = event
        self.ax = event.inaxes
        for jj, axj in enumerate(self.ax_list):
            for ll, axl in enumerate(axj):
                if self.ax == axl:
                    self.ax_num = ll
        self.line = self.line_list[self.fig_num][self.ax_num]

    # get the figure number that the mouse is in
    def inFigure(self, event):
        """
        gets the figure number that the mouse is in
        
        Arguments:
        ----------
            **event** : figure_enter_event
            
        Returns:
        --------
            **OccamPointPicker.fig_num** : figure number that corresponds to the
                                          index in the ax_list, datalist, errorlist
                                          and line_list.
                        
        """
        self.event3 = event
        self.fig_num = self.fndict['{0}'.format(event.canvas.figure.number)]

    # type the q key to quit the figure and disconnect event handling
    def on_close(self, event):
        """
        close the figure with a 'q' key event and disconnect the event ids
        
        Arguments:
        ----------
            **event** : key_press_event
               
        Returns:
        --------
            print statement saying the figure is closed
        """
        self.event3 = event
        if self.event3.key == 'q':
            for cid in self.cidlist[self.fig_num]:
                event.canvas.mpl_disconnect(cid)
            plt.close(event.canvas.figure)
            print('Closed figure ', self.fig_num)


class Run():
    """
    Run Occam2D by system call.

    Future plan: implement Occam in Python and call it from here directly.
    """


class Mask(Data):
    """
    Allow masking of points from data file (effectively commenting them out, 
    so the process is reversable). Inheriting from Data class.
    """


class OccamInputError(Exception):
    pass
