"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
import os
import time

import numpy as np
from matplotlib import pyplot as plt, gridspec as gridspec, colorbar as mcb
from matplotlib.colors import Normalize
from matplotlib.widgets import Button, RadioButtons, SpanSelector

from mtpy.modeling.modem.data import Data
from mtpy.modeling.modem.data import Model
from mtpy.utils import exceptions as mtex,basemap_tools
from mtpy.utils.gis_tools import get_epsg,epsg_project
from mtpy.utils.calculator import nearest_index
from mtpy.utils.mesh_tools import rotate_mesh

from mtpy.imaging.seismic import Segy, VelocityModel

from scipy.spatial import cKDTree
from scipy.interpolate import interp1d, UnivariateSpline
from matplotlib import colors,cm
from matplotlib.ticker import LogLocator

__all__ = ['PlotSlices']


class PlotSlices(object):
    """
    * Plot all cartesian axis-aligned slices and be able to scroll through the model
    * Extract arbitrary profiles (e.g. along a seismic line) from a model

    :Example: ::

        >>> import mtpy.modeling.modem as modem
        >>> mfn = r"/home/modem/Inv1/Modular_NLCG_100.rho"
        >>> dfn = r"/home/modem/Inv1/ModEM_data.dat"
        >>> pds = ws.PlotSlices(model_fn=mfn, data_fn=dfn)

    ======================= ===================================================
    Buttons                  Description
    ======================= ===================================================
    'e'                     moves n-s slice east by one model block
    'w'                     moves n-s slice west by one model block
    'n'                     moves e-w slice north by one model block
    'm'                     moves e-w slice south by one model block
    'd'                     moves depth slice down by one model block
    'u'                     moves depth slice up by one model block
    ======================= ===================================================


    ======================= ===================================================
    Attributes              Description
    ======================= ===================================================
    ax_en                   matplotlib.axes instance for depth slice  map view
    ax_ez                   matplotlib.axes instance for e-w slice
    ax_map                  matplotlib.axes instance for location map
    ax_nz                   matplotlib.axes instance for n-s slice
    climits                 (min , max) color limits on resistivity in log
                            scale. *default* is (0, 4)
    cmap                    name of color map for resisitiviy.
                            *default* is 'jet_r'
    data_fn                 full path to data file name
    draw_colorbar           show colorbar on exported plot; default True
    dscale                  scaling parameter depending on map_scale
    east_line_xlist         list of line nodes of east grid for faster plotting
    east_line_ylist         list of line nodes of east grid for faster plotting
    ew_limits               (min, max) limits of e-w in map_scale units
                            *default* is None and scales to station area
    fig                     matplotlib.figure instance for figure
    fig_aspect              aspect ratio of plots. *default* is 1
    fig_dpi                 resolution of figure in dots-per-inch
                            *default* is 300
    fig_num                 figure instance number
    fig_size                [width, height] of figure window.
                            *default* is [6,6]
    font_dict               dictionary of font keywords, internally created
    font_size               size of ticklables in points, axes labes are
                            font_size+2. *default* is 4
    grid_east               relative location of grid nodes in e-w direction
                            in map_scale units
    grid_north              relative location of grid nodes in n-s direction
                            in map_scale units
    grid_z                  relative location of grid nodes in z direction
                            in map_scale units
    index_east              index value of grid_east being plotted
    index_north             index value of grid_north being plotted
    index_vertical          index value of grid_z being plotted
    initial_fn              full path to initial file
    key_press               matplotlib.canvas.connect instance
    map_scale               [ 'm' | 'km' ] scale of map. *default* is km
    mesh_east               np.meshgrid(grid_east, grid_north)[0]
    mesh_en_east            np.meshgrid(grid_east, grid_north)[0]
    mesh_en_north           np.meshgrid(grid_east, grid_north)[1]
    mesh_ez_east            np.meshgrid(grid_east, grid_z)[0]
    mesh_ez_vertical        np.meshgrid(grid_east, grid_z)[1]
    mesh_north              np.meshgrid(grid_east, grid_north)[1]
    mesh_nz_north           np.meshgrid(grid_north, grid_z)[0]
    mesh_nz_vertical        np.meshgrid(grid_north, grid_z)[1]
    model_fn                full path to model file
    ms                      size of station markers in points. *default* is 2
    nodes_east              relative distance betwen nodes in e-w direction
                            in map_scale units
    nodes_north             relative distance betwen nodes in n-s direction
                            in map_scale units
    nodes_z                 relative distance betwen nodes in z direction
                            in map_scale units
    north_line_xlist        list of line nodes north grid for faster plotting
    north_line_ylist        list of line nodes north grid for faster plotting
    ns_limits               (min, max) limits of plots in n-s direction
                            *default* is None, set veiwing area to station area
    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
                            *default* is 'y'
    plot_stations           default False
    plot_grid               show grid on exported plot; default False
    res_model               np.ndarray(n_north, n_east, n_vertical) of
                            model resistivity values in linear scale
    save_format             exported format; default png
    save_path               path to save exported plots to; default current working folder
    station_color           color of station marker. *default* is black
    station_dict_east       location of stations for each east grid row
    station_dict_north      location of stations for each north grid row
    station_east            location of stations in east direction
    station_fn              full path to station file
    station_font_color      color of station label
    station_font_pad        padding between station marker and label
    station_font_rotation   angle of station label
    station_font_size       font size of station label
    station_font_weight     weight of font for station label
    station_id              [min, max] index values for station labels
    station_marker          station marker
    station_names           name of stations
    station_north           location of stations in north direction
    subplot_bottom          distance between axes and bottom of figure window
    subplot_hspace          distance between subplots in vertical direction
    subplot_left            distance between axes and left of figure window
    subplot_right           distance between axes and right of figure window
    subplot_top             distance between axes and top of figure window
    subplot_wspace          distance between subplots in horizontal direction
    title                   title of plot
    xminorticks             location of xminorticks
    yminorticks             location of yminorticks
    z_limits                (min, max) limits in vertical direction,
    ======================= ===================================================

    """

    def __init__(self, model_fn, data_fn=None, **kwargs):
        self.model_fn = model_fn
        self.data_fn = data_fn

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.font_size = kwargs.pop('font_size', 4)

        self.subplot_wspace = .20
        self.subplot_hspace = .30
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.index_vertical = kwargs.pop('index_vertical', 0)
        self.index_east = kwargs.pop('index_east', 0)
        self.index_north = kwargs.pop('index_north', 0)

        self.cmap = kwargs.pop('cmap', 'jet_r')
        self.climits = kwargs.pop('climits', (0, 4))

        self.map_scale = kwargs.pop('map_scale', 'km')
        # make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
        elif self.map_scale == 'm':
            self.dscale = 1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)
        self.z_limits = kwargs.pop('z_limits', None)

        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None
        self.model_epsg = kwargs.pop('model_epsg',None)

        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        self.mesh_east = None
        self.mesh_north = None

        self.station_east = None
        self.station_north = None
        self.station_names = None

        self.station_id = kwargs.pop('station_id', None)
        self.station_font_size = kwargs.pop('station_font_size', 4)
        self.station_font_pad = kwargs.pop('station_font_pad', 1.0)
        self.station_font_weight = kwargs.pop('station_font_weight', 'bold')
        self.station_font_rotation = kwargs.pop('station_font_rotation', 60)
        self.station_font_color = kwargs.pop('station_font_color', 'k')
        self.station_marker = kwargs.pop('station_marker',
                                         r"$\blacktriangledown$")
        self.station_color = kwargs.pop('station_color', 'k')
        self.ms = kwargs.pop('ms', 10)

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.plot_stations = kwargs.pop('plot_stations', False)
        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)
        self.plot_grid = kwargs.pop('plot_grid', False)
        self.draw_colorbar = kwargs.pop('draw_colorbar', True)
        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.save_format = kwargs.pop('save_format', 'png')
        

        # read data
        self.read_files()
        self.get_station_grid_locations()


        # set up kd-tree for interpolation on to arbitrary surfaces
        # intersecting the model
        self._initialize_interpolation()
        
        
        if self.plot_yn == 'y': 
            self.plot()

    def _initialize_interpolation(self):
        self._mx = np.array(self.grid_east)
        self._my = np.array(self.grid_north)
        self._mz = np.array(self.grid_z)

        # Compute cell-centre coordinates
        self._mcx = (self._mx[1:] + self._mx[:-1]) / 2.
        self._mcy = (self._my[1:] + self._my[:-1]) / 2.
        self._mcz = (self._mz[1:] + self._mz[:-1]) / 2.

        # Create mesh-grid based on cell-centre coordinates
        self._mgx, self._mgy, self._mgz = np.meshgrid(self._mcx,
                                                      self._mcy,
                                                      self._mcz)

        # List of xyz coodinates of mesh-grid
        self._mgxyz = np.vstack([self._mgx.flatten(),
                                 self._mgy.flatten(),
                                 self._mgz.flatten()]).T

        # Create Kd-Tree based on mesh-grid coordinates
        self._tree = cKDTree(self._mgxyz)
    # end func

    def get_slice(self, option='STA', coords=[], nsteps=-1, nn=1, p=4,
                  absolute_query_locations = False,
                  extrapolate=True):
        """

        :param option: can be either of 'STA', 'XY' or 'XYZ'. For 'STA' or 'XY', a vertical
                       profile is returned based on station coordinates or arbitrary XY
                       coordinates, respectively. For 'XYZ', interpolated values at those
                       coordinates are returned
        :param coords: Numpy array of shape (np, 2) for option='XY' or of shape (np, 3) for
                       option='XYZ', where np is the number of coordinates. Not used for
                       option='STA', in which case a vertical profile is created based on
                       station locations.
        :param nsteps: When option is set to 'STA' or 'XY', nsteps can be used to create a
                       regular grid along the profile in the horizontal direction. By default,
                       when nsteps=-1, the horizontal grid points are defined by station
                       locations or values in the XY array. This parameter is ignored for
                       option='XYZ'
        :param nn: Number of neighbours to use for interpolation.
                   Nearest neighbour interpolation is returned when nn=1 (default).
                   When nn>1, inverse distance weighted interpolation is returned. See
                   link below for more details:
                   https://en.wikipedia.org/wiki/Inverse_distance_weighting
        :param p: Power parameter, which determines the relative influence of near and far
                  neighbours during interpolation. For p<=3, causes interpolated values to
                  be dominated by points far away. Larger values of p assign greater influence
                  to values near the interpolated point.
        :param absolute_query_locations: if True, query locations are shifted to be centered
               on the center of station locations; default False, in which case the function
               treats query locations as relative coordinates. For option='STA', this parameter
               is ignored, since station locations are internally treated as relative
               coordinates
        :param extrapolate: Extrapolates values (default), which can be particularly useful
                            for extracting values at nodes, since the field values are given
                            for cell-centres.
        :return: 1: when option is 'STA' or 'XY'
                    gd, gz, gv : where gd, gz and gv are 2D grids of distance (along profile),
                    depth and interpolated values, respectively. The shape of the 2D grids
                    depend on the number of stations or the number of xy coordinates provided,
                    for options 'STA' or 'XY', respectively, the number of vertical model grid
                    points and whether regular gridding in the horizontal direction was enabled
                    with nsteps>-1.
                 2: when option is 'XYZ'
                    gv : list of interpolated values of shape (np)
        """

        def distance(P1, P2):
            """
            Compute Euclidean distance
            """

            return ((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2) ** 0.5
        # end func

        def optimized_path(coords, start=None):
            """
            This function was adopted verbatim from the following link:

            https://stackoverflow.com/questions/45829155/sort-points-in-order-to-have-a-continuous-curve-using-python

            Order coordinates such that a continuous line can be formed.
            coords: a list containing coordnates = [ [x1, y1], [x2, y2] , ...]

            """
            if start is None:
                start = coords[0]
            pass_by = coords
            path = [start]
            pass_by.remove(start)
            while pass_by:
                nearest = min(pass_by, key=lambda x: distance(path[-1], x))
                path.append(nearest)
                pass_by.remove(nearest)
            return path
        # end func

        assert option in ['STA', 'XY', 'XYZ'], 'Invalid option; Aborting..'
        if(option == 'STA'):
            if(self.md_data is None):
                print('Station coordinates not available. Aborting..')
                exit(-1)
        elif(option == 'XY'):
            assert type(coords)==np.ndarray and coords.ndim==2 and coords.shape[1]==2, \
                'Shape of coords should be (np, 2); Aborting..'
        elif(option == 'XYZ'):
            assert type(coords)==np.ndarray and coords.ndim==2 and coords.shape[1]==3, \
                'Shape of coords should be (np, 3); Aborting..'

        xyz_list = []
        d = None
        x = None
        y = None
        xmin = 0
        ymin = 0
        if(option == 'STA' or option == 'XY'):
            if(nsteps > -1): assert nsteps > 2, 'Must have more than 2 grid points in the ' \
                                              'horizontal direction. Aborting..'

            x = None
            y = None
            d = None
            if(option == 'STA'):
                x = np.array(self.station_east)
                y = np.array(self.station_north)

                if(nsteps==-1): nsteps = len(x)
            elif(option == 'XY'):
                x = np.array(coords[:,0])
                y = np.array(coords[:,1])

                if(nsteps==-1): nsteps = len(x)
            # end if

            xy = [[a, b] for a, b in zip(x,y)]
            ordered_xy = np.array(optimized_path(xy))
            xx = ordered_xy[:, 0]
            yy = ordered_xy[:, 1]

            dx = xx[:-1] - xx[1:]
            dy = yy[:-1] - yy[1:]
            dst = np.cumsum(np.sqrt(dx ** 2 + dy ** 2))
            dst = np.insert(dst, 0, 0)

            xio = interp1d(dst, xx)
            yio = interp1d(dst, yy)
            # ====

            if(nsteps>-1):
                d = np.linspace(dst.min(), dst.max(), nsteps) # create regular grid
            for zi in self.grid_z:
                for xi,yi in zip(xio(d), yio(d)):
                    xyz_list.append([xi+xmin, yi+ymin, zi])
            xyz_list = np.array(xyz_list)
        elif(option == 'XYZ'):
            xyz_list = coords
        # end if

        gv = self._get_slice_helper(xyz_list, nn, p, absolute_query_locations, extrapolate)

        if(option=='STA' or option=='XY'):
            gz, gd = np.meshgrid(self.grid_z, d, indexing='ij')
            gv = gv.reshape(gd.shape)
            return gd, gz, gv
        elif(option=='XYZ'):
            return gv

        return None
    # end func

    def _get_slice_helper(self, _xyz_list, nn=1, p=4, absolute_query_locations=False,
                          extrapolate=True):
        '''
        Function to retrieve interpolated field values at arbitrary locations

        :param xyz_list: numpy array of shape (np,3), where np in the number of points
        :param nn: as above
        :param p: as above
        :param absolute_query_locations: as above
        :param extrapolate: as above
        :return: numpy array of interpolated values of shape (np)
        '''

        xyz_list = np.array(_xyz_list)
        if(absolute_query_locations):
            if(self.md_data is None):
                print('Station coordinates not available. Aborting..')
                exit(-1)
            #end if

            xyz_list[:, 0] -= self.md_data.center_point['east']
            xyz_list[:, 1] -= self.md_data.center_point['north']
        # end if

        # query Kd-tree instance to retrieve distances and
        # indices of k nearest neighbours
        d, l = self._tree.query(xyz_list, k=nn)

        img = None
        if (nn == 1):
            # extract nearest neighbour values
            img = self.res_model.flatten()[l]
        else:
            vals = self.res_model.flatten()
            img = np.zeros((xyz_list.shape[0]))

            # field values are directly assigned for coincident locations
            coincidentValIndices = d[:, 0] == 0
            img[coincidentValIndices] = vals[l[coincidentValIndices, 0]]

            # perform idw interpolation for non-coincident locations
            idwIndices = d[:, 0] != 0
            w = np.zeros(d.shape)
            w[idwIndices, :] = 1. / np.power(d[idwIndices, :], p)

            img[idwIndices] = np.sum(w[idwIndices, :] * vals[l[idwIndices, :]], axis=1) / \
                              np.sum(w[idwIndices, :], axis=1)
        # end if

        if (extrapolate == False):
            # if extrapolate is false, set interpolation values to NaN for locations
            # outside the model domain
            minX = np.min(self._mgxyz[:, 0])
            maxX = np.max(self._mgxyz[:, 0])

            minY = np.min(self._mgxyz[:, 1])
            maxY = np.max(self._mgxyz[:, 1])

            minZ = np.min(self._mgxyz[:, 2])
            maxZ = np.max(self._mgxyz[:, 2])

            xFilter = np.array(xyz_list[:, 0] < minX) + \
                      np.array(xyz_list[:, 0] > maxX)
            yFilter = np.array(xyz_list[:, 1] < minY) + \
                      np.array(xyz_list[:, 1] > maxY)
            zFilter = np.array(xyz_list[:, 2] < minZ) + \
                      np.array(xyz_list[:, 2] > maxZ)

            img[xFilter] = np.nan
            img[yFilter] = np.nan
            img[zFilter] = np.nan
        # end if

        return img
    # end func

    def read_files(self):
        """
        read in the files to get appropriate information
        """
        # --> read in model file
        if self.model_fn is not None:
            if os.path.isfile(self.model_fn) == True:
                md_model = Model()
                md_model.read_model_file(self.model_fn)
                self.res_model = md_model.res_model
                self.grid_east = md_model.grid_east / self.dscale
                self.grid_north = md_model.grid_north / self.dscale
                self.grid_z = md_model.grid_z / self.dscale
                self.nodes_east = md_model.nodes_east / self.dscale
                self.nodes_north = md_model.nodes_north / self.dscale
                self.nodes_z = md_model.nodes_z / self.dscale

                self.md_model = md_model
            else:
                raise mtex.MTpyError_file_handling(
                    '{0} does not exist, check path'.format(self.model_fn))

        # --> read in data file to get station locations
        if self.data_fn is not None:
            if os.path.isfile(self.data_fn) == True:
                md_data = Data(model_epsg=self.model_epsg)
                md_data.read_data_file(self.data_fn)
                
                self.station_east = md_data.station_locations.rel_east / self.dscale
                self.station_north = md_data.station_locations.rel_north / self.dscale
                self.station_names = md_data.station_locations.station
                self.station_elev = md_data.station_locations.elev / self.dscale

                self.md_data = md_data
            else:
                print('Could not find data file {0}'.format(self.data_fn))



                
                
    def basemap_plot(self, depth, basemap = None,tick_interval=None, save=False, 
                     save_path=None, new_figure=True,mesh_rotation_angle=0.,
                     overlay=False,clip=[0,0],**basemap_kwargs):
        """
        plot model depth slice on a basemap using basemap modules in matplotlib
        
        :param depth: depth in model to plot
        :param tick_interval: tick interval on map in degrees, if None it is 
                              calculated from the data extent
        :param save: True/False, whether or not to save and close figure
        :param savepath: full path of file to save to, if None, saves to 
                         self.save_path
        :new_figure: True/False, whether or not to initiate a new figure for
                     the plot
        :param mesh_rotation_angle: rotation angle of mesh, in degrees 
                                    clockwise from north
        :param **basemap_kwargs: provide any valid arguments to Basemap
                                 instance and these will be passed to the map.
        
        """
        
        if self.model_epsg is None:
            print("No projection information provided, please provide the model epsg code relevant to your model")
            return
        
        # initialise plot parameters
        mpldict={}
        mpldict['cmap'] = cm.get_cmap(self.cmap)
        mpldict['norm'] = colors.LogNorm()
        mpldict['vmin'] = 10**self.climits[0]
        mpldict['vmax'] = 10**self.climits[1]
        
        # find nearest depth index
        depthIdx = nearest_index(depth,self._mcz*self.dscale)
        
        
        
        if new_figure:
            plt.figure()

        
        # get eastings/northings of mesh
        ge,gn = self.md_model.grid_east, self.md_model.grid_north
        e0,n0 = self.md_data.center_point['east'],self.md_data.center_point['north']
        
        if mesh_rotation_angle != 0:
            if hasattr(self,'mesh_rotation_angle'):
                angle_to_rotate_stations = self.mesh_rotation_angle - mesh_rotation_angle
            else:
                angle_to_rotate_stations = -mesh_rotation_angle
                
            self.mesh_rotation_angle = mesh_rotation_angle


            mgx, mgy = rotate_mesh(ge,gn,[e0,n0],
                                   -mesh_rotation_angle)
        else:
            mgx, mgy = np.meshgrid(ge + e0, gn + n0)

        # rotate stations if necessary
        if mesh_rotation_angle != 0:
            self.md_data.station_locations.rotate_stations(angle_to_rotate_stations)
                        
            # get relative locations
            seast,snorth = self.md_data.station_locations.rel_east + self.md_data.station_locations.center_point['east'],\
                           self.md_data.station_locations.rel_north + self.md_data.station_locations.center_point['north']
            
            # project station location eastings and northings to lat/long
            slon,slat = epsg_project(seast,snorth,self.model_epsg,4326)
            self.md_data.station_locations.station_locations['lon'] = slon
            self.md_data.station_locations.station_locations['lat'] = slat    
        
        
        if basemap is None:
            # initialise a basemap with extents, projection etc calculated from data 
            # if not provided in basemap_kwargs
            self.bm = basemap_tools.initialise_basemap(self.md_data.station_locations,**basemap_kwargs)
            # add frame to basemap and plot data
            basemap_tools.add_basemap_frame(self.bm,tick_interval=tick_interval)
        else:
            self.bm = basemap

        
        # lat/lon coordinates of resistivity model values
        loncg,latcg = epsg_project(mgx,mgy,
                                   self.model_epsg,
                                   4326)
               
        # get x and y projected positions on the basemap                    
        xcg,ycg = self.bm(loncg,latcg)
        
        # get clip extents
        rx0,rx1 = clip[0],xcg.shape[1]-clip[0]
        ry0,ry1 = clip[1],ycg.shape[0]-clip[1]
        
        # plot model on basemap, applying clip
        basemap_tools.plot_data(xcg[ry0:ry1,rx0:rx1],
                                ycg[ry0:ry1,rx0:rx1],
                                self.res_model[ry0:ry1,rx0:rx1,depthIdx],
                                basemap=self.bm,
                                **mpldict)
                                   
        # plot stations
        if self.plot_stations:
            # rotate stations
            seast,snorth = self.md_data.station_locations.rel_east + self.md_data.center_point['east'],\
                           self.md_data.station_locations.rel_north + self.md_data.center_point['north']

            # reproject station location eastings and northings
            slon,slat = epsg_project(seast,snorth,self.model_epsg,4326)
            
            sx,sy = self.bm(slon,slat)
            self.bm.plot(sx,sy,'k.')
            
            
        # draw colorbar
        if self.draw_colorbar:
            plt.colorbar(shrink=0.5)
            
        if save:
            if save_path is not None:
                self.save_path = save_path
            plt.savefig(os.path.join(self.save_path,'DepthSlice%1i%1s.png'%(depth/self.dscale,self.map_scale)),
                        dpi=self.fig_dpi)
            plt.close()
    
        

    def plot(self):
        """
        plot:
            east vs. vertical,
            north vs. vertical,
            east vs. north


        """
        print("=============== ===============================================")
        print("    Buttons                  Description                       ")
        print("=============== ===============================================")
        print("     'e'          moves n-s slice east by one model block")
        print("     'w'          moves n-s slice west by one model block")
        print("     'n'          moves e-w slice north by one model block")
        print("     'm'          moves e-w slice south by one model block")
        print("     'd'          moves depth slice down by one model block")
        print("     'u'          moves depth slice up by one model block")
        print("=============== ===============================================")

        self.font_dict = {'size': self.font_size*0.75, 'weight': 'bold'}

        # --> set default font size
        plt.rcParams['font.size'] = self.font_size*0.75
        plt.rcParams['xtick.major.pad'] = '1'
        plt.rcParams['ytick.major.pad'] = '1'
        plt.rcParams['ytick.major.pad'] = '1'

        # set the limits of the plot
        if self.ew_limits == None:
            if self.station_east is not None:
                self.ew_limits = (np.floor(self.station_east.min()),
                                  np.ceil(self.station_east.max()))
            else:
                self.ew_limits = (self.grid_east[5], self.grid_east[-5])

        if self.ns_limits == None:
            if self.station_north is not None:
                self.ns_limits = (np.floor(self.station_north.min()),
                                  np.ceil(self.station_north.max()))
            else:
                self.ns_limits = (self.grid_north[5], self.grid_north[-5])

        if self.z_limits == None:
            depth_limit = max([(abs(self.ew_limits[0]) + abs(self.ew_limits[1])),
                               (abs(self.ns_limits[0]) + abs(self.ns_limits[1]))])
            self.z_limits = (-5000 / self.dscale, depth_limit)


        self.fig = plt.figure(self.fig_num, figsize=self.fig_size,
                              dpi=self.fig_dpi,frameon=False)
        plt.clf()

        # annotations
        self.ax_border = plt.axes([0.01, 0.01, 0.98, 0.3])
        self.ax_border.set_xticks([])
        self.ax_border.set_yticks([])
        self.ax_border.set_title('Select/Export Slices',
                                 y=-0.01, fontdict={'size': self.font_size*2,
                                                     'weight': 'bold'})

        # set up plot axes
        gs = gridspec.GridSpec(3, 2,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom,
                               right=self.subplot_right,
                               hspace=self.subplot_hspace)

        # make subplots
        self.ax_ez = self.fig.add_subplot(gs[0, 0], aspect=self.fig_aspect)
        self.ax_nz = self.fig.add_subplot(gs[1, 1], aspect=self.fig_aspect)
        self.ax_en = self.fig.add_subplot(gs[1, 0], aspect=self.fig_aspect)
        self.ax_map = self.fig.add_subplot(gs[0, 1])
        self.ax_radio = plt.axes([0.1, 0.05, 0.1, 0.2])
        self.ax_span =  plt.axes([0.3, 0.15, 0.6, 0.1])
        self.ax_button = plt.axes([0.57, 0.075, 0.06, 0.03])

        # set tick sizes
        axList = [self.ax_ez, self.ax_nz, self.ax_en, self.ax_map]
        for ax in axList: ax.tick_params(axis='both', length=2)

        # make grid meshes being sure the indexing is correct
        self.mesh_ez_east, self.mesh_ez_vertical = np.meshgrid(self.grid_east,
                                                               self.grid_z,
                                                               indexing='ij')
        self.mesh_nz_north, self.mesh_nz_vertical = np.meshgrid(self.grid_north,
                                                                self.grid_z,
                                                                indexing='ij')
        self.mesh_en_east, self.mesh_en_north = np.meshgrid(self.grid_east,
                                                            self.grid_north,
                                                            indexing='ij')

        # --> plot east vs vertical
        self._update_ax_ez()

        # --> plot north vs vertical
        self._update_ax_nz()

        # --> plot east vs north
        self._update_ax_en()

        # --> plot the grid as a map view
        self._update_map()

        # plot color bar
        cbx = mcb.make_axes(self.ax_map, fraction=.15, shrink=.75, pad=.15)

        if type(self.cmap) == str:
            self.cmap=cm.get_cmap(self.cmap)
        cb = mcb.ColorbarBase(cbx[0],
                              cmap=self.cmap,
                              norm=Normalize(vmin=self.climits[0],
                                             vmax=self.climits[1]))

        cb.ax.yaxis.set_label_position('right')
        cb.ax.yaxis.set_label_coords(1.25, .5)
        cb.ax.yaxis.tick_left()
        cb.ax.tick_params(axis='y', direction='in')

        cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size': self.font_size}, x=2)

        cb.set_ticks(np.arange(np.ceil(self.climits[0]),
                               np.floor(self.climits[1] + 1)))
        cblabeldict = {-2: '$10^{-3}$', -1: '$10^{-1}$', 0: '$10^{0}$', 1: '$10^{1}$',
                       2: '$10^{2}$', 3: '$10^{3}$', 4: '$10^{4}$', 5: '$10^{5}$',
                       6: '$10^{6}$', 7: '$10^{7}$', 8: '$10^{8}$'}
        cb.set_ticklabels([cblabeldict[cc]
                           for cc in np.arange(np.ceil(self.climits[0]),
                                               np.floor(self.climits[1] + 1))])

        self.key_press = self.fig.canvas.mpl_connect('key_press_event',
                                                     self.on_key_press)

        # Interactive widgets ==========================================
        def getCursorValue():
            if(self.current_label == 'N-E'):
                return self.grid_z[self.index_vertical]
            elif(self.current_label == 'N-Z'):
                return self.grid_east[self.index_east]
            elif(self.current_label == 'E-Z'):
                return self.grid_north[self.index_north]
        # end func

        self.current_range = self.z_limits
        self.current_label = 'N-E'
        self.current_label_desc = {'N-E': 'Depth',
                                   'N-Z': 'Easting',
                                   'E-Z': 'Northing'}
        self.axis_values = {'N-E':self.grid_z,
                            'N-Z':self.grid_east,
                            'E-Z':self.grid_north}
        self.axis_cursor_colors = {'N-E':'r',
                                   'N-Z':'b',
                                   'E-Z':'g'}
        self.selected_indices = []

        self.ax_span.scatter(self.axis_values[self.current_label],
                             np.ones(self.axis_values[self.current_label].shape[0]) *
                             (self.current_range[0] + self.current_range[1]) / 2.,
                             0.5, zorder=100, marker='o', color='k')

        self.ax_span.fill_between(self.current_range,
                                  self.current_range[0] * np.ones(len(self.current_range)),
                                  self.current_range[1] * np.ones(len(self.current_range)),
                                  alpha=0.5, facecolor='b')
        self.ax_span.plot(np.ones(2)*getCursorValue(),
                          np.array(self.current_range),
                          c=self.axis_cursor_colors[self.current_label], lw=1)

        self.ax_span.set_xlim(self.current_range)
        self.ax_span.set_ylim(self.current_range)
        self.ax_span.set_yticks([])
        self.ax_span.set_aspect(0.05)

        self.ax_span.set_title('Depth Extent: Click+Drag to Select Sub-range')
        def updateRange(label):
            self.current_label = label
            if(label == 'N-E'):
                self.current_range = self.z_limits
            elif(label == 'N-Z'):
                self.current_range = self.ew_limits
            else:
                self.current_range = self.ns_limits

            self.ax_span.cla()

            self.ax_span.scatter(self.axis_values[self.current_label],
                                 np.ones(self.axis_values[self.current_label].shape[0])*
                                 (self.current_range[0]+self.current_range[1])/2.,
                                 0.5, zorder=100, marker='o', color='k')

            self.ax_span.fill_between(self.current_range,
                                      self.current_range[0] * np.ones(len(self.current_range)),
                                      self.current_range[1] * np.ones(len(self.current_range)),
                                      alpha=0.5, facecolor='b')
            self.ax_span.plot(np.ones(2) * getCursorValue(),
                              np.array(self.current_range),
                              c=self.axis_cursor_colors[self.current_label], lw=1)

            self.ax_span.set_yticks([])
            self.ax_span.set_title('%s Extent: Click+Drag to Select Sub-range'%
                                   (self.current_label_desc[label]))
            self.ax_span.set_xlim(self.current_range)
            self.ax_span.set_ylim(self.current_range)
            self.ax_span.set_aspect(0.05)
            self.fig.canvas.draw_idle()

            self.selected_indices = []
        # end func

        def onSelect(xmin, xmax):
            updateRange(self.current_label)
            indmin, indmax = np.searchsorted(self.axis_values[self.current_label], (xmin, xmax))

            self.ax_span.fill_between(np.linspace(xmin, xmax, 100),
                                      self.current_range[0] * np.ones(100),
                                      self.current_range[1] * np.ones(100),
                                      facecolor='red',
                                      alpha=0.4,
                                      edgecolor='none')

            self.selected_indices = np.arange(indmin, indmax)
            print('Selected indices: ' + str(self.selected_indices))

            self.ax_span.set_yticks([])
            self.ax_span.set_title('%s Extent: Click+Drag to Select Sub-range'%
                                   (self.current_label_desc[self.current_label]))
            self.ax_span.set_xlim(self.current_range)
            self.ax_span.set_ylim(self.current_range)
            self.ax_span.set_aspect(0.05)
            self.fig.canvas.draw_idle()
        #end func

        def buttonClicked(event):
            self.export_slices(self.current_label, self.selected_indices)
        # end func

        radio = RadioButtons(self.ax_radio, ('N-E', # (Depth Slice)
                                   'N-Z', # (North-south-aligned vertical profile)
                                   'E-Z'), #(East-west-aligned vertical profile)
                             active=0)
        self.ax_radio.set_title('Plane')

        radio.on_clicked(updateRange)

        span = SpanSelector(self.ax_span, onSelect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red', edgecolor='none'))

        button = Button(self.ax_button, 'Export', color='lightgoldenrodyellow',
                        hovercolor='orange')
        button.on_clicked(buttonClicked)
        self.update_range_func = updateRange

        # Only show interactive popout if plot_yn is set to True; otherwise hide
        # popout
        if(self.plot_yn == 'y'):
            plt.show()
        else:
            self.fig.set_visible(False)
            plt.close()
#            plt.draw()
    # end func

    def export_slices(self, plane='N-E', indexlist=[], station_buffer=200, save=True):
        """
        Plot Slices

        :param plane: must be either 'N-E', 'N-Z' or 'E-Z'
        :param indexlist: must be a list or 1d numpy array of indices
        :param station_buffer: spatial buffer (in metres) used around grid locations for
                               selecting stations to be projected and plotted on profiles.
                               Ignored if .plot_stations is set to False.
        :return: [figlist, savepaths]. A list containing (1) lists of Figure objects,
                 for further manipulation (2) corresponding paths for saving them to disk
        """

        station_buffer /= self.dscale
        assert plane in ['N-E', 'N-Z', 'E-Z'], 'Invalid plane; Aborting..'
        assert type(indexlist) == list or type(indexlist) == np.ndarray, \
                'Index list must be of type list or a 1d numpy array. Aborting..'

        fdict = {'size': self.font_size, 'weight': 'bold'}

        cblabeldict = {-2: '$10^{-3}$', -1: '$10^{-1}$', 0: '$10^{0}$', 1: '$10^{1}$',
                       2: '$10^{2}$', 3: '$10^{3}$', 4: '$10^{4}$', 5: '$10^{5}$',
                       6: '$10^{6}$', 7: '$10^{7}$', 8: '$10^{8}$'}

        # make a mesh grid of nodes
        xg, yg = None, None
        if(plane == 'N-E'):
            xg, yg = self.mesh_en_east, self.mesh_en_north
        elif(plane == 'N-Z'):
            xg, yg = self.mesh_nz_north, self.mesh_nz_vertical
        elif(plane == 'E-Z'):
            xg, yg = self.mesh_ez_east, self.mesh_ez_vertical

        plt.rcParams['font.size'] = self.font_size

        figlist = []
        fnlist = []
        # --> plot slices into individual figures
        for ii in indexlist:
            #depth = '{0:.3f} ({1})'.format(self.grid_z[ii],
            #                               self.map_scale)

            fig = plt.figure(figsize=self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            ax1 = fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)
            ax1.tick_params(axis='both', length=2)

            if (plane == 'N-E'):
                plot_res = np.log10(self.res_model[:, :, ii].T)
                ax1.set_xlim(self.ew_limits)
                ax1.set_ylim(self.ns_limits)
                ax1.set_ylabel('Northing (' + self.map_scale + ')', fontdict=fdict)
                ax1.set_xlabel('Easting (' + self.map_scale + ')', fontdict=fdict)
            elif (plane == 'N-Z'):
                plot_res = np.log10(self.res_model[:, ii, :])
                ax1.set_xlim(self.ns_limits)
                ax1.set_ylim(self.z_limits)
                ax1.invert_yaxis()
                ax1.set_ylabel('Depth (' + self.map_scale + ')', fontdict=fdict)
                ax1.set_xlabel('Northing (' + self.map_scale + ')', fontdict=fdict)
            elif (plane == 'E-Z'):
                plot_res = np.log10(self.res_model[ii, :, :])
                ax1.set_xlim(self.ew_limits)
                ax1.set_ylim(self.z_limits)
                ax1.invert_yaxis()
                ax1.set_ylabel('Depth (' + self.map_scale + ')', fontdict=fdict)
                ax1.set_xlabel('Easting (' + self.map_scale + ')', fontdict=fdict)
            # end if

            mesh_plot = ax1.pcolormesh(xg,
                                       yg,
                                       plot_res,
                                       cmap=self.cmap,
                                       vmin=self.climits[0],
                                       vmax=self.climits[1])
            # plot the stations
            if (self.station_east is not None \
                    and self.plot_stations):

                if(plane == 'N-E'):
                    for ee, nn, slabel in zip(self.station_east, self.station_north, self.station_names):
                        if self.station_id is not None:
                            slabel = slabel[self.station_id[0]:self.station_id[1]]
                        # plot marker
#                        ax1.plot(ee, nn, 'k.')
                        ax1.text(ee, nn, '*', 
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 fontdict={'size': 3, 'weight': 'bold'})
                        # plot label
                        ax1.text(ee, nn, slabel,
                                 rotation=45,
                                 verticalalignment='bottom',
                                 horizontalalignment='left',
                                 fontdict={'size': 3, 'weight': 'bold'})
                elif(plane == 'N-Z'):
                    sids = np.fabs(self.grid_east[ii] - self.station_east) < station_buffer
                    nvals = self.station_north[sids]
                    for x in nvals:
                        ax1.text(x,
                                 0,
                                 self.station_marker,
                                 horizontalalignment='center',
                                 verticalalignment='baseline',
                                 fontdict={'size': self.ms,
                                 'color': self.station_color})
                elif (plane == 'E-Z'):
                    sids = np.fabs(self.grid_north[ii] - self.station_north) < station_buffer
                    evals = self.station_east[sids]
                    for x in evals:
                        ax1.text(x,
                                 0,
                                 self.station_marker,
                                 horizontalalignment='center',
                                 verticalalignment='baseline',
                                 fontdict={'size': self.ms,
                                 'color': self.station_color})
                # end if

                # plot the grid if desired
            if self.plot_grid == 'y':
                x_line_xlist = []
                x_line_ylist = []
                for xx in xg[:,0]:
                    x_line_xlist.extend([xx, xx])
                    x_line_xlist.append(None)
                    x_line_ylist.extend([yg[0,:].min(),
                                         yg[0,:].max()])
                    x_line_ylist.append(None)
                ax1.plot(x_line_xlist,
                         x_line_ylist,
                         lw=.25,
                         color='k')

                y_line_xlist = []
                y_line_ylist = []
                for yy in yg[0,:]:
                    y_line_xlist.extend([xg[:,0].min(),
                                         xg[:,0].max()])
                    y_line_xlist.append(None)
                    y_line_ylist.extend([yy, yy])
                    y_line_ylist.append(None)
                ax1.plot(y_line_xlist,
                         y_line_ylist,
                         lw=.25,
                         color='k')

            # plot the colorbar
            if self.draw_colorbar:
                
                cbx = mcb.make_axes(ax1, fraction=.15, shrink=.75, pad=.15)
                cb = mcb.ColorbarBase(cbx[0],
                                      cmap=self.cmap,
                                      norm=Normalize(vmin=self.climits[0],
                                                     vmax=self.climits[1]))

                cb.ax.yaxis.set_label_position('right')
                cb.ax.yaxis.set_label_coords(1.25, .5)
                cb.ax.yaxis.tick_left()
                cb.ax.tick_params(axis='y', direction='in')

                cb.set_label('Resistivity ($\Omega \cdot$m)',
                             fontdict={'size': self.font_size}, x=2)

                cb.set_ticks(np.arange(np.ceil(self.climits[0]),
                                       np.floor(self.climits[1] + 1)))
                cblabeldict = {-2: '$10^{-3}$', -1: '$10^{-1}$', 0: '$10^{0}$', 1: '$10^{1}$',
                               2: '$10^{2}$', 3: '$10^{3}$', 4: '$10^{4}$', 5: '$10^{5}$',
                               6: '$10^{6}$', 7: '$10^{7}$', 8: '$10^{8}$'}
                cb.set_ticklabels([cblabeldict[cc]
                                   for cc in np.arange(np.ceil(self.climits[0]),
                                                       np.floor(self.climits[1] + 1))])
            # end if

            #plt.show()
            figlist.append(fig)

            if self.title == 'on':
                fig.suptitle('%s Plane at %s: %0.4f %s'%(plane,
                                               self.current_label_desc[plane],
                                               self.axis_values[plane][ii],
                                               self.map_scale))
            if save:
                # --> save plots to a common folder
                fn = '%s-plane-at-%s%03i.%0.3f.%s.%s'%(plane,
                                               self.current_label_desc[plane],
                                               ii,
                                               self.axis_values[plane][ii],
                                               self.map_scale,
                                               self.save_format)
    

                fpath = os.path.join(self.save_path, fn)
                print(('Exporting %s..'%(fpath)))
                fig.savefig(fpath, dpi=self.fig_dpi, bbox_inches='tight')
                fnlist.append(fpath)
    
                #fig.clear()
#                plt.close()
        # end for
        return figlist, fnlist
    #end func

    def on_key_press(self, event):
        """
        on a key press change the slices

        """

        key_press = event.key

        if key_press == 'n':
            if self.index_north == self.grid_north.size:
                print('Already at northern most grid cell')
            else:
                self.index_north += 1
                if self.index_north > self.grid_north.size:
                    self.index_north = self.grid_north.size
            self._update_ax_ez()
            self._update_map()

        if key_press == 'm':
            if self.index_north == 0:
                print('Already at southern most grid cell')
            else:
                self.index_north -= 1
                if self.index_north < 0:
                    self.index_north = 0
            self._update_ax_ez()
            self._update_map()

        if key_press == 'e':
            if self.index_east == self.grid_east.size:
                print('Already at eastern most grid cell')
            else:
                self.index_east += 1
                if self.index_east > self.grid_east.size:
                    self.index_east = self.grid_east.size
            self._update_ax_nz()
            self._update_map()

        if key_press == 'w':
            if self.index_east == 0:
                print('Already at western most grid cell')
            else:
                self.index_east -= 1
                if self.index_east < 0:
                    self.index_east = 0
            self._update_ax_nz()
            self._update_map()

        if key_press == 'd':
            if self.index_vertical == self.grid_z.size:
                print('Already at deepest grid cell')
            else:
                self.index_vertical += 1
                if self.index_vertical > self.grid_z.size:
                    self.index_vertical = self.grid_z.size
            self._update_ax_en()
            self._update_ax_nz()
            print('Depth = {0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                                 self.map_scale))

        if key_press == 'u':
            if self.index_vertical == 0:
                print('Already at surface grid cell')
            else:
                self.index_vertical -= 1
                if self.index_vertical < 0:
                    self.index_vertical = 0
            self._update_ax_en()
            self._update_ax_nz()
            print('Depth = {0:.5gf} ({1})'.format(self.grid_z[self.index_vertical],
                                                  self.map_scale))
        self.update_range_func(self.current_label)
    # end func

    def _update_ax_ez(self):
        """
        update east vs vertical plot
        """
        self.ax_ez.cla()
        plot_ez = np.log10(self.res_model[self.index_north, :, :])
        self.ax_ez.pcolormesh(self.mesh_ez_east,
                              self.mesh_ez_vertical,
                              plot_ez,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        # plot stations
        for sx in self.station_dict_north[self.grid_north[self.index_north]]:
            self.ax_ez.text(sx,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size': self.ms,
                                      'color': self.station_color})

        self.ax_ez.set_xlim(self.ew_limits)
        self.ax_ez.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_ez.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_ez.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()

    def _update_ax_nz(self):
        """
        update east vs vertical plot
        """
        self.ax_nz.cla()
        plot_nz = np.log10(self.res_model[:, self.index_east, :])
        self.ax_nz.pcolormesh(self.mesh_nz_north,
                              self.mesh_nz_vertical,
                              plot_nz,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])

        # --> depth indication line
        self.ax_nz.plot([self.grid_north.min(),
                         self.grid_north.max()],
                        [self.grid_z[self.index_vertical],
                         self.grid_z[self.index_vertical]],
                         lw=1,
                         color='r')

        # plot stations
        for sy in self.station_dict_east[self.grid_east[self.index_east]]:
            self.ax_nz.text(sy,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size': self.ms,
                                      'color': self.station_color})
        self.ax_nz.set_xlim(self.ns_limits)
        self.ax_nz.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_nz.set_xlabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_nz.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()

    def _update_ax_en(self):
        """
        update east vs vertical plot
        """

        self.ax_en.cla()
        plot_en = np.log10(self.res_model[:, :, self.index_vertical].T)
        self.ax_en.pcolormesh(self.mesh_en_east,
                              self.mesh_en_north,
                              plot_en,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        self.ax_en.set_xlim(self.ew_limits)
        self.ax_en.set_ylim(self.ns_limits)
        self.ax_en.set_ylabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_en.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        # --> plot the stations
        if self.station_east is not None and self.plot_stations:
            for ee, nn, elev, name in zip(self.station_east,
                                          self.station_north,
                                          self.station_elev,
                                          self.station_names):
                if elev <= self.grid_z[self.index_vertical]:
                    self.ax_en.text(ee, nn, '+',
                                    verticalalignment='center',
                                    horizontalalignment='center',
                                    fontdict={'size': 1, 'weight': 'bold',
                                              'color': (.75, 0, 0)})
                    self.ax_en.text(ee, nn, name[2:],
                                    verticalalignment='center',
                                    horizontalalignment='center',
                                    fontdict={'size': 1, 'weight': 'bold',
                                              'color': (.75, 0, 0)})

        self.fig.canvas.draw()
        self._update_map()

    def _update_map(self):
        self.ax_map.cla()
        self.east_line_xlist = []
        self.east_line_ylist = []
        for xx in self.grid_east:
            self.east_line_xlist.extend([xx, xx])
            self.east_line_xlist.append(None)
            self.east_line_ylist.extend([self.grid_north.min(),
                                         self.grid_north.max()])
            self.east_line_ylist.append(None)
        self.ax_map.plot(self.east_line_xlist,
                         self.east_line_ylist,
                         lw=.25,
                         color='k')

        self.north_line_xlist = []
        self.north_line_ylist = []
        for yy in self.grid_north:
            self.north_line_xlist.extend([self.grid_east.min(),
                                          self.grid_east.max()])
            self.north_line_xlist.append(None)
            self.north_line_ylist.extend([yy, yy])
            self.north_line_ylist.append(None)
        self.ax_map.plot(self.north_line_xlist,
                         self.north_line_ylist,
                         lw=.25,
                         color='k')
        # --> e-w indication line
        self.ax_map.plot([self.grid_east.min(),
                          self.grid_east.max()],
                         [self.grid_north[self.index_north],
                          self.grid_north[self.index_north]],
                         lw=1,
                         color='g')

        # --> e-w indication line
        self.ax_map.plot([self.grid_east[self.index_east],
                          self.grid_east[self.index_east]],
                         [self.grid_north.min(),
                          self.grid_north.max()],
                         lw=1,
                         color='b')
        # --> plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax_map.text(ee, nn, '*',
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 fontdict={'size': 5, 'weight': 'bold'})

        self.ax_map.set_xlim(self.ew_limits)
        self.ax_map.set_ylim(self.ns_limits)
        self.ax_map.set_ylabel('Northing ({0})'.format(self.map_scale),
                               fontdict=self.font_dict)
        self.ax_map.set_xlabel('Easting ({0})'.format(self.map_scale),
                               fontdict=self.font_dict)

        # plot stations
        self.ax_map.text(self.ew_limits[0] * .95, self.ns_limits[1] * .95,
                         '{0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                                self.map_scale),
                         horizontalalignment='left',
                         verticalalignment='top',
                         bbox={'facecolor': 'white'},
                         fontdict=self.font_dict)

        self.fig.canvas.draw()

    def get_station_grid_locations(self):
        """
        get the grid line on which a station resides for plotting

        """
        self.station_dict_east = dict([(gx, []) for gx in self.grid_east])
        self.station_dict_north = dict([(gy, []) for gy in self.grid_north])
        if self.station_east is not None:
            for ss, sx in enumerate(self.station_east):
                gx = np.where(self.grid_east <= sx)[0][-1]
                self.station_dict_east[self.grid_east[gx]].append(self.station_north[ss])

            for ss, sy in enumerate(self.station_north):
                gy = np.where(self.grid_north <= sy)[0][-1]
                self.station_dict_north[self.grid_north[gy]].append(self.station_east[ss])
        else:
            return

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

    def plot_resistivity_on_seismic(self, segy_fn, velocity_model=6000, pick_every=10, ax=None, cb_ax=None, percent_clip=99, alpha=0.5, **kwargs):
        """
        :param segy_fn: SegY file name
        :param velocity_model: can be either the name of a velocity-model file containing stacking velocities
                               for the given 2D seismic line, or a floating point value representing a
                               constant velocity (m/s)
        :param pick_every: this parameter controls the decimation factor for the SegY file; e.g. if pick_every=10,
                           every 10th trace from the SegY file is read in. This significantly speeds up plotting
                           routines.
        :param ax: figure axes
        :param cb_ax: colorbar axes
        :param percent_clip: percentile value used for filtering out seismic amplitudes from plot; e.g. for a value
                             of 99, only seismic amplitudes above the 99th percentile are plotted. The parameter is
                             tuned to plot only the required level of seismic detail.
        :param alpha: alpha value used while resistivity and seismic values
        :param kwargs:

        max_depth : maximum depth extent of plots
        time_shift : time shift in ms to remove topography

        :return: fig, ax : a figure and an plot axes object are returned when the parameter ax is not provided
        """

        # process keyword arguments
        max_depth = kwargs.pop('max_depth', 60e3) # 60 km default
        time_shift = kwargs.pop('time_shift', 225) # 225 ms

        sl = Segy(segy_fn, pick_every=pick_every)
        vm = None

        if(type(velocity_model) == str):
            vm = VelocityModel(velocity_model)
        else:
            try:
                vm = np.float_(velocity_model)
            except:
                raise ValueError('Invalid velocity model')
            # end try
        # end if

        # fetch a depth-migrated image (done using velocity model
        mdepth, mdist, svals, xy_list = sl.getMigratedProfile(vm, max_depth=max_depth, time_shift=time_shift)

        # fetch resistivity along seismic line
        gd, gz, mvals = self.get_slice(option='XY', coords=xy_list, nn=1, absolute_query_locations=True)

        # create plot axes, if not provided
        fig = None
        if(ax is None):
            fig, ax = plt.subplots(1, 1)
            fig.set_size_inches(10,5)

            # if a new figure object is created, user-provided colorbar axes is ignored
            cb_ax = fig.add_axes([0.25, 0.25, 0.5, 0.025])
        # end if

        # plot resistivity ==================================================
        ci = ax.pcolor(gd, gz, mvals,
                       norm=colors.LogNorm(),
                       vmin=np.power(10, 0.5), vmax=np.power(10, 7),
                       cmap='nipy_spectral',
                       alpha=alpha, linewidth=0,
                       edgecolors='None',
                       rasterized=True)

        # deal with white stripes
        ci.set_antialiaseds(True)
        ci.set_rasterized(True)

        # plot colorbar
        if(cb_ax):
            cb = plt.colorbar(ci, cax=cb_ax, ticks=LogLocator(subs=range(10)),
                              orientation="horizontal")
            cb.solids.set_edgecolor('none')
            cb.solids.set_antialiased(True)
            cb.solids.set_rasterized(True)
        # end if

        # Plot seismic ======================================================
        # compute the 99th percentile and zero out all values below that. This can
        # be tweaked to plot the required amount of detail without cluttering the
        # plot

        vmm = np.percentile(svals.flatten(), percent_clip)
        svalsClipped = np.array(svals)
        svalsClipped[svalsClipped < vmm] = 0
        ci = ax.contourf(mdist, mdepth, svalsClipped, 50,
                         vmin=-vmm, vmax=vmm,
                         cmap='Greys', alpha=alpha/3., rasterized=True) # use a lower alpha value for seismic amplitudes

        ax.invert_yaxis()
        ax.set_aspect(1)
        ax.set_ylim(max_depth)

        if(fig):
            return fig, ax
        # end if
    # end func

    def save_figure(self, save_fn=None, fig_dpi=None, file_format='pdf',
                    orientation='landscape', close_fig='y'):
        """
        save_figure will save the figure to save_fn.

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
            save_fn = os.path.join(save_fn, '_E{0}_N{1}_Z{2}.{3}'.format(
                self.index_east, self.index_north,
                self.index_vertical, file_format))
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)


if __name__=='__main__':
    modem = os.path.dirname(__file__)
    modeling = os.path.dirname(modem)
    mtpy = os.path.dirname(modeling)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    ModEM_files = os.path.join(data, 'ModEM_files')

    mfn = os.path.join(ModEM_files, 'Modular_MPI_NLCG_056_im2.rho')
    dfn = os.path.join(ModEM_files, 'ModEM_Data_im2.dat')
    ps = PlotSlices(model_fn=mfn, data_fn=dfn,
                    save_path='/tmp',
                    plot_stations=True,
                    plot_yn='n')
    figs, fpaths = ps.export_slices('E-Z', [20], station_buffer=2000)

    # Updating cb-axis location. This first axis in each fig object is the
    # plot axis and the second being the colorbar axis.
    for f,fp in zip(figs, fpaths):
        cbax = f.axes[1]
        oldPos = cbax.get_position()  # get the original position
        newPos = [oldPos.x0, oldPos.y0, oldPos.width / 2.0, oldPos.height / 2.0]
        cbax.set_position(newPos)
        f.savefig(fp, dpi=ps.fig_dpi)


    # Exporting slices without saving
    figs, fpaths = ps.export_slices('E-Z', [20], station_buffer=2000, save=False)
    figs[0].savefig('/tmp/f.png', dpi=600)


    # Fetch a profile along station locations
    gd, gz, gv = ps.get_slice("STA", nsteps=1000)
