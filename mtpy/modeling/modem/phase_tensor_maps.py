#!/usr/bin/env python
"""
    Create Phase Tensor Map from the ModEM's output Resistivity model
"""
import os
import os.path as op

import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse

import mtpy.analysis.pt as mtpt
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtplottools
import mtpy.modeling.ws3dinv as ws
import mtpy.utils.exceptions as mtex
from mtpy.utils.calculator import nearest_index
from mtpy.utils.gis_tools import epsg_project
from mtpy.utils import basemap_tools
from mtpy.modeling.modem import Data, Model
import logging, traceback

try:
    from pyevtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need to pip install PyEVTK:'
           ' https://bitbucket.org/pauloh/pyevtk')

    print ('Note: if you are using Windows you should build evtk first with'
           'either MinGW or cygwin using the command: \n'
           '    python setup.py build -compiler=mingw32  or \n'
           '    python setup.py build -compiler=cygwin')


# ==============================================================================
# plot phase tensors
# ==============================================================================
class PlotPTMaps(mtplottools.MTEllipse):
    """
    Plot phase tensor maps including residual pt if response file is input.

    :Plot only data for one period: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, plot_period_list=[0])

    :Plot data and model response: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> rfn = r"/home/MT/ws3dinv/Inv1/Test_resp.00"
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, resp_fn=rfn, model_fn=mfn,
        >>> ...                 plot_period_list=[0])
        >>> # adjust colorbar
        >>> ptm.cb_res_pad = 1.25
        >>> ptm.redraw_plot()


    ========================== ================================================
    Attributes                 Description
    ========================== ================================================
    cb_pt_pad                  percentage from top of axes to place pt
                               color bar. *default* is 1.2
    cb_res_pad                 percentage from bottom of axes to place
                               resistivity color bar. *default* is 0.5
    cb_residual_tick_step      tick step for residual pt. *default* is 3
    cb_tick_step               tick step for phase tensor color bar,
                               *default* is 45
    data_obj                   data object (read in from ModEM data file)
    data_fn                    full path to data fle
    dscale                     scaling parameter depending on map_scale
    ellipse_cmap               color map for pt ellipses. *default* is
                               mt_bl2gr2rd
    ellipse_colorby            [ 'skew' | 'skew_seg' | 'phimin' | 'phimax'|
                                 'phidet' | 'ellipticity' ] parameter to color
                                 ellipses by. *default* is 'phimin'
    ellipse_range              (min, max, step) min and max of colormap, need
                               to input step if plotting skew_seg
    ellipse_size               relative size of ellipses in map_scale
    ew_limits                  limits of plot in e-w direction in map_scale
                               units.  *default* is None, scales to station
                               area
    fig_aspect                 aspect of figure. *default* is 1
    fig_dpi                    resolution in dots-per-inch. *default* is 300
    fig_list                   list of matplotlib.figure instances for each
                               figure plotted.
    fig_size                   [width, height] in inches of figure window
                               *default* is [6, 6]
    font_size                  font size of ticklabels, axes labels are
                               font_size+2. *default* is 7
    grid_east                  relative location of grid nodes in e-w direction
                               in map_scale units
    grid_north                 relative location of grid nodes in n-s direction
                               in map_scale units
    grid_z                     relative location of grid nodes in z direction
                               in map_scale units
    model_fn                 full path to initial file
    map_scale                  [ 'km' | 'm' ] distance units of map.
                               *default* is km
    mesh_east                  np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north                 np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                   full path to model file
    nodes_east                 relative distance betwen nodes in e-w direction
                               in map_scale units
    nodes_north                relative distance betwen nodes in n-s direction
                               in map_scale units
    nodes_z                    relative distance betwen nodes in z direction
                               in map_scale units
    ns_limits                  (min, max) limits of plot in n-s direction
                               *default* is None, viewing area is station area
    pad_east                   padding from extreme stations in east direction
    pad_north                  padding from extreme stations in north direction
    period_list                list of periods from data
    plot_grid                  [ 'y' | 'n' ] 'y' to plot grid lines
                               *default* is 'n'
    plot_period_list           list of period index values to plot
                               *default* is None
    plot_yn                    ['y' | 'n' ] 'y' to plot on instantiation
                               *default* is 'y'
    res_cmap                   colormap for resisitivity values.
                               *default* is 'jet_r'
    res_limits                 (min, max) resistivity limits in log scale
                               *default* is (0, 4)
    res_model                  np.ndarray(n_north, n_east, n_vertical) of
                               model resistivity values in linear scale
    residual_cmap              color map for pt residuals.
                               *default* is 'mt_wh2or'
    resp                       np.ndarray(n_stations, n_periods, 2, 2)
                               impedance tensors for model response
    resp_fn                    full path to response file
    save_path                  directory to save figures to
    save_plots                 [ 'y' | 'n' ] 'y' to save plots to save_path
    station_east               location of stations in east direction in
                               map_scale units
    station_fn                 full path to station locations file
    station_names              station names
    station_north              location of station in north direction in
                               map_scale units
    subplot_bottom             distance between axes and bottom of figure window
    subplot_left               distance between axes and left of figure window
    subplot_right              distance between axes and right of figure window
    subplot_top                distance between axes and top of figure window
    title                      titiel of plot *default* is depth of slice
    xminorticks                location of xminorticks
    yminorticks                location of yminorticks
    ========================== ================================================
    """

    def __init__(self, data_fn=None, resp_fn=None, model_fn=None, **kwargs):
        # MTEllipse.__init__(self, **kwargs)
        super(PlotPTMaps, self).__init__(**kwargs)

        self.model_fn = model_fn
        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)

        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

        self.save_plots = kwargs.pop('save_plots', 'y')
        self.plot_period_list = kwargs.pop('plot_period_list', None)
        self.period_dict = None
        self.d_index = kwargs.pop('d_index',None)

        self.map_scale = kwargs.pop('map_scale', 'km')
        # make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
        elif self.map_scale == 'm':
            self.dscale = 1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)

        self.pad_east = kwargs.pop('pad_east', 2000)
        self.pad_north = kwargs.pop('pad_north', 2000)

        self.plot_grid = kwargs.pop('plot_grid', 'n')

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.fig_list = []

        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)

        self.residual_cmap = kwargs.pop('residual_cmap', 'mt_wh2or')
        self.font_size = kwargs.pop('font_size', 7)

        self.cb_pt_pad = kwargs.pop('cb_pt_pad', 1.2)
        self.cb_res_pad = kwargs.pop('cb_res_pad', .5)

        self.res_limits = kwargs.pop('res_limits', (0, 4))
        self.res_cmap = kwargs.pop('res_cmap', 'jet_r')

        # --> set the ellipse properties -------------------
        self._ellipse_dict = kwargs.pop('ellipse_dict',
                                        {'size': 2,
                                         'ellipse_range':[0,0],
                                         'ellipse_colorby':'phimin',
                                         'ellipse_cmap':'mt_bl2gr2rd',
                                         'normalise':False})

        self._read_ellipse_dict(self._ellipse_dict)

        self.ellipse_size = kwargs.pop(
            'ellipse_size', self._ellipse_dict['size'])
            
        self.normalise_ellipses = kwargs.pop('normalise_ellipses',False)

        self.cb_tick_step = kwargs.pop('cb_tick_step', None)
        # update default colorbar tick step based on ellipse_range
        if self.cb_tick_step is None:
            self.cb_tick_step = int((self.ellipse_range[1] - self.ellipse_range[0])/2.)
            
        self.cb_residual_tick_step = kwargs.pop('cb_residual_tick_step', 3)        

        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .2
        self.subplot_wspace = .05

        self.data_obj = None
        self.resp_obj = None
        self.model_obj = None
        self.period_list = None

        self.pt_data_arr = None
        self.pt_resp_arr = None
        self.pt_resid_arr = None
        
        self.residual_pt_type= kwargs.pop('residual_pt_type','heise')

        # FZ: do not call plot in the constructor! it's not pythonic
        self.plot_yn = kwargs.pop('plot_yn', 'n')
        if self.plot_yn == 'y':
            self.plot()
        else:
            self._read_files()

    def _read_files(self):
        """
        get information from files
        """

        # --> read in data file
        self.data_obj = Data()
        self.data_obj.read_data_file(self.data_fn)

        # --> read response file
        if self.resp_fn is not None:
            self.resp_obj = Data()
            self.resp_obj.read_data_file(self.resp_fn)

        # --> read mode file
        if self.model_fn is not None:
            self.model_obj = Model()
            self.model_obj.read_model_file(self.model_fn)

        self._get_plot_period_list()
        self._get_pt()

    def _get_plot_period_list(self):
        """
        get periods to plot from input or data file
        """
        # --> get period list to plot
        if self.plot_period_list is None:
            self.plot_period_list = self.data_obj.period_list
        else:
            if isinstance(self.plot_period_list, list):
                # check if entries are index values or actual periods
                if isinstance(self.plot_period_list[0], int):
                    self.plot_period_list = [self.period_list[ii]
                                             for ii in self.plot_period_list]
                else:
                    pass
            elif isinstance(self.plot_period_list, int):
                self.plot_period_list = self.period_list[self.plot_period_list]
            elif isinstance(self.plot_period_list, float):
                self.plot_period_list = [self.plot_period_list]

        self.period_dict = dict([(key, value) for value, key in
                                 enumerate(self.data_obj.period_list)])

    def _get_pt(self):
        """
        put pt parameters into something useful for plotting
        """

        ns = len(list(self.data_obj.mt_dict.keys()))
        nf = len(self.data_obj.period_list)

        data_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                ('phimax', np.float),
                                                ('skew', np.float),
                                                ('azimuth', np.float),
                                                ('east', np.float),
                                                ('north', np.float),
                                                ('lon', np.float),
                                                ('lat', np.float),
                                                ('station', 'S10')])
        if self.resp_fn is not None:
            model_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                     ('phimax', np.float),
                                                     ('skew', np.float),
                                                     ('azimuth', np.float),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('lon', np.float),
                                                     ('lat', np.float),
                                                     ('station', 'S10')])

            res_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                   ('phimax', np.float),
                                                   ('skew', np.float),
                                                   ('azimuth', np.float),
                                                   ('east', np.float),
                                                   ('north', np.float),
                                                   ('lon', np.float),
                                                   ('lat', np.float),
                                                   ('geometric_mean', np.float),
                                                   ('station', 'S10')])

        for ii, key in enumerate(self.data_obj.mt_dict.keys()):
            east = self.data_obj.mt_dict[key].grid_east / self.dscale
            north = self.data_obj.mt_dict[key].grid_north / self.dscale
            lon = self.data_obj.mt_dict[key].lon
            lat = self.data_obj.mt_dict[key].lat
            dpt = self.data_obj.mt_dict[key].pt
            data_pt_arr[:, ii]['east'] = east
            data_pt_arr[:, ii]['north'] = north
            data_pt_arr[:, ii]['lon'] = lon
            data_pt_arr[:, ii]['lat'] = lat
            data_pt_arr[:, ii]['phimin'] = dpt.phimin
            data_pt_arr[:, ii]['phimax'] = dpt.phimax
            data_pt_arr[:, ii]['azimuth'] = dpt.azimuth
            data_pt_arr[:, ii]['skew'] = dpt.beta
            data_pt_arr[:, ii]['station'] = self.data_obj.mt_dict[key].station
            if self.resp_fn is not None:
                mpt = self.resp_obj.mt_dict[key].pt
                try:
                    rpt = mtpt.ResidualPhaseTensor(pt_object1=dpt,
                                                   pt_object2=mpt,
                                                   residualtype=self.residual_pt_type)
                    rpt = rpt.residual_pt
                    res_pt_arr[:, ii]['east'] = east
                    res_pt_arr[:, ii]['north'] = north
                    res_pt_arr[:, ii]['lon'] = lon
                    res_pt_arr[:, ii]['lat'] = lat
                    res_pt_arr[:, ii]['phimin'] = rpt.phimin
                    res_pt_arr[:, ii]['phimax'] = rpt.phimax
                    res_pt_arr[:, ii]['azimuth'] = rpt.azimuth
                    res_pt_arr[:, ii]['skew'] = rpt.beta
                    res_pt_arr[:, ii]['station'] = self.data_obj.mt_dict[key].station
                    res_pt_arr[:, ii]['geometric_mean'] = np.sqrt(np.abs(rpt.phimin) *
                                                                  np.abs(rpt.phimax))
                except mtex.MTpyError_PT:
                    print(key, dpt.pt.shape, mpt.pt.shape)

                model_pt_arr[:, ii]['east'] = east
                model_pt_arr[:, ii]['north'] = north
                model_pt_arr[:, ii]['lon'] = lon
                model_pt_arr[:, ii]['lat'] = lat
                model_pt_arr[:, ii]['phimin'] = mpt.phimin
                model_pt_arr[:, ii]['phimax'] = mpt.phimax
                model_pt_arr[:, ii]['azimuth'] = mpt.azimuth
                model_pt_arr[:, ii]['skew'] = mpt.beta
                model_pt_arr[
                    :, ii]['station'] = self.data_obj.mt_dict[key].station

        # make these attributes
        self.pt_data_arr = data_pt_arr
        if self.resp_fn is not None:
            self.pt_resp_arr = model_pt_arr
            self.pt_resid_arr = res_pt_arr

    def plot_on_axes(self, ax, m, periodIdx, ptarray='data', ellipse_size_factor=10000,
                     cvals=None, map_scale='m', centre_shift=[0, 0], plot_tipper='n',
                     tipper_size_factor=1e5, **kwargs):

        '''
        Plots phase tensors for a given period index.

        :param ax: plot axis
        :param m: basemap instance
        :param periodIdx: period index
        :param ptarray: name of data-array to access for retrieving attributes;
                        can be either 'data', 'resp' or 'resid'
        :param ellipse_size_factor: factor to control ellipse size
        :param cvals: list of colour values for colouring each ellipse; must be of
                      the same length as the number of tuples for each period
        :param map_scale: map length scale
        :param kwargs: list of relevant matplotlib arguments (e.g. zorder, alpha, etc.)
        :param plot_tipper: string ('n', 'yr', 'yi', or 'yri') to plot
                            no tipper, real only, imaginary only, or both
        :param tipper_size_factor: scaling factor for tipper vectors
        '''

        assert (periodIdx >= 0 and periodIdx < len(self.data_obj.period_list)), \
            'Error: Index for plot-period out of bounds.'

        k = periodIdx
        pt_array = getattr(self, 'pt_' + ptarray + '_arr')

        for i in range(len(pt_array[k])):
            lon = pt_array[k]['lon'][i]
            lat = pt_array[k]['lat'][i]
            if self.normalise_ellipses:
                phimax = pt_array[k]['phimax'][i] / pt_array[k]['phimax'][i]
                phimin = pt_array[k]['phimin'][i] / pt_array[k]['phimax'][i]
            else:
                phimax = pt_array[k]['phimax'][i] / pt_array[k]['phimax'].max()
                phimin = pt_array[k]['phimin'][i] / pt_array[k]['phimax'].max()
            az = pt_array[k]['azimuth'][i]
            if ptarray == 'resid':
                phimin = np.abs(phimin)
            nskew = pt_array[k]['skew'][i]

            # print az
            if (phimax > 0 and phimin > 0):
                c = None
                if (cvals is not None): c = cvals[i]
                if (c is not None): kwargs['facecolor'] = c

                if m is None:
                    x = pt_array[k]['east'][i]
                    y = pt_array[k]['north'][i]
                    if map_scale == 'km':
                        x /= 1e3
                        y /= 1e3
                else:
                    x, y = m(lon, lat)

                # matplotlib angles are defined as degrees anticlockwise from positive x direction.
                # therefore we need to adjust az accordingly
                e = Ellipse([x, y],
                            phimax * ellipse_size_factor,
                            phimin * ellipse_size_factor,
                            90. - az, **kwargs)
                ax.add_artist(e)
            # end if
        # end for
        if 'y' in plot_tipper:
            # if neither r or i provided, assume that we want to plot both
            if plot_tipper == 'y':
                plot_tipper = 'yri'
            self._plot_induction_vectors(ax, m, periodIdx, 
                                         ptarray=ptarray, size_factor=tipper_size_factor,
                                         map_scale=map_scale, centre_shift=centre_shift,
                                         plot_tipper=plot_tipper, **kwargs)
    # end func

        


    def plot(self, period = None, periodIdx = 0, save2file=None, **kwargs):
        """ Plot phase tensor maps for data and or response, each figure is of a
        different period.  If response is input a third column is added which is
        the residual phase tensor showing where the model is not fitting the data
        well.  The data is plotted in km.

        Args:
            period: the period index to plot, default=0

        Returns:

        """

        print(("The input parameter period is", period))

        # --> read in data first
        if self.data_obj is None:
            self._read_files()

        # set plot properties
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        font_dict = {'size': self.font_size + 2, 'weight': 'bold'}

        # make a grid of subplots
        gs = gridspec.GridSpec(1, 3, hspace=self.subplot_hspace,
                               wspace=self.subplot_wspace)

        # set some parameters for the colorbar
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                ckstep = 3
        bounds = np.arange(ckmin, ckmax + ckstep, ckstep)

        # set plot limits to be the station area
        if self.ew_limits is None:
            east_min = self.data_obj.data_array['rel_east'].min() - \
                self.pad_east
            east_max = self.data_obj.data_array['rel_east'].max() + \
                self.pad_east
            self.ew_limits = (east_min / self.dscale, east_max / self.dscale)

        if self.ns_limits is None:
            north_min = self.data_obj.data_array['rel_north'].min() - \
                self.pad_north
            north_max = self.data_obj.data_array['rel_north'].max() + \
                self.pad_north
            self.ns_limits = (north_min / self.dscale, north_max / self.dscale)

        # -------------plot phase tensors------------------------------------
        if period > len(self.plot_period_list) - 1:
            print((
                "Error: the period exceeds the max value:", len(
                    self.plot_period_list) - 1))

        # FZ: changed below to plot a given period index
        # for ff, per in enumerate(self.plot_period_list):
        # first, reset fig list
        self.fig_list = []
        for ff, per in enumerate(self.plot_period_list[period:period + 1]):
            
            data_ii = self.period_dict[per]

            print('Plotting Period: {0:.5g}'.format(per))
            fig = plt.figure('{0:.5g}'.format(per), figsize=self.fig_size,
                             dpi=self.fig_dpi)
            fig.clf()

            if self.resp_fn is not None:
                axd = fig.add_subplot(gs[0, 0], aspect='equal')
                axm = fig.add_subplot(gs[0, 1], aspect='equal')
                axr = fig.add_subplot(gs[0, 2], aspect='equal')
                ax_list = [axd, axm, axr]

            else:
                axd = fig.add_subplot(gs[0, :], aspect='equal')
                ax_list = [axd]

            # plot model below the phase tensors
            if self.model_fn is not None:
                gridzcentre = np.mean(
                    [self.model_obj.grid_z[1:], self.model_obj.grid_z[:-1]], axis=0)
                if self.d_index is not None:
                    approx_depth, d_index = ws.estimate_skin_depth(self.model_obj.res_model.copy(),
                                                                   gridzcentre / self.dscale,
                                                                   per,
                                                                   dscale=self.dscale)
                else:
                    d_index = self.d_index
                    approx_depth = self.model_obj.grid_z[d_index]

                # need to add an extra row and column to east and north to make sure
                # all is plotted see pcolor for details.
                plot_east = np.append(self.model_obj.grid_east,
                                      self.model_obj.grid_east[-1] * 1.25) / \
                    self.dscale
                plot_north = np.append(self.model_obj.grid_north,
                                       self.model_obj.grid_north[-1] * 1.25) / \
                    self.dscale

                # make a mesh grid for plotting
                # the 'ij' makes sure the resulting grid is in east, north
                try:
                    self.mesh_east, self.mesh_north = np.meshgrid(plot_east,
                                                                  plot_north,
                                                                  indexing='ij')
                except TypeError:
                    self.mesh_east, self.mesh_north = [arr.T for arr in np.meshgrid(plot_east,
                                                                                    plot_north)]

                for ax in ax_list:
                    plot_res = np.log10(
                        self.model_obj.res_model[
                            :, :, d_index].T)
                    ax.pcolormesh(self.mesh_east,
                                  self.mesh_north,
                                  plot_res,
                                  cmap=self.res_cmap,
                                  vmin=self.res_limits[0],
                                  vmax=self.res_limits[1])

            # --> plot data phase tensors
            print(kwargs)
            for pt in self.pt_data_arr[data_ii]:
                eheight = pt['phimin'] / \
                    self.pt_data_arr[data_ii]['phimax'].max() * \
                    self.ellipse_size
                ewidth = pt['phimax'] / \
                    self.pt_data_arr[data_ii]['phimax'].max() * \
                    self.ellipse_size

                ellipse = Ellipse((pt['east'],
                                   pt['north']),
                                  width=ewidth,
                                  height=eheight,
                                  angle=90 - pt['azimuth'],
                                  **kwargs)

                # get ellipse color
                if self.ellipse_cmap.find('seg') > 0:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                              self.ellipse_colorby,
                                                              self.ellipse_cmap,
                                                              ckmin,
                                                              ckmax,
                                                              bounds=bounds))
                else:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                              self.ellipse_colorby,
                                                              self.ellipse_cmap,
                                                              ckmin,
                                                              ckmax))

                axd.add_artist(ellipse)

            # -----------plot response phase tensors---------------
            if self.resp_fn is not None:
                rcmin = np.floor(self.pt_resid_arr['geometric_mean'].min())
                rcmax = np.floor(self.pt_resid_arr['geometric_mean'].max())
                for mpt, rpt in zip(self.pt_resp_arr[data_ii],
                                    self.pt_resid_arr[data_ii]):
                    eheight = mpt['phimin'] / \
                        self.pt_resp_arr[data_ii]['phimax'].max() * \
                        self.ellipse_size
                    ewidth = mpt['phimax'] / \
                        self.pt_resp_arr[data_ii]['phimax'].max() * \
                        self.ellipse_size

                    ellipsem = Ellipse((mpt['east'],
                                        mpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=90 - mpt['azimuth'],
                                       **kwargs)

                    # get ellipse color
                    if self.ellipse_cmap.find('seg') > 0:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                                   self.ellipse_colorby,
                                                                   self.ellipse_cmap,
                                                                   ckmin,
                                                                   ckmax,
                                                                   bounds=bounds))
                    else:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                                   self.ellipse_colorby,
                                                                   self.ellipse_cmap,
                                                                   ckmin,
                                                                   ckmax))

                    axm.add_artist(ellipsem)

                    # -----------plot residual phase tensors---------------
                    eheight = rpt['phimin'] / \
                        self.pt_resid_arr[data_ii]['phimax'].max() * \
                        self.ellipse_size
                    ewidth = rpt['phimax'] / \
                        self.pt_resid_arr[data_ii]['phimax'].max() * \
                        self.ellipse_size

                    ellipser = Ellipse((rpt['east'],
                                        rpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=rpt['azimuth'],
                                       **kwargs)

                    # get ellipse color
                    rpt_color = np.sqrt(abs(rpt['phimin'] * rpt['phimax']))
                    if self.ellipse_cmap.find('seg') > 0:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                                   'geometric_mean',
                                                                   self.residual_cmap,
                                                                   ckmin,
                                                                   ckmax,
                                                                   bounds=bounds))
                    else:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                                   'geometric_mean',
                                                                   self.residual_cmap,
                                                                   ckmin,
                                                                   ckmax))

                    axr.add_artist(ellipser)

            # --> set axes properties
            # data
            axd.set_xlim(self.ew_limits)
            axd.set_ylim(self.ns_limits)
            axd.set_xlabel('Easting ({0})'.format(self.map_scale),
                           fontdict=font_dict)
            axd.set_ylabel('Northing ({0})'.format(self.map_scale),
                           fontdict=font_dict)
            # make a colorbar for phase tensors
            # bb = axd.axes.get_position().bounds
            bb = axd.get_position().bounds
            y1 = .25 * (2 + (self.ns_limits[1] - self.ns_limits[0]) /
                        (self.ew_limits[1] - self.ew_limits[0]))
            cb_location = (3.35 * bb[2] / 5 + bb[0],
                           y1 * self.cb_pt_pad, .295 * bb[2], .02)
            cbaxd = fig.add_axes(cb_location)
            if self.ellipse_cmap in list(mtcl.cmapdict.keys()):
                ecmap = mtcl.cmapdict[self.ellipse_cmap]
            else:
                ecmap = self.ellipse_cmap
            cbd = mcb.ColorbarBase(cbaxd,
                                   cmap=ecmap,
                                   norm=Normalize(vmin=ckmin,
                                                  vmax=ckmax),
                                   orientation='horizontal')
            cbd.ax.xaxis.set_label_position('top')
            cbd.ax.xaxis.set_label_coords(.5, 1.75)
            cbd.set_label(mtplottools.ckdict[self.ellipse_colorby])
            cbd.set_ticks(np.arange(ckmin, ckmax + self.cb_tick_step,
                                    self.cb_tick_step))

            axd.text(self.ew_limits[0] * .95,
                     self.ns_limits[1] * .95,
                     'Data',
                     horizontalalignment='left',
                     verticalalignment='top',
                     bbox={'facecolor': 'white'},
                     fontdict={'size': self.font_size + 1})

            # Model and residual
            if self.resp_fn is not None:
                for aa, ax in enumerate([axm, axr]):
                    ax.set_xlim(self.ew_limits)
                    ax.set_ylim(self.ns_limits)
                    ax.set_xlabel('Easting ({0})'.format(self.map_scale),
                                  fontdict=font_dict)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    # make a colorbar ontop of axis
                    bb = ax.axes.get_position().bounds
                    y1 = .25 * (2 + (self.ns_limits[1] - self.ns_limits[0]) /
                                (self.ew_limits[1] - self.ew_limits[0]))
                    cb_location = (3.35 * bb[2] / 5 + bb[0],
                                   y1 * self.cb_pt_pad, .295 * bb[2], .02)
                    cbax = fig.add_axes(cb_location)
                    if aa == 0:
                        cb = mcb.ColorbarBase(cbax,
                                              cmap=mtcl.cmapdict[
                                                  self.ellipse_cmap],
                                              norm=Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                              orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(mtplottools.ckdict[self.ellipse_colorby])
                        cb.set_ticks(np.arange(ckmin, ckmax + self.cb_tick_step,
                                               self.cb_tick_step))
                        ax.text(self.ew_limits[0] * .95,
                                self.ns_limits[1] * .95,
                                'Model',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor': 'white'},
                                fontdict={'size': self.font_size + 1})
                    else:
                        cb = mcb.ColorbarBase(cbax,
                                              cmap=mtcl.cmapdict[
                                                  self.residual_cmap],
                                              norm=Normalize(vmin=rcmin,
                                                             vmax=rcmax),
                                              orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(r"$\sqrt{\Phi_{min} \Phi_{max}}$")
                        cb_ticks = [rcmin, (rcmax - rcmin) / 2, rcmax]
                        cb.set_ticks(cb_ticks)
                        ax.text(self.ew_limits[0] * .95,
                                self.ns_limits[1] * .95,
                                'Residual',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor': 'white'},
                                fontdict={'size': self.font_size + 1})

            if self.model_fn is not None:
                for ax in ax_list:
                    ax.tick_params(direction='out')
                    bb = ax.axes.get_position().bounds
                    y1 = .25 * (2 - (self.ns_limits[1] - self.ns_limits[0]) /
                                (self.ew_limits[1] - self.ew_limits[0]))
                    cb_position = (3.0 * bb[2] / 5 + bb[0],
                                   y1 * self.cb_res_pad, .35 * bb[2], .02)
                    cbax = fig.add_axes(cb_position)
                    cb = mcb.ColorbarBase(cbax,
                                          cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits[0],
                                                         vmax=self.res_limits[1]),
                                          orientation='horizontal')
                    cb.ax.xaxis.set_label_position('top')
                    cb.ax.xaxis.set_label_coords(.5, 1.5)
                    cb.set_label('Resistivity ($\Omega \cdot$m)')
                    cb_ticks = np.arange(np.floor(self.res_limits[0]),
                                         np.ceil(self.res_limits[1] + 1), 1)
                    cb.set_ticks(cb_ticks)
                    cb.set_ticklabels([mtplottools.labeldict[ctk]
                                       for ctk in cb_ticks])

            if save2file is not None:
                fig.savefig(save2file, dpi=self.fig_dpi, bbox_inches='tight')

            plt.show()
            self.fig_list.append(fig)

            return fig

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
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()

    
    def _plot_induction_vectors(self, ax, m, periodIdx, ptarray='data', size_factor=10000,
                                 map_scale='m', centre_shift=[0, 0], plot_tipper='yri', **kwargs):
                                                  
        if ptarray == 'data':
            data_array = self.data_obj.data_array
        elif ptarray == 'resp':
            data_array = self.resp_obj.data_array
        rx = data_array['tip'].real[:,periodIdx,0,0]
        ry = data_array['tip'].real[:,periodIdx,0,1]
        ix = data_array['tip'].imag[:,periodIdx,0,0]
        iy = data_array['tip'].imag[:,periodIdx,0,1]
        
        lon,lat = self.data_obj.station_locations.lon, self.data_obj.station_locations.lat
        x,y = m(lon,lat)

        kwargs_tip = {'length_includes_head':True, 
                      'head_width':size_factor*0.07, 
                      'head_length': size_factor*0.1}
        kwargs_tip.update(kwargs)
        
        for sidx in range(len(self.data_obj.data_array)):
            if 'r' in plot_tipper:
                ax.arrow(x[sidx],y[sidx],size_factor*rx[sidx],size_factor*ry[sidx],color='k',**kwargs_tip)
            if 'i' in plot_tipper:
                ax.arrow(x[sidx],y[sidx],size_factor*ix[sidx],size_factor*iy[sidx],color='b',**kwargs_tip)


    def _get_pt_data_list(self, attribute, xykeys=['east', 'north']):

        headerlist = ['period', 'station'] + xykeys + \
            ['azimuth', 'phimin', 'phimax', 'skew']
        data = getattr(self, attribute).T.copy()
        indices = np.argsort(data['station'][:, 0])

        data = data[indices].T
        dtype = []
        for val in headerlist:
            if val == 'station':
                dtype.append((val, 'S10'))
            else:
                dtype.append((val, np.float))

        data_to_write = np.zeros(np.product(data.shape), dtype=dtype)
        data_to_write['period'] = np.vstack(
            [self.plot_period_list] * data.shape[1]).T.flatten()

        for val in headerlist[1:]:
            if val in ['east', 'north']:
                data[val] *= self.dscale
            data_to_write[val] = data[val].flatten()

        return data_to_write, headerlist

    def get_period_attributes(self, periodIdx, key, ptarray='data'):
        '''
        Returns, for a given period, a list of attribute values for key
        (e.g. skew, phimax, etc.).

        :param periodIdx: index of period; print out _plot_period for periods available
        :param key: attribute key
        :param ptarray: name of data-array to access for retrieving attributes;
                        can be either 'data', 'resp' or 'resid'
        :return: numpy array of attribute values
        '''

        # load data if necessary
        if self.data_obj is None:
            self._read_files()

        assert (periodIdx >= 0 and periodIdx < len(self.plot_period_list)), \
            'Error: Index for plot-period out of bounds.'

        pk = periodIdx
        try:
            print("getting", key)
            if key == 'phimean':
                vals = np.mean([getattr(self, 'pt_' + ptarray + '_arr')[pk]['phimin'],
                                getattr(self, 'pt_' + ptarray + '_arr')[pk]['phimax']],axis=0)
            else:
                vals = getattr(self, 'pt_' + ptarray + '_arr')[pk][key]
            return vals
        except:

            print('Attribute %s not found' % ('pt_' + ptarray + '_arr'))
            logging.error(traceback.format_exc())
            exit(-1)

        return None
    # end func

    def write_pt_data_to_text(self, savepath='.'):

        if self.pt_data_arr is None:
            self._read_files()

        for att in ['pt_data_arr', 'pt_resp_arr', 'pt_resid_arr']:
            if hasattr(self, att):
                data_to_write, headerlist = self._get_pt_data_list(att)
                header = ' '.join(headerlist)

                filename = op.join(savepath, att[:-4] + '.txt')
                if att == 'pt_resid_arr':
                    data_to_write['azimuth'] = 90. - data_to_write['azimuth']
                np.savetxt(filename, data_to_write, header=header,
                           fmt=['%.4e', '%s', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.3f'])


    def write_pt_data_to_gmt(self, period=None, epsg=None, savepath='.', center_utm=None,
                             colorby='phimin', attribute='data', clim=None):
        """
        write data to plot phase tensor ellipses in gmt.
        saves a gmt script and text file containing ellipse data

        provide:
        period to plot (seconds)
        epsg for the projection the model was projected to
        (google "epsg your_projection_name" and you will find it)
        centre_utm - utm coordinates for centre position of model, if not
                     provided, script will try and extract it from data file
        colorby - what to colour the ellipses by, 'phimin', 'phimax', or 'skew'
        attribute - attribute to plot 'data', 'resp', or 'resid' for data,
                    response or residuals

        """

        att = 'pt_{}_arr'.format(attribute)

        # if centre utm not provided, get station locations from the data
        # object
        project = False
        xykeys = ['lon', 'lat']

        if epsg is not None:
            if center_utm is not None:
                project = True
            else:
                if hasattr(self.data_obj, 'center_position'):
                    if np.all(np.array(self.data_obj.center_position) > 0):
                        project = True
                        center_utm = self.data_obj.project_xy(self.data_obj.center_position[0],
                                                              self.data_obj.center_position[
                                                                  1],
                                                              epsg_from=4326, epsg_to=epsg)
        if project:
            xykeys = ['east', 'north']

        # get text data list
        data, headerlist = self._get_pt_data_list(att, xykeys=xykeys)

        # extract relevant columns in correct order
        periodlist = data['period']

        columns = xykeys + [colorby, 'azimuth', 'phimax', 'phimin']
        gmtdata = np.vstack([data[i] for i in columns]).T

        # make a filename based on period
        if period >= 1.:
            suffix = '%1i' % round(period)
        else:
            nzeros = np.abs(np.int(np.floor(np.log10(period))))
            fmt = '%0' + str(nzeros + 1) + 'i'
            suffix = fmt % (period * 10 ** nzeros)

        filename = 'ellipse_' + attribute + '.' + suffix

        if period is not None:
            # extract relevant period
            unique_periods = np.unique(periodlist)
            closest_period = unique_periods[np.abs(unique_periods - period) ==
                                            np.amin(np.abs(unique_periods - period))]
            # indices to select all occurrances of relevant period (to nearest
            # 10^-8 s)
            pind = np.where(np.abs(closest_period - periodlist) < 1e-8)[0]
        else:
            # take the first period
            pind = 0

        # select relevant periods
        periodlist, gmtdata = periodlist[pind], gmtdata[pind]

        if project:
            gmtdata[:, 0] += center_utm[0]
            gmtdata[:, 1] += center_utm[1]

            # now that x y coordinates are in utm, project to lon/lat
            self.data_obj.epsg = epsg
            gmtdata[
                :, 0], gmtdata[
                :, 1] = self.data_obj.project_xy(
                gmtdata[
                    :, 0], gmtdata[
                    :, 1])

        if self.normalise_ellipses:
            norm = gmtdata[:,4]
        else:
            # normalise by maximum value of phimax
            norm = np.amax(gmtdata[:, 4])
        gmtdata[:, 5] /= norm
        gmtdata[:, 4] /= norm
        if attribute != 'resid':
            gmtdata[:, 3] = 90. - gmtdata[:, 3]

        # write to text file in correct format
        fmt = ['%+11.6f', '%+10.6f'] + ['%+9.4f'] * 2 + ['%8.4f'] * 2
        np.savetxt(op.join(savepath, filename), gmtdata, fmt)

        # write gmt script
        xmin, xmax = gmtdata[:, 0].min(), gmtdata[:, 0].max()
        ymin, ymax = gmtdata[:, 1].min(), gmtdata[:, 1].max()

        pad = min(ymax - ymin, xmax - xmin) / 10.
        tr = -int(np.log10(20. * (xmax - xmin)))
        tickspacing = int(np.round(20. * (xmax - xmin), tr))
        scalebarlat = int(round(ymax + ymin) / 2.)
        if clim is None:
            cr = int(np.ceil(-np.log10(np.amax(gmtdata[:, 2]))))
            clim = np.round([gmtdata[:, 2].min(), gmtdata[
                            :, 2].max()], cr).astype(int)

        gmtlines = [line + '\n' for line in ['w={}'.format(xmin - pad),
                                             'e={}'.format(xmax + pad),
                                             's={}'.format(ymin - pad),
                                             'n={}'.format(ymax + pad),
                                             r"wesn=$w/$s/$e/$n'r'",
                                             '',
                                             '# define output file and remove it if it exists',
                                             'PS={}.ps'.format(
                                                 filename.replace('.', '')),
                                             'rm $PS',
                                             '',
                                             '# set gmt parameters',
                                             'gmtset FORMAT_GEO_MAP ddd:mm:ss',
                                             'gmtset FONT_ANNOT_PRIMARY 9p,Helvetica,black',
                                             'gmtset MAP_FRAME_TYPE fancy',
                                             '',
                                             '# make colour palette',
                                             'makecpt -Cpolar -T{}/{} -Z > {}.cpt'.format(
                                                 clim[0], clim[1], colorby),
                                             '',
                                             '# draw coastline',
                                             'pscoast -R$wesn -JM18c -W0.5p -Ba1f1/a1f1WSen -Gwhite -Slightgrey -Lfx14c/1c/{}/{}+u -Df -P -K >> $PS'.format(
                                                 scalebarlat, tickspacing),
                                             '',
                                             '# draw ellipses',
                                             'psxy {} -R -J -P -Se -C{}.cpt -W0.01p -O >> $PS'.format(filename,
                                                                                                      colorby),
                                             '',
                                             '# save to png',
                                             'ps2raster -Tg -A -E400 $PS']]

        with open(op.join(savepath, 'gmtscript_{}.gmt'.format(attribute)), 'wb') as scriptfile:
            scriptfile.writelines(gmtlines)




    def save_all_figures(self, save_path=None, fig_dpi=None, file_format='pdf',
                    orientation='landscape', close_fig='y'):
        """
        save_figure will save all figures in fig_list to save_fn.

        Arguments:
        -----------

            **save_path** : string
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

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_path) == False:
            try:
                os.mkdir(save_path)
            except:
                raise IOError('Need to input a correct directory path')

        for fig in self.fig_list:
            per = fig.canvas.get_window_title()
            save_fn = os.path.join(save_path, 'PT_DepthSlice_{0}s.{1}'.format(
                per, file_format))
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

            if close_fig == 'y':
                plt.close(fig)

            else:
                pass

            self.fig_fn = save_fn
            print('Saved figure to: ' + self.fig_fn)
            