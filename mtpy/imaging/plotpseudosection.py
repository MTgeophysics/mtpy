# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:39:58 2013

@author: jpeacock-pr
"""

#==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import mtpy.imaging.mtplottools as mtpl
import mtpy.imaging.mtcolors as mtcl
import matplotlib.gridspec as gridspec

#==============================================================================

class PlotResPhasePseudoSection(object):
    """
    plot a resistivity and phase pseudo section for different components

    Need to input one of the following lists:

    Arguments:
    ----------
        **fn_list** : list of strings
                     full paths to .edi files to plot. *default* is None

        **z_object** : list of class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.frequency is filled.  *default* is None

        **res_object_list** : list of class mtpy.mtplottools.ResPhase
                             list of ResPhase objects. *default* is None

        **mt_object** : class mtpy.imaging.mtplot.MTplot
                        object of mtpy.imaging.mtplot.MTplot
                        *default* is None
    Optional Key Words:
    -------------------

        *aspect*: [ 'equal' | 'auto' | float ]
                  aspect ratio of each subplot (height/width),
                  *default* is 'auto'

        *cb_orientation*: [ 'vertical' | 'horizontal' ]
                          orientation of colorbars, *default* is 'vertical'

        *cb_pad*: float
                  padding between edge of plot and edge of colorbar,
                  *default* is .0375

        *cb_position*: (res_position (x, y, ds, dy),
                        phase_position(x, y, dx, dy))
                        position of colorbars on the plot, need to input both
                        resistivity and phase positions.  *default* is None
                        and will automatically place the colorbars

        *cb_shrink*: float
                     factor to shrink the colorbar relative to the height of
                     the y-axis. *default* is 0.75

        *fig_dpi*: float
                   dots-per-inch resolution of figure, *default* is 300

        *fig_num*: int
                  number of figure instance, *default* is 1

        *fig_size*: (x, y) in inches
                    figure size in inches, *default* is (8, 4)

        *font_size*: float
                     size of font for axes tick labels, note label are +1.
                     *default* is 7

        *ftol*: float
                tolerance to extract periods relative to plot_period.
                *default* is 0.1

        *imshow_interp*: [ 'none' | 'nearest' | 'bilinear' | 'bicubic' |
                           'spline16' | 'spline36' | 'hanning' | 'hamming' |
                           'hermite' | 'kaiser' | 'quadric' | 'catrom' |
                           'gaussian' | 'bessel' | 'mitchell' | 'sinc' |
                           'lanczos']
                         defines the interpolation method if plot_style is
                         'imshow', if you use 'nearest' same as pcolormesh
                         except the lateral boxes are equal size instead of
                         set in a grid like pcolormesh.  Imshow just gives
                         a smoother interpretation of the pseudosection.
                         *default* is 'bicubic'

        *linedir*: [ 'ew' | 'ns' ]
                   predominant direction of profile line. *default* is 'ew'

        *period_limits*: (min_period, max_period)
                         limits on the plotting period in a linear scale.
                         *default* is None and extracts limits from data

        *phase_cmap*: [ 'mt_yl2rd' | 'mt_bl2yl2rd' | 'mt_wh2bl' |
                        'mt_rd2bl' | 'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                        'mt_rd2gr2bl' ]
                      color map for phase plots, at the moment only supports
                      colormaps from mtpy.imaging.mtcolors, which are:

                           - 'mt_yl2rd' -> yellow to red
                           - 'mt_bl2yl2rd' -> blue to yellow to red
                           - 'mt_wh2bl' -> white to blue
                           - 'mt_rd2bl' -> red to blue
                           - 'mt_bl2wh2rd' -> blue to white to red
                           - 'mt_bl2gr2rd' -> blue to green to red *default*
                           - 'mt_rd2gr2bl' -> red to green to blue
                           - 'mt_seg_bl2wh2rd' -> discrete blue to
                                                 white to red

        *phase_limits*: (min_phase, max_phase)
                        minimum and maximum phase in degrees for coloring
                        *default* is (0, 90)

        *plot_period*: np.ndarray(periods)
                       array of periods to plot.  *default* is None, which
                       extracts the longest period from the input data
                       assuming that the longest is the most complete, if it
                       is not input manually.  If there are same lengths, it
                       picks the first one it finds.

        *plot_style*: [ 'imshow' | 'pcolormesh' ]
                      type of gridding for the plot. 'imshow' plots the data
                      as an image and can be interpolated, though the image
                      is stretched to the station spacing and plot_period, the
                      cells remain of equal size, so the interpolation might be
                      a little skewed.  For an accurate location of resistivity
                      values use pcolormesh, which can plot the data on an
                      irregular grid, but with no interpolation.
                      *default* is 'imshow'

        *plot_xx*: [ 'y' | 'n' ]
                  boolean to plot Z_xx, *default* is 'n'

        *plot_xy*: [ 'y' | 'n' ]
                  boolean to plot Z_xy, *default* is 'y'

        *plot_yx*: [ 'y' | 'n' ]
                  boolean to plot Z_yx, *default* is 'y'

        *plot_yy*: [ 'y' | 'n' ]
                  boolean to plot Z_yy, *default* is 'n'

        *plot_yn*: [ 'y' | 'n' ]
                  boolean to plot on instance creation, *default* is 'y'

        *res_cmap*: [ 'mt_yl2rd' | 'mt_bl2yl2rd' | 'mt_wh2bl' |
                        'mt_rd2bl' | 'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                        'mt_rd2gr2bl' ]
                  color map for phase plots, at the moment only supports
                  colormaps from mtpy.imaging.mtcolors, which are:

                       - 'mt_yl2rd' -> yellow to red
                       - 'mt_bl2yl2rd' -> blue to yellow to red
                       - 'mt_wh2bl' -> white to blue
                       - 'mt_rd2bl' -> red to blue
                       - 'mt_bl2wh2rd' -> blue to white to red
                       - 'mt_bl2gr2rd' -> blue to green to red
                       - 'mt_rd2gr2bl' -> red to green to blue *default*
                       - 'mt_seg_bl2wh2rd' -> discrete blue to
                                             white to red
        *res_limits*: (min_resistivity, max_resistivity)
                      limits on resistivity in log scale, *default* is (0,3)

        *stationid*: (min, max)
                     min and max indicies to extract from each station name.
                     *default* is (0,4)

        *text_location*: (x, y)
                         location for text label for each resistivity subplot.
                         location is in relative coordinates of the data.
                         *default* is None, which locates the label in the
                         upper left hand corner

        *text_size*: [ size in points | 'xx-small' | 'x-small' | 'small' |
                      'medium' | 'large' | 'x-large' | 'xx-large' ]
                     size of text for subplot label, *default* is font_size

        *text_weight*: [ a numeric value in range 0-1000 | 'ultralight' |
                        'light' | 'normal' | 'regular' | 'book' | 'medium' |
                        'roman' | 'semibold' | 'demibold' | 'demi' | 'bold' |
                        'heavy' | 'extra bold' | 'black' ]
                       weight of text label font

        *text_xpad*: float
                     padding from x-axis, as a percentage of the xmin,
                     *default* is 0.95

        *text_ypad*: float
                     padding from y-axis, as a percentage of the ymin,
                     *default* is 0.95

        *xtickspace*: int
                      integer telling at what interval to place station names
                      as the tick labels, in case they are closely spaced.
                      *default* is 1

    ================ ==========================================================
    Attributes        Description
    ================ ==========================================================
        ax_pxx       axes instance for phase_xx
        ax_rxx       axes instance for res_xx
        ax_pxy       axes instance for phase_xy
        ax_rxy       axes instance for res_xy
        ax_pyx       axes instance for phase_yx
        ax_ryx       axes instance for res_yx
        ax_pyy       axes instance for phase_yy
        ax_ryy       axes instance for res_yx
        cbaxp        axes instance for phase colorbar
        cbaxr        axes instance for resistivity colorbar
        cbp          colorbar instance for phase
        cbr          colorbar instance for resistivity
        fig          figure instance for the plot
        mt_list       list of mtpy.mtplottools.MTplot instances
        mt_list_sort  same as mt_list but sorted by location assuming linedir
        offset_list   list of station offsetes assuming linedir
        phasexx      np.ndarray of phase_xx values
        phasexy      np.ndarray of phase_xy values
        phaseyx      np.ndarray of phase_yx values
        phaseyy      np.ndarray of phase_yy values

        resxx        np.ndarray of res_xx values
        resxy        np.ndarray of res_xy values
        resyx        np.ndarray of res_yx values
        resyy        np.ndarray of res_yy values
        station_list  list of stations corresponding to offset_list
        text         text instance for subplot label
    ================ ==========================================================

    Methods:
    ---------
        * *plot*: plots the pseudosection according to keywords
        * *redraw_plot*: redraws the plot, use if you change some of the
                         attributes.
        * *update_plot*: updates the plot, use if you change some of the
                         axes attributes, figure needs to be open to update.
        * *save_plot*: saves the plot to given filepath.

    """

    def __init__(self, **kwargs):
        """
        Initialize parameters
        """

        #read in key word arguments and set defaults if none given
        fn_list = kwargs.pop('fn_list', None)
        res_object_list = kwargs.pop('res_object_list', None)
        z_object_list = kwargs.pop('z_object_list', None)
        mt_object_list = kwargs.pop('mt_object_list', None)

        #--> get the inputs into a list of mt objects
        self.mt_list = mtpl.get_mtlist(fn_list=fn_list,
                                     res_object_list=res_object_list,
                                     z_object_list=z_object_list,
                                     mt_object_list=mt_object_list)

        #--> set figure parameters
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [8, 4])
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.font_size = kwargs.pop('font_size', 7)
        self.aspect = kwargs.pop('aspect', 'auto')

        self.xtickspace = kwargs.pop('xtickspace', 1)
        self.ftol = kwargs.pop('ftol', 0.1)
        self.stationid = kwargs.pop('stationid', [0,4])
        self.linedir = kwargs.pop('linedir', 'ew')

        #--> set plots to plot and how to plot them
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.plot_xx = kwargs.pop('plot_xx', 'n')
        self.plot_xy = kwargs.pop('plot_xy', 'y')
        self.plot_yx = kwargs.pop('plot_yx', 'y')
        self.plot_yy = kwargs.pop('plot_yy', 'n')
        self.plot_style = kwargs.pop('plot_style', 'imshow')
        self.imshow_interp = kwargs.pop('imshow_interp', 'bicubic')
        self.plot_period = kwargs.pop('plot_period', None)
        self.shift_yx_phase = kwargs.pop('shift_yx_phase',False)

        #--> set plot limits
        self.res_limits = kwargs.pop('res_limits', (0, 3))
        self.phase_limits = kwargs.pop('phase_limits', (0, 90))
        self.period_limits = kwargs.pop('period_limits', None)

        #--> set colorbar properties
        self.cb_pad = kwargs.pop('cb_pad', .0375)
        self.cb_orientation = kwargs.pop('cb_orientation', 'vertical')
        self.cb_shrink = kwargs.pop('cb_shrink', .75)
        self.cb_position = kwargs.pop('cb_position', None)

        #--> set text box parameters
        self.text_location = kwargs.pop('text_location', None)
        self.text_xpad = kwargs.pop('text_xpad', .95)
        self.text_ypad = kwargs.pop('text_ypad', .95)
        self.text_size = kwargs.pop('text_size', self.font_size)
        self.text_weight = kwargs.pop('text_weight', 'bold')
        self.station_label_rotation = kwargs.pop('station_label_rotation',0)
        self.show_grid = kwargs.pop('show_grid',True)
    
        #--> set colormaps Note only mtcolors is supported
        self.res_cmap = kwargs.pop('res_cmap', mtcl.cmapdict['mt_rd2gr2bl'])
        self.phase_cmap = kwargs.pop('phase_cmap', mtcl.cmapdict['mt_bl2gr2rd'])

        #create empty lists to put things into
        self.stationlist = []
        self.offsetlist = []

        #make a list of periods from each station assuming the longest one
        #is the most complete ==> largest range.
        period_list = np.array([len(mt.period) for mt in self.mt_list])

        #find index where the longest period is if multiple pick the first one
        max_find = np.where(period_list==period_list.max())[0]
        if len(max_find)>0:
            max_find = max_find[0]

        if self.plot_period is None:
            self.plot_period = \
                  self.mt_list[max_find].period

        #create empty arrays to put data into
        ns = len(self.mt_list)
        nt = len(self.plot_period)

        self.resxx = np.zeros((nt, ns))
        self.resxy = np.zeros((nt, ns))
        self.resyx = np.zeros((nt, ns))
        self.resyy = np.zeros((nt, ns))

        self.phasexx = np.zeros((nt, ns))
        self.phasexy = np.zeros((nt, ns))
        self.phaseyx = np.zeros((nt, ns))
        self.phaseyy = np.zeros((nt, ns))

        rot_z = kwargs.pop('rot_z', 0)
        # if rotation angle is an int or float make an array the length of
        # mt_list for plotting purposes
        if isinstance(rot_z, float) or isinstance(rot_z, int):
            self.rot_z = np.array([rot_z] * len(self.mt_list))

        # if the rotation angle is an array for rotation of different
        # freq than repeat that rotation array to the len(mt_list)
        elif isinstance(rot_z, np.ndarray):
            if rot_z.shape[0] != len(self.mt_list):
                self.rot_z = np.repeat(rot_z, len(self.mt_list))

        else:
            self.rot_z = rot_z

        if self.plot_yn == 'y':
            self.plot()

    #---need to rotate data on setting rotz
    def _set_rot_z(self, rot_z):
        """
        need to rotate data when setting z
        """

        # if rotation angle is an int or float make an array the length of
        # mt_list for plotting purposes
        if isinstance(rot_z, float) or isinstance(rot_z, int):
            rot_z = np.array([rot_z] * len(self.mt_list))

        # if the rotation angle is an array for rotation of different
        # freq than repeat that rotation array to the len(mt_list)
        elif isinstance(rot_z, np.ndarray):
            if rot_z.shape[0] != len(self.mt_list):
                rot_z = np.repeat(rot_z, len(self.mt_list))

        else:
            pass

        #rotate the data
        for ii,mt in enumerate(self.mt_list):
            mt.rot_z = rot_z[ii]

    rot_z = property(fset=_set_rot_z, doc="rotation angle(s)")


    def sort_by_offsets(self):
        """
        get list of offsets to sort the mt list

        """

        dtype = [('station', 'S10'), ('offset', float), ('spot', int)]
        slist = []
        #get offsets
        for ii, mt in enumerate(self.mt_list):
            #get offsets between stations
            if ii == 0:
                east0 = mt.lon
                north0 = mt.lat
                offset = 0.0
            else:
                east = mt.lon
                north = mt.lat
                #if line is predominantly e-w
                if self.linedir=='ew':
                    if east0 < east:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif east0 > east:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset = 0
                #if line is predominantly n-s
                elif self.linedir == 'ns':
                    if north0 < north:
                        offset = np.sqrt((east0-east)**2+(north0-north)**2)
                    elif north0 > north:
                        offset = -1*np.sqrt((east0-east)**2+(north0-north)**2)
                    else:
                        offset=0
            #append values to list for sorting
            slist.append((mt.station, offset, ii))

        #create a structured array according to the data type and values
        v_array = np.array(slist, dtype=dtype)

        #sort the structured array by offsets
        sorted_array = np.sort(v_array, order=['offset'])

        #create an offset list as an attribute
        self.offset_list = np.array([ss[1] for ss in sorted_array])

        #create a station list as an attribute
        self.station_list = np.array([ss[0][self.stationid[0]:self.stationid[1]]
                                     for ss in sorted_array])

        #create an index list of the sorted index values
        index_list = [ss[2] for ss in sorted_array]

        #create a new mt_list according to the offsets from the new index_list
        new_mt_list = [self.mt_list[ii] for ii in index_list]

        #set the mt_list attribute as the new sorted mt_list
        self.mt_list_sort = new_mt_list

    def get_rp_arrays(self):
        """
        get resistivity and phase values in the correct order according to
        offsets and periods.

        """

        self.sort_by_offsets()

        #create empty arrays to put data into need to reset to zero in case
        #something has changed
        ns = len(self.mt_list)
        nt = len(self.plot_period)

        self.resxx = np.zeros((nt, ns))
        self.resxy = np.zeros((nt, ns))
        self.resyx = np.zeros((nt, ns))
        self.resyy = np.zeros((nt, ns))

        self.phasexx = np.zeros((nt, ns))
        self.phasexy = np.zeros((nt, ns))
        self.phaseyx = np.zeros((nt, ns))
        self.phaseyy = np.zeros((nt, ns))

        #make a dictionary of the periods to plot for a reference
        period_dict = dict([(key, vv)
                             for vv, key in enumerate(self.plot_period)])

        for ii, mt in enumerate(self.mt_list_sort):
            #get resisitivity and phase in a dictionary and append to a list
#            rp = mt.get_ResPhase()
            rp = mt.Z

            for rr, rper in enumerate(self.plot_period):
                jj = None
                for kk, iper in enumerate(mt.period):
                    if iper == rper:
                        jj = period_dict[rper]
                        self.resxx[jj, ii] = np.log10(rp.res_xx[kk])
                        self.resxy[jj, ii] = np.log10(rp.res_xy[kk])
                        self.resyx[jj, ii] = np.log10(rp.res_yx[kk])
                        self.resyy[jj, ii] = np.log10(rp.res_yy[kk])

                        self.phasexx[jj, ii] = rp.phase_xx[kk]
                        self.phasexy[jj, ii] = rp.phase_xy[kk]
                        self.phaseyx[jj, ii] = rp.phase_yx[kk]
                        self.phaseyy[jj, ii] = rp.phase_yy[kk]

                        break

                    elif rper*(1-self.ftol) <= iper and \
                         iper <= rper*(1+self.ftol):
                             jj = period_dict[rper]
                             self.resxx[jj, ii] = np.log10(rp.res_xx[kk])
                             self.resxy[jj, ii] = np.log10(rp.res_xy[kk])
                             self.resyx[jj, ii] = np.log10(rp.res_yx[kk])
                             self.resyy[jj, ii] = np.log10(rp.res_yy[kk])

                             self.phasexx[jj, ii] = rp.phase_xx[kk]
                             self.phasexy[jj, ii] = rp.phase_xy[kk]
                             self.phaseyx[jj, ii] = rp.phase_yx[kk]
                             self.phaseyy[jj, ii] = rp.phase_yy[kk]

                             break
                    else:
                        pass

                if jj is None:
                    print('did not find period {0:.6g} (s) for {1}'.format(
                               rper, self.station_list[ii]))

    def plot(self, show=True, get_rp_arrays=True):

        #--> set subplot spacing
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .12
        plt.rcParams['figure.subplot.right'] = .90
        plt.rcParams['figure.subplot.bottom'] = .09
        plt.rcParams['figure.subplot.top'] = .98

        #get apparent resistivity and phase
        if get_rp_arrays:
            self.get_rp_arrays()
            if self.shift_yx_phase:
                self.phaseyx = self.phaseyx + 180

        #make a list of tuples to see how many subplots are needed
        ynlist = [self.plot_xx+'xx', self.plot_xy+'xy', self.plot_yx+'yx',
                 self.plot_yy+'yy']
        reslist = [self.resxx, self.resxy, self.resyx, self.resyy]
        phaselist = [self.phasexx, self.phasexy, self.phaseyx, self.phaseyy]
        plist = [(yn[1:], res, phase) for yn, res, phase in zip(ynlist,
                                                           reslist,
                                                           phaselist)
                                                           if yn[0]=='y']

        #make a general subplot array
        gs = gridspec.GridSpec(2, len(plist),
                               height_ratios=[1, 1],
                               hspace=.00,
                               wspace=.025)

        # get ylimits for plot
        if self.period_limits is None:
            self.period_limits = (self.plot_period.min(),
                                 self.plot_period.max())

        font_dict = {'size':self.font_size+2, 'weight':'bold'}
        ns = len(self.station_list)
        #--> plot data
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        #plot as a mesh where the data are left as blocks
        if self.plot_style == 'pcolormesh':
            #need to add another element at the end of the array so pcolor
            #will plot the full array
            # first, get median station spacing
            mss = np.median(np.abs(self.offset_list[1:] - self.offset_list[:-1]))
            xgrid_edges = np.mean([self.offset_list[1:],self.offset_list[:-1]],axis=0)
            xgrid = np.hstack([self.offset_list[:1],xgrid_edges,self.offset_list[-1:]])
            
            ygrid_edges = 10**np.mean([np.log10(self.plot_period[1:]),np.log10(self.plot_period[:-1])],axis=0)
            ygrid = np.hstack([self.plot_period[:1],ygrid_edges,self.plot_period[-1:]])
            xgrid, ygrid = np.meshgrid(xgrid,ygrid)
            
#            xgrid, ygrid = np.meshgrid(np.append(self.offset_list,
#                                                 self.offset_list[-1]*1.1),
#                                       np.append(self.plot_period,
#                                                 self.plot_period[-1]*1.1))

            for ii, tt in enumerate(plist):
                axr = self.fig.add_subplot(gs[0, ii])
                axp = self.fig.add_subplot(gs[1, ii])

                #plot apparent resistivity
                axr.pcolormesh(xgrid, ygrid, tt[1],#np.flipud(tt[1]),
                               cmap=self.res_cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
                
                axr.set_aspect(self.aspect)
                axp.set_aspect(self.aspect)

                axr.set_xticks(self.offset_list[list(range(0,ns,self.xtickspace))])
                if self.xtickspace != 1:
                    axr.set_xticks(self.offset_list)#, minor=True)
                plt.setp(axr.get_xticklabels(), visible=False)
                if self.show_grid:
                    axr.grid(which='major', alpha=.25)
                axr.set_yscale('log', nonposy='clip')
                axr.set_xlim(
                    self.offset_list.min(),
                    self.offset_list.max() * 1.1)
                axr.set_ylim(self.period_limits)

                #label the plot with a text box
                if self.text_location is None:
                    txloc = self.offset_list.min()*self.text_xpad
                    tyloc = self.period_limits[1]*self.text_ypad
                else:
                    txloc = self.text_location[0]
                    tyloc = self.text_location[1]

                self.text = axr.text(txloc,
                                     tyloc,
                                     '$Z_{'+tt[0]+'}$',
                                     fontdict={'size':self.text_size,
                                               'weight':self.text_weight},
                                     verticalalignment='top',
                                     horizontalalignment='left',
                                     bbox={'facecolor':'white', 'alpha':.5}
                                     )

                #plot phase
                axp.pcolormesh(xgrid, ygrid, tt[2],#np.flipud(tt[2]),
                               cmap=self.phase_cmap,
                               vmin=self.phase_limits[0],
                               vmax=self.phase_limits[1])
                if self.show_grid:
                    axp.grid(which='major', alpha=.25)
                axp.set_xticks(self.offset_list[list(range(0,ns,self.xtickspace))])
                axp.set_xticklabels([self.station_list[st]
                                    for st in range(0,ns,self.xtickspace)],
                                    rotation=self.station_label_rotation,
                                    fontsize=self.text_size)
                if self.xtickspace != 1:
                    axp.set_xticks(self.offset_list, minor=True)
                axp.set_yscale('log', nonposy='clip')
                axp.set_xlim(
                    self.offset_list.min(),
                    self.offset_list.max() * 1.1)
                axp.set_ylim(self.period_limits)
                if ii == 0:
                    axp.set_ylabel('Period (s)', font_dict)
                    axr.set_ylabel('Period (s)', font_dict)

                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)

                #add colorbars
                if ii == len(plist)-1:
                    cminr = self.res_limits[0]
                    cmaxr = self.res_limits[1]
                    #add colorbar for res
                    axrpos = axr.get_position()

                    #set position just to the right of the figure
                    if self.cb_position is None:
                        cbr_position = (axrpos.bounds[0]+axrpos.bounds[2]+\
                                        self.cb_pad,
                                        axrpos.bounds[1]+.05,
                                        .015,
                                        axrpos.bounds[3]*self.cb_shrink)
                    else:
                        cbr_position = self.cb_position[0]

                    self.cbaxr = self.fig.add_axes(cbr_position)
                    self.cbr = mcb.ColorbarBase(self.cbaxr,
                                            cmap=self.res_cmap,
                                            norm=colors.Normalize(vmin=cminr,
                                                                  vmax=cmaxr),
                                            orientation=self.cb_orientation)
                    tkrmin = np.ceil(cminr)
                    tkrmax = np.floor(cmaxr)

                    self.cbr.set_ticks(np.arange(tkrmin, tkrmax+1))
                    cbr_ticklabels = [mtpl.labeldict[ll]
                                      for ll in np.arange(tkrmin, tkrmax+1)]

                    self.cbr.set_ticklabels(cbr_ticklabels)
                    self.cbr.ax.yaxis.set_label_position('right')
                    self.cbr.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbr.ax.yaxis.tick_left()
                    self.cbr.ax.tick_params(axis='y', direction='in', pad=1)
                    self.cbr.set_label('App. Res ($\Omega \cdot$m)',
                                        fontdict={'size':self.font_size})

                    #--> add colorbar for phase
                    cminp = self.phase_limits[0]
                    cmaxp = self.phase_limits[1]

                    axppos = axp.get_position()

                    #set position just to the right of the figure
                    if self.cb_position is None:
                        cbp_position = (axppos.bounds[0]+axppos.bounds[2]+\
                                        self.cb_pad,
                                        axppos.bounds[1]+.05,
                                        .015,
                                        axppos.bounds[3]*self.cb_shrink)
                    else:
                        cbp_position = self.cb_position[1]

                    self.cbaxp = self.fig.add_axes(cbp_position)
                    self.cbp = mcb.ColorbarBase(self.cbaxp,
                                            cmap=self.phase_cmap,
                                            norm=colors.Normalize(vmin=cminp,
                                                                  vmax=cmaxp),
                                            orientation=self.cb_orientation)
                    self.cbp.set_ticks([cminp, (cmaxp-cminp)/2, cmaxp])
                    self.cbp.set_ticklabels(['{0:.0f}'.format(cminp),
                                             '{0:.0f}'.format((cmaxp-cminp)/2),
                                             '{0:.0f}'.format(cmaxp)])
                    self.cbp.ax.yaxis.set_label_position('right')
                    self.cbp.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbp.ax.yaxis.tick_left()
                    self.cbp.ax.tick_params(axis='y', direction='in', pad=.5)
                    self.cbp.set_label('Phase (deg)',
                                       fontdict={'size':self.font_size})

                #make axes attributes for user editing
                if tt == 'xx':
                    self.ax_rxx = axr
                    self.ax_pxx = axp
                elif tt == 'xy':
                    self.ax_rxy = axr
                    self.ax_pxy = axp
                elif tt == 'yx':
                    self.ax_ryx = axr
                    self.ax_pyx = axp
                elif tt == 'yy':
                    self.ax_ryy = axr
                    self.ax_pyy = axp
            if show:
                plt.show()

        #plot data as an image which can have interpolation
        elif self.plot_style == 'imshow':
            #make ticks simulate a log scale in the y direction
            #--> set major and minor ticks with appropriate labels
            major_yticks = np.arange(np.ceil(np.log10(self.period_limits[0]) if self.period_limits[0]!=0 else 0),
                                     np.floor(np.log10(self.period_limits[1]) if self.period_limits[0]!=0 else 0)+1)

            #make minor ticks look like they are on a log scale
            minor_yticks = []
            for ll in major_yticks:
                minor_yticks += [np.arange(1,10)*10**ll]
            minor_yticks = np.array(minor_yticks)
            minor_yticks = np.log10(minor_yticks.flatten())

            #set ticklabels as 10**
            yticklabels = [mtpl.labeldict[ll] for ll in major_yticks]

            for ii, tt in enumerate(plist):
                axr = self.fig.add_subplot(gs[0, ii])
                axp = self.fig.add_subplot(gs[1, ii])

                #plot apparent resistivity
                axr.imshow(tt[1],
                           cmap=self.res_cmap,
                           vmin=self.res_limits[0],
                           vmax=self.res_limits[1],
                           aspect=self.aspect,
                           interpolation=self.imshow_interp,
                           extent=(self.offset_list.min(),
                                   self.offset_list.max(),
                                   np.log10(self.plot_period.min()),
                                   np.log10(self.plot_period.max())))

                #set x ticks but remove labels
                axr.set_xticks(self.offset_list[list(range(0,ns,self.xtickspace))])
                if self.xtickspace != 1:
                    axr.set_xticks(self.offset_list, minor=True)

                plt.setp(axr.get_xticklabels(), visible=False)

                #set y-axis ticks
                axr.yaxis.set_ticks(major_yticks)
                axr.yaxis.set_ticks(minor_yticks, minor=True)
                axr.set_yticklabels(yticklabels[::-1])

                axr.grid(which='major', alpha=.25)
                axr.set_xlim(self.offset_list.min(), self.offset_list.max())
                axr.set_ylim(np.log10(self.period_limits[0]) if self.period_limits[0]!=0 else 0,
                             np.log10(self.period_limits[1]) if self.period_limits[0]!=0 else 0)

                #label the plot with a text box
                if self.text_location is None:
                    txloc = self.offset_list.min()*self.text_xpad
                    tyloc = np.log10(self.period_limits[1]*self.text_ypad) if self.period_limits[0]!=0 else 0
                else:
                    txloc = self.text_location[0]
                    tyloc = self.text_location[1]

                self.text = axr.text(txloc,
                                     tyloc,
                                     '$Z_{'+tt[0]+'}$',
                                     fontdict={'size':self.text_size,
                                               'weight':self.text_weight},
                                     verticalalignment='top',
                                     horizontalalignment='left',
                                     bbox={'facecolor':'white', 'alpha':.5})

                if ii == 0:
                    axr.set_ylabel('Period (s)', font_dict)

                #plot phase
                axp.imshow(tt[2],
                           cmap=self.phase_cmap,
                           vmin=self.phase_limits[0],
                           vmax=self.phase_limits[1],
                           aspect=self.aspect,
                           interpolation=self.imshow_interp,
                           extent=(self.offset_list.min(),
                                   self.offset_list.max(),
                                   np.log10(self.plot_period.min()),
                                   np.log10(self.plot_period.max())))

                axp.grid(which='major', alpha=.25)
                axp.set_xticks(self.offset_list[list(range(0,ns,self.xtickspace))])
                axp.set_xticklabels([self.station_list[st]
                                    for st in range(0,ns,self.xtickspace)])
                if self.xtickspace != 1:
                    axp.set_xticks(self.offset_list, minor=True)

                #remove tick labels if not the first subplot
                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)

                #set y-axis ticks
                axp.yaxis.set_ticks(major_yticks)
                axp.yaxis.set_ticks(minor_yticks, minor=True)
                axp.set_yticklabels(yticklabels[::-1])

                axp.set_xlim(self.offset_list.min(), self.offset_list.max())
                axp.set_ylim(np.log10(self.period_limits[0] if self.period_limits[0]!=0 else 0),
                             np.log10(self.period_limits[1]) if self.period_limits[0]!=0 else 0)

                if ii == 0:
                    axp.set_ylabel('Period (s)', font_dict)

                #add colorbars
                if ii == len(plist)-1:
                    cminr = self.res_limits[0]
                    cmaxr = self.res_limits[1]
                    #add colorbar for res
                    axrpos = axr.get_position()

                    #set position just to the right of the figure
                    if self.cb_position is None:
                        cbr_position = (axrpos.bounds[0]+axrpos.bounds[2]+\
                                        self.cb_pad,
                                        axrpos.bounds[1]+.05,
                                        .015,
                                        axrpos.bounds[3]*self.cb_shrink)
                    else:
                        cbr_position = self.cb_position[0]

                    self.cbaxr = self.fig.add_axes(cbr_position)
                    self.cbr = mcb.ColorbarBase(self.cbaxr,
                                            cmap=self.res_cmap,
                                            norm=colors.Normalize(vmin=cminr,
                                                                  vmax=cmaxr),
                                            orientation=self.cb_orientation)
                    tkrmin = np.ceil(cminr)
                    tkrmax = np.floor(cmaxr)

                    self.cbr.set_ticks(np.arange(tkrmin, tkrmax+1))
                    cbr_ticklabels = [mtpl.labeldict[ll]
                                      for ll in np.arange(tkrmin, tkrmax+1)]

                    self.cbr.set_ticklabels(cbr_ticklabels)
                    self.cbr.ax.yaxis.set_label_position('right')
                    self.cbr.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbr.ax.yaxis.tick_left()
                    self.cbr.ax.tick_params(axis='y', direction='in', pad=1)
                    self.cbr.set_label('App. Res ($\Omega \cdot$m)',
                                        fontdict={'size':self.font_size})

                    #--> add colorbar for phase
                    cminp = self.phase_limits[0]
                    cmaxp = self.phase_limits[1]

                    axppos = axp.get_position()

                    #set position just to the right of the figure
                    if self.cb_position is None:
                        cbp_position = (axppos.bounds[0]+axppos.bounds[2]+\
                                        self.cb_pad,
                                        axppos.bounds[1]+.05,
                                        .015,
                                        axppos.bounds[3]*self.cb_shrink)
                    else:
                        cbp_position = self.cb_position[1]

                    self.cbaxp = self.fig.add_axes(cbp_position)
                    self.cbp = mcb.ColorbarBase(self.cbaxp,
                                            cmap=self.phase_cmap,
                                            norm=colors.Normalize(vmin=cminp,
                                                                  vmax=cmaxp),
                                            orientation=self.cb_orientation)
                    self.cbp.set_ticks([cminp, (cmaxp-cminp)/2, cmaxp])
                    self.cbp.set_ticklabels(['{0:.0f}'.format(cminp),
                                             '{0:.0f}'.format((cmaxp-cminp)/2),
                                             '{0:.0f}'.format(cmaxp)])
                    self.cbp.ax.yaxis.set_label_position('right')
                    self.cbp.ax.yaxis.set_label_coords(1.35, .5)
                    self.cbp.ax.yaxis.tick_left()
                    self.cbp.ax.tick_params(axis='y', direction='in', pad=.5)
                    self.cbp.set_label('Phase (deg)',
                                       fontdict={'size':self.font_size})

                if tt[0] == 'xx':
                    self.ax_rxx = axr
                    self.ax_pxx = axp
                elif tt[0] == 'xy':
                    self.ax_rxy = axr
                    self.ax_pxy = axp
                elif tt[0] == 'yx':
                    self.ax_ryx = axr
                    self.ax_pyx = axp
                elif tt[0] == 'yy':
                    self.ax_ryy = axr
                    self.ax_pyy = axp
            if show:
                plt.show()

    def save_plot(self, save_fn, file_format='pdf', orientation='portrait',
                  fig_dpi=None, close_plot='y'):
        """
        save_plot will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as
                            save_fn/station_name_ResPhase.file_format

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
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation)
            # plt.clf()
            # plt.close(self.fig)

        else:
            save_fn = os.path.join(save_fn,
                                   self._mt.station+'_ResPhasePseudoSection.'+
                                    file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation)

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: '+self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> [ax.grid(True, which='major') for ax in [p1.axr,p1.axp]]
            >>> p1.update_plot()

        """

        self.fig.canvas.draw()

    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.xy_color = (.5,.5,.9)
            >>> p1.xy_marker = '*'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()


    def writeTextFiles(self, save_path=None, ptol=0.10):
        """
        This will write text files for all the phase tensor parameters
        """

        if save_path is None:
            try:
                svpath = os.path.dirname(self.mt_list[0].fn)
            except TypeError:
                raise IOError('Need to input save_path, could not find a path')
        else:
            svpath = save_path

        if self.resxy.mean() == 0 :
            self.get_rp_arrays()

        header_list = ['{0:^10}'.format('period(s)')]+\
                     ['{0:^8}'.format(ss) for ss in self.station_list]+['\n']

        fn_dict = {'resxx':self.resxx,
                   'resxy':self.resxy,
                   'resyx':self.resyx,
                   'resyy':self.resyy,
                   'phasexx':self.phasexx,
                   'phasexy':self.phasexy,
                   'phaseyx':self.phaseyx,
                   'phaseyy':self.phaseyy}

        #write the arrays into lines properly formatted
        t1_kwargs = {'spacing':'{0:^10} ', 'value_format':'{0:.2e}',
                     'append':False, 'add':False}

        tr_kwargs = {'spacing':'{0:^8}', 'value_format':'{0: .2f}',
                     'append':False, 'add':False}

        tp_kwargs = {'spacing':'{0:^8}', 'value_format':'{0: .2f}',
                     'append':False, 'add':False}


        for key in list(fn_dict.keys()):
            fid = file(os.path.join(svpath, 'PseudoSection.'+key), 'w')
            fid.write(''.join(header_list))
            for ii, per in enumerate(self.plot_period):
                if key[0] == 'r':
                    line = [mtpl.make_value_str(per, **t1_kwargs)] + \
                           [mtpl.make_value_str(rr, **tr_kwargs)
                            for rr in fn_dict[key][ii]]+['\n']
                elif key[0] == 'p':
                    line = [mtpl.make_value_str(per, **t1_kwargs)] + \
                           [mtpl.make_value_str(rr, **tp_kwargs)
                            for rr in fn_dict[key][ii]]+['\n']
                fid.write(''.join(line))
            fid.close()

        print('Wrote files to: '+\
                        os.path.join(svpath, 'PseudoSection.component'))

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots Resistivity and phase as a pseudo section."


