#!/bin/env python
"""
Description:
    Plots resistivity and phase maps for a given frequency
   
References:
 
CreationDate:   4/19/18
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     4/19/18   RH

"""

import matplotlib.pyplot as plt
import numpy as np
import os, glob
from matplotlib.ticker import FormatStrFormatter
import mtpy.utils.gis_tools as gis_tools
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import mtpy.imaging.mtcolors as mtcl
import mtpy.imaging.mtplottools as mtpl
from mtpy.utils.mtpylog import MtPyLog
import mtpy.analysis.pt as MTpt
import matplotlib.tri as tri
from scipy.spatial import cKDTree
from scipy.spatial import Delaunay
from matplotlib import ticker
from matplotlib import colors


class PlotResPhaseMaps(mtpl.PlotSettings):
    """
    Plots apparent resistivity and phase in map view from a list of edi files

    Arguments:
    -------------

        **fn_list** : list of strings
                          full paths to .edi files to plot

        **fig_size** : tuple or list (x, y) in inches
                      dimensions of the figure box in inches, this is a default
                      unit of matplotlib.  You can use this so make the plot
                      fit the figure box to minimize spaces from the plot axes
                      to the figure box.  *default* is [8, 8]

        **mapscale** : [ 'deg' | 'm' | 'km' ]
                       Scale of the map coordinates.

                       * 'deg' --> degrees in latitude and longitude

                       * 'm' --> meters for easting and northing

                       * 'km' --> kilometers for easting and northing

        **plot_yn** : [ 'y' | 'n' ]
                      *'y' to plot on creating an instance

                      *'n' to not plot on creating an instance

        **title** : string
                    figure title

        **dpi** : int
                  dots per inch of the resolution. *default* is 300

        **font_size** : float
                        size of the font that labels the plot, 2 will be added
                        to this number for the axis labels.
    """

    def __init__(self, **kwargs):
        """
        Initialise the object
        :param kwargs: keyword-value pairs
        """
        super(PlotResPhaseMaps, self).__init__(**kwargs)

        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)

        fn_list = kwargs.pop('fn_list', [])

        if(len(fn_list)==0): raise NameError('File list is empty.')

        # ----set attributes for the class-------------------------
        self.mt_list = mtpl.get_mtlist(fn_list=fn_list)

        # read in map scale
        self.mapscale = kwargs.pop('mapscale', 'deg')
        if self.mapscale == 'km':
            self.dscale = 1000.
        elif self.mapscale == 'm':
            self.dscale = 1.
        # end if

        self.plot_title = kwargs.pop('plot_title', None)
        self.fig_dpi = kwargs.pop('fig_dpi', 100)

        self.fig_size = kwargs.pop('fig_size', [8, 6])

        self.font_size = kwargs.pop('font_size', 7)

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.save_fn = kwargs.pop('save_fn', "/c/tmp")

        # By this stage all keyword arguments meant to be set as class properties will have
        # been processed. Popping all class properties that still exist in kwargs
        self.kwargs = kwargs
        for key in vars(self):
            self.kwargs.pop(key, None)

        self.axesList = []
    # end func

    # -----------------------------------------------
    # The main plot method for this module
    # -------------------------------------------------
    def plot(self, freq, type, vmin, vmax,
             extrapolation_buffer_degrees=1,
             regular_grid_nx=100, regular_grid_ny=100,
             nn=7,
             p = 4,
             show_stations = True,
             show_station_names = False,
             save_path = os.getcwd(),
             file_ext = 'png',
             cmap='rainbow',
             show = True):
        """
        :param freq: plot frequency
        :param type: plot type; can be either 'res' or 'phase'
        :param vmin: minimum value used in color-mapping
        :param vmax: maximum value used in color-mapping
        :param extrapolation_buffer_degrees: extrapolation buffer in degrees
        :param regular_grid_nx: number of longitudinal grid points to use during interpolation
        :param regular_grid_ny: number of latitudinal grid points to use during interpolation
        :param nn: number of nearest neighbours to use in inverse distance weighted interpolation
        :param p: power parameter in inverse distance weighted interpolation
        :param save_path: path where plot is saved
        :param file_ext: file extension
        :param show: boolean to toggle display of plot
        :return: fig object
        """

        if(type not in ['res', 'phase']): raise NameError("type must be 'res' or 'phase'")
        if(not os.path.isdir(save_path)): raise NameError("Invalid save_path")

        def in_hull(p, hull):
            """
            Test if points in p are within the convex hull
            """

            try:
                if not isinstance(hull, Delaunay):
                    hull = Delaunay(hull)

                return hull.find_simplex(p)>=0
            except:
                from scipy.optimize import linprog

                # Delaunay triangulation will fail if there are collinear points; in those instances
                # use linear programming (much slower) to define a convex hull.
                def in_hull_lp(points, x):
                    """
                    :param points:
                    :param x:
                    :return:
                    """
                    n_points = len(points)
                    n_dim = len(x)
                    c = np.zeros(n_points)
                    A = np.r_[points.T, np.ones((1, n_points))]
                    b = np.r_[x, np.ones(1)]
                    lp = linprog(c, A_eq=A, b_eq=b)
                    return not lp.success
                # end func

                result = []
                for cp in p:
                    result.append(in_hull_lp(hull, cp))
                # end for

                return np.array(result)
            # end try
        # end func
        
        # change vmin, vmax to 2x2 array
        if not np.iterable(vmin):
            vmin = np.ones((2,2))*vmin
        if not np.iterable(vmax):
            vmax = np.ones((2,2))*vmax       
        

        # set position properties for the plot
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = .1
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .93
        plt.rcParams['figure.subplot.wspace'] = .55
        plt.rcParams['figure.subplot.hspace'] = .70

        # make figure instance
        self.fig = plt.figure(1, figsize=self.fig_size, dpi=self.fig_dpi)

        # clear the figure if there is already one up
        plt.clf()

        # interpolate data
        res, phase, lat, lon = [], [], [], []
        for mt_obj in self.mt_list:
            z_obj_i, tipper_obj_i = mt_obj.interpolate([freq], bounds_error=False)
            z_obj_i.compute_resistivity_phase()
            res.append(z_obj_i.resistivity[0])
            phase.append(z_obj_i.phase[0])
            lat = np.append(lat, mt_obj.lat)
            lon = np.append(lon, mt_obj.lon)
        # end for

        lon = np.array(lon)
        lat = np.array(lat)
        res = np.array(res)
        phase = np.array(phase)
        phase[:,1,0] += 180

        elon = np.array(lon)
        elat = np.array(lat)

        elon[np.argmin(elon)] -= extrapolation_buffer_degrees
        elon[np.argmax(elon)] += extrapolation_buffer_degrees
        elat[np.argmin(elat)] -= extrapolation_buffer_degrees
        elat[np.argmax(elat)] += extrapolation_buffer_degrees

        x = np.zeros(lon.shape)
        y = np.zeros(lon.shape)
        ex = np.zeros(lon.shape)
        ey = np.zeros(lon.shape)

        # plot results
        insideIndices = []
        tree = None
        triangulation = None
        foundCoordinates = False
        plotIdx = 1
        for i in range(2):
            for j in range(2):
                ax = self.fig.add_subplot(2, 2, plotIdx)
                self.axesList.append(ax)

                nx = regular_grid_nx
                ny = regular_grid_ny
                if (not foundCoordinates):
                    # transform coordinates if necessary
                    if (self.mapscale == 'm' or self.mapscale=='km'):
                        zl = []
                        zle = []
                        for k in range(len(lon)):
                            east, north, zone = gis_tools.project_point_ll2utm(lat[k],
                                                                               lon[k])
                            x[k] = east / self.dscale
                            y[k] = north / self.dscale
                            zl.append(zone)

                            east, north, zone = gis_tools.project_point_ll2utm(elat[k],
                                                                               elon[k])
                            ex[k] = east / self.dscale
                            ey[k] = north / self.dscale
                            zle.append(zone)
                        # end for

                        if (len(set(zl)) > 1 or len(set(zle)) > 1):
                            print('Warning: multiple UTM zones detected. ' \
                                  'Using geographical coordinates instead')
                            x = lon
                            y = lat
                            ex = elon
                            ey = elat
                        # end if
                    else:
                        x = lon
                        y = lat
                        ex = elon
                        ey = elat
                    # end if

                    rx = np.linspace(ex.min(), ex.max(), nx)
                    ry = np.linspace(ey.min(), ey.max(), ny)
                    rx, ry = np.meshgrid(rx, ry)
                    rx = rx.flatten()
                    ry = ry.flatten()

                    triangulation = tri.Triangulation(rx, ry)

                    mx = rx[triangulation.triangles].mean(axis=1)
                    my = ry[triangulation.triangles].mean(axis=1)

                    mxmy = np.array([mx, my]).T
                    exey = np.array([ex, ey]).T

                    insideIndices = in_hull(mxmy, exey)
                    insideIndices = np.bool_(insideIndices)
                    triangulation.set_mask(~insideIndices)

                    foundCoordinates = True

                    tree = cKDTree(np.array([x, y]).T)
                # end if

                xy = np.array([rx, ry]).T
                d, l = tree.query(xy, k=nn)

                img = None
                vs = res if type == 'res' else phase

                if (nn == 1):
                    # extract nearest neighbour values
                    img = vs[:, i, j][l]
                else:
                    vals = vs[:, i, j]
                    img = np.zeros((xy.shape[0]))

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

                if(isinstance(cmap, str)):
                    cmap = plt.get_cmap(cmap)
                # set cmap values for over and under
                norm = colors.Normalize(vmin=vmin[i,j], vmax=vmax[i,j])
                cmap.set_over(cmap(norm(vmax[i,j])))
                cmap.set_under(cmap(norm(vmin[i,j])))

                if (type == 'res'):
                    # Log-normalized contour plots do not support the 'extend' keyword which
                    # can be used to clip data values above/below the given range to their
                    # corresponding colors. We do the following to get around this issue.
                    cbinfo = ax.tricontourf(triangulation, np.log10(img), mask=insideIndices,
                                            levels=np.linspace(np.log10(vmin[i,j]), np.log10(vmax[i,j]), 50),
                                            extend='both',
                                            cmap=cmap)

                    cb = self.fig.colorbar(cbinfo,
                                           ticks=ticker.FixedLocator(
                                               np.arange(int(np.round(np.log10(vmin[i,j]))),
                                                         int(np.round(np.log10(vmax[i,j])))+1)))

                    labels = ['$10^{%d}$'%l for l in
                              np.arange(int(np.round(np.log10(vmin[i,j]))), int(np.round(np.log10(vmax[i,j])))+1)]
                    cb.ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
                elif (type == 'phase'):
                    cbinfo = ax.tricontourf(triangulation, img, mask=insideIndices,
                                             levels=np.linspace(vmin[i,j], vmax[i,j], 50),
                                             norm=colors.Normalize(vmin=vmin[i,j], vmax=vmax[i,j]),
                                             extend='both',
                                             cmap=cmap)

                    cb = self.fig.colorbar(cbinfo, ticks=np.linspace(vmin[i,j], vmax[i,j], 12))
                # end if

                ax.tick_params(axis='both', which='major', labelsize=self.font_size-2)
                ax.tick_params(axis='both', which='minor', labelsize=self.font_size-2)

                cb.ax.tick_params(axis='both', which='major', labelsize=self.font_size-1)
                cb.ax.tick_params(axis='both', which='minor', labelsize=self.font_size-1)

                # show stations
                if (show_stations): 
                    ax.scatter(x, y, 2, marker='v', c='k', edgecolor='none')
                    if show_station_names:
                        for isn, mt_obj in enumerate(self.mt_list):
                            plt.text(lon[isn],lat[isn],mt_obj.station,fontsize=self.font_size-2)

                # Label plots
                label = ''
                if(i==0 and j==0):
                    if(type=='res'):
                        label = '$\\rho_{xx}  \\mathrm{[\Omega m]}$'
                    else:
                        label = '$\\phi_{xx} \\mathrm{[^\circ]}$'
                elif(i==0 and j==1):
                    if(type=='res'):
                        label = '$\\rho_{xy}$'
                    else:
                        label = '$\\phi_{xy}$'
                elif(i==1 and j==0):
                    if(type=='res'):
                        label = '$\\rho_{yx}$'
                    else:
                        label = '$\\phi_{yx}$'
                elif(i==1 and j==1):
                    if(type=='res'):
                        label = '$\\rho_{yy}$'
                    else:
                        label = '$\\phi_{yy}$'


                ax.text(0.8, 0.9, label, fontdict={'size': self.font_size + 3},
                        transform=ax.transAxes)
                plotIdx += 1
            # end for
        # end for

        # Plot title
        suffix = ' %0.2f Hz'%(freq) if (freq>=1) else ' %0.2f s'%(1./freq)
        if(type=='res'):
            self.fig.suptitle('Apparent Resistivity Maps for'+suffix, y=0.985)
        else:
            self.fig.suptitle('Phase Maps for'+suffix, y=0.985)

        plt.tight_layout(rect=[0, 0.025, 1, 0.975])
        if (show): plt.show()

        fn = os.path.join(save_path, '%s.%0.2f.%s'%(type, freq, file_ext))
        self.fig.savefig(fn, dpi=self.fig_dpi)

        return self.fig
    # end func
# end class


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    imaging = os.path.dirname(__file__)
    mtpy = os.path.dirname(imaging)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    edidir = os.path.join(data, 'edi2')

    edi_file_list = glob.glob(edidir + '/*.edi')

    prp = PlotResPhaseMaps(fn_list=edi_file_list,
                           fig_dpi=200, mapscale='m')

    plot_type = 'res'

    if(plot_type=='res'):
        f = prp.plot(0.02, plot_type, 0.005, 1e2,
                     extrapolation_buffer_degrees=0.1,
                     regular_grid_nx=100,
                     regular_grid_ny=100,
                     show=True, save_path='/tmp')
    else:
        f = prp.plot(0.02, plot_type, -180, 180,
                     extrapolation_buffer_degrees=0.1,
                     regular_grid_nx=100,
                     regular_grid_ny=100,
                     show=True, save_path='/tmp')