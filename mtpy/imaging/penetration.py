"""
    Description:
        This file defines imaging functions for penetration.
        The plotting function are extracted and implemented in plot() of each class from penetration_depth1D.py,
        penetration_depth2D.py and penetration_depth3D.py

    Usage:
        see descriptions of each clases

    Author: YingzhiGou
    Date: 20/06/2017
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

import mtpy
import mtpy.modeling.occam2d_rewrite as occam2d
from imaging_base import ImagingBase, ParameterError, ImagingError
from mtpy.core import mt as mt
from mtpy.utils.decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
logger = MtPyLog().get_mtpy_logger(__name__)
# default contains of rholist
DEFAULT_RHOLIST = {'zxy', 'zyx', 'det'}


class Depth1D(ImagingBase):
    """
        Description:
        For a given input MT object, plot the Penetration Depth vs all the periods (1/freq).
    """
    def set_data(self, data):
        # this plot only use one edi each time
        self._set_edi(data)

    def __init__(self, edis=None, rholist=DEFAULT_RHOLIST):
        super(Depth1D, self).__init__()
        self._rholist = None
        self.set_data(edis)
        self.set_rholist(rholist)

    def set_rholist(self, rholist=DEFAULT_RHOLIST):
        if not isinstance(rholist, set):
            rholist = set(rholist)
            self._reset_fig()

        if rholist.difference(DEFAULT_RHOLIST):
            # there are unsupported values
            # previous behaviour: ignore
            pass

        self._rholist = rholist
        self._reset_fig()

    def plot(self):
        if self._data is None or not self._data:
            # todo: raise an exception
            raise NotImplemented
        elif self._rholist is None or not self._rholist:
            raise ZComponentError
        elif self._fig is not None:
            # nothing to plot
            return

        self._logger.info("Plotting the edi file %s", self._data.fn)
        self._fig = plt.figure(figsize=(8, 6), dpi=80)
        self._fig.set_tight_layout(True)
        plt.grid(True)

        zeta = self._data.Z  # the attribute Z represent the impedance tensor 2X2 matrix
        freqs = zeta.freq  # frequencies
        scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

        # The periods array
        periods = 1.0 / freqs
        legendh = []

        if 'zxy' in self._rholist:
            # One of the 4-components: XY
            penetration_depth = scale_param * \
                                np.sqrt(zeta.resistivity[:, 0, 1] * periods)

            # pen_zxy, = plt.semilogx(periods, -penetration_depth, '-*',label='Zxy')
            pen_zxy, = plt.semilogx(
                periods, -penetration_depth, color='#000000', marker='*', label='Zxy')
            # See
            # http://matplotlib.org/1.3.1/examples/pylab_examples/line_styles.html

            legendh.append(pen_zxy)

        if 'zyx' in self._rholist:
            penetration_depth = scale_param * \
                                np.sqrt(zeta.resistivity[:, 1, 0] * periods)

            pen_zyx, = plt.semilogx(
                periods, -penetration_depth, color='g', marker='o', label='Zyx')
            legendh.append(pen_zyx)

        if 'det' in self._rholist:
            # determinant
            det2 = np.abs(zeta.det[0])
            det_penetration_depth = scale_param * np.sqrt(0.2 * periods * det2 * periods)

            # pen_det, = plt.semilogx(periods, -det_penetration_depth, '-^', label='Determinant')
            pen_det, = plt.semilogx(
                periods, -det_penetration_depth, color='b', marker='^', label='Determinant')
            legendh.append(pen_det)

        plt.legend(
            handles=legendh,
            bbox_to_anchor=(
                0.1,
                0.5),
            loc=3,
            ncol=1,
            borderaxespad=0.)

        title = "Penetration Depth for file %s" % self._data.fn
        plt.title(title)
        plt.xlabel("Log Period (seconds)", fontsize=16)
        plt.ylabel("Penetration Depth (meters)", fontsize=16)
        # set window title
        self._fig.canvas.set_window_title(title)


class Depth2D(ImagingBase):
    """
    With a list of MT object and a list of period index,
    generate a profile using occam2d module,
    then plot the Penetration Depth profile at the given periods vs the stations locations.
    """
    def plot(self, **kwargs):
        if self._rho is None:
            raise ZComponentError
        elif self._period_indexes is None or not self._period_indexes:
            raise Exception("please provide a period index list")
        elif self._fig is not None:
            return  # nothing to plot

        pr = occam2d.Profile(edi_list=self._data)
        pr.generate_profile()
        # pr.plot_profile(station_id=[0, 4])

        self._fig = plt.figure(figsize=(8, 6), dpi=80)
        self._fig.set_tight_layout(True)
        for period_index in self._period_indexes:
            self._logger.debug("doing period index %s", period_index)
            (stations, periods, pen, _) = get_penetration_depth(pr.edi_list, int(period_index), whichrho=self._rho)
            line_label = "Period=%s" % periods[0]

            plt.plot(
                pr.station_locations,
                pen,
                "--",
                marker="o",
                markersize=12,
                linewidth=2,
                label=line_label)
            plt.legend()

        plt.ylabel(
            'Penetration Depth (Metres) Computed by %s' %
            self._rho, fontsize=16)
        plt.yticks(fontsize=16)

        plt.xlabel('MT Penetration Depth Profile Over Stations.', fontsize=16)
        self._logger.debug("stations= %s", stations)
        self._logger.debug("station locations: %s", pr.station_locations)
        if pr.station_list is not None:
            plt.xticks(
                pr.station_locations,
                pr.station_list,
                rotation='horizontal',
                fontsize=16)
        else:  # Are the stations in the same order as the profile generated pr.station_list????
            plt.xticks(
                pr.station_locations,
                stations,
                rotation='horizontal',
                fontsize=16)

        # plt.tight_layout()
        plt.gca().xaxis.tick_top()
        self._fig.canvas.set_window_title("MT Penetration Depth Profile by %s" % self._rho)
        plt.legend(loc="best")

    def set_data(self, data):
        # this plot require multiple edi files
        self._set_edis(data)

    def __init__(self, data=None, period_index_list=None, rho='det'):
        super(Depth2D, self).__init__()
        self._period_indexes = None
        self._rho = None
        self.set_data(data)
        self.set_rho(rho)
        self.set_period_index_list(period_index_list)

    def set_rho(self, rho):
        if rho is None or rho in DEFAULT_RHOLIST:
            self._rho = rho
            self._reset_fig()
        else:
            self._logger.critical("unsupported method to compute penetratoin depth: %s", rho)
            # raise Exception("unsupported method to compute penetratoin depth: %s" % rho)

    def set_period_index_list(self, period_index_list):
        if period_index_list:
            self._period_indexes = period_index_list
        else:
            self._logger.error("Please provide a period index list like [1,2,3,4]")


class Depth3D(ImagingBase):
    """
    For a batch of MT_stations (input as a list of MT objects),
    plot the Penetration Depth vs the station_location,
    for a given period value or index (1/freq)-
    Note that the values of periods within10% tolerance (ptol=0.1) are considered as equal.
    Setting a smaller value for ptol(=0.05)may result less MT sites data included.
    """
    def __init__(self, data=None, period=None, rho='det', ptol=0.1):
        super(Depth3D, self).__init__()
        self._rho = None
        self._period = None
        self._period_fmt = None
        self._ptol = ptol
        self.set_data(data)
        self.set_rho(rho)
        self.set_period(period)

    def plot(self, **kwargs):
        if self._rho is None:
            raise ZComponentError
        elif self._period is None:
            raise ParameterError("please set period")
        elif self._fig is not None:
            # nothing to plot
            return

        period_by_index = kwargs.pop("period_by_index", False)  # search periods by its index in the data file
        plot_station_id = kwargs.pop("plot_station_id", False)  # plot the id of each station on the image
        z_unit = kwargs.pop("z_unit", 'km')

        if z_unit not in ('m', 'km'):
            raise ParameterError("z_unit has to be m or km.")

        if period_by_index:  # self._period is considered as an index
            if not isinstance(self._period, int):
                self._logger.warning("period value is not integer but used as an index.")
            (stations, periods, pendep, latlons) = get_penetration_depth(self._data,
                                                                         int(self._period),
                                                                         whichrho=self._rho)
        else:
            (stations, periods, pendep, latlons) = get_penetration_depth_generic(self._data,
                                                                                 self._period,
                                                                                 whichrho=self._rho, ptol=self._ptol)

        # create figure
        self._fig = plt.figure(figsize=(8, 6), dpi=80)
        self._fig.set_tight_layout(True)

        if check_period_values(periods) is False:
            self._logger.error("The period values are NOT equal - Please check!!! %s", periods)
            plt.plot(periods, "-^")
            title = "ERROR: Periods are NOT equal !!!"
            plt.title(title, )
            self._fig.canvas.set_window_title(title)
            raise ImagingError("Period values NOT equal across the EDI files. Please check!!!")
        else:
            # good case
            period0 = periods[0]

            if period0 < 1.0:
                # kept 4 signifiant digits - nonzero digits
                self._period_fmt = str(mtpy.utils.calculator.roundsf(period0, 4))
            else:
                self._period_fmt = "%.2f" % period0
            bbox = get_bounding_box(latlons)

            self._logger.debug("Bounding Box %s", bbox)

            xgrids = bbox[0][1] - bbox[0][0]
            ygrids = bbox[1][1] - bbox[1][0]

            self._logger.debug("xy grids: %s %s", xgrids, ygrids)

            minlat = bbox[1][0]
            minlon = bbox[0][0]

            # Pixel size in Degree:  0.001=100meters, 0.01=1KM 1deg=100KM
            pixelsize = 0.002  # Degree 0.002=200meters, 0.01=1KM 1deg=100KM

            nx = int(np.ceil(xgrids / pixelsize))
            ny = int(np.ceil(ygrids / pixelsize))

            self._logger.debug("number of grids xy: %s %s", nx, ny)

            # make the image slightly bigger than the (nx, ny) to contain all points
            # avoid index out of bound
            pad = 1  # pad = 1 affect the top and right of the plot. it is linked to get_index offset?
            # todo change this part to use xy bound offset? (0.5 gride on each side?)
            nx_padded = nx + pad
            ny_padded = ny + pad

            # fast initialization
            zdep = np.empty((ny_padded, nx_padded))
            zdep.fill(np.nan)  # initialize all pixel value as np.nan

            points = np.zeros((len(latlons), 2))
            values = np.zeros(len(latlons))

            station_points = np.zeros((len(latlons), 2))

            self._logger.debug("zdep shape %s", zdep.shape)

            for iter, (lat, lon) in enumerate(latlons):
                # self._logger.debug(iter, pair)
                (xi, yi) = get_index(lat, lon, minlat, minlon, pixelsize)
                zdep[ny_padded - yi - 1, xi] = np.abs(pendep[iter])
                # print pair
                points[iter, 0] = ny_padded - yi - 1
                points[iter, 1] = xi
                values[iter] = np.abs(pendep[iter])
                station_points[iter, 0] = ny_padded - yi - 1
                station_points[iter, 1] = xi

            # plt.imshow(zdep, interpolation='none') #OR plt.imshow(zdep,  interpolation='spline36')
            # plt.colorbar()
            # plt.show() # without this show(), the 2 figure will be plotted in one
            # canvas, overlay and compare

            # griddata interpolation of the zdep sample MT points.
            # self._logger.debug(zdep.shape)

            # grid_x, grid_y = np.mgrid[0:95:96j, 0:83:84j]  # this syntax with
            # complex step 96j has different meaning
            # this is more straight forward.
            grid_x, grid_y = np.mgrid[0:ny_padded:1, 0:nx_padded:1]

            # grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
            grid_z = griddata(points, values, (grid_x, grid_y), method='linear')
            # grid_z = griddata(points, values, (grid_x, grid_y), method='cubic')

            # method='cubic' may cause negative interp values; set them nan to make
            # empty
            grid_z[grid_z < 0] = np.nan
            if z_unit == 'km':  # change to km
                grid_z = grid_z / 1000.0

            # use reverse color map in imshow and the colorbar
            my_cmap = matplotlib.cm.jet_r
            # my_cmap = matplotlib.cm.jet
            # from mtpy.imaging.penetration_depth3d import reverse_colourmap
            # my_cmap = reverse_colourmap(my_cmap)

            # since matplotlib v2.0 the default interpolation is changed to nearest, use "bilinear" to restore the default behaviour in 1.5.3
            # imgplot = plt.imshow(grid_z, origin='upper', cmap=my_cmap, interpolation='bilinear', resample=False)
            imgplot = plt.imshow(grid_z, origin='upper', cmap=my_cmap)

            # the stations sample point 1-lon-j, 0-lat-i
            plt.plot(station_points[:, 1], station_points[:, 0], 'kv', markersize=6, )
            # add station id
            if plot_station_id:
                for label, x, y in zip(stations, station_points[:, 1], station_points[:, 0]):
                    plt.annotate(label, xy=(x, y), fontsize=9)

            # set the axix limit to control white margins
            padx = int(nx * 0.01)
            pady = int(ny * 0.01)
            min_margin = 4

            # adjusted if necessary, the number of grids extended out of the sample points area
            margin = max(padx, pady, min_margin)
            self._logger.debug("**** station_points shape ***** %s", station_points.shape)
            self._logger.debug("**** grid_z shape ***** %s", grid_z.shape)
            self._logger.debug("margin = %s" % margin)

            # horizontal axis 0-> the second index (i,j) of the matrix
            plt.xlim(-margin, grid_z.shape[1] + margin)
            # vertical axis origin at upper corner, not the lower corner.
            plt.ylim(grid_z.shape[0] + margin, -margin)

            ax = plt.gca()
            plt.gcf().set_size_inches(6, 6)

            ftsize = 14
            numticks = 5  # number of ticks to draw 5,10?
            stepx = int(zdep.shape[1] / numticks)
            stepy = int(zdep.shape[0] / numticks)
            xticks = np.arange(0, zdep.shape[1], stepx)  # 10, 100
            yticks = np.arange(0, zdep.shape[0], stepy)

            xticks_label = ['%.2f' % (bbox[0][0] + pixelsize * xtick)
                            for xtick in xticks]  # formatted float numbers
            yticks_label = ['%.2f' % (bbox[1][0] + pixelsize * ytick)
                            for ytick in yticks]

            self._logger.debug("xticks_labels= %s", xticks_label)
            self._logger.debug("yticks_labels= %s", yticks_label)
            # make sure the latitude y-axis is correctly labeled.
            yticks_label.reverse()

            plt.xticks(xticks, xticks_label, rotation='0', fontsize=ftsize)
            plt.yticks(yticks, yticks_label, rotation='horizontal', fontsize=ftsize)
            ax.set_ylabel('Latitude', fontsize=ftsize)
            ax.set_xlabel('Longitude', fontsize=ftsize)
            ax.tick_params(
                axis='both',
                which='major',
                width=2,
                length=5,
                labelsize=ftsize)
            # plt.title('Penetration Depth at the Period=%.6f (Cubic Interpolation)\n'
            # % period_fmt)  # Cubic
            title = "Penetration Depth at the Period=%s seconds \n" % self._period_fmt  # todo is \n necessory?
            plt.title(title)  # Cubic
            self._fig.canvas.set_window_title(title)

            # method-1. this is the simplest colorbar, but cannot take cmap_r
            # plt.colorbar(cmap=my_cmap_r).set_label(label='Penetration Depth
            # (Meters)', size=ftsize)  # ,weight='bold')

            # method-2 A more controlled colorbar:
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(ax)
            # pad = separation from figure to colorbar
            cax = divider.append_axes("right", size="3%", pad=0.2)
            mycb = plt.colorbar(imgplot, cax=cax, use_gridspec=True)  # cmap=my_cmap_r, does not work!!
            mycb.outline.set_linewidth(2)
            mycb.set_label(label='Penetration Depth ({})'.format(z_unit), size=ftsize)
            mycb.set_cmap(my_cmap)

    def set_data(self, data):
        # this plot need a list of edi files
        self._set_edis(data)

    def set_rho(self, rho):
        if rho is None or rho in DEFAULT_RHOLIST:
            self._rho = rho
            self._reset_fig()
        else:
            self._logger.critical("unsupported method to compute penetratoin depth: %s", rho)
            # raise Exception("unsupported method to compute penetratoin depth: %s" % rho)

    def set_period(self, period):
        if period is not None:
            self._period = period
            self._reset_fig()

    @deprecated("this function is added only to compatible with the old script penetration_depth3d.py")
    def get_period_fmt(self):
        if self._fig is not None:
            return self._period_fmt
        else:
            return None


# Utility functions (may need to move to utility module

def get_penetration_depth(mt_obj_list, per_index, whichrho='det'):
    """
    compute the penetration depth of mt_obj at the given period_index, and using whichrho option
    :param per_index: the index of periods 0, 1, ...
    :param mt_obj_list: list of edi file paths or mt objects
    :param whichrho: det, zxy, or zyx
    :return:
    """

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    # per_index=0,1,2,....
    periods = []
    pen_depth = []
    stations = []
    latlons = []
    for mt_obj in mt_obj_list:
        if isinstance(mt_obj, str) and os.path.isfile(mt_obj):
            mt_obj = mt.MT(mt_obj)
        elif not isinstance(mt_obj, mt.MT):
            raise Exception("Unsupported list of objects %s" % type(mt_obj))
        # station id
        stations.append(mt_obj.station)
        # latlons
        latlons.append((mt_obj.lat, mt_obj.lon))
        # the attribute Z
        zeta = mt_obj.Z

        if per_index >= len(zeta.freq):
            logger.debug(
                "Number of frequecies (Max per_index)= %s", len(
                    zeta.freq))
            raise Exception(
                "Index out_of_range Error: period index must be less than number of periods in zeta.freq")

        per = 1.0 / zeta.freq[per_index]
        periods.append(per)

        if whichrho == 'zxy':
            penetration_depth = - scale_param * \
                                np.sqrt(zeta.resistivity[per_index, 0, 1] * per)
        elif whichrho == 'zyx':
            penetration_depth = - scale_param * \
                                np.sqrt(zeta.resistivity[per_index, 1, 0] * per)
        elif whichrho == 'det':  # the 2X2 complex Z-matrix's determinant abs value
            # determinant value at the given period index
            det2 = np.abs(zeta.det[0][per_index])
            penetration_depth = -scale_param * np.sqrt(0.2 * per * det2 * per)
        else:
            logger.critical(
                "unsupported method to compute penetration depth: %s",
                whichrho)
            # sys.exit(100)
            raise Exception("unsupported method to compute penetratoin depth: %s" % whichrho)

        pen_depth.append(penetration_depth)

    return stations, periods, pen_depth, latlons


def load_edi_files(edi_path):
    edi_list = []
    if edi_path is not None:
        edi_list = [mt.MT(os.path.join(edi_path, edi)) for edi in os.listdir(edi_path) if edi.endswith("edi")]
    return edi_list


def get_penetration_depth_generic(edi_file_list, period_sec, whichrho='det', ptol=0.1):
    """
    This is a more generic and useful function to compute the penetration depths
    of a list of edi files at given period_sec (seconds).
    No assumption is made about the edi files period list.
    A tolerance of 10% is used to identify the relevant edi files which contain the period of interest.

    :param ptol: freq error/tolerance, need to be consistent with phase_tensor_map.py, default is 0.1
    :param edi_file_list: edi file list of mt object list
    :param period_sec: the float number value of the period in second: 0.1, ...20.0
    :param whichrho:
    :return: tuple of (stations, periods, penetrationdepth, lat-lons-pairs)
    """

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))
    logger.debug("The scaling parameter=%.6f" % scale_param)

    # per_index=0,1,2,....
    periods = []
    pendep = []
    stations = []
    latlons = []

    all_freqs = []  # gather all freqs

    for afile in edi_file_list:
        if isinstance(afile, str) and os.path.isfile(afile):
            mt_obj = mt.MT(afile)
        elif isinstance(afile, mt.MT):
            mt_obj = afile
        else:
            raise Exception("Unsupported list of objects %s" % type(afile))

        all_freqs.extend(list(mt_obj.Z.freq))

        # all stations positions included
        # stations.append((mt_obj.lat, mt_obj.lon))
        stations.append(mt_obj.station)

        p_index = [ff for ff, f2 in enumerate(1.0 / mt_obj.Z.freq)
                   if (f2 > period_sec * (1 - ptol)) and (f2 < period_sec * (1 + ptol))]

        logger.debug("Period index found: %s", p_index)

        if len(p_index) >= 1:  # this edi can be included
            per_index = p_index[0]

            # stations.append(mt_obj.station)
            latlons.append((mt_obj.lat, mt_obj.lon))

            # the attribute Z
            zeta = mt_obj.Z

            if per_index >= len(zeta.freq):
                logger.debug(
                    "Number of frequecies (Max per_index)= %s", len(
                        zeta.freq))
                raise Exception(
                    "Index out_of_range Error: period index must be less than number of periods in zeta.freq")

            per = 1.0 / zeta.freq[per_index]
            periods.append(per)

            if whichrho == 'det':  # the 2X2 complex Z-matrix's determinant abs value
                # determinant value at the given period index
                det2 = np.abs(zeta.det[0][per_index])
                penetration_depth = -scale_param * np.sqrt(0.2 * per * det2 * per)
            elif whichrho == 'zxy':
                penetration_depth = -scale_param * np.sqrt(zeta.resistivity[per_index, 0, 1] * per)
            elif whichrho == 'zyx':
                penetration_depth = -scale_param * np.sqrt(zeta.resistivity[per_index, 1, 0] * per)

            else:
                logger.critical(
                    "un-supported method to compute penetration depth: %s", whichrho)
                sys.exit(100)

            pendep.append(penetration_depth)

        else:
            logger.warn(
                '%s was not used in the 3d profile, because it has no required period.',
                afile)
            pass

    # sort all frequencies so that they are in descending order,
    # use set to remove repeats and make an array.
    all_periods = 1.0 / np.array(sorted(list(set(all_freqs)), reverse=True))
    print("Here is a list of ALL the periods in your edi files:\t ", all_periods)

    return stations, periods, pendep, latlons


class ZComponentError(ParameterError):
    def __init__(self, *args, **kwargs):
        if args is None:
            ParameterError.__init__(self, "please set zcomponent (rho) to either \"zxy\", \"zyx\" or \"det\"", **kwargs)
        else:
            ParameterError.__init__(self, *args, **kwargs)


def check_period_values(period_list, ptol=0.1):
    """
    check if all the values are equal in the input list
    :param period_list: a list of period
    :param ptol=0.1 # 1% percentage tolerance of period values considered as equal
    :return: True/False
    """

    logger.debug("The Periods List to be checked : %s", period_list)

    if not period_list:
        logger.error("The MT periods list is empty - No relevant data found in the EDI files.")
        return False  # what is the sensible default value for this ?

    p0 = period_list[0]  # the first value as a ref

    upper_bound = p0 * (1 + ptol)
    lower_bound = p0 * (1 - ptol)
    if all((per > lower_bound and (per < upper_bound)) for per in period_list[1:]):
        return True
    else:
        return False

        # pcounter = 0
        #
        # for per in period_list:
        #     if (per > p0 * (1 - ptol)) and (per < p0 * (1 + ptol)):  # approximately equal by <5% error
        #         pcounter = pcounter + 1
        #     else:
        #         logger.warn("Periods NOT Equal!!!  %s != %s", p0, per)
        #         return False
        #
        # return True


def get_bounding_box(latlons):
    """ get min max lat lon from the list of lat-lon-pairs points"""
    minlat = minlon = float("inf")
    maxlat = maxlon = float("-inf")

    for (lat, lon) in latlons:
        minlat = min(minlat, lat)
        maxlat = max(maxlat, lat)
        minlon = min(minlon, lon)
        maxlon = max(maxlon, lon)

    logger.debug("Latitude Range: [%.5f, %.5f]", minlat, maxlat)
    logger.debug("Longitude Range: [%.5f, %.5f]", minlon, maxlon)

    # lats = [tup[0] for tup in latlons]
    # lons = [tup[1] for tup in latlons]
    #
    # minlat = min(lats)
    # maxlat = max(lats)
    #
    # print("Latitude Range:", minlat, maxlat)
    #
    # minlon = min(lons)
    # maxlon = max(lons)
    #
    # print("Longitude Range:", minlon, maxlon)

    return (minlon, maxlon), (minlat, maxlat)


def get_index(lat, lon, minlat, minlon, pixelsize, offset=0):
    """
    compute the grid index from the lat lon float value
    :param lat: float lat
    :param lon: float lon
    :param minlat: min lat at low left corner
    :param minlon: min long at left
    :param pixelsize: pixel size in lat long degree
    :param offset: a shift of grid index. should be =0.
    :return: a paire of integer
    """
    index_x = (lon - minlon) / pixelsize
    index_y = (lat - minlat) / pixelsize

    ix = int(round(index_x))
    iy = int(round(index_y))

    # any negative values, out-of-bound?
    logger.debug("Grid index: (%s, %s)", ix, iy)

    return ix + offset, iy + offset