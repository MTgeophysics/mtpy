"""
    Description:
        todo: write description

    Usage:
        todo: write usage

    Author: YingzhiGou
    Date: 20/06/2017
"""

import os
import sys

import mtpy
from imaging_base import ImagingBase, ParameterError, ImagingError
import numpy as np
import matplotlib.pyplot as plt

import mtpy.modeling.occam2d_rewrite as occam2d
from mtpy.core import mt as mt
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
logger = MtPyLog().get_mtpy_logger(__name__)
# default contains of rholist
DEFAULT_RHOLIST = {'zxy', 'zyx', 'det'}


class Depth1D(ImagingBase):
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
        elif self._rholist is None or self._rholist:
            raise ZComponentError
        elif self._fig is not None:
            # nothing to plot
            return

        self._logger.info("Plotting the edi file %s", self._data.fn)
        self._fig = plt.figure()
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

        title = "Penetration Depth for file %s" % self._data.fn()
        plt.title(title)
        plt.xlabel("Log Period (seconds)", fontsize=16)
        plt.ylabel("Penetration Depth (meters)", fontsize=16)
        # set window title
        self._fig.canvas.set_window_title(title)


class Depth2D(ImagingBase):
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

        self._fig = plt.figure()
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
    def __init__(self, data=None, period=None, rho='det'):
        super(Depth3D, self).__init__()
        self._rho = None
        self._period = None
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
        if period_by_index:  # self._period is considered as an index
            if not isinstance(self._period, int):
                self._logger.warning("period value is not integer but used as an index.")
            (stations, periods, pendep, latlons) = get_penetration_depth(self._data,
                                                                         int(self._period),
                                                                         whichrho=self._rho)
        else:
            (stations, periods, pendep, latlons) = get_penetration_depth_generic(self._data,
                                                                                 self._period,
                                                                                 whichrho=self._rho)

        # create figure
        self._fig = plt.figure()

        if check_period_values(periods) is False:
            logger.error("The period values are NOT equal - Please check!!! %s", periods)
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
                period_fmt = str(mtpy.utils.calculator.roundsf(period0, 4))
            else:
                period_fmt = "%.2f" % period0
            bbox = get_bounding_box(latlons)

            logger.debug("Bounding Box %s", bbox)

            xgrids = bbox[0][1] - bbox[0][0]
            ygrids = bbox[1][1] - bbox[1][0]

            logger.debug("xy grids: %s %s", xgrids, ygrids)

            minlat = bbox[1][0]
            minlon = bbox[0][0]

            # Pixel size in Degree:  0.001=100meters, 0.01=1KM 1deg=100KM
            pixelsize = 0.002  # Degree 0.002=200meters, 0.01=1KM 1deg=100KM

            nx = int(np.ceil(xgrids/pixelsize))
            ny = int(np.ceil(ygrids/pixelsize))

            logger.debug("number of grids xy: %s %s", nx, ny)

            # make the image slightly bigger than the (nx, ny) to contain all points
            # avoid index out of bound
            pad = 1  # pad = 1 affect the top and right of the plot. it is linked to get_index offset?
            # todo change this part to use xy bound offset? (0.5 gride on each side?)
            nx2 = nx + pad
            ny2 = ny + pad

            # fast initialization
            zdep = np.empty((ny2, nx2))
            zdep.fill(np.nan)  # initialize all pixel value as np.nan

            logger.debug("zdep shape %s", zdep.shape)

            for iter, pair in enumerate(latlons):
                # logger.debug(iter, pair)
                pass

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


def get_penetration_depth_generic(edi_file_list, period_sec, whichrho='det'):
    """
    This is a more generic and useful function to compute the penetration depths
    of a list of edi files at given period_sec (seconds).
    No assumption is made about the edi files period list.
    A tolerance of 10% is used to identify the relevant edi files which contain the period of interest.

    :param edi_file_list: edi file list of mt object list
    :param period_sec: the float number value of the period in second: 0.1, ...20.0
    :param whichrho:
    :return: tuple of (stations, periods, penetrationdepth, lat-lons-pairs)
    """
    ptol = 0.1  # 10% freq error, need to be consistent with phase_tensor_map.py

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
            ParameterError.__init__("please set zcomponent (rho) to either \"zxy\", \"zyx\" or \"det\"", **kwargs)
        else:
            ParameterError.__init__(*args, **kwargs)


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

    upper_bound = p0 * (1+ptol)
    lower_bound = p0 * (1-ptol)
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

    logger.debug("Latitude Range:", minlat, maxlat)
    logger.debug("Longitude Range:", minlon, maxlon)

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
