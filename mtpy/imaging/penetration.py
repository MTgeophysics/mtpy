"""
    Description:
        todo: write description

    Usage:
        todo: write usage

    Author: YingzhiGou
    Date: 20/06/2017
"""

import os
from imaging_base import ImagingBase
import numpy as np
import matplotlib.pyplot as plt

import mtpy.core.mt as mt
import mtpy.modeling.occam2d_rewrite as occam2d
from mtpy.utils.mtpylog import MtPyLog

# get a logger object for this module, using the utility class MtPyLog to
# config the logger
logger = MtPyLog().get_mtpy_logger(__name__)
# default contains of rholist
DEFAULT_RHOLIST = set(['zxy', 'zyx', 'det'])


class Depth1D(ImagingBase):
    def set_data(self, data):
        # this plot only use one edi each time
        self._set_edi(data)

    def __init__(self, edis=None, rholist=DEFAULT_RHOLIST):
        super(Depth1D, self).__init__()
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

    def plot(self):
        if self._data == None or not self._data:
            # todo: raise an exception
            raise NotImplemented
        elif self._fig is not None:
            # nothing to plot
            pass

        self._logger.info("Plotting the edi file %s", self._data._get_fn())
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
            det_penetration_depth = scale_param * \
                                    np.sqrt(0.2 * periods * det2 * periods)

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

        title = "Penetration Depth for file %s" % self._data._get_fn()
        plt.title(title)
        plt.xlabel("Log Period (seconds)", fontsize=16)
        plt.ylabel("Penetration Depth (meters)", fontsize=16)
        # set window title
        self._fig.canvas.set_window_title(title)


class Depth2D(ImagingBase):
    def plot(self, **kwargs):
        if self._rho is None:
            raise Exception("please set zcomponent (rho) to either \"zxy\", \"zyx\" or \"det\"")
        elif self._period_indexes is None or not self._period_indexes:
            raise Exception("please provide a period index list")
        else:

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
        if (pr.station_list is not None):
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
        plt.title("MT Penetration Depth Profile by %s" % self._rho)

    def set_data(self, data):
        # this plot require multiple edi files
        self._set_edis(data)

    def __init__(self, data=None, period_index_list=None, rho='det'):
        super(Depth2D, self).__init__()
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
        self.set_rho(rho)
        self.set_period(period)

    def plot(self, **kwargs):
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
        pass



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
        elif whichrho == 'det':   # the 2X2 complex Z-matrix's determinant abs value
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
