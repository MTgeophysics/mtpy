#!/bin/env python
"""
Description:
    Class for loading MODEM data and plotting phase tensors.

References:
 
CreationDate:   11/24/17
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     11/24/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""
import os
import mtpy.modeling.modem as modem
import numpy as np
import logging, traceback

from matplotlib.patches import Polygon as MPLPolygon
from mpl_toolkits.basemap import pyproj

class plot_modem_ptensors:
    def __init__(self, data_fn):
        '''
        Class for loading MODEM data and plotting phase tensors.

        :param data_fn: MODEM data file name
        '''

        self._data_fn = data_fn

        try:
            self._modem_obj = modem.Data()
            self._modem_obj.read_data_file(self._data_fn)
        except Exception as err:
            print 'Failed to read %s' % (self._data_fn)
            logging.error(traceback.format_exc())
            exit(-1)

        self._plot_period = self._modem_obj.period_list.copy()
        self._mt_obj_list = [self._modem_obj.mt_dict[key]
                             for key in self._modem_obj.mt_dict.keys()]

        # Read data
        self._pt_dict = {}
        self._ptol = 0.05
        for plot_per in self._plot_period:
            self._pt_dict[plot_per] = []
            for mt_obj in self._mt_obj_list:
                p_index = [ff for ff, f2 in enumerate(1. / mt_obj.Z.freq)
                           if (f2 > plot_per * (1 - self._ptol)) and
                           (f2 < plot_per * (1 + self._ptol))][0]

                pt_tuple = (mt_obj.station, mt_obj.lon, mt_obj.lat,
                            mt_obj.pt.phimin[p_index],
                            mt_obj.pt.phimax[p_index],
                            mt_obj.pt.azimuth[p_index],
                            mt_obj.pt.beta[p_index],
                            2 * mt_obj.pt.beta[p_index],
                            mt_obj.pt.ellipticity[p_index])

                self._pt_dict[plot_per].append(pt_tuple)
            # end for
            self._pt_dict[plot_per] = np.array(self._pt_dict[plot_per],
                                               dtype=[('station', '|S15'),
                                                      ('lon', np.float),
                                                      ('lat', np.float),
                                                      ('phimin', np.float),
                                                      ('phimax', np.float),
                                                      ('azimuth', np.float),
                                                      ('skew', np.float),
                                                      ('n_skew', np.float),
                                                      ('ellipticity', np.float)])
        # end for
    # end func

    def ellipse(self, x0, y0, a, b, azimuth, n, ax, m, **kwargs):
        '''
        This function has been adapted from:
        https://stackoverflow.com/questions/8161144/drawing-ellipses-on-matplotlib-basemap-projections

        Draws a polygon centered at ``x0, y0``. The polygon approximates an
        ellipse on the surface of the Earth with semi-major-axis ``a`` and
        semi-minor axis ``b`` degrees longitude and latitude, made up of
        ``n`` vertices.

        For a description of the properties of ellipsis, please refer to [1].

        The polygon is based upon code written do plot Tissot's indicatrix
        found on the matplotlib mailing list.

        :param x0: latitude
        :param y0: longitude
        :param a: semi-major axis of ellipse
        :param b: semi-minor axis of ellipse
        :param azimuth: azimuth of rotated ellipse
        :param n: number of vertices in ellipse
        :param ax: plot axis
        :param m: instance of Basemap
        :param kwargs: list of matplotlib relevant keywords, e.g. alpha, zorder, color, etc. 
        :return: returns a matplotlib Polygon instance
        '''

        g = pyproj.Geod(a=m.rmajor, b=m.rminor)
        # Gets forward and back azimuths, plus distances between initial
        # points (x0, y0)
        azf, azb, dist = g.inv([x0, x0], [y0, y0], [x0 + a, x0], [y0, y0 + b])
        tsid = dist[0] * dist[1]  # a * b

        # Initializes list of segments, calculates \del azimuth, and goes on
        # for every vertex
        seg = [m(x0 + a, y0)]
        AZ = np.linspace(azf[0], 360. + azf[0], n)
        for i, az in enumerate(AZ):
            # Skips segments along equator (Geod can't handle equatorial arcs).
            if np.allclose(0., y0) and (np.allclose(90., az) or
                                               np.allclose(270., az)):
                continue

            # In polar coordinates, with the origin at the center of the
            # ellipse and with the angular coordinate ``az`` measured from the
            # major axis, the ellipse's equation  is [1]:
            #
            #                           a * b
            # r(az) = ------------------------------------------
            #         ((b * cos(az))**2 + (a * sin(az))**2)**0.5
            #
            # Azymuth angle in radial coordinates and corrected for reference
            # angle.
            rot = np.deg2rad(azimuth)
            azr = 2. * np.pi / 360. * (az + 90.)
            # A = dist[0] * np.sin(azr)
            # B = dist[1] * np.cos(azr)

            A = dist[0] * np.sin(azr) * np.cos(rot) - dist[1] * np.cos(azr) * np.sin(rot)
            B = dist[0] * np.sin(azr) * np.sin(rot) + dist[1] * np.cos(azr) * np.cos(rot)

            r = tsid / (B ** 2. + A ** 2.) ** 0.5
            lon, lat, azb = g.fwd(x0, y0, az - azimuth, r)
            x, y = m(lon, lat)

            # Add segment if it is in the map projection region.
            if x < 1e20 and y < 1e20:
                seg.append((x, y))

        poly = MPLPolygon(seg, **kwargs)
        ax.add_patch(poly)

        # Set axes limits to fit map region.
        m.set_axes_limits(ax=ax)

        return poly
    # end func

    def get_period_attributes(self, periodIdx, key):
        '''
        Returns, for a given period, a list of attribute values for key
        (e.g. skew, phimax, etc.).

        :param periodIdx: index of period; print out _plot_period for periods available
        :param key: attribute key
        :return: numpy array of attribute values
        '''
        assert (periodIdx > 0 and periodIdx < len(self._plot_period)), \
            'Error: Index for plot-period out of bounds.'

        pk = self._pt_dict.keys()[periodIdx]

        try:
            vals = self._pt_dict[pk][key]
            return vals
        except Exception as err:
            print 'Key error: key %s not found' % (key)
            logging.error(traceback.format_exc())
            exit(-1)
    # end func

    def plot(self, ax, m, periodIdx, ellipse_size_factor=0.05,
             nverts=100, cvals=None, **kwargs):

        '''
        Plots phase tensors for a given period index.

        :param ax: plot axis
        :param m: basemap instance
        :param periodIdx: period index
        :param ellipse_size_factor: factor to control ellipse size
        :param nverts: number of vertices in each ellipse
        :param cvals: list of colour values for colouring each ellipse
        :param kwargs: list of relevant matplotlib arguments (e.g. zorder, alpha, etc.)
        '''

        assert (periodIdx > 0 and periodIdx < len(self._plot_period)), \
            'Error: Index for plot-period out of bounds.'

        k = self._pt_dict.keys()[periodIdx]
        for i in range(len(self._pt_dict[k])):
            lon = self._pt_dict[k]['lon'][i]
            lat = self._pt_dict[k]['lat'][i]
            phimax = self._pt_dict[k]['phimax'][i] / self.pt_dict[k]['phimax'].max()
            phimin = self._pt_dict[k]['phimin'][i] / self.pt_dict[k]['phimax'].max()
            az = self._pt_dict[k]['azimuth'][i]
            nskew = self._pt_dict[k]['n_skew'][i]

            # print az
            if (phimax > 0 and phimin > 0):
                c = None
                if (cvals is not None): c = cvals[i]
                if (c is not None): kwargs['facecolor'] = c
                self.ellipse(lon, lat, phimax * ellipse_size_factor,
                             phimin * ellipse_size_factor, az,
                             nverts, ax, m, **kwargs)
            # end if
        # end for
    # end func
# end class

def main():
    """
    define main function
    :return:
    """

    imaging2 = os.path.dirname(__file__)
    mtpy = os.path.dirname(imaging2)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    ModEM_files = os.path.join(data, 'ModEM_files')
    data_fn = os.path.join(ModEM_files, 'ModEM_Data_im2.dat')

    pmp = plot_modem_ptensors(data_fn=data_fn)

    return
# end


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    main()