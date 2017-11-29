#!/bin/env python
"""
Description:
    Class for extracting field values at arbitrary locations from a 3D MODEM model.

References:
 
CreationDate:   11/27/17
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     11/27/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import mtpy.modeling.modem as modem
import numpy as np
import logging, traceback
from scipy.spatial import cKDTree

class MODEM_slice:
    def __init__(self, model_fn, data_fn):
        '''
        Class for extracting field values at arbitrary locations from a 3D MODEM model.

        :param model_fn: name of model file
        :param data_fn: name of data file
        '''

        # Read data and model files
        self._model_fn = model_fn
        self._data_fn = data_fn
        try:
            self._d = modem.Data()
            self._d.read_data_file(data_fn=self._data_fn)
        except Exception as err:
            print 'Failed to read %s' % (self._data_fn)
            logging.error(traceback.format_exc())
            exit(-1)

        try:
            self._m = modem.Model(model_fn=self._model_fn,
                                  station_object=self._d.station_locations)
            self._m.read_model_file()
        except Exception as err:
            print 'Failed to read %s' % (self._model_fn)
            logging.error(traceback.format_exc())
            exit(-1)

        # Re-orient model coordinates based on the centre of data locations
        self._mx = self._m.grid_east + self._d.center_point['east']
        self._my = self._m.grid_north + self._d.center_point['north']
        self._mz = self._m.grid_z

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

    def get_slice(self, xyz_list, nn=1, p=4, extrapolate=True):
        '''
        Function to retrieve interpolated field values at arbitrary locations

        :param xyz_list: numpy array of shape (np,3), where np in the number of points
        :param nn: Number of neighbours to use for interpolation.
                   Nearest neighbour interpolation is returned when nn=1 (default).
                   When nn>1, inverse distance weighted interpolation is returned. See
                   link below for more details:
                   https://en.wikipedia.org/wiki/Inverse_distance_weighting
        :param p: Power parameter, which determines the relative influence of near and far
                  neighbours during interpolation. For p<=3, causes interpolated values to
                  be dominated by points far away. Larger values of p assign greater influence
                  to values near the interpolated point.
        :param extrapolate: Extrapolates values (default), which can be particularly useful
                            for extracting values at nodes, since the field values are given
                            for cell-centres.
        :return: numpy array of interpolated values of shape (np)
        '''

        # query Kd-tree instance to retrieve distances and
        # indices of k nearest neighbours
        d, l = self._tree.query(xyz_list, k=nn)

        img = None
        if (nn == 1):
            # extract nearest neighbour values
            img = self._m.res_model.flatten()[l]
        else:
            vals = self._m.res_model.flatten()
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
# end class

def main():
    """
    define main function
    :return:
    """

    utils = os.path.dirname(__file__)
    mtpy = os.path.dirname(utils)
    base = os.path.dirname(mtpy)
    examples = os.path.join(base, 'examples')
    data = os.path.join(examples, 'data')
    ModEM_files = os.path.join(data, 'ModEM_files')
    d_fn = os.path.join(ModEM_files, 'ModEM_Data_im2.dat')
    m_fn = os.path.join(ModEM_files, 'Modular_MPI_NLCG_056_im2.rho')

    ms = MODEM_slice(model_fn=m_fn, data_fn=d_fn)
    return
# end


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    main()