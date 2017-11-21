#!/usr/bin/env python
"""
Description:
   Define the Covariance class.
   This module is refactored from modem.py which is too big to manage and edit

Author: fei.zhang@ga.gov.au

Date: 2017-06-30
"""

__author__ = 'fei.zhang@ga.gov.au'

# from __future__ import print_function
import os

import numpy as np

from mtpy.modeling.modem_model import Model
from mtpy.utils.mtpylog import MtPyLog

try:
    from evtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')

    print ('Note: if you are using Windows you should build evtk first with'
           'either MinGW or cygwin using the command: \n'
           '    python setup.py build -compiler=mingw32  or \n'
           '    python setup.py build -compiler=cygwin')

logger = MtPyLog.get_mtpy_logger(__name__)


class MTException(Exception):
    """Raise for MTPy specific exception"""
    pass


# ==============================================================================
# covariance
# ==============================================================================
class Covariance(object):
    """ Read and write covariance files
    Arguments
    -----------

    ====================== ====================================================
    Attributes/Key Words   Description
    ====================== ====================================================
    grid_dimensions        Grid dimensions excluding air layers (Nx, Ny, NzEarth)
    smoothing_east         Smoothing in the X direction (NzEarth real values)
    smoothing_north        Smoothing in the Y direction (NzEarth real values)
    smoothing_z            Vertical smoothing (1 real value)
    smoothing_num          Number of times the smoothing should be applied (1 integer >= 0)

    exception_list         Exceptions in the for e.g. 2 3 0. (to turn off smoothing between 3 & 4)
    mask_arr

    save_path              path to save data file to
    cov_fn_basename
    cov_fn
    """

    def __init__(self, grid_dimensions=None, **kwargs):

        self.grid_dimensions = grid_dimensions
        self.smoothing_east = kwargs.pop('smoothing_east', 0.3)
        self.smoothing_north = kwargs.pop('smoothing_north', 0.3)
        self.smoothing_z = kwargs.pop('smoothing_z', 0.3)
        self.smoothing_num = kwargs.pop('smoothing_num', 1)

        self.exception_list = kwargs.pop('exception_list', [])
        self.mask_arr = kwargs.pop('mask_arr', None)

        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.cov_fn_basename = kwargs.pop('cov_fn_basename', 'covariance.cov')

        self.cov_fn = kwargs.pop('cov_fn', None)

        self._header_str = '\n'.join(['+{0}+'.format('-' * 77),
                                      '| This file defines model covariance for a recursive autoregression scheme.   |',
                                      '| The model space may be divided into distinct areas using integer masks.     |',
                                      '| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |',
                                      '| air, ocean and the rest of the model is turned off automatically. You can   |',
                                      '| also define exceptions to override smoothing between any two model areas.   |',
                                      '| To turn off smoothing set it to zero.  This header is 16 lines long.        |',
                                      '| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |',
                                      '| 2. Smoothing in the X direction (NzEarth real values)                       |',
                                      '| 3. Smoothing in the Y direction (NzEarth real values)                       |',
                                      '| 4. Vertical smoothing (1 real value)                                        |',
                                      '| 5. Number of times the smoothing should be applied (1 integer >= 0)         |',
                                      '| 6. Number of exceptions (1 integer >= 0)                                    |',
                                      '| 7. Exceptions in the for e.g. 2 3 0. (to turn off smoothing between 3 & 4)  |',
                                      '| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|',
                                      '+{0}+'.format('-' * 77)])


    def write_covariance_file(self, cov_fn=None, save_path=None,
                              cov_fn_basename=None, model_fn=None,
                              sea_water=0.3, air=1e17):
        """
        write a covariance file
        """

        if model_fn is not None:
            mod_obj = Model()
            mod_obj.read_model_file(model_fn)
            print 'Done Reading {0}'.format(model_fn)

            self.grid_dimensions = mod_obj.res_model.shape
            if self.mask_arr is None:
                self.mask_arr = np.ones_like(mod_obj.res_model)
                self.mask_arr[np.where(mod_obj.res_model > air * .9)] = 0
                self.mask_arr[np.where((mod_obj.res_model < sea_water * 1.1) &
                                       (mod_obj.res_model > sea_water * .9))] = 9
                # flip mask arr as it needs to be in opposite order
                self.mask_arr = self.mask_arr[::-1]
        else:
            if self.mask_arr is None:
                self.mask_arr = np.ones((self.grid_dimensions[0],
                                         self.grid_dimensions[1],
                                         self.grid_dimensions[2]))

        if self.grid_dimensions is None:
            raise MTException('Grid dimensions are None, input as (Nx, Ny, Nz)')

        if cov_fn is not None:
            self.cov_fn = cov_fn
        else:
            if save_path is not None:
                self.save_path = save_path
            if cov_fn_basename is not None:
                self.cov_fn_basename = cov_fn_basename
            self.cov_fn = os.path.join(self.save_path, self.cov_fn_basename)

        clines = [self._header_str]
        clines.append('\n\n')

        # --> grid dimensions
        clines.append(' {0:<10}{1:<10}{2:<10}\n'.format(self.grid_dimensions[0],
                                                        self.grid_dimensions[
                                                            1],
                                                        self.grid_dimensions[2]))
        clines.append('\n')

        # --> smoothing in north direction
        n_smooth_line = ''
        for zz in range(self.grid_dimensions[2]):
            n_smooth_line += ' {0:<5.1f}'.format(self.smoothing_north)
        clines.append(n_smooth_line + '\n')

        # --> smoothing in east direction
        e_smooth_line = ''
        for zz in range(self.grid_dimensions[2]):
            e_smooth_line += ' {0:<5.1f}'.format(self.smoothing_east)
        clines.append(e_smooth_line + '\n')

        # --> smoothing in vertical direction
        clines.append(' {0:<5.1f}\n'.format(self.smoothing_z))
        clines.append('\n')

        # --> number of times to apply smoothing
        clines.append(' {0:<2.0f}\n'.format(self.smoothing_num))
        clines.append('\n')

        # --> exceptions
        clines.append(' {0:<.0f}\n'.format(len(self.exception_list)))
        for exc in self.exception_list:
            clines.append('{0:<5.0f}{1:<5.0f}{2:<5.0f}\n'.format(exc[0],
                                                                 exc[1],
                                                                 exc[2]))
        clines.append('\n')
        clines.append('\n')
        # --> mask array
        # self.mask_arr was constructed in the Model.add_topography()
        # and passed to there through constructor param mask_arr=model.covariance_mask
        for zz in range(self.mask_arr.shape[2]):
            clines.append(' {0:<8.0f}{0:<8.0f}\n'.format(zz + 1))

            for nn in range(self.mask_arr.shape[0]):
                cline = ''
                for ee in range(self.mask_arr.shape[1]):
                    cline += '{0:^3.0f}'.format(self.mask_arr[nn, ee, zz])
                clines.append(cline + '\n')

        cfid = file(self.cov_fn, 'w')
        cfid.writelines(clines)
        cfid.close()

        print 'Wrote covariance file to {0}'.format(self.cov_fn)

        return self.cov_fn

# ======================================
# example usage
# ======================================
if __name__ == "__name__":

    # make covariance file

    model=None  # define modem_model

    cov = Covariance(mask_arr=model.covariance_mask,
                     save_path="/outputdir",
                     smoothing_east=0.3,
                     smoothing_north=0.4,
                     smoothing_z=0.5)

    cov.write_covariance_file(model_fn=model.model_fn)
