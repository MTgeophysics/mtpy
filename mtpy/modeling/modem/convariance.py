"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
import os

import numpy as np

from mtpy.utils.mtpylog import MtPyLog
from .exception import CovarianceError
from .model import Model

try:
    from evtk.hl import gridToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')

__all__ = ['Covariance']


class Covariance(object):
    """
    read and write covariance files

    """

    def __init__(self, grid_dimensions=None, **kwargs):
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)

        self.grid_dimensions = grid_dimensions
        self.smoothing_east = 0.3
        self.smoothing_north = 0.3
        self.smoothing_z = 0.3
        self.smoothing_num = 1

        self.exception_list = []
        self.mask_arr = None

        self.save_path = os.getcwd()
        self.cov_fn_basename = 'covariance.cov'

        self.cov_fn = None

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

        for key in list(kwargs.keys()):
            if hasattr(self, key):
                setattr(self, key, kwargs[key])
            else:
                self._logger.warn("Argument {}={} is not supportted thus not been set.".format(key, kwargs[key]))



    def write_covariance_file(self, cov_fn=None, save_path=None,
                              cov_fn_basename=None, model_fn=None,
                              sea_water=0.3, air=1e12):  #
        """
        write a covariance file
        """

        if model_fn is not None:
            mod_obj = Model()
            mod_obj.read_model_file(model_fn)

            # update save_path from model path if not provided separately
            if save_path is None:
                save_path = os.path.dirname(model_fn)

            print('Reading {0}'.format(model_fn))
            self.grid_dimensions = mod_obj.res_model.shape
            if self.mask_arr is None:
                self.mask_arr = np.ones_like(mod_obj.res_model)
            self.mask_arr[np.where(mod_obj.res_model >= air * .9)] = 0
            self.mask_arr[np.where((mod_obj.res_model <= sea_water * 1.1) &
                                   (mod_obj.res_model >= sea_water * .9))] = 9

        if self.grid_dimensions is None:
            raise CovarianceError('Grid dimensions are None, input as (Nx, Ny, Nz)')

        if cov_fn is not None:
            self.cov_fn = cov_fn
        else:
            if save_path is not None:
                self.save_path = save_path
            if cov_fn_basename is not None:
                self.cov_fn_basename = cov_fn_basename
            self.cov_fn = os.path.join(self.save_path, self.cov_fn_basename)

        clines = [self._header_str, '\n\n', ' {0:<10}{1:<10}{2:<10}\n'.format(self.grid_dimensions[0],
                                                                              self.grid_dimensions[1],
                                                                              self.grid_dimensions[2]), '\n']

        # --> grid dimensions


        # --> smoothing in north direction
        n_smooth_line = ''
        for zz in range(self.grid_dimensions[2]):
            if not np.iterable(self.smoothing_north):
                n_smooth_line += ' {0:<5.2f}'.format(self.smoothing_north)
            else:
                n_smooth_line += ' {0:<5.2f}'.format(self.smoothing_north[zz])
        clines.append(n_smooth_line + '\n')

        # --> smoothing in east direction
        e_smooth_line = ''
        for zz in range(self.grid_dimensions[2]):
            if not np.iterable(self.smoothing_east):
                e_smooth_line += ' {0:<5.2f}'.format(self.smoothing_east)
            else:
                e_smooth_line += ' {0:<5.2f}'.format(self.smoothing_east[zz])
        clines.append(e_smooth_line + '\n')

        # --> smoothing in vertical direction
        clines.append(' {0:<5.2f}\n'.format(self.smoothing_z))
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
        if self.mask_arr is None:
            self.mask_arr = np.ones((self.grid_dimensions[0],
                                     self.grid_dimensions[1],
                                     self.grid_dimensions[2]))

        # need to flip north and south.
        write_mask_arr = self.mask_arr[::-1, :, :].copy()
        for zz in range(self.mask_arr.shape[2]):
            clines.append(' {0:<8.0f}{0:<8.0f}\n'.format(zz + 1))
            for nn in range(self.mask_arr.shape[0]):
                cline = ''
                for ee in range(self.mask_arr.shape[1]):
                    cline += '{0:^3.0f}'.format(write_mask_arr[nn, ee, zz])
                clines.append(cline + '\n')

        with open(self.cov_fn, 'w') as cfid:
            cfid.writelines(clines)

        # not needed cfid.close()

        self._logger.info('Wrote covariance file to {0}'.format(self.cov_fn))

    def read_cov_file(self, cov_fn):
        """
        read a covariance file
        """
        if not os.path.isfile(cov_fn):
            raise CovarianceError('{0} not found, check path'.format(cov_fn))

        self.cov_fn = cov_fn
        self.save_path = os.path.dirname(self.cov_fn)
        self.cov_fn_basename = os.path.basename(self.cov_fn)

        with open(cov_fn, 'r') as fid:
            lines = fid.readlines()

        num_find = False
        east_find = False
        north_find = False
        count = 0

        for line in lines:
            if line.find('+') >= 0 or line.find('|') >= 0:
                continue
            else:
                line_list = line.strip().split()
                if len(line_list) == 0:
                    continue
                elif len(line_list) == 1 and not num_find and line_list[0].find('.') == -1:
                    self.smoothing_num = int(line_list[0])
                    num_find = True
                elif len(line_list) == 1 and num_find and line_list[0].find('.') == -1:
                    self.exceptions_num = int(line_list[0])
                elif len(line_list) == 1 and line_list[0].find('.') >= 0:
                    self.smoothing_z = float(line_list[0])
                elif len(line_list) == 3:
                    nx, ny, nz = [int(ii) for ii in line_list]
                    self.grid_dimensions = (nx, ny, nz)
                    self.mask_arr = np.ones((nx, ny, nz), dtype=np.int)
                    self.smoothing_east = np.zeros(ny)
                    self.smoothing_north = np.zeros(nx)
                elif len(line_list) == 2:
                    # starts at 1 but python starts at 0
                    index_00, index_01 = [int(ii) - 1 for ii in line_list]
                    count = 0
                elif line_list[0].find('.') >= 0 and north_find == False:
                    self.smoothing_north = np.array(line_list, dtype=np.float)
                    north_find = True
                elif line_list[0].find('.') >= 0 and north_find == True:
                    self.smoothing_east = np.array(line_list, dtype=np.float)
                    east_find = True
                elif north_find and east_find:
                    line_list = np.array(line_list, dtype=np.int)
                    line_list = line_list.reshape((ny, 1))

                    self.mask_arr[count, :, index_00:index_01 + 1] = line_list
                    count += 1

    def get_parameters(self):

        parameter_list = ['smoothing_north',
                          'smoothing_east',
                          'smoothing_z',
                          'smoothing_num']

        parameter_dict = {}
        for parameter in parameter_list:
            key = 'covariance.{0}'.format(parameter)
            parameter_dict[key] = getattr(self, parameter)

        return parameter_dict

    def write_cov_vtk_file(self, cov_vtk_fn, model_fn=None, grid_east=None,
                           grid_north=None, grid_z=None):
        """
        write a vtk file of the covariance to match things up
        """

        if model_fn is not None:
            m_obj = Model()
            m_obj.read_model_file(model_fn)
            grid_east = m_obj.grid_east
            grid_north = m_obj.grid_north
            grid_z = m_obj.grid_z

        if grid_east is not None:
            grid_east = grid_east
        if grid_north is not None:
            grid_north = grid_north
        if grid_z is not None:
            grid_z = grid_z

        # use cellData, this makes the grid properly as grid is n+1
        gridToVTK(cov_vtk_fn,
                  grid_north / 1000.,
                  grid_east / 1000.,
                  grid_z / 1000.,
                  cellData={'covariance_mask': self.mask_arr})

        self._logger.info('Wrote covariance file to {0}\n'.format(cov_vtk_fn))
