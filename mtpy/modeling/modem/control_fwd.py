"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""

import os

from mtpy.utils import exceptions as mtex

__all__ = ['ControlFwd']


class ControlFwd(object):
    """
    read and write control file for

    This file controls how the inversion starts and how it is run

    """

    def __init__(self, **kwargs):

        self.num_qmr_iter = kwargs.pop('num_qmr_iter', 40)
        self.max_num_div_calls = kwargs.pop('max_num_div_calls', 20)
        self.max_num_div_iters = kwargs.pop('max_num_div_iters', 100)
        self.misfit_tol_fwd = kwargs.pop('misfit_tol_fwd', 1.0e-7)
        self.misfit_tol_adj = kwargs.pop('misfit_tol_adj', 1.0e-7)
        self.misfit_tol_div = kwargs.pop('misfit_tol_div', 1.0e-5)

        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.fn_basename = kwargs.pop('fn_basename', 'control.fwd')
        self.control_fn = kwargs.pop('control_fn', os.path.join(self.save_path,
                                                                self.fn_basename))

        self._control_keys = ['Number of QMR iters per divergence correction',
                              'Maximum number of divergence correction calls',
                              'Maximum number of divergence correction iters',
                              'Misfit tolerance for EM forward solver',
                              'Misfit tolerance for EM adjoint solver',
                              'Misfit tolerance for divergence correction']

        self._control_dict = dict([(key, value)
                                   for key, value in zip(self._control_keys,
                                                         [self.num_qmr_iter,
                                                          self.max_num_div_calls,
                                                          self.max_num_div_iters,
                                                          self.misfit_tol_fwd,
                                                          self.misfit_tol_adj,
                                                          self.misfit_tol_div])])
        self._string_fmt_dict = dict([(key, value)
                                      for key, value in zip(self._control_keys,
                                                            ['<.0f', '<.0f', '<.0f', '<.1e', '<.1e',
                                                             '<.1e'])])

    def write_control_file(self, control_fn=None, save_path=None,
                           fn_basename=None):
        """
        write control file

        Arguments:
        ------------
            **control_fn** : string
                             full path to save control file to
                             *default* is save_path/fn_basename

            **save_path** : string
                            directory path to save control file to
                            *default* is cwd

            **fn_basename** : string
                              basename of control file
                              *default* is control.inv

        """

        if control_fn is not None:
            self.save_path = os.path.dirname(control_fn)
            self.fn_basename = os.path.basename(control_fn)

        if save_path is not None:
            self.save_path = save_path

        if fn_basename is not None:
            self.fn_basename = fn_basename

        self.control_fn = os.path.join(self.save_path, self.fn_basename)

        self._control_dict = dict([(key, value)
                                   for key, value in zip(self._control_keys,
                                                         [self.num_qmr_iter,
                                                          self.max_num_div_calls,
                                                          self.max_num_div_iters,
                                                          self.misfit_tol_fwd,
                                                          self.misfit_tol_adj,
                                                          self.misfit_tol_div])])

        clines = []
        for key in self._control_keys:
            value = self._control_dict[key]
            str_fmt = self._string_fmt_dict[key]
            clines.append('{0:<47}: {1:{2}}\n'.format(key, value, str_fmt))

        cfid = file(self.control_fn, 'w')
        cfid.writelines(clines)
        cfid.close()

        print('Wrote ModEM control file to {0}'.format(self.control_fn))

    def read_control_file(self, control_fn=None):
        """
        read in a control file
        """

        if control_fn is not None:
            self.control_fn = control_fn

        if self.control_fn is None:
            raise mtex.MTpyError_file_handling('control_fn is None, input '
                                               'control file')

        if os.path.isfile(self.control_fn) is False:
            raise mtex.MTpyError_file_handling('Could not find {0}'.format(
                self.control_fn))

        self.save_path = os.path.dirname(self.control_fn)
        self.fn_basename = os.path.basename(self.control_fn)

        cfid = file(self.control_fn, 'r')
        clines = cfid.readlines()
        cfid.close()
        for cline in clines:
            clist = cline.strip().split(':')
            if len(clist) == 2:

                try:
                    self._control_dict[clist[0].strip()] = float(clist[1])
                except ValueError:
                    self._control_dict[clist[0].strip()] = clist[1]

        # set attributes
        attr_list = ['num_qmr_iter', 'max_num_div_calls', 'max_num_div_iters',
                     'misfit_tol_fwd', 'misfit_tol_adj', 'misfit_tol_div']
        for key, kattr in zip(self._control_keys, attr_list):
            setattr(self, kattr, self._control_dict[key])


