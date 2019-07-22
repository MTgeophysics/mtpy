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

__all__ = ['ControlInv']


class ControlInv(object):
    """
    read and write control file for how the inversion starts and how it is run

    """

    def __init__(self, **kwargs):

        self.output_fn = kwargs.pop('output_fn', 'MODULAR_NLCG')
        self.lambda_initial = kwargs.pop('lambda_initial', 10)
        self.lambda_step = kwargs.pop('lambda_step', 10)
        self.model_search_step = kwargs.pop('model_search_step', 1)
        self.rms_reset_search = kwargs.pop('rms_reset_search', 2.0e-3)
        self.rms_target = kwargs.pop('rms_target', 1.05)
        self.lambda_exit = kwargs.pop('lambda_exit', 1.0e-4)
        self.max_iterations = kwargs.pop('max_iterations', 100)
        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.fn_basename = kwargs.pop('fn_basename', 'control.inv')
        self.control_fn = kwargs.pop('control_fn', os.path.join(self.save_path,
                                                                self.fn_basename))

        self._control_keys = ['Model and data output file name',
                              'Initial damping factor lambda',
                              'To update lambda divide by',
                              'Initial search step in model units',
                              'Restart when rms diff is less than',
                              'Exit search when rms is less than',
                              'Exit when lambda is less than',
                              'Maximum number of iterations']

        self._control_dict = dict([(key, value)
                                   for key, value in zip(self._control_keys,
                                                         [self.output_fn, self.lambda_initial,
                                                          self.lambda_step, self.model_search_step,
                                                          self.rms_reset_search, self.rms_target,
                                                          self.lambda_exit, self.max_iterations])])
        self._string_fmt_dict = dict([(key, value)
                                      for key, value in zip(self._control_keys,
                                                            ['<', '<.1f', '<.1f', '<.1f', '<.1e',
                                                             '<.2f', '<.1e', '<.0f'])])

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
                                                         [self.output_fn, self.lambda_initial,
                                                          self.lambda_step, self.model_search_step,
                                                          self.rms_reset_search, self.rms_target,
                                                          self.lambda_exit, self.max_iterations])])

        clines = []
        for key in self._control_keys:
            value = self._control_dict[key]
            str_fmt = self._string_fmt_dict[key]
            clines.append('{0:<35}: {1:{2}}\n'.format(key, value, str_fmt))

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
        attr_list = ['output_fn', 'lambda_initial', 'lambda_step',
                     'model_search_step', 'rms_reset_search', 'rms_target',
                     'lambda_exit', 'max_iterations']
        for key, kattr in zip(self._control_keys, attr_list):
            setattr(self, kattr, self._control_dict[key])


