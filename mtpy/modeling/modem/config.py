"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
import os

from mtpy.utils import configfile as mtcfg
from .exception import ModEMError
from .data import Data
from .model import Model
from .convariance import Covariance

__all__ = ['ModEMConfig']


class ModEMConfig(object):
    """
    read and write configuration files for how each inversion is run
    """

    def __init__(self, **kwargs):
        self.cfg_dict = {'ModEM_Inversion_Parameters': {}}

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    def write_config_file(self, save_dir=None,
                          config_fn_basename='ModEM_inv.cfg'):
        """
        write a config file based on provided information
        """

        if save_dir is None:
            save_dir = os.getcwd()

        cfg_fn = os.path.join(save_dir, config_fn_basename)

        if self.cfg_dict is not None:
            mtcfg.write_dict_to_configfile(self.cfg_dict,
                                           cfg_fn)

    def add_dict(self, fn=None, obj=None):
        """
        add dictionary based on file name or object
        """

        if fn is not None:
            if fn.endswith('.rho'):
                m_obj = Model()
                m_obj.read_model_file(fn)
            elif fn.endswith('.dat'):
                m_obj = Data()
                m_obj.read_data_file(fn)
            elif fn.endswith('.cov'):
                m_obj = Covariance()
                m_obj.read_cov_fn(fn)
        elif obj is not None:
            m_obj = obj

        else:
            raise ModEMError('Need to input a file name or object')

        add_dict = m_obj.get_parameters()

        for key in list(add_dict.keys()):
            self.cfg_dict['ModEM_Inversion_Parameters'][key] = add_dict[key]


