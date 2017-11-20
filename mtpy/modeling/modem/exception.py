"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""

__all__ = ['ModEMError', 'DataError']


class ModEMError(Exception):
    pass


class DataError(ModEMError):
    """Raise for ModEM Data class specific exceptions"""
    pass


class ModelError(ModEMError):
    """ Raise for ModEM Model class specific exceptions"""
    pass


class CovarianceError(ModEMError):
    """ Raise for Covariance class specific exceptions"""
    pass
