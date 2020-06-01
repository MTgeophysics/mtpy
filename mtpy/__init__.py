
# define mtpy release version below
# see https://packaging.python.org/guides/single-sourcing-package-version/
__version__ = "1.1.3"

import logging

# load mtpy default logging config
from mtpy.utils.mtpylog import MtPyLog

MtPyLog.load_configure()

logging.getLogger('matplotlib').setLevel(logging.WARNING)
