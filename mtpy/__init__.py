
# define mtpy release version through the variable __version__
# see https://packaging.python.org/guides/single-sourcing-package-version/
__version__ = "1.1.5"

import logging

# load mtpy default logging config
from mtpy.utils.mtpylog import MtPyLog

MtPyLog.load_configure()

logging.getLogger('matplotlib').setLevel(logging.WARNING)
