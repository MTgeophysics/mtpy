
# define mtpy release version through the variable __version__
# see https://packaging.python.org/guides/single-sourcing-package-version/
__version__ = "1.1.5"

import logging
import os
from pathlib import Path

# load mtpy default logging config
from mtpy.utils.mtpylog import MtPyLog
CONFIG_PATH = Path(os.path.dirname(os.path.abspath(__file__))).joinpath("utils")
CONFIG_FILE = Path.joinpath(CONFIG_PATH, "logging.yml")

MtPyLog.load_configure(CONFIG_FILE)

logging.getLogger('matplotlib').setLevel(logging.WARNING)
