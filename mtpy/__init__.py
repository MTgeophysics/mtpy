# -*- coding: utf-8 -*-
"""
==================
MTpy
==================

"""

# define mtpy release version through the variable __version__
# see https://packaging.python.org/guides/single-sourcing-package-version/
__version__ = "1.1.5"

# load mtpy default logging config
from mtpy.utils.mtpy_logger import load_configure, get_mtpy_logger

load_configure()

debug_logger = get_mtpy_logger(__name__, fn="mtpy_debug", level="debug")
debug_logger.debug("Starting MTpy Debug Log File")

error_logger = get_mtpy_logger("error", fn="mtpy_error", level="error")
matplotlib_logger = get_mtpy_logger("matplotlib", fn="matplotlib_warn", level="warning")
