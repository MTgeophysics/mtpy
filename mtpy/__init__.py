# -*- coding: utf-8 -*-
"""
==================
MTpy
==================

"""
# =============================================================================
# Commonly used objects
# =============================================================================
import sys
from loguru import logger

# common objects
from mtpy.core.mt import MT
from mtpy.core.mt_data import MTData
from mtpy.core.mt_collection import MTCollection
from mtpy.imaging.mtcolors import MT_CMAP_DICT, register_cmaps


__version__ = "2.0.0"
__all__ = ["MT", "MTData", "MTCollection"]

# =============================================================================
# Initiate loggers
# =============================================================================
config = {
    "handlers": [
        {
            "sink": sys.stdout,
            "level": "INFO",
            "colorize": True,
            "format": (
                "<level>{time:YY:MM:DDTHH:mm:ss} | {level: <3} | line:{line} |"
                "{name} | {function} | {message}</level>"
            ),
        },
    ],
    "extra": {"user": "someone"},
}
logger.configure(**config)
# logger.disable("mt_metadata")

# register custom colormaps
register_cmaps(MT_CMAP_DICT)
