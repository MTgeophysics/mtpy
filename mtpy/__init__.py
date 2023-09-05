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
            "format": "<level>{time} | {level: <3} | {name} | {function} | {message}</level>",
        },
    ],
    "extra": {"user": "someone"},
}
logger.configure(**config)
# logger.disable("mt_metadata")
