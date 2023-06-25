# -*- coding: utf-8 -*-
"""
==================
MTpy
==================

"""

__version__ = "2.0.0"

# load mtpy default logging config
from mtpy.utils.mtpy_logger import get_mtpy_logger

debug_logger = get_mtpy_logger(__name__, level="debug")
error_logger = get_mtpy_logger("error", level="error")
matplotlib_logger = get_mtpy_logger("matplotlib", level="warning")

# =============================================================================
# Commonly used objects
# =============================================================================
from mtpy.core.mt import MT
from mtpy.core.mt_data import MTData
from mtpy.core.mt_collection import MTCollection


__all__ = ["MT", "MTData", "MTCollection"]
