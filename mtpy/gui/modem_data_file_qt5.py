# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 12:42:36 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# ==============================================================================
# Imports
# ==============================================================================
import sys
from pathlib import Path

try:
    from PyQt5 import QtCore, QtWidgets
except ImportError:
    raise ImportError("This version needs PyQt5")


class ModEMDataFile(QtWidgets.QWidget):
    """
    Widget to build a data file
    
    """

    pass
