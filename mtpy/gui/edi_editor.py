# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 15:04:27 2016

@author: jpeacock
"""

from PyQt4 import QtCore, QtGui
import mtpy.modeling.modem_new as modem
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as plt
import os


class Edit_EDI_Window(object):
    """
    This is the main window for editing an edi file.
    
    Includes:
        * Editing points by removing them
        * Adjusting for static shift
        * Removing distortion
    """