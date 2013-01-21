#!/usr/bin/env python

"""
This modules contains helper functions for the calibration of raw time series. 

The various functions deal with the calibration of data from fluxgates, coils, dipoles,...
The calibration depends on  the instrument as well as on the respective data logger. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import calendar
import time


from mtpy.utils.exceptions import *
#=================================================================




#=================================================================
# section for amplification factors



#=================================================================


