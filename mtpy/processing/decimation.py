#!/usr/bin/env python

"""
mtpy/processing/decimation.py

Functions for the decimation of raw time series. 



For calling a batch decimation rather than just one file, use the appropriate scripts from the mtpy.utils subpackage. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import sys, os
import os.path as op
import time
import copy


import  mtpy.utils.exceptions as MTex

#=================================================================