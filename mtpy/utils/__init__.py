# -*- coding: utf-8 -*-
"""
Got rid of the GDAL check because have moved geographic operations to use 
pyproj.  There are still some functions that use GDAL, but those are for
raster and shapefile making.  Therefore the check is not needed.

Created on Tue Sep  5 14:35:54 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np

# =============================================================================

### constants magnetic permeability
MU0 = 4e-7 * np.pi

EPSILON = epsilon = 1e-9
