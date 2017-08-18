# -*- coding: utf-8 -*-
"""
Created on Thu Jun 05 14:52:33 2014

@author: a1655681
"""

import numpy as np
import scipy.interpolate as si
import os
import mtpy.utils.exceptions as MTex

def get_elevation(x,y,elevfn,skiprows = 1):
    """
    get elevation at an arbitrary point, interpolated from an xyz file containing elevations
    x = x location of point, can be a numpy array to return multiple z values
    y = y location of point, can be a numpy array to return multiple z values
    elevfn = full path to elevation filename
    """
    elev = np.loadtxt(elevfn)
    f = si.LinearNDInterpolator(elev[:,0:2],elev[:,2])
    return f(x,y)
    
    
def project_interface(interface,epsg_from,epsg_to,suffix,skiprows=1):
    """
    project interface, save into a new file with suffix
    interface = full path to xyz file containing x,y,z values
    epsg_from = epsg from
    epsg_to = epsg to
    suffix = suffix to add to new files
    
    useful epsg numbers:
    28354 - GDA 94 mga zone 54
    4326 - WGS 84
    """
    try:
        import pyproj
    except:
        raise MTex.MTpyError_module_import('pyproj module not found; cannot continue')
        
    coord_from = pyproj.Proj("+init=EPSG:%i"%epsg_from)
    coord_to = pyproj.Proj("+init=EPSG:%i"%epsg_to)
    
    data = np.loadtxt(interface,skiprows=skiprows)
    xp,yp = pyproj.transform(coord_from,coord_to,data[:,0],data[:,1])
    
    filename,extension = os.path.splitext(interface)
    outfile = filename + suffix + extension
    np.savetxt(outfile,np.vstack([xp,yp,data[:,2]]).T,fmt = ['%12.6f','%12.6f','%10.2f'])