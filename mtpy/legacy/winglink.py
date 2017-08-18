# -*- coding: utf-8 -*-
"""
===================
Winglink Tools
===================

    * module to deal with WingLink output files
    
    
Created on Mon Aug 26 17:44:20 2013

@author: jpeacock-pr
"""
#==============================================================================
import numpy as np
#==============================================================================

#==============================================================================
# read an out file
#==============================================================================
def read_out_file(out_fn, ncol=5):
    """
    read .out file from winglink
    
    Arguments:
    -----------
        **out_fn** : full path to .out file from winglink
        
    Returns:
    ---------
        **dx, dy, dz** : np.ndarrays
                        cell nodes in x, y, z directions
                        (note x is to the East here and
                          y is to the north in meters.)
    """
    
    wl_ofid = file(out_fn,'r')
    raw_data = wl_ofid.read().strip().split()
    
    nx = int(raw_data[0])
    ny = int(raw_data[1])
    nz = int(raw_data[2])


    dx = np.zeros(nx)
    dy = np.zeros(ny)
    dz = np.zeros(nz)
    
    for xx in range(nx):
      dx[xx] = raw_data[xx + ncol]
    for yy in range(ny):
      dy[yy] = raw_data[yy+ncol+nx]
    for zz in range(nz):
      dz[zz] = raw_data[zz+ncol+nx+ny]
            
    return dx, dy, dz

#==============================================================================
# read a sites file
#==============================================================================
def read_sites_file(sites_fn):
    """
    read sites_ file output from winglink
    
    Arguments:
    -----------
        **sites_fn** : string
                       full path to the sites file output by winglink
        
    Returns:
    ----------
        **slst** : list of dictionaries for each station.  
                   Keys include:
                       * station = station name
                       * dx = number of blocks from center of grid in 
                              East-West direction
                       * dy = number of blocks from center of grid in
                              North-South direction
                       * dz = number of blocks from center of grid 
                              vertically
                       * number = block number in the grid
        
        **site_list** : list of station names 
    """
    
    sfid = file(sites_fn,'r')
    slines = sfid.readlines()
    
    slst = []
    site_list = []
    for ss in slines:
        sdict = {}
        sline = ss.strip().split()
        sdict['station'] = sline[0][0:-4]
        sdict['dx'] = int(sline[1])-1
        sdict['dy'] = int(sline[2])-1
        sdict['dz'] = int(sline[3])-1
        sdict['something'] = int(sline[4])
        sdict['number'] = int(sline[5])
        slst.append(sdict)
        site_list.append(sline[0][0:-4])
    return slst, site_list
    
#==============================================================================
# get station locations from sites file
#==============================================================================
def get_station_locations(sites_fn, out_fn, ncol=5):
    """
    get x (e-w) and y (n-s) position of station and put in middle of cell for 
    a 3D model.
    
    Arguments:
    -----------
        **sites_fn** : string
                       full path to sites file output from winglink
        
        **out_fn** : string
                     full path to .out file output from winglink
        
        **ncol** : int
                   number of columns the data is in
                   *default* is 5
        
    Returns:
    ---------
        **xarr** : np.ndarray()
                  array of relative distance for each station from center of
                  the grid.  Note this is E-W direction
        **yarr** : np.ndarray()
                   array of relative distance for each station from center of 
                   the grid.  Note this is N-S direction
                
    """
    
    slst, sitelst = read_sites_file(sites_fn)
    
    dx, dy, dz = read_out_file(out_fn, ncol=ncol)
    
    ns = len(slst)
    nxh = len(dx)/2
    nyh = len(dy)/2
    xarr = np.zeros(ns)
    yarr = np.zeros(ns)
    
    
    for ii,sdict in enumerate(slst):
        xx = sdict['dx']
        yy = sdict['dy']
        if xx < nxh:
            xarr[ii] = dx[xx:nxh].sum()-dx[xx]/2
        else:
            xarr[ii] = dx[nxh:xx].sum()+dx[xx]/2                    
        if yy < nyh:
            yarr[ii] = -1*(dy[yy:nyh].sum()-dy[yy]/2)
        else:
            yarr[ii] = -1*(dy[nyh:yy].sum()+dy[yy]/2)   

    return xarr, yarr  
