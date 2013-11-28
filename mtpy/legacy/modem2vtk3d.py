#!/usr/bin/env python
# -*- coding: utf-8 -*-
######################################################################
#
# Create 3D vtk format file from ModEM output model file
#   - changes units from meters to kilometers
#
# 16.12.2012
# EES
# LK
#
######################################################################
#

from evtk.hl import gridToVTK, pointsToVTK
import numpy as np
import os
import sys
import mtpy.modeling.modemtools as mdt
import numpy as np


######################################################################

## Orientation convention:
#
# coordinate system NED is used! First component is positive from South to North,
# second component is positve from West to East, third component is positive Downwards#


######################################################################


def model2vtkgrid(ModEMmodelfn,VTKfn='VTKresistivitymodel' ):
    """
    Convert ModEM output files (model and responses) into 3D VTK resistivity grid

    Input:
    - ModEM model name
    - [optional] VTK resistivity grid file - output file name
    """

    if not os.path.isfile(os.path.abspath(os.path.realpath(ModEMmodelfn))):
        sys.exit('ERROR - could not find file:\n%s'%(os.path.abspath(os.path.realpath(ModEMmodelfn))))

    if not os.path.isfile(os.path.abspath(os.path.realpath(VTKfn))):
        VTKfn = os.path.abspath(os.path.realpath('VTKresistivitymodel'))



    F = open(ModEMmodelfn, 'r')
    raw_data = F.readlines()
    F.close()

    coords_list = mdt.getmeshblockcoordinates(ModEMmodelfn)
    n_north_blocks = len(coords_list[0])
    n_east_blocks  = len(coords_list[1])
    n_depth_blocks = len(coords_list[2])

    res_model = np.zeros((n_north_blocks,n_east_blocks,n_depth_blocks))


    current_line_idx = 5
    for idx_depth in range(n_depth_blocks):
        #catering for the empty line between depths
        current_line_idx += 1
        for idx_east in range(n_east_blocks):
            #use North-to-South convention in the grid!!
            current_line = raw_data[current_line_idx]
            lo_data_tmp = current_line.strip().split()
            for idx_north in range(n_north_blocks):
                res_model[idx_north,idx_east, idx_depth] = np.exp(float(lo_data_tmp[idx_east]) )
            current_line_idx +=1

    #transfer grid to km  instead of m
    #use North-to-South convention in the grid!!
    coords_list[0].reverse()
    N = np.array(coords_list[0])/1000.
    E = np.array(coords_list[1])/1000.
    D = np.array(coords_list[2])/1000.


    gridToVTK(VTKfn, N, E, D, cellData = {'resistivity (in Ohm)' : res_model})

    print 'Created Resistivity File: ',VTKfn

    return VTKfn


def data2vtkstationsgrid(ModEMdatafn, VTKfn='VTKstations'):
    """
    Convert ModEM data file into 2D VTK station set (unstructured grid)

    Input:
    - ModEM data file name
    - [optional] VTK station grid file - output file name
    """

    if not os.path.isfile(os.path.abspath(os.path.realpath(ModEMdatafn))):
        sys.exit('ERROR - could not find file:\n%s'%(os.path.abspath(os.path.realpath(ModEMdatafn))))

    if not os.path.isfile(os.path.abspath(os.path.realpath(VTKfn))):
        VTKfn = os.path.abspath(os.path.realpath('VTKstations'))


    F = open(ModEMdatafn, 'r')
    raw_data = F.readlines()
    F.close()

    lo_stations   = []
    stations_dict = {}
    lo_norths= []
    lo_easts = []


    for tmp_line in raw_data:
        linelist = tmp_line.strip().split()
        if len(linelist) == 11:
            tmp_station = linelist[1]
            tmp_north  = float(linelist[4])
            tmp_east   = float(linelist[5])
            if not tmp_station in lo_stations:
                lo_stations.append(tmp_station)
                stations_dict[tmp_station] = [tmp_north,tmp_east]
                lo_norths.append(tmp_north)
                lo_easts.append(tmp_east)

    #convert m to km
    N = np.array(lo_norths)/1000.
    E = np.array(lo_easts)/1000.
    D = np.zeros((len(lo_norths)))

    #dummy scalar values
    dummy = np.ones((len(lo_norths)))

    pointsToVTK(VTKfn, N, E, D, data = {"value" : dummy})


    print 'Created Station File: ',VTKfn


