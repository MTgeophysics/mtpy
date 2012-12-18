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


def model2vtkgrid():
    """
    Convert ModEM output files (model and responses) into 3D VTK resistivity grid

    Input:
    - ModEM model name
    - [optional] VTK resistivity grid file - output file name
    """

    arguments = sys.argv

    if len(arguments) < 2:
        sys.exit('ERROR - provide at least 1 file name: <modeldata file> ')

    try:
        ModEMmodelfn = os.path.abspath(os.path.realpath(arguments[1]))

        try:
            VTKfn = os.path.abspath(os.path.realpath(arguments[2]))
        except:
            VTKfn = os.path.abspath(os.path.realpath('VTKresistivitymodel.vtk'))


    except:
        sys.exit('ERROR - could not find file:\n%s'%(os.path.abspath(os.path.realpath(arguments[1]))))


    F = open(ModEMmodelfn, 'r')
    raw_data = f.readlines()
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
        for idx_north in range(n_north_blocks):
            #use North-to-South convention in the grid!!
            current_line = raw_data[current_line_idx]
            lo_data_tmp = current_line.strip().split()
            for idx_east in range(n_east_blocks):
                res_model[idx_north,idx_east, idx_depth] = exp(float(lo_data_tmp[idx_east]) )
            current_line_idx +=1

    #transfer grid to km  instead of m
    #use North-to-South convention in the grid!!
    N = np.array(coords_list[0].reverse())/1000.
    E = np.array(coords_list[1])/1000.
    D = np.array(coords_list[2])/1000.


    gridToVTK(VTKfn, N, E, D, cellData = {'resistivity' : res_model})
    print 'Created Resistivity File: ',VTKfn

    return VTKfn


def stations2vtkgrid():
    """
    Convert ModEM data file into 2D VTK station set (unstructured grid)

    Input:
    - ModEM data file name
    - [optional] VTK station grid file - output file name
    """

    arguments = sys.argv

    if len(arguments) < 2:
        sys.exit('ERROR - provide at least 1 file name: <data file> ')

        try:
            ModEMdatafn = os.path.abspath(os.path.realpath(arguments[1]))

            try:
                VTKfn = os.path.abspath(os.path.realpath(arguments[2]))
            except:
                VTKfn = os.path.abspath(os.path.realpath('VTKstations.vtk'))


        except:
            sys.exit('ERROR - could not find file:\n%s'%(os.path.abspath(os.path.realpath(arguments[1]))))


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


#if __name__ == '__main__':
    #arglist = sys.argv

    #main()
