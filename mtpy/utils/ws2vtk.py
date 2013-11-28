#!/usr/bin/env python
# -*- coding: utf-8 -*-
######################################################################
#
# Create vtk format file from WS3dinv output model file
#   - changes units from meters to kilometers
#

# LK & JP UofA
# Geophysics
# 2013

# adapted from
# 15.05.2012
# Institute of Earth Science and Engineering
# J Rugis
#
######################################################################
#

from evtk.hl import gridToVTK, pointsToVTK
import numpy as np
import os
import sys

######################################################################

## Orientation convention:
#
# coordinate system NED is used! First component is positive from South to North,
# second component is positve from West to East, third component is positive Downwards#

#Important note:
# The output of wsinv3D is not sorted straight forward!
# The x-components (North) of the model have to be read in opposite order.

# This holds for the model as well as the station locations.

######################################################################


def main():
    """
    Convert ws3Dinv output files (model and responses) into 3D VTK resistivity grid
    and unstructured VTKGrid containing station locations.

    Input:
    - ws3dInv model name
    - ws3DInv response file name
    - [optional] VTK resistivity grid file - output file name
    - [optional] VTK stations grid file - output file name
    """

    arguments = sys.argv

    if len(arguments) < 3:
        sys.exit('ERROR - provide at least 2 file names: <modeldata file>  <responses file>')

    try:
        WSMTmodel = os.path.abspath(os.path.realpath(arguments[1]))
        WSMTresp  = os.path.abspath(os.path.realpath(arguments[2]))

        try:
            VTKresist = os.path.abspath(os.path.realpath(arguments[3]))
        except:
            VTKresist = os.path.abspath(os.path.realpath('VTKResistivityGrid'))

        try:
            VTKstations = os.path.abspath(os.path.realpath(arguments[4]))
        except:
            VTKstations = os.path.abspath(os.path.realpath('VTKStationGrid'))


    except:
        sys.exit('ERROR - could not find file(s)')


    f = open(WSMTmodel, 'r')

    # skip first line in file
    f.readline()

    # read N,E,D mesh dimensions
    dims = []
    modeldata_firstline = f.readline().split()
    for n in range(3):
        dims.append(int(modeldata_firstline[n]))
    size = dims[0]*dims[1]*dims[2]
    print 'Mesh     ', dims
    print 'Data     ', size

    # read N,E,D spacing
    #  (depends on line break only after final value)
    spacing = []
    for n in range(3):
        i=0
        while i < dims[n]:
            modeldata_nextlines = f.readline().split()
            for j in modeldata_nextlines:
                spacing.append(float(j)/1000.0)
                i += 1

    # read mt data
    #  (depends on line break only after final value)
    mt = np.zeros(size)
    i=0
    while i < size:
        modeldata_morelines = f.readline().split()
        for j in modeldata_morelines:
            mt[i] = float(j)
            i += 1

    # calc North coordinates of vtk mesh
    Ndist = 0 # calculate total North distance
    for i in range(dims[0]):
        Ndist += spacing[i]
    N = np.zeros(dims[0]+1)
    N[0] = -0.5 * Ndist # zero center of model
    for i in range(dims[0]):
        N[i+1] = N[i] + spacing[i]

    # calc East coordinates of vtk mesh
    Edist = 0 # calculate total y distance
    for i in range(dims[1]):
        Edist += spacing[dims[0] + i]
    E = np.zeros(dims[1]+1)
    E[0] = -0.5 * Edist # zero center of model
    for i in range(dims[1]):
        E[i+1] = E[i] + spacing[dims[0] + i]

    # calc Down coordinates of vtk mesh
    D = np.zeros(dims[2]+1)
    D[0] = 0.0
    for i in range(dims[2]):
        D[i+1] = D[i] + spacing[dims[0] + dims[1] + i]

    # output to vtk format
    # first components read in reverse order!!
    mtNS = np.zeros((dims[0],dims[1],dims[2])) # North-to-South conversion
    n=0
    for idx_D in range(dims[2]):
        for idx_E in range(dims[1]):
            for idx_S in range(dims[0]):
                mtNS[(dims[0]-1)-idx_S,idx_E,idx_D] = mt[n]
                n += 1
    gridToVTK(VTKresist, N, E, D, cellData = {'resistivity' : mtNS})

    f.close()

    f = open(WSMTresp, 'r')

    # get station count
    respdata_firstline = f.readline().split()
    nstations = int(respdata_firstline[0])
    print 'Stations ', nstations

    # read North locations
    f.readline() #skip line
    N = np.zeros(nstations)
    i=0
    while i < nstations:
        respdata_nextlines = f.readline().split()
        for j in respdata_nextlines:
            N[i] = float(j)/1000.0
            i += 1

    # read East locations
    f.readline() #skip line
    E = np.zeros(nstations)
    i=0
    while i < nstations:
        respdata_morelines = f.readline().split()
        for j in respdata_morelines:
            E[i] = float(j)/1000.0
            i += 1
    f.close()


    # set Depths -- all stations at the surface!!
    D = np.zeros(nstations)

    # output to vtk format - dummy value for scalar field needed
    dummy = np.ones(nstations)
    #for j in range(nstations):
    #    dummy[j] = 1.0

    print np.shape(dummy),np.shape(N),np.shape(E),np.shape(D)

    pointsToVTK(VTKstations, N, E, D, data = {"dummyvalue" : dummy})



    print 'Created Resistivity File: {0}.vtr '.format(VTKresist)
    print 'Created Station File: {0}.vtu '.format(VTKstations)


if __name__ == '__main__':

    main()
