#!/usr/bin/env python
# -*- coding: utf-8 -*-
######################################################################
#
# Create vtk format file from WSMT output model file
#   - changes units from meters to kilometers
#
# 15.05.2012
# Institute of Earth Science and Engineering
# J Rugis
#
######################################################################
#

from evtk.hl import gridToVTK, pointsToVTK
import numpy as np
import os

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
    #WSMTmodel = './ParwQuantec_01_model.01' # WSMT output model file
    #WSMTresp = './ParwQuantec_01_resp.01'   # WSMT initial response file
    #VTKresist = './ParwQuantec_01_01_res'            # VTK file to create
    #VTKstations = './ParwQuantec_01_0_sta'          # VTK file to create

    #==============================================================================
    # Paralana Regional
    #==============================================================================
    #dirpath=r"/home/mt/Documents/wsinv/Paralana/Big3D/Inv1rough"
    #savepath=r"/home/mt/Documents/ParaviewFiles/Paralana"
    ##
    ## WSMT output model file
    #WSMTmodel = os.path.join(dirpath,'ParalanaRegional_rough_model.05')
    #
    ## WSMT initial response file
    #WSMTresp = os.path.join(dirpath,'ParalanaRegional_rough_resp.05')
    #
    ## VTK file to create
    #VTKresist = os.path.join(savepath,'ParalanaRegionalR5_res')
    #
    ## VTK file to create
    #VTKstations = os.path.join(savepath,'ParalanaRegionalR5_sta')

    #==============================================================================
    # Paralana Dense
    #==============================================================================
    dirpath='.'#r"/home/mt/Documents/wsinv/od/inv0503rough"
    savepath='.'#r"/home/mt/Documents/ParaviewFiles/od"

    # WSMT output model file
    WSMTmodel = os.path.join(dirpath,'olympic_model.03')

    # WSMT initial response file
    WSMTresp = os.path.join(dirpath,'olympic_resp.03')

    # VTK file to create
    VTKresist = os.path.join(savepath,'olympicR3_res')

    # VTK file to create
    VTKstations = os.path.join(savepath,'olympicR3_sta')

    #==============================================================================
    # Olympic Dam
    #==============================================================================
    #dirpath=r"/home/mt/Documents/wsinv/od/inv0503"
    #savepath=r"/home/mt/Documents/ParaviewFiles/od"
    #
    #if not os.path.exists(savepath):
    #    os.mkdir(savepath)
    #
    ## WSMT output model file
    #WSMTmodel = os.path.join(dirpath,'olympic_model.05_01')
    #
    ## WSMT initial response file
    #WSMTresp = os.path.join(dirpath,'olympic_resp.05_01')
    #
    ## VTK file to create for resistivity blocks
    #VTKresist = os.path.join(savepath,'olympic_res')
    #
    # # VTK file to create for station locations
    #VTKstations = os.path.join(savepath,'olympic_sta')
    #####################################################################

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
            for j in range(len(list)):
                spacing.append(float(modeldata_nextlines[j])/1000.0)
                i += 1

    # read mt data
    #  (depends on line break only after final value)
    mt = np.zeros(size)
    i=0
    while i < size:
        modeldata_morelines = f.readline().split()
        for j in range(len(modeldata_morelines)):
            mt[i] = float(modeldata_morelines[j])
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
    gridToVTK(VTKresist, V, E, D, cellData = {'resistivity' : mtNS})

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

    # set Depths -- all stations at the surface!!
    D = np.zeros(nstations)

    # output to vtk format - dummy value for scalar field needed
    dummy = np.zeros(nstations)
    for j in range(nstations):
        dummy[j] = 1.0
    pointsToVTK(VTKstations, x, y, z, data = {"value" : dummy})

    f.close()


    print 'Created Resistivity File: ',VTKresist
    print 'Created Station File: ',VTKstations


if __name__ == '__main__':
  main()