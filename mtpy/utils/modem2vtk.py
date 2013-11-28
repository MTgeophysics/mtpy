#!/usr/bin/env python
# -*- coding: utf-8 -*-
######################################################################
#
# Create vtk format file from MODEM output model file
#   - changes units from meters to kilometers
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
    Convert ModEM output files (model) into 3D VTK resistivity grid

    Input:
    - ModEM model file 
    - ModEM data file 

    - [optional] VTK resistivity grid file - output file name
    - [optional] VTK station grid file - output file name

    """

    arguments = sys.argv

    if len(arguments) < 2:
        sys.exit('ERROR - provide at least 1 file name: <modeldata file>')

    try:
        Mmodel = os.path.abspath(os.path.realpath(arguments[1]))

        try:
            VTKresist = os.path.abspath(os.path.realpath(arguments[3]))
        except:
            VTKresist = os.path.abspath(os.path.realpath('VTKResistivityGrid'))

        try:
            Mdata = os.path.abspath(os.path.realpath(arguments[2]))
        except:
            pass

        try:
            VTKstations = os.path.abspath(os.path.realpath(arguments[4]))
        except:
            VTKstations = os.path.abspath(os.path.realpath('VTKStationGrid'))



    except:
        sys.exit('ERROR - could not find file(s)')


    f = open(Mmodel, 'r')

    # skip first line in file
    f.readline()

    # read N,E,D mesh dimensions
    dims = []
    modeldata_firstline = f.readline().split()

    try:
        function =  modeldata_firstline[4].lower()
    except:
        function = None

    for n in range(3):
        dims.append(int(modeldata_firstline[n]))
    size = dims[0]*dims[1]*dims[2]
    print 'Mesh:     ', dims
    print 'Datapoints:     ', size

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
        if len(modeldata_morelines) == 0:
            continue
        for j in modeldata_morelines:
            if function is None:
                mt[i] = float(j)
            elif function in ['loge']:
                mt[i] = np.e**(float(j))
            elif function in ['log','log10']:
                mt[i] = 10**(float(j))
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

    print 'Created Resistivity File: {0}.vtr'.format(VTKresist)

    try:
        f = open(Mdata, 'r')

        # get stations by parsing all lines that do NOT start with > or #
        rawdata = f.readlines()
        f.close()

        lo_datalines = []
        for l in rawdata:
            if l.strip()[0] not in ['#','>']: 
                lo_datalines.append(l.strip())
    
        lo_coords = []
        for line in lo_datalines:
            line = line.split()
            x = float(line[4])/1000.
            y = float(line[5])/1000.
            z = float(line[6])/1000.
            point = (x,y,z)
            if point not in lo_coords:
                lo_coords.append(point)

        all_coords = np.array(lo_coords)[:]
        print 'Stations:     ', len(all_coords)

        #sorting N->S, W->E
        all_coords = all_coords[np.lexsort((all_coords[:, 1], all_coords[:, 0]))]

        N = np.array(all_coords[:,0])
        E = np.array(all_coords[:,1])
        D = np.array(all_coords[:,2])


        dummy = np.ones(len(all_coords))

        pointsToVTK(VTKstations, N, E, D, data = {"dummyvalue" : dummy})



        print 'Created Station File: {0}.vtu'.format(VTKstations)
    except:
        print 'Could not create Station File '


if __name__ == '__main__':

    main()
