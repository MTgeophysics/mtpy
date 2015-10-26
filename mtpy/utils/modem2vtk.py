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
# The output of ModEM is not sorted straight forward!
# The x-components (North) of the model have to be read in opposite order.

# This holds for the model as well as the station locations.

######################################################################


def main():
    """
    Convert ModEM output files (model) into 3D VTK resistivity grid
    
    (distance is given in kilometers)


    Input:
    - ModEM model file 
    - ModEM data file 

    - [optional] VTK resistivity grid file - output file name
    - [optional] VTK station grid file - output file name

    """

    arguments = sys.argv

    if len(arguments) < 2:
        sys.exit('\nERROR - provide at least 1 file name: <model file> [<data>]'\
            ' [out:rho] [out:stations]\n')

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

    # read X,Y,Z mesh dimensions
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

    # read model values
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
    #reading the mesh blocks and calculate coordinates.
    #Coords are taken at the center points of the blocks

    #read the lines following the model:
    appendix = f.readlines()
    
    # the second last non empty line contains a triple that defines the 
    # lower left corner coordinates of the model 
    x0 = None
    for line in appendix[::-1]:
        line = line.strip().split()
        if len(line) != 3:
            continue
        else:
            try:
                x0 = float(line[0])/1000.
                y0 = float(line[1])/1000.
                z0 = float(line[2])/1000.
                break
            except:
                continue
    if x0 is None:
        x0 = 0
        y0 = 0
        z0 = 0 
        print 'Warning - no reference point found - lower left corner of model'\
        ' set to 0,0,0'

    Xspacing = spacing[:dims[0]]

    # calc North coordinates of vtk mesh
    Xdist = 0 # calculate total North distance by summing all blocks
    for i in range(dims[0]):
        Xdist += Xspacing[i]

    X = np.zeros(dims[0])
    #X[0] = -0.5 * Xdist + 0.5*(Xspacing[0])# define zero at the center of the model
    X[0] = 0.5*(Xspacing[0])# define zero at the lower left corner of the model
    for i in range(dims[0]-1):
        local_spacing = 0.5*(Xspacing[i]+Xspacing[i+1])
        X[i+1] = X[i] + local_spacing
    
    #add reference point coordinate
    X += x0

    Yspacing = spacing[dims[0]:dims[0]+dims[1]]
    # calc Yast coordinates of vtk mesh
    Ydist = 0 # calculate total Yast distance by summing all blocks
    for i in range(dims[1]):
        Ydist += Yspacing[i]

    Y = np.zeros(dims[1])
    #Y[0] = -0.5 * Ydist + 0.5*(Yspacing[0])# define zero at the center of the model
    Y[0] = 0.5*(Yspacing[0])# define zero at the lower left corner of the model
    for i in range(dims[1]-1):
        local_spacing = 0.5*(Yspacing[i]+Yspacing[i+1])
        Y[i+1] = Y[i] + local_spacing
    
    #add reference point coordinate
    Y += y0

    Dspacing=spacing[dims[0]+dims[1]:]
    # calc Down coordinates of vtk mesh
    D = np.zeros(dims[2])
    D[0] = 0.5*Dspacing[0]
    for i in range(dims[2]-1):
        local_spacing = 0.5*(Dspacing[i]+Dspacing[i+1])
        D[i+1] = D[i] + local_spacing

    
    #add reference point coordinate
    D+= z0

    # output to vtk format
    # first components read in reverse order!!
    mtNS = np.zeros((dims[0],dims[1],dims[2])) 

    n=0

    #define array for alternative output as csv file, which
    #can be read by paraview directly 
    arr = np.zeros((len(X)*len(Y)*len(D),4))
    
    for idx_D in range(dims[2]):
        for idx_Y in range(dims[1]):
            for idx_X in range(dims[0]):
                #switch North/South by convention
                mtNS[-(idx_X+1),idx_Y,idx_D] = mt[n]
                arr[n,0] = X[idx_X]
                arr[n,1] = Y[idx_Y]
                arr[n,2] = D[idx_D]
                arr[n,3] = mt[n]
                n += 1
    gridToVTK(VTKresist, X, Y, D, pointData = {'resistivity' : mtNS})

    f.close()

    #fn2 = VTKresist+'.csv'
    #F = open(fn2,'w')
    #F.write('X , Y , Down , Rho \n')
    #np.savetxt(F,arr,delimiter=',')
    #F.close()

    #print 'Created Resistivity Array: {0}.csv'.format(VTKresist)
    print 'Created Resistivity VTK File: {0}.vtr'.format(VTKresist)

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



        print 'Created Station VTK File: {0}.vtu'.format(VTKstations)
    except:
        print 'Could not create Station File '


if __name__ == '__main__':

    main()
