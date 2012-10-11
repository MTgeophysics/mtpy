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
  dirpath=r"/home/mt/Documents/wsinv/od/inv0503rough"
  savepath=r"/home/mt/Documents/ParaviewFiles/od"
  
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
  
  # read x,y,z mesh dimensions
  dims = []
  modeldata_firstline = f.readline().split()
  for n in range(3):
    dims.append(int(modeldata_firstline[n]))
  size = dims[0]*dims[1]*dims[2]
  print 'Mesh     ', dims
  print 'Data     ', size

  # read x,y,z spacing
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
  
  # calc x coordinates of vtk mesh
  xdist = 0 # calculate total x distance
  for i in range(dims[0]):
    xdist += spacing[i]
  x = np.zeros(dims[0]+1)
  x[0] = -0.5 * xdist # zero center of model
  for i in range(dims[0]):
    x[i+1] = x[i] + spacing[i]
  
  # calc y coordinates of vtk mesh
  ydist = 0 # calculate total y distance
  for i in range(dims[1]):
    ydist += spacing[dims[0] + i]
  y = np.zeros(dims[1]+1)
  y[0] = -0.5 * ydist # zero center of model
  for i in range(dims[1]):
    y[i+1] = y[i] + spacing[dims[0] + i]
  
  # calc z coordinates of vtk mesh
  z = np.zeros(dims[2]+1)
  z[0] = 0.0
  for i in range(dims[2]):
    z[i+1] = z[i] + spacing[dims[0] + dims[1] + i]
  
  # output to vtk format
  mtNS = np.zeros((dims[0],dims[1],dims[2])) # North-to-South conversion
  n=0
  for k in range(dims[2]):
    for j in range(dims[1]):
      for i in range(dims[0]):
        mtNS[(dims[0]-1)-i,j,k] = mt[n]
        n += 1
  gridToVTK(VTKresist, x, y, z, cellData = {'resistivity' : mtNS})
  
  f.close()
  
  f = open(WSMTresp, 'r')
  
  # get station count
  respdata_firstline = f.readline().split()
  nstations = int(respdata_firstline[0])
  print 'Stations ', nstations
  
  # read x locations
  f.readline() #skip line
  x = np.zeros(nstations)
  i=0
  while i < nstations:
    respdata_nextlines = f.readline().split()
    for j in range(len(respdata_nextlines)):
      x[i] = float(respdata_nextlines[j])/1000.0
      i += 1
  
  # read y locations
  f.readline() #skip line
  y = np.zeros(nstations)
  i=0
  while i < nstations:
    respdata_morelines = f.readline().split()
    for j in range(len(respdata_morelines)):
      y[i] = float(respdata_morelines[j])/1000.0
      i += 1
  
  # set z locations
  z = np.zeros(nstations)
  
  # output to vtk format
  dummy = np.zeros(nstations)
  for j in range(nstations):
    dummy[j] = 1.0
  pointsToVTK(VTKstations, x, y, z, data = {"value" : dummy})
  
  f.close()
  
  
  print 'Created Resistivity File: ',VTKresist
  print 'Created Station File: ',VTKstations
  

if __name__ == '__main__':
  main()