# -*- coding: utf-8 -*-
"""
Created on 12.12.2012

@author: LK
"""

import os
import os.path as op
import sys
import numpy as np
import glob

import mtpy.core.mttools as mtt
import mtpy.modeling.winglinktools as wlt



def winglinkmesh2modelfile(WLoutputfile, modelfilename= 'ModEM_initmodel', res_value=100):
    """return init3d file/start model
    
    mainly copied from ws3dtools...
    
    Inputs:
        WLoutputfile -  *.out-file from winglink
        modelfilename - name of the modelfile 
        res_value - starting homogeneous half space in Ohm-meters (default: 100)
        
        
    Output:
        init file path
        
    """
    
     #create the output filename
    model_fn = op.abspath(modelfilename)

        
    dx,dy,dz=wlt.readWLOutFile(WLoutputfile,ncol=5)
    
    nx=len(dx)
    ny=len(dy)
    nz=len(dz)
    
    init_modelFH = open(model_fn,'w')
    init_modelFH.write('#Initial halfspace model, based on WingLink generated mesh \n')
    init_modelFH.write('%i %i %i 1 \n'%(ny,nx,nz))
        
    #write y locations
    y_string=''
    y_counter=0 
    for y_idx in range(ny):
        y_string += '%.3e  '%(dy[y_idx])
        y_counter+=1
        if y_counter == 8:
            y_string += '\n'
            y_counter = 0
    if ny%8:
        y_string +='\n'
    init_modelFH.write(y_string)
    
    #write x locations
    x_string=''
    x_counter=0 
    for x_idx in range(nx):
        x_string += '%.3e  '%(dx[x_idx])
        x_counter+=1
        if x_counter == 8:
            x_string += '\n'
            x_counter = 0
    if nx%8:		    
	x_string +='\n'
    init_modelFH.write(x_string)

    #write z locations
    z_string=''
    z_counter=0 
    for z_idx in range(nz):
        z_string += '%.3e  '%(dz[z_idx])
        z_counter+=1
        if z_counter == 8:
            z_string += '\n'
            z_counter = 0   
    if nz%8:
        z_string +='\n'
    init_modelFH.write(z_string)
           
    init_modelFH.write('%i \n'%int(rhostart))
    
    init_modelFH.close()
    
    print 'Wrote initial halfspace model to file: %s '%(model_fn)

    
    return ifile
        

def  edis2datafile(edilist, comment):

    datafilename = op.abspath('ModEM_inputdata')
    if len(comment)>100:
        sys.exit('comment string is too long (cannot exceed 100 characters)\n')

    for idx_edi,edifile in enumerate(edilist):
        #generate dictionary containing all info
        
        
        
    #check dictionary for values n_stations, n_periods:
    
    #write header info
    F = open(datafilename,'w')
    F.write('#%s\n'%(comment))
    F.write('#Period Station Lat Lon X Y Z Component Real Imag Err\n')
    F.write('>Full_Impedance\n')
    F.write('>exp(-i\omega t\n')    
    F.write('>[mV/km]/[nT]\n')    
    F.write('>0.00\n')    
    F.write('>0 0 \n')    
    F.write('>%i %i\n'%(no_periods, no_stations))        
        
    #iterate over dictionary and write data file lines successively:   
    
    
    
 
    
    return datafilename
    
    
    
    
def generate_edilist(edifolder):
    
    lo_edifiles = [op.abspath(i) for i in glob.glob('*.[eE][dD][iI]')]
    
    return edilist 

    


def winglink2modem(edifolder, winglinkoutput, resistivity):
    
    if not op.isdir(edifolder):
        sys.exit('cannot find EDI files: no such directory: \n%s'%(edifolder))   

    edilist = generate_edilist(edifolder)       

    if len(edilist)==0:
        sys.exit('cannot find EDI files in given directory: \n%s'%(edifolder))          
    
    try:
        WL_outfile = op.abspath(winglinkoutput)
    except:
        sys.exit('cannot find specified WingLink output file: \n%s'%(winglinkoutput))
    
    try:
        HS_rho_value = float(resistivity)
    except:
        sys.exit('provided resistivity value is not a proper number')


    datafn  = edis2datafile(edilist)
    
    modelfn = winglinkmesh2modelfile(WL_outfile, modelfilename=chosen_model, res_value=HS_rho_value)

    return datafn, modelfn

