# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 15:50:53 2016

@author: Alison Kirkby

read and write gocad objects

"""

import numpy as np
import os.path as op

class Sgrid():
    """
    create a gocad sgrid file
    
    need to provide:
    workdir = working directory
    fn = filename for the sgrid
    resistivity = 3d numpy array containing resistivity values, shape (ny,nx,nz)
    grid_xyz = tuple containing x,y,z locations of edges of cells for each 
               resistivity value. Each item in tuple has shape (ny+1,nx+1,nz+1)
    
    """
    
    def __init__(self, workdir='.', **kwargs):
        self.workdir = workdir
        self.fn = kwargs.pop('fn','model')
        self.property_fn = self.fn + '__ascii@@'
        self.property_name = 'Resistivity'
        self.grid_xyz = kwargs.pop('grid_xyz',None)
        self.resistivity = kwargs.pop('resistivity',None)
        self.no_data_value = -99999
        
    def _write_header(self):
                
        ny,nx,nz = np.array(self.resistivity.shape) + 1

        headerlines = [r''+ item + '\n' for item in ['GOCAD SGrid 1 ',
        'HEADER {',
        'name:{}'.format(op.basename(self.fn)),
        'ascii:on',
        'double_precision_binary:off',
        '}',
        'GOCAD_ORIGINAL_COORDINATE_SYSTEM',
        'NAME Default',
        'AXIS_NAME "X" "Y" "Z"',
        'AXIS_UNIT "m" "m" "m"',
        'ZPOSITIVE Elevation',
        'END_ORIGINAL_COORDINATE_SYSTEM',
        'AXIS_N {} {} {} '.format(ny,nx,nz),
        'PROP_ALIGNMENT CELLS',
        'ASCII_DATA_FILE {}'.format(op.basename(self.property_fn)),
        '',
        '',
        'PROPERTY 1 "{}"'.format(self.property_name),
        'PROPERTY_CLASS 1 "{}"'.format(self.property_name),
        'PROPERTY_KIND 1 "Resistivity"',
        'PROPERTY_CLASS_HEADER 1 "{}" '.format(str.lower(self.property_name))+'{',
        'low_clip:1',
        'high_clip:10000',
        'pclip:99',
        'colormap:flag',
        'last_selected_folder:Property',
        'scale_function:log10',
        '*colormap*reverse:true',
        '}',
        'PROPERTY_SUBCLASS 1 QUANTITY Float',
        'PROP_ORIGINAL_UNIT 1 ohm*m',
        'PROP_UNIT 1 ohm*m',
        'PROP_NO_DATA_VALUE 1 {}'.format(self.no_data_value),
        'PROP_ESIZE 1 4',
        'END']]        
        
        hdrfn = self.fn
        if not hdrfn.endswith('.sg'):
            hdrfn += '.sg'
            
        with open(hdrfn,'w') as hdrfile:
            hdrfile.writelines(headerlines)
  
          
    def _write_data(self):
        
        resmodel = np.ones(np.array(self.resistivity.shape)+1)*self.no_data_value
        resmodel[:-1,:-1,:-1] = self.resistivity        
        resvals = resmodel.transpose(2,1,0).flatten()

        x,y,z = [arr.transpose(2,1,0).flatten() for arr in self.grid_xyz]
        ny,nx,nz = self.grid_xyz[0].shape
        j,i,k = [arr.transpose(2,1,0).flatten() for arr in \
        np.meshgrid(*[np.arange(ll) for ll in [nx,ny,nz]])]

        # make an array containing the data
        data = np.vstack([x,y,z,resvals,i,j,k]).T
        
        # make data header
        datahdr = '\n X Y Z {} I J K\n'.format(self.property_name)
        
        # write property values
        np.savetxt(self.property_fn,data,header=datahdr,comments='*',fmt=['%10.3f']*4 + ['%10i']*3)

        
    def write_sgrid_file(self):
        self._write_header()
        self._write_data()
                