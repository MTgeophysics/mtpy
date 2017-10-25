# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 09:35:31 2017

@author: Alison Kirkby

functions to assist with mesh generation

"""

import numpy as np



def make_log_increasing_array(z1_layer, target_depth, n_layers, increment_factor=0.9):
    """
    create depth array with log increasing cells, down to target depth,
    inputs are z1_layer thickness, target depth, number of layers (n_layers)
    """        
    
    # make initial guess for maximum cell thickness
    max_cell_thickness = target_depth
    # make initial guess for log_z
    log_z = np.logspace(np.log10(z1_layer), 
                        np.log10(max_cell_thickness),
                        num=n_layers)
    counter = 0
    
    while np.sum(log_z) > target_depth:
        max_cell_thickness *= increment_factor
        log_z = np.logspace(np.log10(z1_layer), 
                            np.log10(max_cell_thickness),
                            num=n_layers) 
        counter += 1
        if counter > 1e6:
            break        

    return log_z


def get_padding_cells(cell_width, max_distance, num_cells, stretch):
    """
    get padding cells, which are exponentially increasing to a given 
    distance.  Make sure that each cell is larger than the one previously.
    
    Arguments
    -------------
    
        **cell_width** : float
                         width of grid cell (m)
                         
        **max_distance** : float
                           maximum distance the grid will extend (m)
                           
        **num_cells** : int
                        number of padding cells
                        
        **stretch** : float
                      base geometric factor
                        
    Returns
    ----------------
    
        **padding** : np.ndarray
                      array of padding cells for one side
    
    """

    # compute scaling factor
    scaling = ((max_distance)/(cell_width*stretch))**(1./(num_cells-1)) 
    
    # make padding cell
    padding = np.zeros(num_cells)
    for ii in range(num_cells):
        # calculate the cell width for an exponential increase
        exp_pad = np.round((cell_width*stretch)*scaling**ii, -2)
        
        # calculate the cell width for a geometric increase by 1.2
        mult_pad = np.round((cell_width*stretch)*((1-stretch**(ii+1))/(1-stretch)), -2)
        
        # take the maximum width for padding
        padding[ii] = max([exp_pad, mult_pad])

    return padding


def get_padding_from_stretch(cell_width, pad_stretch, num_cells):
    """
    get padding cells using pad stretch factor
    
    """
    nodes = np.around(cell_width * (np.ones(num_cells)*pad_stretch)**np.arange(num_cells),-2)
    
    return np.array([nodes[:i].sum() for i in range(1,len(nodes)+1)])
    
    

def get_padding_cells2(cell_width, core_max, max_distance, num_cells):
    """
    get padding cells, which are exponentially increasing to a given 
    distance.  Make sure that each cell is larger than the one previously.
    """
    # check max distance is large enough to accommodate padding
    max_distance = max(cell_width*num_cells, max_distance)

    cells = np.around(np.logspace(np.log10(core_max),np.log10(max_distance),num_cells), -2)
    cells -= core_max
    
    # check if first padding cell is at least as big as cell width
    if cells[1] - cells[0] < cell_width:
        pad_stretch = (core_max + cell_width)/core_max
        cells = get_padding_from_stretch(cell_width, pad_stretch, num_cells)
        print "Provided model extent not wide enough to contain padding, "+\
        "expanding model to {} m".format((cells[-1] + core_max)*2)
        
    return cells