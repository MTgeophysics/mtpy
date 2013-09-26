# -*- coding: utf-8 -*-
"""
=============
Static Shift
=============

    * module for estimating static shift
    
Created on Mon Aug 19 10:06:21 2013

@author: jpeacock
"""

#==============================================================================
import mtpy.core.edi as mtedi
import os
import numpy as np
#==============================================================================

def compute_spatial_median_ss(edi_file, ss_tol=0.2, freq_tol=0.15, 
                              distance_radius=1000, n_freq=20, 
                              write_new_edi='y'):
    """
    Compute the median of all stations within a given radius (distance_radius).  
    For each station only resistivities up to  n_freq will be used to
    estimate the median resistivity.  The input station resistivity is then
    divided by the median.  The median is only calculated for the off-diagonal
    terms.  This works well if your station spacing is small <= 500m.
    
    Arguments:
    -----------
        **edi_file**  : string
                        full path to edi file to compute static shift for
                        Note that the other stations are assumed to be in 
                        the same directory as edi_file.
        
        **ss_tol** : float
                     tolerance or error about 1.0 to give to static shift 
                     estimation.  if 1-ss_tol < ss_estimation > 1+ss_tol then
                     the returned value is the ss_estimation, otherwise 1.0 
                     is returned.  *default* is 0.2
                     
        **freq_tol** : float
                       tolerance to look for frequencies in other edi files
                       that were found in edi_file.  This is important if
                       the frequencies for each edi file are different.
                       frequency match occurs when 
                       edi_file(freq)*(1-freq_tol) < edi_file2(freq) and 
                       edi_file2(freq) > edi_file(freq)*(1+freq_tol) 
                       *default* is 0.15
                       
        **distance_radius** : float
                              radius to search for other stations.  Any
                              station that is within this radius is used
                              to calculate the median static shift.
                              *default* is 1000.
                              
        **n_freq** : int
                     number of frequencies to use to calculate median, assuming
                     the first frequency is the highest frequency.
                     *default* is 20.
        
        **write_new_edi** : [ 'y' | 'n' ]
                            * 'y' will write a new edi file if there is a 
                              static shift
                            * 'n' will not write a new edi_file
                            
                            File will be saved to edi_path\SS\station.edi
                     
    Returns:
    --------
        **static_shift_x** : float
                             estimated median static shift in x direction
                             
        **static_shift_y** : float
                             estimated median static shift in y direction
        
        **new_edi_fn** : string
                      full path to new edi file if write_new_edi == 'y'
                      otherwise None is returned.
              
    """
    
    #convert meters to decimal degrees so we don't have to deal with zone 
    #changes
    dm_deg = distance_radius*8.994423457456377e-06
    
    #get path to all the other edi_files
    edi_path = os.path.dirname(edi_file)
    
    #make a list of edi files in the directory
    edi_list = [os.path.join(edi_path, edi) for edi in os.listdir(edi_path)
                if edi.find('.edi') > 0]
   
   #remove the edi_file to be estimated for
    edi_list.remove(edi_file)
    n_edi = len(edi_list)
    
    #read the edi file
    edi1 = mtedi.Edi(filename=edi_file)

    #check to make sure the first frequency is the highest, if not flip 
    #things around    
    if edi1.freq[0] < edi1.freq[-1]:
        edi1.freq = edi1.freq[::-1]
        edi1.Z.z = edi1.Z.z[::-1]
    
    #make a dictionary of frequencies
    freq_dict = dict([('{:.6g}'.format(ff), ii) 
                        for ii, ff in enumerate(edi1.freq[0:n_freq])])                
    
    #get resistivity and phase
    res1, phase1, reserr1, phaseerr1 = edi1.Z.res_phase
    
    #make a resistivity array invovling all stations that are within dm
    res_array = np.zeros((n_edi, n_freq, 2, 2))
    
    kk = 0
    for kedi in edi_list:
        edi2 = mtedi.Edi(kedi)
        
        #check to make sure the first frequency is the highest, if not flip 
        #things around    
        if edi2.freq[0] < edi2.freq[-1]:
            edi2.freq = edi2.freq[::-1]
            edi2.Z.z = edi2.Z.z[::-1]
            
        delta_d = np.sqrt((edi1.lat-edi2.lat)**2+(edi1.lon-edi2.lon)**2)
        if delta_d <= dm_deg:
            print 'Station {0} is within dm'.format(edi2.station)
            res2, phase2, reserr2, phaseerr2 = edi2.Z.res_phase
            for jj, ff in enumerate(edi2.freq[0:n_freq]):
                try:
                    ii = freq_dict['{0:.6g}'.format(ff)]
                    res_array[kk, ii, :, :] = res2[jj, :, :]
                except KeyError:
                    freq_find = False
                    for fkey in freq_dict.keys():
                        if float(fkey)*(1-freq_tol) < ff < float(fkey)*(1+freq_tol):
                            ii = freq_dict[fkey]
                            res_array[kk, ii, :, :] = res2[jj, :, :]
                            freq_find = True
                            break
                        else:
                            pass
                    if freq_find == False:
                        print 'Did not find frequency {0:.6g} for {1}'.format(
                                ff, edi2.station)
            kk += 1 
    
    if kk == 0:
        print '**** No stations with in {0} m'.format(distance_radius)
        return 1.0, 1.0, None

    #convert the res array into a masked array for averaging
    res_array = np.ma.masked_equal(res_array,0)
    
    #compute the static shift of x-components
    
    static_shift_x = res1[0:n_freq, 0, 1]/np.ma.extras.median(res_array[:, :, 0, 1], 
                                                    axis=0)                                
    static_shift_x = np.median(static_shift_x[np.where(static_shift_x!=np.inf)])
    
    #check to see if the estimated static shift is within given tolerance
    if 1.0-ss_tol < static_shift_x < 1.0+ss_tol:
        static_shift_x = 1.0
    
    #compute the static shift of y-components
    static_shift_y = res1[0:n_freq, 1, 0]/np.ma.extras.median(res_array[:, :, 1, 0], 
                                                    axis=0)
    static_shift_y = np.median(static_shift_y[np.where(static_shift_y!=np.inf)])
  
    #check to see if the estimated static shift is within given tolerance
    if 1.0-ss_tol < static_shift_y < 1.0+ss_tol:
        static_shift_y = 1.0
    
    print 'Static shift in x-direction = {0:.2f}'.format(static_shift_x)
    print 'Static shift in y-direction = {0:.2f}'.format(static_shift_y)
    
    #write new edi file if there is a static shift correction
    if write_new_edi == 'y':
        svpath = os.path.join(edi_path, 'SS')
        if not os.path.exists(svpath):
            os.mkdir(svpath)
        
        new_edi_fn = os.path.join(svpath, 
                                  os.path.basename(edi_file)[:-4]+'_ss.edi')
        
        static_shift, new_z = edi1.Z.no_ss(static_shift_x, 
                                           static_shift_y)
        
        edi1.Z.z = new_z

        edi1.writefile(new_edi_fn)
        
        return static_shift_x, static_shift_y, new_edi_fn
    else:
        
        return static_shift_x, static_shift_y, None

                            
            






