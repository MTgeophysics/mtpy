#!/usr/bin/env python
"""
==================
ModEM
==================

residuals class to contain RMS information

"""

import os.path as op
import numpy as np
from numpy.lib import recfunctions
from mtpy.modeling.modem_data import Data

import mtpy.utils.gocad as mtgocad
reload(mtgocad)
try:
    from evtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')
           
    print ('Note: if you are using Windows you should build evtk first with'
           'either MinGW or cygwin using the command: \n'
           '    python setup.py build -compiler=mingw32  or \n'
           '    python setup.py build -compiler=cygwin')



        


class Residual():
    """
    class to contain residuals for each data point, and rms values for each
    station
    
    ====================== ====================================================
    Attributes/Key Words   Description    
    ====================== ====================================================

    center_position_EN     (east, north, evel) for center point of station 
                           array.  All stations are relative to this location
                           for plotting purposes.
    rms_array              numpy.ndarray structured to store station 
                           location values and rms.  Keys are:
                               * station --> station name
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * zone --> UTM zone
                               * rel_east -- > relative east location to 
                                               center_position (m)
                               * rel_north --> relative north location to 
                                               center_position (m)
                               * rms --> root-mean-square residual for each
                                         station
    residual_array         numpy.ndarray (num_stations) structured to store
                           data.  keys are:
                               * station --> station name
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * rel_east -- > relative east location to 
                                               center_position (m)
                               * rel_north --> relative north location to 
                                               center_position (m)
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * zone --> UTM zone
                               * z --> impedance tensor residual (measured - modelled)
                                       (num_freq, 2, 2)
                               * z_err --> impedance tensor error array with
                                       shape (num_freq, 2, 2)
                               * tip --> Tipper residual (measured - modelled)
                                       (num_freq, 1, 2)
                               * tipperr --> Tipper array with shape
                                       (num_freq, 1, 2)
    residual_fn            full path to data file 
    data_period_list       period list from all the data

    fn_basename            basename of residual file
    header_strings         strings for header of data file following the format
                           outlined in the ModEM documentation
    inv_comp_dict          dictionary of inversion componets
    inv_mode               inversion mode, options are: *default* is '1'
                               * '1' --> for 'Full_Impedance' and 
                                             'Full_Vertical_Components'
                               * '2' --> 'Full_Impedance'
                               * '3' --> 'Off_Diagonal_Impedance' and 
                                         'Full_Vertical_Components'
                               * '4' --> 'Off_Diagonal_Impedance'
                               * '5' --> 'Full_Vertical_Components'
                               * '6' --> 'Full_Interstation_TF'
                               * '7' --> 'Off_Diagonal_Rho_Phase' 

    inv_mode_dict          dictionary for inversion modes
    mt_dict                dictionary of mtpy.core.mt.MT objects with keys 
                           being station names
    units                  [ [V/m]/[T] | [mV/km]/[nT] | Ohm ] units of Z
                           *default* is [mV/km]/[nT]
    wave_sign              [ + | - ] sign of time dependent wave.  
                           *default* is '+' as positive downwards. 
    ====================== ====================================================    
    
    """      

    def __init__(self, **kwargs):
        
        self.workdir = kwargs.pop('workdir','.')
        self.residual_fn = kwargs.pop('residual_fn',None)
        
        
        return
    
    def read_residual_file(self,residual_fn=None):
        
        if residual_fn is not None:
            self.residual_fn = residual_fn
            resObj = Data()
            resObj.read_data_file(self.residual_fn)
        else:
            print "Cannot read residuals, please provide residual_fn"
            return
        
        # pass relevant arguments through residual object
        for att in ['center_position_EN','data_period_list',
                    'wave_sign_impedance','wave_sign_tipper']:
            if hasattr(resObj,att):
                setattr(self,att,getattr(resObj,att))
        
        # define new data types for residual arrays by copying/modifying dtype from data object
        self.residual_array = resObj.data_array.copy()
        
        # append some new fields to contain rms values
        self.rms_array = resObj.station_locations.copy()
        for fieldname in ['rms','rms_z','rms_tip']:
            self.rms_array = recfunctions.append_fields(self.rms_array.copy(),
                                                          fieldname,
                                                          np.zeros(len(resObj.station_locations)),
                                                          usemask=False)

    def calculate_residual_from_data(self,data_fn=None,resp_fn=None):
        
        dataObj = self._read_data_file(data_fn=data_fn)
        respObj = self._read_resp_file(resp_fn=resp_fn)
        
        self.residual_array = dataObj.data_array
        for comp in ['z','tip']:
            self.residual_array[comp] = self.residual_array[comp] - respObj.data_array[comp]
        
        dataObj.fn_basename = respObj.fn_basename[:-3]+'res'
        
        dataObj.write_data_file(fill=False,compute_error=False)
        
        
    def _read_data_file(self,data_fn=None):
        if data_fn is not None:
            self.data_fn = data_fn
            dataObj = Data()
            dataObj.read_data_file(self.data_fn)
        else:
            print "Cannot read data, please provide data_fn"
            return
        
        # pass relevant arguments through residual object
        for att in ['center_position_EN','data_period_list',
                    'wave_sign_impedance','wave_sign_tipper']:
            if hasattr(dataObj,att):
                setattr(self,att,getattr(dataObj,att))        
        
        return dataObj

    def _read_resp_file(self,resp_fn=None):
        if resp_fn is not None:
            self.resp_fn = resp_fn
            respObj = Data()
            respObj.read_data_file(self.resp_fn)
        else:
            print "Cannot read data, please provide data_fn"
            return
        
        # pass relevant arguments through residual object
        for att in ['center_position_EN','data_period_list',
                    'wave_sign_impedance','wave_sign_tipper']:
            if hasattr(respObj,att):
                setattr(self,att,getattr(respObj,att))

        return respObj       
        
    def get_rms(self,residual_fn=None):
        
        if self.residual_array is None:
            self._read_residual_fn()
        if self.residual_array is None:
            return
            
        rms_z_comp = np.zeros((len(self.rms_array),2,2))
        rms_tip_comp = np.zeros((len(self.rms_array),2))
        rms_valuelist_all = np.zeros(0)
        rms_valuelist_z = np.zeros(0)
        rms_valuelist_tip = np.zeros(0)
        
        for stname in self.rms_array['station']:
            rms_valuelist = []
            sta_ind = np.where(self.rms_array['station']==stname)[0][0]
            sta_indd = np.where(self.residual_array['station']==stname)[0][0]
            resvals = self.residual_array[sta_indd]
            znorm,tipnorm = None,None
            if np.amax(np.abs(resvals['z'])) > 0:

                # sum over absolute value of z
                # need to divide by sqrt(2) to normalise (code applies same error to real and imag components)
                znorm = np.abs(resvals['z'])/(np.real(resvals['z_err'])*2.**0.5)
                znorm = znorm[np.all(np.isfinite(znorm),axis=(1,2))]
                
                # append individual normalised errors to a master list for all stations
                rms_valuelist_all = np.append(rms_valuelist_all,znorm.flatten())
                rms_valuelist_z = np.append(rms_valuelist_z,znorm.flatten())
                
                # normalised error for separate components
                rms_z_comp[sta_ind] = (((znorm**2.).sum(axis=0))/(znorm.shape[0]))**0.5
                rms_valuelist.append(rms_z_comp[sta_ind])
                
            if np.amax(np.abs(resvals['tip'])) > 0:
                # sum over absolute value of tipper
                # need to divide by sqrt(2) to normalise (code applies same error to real and imag components)
                tipnorm = np.abs(resvals['tip'])/(np.real(resvals['tip_err'])*2.**0.5)
                tipnorm = tipnorm[np.all(np.isfinite(tipnorm),axis=(1,2))]
                
                # append individual normalised errors to a master list for all stations
                rms_valuelist_all = np.append(rms_valuelist_all,tipnorm.flatten())
                rms_valuelist_tip = np.append(rms_valuelist_tip,tipnorm.flatten())
                
                # normalised error for separate components
                rms_tip_comp[sta_ind] = (((tipnorm**2.).sum(axis=0))/len(tipnorm))**0.5
                rms_valuelist.append(rms_tip_comp[sta_ind])

            rms_valuelist = np.vstack(rms_valuelist).flatten()
            
            rms_value = ((rms_valuelist**2.).sum()/rms_valuelist.size)**0.5

            self.rms_array[sta_ind]['rms'] = rms_value
            
            if znorm is not None:
                self.rms_array[sta_ind]['rms_z'] = ((rms_z_comp[sta_ind]**2.).sum()/rms_z_comp[sta_ind].size)**0.5
            if tipnorm is not None:
                self.rms_array[sta_ind]['rms_tip'] = ((rms_tip_comp[sta_ind]**2.).sum()/rms_z_comp[sta_ind].size)**0.5
            
        self.rms = np.mean(rms_valuelist_all**2.)**0.5
        self.rms_z = np.mean(rms_valuelist_z**2.)**0.5
        self.rms_tip = np.mean(rms_valuelist_tip**2.)**0.5


    def write_rms_to_file(self):
        """
        write rms station data to file
        """
        
        fn = op.join(self.workdir,'rms_values.dat')
        
        if not hasattr(self,'rms'):
            self.get_rms()

        headerlist = ['station','lon','lat','rel_east','rel_north','rms','rms_z','rms_tip']
        
        dtype = []
        for val in headerlist:
            if val == 'station':
                dtype.append((val,'S10'))
            else:
                dtype.append((val,np.float))        
        
        savelist = np.zeros(len(self.rms_array),dtype=dtype)
        for val in headerlist:
            savelist[val] = self.rms_array[val]
        
        header = ' '.join(headerlist)
        
        np.savetxt(fn,savelist,header=header,fmt=['%s','%.6f','%.6f','%.1f','%.1f','%.3f','%.3f','%.3f'])

