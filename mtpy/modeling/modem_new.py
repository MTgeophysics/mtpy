#!/usr/bin/env python
"""
==================
ModEM
==================


# Generate data file for ModEM
# by Paul Soeffky 2013
# revised by LK 2014
# revised by JP 2014

"""

import os
import mtpy.core.edi as mtedi
import mtpy.core.z as mtz
import mtpy.core.mt as mt
import numpy as np
import mtpy.utils.latlongutmconversion as utm2ll
import mtpy.modeling.ws3dinv as ws
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.widgets as widgets
import matplotlib.colors as colors
import matplotlib.cm as cm
import mtpy.modeling.winglink as wl
import mtpy.utils.exceptions as mtex
import mtpy.analysis.pt as mtpt
import mtpy.imaging.mtcolors as mtcl

#==============================================================================

class Data(object):
    """
    Data will read and write .dat files for ModEM and convert a WS data file 
    to ModEM format.  
    
    Arguments:
    ------------
        **edi_list** : list 
                       list of full paths to .edi files you want to invert for
        
    ====================== ====================================================
    Attributes/Key Words   Description    
    ====================== ====================================================
    _dtype                 internal variable defining the data type of 
                           data_array 
    _t_shape               internal variable defining shape of tipper array in
                           _dtype  
    _z_shape               internal variable defining shape of Z array in 
                           _dtype
    center_position        (east, north, evel) for center point of station 
                           array.  All stations are relative to this location
                           for plotting purposes.
    comp_index_dict        dictionary for index values of component of Z and T
    coord_array            numpy.ndarray structured to store station 
                           location values.  Keys are:
                               * station --> station name
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * rel_east -- > relative east location to 
                                               center_position (m)
                               * rel_north --> relative north location to 
                                               center_position (m)
    data_array             numpy.ndarray (num_stations) structured to store
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
                               * z --> impedance tensor array with shape
                                       (num_freq, 2, 2)
                               * z_err --> impedance tensor error array with
                                       shape (num_freq, 2, 2)
                               * tip --> Tipper array with shape
                                       (num_freq, 1, 2)
                               * tipperr --> Tipper array with shape
                                       (num_freq, 1, 2)
    data_fn                full path to data file 
    data_period_list       period list from all the data
    edi_list               list of full paths to edi files
    error_egbert           percentage to multiply sqrt(Z_xy*Zyx) by. 
                           *default* is 3 as prescribed by Egbert & Kelbert
    error_floor            percentage to set the error floor at, anything below
                           this number will be set to error_floor. 
                           *default* is 10
    error_tipper           absolute tipper error, all tipper error will be 
                           set to this value. *default* is .05 for 5%
    error_type             [ 'floor' | 'value' | 'egbert' ]
                           *default* is 'egbert'
                                * 'floor' sets the error floor to error_floor
                                * 'value' sets error to error_value
                                * 'egbert' sets error to  
                                           error_egbert * sqrt(abs(zxy*zyx))
                                           
    error_value            percentage to multiply Z by to set error
                           *default* is 5 for 5% of Z as error
    fn_basename            basename of data file. *default* is 'ModEM_Data.dat'
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
    max_num_periods        maximum number of periods
    mt_dict                dictionary of mtpy.core.mt.MT objects with keys 
                           being station names
    period_dict            dictionary of period index for period_list
    period_list            list of periods to invert for
    period_max             maximum value of period to invert for
    period_min             minimum value of period to invert for           
    save_path              path to save data file to
    units                  [ [V/m]/[T] | [mV/km]/[nT] | Ohm ] units of Z
                           *default* is [mV/km]/[nT]
    wave_sign              [ + | - ] sign of time dependent wave.  
                           *default* is '+' as positive downwards. 
    ====================== ====================================================
   
    ========================== ================================================
    Methods                    Description 
    ========================== ================================================
    convert_ws3dinv_data_file  convert a ws3dinv file to ModEM fomrat, 
                               **Note** this doesn't include tipper data and 
                               you need a station location file like the one
                               output by mtpy.modeling.ws3dinv
    get_data_from_edi          get data from given .edi files and fill 
                               attributes accordingly
    get_mt_dict                get a dictionary of mtpy.core.mt.MT objects 
                               with keys being station names
    get_period_list            get a list of periods to invert for
    get_station_locations      get station locations and relative locations
                               filling in coord_array
    read_data_file             read in a ModEM data file and fill attributes
                               data_array, coord_array, period_list, mt_dict
    write_data_file            write a ModEM data file 
    ========================== ================================================
    
    
    :Example 1 --> create inversion period list: ::
    
        >>> import os
        >>> import mtpy.modeling.modem as modem
        >>> edi_path = r"/home/mt/edi_files"
        >>> edi_list = [os.path.join(edi_path, edi) \
                        for edi in os.listdir(edi_path)\
                        if edi.find('.edi') > 0]
        >>> md = modem.Data(edi_list, period_min=.1, period_max=300,\
                            max_num_periods=12)
        >>> md.write_data_file(save_path=r"/home/modem/inv1")
        
    :Example 2 --> set inverions period list from data: ::
    
        >>> import os
        >>> import mtpy.modeling.modem as modem
        >>> edi_path = r"/home/mt/edi_files"
        >>> edi_list = [os.path.join(edi_path, edi) \
                        for edi in os.listdir(edi_path)\
                        if edi.find('.edi') > 0]
        >>> md = modem.Data(edi_list)
        >>> #get period list from an .edi file
        >>> mt_obj1 = modem.mt.MT(edi_list[0])
        >>> inv_period_list = 1./mt_obj1.Z.freq
        >>> #invert for every third period in inv_period_list
        >>> inv_period_list = inv_period_list[np.arange(0, len(inv_period_list, 3))]
        >>> md.period_list = inv_period_list        
        >>> md.write_data_file(save_path=r"/home/modem/inv1")
                
    :Example 3 --> change error values: ::
        
        >>> import mtpy.modeling.modem as modem
        >>> mdr = modem.Data()
        >>> mdr.read_data_file(r"/home/modem/inv1/ModEM_Data.dat")
        >>> mdr.error_type = 'floor'
        >>> mdr.error_floor = 10
        >>> mdr.error_tipper = .03
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")
        
    :Example 4 --> change inversion type: ::
        
        >>> import mtpy.modeling.modem as modem
        >>> mdr = modem.Data()
        >>> mdr.read_data_file(r"/home/modem/inv1/ModEM_Data.dat")
        >>> mdr.inv_mode = '3'
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")
        
    :Example 5 --> create mesh first then data file: ::
    
        >>> import mtpy.modeling.modem as modem
        >>> import os
        >>> #1) make a list of all .edi files that will be inverted for 
        >>> edi_path = r"/home/EDI_Files"
        >>> edi_list = [os.path.join(edi_path, edi) 
                        for edi in os.listdir(edi_path) 
        >>> ...         if edi.find('.edi') > 0]
        >>> #2) make a grid from the stations themselves with 200m cell spacing
        >>> mmesh = modem.Model(edi_list=edi_list, cell_size_east=200, 
        >>> ...                cell_size_north=200)
        >>> mmesh.make_mesh()
        >>> # check to see if the mesh is what you think it should be
        >>> msmesh.plot_mesh()
        >>> # all is good write the mesh file
        >>> msmesh.write_model_file(save_path=r"/home/modem/Inv1")
        >>> # create data file
        >>> md = modem.Data(edi_list, coord_array=mmesh.station_locations)
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")
    
                  
    """
    
    def __init__(self, edi_list=None, **kwargs):
        self.edi_list = edi_list

        self.error_type = kwargs.pop('error_type', 'egbert')
        self.error_floor = kwargs.pop('error_floor', 5.0)
        self.error_value = kwargs.pop('error_value', 5.0)
        self.error_egbert = kwargs.pop('error_egbert', 3.0)
        self.error_tipper = kwargs.pop('error_tipper', .05)
        
        self.wave_sign = kwargs.pop('wave_sign', '+')
        self.units = kwargs.pop('units', '[mV/km]/[nT]')
        self.inv_mode = kwargs.pop('inv_mode', '1')
        self.period_list = kwargs.pop('period_list', None)
        self.period_step = kwargs.pop('period_step', 1)
        self.period_min = kwargs.pop('period_min', None)
        self.period_max = kwargs.pop('period_max', None)
        self.max_num_periods = kwargs.pop('max_num_periods', None)
        self.period_dict = None
        self.data_period_list = None
        
        self.fn_basename = kwargs.pop('fn_basename', 'ModEM_Data.dat')
        self.save_path = kwargs.pop('save_path', os.getcwd())
        
        self.coord_array = None
        self.center_position = (0.0, 0.0, 0.0)
        self.data_array = None
        self.mt_dict = None
        
        self._z_shape = (1, 2, 2)
        self._t_shape = (1, 1, 2)
        self._dtype = [('station', '|S10'),
                       ('lat', np.float),
                       ('lon', np.float),
                       ('elev', np.float),
                       ('rel_east', np.float), 
                       ('rel_north', np.float), 
                       ('east', np.float),
                       ('north', np.float),
                       ('z', (np.complex, self._z_shape)),
                       ('z_err', (np.complex, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.complex, self._t_shape))]
        
        self.inv_mode_dict = {'1':['Full_Impedance', 'Full_Vertical_Components'],
                              '2':['Full_Impedance'],
                              '3':['Off_Diagonal_Impedance', 
                                   'Full_Vertical_Components'],
                              '4':['Off_Diagonal_Impedance'],
                              '5':['Full_Vertical_Components'],
                              '6':['Full_Interstation_TF'],
                              '7':['Off_Diagonal_Rho_Phase']}
        self.inv_comp_dict = {'Full_Impedance':['zxx', 'zxy', 'zyx', 'zyy'],
                              'Off_Diagonal_Impedance':['zxy', 'zyx'],
                              'Full_Vertical_Components':['tx', 'ty']}
                              
        self.comp_index_dict = {'zxx': (0, 0), 'zxy':(0, 1), 'zyx':(1, 0), 
                                'zyy':(1, 1), 'tx':(0, 0), 'ty':(0, 1)}
        
                              
        self.header_strings = \
        ['# Created using MTpy error {0} of {1:.0f}%\n'.format(self.error_type, self.error_floor), 
        '# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n']

    def _set_dtype(self, z_shape, t_shape):
        """
        reset dtype
        """
        
        self._z_shape = z_shape
        self._t_shape = t_shape
        
        self._dtype = [('station', '|S10'),
                       ('lat', np.float),
                       ('lon', np.float),
                       ('elev', np.float),
                       ('rel_east', np.float), 
                       ('rel_north', np.float), 
                       ('east', np.float),
                       ('north', np.float),
                       ('z', (np.complex, self._z_shape)),
                       ('z_err', (np.complex, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.complex, self._t_shape))]
                       
    def _set_header_string(self):
        """
        reset the header sring for file
        """
        self.header_strings = \
        ['# Created using MTpy error {0} of {1:.0f}%\n'.format(self.error_type, 
                                                             self.error_floor), 
        '# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n']
   
    def get_mt_dict(self):
        """
        get mt_dict from edi file list
        """
        
        if self.edi_list is None:
            raise ModEMError('edi_list is None, please input a list of '
                             '.edi files containing the full path')
                             
        if len(self.edi_list) == 0:
            raise ModEMError('edi_list is empty, please input a list of '
                             '.edi files containing the full path' )
                             
        self.mt_dict = {}
        for edi in self.edi_list:
            mt_obj = mt.MT(edi)
            self.mt_dict[mt_obj.station] = mt_obj
        
    def get_station_locations(self):
        """
        get station locations from edi files
        """
        
        if self.mt_dict is None:
            self.get_mt_dict()

        #--> read in .edi files and get position information as well as center
        #    station. 
        ns = len(self.mt_dict.keys())
        self.coord_array = np.zeros(ns, dtype=[('station','|S10'),
                                               ('east', np.float),
                                               ('north', np.float),
                                               ('lat', np.float),
                                               ('lon', np.float),
                                               ('elev', np.float),
                                               ('rel_east', np.float),
                                               ('rel_north', np.float)])
                                        
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]
            self.coord_array[ii]['station'] = mt_obj.station
            self.coord_array[ii]['lat'] = float(mt_obj.lat)
            self.coord_array[ii]['lon'] = float(mt_obj.lon)
            self.coord_array[ii]['east'] = float(mt_obj.east)
            self.coord_array[ii]['north'] = float(mt_obj.north)
            self.coord_array[ii]['elev'] = float(mt_obj.elev)
            
        #--> get center of the grid
        east_0 = self.coord_array['east'].mean()
        north_0 = self.coord_array['north'].mean()
        
        #set the center of the grid, for now leve elevation at 0, but later
        #add in elevation.  Also should find the closest station to center
        #of the grid.
        self.center_position = (east_0, north_0, 0.0)
        
        self.coord_array['rel_east'] = self.coord_array['east']-east_0
        self.coord_array['rel_north'] = self.coord_array['north']-north_0
        
        #fill in value for relative location in mt_obj
        for cc in self.coord_array:
            self.mt_dict[cc['station']].grid_east = cc['rel_east']
            self.mt_dict[cc['station']].grid_north = cc['rel_north']
        
    def get_period_list(self):
        """
        make a period list to invert for
        
        """
        
        if self.period_list is not None:
            print '-'*50
            print 'Inverting for periods:'
            for per in self.period_list:
                print '     {0:<12.6f}'.format(per)
            print '-'*50
            return
        
        if self.mt_dict is None:
            self.get_mt_dict()
            
        data_period_list = []
        for s_key in sorted(self.mt_dict.keys()):
            mt_obj = self.mt_dict[s_key]
            data_period_list.extend(list(1./mt_obj.Z.freq))
            
        self.data_period_list = np.array(sorted(list(set(data_period_list)), 
                                       reverse=False))
                                       
        if self.period_min is not None:
            if self.period_max is None:
                raise ModEMError('Need to input period_max')
        if self.period_max is not None:
            if self.period_min is None:
                raise ModEMError('Need to input period_min')
        if self.period_min is not None and self.period_max is not None:
            if self.max_num_periods is None:
                raise ModEMError('Need to input number of periods to use')
                
            min_index = np.where(self.data_period_list >= self.period_min)[0][0]
            max_index = np.where(self.data_period_list <= self.period_max)[0][-1]
            
            pmin = np.log10(self.data_period_list[min_index])
            pmax = np.log10(self.data_period_list[max_index])
            self.period_list = np.logspace(pmin, pmax, num=self.max_num_periods)
            
            print '-'*50
            print 'Inverting for periods:'
            for per in self.period_list:
                print '     {0:<12.6f}'.format(per)
            print '-'*50
                
        if self.period_list is None:
            raise ModEMError('Need to input period_min, period_max, '
                             'max_num_periods or a period_list')
        
        
    
    def get_data_from_edi(self):
        """
        get data from edi files and put into an array for easy manipulation 
        later, this will be handy if you want to rewrite the data file from
        an existing file
        
        """

        self.get_station_locations()
        self.get_period_list()
             
        ns = len(self.mt_dict.keys())
        nf = len(self.period_list)
        
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)
                                              
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]
            self.data_array[ii]['station'] = mt_obj.station
            self.data_array[ii]['lat'] = mt_obj.lat
            self.data_array[ii]['lon'] = mt_obj.lon
            self.data_array[ii]['east'] = mt_obj.east
            self.data_array[ii]['north'] = mt_obj.north
            self.data_array[ii]['elev'] = mt_obj.elev
            self.data_array[ii]['rel_east'] = mt_obj.grid_east
            self.data_array[ii]['rel_north'] = mt_obj.grid_north
            
            #make a dictionary for period as keys and values as index within
            #each mt_obj
            p_dict = dict([(np.round(per, 5), kk) for kk, per in 
                            enumerate(1./mt_obj.Z.freq)])
                            
            #search for period in given period list
            for ff, per in enumerate(self.period_list):
                per = np.round(per, 5)
                jj = None
                try:
                    jj = p_dict[per]
                except KeyError:
                    try:
                        jj = np.where(((1./mt_obj.Z.freq)*.95 <= per) & 
                                      ((1./mt_obj.Z.freq) >= per))[0][0]
                    except IndexError:
                        print 'Could not find {0:<12.6f} in {1}'.format(per,
                                                                mt_obj.station)
                if jj is not None:
                    self.data_array[ii]['z'][ff] = mt_obj.Z.z[jj, :, :]
                    self.data_array[ii]['z_err'][ff] = mt_obj.Z.zerr[jj, :, :]

                    if mt_obj.Tipper.tipper is not None:
                        self.data_array[ii]['tip'][ff] = \
                                        mt_obj.Tipper.tipper[jj, :, :]
                        
                        self.data_array[ii]['tip_err'][ff] = \
                                        mt_obj.Tipper.tippererr[jj, :, :]                                            
                    
    def write_data_file(self, save_path=None, fn_basename=None):
        """
        write data file for ModEM
        
        will save file as save_path/fn_basename
        
        Arguments:
        ------------
            **save_path** : string
                            directory path to save data file to.
                            *default* is cwd
                            
            **fn_basename** : string
                              basename to save data file as
                              *default* is 'ModEM_Data.dat'
                              
        Outputs:
        ----------
            **data_fn** : string
                          full path to created data file
                          
        :Example: :: 
                            
            >>> import os
            >>> import mtpy.modeling.modem as modem
            >>> edi_path = r"/home/mt/edi_files"
            >>> edi_list = [os.path.join(edi_path, edi) \
                            for edi in os.listdir(edi_path)\
                            if edi.find('.edi') > 0]
            >>> md = modem.Data(edi_list, period_min=.1, period_max=300,\
                                max_num_periods=12)
            >>> md.write_data_file(save_path=r"/home/modem/inv1")
        """
        
        if save_path is not None:
            self.save_path = save_path
        if fn_basename is not None:
            self.fn_basename = fn_basename
            
        self.data_fn = os.path.join(self.save_path, self.fn_basename)
        
        if self.coord_array is None:
            self.get_station_locations()
        
        if self.period_list is None:
            self.get_period_list()
        
        if self.data_array is None:
            self.get_data_from_edi()
            
        self._set_header_string()

        dlines = []        
        for inv_mode in self.inv_mode_dict[self.inv_mode]:
            dlines.append(self.header_strings[0])
            dlines.append(self.header_strings[1])
            dlines.append('> {0}\n'.format(inv_mode))
            dlines.append('> exp({0}i\omega t)\n'.format(self.wave_sign))
            if inv_mode.find('Impedance') > 0:
                dlines.append('> {0}\n'.format(self.units))
            elif inv_mode.find('Vertical') >=0:
                dlines.append('> []\n')
            dlines.append('> 0\n') #oriention, need to add at some point
            dlines.append('> {0: >7.3f} {1: >7.3f}\n'.format(
                          self.center_position[0], self.center_position[1]))
            dlines.append('> {0} {1}\n'.format(self.data_array['z'].shape[1],
                                               self.data_array['z'].shape[0]))
                                               
            for ss in range(self.data_array['z'].shape[0]):
                for ff in range(self.data_array['z'].shape[1]):
                    for comp in self.inv_comp_dict[inv_mode]:
                        #index values for component with in the matrix
                        z_ii, z_jj = self.comp_index_dict[comp]
                        
                        #get the correct key for data array according to comp
                        if comp.find('z') == 0:
                            c_key = 'z'
                        elif comp.find('t') == 0:
                            c_key = 'tip'
                        
                        #get the value for that compenent at that frequency
                        zz = self.data_array[ss][c_key][ff, z_ii, z_jj]
                        if zz != 0.0+0.0j:
                            per = '{0:<12.5e}'.format(self.period_list[ff])
                            sta = '{0:>7}'.format(self.data_array[ss]['station'])
                            lat = '{0:> 9.3f}'.format(self.data_array[ss]['lat'])
                            lon = '{0:> 9.3f}'.format(self.data_array[ss]['lon'])
                            eas = '{0:> 12.3f}'.format(self.data_array[ss]['rel_east'])
                            nor = '{0:> 12.3f}'.format(self.data_array[ss]['rel_north'])
                            ele = '{0:> 12.3f}'.format(0)
                            com = '{0:>4}'.format(comp.upper())
                            rea = '{0:> 14.6e}'.format(zz.real)
                            ima = '{0:> 14.6e}'.format(zz.imag)
                            
                            #compute relative error
                            if comp.find('t') == 0:
                                rel_err = self.error_tipper
                            elif comp.find('z') == 0:
                                if self.error_type == 'floor':
                                    rel_err = self.data_array[ss][c_key+'_err'][ff, z_ii, z_jj]/\
                                              abs(zz)
                                    if rel_err < self.error_floor/100.:
                                        rel_err = self.error_floor/100.*abs(zz)
                                
                                elif self.error_type == 'value':
                                    rel_err = abs(zz)*self.error_value/100.
                                
                                elif self.error_type == 'egbert':
                                    rel_err = np.sqrt(abs(self.data_array[ss][c_key][ff, 0, 1]*\
                                              self.data_array[ss][c_key][ff, 1, 0]))*\
                                              self.error_egbert/100.
                            
                            rel_err = '{0:> 14.6e}'.format(abs(rel_err))
                            #make sure that x==north, y==east, z==+down
                            dline = ''.join([per, sta, lat, lon, nor, eas, ele, 
                                             com, rea, ima, rel_err, '\n'])
                            dlines.append(dline)
        
        dfid = file(self.data_fn, 'w')
        dfid.writelines(dlines)
        dfid.close()
        
        print 'Wrote ModEM data file to {0}'.format(self.data_fn)
        
    def convert_ws3dinv_data_file(self, ws_data_fn, station_fn=None, 
                                  save_path=None, fn_basename=None):
        """ 
        convert a ws3dinv data file into ModEM format
        
        Arguments:
        ------------
            **ws_data_fn** : string
                             full path to WS data file
            
            **station_fn** : string
                             full path to station info file output by
                             mtpy.modeling.ws3dinv. Or you can create one using
                             mtpy.modeling.ws3dinv.WSStation
            
            **save_path** : string
                            directory path to save data file to.
                            *default* is cwd
                            
            **fn_basename** : string
                              basename to save data file as
                              *default* is 'ModEM_Data.dat'
                              
        Outputs:
        -----------
            **data_fn** : string
                          full path to created data file
        
        :Example: ::
            
            >>> import mtpy.modeling.modem as modem
            >>> mdr = modem.Data()
            >>> mdr.convert_ws3dinv_data_file(r"/home/ws3dinv/inv1/WSData.dat",
                    station_fn=r"/home/ws3dinv/inv1/WS_Station_Locations.txt")
        """
        
        if os.path.isfile(ws_data_fn) == False:
            raise ws.WSInputError('Did not find {0}, check path'.format(ws_data_fn))
        
        if save_path is not None:
            self.save_path = save_path
        else:
            self.save_path = os.path.dirname(ws_data_fn)
    
        if fn_basename is not None:
            self.fn_basename = fn_basename
            
        #--> get data from data file
        wsd = ws.WSData()
        wsd.read_data_file(ws_data_fn, station_fn=station_fn)
        
        ns = wsd.data['station'].shape[0]
        nf = wsd.period_list.shape[0]
        
        self.period_list = wsd.period_list.copy()
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)
        
        #--> fill coord array
        self.coord_array = np.zeros(ns, dtype=[('station','|S10'),
                                               ('east', np.float),
                                               ('north', np.float),
                                               ('lat', np.float),
                                               ('lon', np.float),
                                               ('elev', np.float),
                                               ('rel_east', np.float),
                                               ('rel_north', np.float)])
                                                    
        #--> fill data array
        for ii, d_arr in enumerate(wsd.data):
            self.data_array[ii]['station'] = d_arr['station']
            self.data_array[ii]['rel_east'] = d_arr['east']
            self.data_array[ii]['rel_north'] = d_arr['north']
            self.data_array[ii]['z'][:] = d_arr['z_data']
            self.data_array[ii]['z_err'][:] = d_arr['z_data_err'].real*\
                                                d_arr['z_err_map'].real
            self.coord_array[ii]['station'] = d_arr['station']
            self.coord_array[ii]['lat'] = 0.0
            self.coord_array[ii]['lon'] = 0.0
            self.coord_array[ii]['rel_east'] = d_arr['east']
            self.coord_array[ii]['rel_north'] = d_arr['north']
            self.coord_array[ii]['elev'] = 0.0
            self.data_array[ii]['station'] = d_arr['station']
            self.data_array[ii]['lat'] = 0.0
            self.data_array[ii]['lon'] = 0.0
            self.data_array[ii]['rel_east'] = d_arr['east']
            self.data_array[ii]['rel_north'] = d_arr['north']
            self.data_array[ii]['elev'] = 0.0
      
        #need to change the inversion mode to be the same as the ws_data file
        if self.data_array['z'].all() == 0.0:
            if self.data_array['tip'].all() == 0.0:
                self.inv_mode = '4'
            else:
                self.inv_mode = '3'
        else:
            if self.data_array['tip'].all() == 0.0:
                self.inv_mode = '2'
            else:
                self.inv_mode = '1'
        
        #-->write file
        self.write_data_file()
        
    def read_data_file(self, data_fn=None):
        """
        read ModEM data file
        
        Fills attributes: 
            * data_array
            * period_list
            * mt_dict
        
        """
        
        if data_fn is not None:
            self.data_fn = data_fn
            
        if self.data_fn is None:
            raise ModEMError('data_fn is None, enter a data file to read.')
        elif os.path.isfile(self.data_fn) is False:
            raise ModEMError('Could not find {0}, check path'.format(self.data_fn))
            
        dfid = file(self.data_fn, 'r')
        dlines = dfid.readlines()
        dfid.close()
        
        header_list = []
        metadata_list = []
        data_list = []
        period_list = []
        station_list = []
        for dline in dlines:
            if dline.find('#') == 0:
                header_list.append(dline.strip())
            elif dline.find('>') == 0:
                metadata_list.append(dline[1:].strip())
            else:
                dline_list = dline.strip().split()
                if len(dline_list) == 11:
                    for ii, d_str in enumerate(dline_list):
                        try:
                            dline_list[ii] = float(d_str.strip())
                        except ValueError:
                            pass
                    period_list.append(dline_list[0])
                    station_list.append(dline_list[1])
                    
                    data_list.append(dline_list)
            
        self.period_list = np.array(sorted(set(period_list)))
        station_list = sorted(set(station_list))
        
        #make a period dictionary to with key as period and value as index
        period_dict = dict([(per, ii) for ii, per in enumerate(self.period_list)])
        
        #--> need to sort the data into a useful fashion such that each station
        #    is an mt object
        
        data_dict = {}
        z_dummy = np.zeros((len(self.period_list), 2, 2), dtype='complex')
        t_dummy = np.zeros((len(self.period_list), 1, 2), dtype='complex')
        
        index_dict = {'zxx': (0, 0), 'zxy':(0, 1), 'zyx':(1, 0), 'zyy':(1, 1),
                      'tx':(0, 0), 'ty':(0, 1)} 
        
        #dictionary for true false if station data (lat, lon, elev, etc) 
        #has been filled already so we don't rewrite it each time              
        tf_dict = {}
        for station in station_list:
            data_dict[station] = mt.MT()
            data_dict[station].Z = mtz.Z(z_array=z_dummy.copy(), 
                                        zerr_array=z_dummy.copy().real, 
                                        freq=1./self.period_list)
            data_dict[station].Tipper = mtz.Tipper(tipper_array=t_dummy.copy(), 
                                                   tippererr_array=t_dummy.copy().real, 
                                                   freq=1./self.period_list)
            #make sure that the station data starts out with false to fill
            #the data later
            tf_dict[station] = False
            
        #fill in the data for each station                                      
        for dd in data_list:
            #get the period index from the data line
            p_index = period_dict[dd[0]]
            #get the component index from the data line
            ii, jj = index_dict[dd[7].lower()]
            
            #if the station data has not been filled yet, fill it
            if tf_dict[dd[1]] == False:
                data_dict[dd[1]].lat = dd[2]
                data_dict[dd[1]].lon = dd[3]
                data_dict[dd[1]].grid_north = dd[4]
                data_dict[dd[1]].grid_east = dd[5]
                data_dict[dd[1]].grid_elev = dd[6]
                data_dict[dd[1]].station = dd[1]
                tf_dict[dd[1]] = True
            #fill in the impedance tensor with appropriate values
            if dd[7].find('Z') == 0:
                data_dict[dd[1]].Z.z[p_index, ii, jj] = dd[8]+1j*dd[9]
                data_dict[dd[1]].Z.zerr[p_index, ii, jj] = dd[10]
            #fill in tipper with appropriate values
            elif dd[7].find('T') == 0:
                data_dict[dd[1]].Tipper.tipper[p_index, ii, jj] = dd[8]+1j*dd[9]
                data_dict[dd[1]].Tipper.tippererr[p_index, ii, jj] = dd[10]
       
        #make mt_dict an attribute for easier manipulation later
        self.mt_dict = data_dict
        
        #--> fill coord array
        ns = len(self.mt_dict.keys())
        nf = len(self.period_list)
        self.coord_array = np.zeros(ns, dtype=[('station','|S10'),
                                               ('east', np.float),
                                               ('north', np.float),
                                               ('lat', np.float),
                                               ('lon', np.float),
                                               ('elev', np.float),
                                               ('rel_east', np.float),
                                               ('rel_north', np.float)])
        
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)

        #Be sure to caclulate invariants and phase tensor for each station
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]
            
            self.mt_dict[s_key].zinv.compute_invariants()
            self.mt_dict[s_key].pt.set_z_object(mt_obj.Z)
            self.coord_array[ii]['station'] = mt_obj.station
            self.coord_array[ii]['lat'] = mt_obj.lat
            self.coord_array[ii]['lon'] = mt_obj.lon
            self.coord_array[ii]['elev'] = mt_obj.elev
            self.coord_array[ii]['rel_east'] = mt_obj.grid_east
            self.coord_array[ii]['rel_north'] = mt_obj.grid_north
            
            
            self.data_array[ii]['station'] = mt_obj.station
            self.data_array[ii]['lat'] = mt_obj.lat
            self.data_array[ii]['lon'] = mt_obj.lon
            self.data_array[ii]['east'] = mt_obj.east
            self.data_array[ii]['north'] = mt_obj.north
            self.data_array[ii]['elev'] = mt_obj.elev
            self.data_array[ii]['rel_east'] = mt_obj.grid_east
            self.data_array[ii]['rel_north'] = mt_obj.grid_north
            
            
            self.data_array[ii]['z'][:] = mt_obj.Z.z
            self.data_array[ii]['z_err'][:] = mt_obj.Z.zerr

            self.data_array[ii]['tip'][:] =  mt_obj.Tipper.tipper
            self.data_array[ii]['tip_err'][:] = mt_obj.Tipper.tippererr 
            
        
#==============================================================================
# mesh class
#==============================================================================
class Model(object):
    """
    make and read a FE mesh grid
    
    The mesh assumes the coordinate system where:
        x == North
        y == East
        z == + down
        
    All dimensions are in meters.
    
    :Example 1 --> create mesh first then data file: ::
    
        >>> import mtpy.modeling.modem as modem
        >>> import os
        >>> #1) make a list of all .edi files that will be inverted for 
        >>> edi_path = r"/home/EDI_Files"
        >>> edi_list = [os.path.join(edi_path, edi) 
                        for edi in os.listdir(edi_path) 
        >>> ...         if edi.find('.edi') > 0]
        >>> #2) make a grid from the stations themselves with 200m cell spacing
        >>> mmesh = modem.Model(edi_list=edi_list, cell_size_east=200, 
        >>> ...                cell_size_north=200)
        >>> mmesh.make_mesh()
        >>> # check to see if the mesh is what you think it should be
        >>> msmesh.plot_mesh()
        >>> # all is good write the mesh file
        >>> msmesh.write_model_file(save_path=r"/home/modem/Inv1")
        >>> # create data file
        >>> md = modem.Data(edi_list, coord_array=mmesh.station_locations)
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")
    
    :Example 2 --> create data file first then model file: ::
    
        >>> import mtpy.modeling.modem as modem
        >>> import os
        >>> #1) make a list of all .edi files that will be inverted for 
        >>> edi_path = r"/home/EDI_Files"
        >>> edi_list = [os.path.join(edi_path, edi) 
                        for edi in os.listdir(edi_path) 
        >>> ...         if edi.find('.edi') > 0]
        >>> #2) create data file
        >>> md = modem.Data(edi_list, coord_array=mmesh.station_locations)
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")
        >>> #3) make a grid from the stations themselves with 200m cell spacing
        >>> mmesh = modem.Model(edi_list=edi_list, cell_size_east=200, 
                                cell_size_north=200, 
                                station_locations=md.coord_array)
        >>> mmesh.make_mesh()
        >>> # check to see if the mesh is what you think it should be
        >>> msmesh.plot_mesh()
        >>> # all is good write the mesh file
        >>> msmesh.write_model_file(save_path=r"/home/modem/Inv1")
        
        
    
    ==================== ======================================================
    Attributes           Description    
    ==================== ======================================================
    cell_size_east       mesh block width in east direction
                         *default* is 500
    cell_size_north      mesh block width in north direction
                         *default* is 500
    edi_list             list of .edi files to invert for
    grid_east            overall distance of grid nodes in east direction 
    grid_north           overall distance of grid nodes in north direction 
    grid_z               overall distance of grid nodes in z direction 
    model_fn           full path to initial file name
    n_layers             total number of vertical layers in model
    nodes_east           relative distance between nodes in east direction 
    nodes_north          relative distance between nodes in north direction 
    nodes_z              relative distance between nodes in east direction 
    pad_east             number of cells for padding on E and W sides
                         *default* is 5
    pad_north            number of cells for padding on S and N sides
                         *default* is 5
    pad_root_east        padding cells E & W will be pad_root_east**(x)
    pad_root_north       padding cells N & S will be pad_root_north**(x) 
    pad_z                number of cells for padding at bottom
                         *default* is 5
    res_list             list of resistivity values for starting model
    res_model            starting resistivity model
    rotation_angle       Angle to rotate the grid to. Angle is measured
                         positve clockwise assuming North is 0 and east is 90.
                         *default* is None
    save_path            path to save file to  
    station_fn           full path to station file
    station_locations    location of stations
    title                title in initial file
    z1_layer             first layer thickness
    z_bottom             absolute bottom of the model *default* is 300,000 
    z_target_depth       Depth of deepest target, *default* is 50,000
    ==================== ======================================================
    
    ==================== ======================================================
    Methods              Description
    ==================== ======================================================
    make_mesh            makes a mesh from the given specifications
    plot_mesh            plots mesh to make sure everything is good
    write_initial_file   writes an initial model file that includes the mesh
    ==================== ======================================================
    
    
    """
    
    def __init__(self, edi_list=None, **kwargs):
        
        self.edi_list = edi_list
        
        # size of cells within station area in meters
        self.cell_size_east = kwargs.pop('cell_size_east', 500)
        self.cell_size_north = kwargs.pop('cell_size_north', 500)
        
        #padding cells on either side
        self.pad_east = kwargs.pop('pad_east', 7)
        self.pad_north = kwargs.pop('pad_north', 7)
        self.pad_z = kwargs.pop('pad_z', 4)
        
        #root of padding cells
        self.pad_stretch_h= kwargs.pop('pad_stretch_h', 1.2)
        self.pad_stretch_v= kwargs.pop('pad_stretch_v', 1.2)
        
        self.z1_layer = kwargs.pop('z1_layer', 10)
        self.z_target_depth = kwargs.pop('z_target_depth', 50000)
        self.z_bottom = kwargs.pop('z_bottom', 300000)
        
        #number of vertical layers
        self.n_layers = kwargs.pop('n_layers', 30)
        
        #strike angle to rotate grid to
        self.strike_angle = kwargs.pop('strike_angle', None)
        
        #--> attributes to be calculated
        #station information
        self.station_locations = kwargs.pop('station_locations', None)
       
       #grid nodes
        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None
        
        #grid locations
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None
        
        #resistivity model
        self.res_model = None
        
        self.grid_center = None
        
        #inital file stuff
        self.model_fn = kwargs.pop('model_fn', None)
        self.save_path = kwargs.pop('save_path', None)
        self.model_fn_basename = kwargs.pop('model_fn_basename', 
                                            'ModEM_Model.ws')
        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
            self.model_fn_basename = os.path.basename(self.model_fn)
        
        self.title = 'Model File written by MTpy.modeling.modem'
        self.res_scale = kwargs.pop('res_scale', 'loge')
        
        
    def make_mesh(self):
        """ 
        create finite element mesh according to parameters set.
        
        The mesh is built by first finding the center of the station area.  
        Then cells are added in the north and east direction with width
        cell_size_east and cell_size_north to the extremeties of the station 
        area.  Padding cells are then added to extend the model to reduce 
        edge effects.  The number of cells are pad_east and pad_north and the
        increase in size is by pad_root_east and pad_root_north.  The station
        locations are then computed as the center of the nearest cell as 
        required by the code.
        
        The vertical cells are built to increase in size exponentially with
        depth.  The first cell depth is first_layer_thickness and should be
        about 1/10th the shortest skin depth.  The layers then increase
        on a log scale to z_target_depth.  Then the model is
        padded with pad_z number of cells to extend the depth of the model.
        
        padding = np.round(cell_size_east*pad_root_east**np.arange(start=.5,
                           stop=3, step=3./pad_east))+west 
        
        """
        #if station locations are not input read from the edi files
        if self.station_locations is None:
            if self.edi_list is None:
                raise AttributeError('edi_list is None, need to input a list of '
                                     'edi files to read in.')
                                     
            n_stations = len(self.edi_list)
            
            if n_stations == 0:
                raise ModEMError('No .edi files in edi_list, please check '
                                 'file locations.')
            
            #make a structured array to put station location information into
            self.station_locations = np.zeros(n_stations,
                                              dtype=[('station','|S10'),
                                                     ('lat', np.float),
                                                     ('lon', np.float),
                                                     ('east', np.float),
                                                     ('north', np.float), 
                                                     ('rel_east', np.float),
                                                     ('rel_north', np.float),
                                                     ('elev', np.float)])
            #get station locations in meters
            for ii, edi in enumerate(self.edi_list):
                zz = mtedi.Edi()
                zz.readfile(edi)
                zone, east, north = utm2ll.LLtoUTM(23, zz.lat, zz.lon)
                self.station_locations[ii]['lat'] = zz.lat
                self.station_locations[ii]['lon'] = zz.lon
                self.station_locations[ii]['station'] = zz.station
                self.station_locations[ii]['east'] = east
                self.station_locations[ii]['north'] = north
                self.station_locations[ii]['elev'] = zz.elev
             
            #remove the average distance to get coordinates in a relative space
            self.station_locations['rel_east'] = self.station_locations['east']-\
                                                 self.station_locations['east'].mean()
            self.station_locations['rel_north'] = self.station_locations['north']-\
                                                  self.station_locations['north'].mean()
         
            #translate the stations so they are relative to 0,0
            east_center = (self.station_locations['rel_east'].max()-
                            np.abs(self.station_locations['rel_east'].min()))/2
            north_center = (self.station_locations['rel_north'].max()-
                            np.abs(self.station_locations['rel_north'].min()))/2
            
            #remove the average distance to get coordinates in a relative space
            self.station_locations['rel_east'] -= east_center
            self.station_locations['rel_north'] -= north_center
        
        #pickout the furtherst south and west locations 
        #and put that station as the bottom left corner of the main grid
        west = self.station_locations['rel_east'].min()-self.cell_size_east/2
        east = self.station_locations['rel_east'].max()+self.cell_size_east/2
        south = self.station_locations['rel_north'].min()-self.cell_size_north/2
        north = self.station_locations['rel_north'].max()+self.cell_size_north/2
        west = np.round(west, -2)
        east= np.round(east, -2)
        south= np.round(south, -2)
        north = np.round(north, -2)


        #make sure the variable n_stations is initialized        
        try:
            n_stations
        except NameError:
            n_stations = self.station_locations.shape[0]

        #-------make a grid around the stations from the parameters above------
        #--> make grid in east-west direction
        #cells within station area
        east_gridr = np.arange(start=west, stop=east+self.cell_size_east,
                               step=self.cell_size_east)
        #padding cells in the east-west direction
        
        for ii in range(1, self.pad_east+1):
            east_0 = float(east_gridr[-1])
            west_0 = float(east_gridr[0])
            add_size = np.round(self.cell_size_east*self.pad_stretch_h*ii, -2)
            pad_w = west_0-add_size
            pad_e = east_0+add_size
            east_gridr = np.insert(east_gridr, 0, pad_w)
            east_gridr = np.append(east_gridr, pad_e)
            
            
        #--> need to make sure none of the stations lie on the nodes
        for s_east in sorted(self.station_locations['rel_east']):
            try:
                node_index = np.where(abs(s_east-east_gridr) < 
                                     .02*self.cell_size_east)[0][0]
                if s_east-east_gridr[node_index] > 0:
                    east_gridr[node_index] -= .02*self.cell_size_east
                elif s_east-east_gridr[node_index] < 0:
                    east_gridr[node_index] += .02*self.cell_size_east
            except IndexError:
                continue
            
        
        #--> make grid in north-south direction 
        #N-S cells with in station area
        north_gridr = np.arange(start=south, stop=north+self.cell_size_north, 
                                step=self.cell_size_north)
        
        #padding cells in the east-west direction
        for ii in range(1, self.pad_north+1):
            south_0 = float(north_gridr[0]) 
            north_0 = float(north_gridr[-1])
            add_size = np.round(self.cell_size_north*self.pad_stretch_h*ii, -2)
            pad_s = south_0-add_size
            pad_n = north_0+add_size
            north_gridr = np.insert(north_gridr, 0, pad_s)
            north_gridr = np.append(north_gridr, pad_n)
            
        #--> need to make sure none of the stations lie on the nodes
        for s_north in sorted(self.station_locations['rel_north']):
            try:
                node_index = np.where(abs(s_north-north_gridr) < 
                                     .02*self.cell_size_north)[0][0]
                if s_north-north_gridr[node_index] > 0:
                    north_gridr[node_index] -= .02*self.cell_size_north
                elif s_north-north_gridr[node_index] < 0:
                    north_gridr[node_index] += .02*self.cell_size_north
            except IndexError:
                continue
            
        #--> make depth grid
        log_z = np.logspace(np.log10(self.z1_layer), 
                            np.log10(self.z_target_depth-np.logspace(np.log10(self.z1_layer), 
                            np.log10(self.z_target_depth), 
                            num=self.n_layers)[-2]), 
                            num=self.n_layers-self.pad_z)
        z_nodes = np.array([zz-zz%10**np.floor(np.log10(zz)) for zz in 
                           log_z])
        #padding cells in the east-west direction
        for ii in range(1, self.pad_z+1):
            z_0 = np.float(z_nodes[-2])
            pad_d = np.round(z_0*self.pad_stretch_v*ii, -2)
            z_nodes = np.append(z_nodes, pad_d)                  
        
        #make an array of absolute values
        z_grid = np.array([z_nodes[:ii+1].sum() for ii in range(z_nodes.shape[0])])
        
        #---Need to make an array of the individual cell dimensions for
        #   modem
        east_nodes = east_gridr.copy()    
        nx = east_gridr.shape[0]
        east_nodes[:nx/2] = np.array([abs(east_gridr[ii]-east_gridr[ii+1]) 
                                          for ii in range(int(nx/2))])
        east_nodes[nx/2:] = np.array([abs(east_gridr[ii]-east_gridr[ii+1]) 
                                          for ii in range(int(nx/2)-1, nx-1)])
    
        north_nodes = north_gridr.copy()
        ny = north_gridr.shape[0]
        north_nodes[:ny/2] = np.array([abs(north_gridr[ii]-north_gridr[ii+1]) 
                                       for ii in range(int(ny/2))])
        north_nodes[ny/2:] = np.array([abs(north_gridr[ii]-north_gridr[ii+1]) 
                                       for ii in range(int(ny/2)-1, ny-1)])
                                
        #--put the grids into coordinates relative to the center of the grid
        east_grid = east_nodes.copy()
        east_grid[:int(nx/2)] = -np.array([east_nodes[ii:int(nx/2)].sum() 
                                           for ii in range(int(nx/2))])
        east_grid[int(nx/2):] = np.array([east_nodes[int(nx/2):ii+1].sum() 
                                         for ii in range(int(nx/2), nx)])-\
                                         east_nodes[int(nx/2)]
                                
        north_grid = north_nodes.copy()
        north_grid[:int(ny/2)] = -np.array([north_nodes[ii:int(ny/2)].sum() 
                                            for ii in range(int(ny/2))])
        north_grid[int(ny/2):] = np.array([north_nodes[int(ny/2):ii+1].sum() 
                                            for ii in range(int(ny/2),ny)])-\
                                            north_nodes[int(ny/2)]
        
        #compute grid center
        center_east = -east_nodes.__abs__().sum()/2
        center_north = -north_nodes.__abs__().sum()/2
        center_z = 0
        self.grid_center = np.array([center_north, center_east, center_z])
        
        #make nodes attributes
        self.nodes_east = east_nodes
        self.nodes_north = north_nodes
        self.nodes_z = z_nodes        
        self.grid_east = east_grid
        self.grid_north = north_grid
        self.grid_z = z_grid
            
        #--> print out useful information                    
        print '-'*15
        print '   Number of stations = {0}'.format(len(self.station_locations))
        print '   Dimensions: '
        print '      e-w = {0}'.format(east_grid.shape[0])
        print '      n-s = {0}'.format(north_grid.shape[0])
        print '       z  = {0} (without 7 air layers)'.format(z_grid.shape[0])
        print '   Extensions: '
        print '      e-w = {0:.1f} (m)'.format(east_nodes.__abs__().sum())
        print '      n-s = {0:.1f} (m)'.format(north_nodes.__abs__().sum())
        print '      0-z = {0:.1f} (m)'.format(self.nodes_z.__abs__().sum())
        print '-'*15

    def plot_mesh(self, east_limits=None, north_limits=None, z_limits=None,
                  **kwargs):
        """
        
        Arguments:
        ----------
            **east_limits** : tuple (xmin,xmax)
                             plot min and max distances in meters for the 
                             E-W direction.  If None, the east_limits
                             will be set to furthest stations east and west.
                             *default* is None
                        
            **north_limits** : tuple (ymin,ymax)
                             plot min and max distances in meters for the 
                             N-S direction.  If None, the north_limits
                             will be set to furthest stations north and south.
                             *default* is None
                        
            **z_limits** : tuple (zmin,zmax)
                            plot min and max distances in meters for the 
                            vertical direction.  If None, the z_limits is
                            set to the number of layers.  Z is positive down
                            *default* is None
        """
        
        fig_size = kwargs.pop('fig_size', [6, 6])
        fig_dpi = kwargs.pop('fig_dpi', 300)
        fig_num = kwargs.pop('fig_num', 1)
        
        station_marker = kwargs.pop('station_marker', 'v')
        marker_color = kwargs.pop('station_color', 'b')
        marker_size = kwargs.pop('marker_size', 2)
        
        line_color = kwargs.pop('line_color', 'k')
        line_width = kwargs.pop('line_width', .5)
        
        plt.rcParams['figure.subplot.hspace'] = .3
        plt.rcParams['figure.subplot.wspace'] = .3
        plt.rcParams['figure.subplot.left'] = .08
        plt.rcParams['font.size'] = 7
        
        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        plt.clf()
        
        #---plot map view    
        ax1 = fig.add_subplot(1, 2, 1, aspect='equal')
        
        #make sure the station is in the center of the cell
        ax1.scatter(self.station_locations['rel_east'],
                    self.station_locations['rel_north'], 
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)
                
        #plot the grid if desired
        east_line_xlist = []
        east_line_ylist = []            
        for xx in self.grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend([self.grid_north.min(), 
                                    self.grid_north.max()])
            east_line_ylist.append(None)
        ax1.plot(east_line_xlist,
                      east_line_ylist,
                      lw=line_width,
                      color=line_color)

        north_line_xlist = []
        north_line_ylist = [] 
        for yy in self.grid_north:
            north_line_xlist.extend([self.grid_east.min(),
                                     self.grid_east.max()])
            north_line_xlist.append(None)
            north_line_ylist.extend([yy, yy])
            north_line_ylist.append(None)
        ax1.plot(north_line_xlist,
                      north_line_ylist,
                      lw=line_width,
                      color=line_color)
        
        if east_limits == None:
            ax1.set_xlim(self.station_locations['rel_east'].min()-\
                            10*self.cell_size_east,
                         self.station_locations['rel_east'].max()+\
                             10*self.cell_size_east)
        else:
            ax1.set_xlim(east_limits)
        
        if north_limits == None:
            ax1.set_ylim(self.station_locations['rel_north'].min()-\
                            10*self.cell_size_north,
                         self.station_locations['rel_north'].max()+\
                             10*self.cell_size_east)
        else:
            ax1.set_ylim(north_limits)
            
        ax1.set_ylabel('Northing (m)', fontdict={'size':9,'weight':'bold'})
        ax1.set_xlabel('Easting (m)', fontdict={'size':9,'weight':'bold'})
        
        ##----plot depth view
        ax2 = fig.add_subplot(1, 2, 2, aspect='auto', sharex=ax1)
        

        #plot the grid if desired
        east_line_xlist = []
        east_line_ylist = []            
        for xx in self.grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend([0, 
                                    self.grid_z.max()])
            east_line_ylist.append(None)
        ax2.plot(east_line_xlist,
                 east_line_ylist,
                 lw=line_width,
                 color=line_color)

        z_line_xlist = []
        z_line_ylist = [] 
        for zz in self.grid_z:
            z_line_xlist.extend([self.grid_east.min(),
                                     self.grid_east.max()])
            z_line_xlist.append(None)
            z_line_ylist.extend([zz, zz])
            z_line_ylist.append(None)
        ax2.plot(z_line_xlist,
                 z_line_ylist,
                 lw=line_width,
                 color=line_color)
                      
        
        #--> plot stations
        ax2.scatter(self.station_locations['rel_east'],
                    [0]*self.station_locations.shape[0],
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        
        if z_limits == None:
            ax2.set_ylim(self.z_target_depth, -200)
        else:
            ax2.set_ylim(z_limits)
            
        if east_limits == None:
            ax1.set_xlim(self.station_locations['rel_east'].min()-\
                            10*self.cell_size_east,
                         self.station_locations['rel_east'].max()+\
                             10*self.cell_size_east)
        else:
            ax1.set_xlim(east_limits)
            
        ax2.set_ylabel('Depth (m)', fontdict={'size':9, 'weight':'bold'})
        ax2.set_xlabel('Easting (m)', fontdict={'size':9, 'weight':'bold'})  
        
        plt.show()
        
    
    def write_model_file(self, **kwargs):
        """
        will write an initial file for ModEM.  
        
        Note that x is assumed to be S --> N, y is assumed to be W --> E and
        z is positive downwards.  This means that index [0, 0, 0] is the 
        southwest corner of the first layer.  Therefore if you build a model
        by hand the layer block will look as it should in map view. 
        
        Also, the xgrid, ygrid and zgrid are assumed to be the relative 
        distance between neighboring nodes.  This is needed because wsinv3d 
        builds the  model from the bottom SW corner assuming the cell width
        from the init file.
        
           
        
        Key Word Arguments:
        ----------------------
        
            **nodes_north** : np.array(nx)
                        block dimensions (m) in the N-S direction. 
                        **Note** that the code reads the grid assuming that
                        index=0 is the southern most point.
            
            **nodes_east** : np.array(ny)
                        block dimensions (m) in the E-W direction.  
                        **Note** that the code reads in the grid assuming that
                        index=0 is the western most point.
                        
            **nodes_z** : np.array(nz)
                        block dimensions (m) in the vertical direction.  
                        This is positive downwards.
                        
            **save_path** : string
                          Path to where the initial file will be saved
                          to savepath/model_fn_basename
                          
            **model_fn_basename** : string
                                    basename to save file to
                                    *default* is ModEM_Model.ws
                                    file is saved at savepath/model_fn_basename

            **title** : string
                        Title that goes into the first line 
                        *default* is Model File written by MTpy.modeling.modem 
                        
            **res_model** : np.array((nx,ny,nz))
                        Starting resistivity model. 
                        
                        .. note:: again that the modeling code 
                        assumes that the first row it reads in is the southern
                        most row and the first column it reads in is the 
                        western most column.  Similarly, the first plane it 
                        reads in is the Earth's surface.
                        
            **res_scale** : [ 'loge' | 'log' | 'log10' | 'linear' ]
                            scale of resistivity.  In the ModEM code it 
                            converts everything to Loge, 
                            *default* is 'loge'

        """
        
        keys = ['nodes_east', 'nodes_north', 'nodes_z', 'title',
                'res_model', 'save_path', 'model_fn', 'model_fn_basename']
        for key in keys:
            try:
                setattr(self, key, kwargs[key])
            except KeyError:
                if self.__dict__[key] is None:
                    pass
        
        if self.save_path is not None:
            self.model_fn = os.path.join(self.save_path, 
                                         self.model_fn_basename)
        if self.model_fn is None:
            if self.save_path is None:
                self.save_path = os.getcwd()
                self.model_fn = os.path.join(self.save_path, 
                                             self.model_fn_basename)
            elif os.path.isdir(self.save_path) == True:
                self.model_fn = os.path.join(self.save_path, 
                                             self.model_fn_basename)
            else:
                self.save_path = os.path.dirname(self.save_path)
                self.model_fn= self.save_path
                
        if self.res_model is None or type(self.res_model) is float or\
            type(self.res_model) is int:
            res_model = np.zeros((self.grid_north.shape[0],
                                  self.grid_east.shape[0],
                                  self.grid_z.shape[0]))
                                          
            if self.res_model is None:
                res_model[:, :, :] = 100.0
                self.res_model = res_model
            else:
                res_model[:, :, :] = self.res_model
                self.res_model = res_model
        

        #--> write file
        ifid = file(self.model_fn, 'w')
        ifid.write('# {0}\n'.format(self.title.upper()))
        ifid.write('{0:>5}{1:>5}{2:>5}{3:>5} {4}\n'.format(self.nodes_north.shape[0],
                                              self.nodes_east.shape[0],
                                              self.nodes_z.shape[0],
                                              0,  
                                              self.res_scale.upper()))
    
        #write S --> N node block
        for ii, nnode in enumerate(self.nodes_north):
            ifid.write('{0:>12.3f}'.format(abs(nnode)))

        ifid.write('\n')
        
        #write W --> E node block        
        for jj, enode in enumerate(self.nodes_east):
            ifid.write('{0:>12.3f}'.format(abs(enode)))
        ifid.write('\n')

    
        #write top --> bottom node block
        for kk, zz in enumerate(self.nodes_z):
            ifid.write('{0:>12.3f}'.format(abs(zz)))
        ifid.write('\n')
    
        #write the resistivity in log e format
        if self.res_scale.lower() == 'loge':
            write_res_model = np.log(self.res_model[::-1, :, :])
        elif self.res_scale.lower() == 'log' or \
             self.res_scale.olower() == 'log10':
            write_res_model = np.log10(self.res_model[::-1, :, :])
            
        #write out the layers from resmodel
        for zz in range(self.nodes_z.shape[0]):
            ifid.write('\n')
            for ee in range(self.nodes_east.shape[0]):
                for nn in range(self.nodes_north.shape[0]):
                    ifid.write('{0:>13.5E}'.format(write_res_model[nn, ee, zz]))
                ifid.write('\n')
                
                
        if self.grid_center is None:
            #compute grid center
            center_east = -self.nodes_east.__abs__().sum()/2
            center_north = -self.nodes_norths.__abs__().sum()/2
            center_z = 0
            self.grid_center = np.array([center_north, center_east, center_z])
            
        ifid.write('\n{0:>16.3f}{1:>16.3f}{2:>16.3f}\n'.format(self.grid_center[0],
                   self.grid_center[1], self.grid_center[2]))
                   
        if self.strike_angle is None:
            ifid.write('{0:>9.3f}\n'.format(0))
        else:
            ifid.write('{0:>9.3f}\n'.format(self.strike_angle))
        ifid.close()
        
        print 'Wrote file to: {0}'.format(self.model_fn)
        
    def read_model_file(self, model_fn=None):
        """
        read an initial file and return the pertinent information including
        grid positions in coordinates relative to the center point (0,0) and 
        starting model.
        
        Note that the way the model file is output, it seems is that the 
        blocks are setup as 
        
        ModEM:                           WS:
        ----------                      ----- 
        0-----> N_north                 0-------->N_east
        |                               |
        |                               |
        V                               V
        N_east                          N_north
        
    
        Arguments:
        ----------
        
            **model_fn** : full path to initializing file.
            
        Outputs:
        --------
            
            **nodes_north** : np.array(nx)
                        array of nodes in S --> N direction
            
            **nodes_east** : np.array(ny) 
                        array of nodes in the W --> E direction
                        
            **nodes_z** : np.array(nz)
                        array of nodes in vertical direction positive downwards
            
            **res_model** : dictionary
                        dictionary of the starting model with keys as layers
                        
            **res_list** : list
                        list of resistivity values in the model
            
            **title** : string
                         title string
                           
        """
        
        if model_fn is not None:
            self.model_fn = model_fn
        
        if self.model_fn is None:
            raise ModEMError('model_fn is None, input a model file name')
        
        if os.path.isfile(self.model_fn) is None:
            raise ModEMError('Cannot find {0}, check path'.format(self.model_fn))
        
        self.save_path = os.path.dirname(self.model_fn)
        
        ifid = file(self.model_fn, 'r')    
        ilines = ifid.readlines()
        ifid.close()
        
        self.title = ilines[0].strip()
    
        #get size of dimensions, remembering that x is N-S, y is E-W, z is + down    
        nsize = ilines[1].strip().split()
        n_north = int(nsize[0])
        n_east = int(nsize[1])
        n_z = int(nsize[2])
        log_yn = nsize[4]
    
        #get nodes
        self.nodes_north = np.array([np.float(nn)
                                     for nn in ilines[2].strip().split()])
        self.nodes_east = np.array([np.float(nn)
                                     for nn in ilines[3].strip().split()])
        self.nodes_z = np.array([np.float(nn)
                                     for nn in ilines[4].strip().split()])

        self.res_model = np.zeros((n_north, n_east, n_z))
        
        #put the grids into coordinates relative to the center of the grid
        self.grid_north = self.nodes_north.copy()
        self.grid_north[:int(n_north/2)] =\
                        -np.array([self.nodes_north[ii:int(n_north/2)].sum() 
                                   for ii in range(int(n_north/2))])
        self.grid_north[int(n_north/2):] = \
                        np.array([self.nodes_north[int(n_north/2):ii+1].sum() 
                                 for ii in range(int(n_north/2), n_north)])-\
                                 self.nodes_north[int(n_north/2)]
                                
        self.grid_east = self.nodes_east.copy()
        self.grid_east[:int(n_east/2)] = \
                            -np.array([self.nodes_east[ii:int(n_east/2)].sum() 
                                       for ii in range(int(n_east/2))])
        self.grid_east[int(n_east/2):] = \
                            np.array([self.nodes_east[int(n_east/2):ii+1].sum() 
                                     for ii in range(int(n_east/2),n_east)])-\
                                     self.nodes_east[int(n_east/2)]
                                
        self.grid_z = np.array([self.nodes_z[:ii+1].sum() for ii in range(n_z)])
        
        #get model
        count_z = 0
        line_index= 6 
        count_e = 0
        while count_z < n_z:
            iline = ilines[line_index].strip().split()
            #blank lines spit the depth blocks, use those as a marker to 
            #set the layer number and start a new block
            if len(iline) == 0:
                count_z += 1
                count_e = 0
                line_index += 1
            #each line in the block is a line of N-->S values for an east value
            else:
                north_line = np.array([float(nres) for nres in 
                                     ilines[line_index].strip().split()])
                                     
                # Need to be sure that the resistivity array matches
                # with the grids, such that the first index is the 
                # furthest south 
                self.res_model[:, count_e, count_z] = north_line[::-1]

                count_e += 1
                line_index += 1
                
        #--> get grid center and rotation angle
        if len(ilines) > line_index:
            for iline in ilines[line_index:]:
                ilist = iline.strip().split()
                #grid center
                if len(ilist) == 3:
                    self.grid_center = np.array(ilist, dtype=np.float)
                #rotation angle
                elif len(ilist) == 1:
                    self.rotation_angle = np.float(ilist[0])
                else:
                    pass
        
        #--> make sure the resistivity units are in linear Ohm-m            
        if log_yn.lower() == 'loge':
            self.res_model = np.e**self.res_model
        elif log_yn.lower() == 'log' or log_yn.lower() == 'log10':
            self.res_model = 10**self.res_model
            
    def read_ws_model_file(self, ws_model_fn):
        """
        reads in a WS3INV3D model file
        """
        
        ws_model_obj = ws.WSModel(ws_model_fn)
        ws_model_obj.read_model_file()
        
        #set similar attributes
        for ws_key in ws_model_obj.__dict__.keys():
            for md_key in self.__dict__.keys():
                if ws_key == md_key:
                    setattr(self, ws_key, ws_model_obj.__dict__[ws_key])
                    
        #compute grid center
        center_east = -self.nodes_east.__abs__().sum()/2
        center_north = -self.nodes_norths.__abs__().sum()/2
        center_z = 0
        self.grid_center = np.array([center_north, center_east, center_z]) 
        
            
    def write_vtk_file(self, vtk_fn=None):
        """
        write a vtk file to view in Paraview or other
        """
        pass
            
#==============================================================================
# Control File
#==============================================================================
class Control(object):
    """
    read and write control file
    
    This file controls how the inversion starts and how it is run
    
    """
    
    def __init__(self, **kwargs):
        
        self.output_fn = kwargs.pop('output_fn', 'MODULAR_NLCG')
        self.lambda_initial = kwargs.pop('lambda_initial', 10)
        self.lambda_step = kwargs.pop('lambda_step', 10)
        self.model_search_step = kwargs.pop('model_search_step', 1)
        self.rms_reset_search = kwargs.pop('rms_reset_search', 2.0e-3)
        self.rms_target = kwargs.pop('rms_target', 1.00)
        self.lambda_exit = kwargs.pop('lambda_exit', 1.0e-4)
        self.max_iterations = kwargs.pop('max_iterations', 100)
        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.fn_basename = kwargs.pop('fn_basename', 'control.inv')
        self.control_fn = kwargs.pop('control_fn', os.path.join(self.save_path,
                                                            self.fn_basename))
                                                            
        self._control_keys = ['Model and data output file name',
                              'Initial damping factor lambda',
                              'To update lambda divide by',
                              'Initial search step in model units',
                              'Restart when rms diff is less than',
                              'Exit search when rms is less than',
                              'Exit when lambda is less than',
                              'Maximum number of iterations']
        
        self._control_dict = dict([(key, value) 
                                    for key, value in zip(self._control_keys,
                                    [self.output_fn, self.lambda_initial,
                                     self.lambda_step, self.model_search_step,
                                     self.rms_reset_search, self.rms_target,
                                     self.lambda_exit, self.max_iterations])])
        self._string_fmt_dict = dict([(key, value) 
                                    for key, value in zip(self._control_keys,
                                    ['<', '<.1f', '<.1f', '<.1f', '<.1e', 
                                     '<.2f', '<.1e', '<.0f'])])
                                     
    def write_control_file(self, control_fn=None, save_path=None, 
                           fn_basename=None):
        """
        write control file
        
        Arguments:
        ------------
            **control_fn** : string
                             full path to save control file to
                             *default* is save_path/fn_basename
            
            **save_path** : string
                            directory path to save control file to
                            *default* is cwd
            
            **fn_basename** : string
                              basename of control file
                              *default* is control.inv
                              
        """
        
        if control_fn is not None:
            self.save_path = os.path.dirname(control_fn)
            self.fn_basename = os.path.basename(control_fn)

        if save_path is not None:
            self.save_path = save_path
            
        if fn_basename is not None:
            self.fn_basename = fn_basename
            
        self.control_fn = os.path.join(self.save_path, self.fn_basename)
        
        self._control_dict = dict([(key, value) 
                                    for key, value in zip(self._control_keys,
                                    [self.output_fn, self.lambda_initial,
                                     self.lambda_step, self.model_search_step,
                                     self.rms_reset_search, self.rms_target,
                                     self.lambda_exit, self.max_iterations])])
        
        clines = []
        for key in self._control_keys:
            value = self._control_dict[key]
            str_fmt = self._string_fmt_dict[key]
            clines.append('{0:<35}: {1:{2}}\n'.format(key, value, str_fmt))
            
        cfid = file(self.control_fn, 'w')
        cfid.writelines(clines)
        cfid.close()
        
        print 'Wrote ModEM control file to {0}'.format(self.control_fn)
        
    def read_control_file(self, control_fn=None):
        """
        read in a control file
        """
        
        if control_fn is not None:
            self.control_fn = control_fn

        if self.control_fn is None:
            raise mtex.MTpyError_file_handling('control_fn is None, input '
                                                'control file')
            
        if os.path.isfile(self.control_fn) is False:
            raise mtex.MTpyError_file_handling('Could not find {0}'.format(
                                                self.control_fn))
                                                
        self.save_path = os.path.dirname(self.control_fn)
        self.fn_basename = os.path.basename(self.control_fn)
        
        cfid = file(self.control_fn, 'r')
        clines = cfid.readlines()
        cfid.close()
        for cline in clines:
            clist = cline.strip().split(':')
            if len(clist) == 2:

                try:
                    self._control_dict[clist[0].strip()] = float(clist[1])
                except ValueError:
                    self._control_dict[clist[0].strip()] = clist[1]

        #set attributes
        attr_list = ['output_fn', 'lambda_initial','lambda_step',
                     'model_search_step','rms_reset_search','rms_target',
                     'lambda_exit','max_iterations']
        for key, kattr in zip(self._control_keys, attr_list):
            setattr(self, kattr, self._control_dict[key])
            
            
#==============================================================================
# covariance 
#==============================================================================
class Covariance(object):
    """
    read and write covariance files
    
    """

    def __init__(self, grid_dimensions=None, **kwargs): 
            
        self.grid_dimensions = grid_dimensions
        self.smoothing_east = kwargs.pop('smoothing_east', 0.3)
        self.smoothing_north = kwargs.pop('smoothing_north', 0.3)
        self.smoothing_z = kwargs.pop('smoothing_z', 0.3)
        self.smoothing_num = kwargs.pop('smoothing_num', 1)
        
        self.exception_list = kwargs.pop('exception_list', [])
        self.mask_arr = kwargs.pop('mask_arr', None)
        
        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.cov_fn_basename = kwargs.pop('cov_fn_basename', 'covariance.cov')
        
        self.cov_fn = kwargs.pop('cov_fn', None)
                                                  
        self._header_str = '\n'.join(['+{0}+'.format('-'*77),
    '| This file defines model covariance for a recursive autoregression scheme.   |',
    '| The model space may be divided into distinct areas using integer masks.     |',
    '| Mask 0 is reserved for air; mask 9 is reserved for ocean. Smoothing between |',
    '| air, ocean and the rest of the model is turned off automatically. You can   |',
    '| also define exceptions to override smoothing between any two model areas.   |',
    '| To turn off smoothing set it to zero.  This header is 16 lines long.        |',
    '| 1. Grid dimensions excluding air layers (Nx, Ny, NzEarth)                   |',
    '| 2. Smoothing in the X direction (NzEarth real values)                       |',
    '| 3. Smoothing in the Y direction (NzEarth real values)                       |',
    '| 4. Vertical smoothing (1 real value)                                        |',
    '| 5. Number of times the smoothing should be applied (1 integer >= 0)         |',
    '| 6. Number of exceptions (1 integer >= 0)                                    |', 
    '| 7. Exceptions in the for e.g. 2 3 0. (to turn off smoothing between 3 & 4)  |',
    '| 8. Two integer layer indices and Nx x Ny block of masks, repeated as needed.|',
    '+{0}+'.format('-'*77)])
    
    def write_covariance_file(self, cov_fn=None, save_path=None, 
                              cov_fn_basename=None):
        """
        write a covariance file
        """
        
        if self.grid_dimensions is None:
            raise ModEMError('Grid dimensions are None, input as (Nx, Ny, Nz)')
        if cov_fn is not None:
            self.cov_fn = cov_fn
        else:
            if save_path is not None:
                self.save_path = save_path
            if cov_fn_basename is not None:
                self.cov_fn_basename = cov_fn_basename
            self.cov_fn = os.path.join(self.save_path, self.cov_fn_basename)
            
        clines = [self._header_str]
        clines.append('\n\n')
        
        #--> grid dimensions
        clines.append(' {0:<10}{1:<10}{2:<10}\n'.format(self.grid_dimensions[0],
                                                       self.grid_dimensions[1],
                                                       self.grid_dimensions[2]))
        clines.append('\n')
        
        #--> smoothing in north direction        
        n_smooth_line = ''
        for zz in range(self.grid_dimensions[2]):
            n_smooth_line += ' {0:<5.1f}'.format(self.smoothing_north)
        clines.append(n_smooth_line+'\n')

        #--> smoothing in east direction
        e_smooth_line = ''
        for zz in range(self.grid_dimensions[2]):
            e_smooth_line += ' {0:<5.1f}'.format(self.smoothing_east)
        clines.append(n_smooth_line+'\n')
        
        #--> smoothing in vertical direction
        clines.append(' {0:<5.1f}\n'.format(self.smoothing_z))
        clines.append('\n')
        
        #--> number of times to apply smoothing
        clines.append(' {0:<2.0f}\n'.format(self.smoothing_num))
        clines.append('\n')
        
        #--> exceptions
        clines.append(' {0:<.0f}\n'.format(len(self.exception_list)))
        for exc in self.exception_list:
            clines.append('{0:<5.0f}{1:<5.0f}{2:<5.0f}\n'.format(exc[0],
                                                                 exc[1],
                                                                 exc[2]))
        clines.append('\n')
        clines.append('\n')
        #--> mask array
        if self.mask_arr is None:
            self.mask_arr = np.ones((self.grid_dimensions[0],
                                     self.grid_dimensions[1],
                                     self.grid_dimensions[2]))
        for zz in range(self.mask_arr.shape[2]):
            clines.append(' {0:<8.0f}{0:<8.0f}\n'.format(zz+1))
            for nn in range(self.mask_arr.shape[0]):
                cline = ''
                for ee in range(self.mask_arr.shape[1]):
                    cline += '{0:^3.0f}'.format(self.mask_arr[nn, ee, zz])
                clines.append(cline+'\n')
        
        cfid = file(self.cov_fn, 'w')
        cfid.writelines(clines)
        cfid.close()
        
        print 'Wrote covariance file to {0}'.format(self.cov_fn)
        
            
        

#==============================================================================
# Manipulate the model to test structures or create a starting model
#==============================================================================
class ModelManipulator(Model):
    """
    will plot a model from wsinv3d or init file so the user can manipulate the 
    resistivity values relatively easily.  At the moment only plotted
    in map view.
    
    
    :Example: ::
        >>> import mtpy.modeling.ws3dinv as ws
        >>> initial_fn = r"/home/MT/ws3dinv/Inv1/WSInitialFile"
        >>> mm = ws.WSModelManipulator(initial_fn=initial_fn)
        
    =================== =======================================================
    Buttons              Description    
    =================== =======================================================
    '='                 increase depth to next vertical node (deeper)
    '-'                 decrease depth to next vertical node (shallower)
    'q'                 quit the plot, rewrites initial file when pressed
    'a'                 copies the above horizontal layer to the present layer
    'b'                 copies the below horizonal layer to present layer
    'u'                 undo previous change
    =================== =======================================================
    
    
    =================== =======================================================
    Attributes          Description
    =================== =======================================================
    ax1                 matplotlib.axes instance for mesh plot of the model 
    ax2                 matplotlib.axes instance of colorbar
    cb                  matplotlib.colorbar instance for colorbar 
    cid_depth           matplotlib.canvas.connect for depth
    cmap                matplotlib.colormap instance
    cmax                maximum value of resistivity for colorbar. (linear)
    cmin                minimum value of resistivity for colorbar (linear)
    data_fn             full path fo data file
    depth_index         integer value of depth slice for plotting
    dpi                 resolution of figure in dots-per-inch
    dscale              depth scaling, computed internally
    east_line_xlist     list of east mesh lines for faster plotting
    east_line_ylist     list of east mesh lines for faster plotting
    fdict               dictionary of font properties
    fig                 matplotlib.figure instance
    fig_num              number of figure instance
    fig_size             size of figure in inches
    font_size           size of font in points
    grid_east           location of east nodes in relative coordinates
    grid_north          location of north nodes in relative coordinates
    grid_z              location of vertical nodes in relative coordinates
    initial_fn          full path to initial file
    m_height            mean height of horizontal cells
    m_width             mean width of horizontal cells
    map_scale            [ 'm' | 'km' ] scale of map
    mesh_east           np.meshgrid of east, north
    mesh_north          np.meshgrid of east, north
    mesh_plot           matplotlib.axes.pcolormesh instance
    model_fn            full path to model file
    new_initial_fn      full path to new initial file
    nodes_east          spacing between east nodes 
    nodes_north         spacing between north nodes 
    nodes_z             spacing between vertical nodes
    north_line_xlist    list of coordinates of north nodes for faster plotting
    north_line_ylist    list of coordinates of north nodes for faster plotting
    plot_yn             [ 'y' | 'n' ] plot on instantiation
    radio_res           matplotlib.widget.radio instance for change resistivity
    rect_selector       matplotlib.widget.rect_selector 
    res                 np.ndarray(nx, ny, nz) for model in linear resistivity
    res_copy            copy of res for undo
    res_dict            dictionary of segmented resistivity values 
    res_list            list of resistivity values for model linear scale
    res_model           np.ndarray(nx, ny, nz) of resistivity values from 
                        res_list (linear scale)
    res_model_int       np.ndarray(nx, ny, nz) of integer values corresponding
                        to res_list for initial model
    res_value           current resistivty value of radio_res
    save_path           path to save initial file to
    station_east        station locations in east direction
    station_north       station locations in north direction
    xlimits             limits of plot in e-w direction
    ylimits             limits of plot in n-s direction
    =================== =======================================================

    """

    def __init__(self, model_fn=None, data_fn=None, **kwargs):
        
        #be sure to initialize Model
        Model.__init__(self, model_fn=model_fn, **kwargs)
        
        self.data_fn = data_fn
        self.model_fn_basename = kwargs.pop('model_fn_basename', 
                                             'ModeEM_Model_rw.ws')
        
        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.data_fn is not None:
            self.save_path = os.path.dirname(self.data_fn)
        else:
            self.save_path = os.getcwd()
        
        
        #station locations in relative coordinates read from data file
        self.station_east = None
        self.station_north = None
        
        #--> set map scale
        self.map_scale = kwargs.pop('map_scale', 'km')
        
        self.m_width = 100
        self.m_height = 100
        
        #--> scale the map coordinates
        if self.map_scale=='km':
            self.dscale = 1000.
        if self.map_scale=='m':
            self.dscale = 1.
            
        #figure attributes
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.cb = None
        self.east_line_xlist = None
        self.east_line_ylist = None
        self.north_line_xlist = None
        self.north_line_ylist = None
        
        #make a default resistivity list to change values
        self._res_sea = 0.3
        self._res_air = 1E12
        self.res_dict = None
        self.res_list = kwargs.pop('res_list', None)
        if self.res_list is None:
            self.set_res_list(np.array([self._res_sea, 1, 10, 50, 100, 500, 
                                        1000, 5000, self._res_air],
                                      dtype=np.float))

        #set initial resistivity value
        self.res_value = self.res_list[0]
        self.cov_arr = None
        
        #--> set map limits
        self.xlimits = kwargs.pop('xlimits', None)
        self.ylimits = kwargs.pop('ylimits', None)

        self.font_size = kwargs.pop('font_size', 7)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.cmap = kwargs.pop('cmap', cm.jet_r)
        self.depth_index = kwargs.pop('depth_index', 0)
        
        self.fdict = {'size':self.font_size+2, 'weight':'bold'}

        #plot on initialization
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn=='y':
            self.get_model()
            self.plot()
            
    def set_res_list(self, res_list):
        """
        on setting res_list also set the res_dict to correspond
        """
        self.res_list = res_list
        #make a dictionary of values to write to file.
        self.res_dict = dict([(res, ii) 
                              for ii, res in enumerate(self.res_list,1)])
        if self.fig is not None:
            plt.close()
            self.plot()
        
    
    #---read files-------------------------------------------------------------    
    def get_model(self):
        """
        reads in initial file or model file and set attributes:
            -resmodel
            -northrid
            -eastrid
            -zgrid
            -res_list if initial file
            
        """
        #--> read in model file        
        self.read_model_file()
        
        self.cov_arr = np.ones_like(self.res_model)
        
        #--> read in data file if given
        if self.data_fn is not None:
            md_data = Data()
            md_data.read_data_file(self.data_fn)
            
            #get station locations
            self.station_east = md_data.coord_array['rel_east']
            self.station_north = md_data.coord_array['rel_north']
            
        #get cell block sizes
        self.m_height = np.median(self.nodes_north[5:-5])/self.dscale
        self.m_width = np.median(self.nodes_east[5:-5])/self.dscale
            
        #make a copy of original in case there are unwanted changes
        self.res_copy = self.res_model.copy()
            
    #---plot model-------------------------------------------------------------    
    def plot(self):
        """
        plots the model with:
            -a radio dial for depth slice 
            -radio dial for resistivity value
            
        """
        #make sure there is a model to plot
        if self.res_model is None:
            self.get_model()
            
        self.cmin = np.floor(np.log10(min(self.res_list)))
        self.cmax = np.ceil(np.log10(max(self.res_list)))
        
        #-->Plot properties
        plt.rcParams['font.size'] = self.font_size
        
        #need to add an extra row and column to east and north to make sure 
        #all is plotted see pcolor for details.
        plot_east = np.append(self.grid_east, self.grid_east[-1]*1.25)/self.dscale
        plot_north = np.append(self.grid_north, self.grid_north[-1]*1.25)/self.dscale
        
        #make a mesh grid for plotting
        #the 'ij' makes sure the resulting grid is in east, north
        self.mesh_east, self.mesh_north = np.meshgrid(plot_east, 
                                                      plot_north,
                                                      indexing='ij')
        
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()
        self.ax1 = self.fig.add_subplot(1, 1, 1, aspect='equal')
        
        #transpose to make x--east and y--north
        plot_res = np.log10(self.res_model[:,:,self.depth_index].T)
        
        self.mesh_plot = self.ax1.pcolormesh(self.mesh_east,
                                             self.mesh_north, 
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)
                                             
        #on plus or minus change depth slice
        self.cid_depth = \
                    self.mesh_plot.figure.canvas.mpl_connect('key_press_event',
                                                        self._on_key_callback)
                                    
                       
        #plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax1.text(ee/self.dscale, nn/self.dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size':self.font_size-2,
                                        'weight':'bold'})

        #set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(xmin=self.grid_east.min()/self.dscale, 
                              xmax=self.grid_east.max()/self.dscale)
        
        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(ymin=self.grid_north.min()/self.dscale,
                              ymax=self.grid_north.max()/self.dscale)
            
        #self.ax1.xaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        #self.ax1.yaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        
        self.ax1.set_ylabel('Northing ('+self.map_scale+')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting ('+self.map_scale+')',
                            fontdict=self.fdict)
        
        depth_title = self.grid_z[self.depth_index]/self.dscale
                                                        
        self.ax1.set_title('Depth = {:.3f} '.format(depth_title)+\
                           '('+self.map_scale+')',
                           fontdict=self.fdict)
        
        #plot the grid if desired  
        self.east_line_xlist = []
        self.east_line_ylist = []            
        for xx in self.grid_east:
            self.east_line_xlist.extend([xx/self.dscale, xx/self.dscale])
            self.east_line_xlist.append(None)
            self.east_line_ylist.extend([self.grid_north.min()/self.dscale, 
                                         self.grid_north.max()/self.dscale])
            self.east_line_ylist.append(None)
        self.ax1.plot(self.east_line_xlist,
                      self.east_line_ylist,
                       lw=.25,
                       color='k')

        self.north_line_xlist = []
        self.north_line_ylist = [] 
        for yy in self.grid_north:
            self.north_line_xlist.extend([self.grid_east.min()/self.dscale,
                                          self.grid_east.max()/self.dscale])
            self.north_line_xlist.append(None)
            self.north_line_ylist.extend([yy/self.dscale, yy/self.dscale])
            self.north_line_ylist.append(None)
        self.ax1.plot(self.north_line_xlist,
                      self.north_line_ylist,
                      lw=.25,
                      color='k')
        
        #plot the colorbar
        self.ax2 = mcb.make_axes(self.ax1, orientation='vertical', shrink=.35)
        seg_cmap = ws.cmap_discretize(self.cmap, len(self.res_list))
        self.cb = mcb.ColorbarBase(self.ax2[0],cmap=seg_cmap,
                                   norm=colors.Normalize(vmin=self.cmin,
                                                         vmax=self.cmax))
                                                         
                            
        self.cb.set_label('Resistivity ($\Omega \cdot$m)',
                          fontdict={'size':self.font_size})
        self.cb.set_ticks(np.arange(self.cmin, self.cmax+1))
        self.cb.set_ticklabels([mtplottools.labeldict[cc] 
                                for cc in np.arange(self.cmin, self.cmax+1)])
                            
        #make a resistivity radio button
        resrb = self.fig.add_axes([.85,.1,.1,.2])
        reslabels = ['{0:.4g}'.format(res) for res in self.res_list]
        self.radio_res = widgets.RadioButtons(resrb, reslabels, 
                                        active=self.res_dict[self.res_value])
        
        #make a rectangular selector
        self.rect_selector = widgets.RectangleSelector(self.ax1, 
                                                       self.rect_onselect,
                                                       drawtype='box',
                                                       useblit=True)

        
        plt.show()
        
        #needs to go after show()
        self.radio_res.on_clicked(self.set_res_value)


    def redraw_plot(self):
        """
        redraws the plot
        """
        
        current_xlimits = self.ax1.get_xlim()
        current_ylimits = self.ax1.get_ylim()
        
        self.ax1.cla()
        
        plot_res = np.log10(self.res_model[:,:,self.depth_index].T)
        
        self.mesh_plot = self.ax1.pcolormesh(self.mesh_east,
                                             self.mesh_north, 
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)
                                             
         #plot the stations
        if self.station_east is not None:
            for ee,nn in zip(self.station_east, self.station_north):
                self.ax1.text(ee/self.dscale, nn/self.dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size':self.font_size-2,
                                        'weight':'bold'})

        #set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(current_xlimits)
        
        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(current_ylimits)
        
        self.ax1.set_ylabel('Northing ('+self.map_scale+')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting ('+self.map_scale+')',
                            fontdict=self.fdict)
        
        depth_title = self.grid_z[self.depth_index]/self.dscale
                                                        
        self.ax1.set_title('Depth = {:.3f} '.format(depth_title)+\
                           '('+self.map_scale+')',
                           fontdict=self.fdict)
                     
        #plot finite element mesh
        self.ax1.plot(self.east_line_xlist,
                      self.east_line_ylist,
                      lw=.25,
                      color='k')

        
        self.ax1.plot(self.north_line_xlist,
                      self.north_line_ylist,
                      lw=.25,
                      color='k')
        
        #be sure to redraw the canvas                  
        self.fig.canvas.draw()
        
    def set_res_value(self, label):
        self.res_value = float(label)
        print 'set resistivity to ', label
        print self.res_value
        
        
    def _on_key_callback(self,event):
        """
        on pressing a key do something
        
        """
        
        self.event_change_depth = event
        
        #go down a layer on push of +/= keys
        if self.event_change_depth.key == '=':
            self.depth_index += 1
            
            if self.depth_index>len(self.grid_z)-1:
                self.depth_index = len(self.grid_z)-1
                print 'already at deepest depth'
                
            print 'Plotting Depth {0:.3f}'.format(self.grid_z[self.depth_index]/\
                    self.dscale)+'('+self.map_scale+')'
            
            self.redraw_plot()
        #go up a layer on push of - key
        elif self.event_change_depth.key == '-':
            self.depth_index -= 1
            
            if self.depth_index < 0:
                self.depth_index = 0
                
            print 'Plotting Depth {0:.3f} '.format(self.grid_z[self.depth_index]/\
                    self.dscale)+'('+self.map_scale+')'
            
            self.redraw_plot()
        
        #exit plot on press of q
        elif self.event_change_depth.key == 'q':
            self.event_change_depth.canvas.mpl_disconnect(self.cid_depth)
            plt.close(self.event_change_depth.canvas.figure)
            self.rewrite_model_file()
            
        #copy the layer above
        elif self.event_change_depth.key == 'a':
            try:
                if self.depth_index == 0:
                    print 'No layers above'
                else:
                    self.res_model[:, :, self.depth_index] = \
                                       self.res_model[:, :, self.depth_index-1]
            except IndexError:
                print 'No layers above'
                
            self.redraw_plot()
        
        #copy the layer below
        elif self.event_change_depth.key == 'b':
            try:
                self.res_model[:, :, self.depth_index] = \
                                    self.res_model[:, :, self.depth_index+1]
            except IndexError:
                print 'No more layers below'
                
            self.redraw_plot() 
            
        #undo
        elif self.event_change_depth.key == 'u':
            if type(self.xchange) is int and type(self.ychange) is int:
                self.res_model[self.ychange, self.xchange, self.depth_index] =\
                self.res_copy[self.ychange, self.xchange, self.depth_index]
            else:
                for xx in self.xchange:
                    for yy in self.ychange:
                        self.res_model[yy, xx, self.depth_index] = \
                        self.res_copy[yy, xx, self.depth_index]
            
            self.redraw_plot()
            
    def change_model_res(self, xchange, ychange):
        """
        change resistivity values of resistivity model
        
        """
        if type(xchange) is int and type(ychange) is int:
            self.res_model[ychange, xchange, self.depth_index] = self.res_value
        else:
            for xx in xchange:
                for yy in ychange:
                    self.res_model[yy, xx, self.depth_index] = self.res_value
        
        self.redraw_plot()            
           
    def rect_onselect(self, eclick, erelease):
        """
        on selecting a rectangle change the colors to the resistivity values
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        
        self.xchange = self._get_east_index(x1, x2)
        self.ychange = self._get_north_index(y1, y2)
        
        #reset values of resistivity
        self.change_model_res(self.xchange, self.ychange)
        
        
    def _get_east_index(self, x1, x2):
        """
        get the index value of the points to be changed
        
        """
        if x1 < x2:
            xchange = np.where((self.grid_east/self.dscale >= x1) & \
                               (self.grid_east/self.dscale <= x2))[0]
            if len(xchange) == 0:
                xchange = np.where(self.grid_east/self.dscale >= x1)[0][0]-1
                return [xchange]
                
        if x1 > x2:
            xchange = np.where((self.grid_east/self.dscale <= x1) & \
                               (self.grid_east/self.dscale >= x2))[0]
            if len(xchange) == 0:
                xchange = np.where(self.grid_east/self.dscale >= x2)[0][0]-1
                return [xchange]

            
        #check the edges to see if the selection should include the square
        xchange = np.append(xchange, xchange[0]-1)
        xchange.sort()

        return xchange
                
    def _get_north_index(self, y1, y2):
        """
        get the index value of the points to be changed in north direction
        
        need to flip the index because the plot is flipped
        
        """
        
        if y1 < y2:
            ychange = np.where((self.grid_north/self.dscale > y1) & \
                               (self.grid_north/self.dscale < y2))[0]
            if len(ychange) == 0:
                ychange = np.where(self.grid_north/self.dscale >= y1)[0][0]-1
                return [ychange]
                
        elif y1 > y2:
            ychange = np.where((self.grid_north/self.dscale < y1) & \
                               (self.grid_north/self.dscale > y2))[0]
            if len(ychange) == 0:
                ychange = np.where(self.grid_north/self.dscale >= y2)[0][0]-1
                return [ychange]
        
        ychange -= 1
        ychange = np.append(ychange, ychange[-1]+1)

        return ychange

                
    def rewrite_model_file(self, model_fn=None, save_path=None, 
                           model_fn_basename=None):
        """
        write an initial file for wsinv3d from the model created.
        """
        if save_path is not None:
            self.save_path = save_path
       
        self.model_fn = model_fn
        
        if model_fn_basename is not None:
            self.model_fn_basename = model_fn_basename
            
        self.write_model_file()
                                                         
#==============================================================================
# plot response       
#==============================================================================
class PlotResponse(object):
    """
    plot data and response 
    
    Plots the real and imaginary impedance and induction vector if present.
    
    :Example: ::
    
        >>> import mtpy.modeling.new_modem as modem
        >>> dfn = r"/home/MT/ModEM/Inv1/DataFile.dat"
        >>> rfn = r"/home/MT/ModEM/Inv1/Test_resp_000.res"
        >>> mrp = modem.PlotResponse(data_fn=dfn, resp_fn=rfn)
        >>> # plot only the TE and TM modes
        >>> mrp.plot_component = 2
        >>> mrp.redraw_plot()
    
    ======================== ==================================================
    Attributes               Description    
    ======================== ==================================================
    color_mode               [ 'color' | 'bw' ] color or black and white plots
    cted                     color for data TE mode
    ctem                     color for data TM mode
    ctmd                     color for model TE mode
    ctmm                     color for model TM mode
    data_fn                  full path to data file
    data_object              WSResponse instance
    e_capsize                cap size of error bars in points (*default* is .5)
    e_capthick               cap thickness of error bars in points (*default*
                             is 1)
    fig_dpi                  resolution of figure in dots-per-inch (300)
    fig_list                 list of matplotlib.figure instances for plots
    fig_size                 size of figure in inches (*default* is [6, 6])
    font_size                size of font for tick labels, axes labels are
                             font_size+2 (*default* is 7)
    legend_border_axes_pad   padding between legend box and axes 
    legend_border_pad        padding between border of legend and symbols
    legend_handle_text_pad   padding between text labels and symbols of legend
    legend_label_spacing     padding between labels
    legend_loc               location of legend 
    legend_marker_scale      scale of symbols in legend
    lw                       line width response curves (*default* is .5)
    ms                       size of markers (*default* is 1.5)
    mted                     marker for data TE mode
    mtem                     marker for data TM mode
    mtmd                     marker for model TE mode
    mtmm                     marker for model TM mode 
    phase_limits             limits of phase
    plot_component           [ 2 | 4 ] 2 for TE and TM or 4 for all components
    plot_style               [ 1 | 2 ] 1 to plot each mode in a seperate
                             subplot and 2 to plot xx, xy and yx, yy in same 
                             plots
    plot_type                [ '1' | list of station name ] '1' to plot all 
                             stations in data file or input a list of station
                             names to plot if station_fn is input, otherwise
                             input a list of integers associated with the 
                             index with in the data file, ie 2 for 2nd station
    plot_z                   [ True | False ] *default* is True to plot 
                             impedance, False for plotting resistivity and 
                             phase
    plot_yn                  [ 'n' | 'y' ] to plot on instantiation
    res_limits               limits of resistivity in linear scale
    resp_fn                  full path to response file
    resp_object              WSResponse object for resp_fn, or list of 
                             WSResponse objects if resp_fn is a list of
                             response files
    station_fn               full path to station file written by WSStation
    subplot_bottom           space between axes and bottom of figure
    subplot_hspace           space between subplots in vertical direction
    subplot_left             space between axes and left of figure
    subplot_right            space between axes and right of figure
    subplot_top              space between axes and top of figure
    subplot_wspace           space between subplots in horizontal direction    
    ======================== ==================================================
    """
    
    def __init__(self, data_fn=None, resp_fn=None, **kwargs):
        self.data_fn = data_fn
        self.resp_fn = resp_fn
        
        self.data_object = None
        self.resp_object = []
        
        self.color_mode = kwargs.pop('color_mode', 'color')
        
        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        #color mode
        if self.color_mode == 'color':
            #color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')
            
            #color for occam2d model
            self.ctem = kwargs.pop('ctem', (0, .6, .3))
            self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')
         
        #black and white mode
        elif self.color_mode == 'bw':
            #color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', '*')
            self.mtmd = kwargs.pop('mtmd', 'v')
            
            #color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')
            
        self.phase_limits = kwargs.pop('phase_limits', None)
        self.res_limits = kwargs.pop('res_limits', None)

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        
        self.subplot_wspace = .2
        self.subplot_hspace = .0
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .93
        self.subplot_bottom = .1
        
        self.legend_loc = 'upper left'
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15

        self.font_size = kwargs.pop('font_size', 6)
        
        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_style = kwargs.pop('plot_style', 1)
        self.plot_component = kwargs.pop('plot_component', 4)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.plot_z = kwargs.pop('plot_z', True)
        self.ylabel_pad = kwargs.pop('ylabel_pad', 1.25)
        
        self.fig_list = []
        
        if self.plot_yn == 'y':
            self.plot()
    
    def plot(self):
        """
        plot
        """
        
        self.data_object = Data()
        self.data_object.read_data_file(self.data_fn)
                                        
        #get shape of impedance tensors
        ns = len(self.data_object.mt_dict.keys())
    
        #read in response files
        if self.resp_fn != None:
            self.resp_object = []
            if type(self.resp_fn) is not list:
                resp_obj = Data()
                resp_obj.read_data_file(self.resp_fn)
                self.resp_object = [resp_obj]
            else:
                for rfile in self.resp_fn:
                    resp_obj = Data()
                    resp_obj.read_data_file(rfile)
                    self.resp_object.append(resp_obj)

        #get number of response files
        nr = len(self.resp_object)
        
        if type(self.plot_type) is list:
            ns = len(self.plot_type)
          
        #--> set default font size                           
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size':self.font_size+2, 'weight':'bold'} 
        if self.plot_z == True:
            h_ratio = [1,1]
        elif self.plot_z == False:
            h_ratio = [2, 1.5]
        
        ax_list = []
        line_list = []
        label_list = []
        
        
        #--> make key word dictionaries for plotting
        kw_xx = {'color':self.cted,
                 'marker':self.mted,
                 'ms':self.ms,
                 'ls':':',
                 'lw':self.lw,
                 'e_capsize':self.e_capsize,
                 'e_capthick':self.e_capthick}        
       
        kw_yy = {'color':self.ctmd,
                 'marker':self.mtmd,
                 'ms':self.ms,
                 'ls':':',
                 'lw':self.lw,
                 'e_capsize':self.e_capsize,
                 'e_capthick':self.e_capthick}        
        
        if self.plot_type != '1':
            pstation_list = []
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]
            for ii, station in enumerate(self.data_object.mt_dict.keys()):
                if type(station) is not int:
                    for pstation in self.plot_type:
                        if station.find(str(pstation)) >= 0:
                            pstation_list.append(station)
                else:
                    for pstation in self.plot_type:
                        if station == int(pstation):
                            pstation_list.append(ii)
        else:
            pstation_list = self.data_object.mt_dict.keys()
        
        for jj, station in enumerate(pstation_list):
            z_obj = self.data_object.mt_dict[station].Z
            t_obj = self.data_object.mt_dict[station].Tipper
            period = self.data_object.period_list
            print 'Plotting: {0}'.format(station)
            
            #convert to apparent resistivity and phase
            rp = mtplottools.ResPhase(z_object=z_obj)
            
            #find locations where points have been masked
            nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
            nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
            nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
            nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
            ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
            nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]
            
            if self.resp_fn != None:
                plotr = True
            else:
                plotr = False
            
            #make figure 
            fig = plt.figure(station, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            fig.suptitle(str(station), fontdict=fontdict)
            
            #set the grid of subplots
            plot_tipper = t_obj.tipper.all() == 0.0
            if plot_tipper == False:
                #makes more sense if plot_tipper is True to plot tipper
                plot_tipper = True
            else:
                plot_tipper = False
                
            if plot_tipper == True:
                gs = gridspec.GridSpec(2, 6,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom, 
                                   right=self.subplot_right, 
                                   hspace=self.subplot_hspace,
                                   height_ratios=h_ratio)
            else:
                gs = gridspec.GridSpec(2, 4,
                                       wspace=self.subplot_wspace,
                                       left=self.subplot_left,
                                       top=self.subplot_top,
                                       bottom=self.subplot_bottom, 
                                       right=self.subplot_right, 
                                       hspace=self.subplot_hspace,
                                       height_ratios=h_ratio)
            #---------plot the apparent resistivity-----------------------------------
            #plot each component in its own subplot
            if self.plot_style == 1:
                #plot xy and yx 
                if self.plot_component == 2:
                    if plot_tipper == False:
                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axryx = fig.add_subplot(gs[0, 2:], sharex=axrxy)
                        
                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)
                        axpyx = fig.add_subplot(gs[1, 2:], sharex=axrxy)
                    else:
                        
                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axryx = fig.add_subplot(gs[0, 2:4], sharex=axrxy)
                        
                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)
                        axpyx = fig.add_subplot(gs[1, 2:4], sharex=axrxy)
                        
                        axtr = fig.add_subplot(gs[0, 4:], sharex=axrxy)
                        axti = fig.add_subplot(gs[1, 4:], sharex=axrxy)
                        
                        
                    if self.plot_z == False: 
                        #plot resistivity
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                         period,
                                                         rp.resxy[nzxy],
                                                         rp.resxy_err[nzxy],
                                                         **kw_xx)   

                        eryx = mtplottools.plot_errorbar(axryx, 
                                                         period[nzyx], 
                                                         rp.resyx[nzyx], 
                                                         rp.resyx_err[nzyx],
                                                         **kw_yy)
                        #plot phase                         
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                         period[nzxy], 
                                                         rp.phasexy[nzxy], 
                                                         rp.phasexy_err[nzxy],
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx, 
                                                         period[nzyx], 
                                                         rp.phaseyx[nzyx], 
                                                         rp.phaseyx_err[nzyx],
                                                         **kw_yy)
                        
                    elif self.plot_z == True:
                        #plot real
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].real, 
                                                  z_obj.zerr[nzxy,0,1].real,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axryx, 
                                                     period[nzyx], 
                                                     z_obj.z[nzyx,1,0].real, 
                                                     z_obj.zerr[nzyx,1,0].real,
                                                     **kw_yy)
                        #plot phase                         
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].imag, 
                                                  z_obj.zerr[nzxy,0,1].imag,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx, 
                                                  period[nzyx], 
                                                  z_obj.z[nzyx,1,0].imag, 
                                                  z_obj.zerr[nzyx,1,0].imag,
                                                  **kw_yy)
                    #plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].real,
                                                 t_obj.tippererr[ntx, 0, 0],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].real,
                                                 t_obj.tippererr[nty, 0, 1],
                                                 **kw_yy)
                                                 
                        ertx = mtplottools.plot_errorbar(axti, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].imag,
                                                 t_obj.tippererr[ntx, 0, 0],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axti, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].imag,
                                                 t_obj.tippererr[nty, 0, 1],
                                                 **kw_yy)
                    
                    if plot_tipper == False:                          
                        ax_list = [axrxy, axryx, axpxy, axpyx]
                        line_list = [[erxy[0]], [eryx[0]]]
                        label_list = [['$Z_{xy}$'], ['$Z_{yx}$']]
                    else:                          
                        ax_list = [axrxy, axryx, axpxy, axpyx, axtr, axti]
                        line_list = [[erxy[0]], [eryx[0]], 
                                     [ertx[0], erty[0]]]
                        label_list = [['$Z_{xy}$'], ['$Z_{yx}$'],
                                       ['$T_{x}$', '$T_{y}$']]
                                                           
                elif self.plot_component == 4:
                    if plot_tipper == False:
                        axrxx = fig.add_subplot(gs[0, 0])
                        axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
                        axryx = fig.add_subplot(gs[0, 2], sharex=axrxx)
                        axryy = fig.add_subplot(gs[0, 3], sharex=axrxx)
                        
                        axpxx = fig.add_subplot(gs[1, 0])
                        axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
                        axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
                        axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)
                    else:
                        axrxx = fig.add_subplot(gs[0, 0])
                        axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
                        axryx = fig.add_subplot(gs[0, 2], sharex=axrxx)
                        axryy = fig.add_subplot(gs[0, 3], sharex=axrxx)
                        
                        axpxx = fig.add_subplot(gs[1, 0])
                        axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
                        axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
                        axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)
                        
                        axtxr = fig.add_subplot(gs[0, 4], sharex=axrxx)
                        axtxi = fig.add_subplot(gs[1, 4], sharex=axrxx)
                        axtyr = fig.add_subplot(gs[0, 5], sharex=axrxx)
                        axtyi = fig.add_subplot(gs[1, 5], sharex=axrxx)
                    
                    if self.plot_z == False:
                        #plot resistivity
                        erxx= mtplottools.plot_errorbar(axrxx, 
                                                  period[nzxx], 
                                                  rp.resxx[nzxx], 
                                                  rp.resxx_err[nzxx],
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  rp.resxy[nzxy], 
                                                  rp.resxy_err[nzxy],
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axryx, 
                                                  period[nzyx], 
                                                  rp.resyx[nzyx], 
                                                  rp.resyx_err[nzyx],
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axryy, 
                                                  period[nzyy], 
                                                  rp.resyy[nzyy], 
                                                  rp.resyy_err[nzyy],
                                                  **kw_yy)
                        #plot phase                         
                        erxx= mtplottools.plot_errorbar(axpxx, 
                                                  period[nzxx], 
                                                  rp.phasexx[nzxx], 
                                                  rp.phasexx_err[nzxx],
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzxy], 
                                                  rp.phasexy[nzxy], 
                                                  rp.phasexy_err[nzxy],
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx, 
                                                  period[nzyx], 
                                                  rp.phaseyx[nzyx], 
                                                  rp.phaseyx_err[nzyx],
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpyy, 
                                                  period[nzyy], 
                                                  rp.phaseyy[nzyy], 
                                                  rp.phaseyy_err[nzyy],
                                                  **kw_yy)
                    elif self.plot_z == True:
                        #plot real
                        erxx = mtplottools.plot_errorbar(axrxx, 
                                                  period[nzxx], 
                                                  z_obj.z[nzxx,0,0].real, 
                                                  z_obj.zerr[nzxx,0,0].real,
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].real, 
                                                  z_obj.zerr[nzxy,0,1].real,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axryx, 
                                                  period[nzyx], 
                                                  z_obj.z[nzyx,1,0].real, 
                                                  z_obj.zerr[nzyx,1,0].real,
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axryy, 
                                                  period[nzyy], 
                                                  z_obj.z[nzyy,1,1].real, 
                                                  z_obj.zerr[nzyy,1,1].real,
                                                  **kw_yy)
                        #plot phase                         
                        erxx = mtplottools.plot_errorbar(axpxx, 
                                                  period[nzxx], 
                                                  z_obj.z[nzxx,0,0].imag, 
                                                  z_obj.zerr[nzxx,0,0].imag,
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].imag, 
                                                  z_obj.zerr[nzxy,0,1].imag,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx, 
                                                  period[nzyx], 
                                                  z_obj.z[nzyx,1,0].imag, 
                                                  z_obj.zerr[nzyx,1,0].imag,
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpyy, 
                                                  period[nzyy], 
                                                  z_obj.z[nzyy,1,1].imag, 
                                                  z_obj.zerr[nzyy,1,1].imag,
                                                  **kw_yy)
                                                  
                    #plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtxr, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].real,
                                                 t_obj.tippererr[ntx, 0, 0],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axtyr, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].real,
                                                 t_obj.tippererr[nty, 0, 0],
                                                 **kw_yy)
                                                 
                        ertx = mtplottools.plot_errorbar(axtxi, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].imag,
                                                 t_obj.tippererr[ntx, 0, 1],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axtyi, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].imag,
                                                 t_obj.tippererr[nty, 0, 1],
                                                 **kw_yy)
                    if plot_tipper == False:                    
                        ax_list = [axrxx, axrxy, axryx, axryy, 
                                   axpxx, axpxy, axpyx, axpyy]
                        line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
                        label_list = [['$Z_{xx}$'], ['$Z_{xy}$'], 
                                      ['$Z_{yx}$'], ['$Z_{yy}$']]
                    else:                    
                        ax_list = [axrxx, axrxy, axryx, axryy, 
                                   axpxx, axpxy, axpyx, axpyy, 
                                   axtxr, axtxi, axtyr, axtyi]
                        line_list = [[erxx[0]], [erxy[0]], 
                                     [eryx[0]], [eryy[0]],
                                     [ertx[0]], [erty[0]]]
                        label_list = [['$Z_{xx}$'], ['$Z_{xy}$'], 
                                      ['$Z_{yx}$'], ['$Z_{yy}$'],
                                      ['$T_{x}$'], ['$T_{y}$']]
                    
                #set axis properties
                for aa, ax in enumerate(ax_list):
                    ax.tick_params(axis='y', pad=self.ylabel_pad)
                    ylabels = ax.get_yticks().tolist()
                    ylabels[-1] = ''
                    ylabels[0] = ''
                    ax.set_yticklabels(ylabels)
                    
                    if len(ax_list) == 4 or len(ax_list) == 6:
                        if aa < 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                            
                        #set axes labels
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Z (mV/km nT)]',
                                              fontdict=fontdict)
                        elif aa == 2:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Z (mV/km nT)]',
                                              fontdict=fontdict)
                            
                    elif len(ax_list) == 8 or len(ax_list) == 12:
                        if aa < 4:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            if aa == 8 or aa == 10:
                                plt.setp(ax.get_xticklabels(), visible=False)
                            else:
                                ax.set_ylim(self.phase_limits)
                                ax.set_xlabel('Period (s)', fontdict=fontdict)

                        #set axes labels
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Z (mV/km nT)]',
                                               fontdict=fontdict)
                        elif aa == 4:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Z (mV/km nT)]',
                                              fontdict=fontdict)

                    ax.set_xscale('log')
                    ax.set_xlim(xmin=10**(np.floor(np.log10(period[0])))*1.01,
                             xmax=10**(np.ceil(np.log10(period[-1])))*.99)
                    ax.grid(True, alpha=.25)
                    
            # plot xy and yx together and xx, yy together
            elif self.plot_style == 2:
                if self.plot_component == 2:
                    if plot_tipper == False:
                        axrxy = fig.add_subplot(gs[0, 0:])
                        axpxy = fig.add_subplot(gs[1, 0:], sharex=axrxy)
                    else:
                        axrxy = fig.add_subplot(gs[0, 0:4])
                        axpxy = fig.add_subplot(gs[1, 0:4], sharex=axrxy)
                        axtr = fig.add_subplot(gs[0, 4:], sharex=axrxy)
                        axti = fig.add_subplot(gs[1, 4:], sharex=axrxy)
                        
                    if self.plot_z == False:
                        #plot resistivity
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  rp.resxy[nzxy], 
                                                  rp.resxy_err[nzxy],
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzyx], 
                                                  rp.resyx[nzyx], 
                                                  rp.resyx_err[nzyx],
                                                  **kw_yy)
                        #plot phase                         
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzxy], 
                                                  rp.phasexy[nzxy], 
                                                  rp.phasexy_err[nzxy],
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzyx], 
                                                  rp.phaseyx[nzyx], 
                                                  rp.phaseyx_err[nzyx],
                                                  **kw_yy)
                    elif self.plot_z == True:
                        #plot real
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].real, 
                                                  z_obj.zerr[nzxy,0,1].real,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,1,0].real, 
                                                  z_obj.zerr[nzxy,1,0].real,
                                                  **kw_yy)
                        #plot phase                         
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].imag, 
                                                  z_obj.zerr[nzxy,0,1].imag,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzyx], 
                                                  z_obj.z[nzyx,1,0].imag, 
                                                  z_obj.zerr[nzyx,1,0].imag,
                                                  **kw_yy)
                    #plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].real,
                                                 t_obj.tippererr[ntx, 0, 0],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].real,
                                                 t_obj.tippererr[nty, 0, 1],
                                                 **kw_yy)
                                                 
                        ertx = mtplottools.plot_errorbar(axti, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].imag,
                                                 t_obj.tippererr[ntx, 0, 0],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axti, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].imag,
                                                 t_obj.tippererr[nty, 0, 1],
                                                 **kw_yy)
                    
                    if plot_tipper == False:    
                        ax_list = [axrxy, axpxy]
                        line_list = [erxy[0], eryx[0]]
                        label_list = ['$Z_{xy}$', '$Z_{yx}$']
                    else:    
                        ax_list = [axrxy, axpxy, axtr, axti]
                        line_list = [[erxy[0], eryx[0]], 
                                     [ertx[0], erty[0]]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'],
                                      ['$T_{x}$', '$T_{y}$']]
                    
                elif self.plot_component == 4:
                    if plot_tipper == False:
                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)
                        
                        axrxx = fig.add_subplot(gs[0, 2:], sharex=axrxy)
                        axpxx = fig.add_subplot(gs[1, 2:], sharex=axrxy)
                    else:
                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)
                        
                        axrxx = fig.add_subplot(gs[0, 2:4], sharex=axrxy)
                        axpxx = fig.add_subplot(gs[1, 2:4], sharex=axrxy)
                        
                        axtr = fig.add_subplot(gs[0, 4:], sharex=axrxy)
                        axti = fig.add_subplot(gs[1, 4:], sharex=axrxy)
                        
                    if self.plot_z == False:
                        #plot resistivity
                        erxx= mtplottools.plot_errorbar(axrxx, 
                                                  period[nzxx], 
                                                  rp.resxx[nzxx], 
                                                  rp.resxx_err[nzxx],
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  rp.resxy[nzxy], 
                                                  rp.resxy_err[nzxy],
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzyx], 
                                                  rp.resyx[nzyx], 
                                                  rp.resyx_err[nzyx],
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axrxx, 
                                                  period[nzyy], 
                                                  rp.resyy[nzyy], 
                                                  rp.resyy_err[nzyy],
                                                  **kw_yy)
                        #plot phase                         
                        erxx= mtplottools.plot_errorbar(axpxx, 
                                                  period[nzxx], 
                                                  rp.phasexx[nzxx], 
                                                  rp.phasexx_err[nzxx],
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzxy], 
                                                  rp.phasexy[nzxy], 
                                                  rp.phasexy_err[nzxy],
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzyx], 
                                                  rp.phaseyx[nzyx], 
                                                  rp.phaseyx_err[nzyx],
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpxx, 
                                                  period[nzyy], 
                                                  rp.phaseyy[nzyy], 
                                                  rp.phaseyy_err[nzyy],
                                                  **kw_yy)
                    elif self.plot_z == True:
                         #plot real
                        erxx = mtplottools.plot_errorbar(axrxx, 
                                                  period[nzxx], 
                                                  z_obj.z[nzxx,0,0].real, 
                                                  z_obj.zerr[nzxx,0,0].real,
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].real, 
                                                  z_obj.zerr[nzxy,0,1].real,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy, 
                                                  period[nzyx], 
                                                  z_obj.z[nzyx,1,0].real, 
                                                  z_obj.zerr[nzyx,1,0].real,
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axrxx, 
                                                  period[nzyy], 
                                                  z_obj.z[nzyy,1,1].real, 
                                                  z_obj.zerr[nzyy,1,1].real,
                                                  **kw_yy)
                        #plot phase                         
                        erxx = mtplottools.plot_errorbar(axpxx, 
                                                  period[nzxx], 
                                                  z_obj.z[nzxx,0,0].imag, 
                                                  z_obj.zerr[nzxx,0,0].imag,
                                                  **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzxy], 
                                                  z_obj.z[nzxy,0,1].imag, 
                                                  z_obj.zerr[nzxy,0,1].imag,
                                                  **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy, 
                                                  period[nzyx], 
                                                  z_obj.z[nzyx,1,0].imag, 
                                                  z_obj.zerr[nzyx,1,0].imag,
                                                  **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpxx, 
                                                  period[nzyy], 
                                                  z_obj.z[nzyy,1,1].imag, 
                                                  z_obj.zerr[nzyy,1,1].imag,
                                                  **kw_yy)
                    #plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].real,
                                                 t_obj.tippererr[ntx, 0, 0],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].real,
                                                 t_obj.tippererr[nty, 0, 1],
                                                 **kw_yy)
                                                 
                        ertx = mtplottools.plot_errorbar(axti, 
                                                 period,
                                                 t_obj.tipper[ntx, 0, 0].imag,
                                                 t_obj.tippererr[ntx, 0, 0],
                                                 **kw_xx)
                        erty = mtplottools.plot_errorbar(axti, 
                                                 period,
                                                 t_obj.tipper[nty, 0, 1].imag,
                                                 t_obj.tippererr[nty, 0, 1],
                                                 **kw_yy)
                    
                    if plot_tipper == False:
                        ax_list = [axrxy, axrxx, axpxy, axpxx]
                        line_list = [[erxy[0], eryx[0]], [erxx[0], eryy[0]]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'], 
                                      ['$Z_{xx}$', '$Z_{yy}$']]
                    else:
                        ax_list = [axrxy, axrxx, axpxy, axpxx, axtr, axti]
                        line_list = [[erxy[0], eryx[0]], [erxx[0], eryy[0]],
                                     [ertx[0]], erty[0]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'], 
                                      ['$Z_{xx}$', '$Z_{yy}$'],
                                      ['$T_x$', '$T_y$']]
                        
                #set axis properties
                for aa, ax in enumerate(ax_list):
                    ax.tick_params(axis='y', pad=self.ylabel_pad)
                    ylabels = ax.get_yticks().tolist()
                    ylabels[-1] = ''
                    ylabels[0] = ''
                    ax.set_yticklabels(ylabels)
                    if len(ax_list) == 2:
                        ax.set_xlabel('Period (s)', fontdict=fontdict)
                        if aa == 0:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Impedance (m/s)]',
                                               fontdict=fontdict)
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Impedance (m/s)]',
                                              fontdict=fontdict)
                    elif len(ax_list) == 4 and plot_tipper == False:
                        if aa < 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            if self.plot_z == False:
                                ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Z (mV/km nT)]',
                                               fontdict=fontdict)
                        elif aa == 2:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Z (mV/km nT)]',
                                              fontdict=fontdict)
                    elif len(ax_list) == 4 and plot_tipper == True:
                        if aa == 0 or aa == 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Z (mV/km nT)]',
                                               fontdict=fontdict)
                        elif aa == 1:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Z (mV/km nT)]',
                                              fontdict=fontdict)
#                        else:
#                            plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    ax.set_xscale('log')
                    ax.set_xlim(xmin=10**(np.floor(np.log10(period[0])))*1.01,
                                xmax=10**(np.ceil(np.log10(period[-1])))*.99)
                    ax.grid(True,alpha=.25)

            if plotr == True:
                for rr in range(nr):
                    if self.color_mode == 'color':   
                        cxy = (0,.4+float(rr)/(3*nr),0)
                        cyx = (.7+float(rr)/(4*nr),.13,.63-float(rr)/(4*nr))
                    elif self.color_mode == 'bw':
                        cxy = (1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))                    
                        cyx = (1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))
                    
                    resp_z_obj = self.resp_object[rr].mt_dict[station].Z
                    resp_z_err = (z_obj.z-resp_z_obj.z)/z_obj.zerr
    
                    resp_t_obj = self.resp_object[rr].mt_dict[station].Tipper
                    
                    rrp = mtplottools.ResPhase(resp_z_obj)
    
                    rms = resp_z_err.std()
                    rms_xx = resp_z_err[:, 0, 0].std()
                    rms_xy = resp_z_err[:, 0, 1].std()
                    rms_yx = resp_z_err[:, 1, 0].std()
                    rms_yy = resp_z_err[:, 1, 1].std()
                    print ' --- response {0} ---'.format(rr)
                    print '  RMS = {:.2f}'.format(rms)
                    print '      RMS_xx = {:.2f}'.format(rms_xx)
                    print '      RMS_xy = {:.2f}'.format(rms_xy)
                    print '      RMS_yx = {:.2f}'.format(rms_yx)
                    print '      RMS_yy = {:.2f}'.format(rms_yy)
                    
                    #--> make key word dictionaries for plotting
                    kw_xx = {'color':cxy,
                             'marker':self.mted,
                             'ms':self.ms,
                             'ls':':',
                             'lw':self.lw,
                             'e_capsize':self.e_capsize,
                             'e_capthick':self.e_capthick}        
                   
                    kw_yy = {'color':cyx,
                             'marker':self.mtmd,
                             'ms':self.ms,
                             'ls':':',
                             'lw':self.lw,
                             'e_capsize':self.e_capsize,
                             'e_capthick':self.e_capthick}
                             
                    if self.plot_style == 1:
                        if self.plot_component == 2:
                            if self.plot_z == False:
                                #plot resistivity
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          rrp.resxy[nzxy], 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axryx, 
                                                          period[nzyx], 
                                                          rrp.resyx[nzyx], 
                                                          **kw_yy)
                                #plot phase                         
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          rrp.phasexy[nzxy], 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpyx, 
                                                          period[nzyx], 
                                                          rrp.phaseyx[nzyx], 
                                                          **kw_yy)
                            elif self.plot_z == True:
                                #plot real
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].real, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axryx, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].real, 
                                                          **kw_yy)
                                #plot phase                         
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].imag, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpyx, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].imag, 
                                                          **kw_yy)
                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].real,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].real,
                                             **kw_yy)
                                                         
                                rertx = mtplottools.plot_errorbar(axti, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].imag,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].imag,
                                             **kw_yy)
                            if plot_tipper == False:
                                line_list[0] += [rerxy[0]]
                                line_list[1] += [reryx[0]]
                                label_list[0] += ['$Z^m_{xy}$ '+
                                                   'rms={0:.2f}'.format(rms_xy)]
                                label_list[1] += ['$Z^m_{yx}$ '+
                                               'rms={0:.2f}'.format(rms_yx)]
                            else:
                                line_list[0] += [rerxy[0]]
                                line_list[1] += [reryx[0]]
                                line_list[2] += [rertx[0], rerty[0]]
                                label_list[0] += ['$Z^m_{xy}$ '+
                                                   'rms={0:.2f}'.format(rms_xy)]
                                label_list[1] += ['$Z^m_{yx}$ '+
                                               'rms={0:.2f}'.format(rms_yx)]
                                label_list[2] += ['$T^m_{x}$', '$T^m_{y}$']
                        elif self.plot_component == 4:
                            if self.plot_z == False:
                                #plot resistivity
                                rerxx= mtplottools.plot_errorbar(axrxx, 
                                                          period[nzxx], 
                                                          rrp.resxx[nzxx], 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          rrp.resxy[nzxy], 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axryx, 
                                                          period[nzyx], 
                                                          rrp.resyx[nzyx], 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axryy, 
                                                          period[nzyy], 
                                                          rrp.resyy[nzyy], 
                                                          **kw_yy)
                                #plot phase                         
                                rerxx= mtplottools.plot_errorbar(axpxx, 
                                                          period[nzxx], 
                                                          rrp.phasexx[nzxx], 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          rrp.phasexy[nzxy], 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpyx, 
                                                          period[nzyx], 
                                                          rrp.phaseyx[nzyx], 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpyy, 
                                                          period[nzyy], 
                                                          rrp.phaseyy[nzyy], 
                                                          **kw_yy)
                            elif self.plot_z == True:
                                #plot real
                                rerxx = mtplottools.plot_errorbar(axrxx, 
                                                          period[nzxx], 
                                                          resp_z_obj.z[nzxx,0,0].real, 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].real, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axryx, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].real, 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axryy, 
                                                          period[nzyy], 
                                                          resp_z_obj.z[nzyy,1,1].real, 
                                                          **kw_yy)
                                #plot phase                         
                                rerxx = mtplottools.plot_errorbar(axpxx, 
                                                          period[nzxx], 
                                                          resp_z_obj.z[nzxx,0,0].imag, 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].imag, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpyx, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].imag, 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpyy, 
                                                          period[nzyy], 
                                                          resp_z_obj.z[nzyy,1,1].imag, 
                                                          **kw_yy)
                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtxr, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].real,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtyr, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].real,
                                             **kw_yy)
                                                         
                                rertx = mtplottools.plot_errorbar(axtxi, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].imag,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtyi, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].imag,
                                             **kw_yy)
                                             
                            if plot_tipper == False:
                                line_list[0] += [rerxx[0]]
                                line_list[1] += [rerxy[0]]
                                line_list[2] += [reryx[0]]
                                line_list[3] += [reryy[0]]
                                label_list[0] += ['$Z^m_{xx}$ '+
                                                   'rms={0:.2f}'.format(rms_xx)]
                                label_list[1] += ['$Z^m_{xy}$ '+
                                               'rms={0:.2f}'.format(rms_xy)]
                                label_list[2] += ['$Z^m_{yx}$ '+
                                               'rms={0:.2f}'.format(rms_yx)]
                                label_list[3] += ['$Z^m_{yy}$ '+
                                               'rms={0:.2f}'.format(rms_yy)]
                            else:
                                line_list[0] += [rerxx[0]]
                                line_list[1] += [rerxy[0]]
                                line_list[2] += [reryx[0]]
                                line_list[3] += [reryy[0]]
                                line_list[4] += [rertx[0]]
                                line_list[5] += [rerty[0]]
                                label_list[0] += ['$Z^m_{xx}$ '+
                                                   'rms={0:.2f}'.format(rms_xx)]
                                label_list[1] += ['$Z^m_{xy}$ '+
                                               'rms={0:.2f}'.format(rms_xy)]
                                label_list[2] += ['$Z^m_{yx}$ '+
                                               'rms={0:.2f}'.format(rms_yx)]
                                label_list[3] += ['$Z^m_{yy}$ '+
                                               'rms={0:.2f}'.format(rms_yy)]
                                label_list[4] += ['$T^m_{x}$ ']
                                label_list[5] += ['$T^m_{y}$']
                                           
                    elif self.plot_style == 2:
                        if self.plot_component == 2:
                            if self.plot_z == False:                            
                                #plot resistivity
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          rrp.resxy[nzxy], 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzyx], 
                                                          rrp.resyx[nzyx], 
                                                          **kw_yy)
                                #plot phase                         
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          rrp.phasexy[nzxy], 
                                                         **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzyx], 
                                                          rrp.phaseyx[nzyx], 
                                                          **kw_yy)
                            elif self.plot_z == True:
                                #plot real
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].real, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].real, 
                                                          **kw_yy)
                                #plot phase                         
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].imag, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].imag, 
                                                          **kw_xx)
                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].real,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].real,
                                             **kw_yy)
                                                         
                                rertx = mtplottools.plot_errorbar(axti, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].imag,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].imag,
                                             **kw_yy)
                                
                            if plot_tipper == False:
                                line_list += [rerxy[0], reryx[0]]
                                label_list += ['$Z^m_{xy}$ '+
                                               'rms={0:.2f}'.format(rms_xy),
                                               '$Z^m_{yx}$ '+
                                               'rms={0:.2f}'.format(rms_yx)]
                            else:
                                line_list[0] += [rerxy[0], reryx[0]]
                                line_list[1] += [rertx[0], rerty[0]]
                                label_list[0] += ['$Z^m_{xy}$ '+
                                               'rms={0:.2f}'.format(rms_xy),
                                               '$Z^m_{yx}$ '+
                                               'rms={0:.2f}'.format(rms_yx)]
                                label_list[1] += ['$T^m_{x}$ ', '$T^m_{y}$']
                                
                        elif self.plot_component == 4:
                            if self.plot_z == False:
                                #plot resistivity
                                rerxx= mtplottools.plot_errorbar(axrxx, 
                                                          period[nzxx], 
                                                          rrp.resxx[nzxx], 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          rrp.resxy[nzxy], 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzyx], 
                                                          rrp.resyx[nzyx], 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axrxx, 
                                                          period[nzyy], 
                                                          rrp.resyy[nzyy], 
                                                          **kw_yy)
                                #plot phase                         
                                rerxx= mtplottools.plot_errorbar(axpxx, 
                                                          period[nzxx], 
                                                          rrp.phasexx[nzxx], 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          rrp.phasexy[nzxy], 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzyx], 
                                                          rrp.phaseyx[nzyx], 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpxx, 
                                                          period[nzyy], 
                                                          rrp.phaseyy[nzyy], 
                                                          **kw_yy)
                            elif self.plot_z == True:
                                #plot real
                                rerxx = mtplottools.plot_errorbar(axrxx, 
                                                          period[nzxx], 
                                                          resp_z_obj.z[nzxx,0,0].real, 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].real, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].real, 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axrxx, 
                                                          period[nzyy], 
                                                          resp_z_obj.z[nzyy,1,1].real, 
                                                          **kw_yy)
                                #plot phase                         
                                rerxx = mtplottools.plot_errorbar(axpxx, 
                                                          period[nzxx], 
                                                          resp_z_obj.z[nzxx,0,0].imag, 
                                                          **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzxy], 
                                                          resp_z_obj.z[nzxy,0,1].imag, 
                                                          **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy, 
                                                          period[nzyx], 
                                                          resp_z_obj.z[nzyx,1,0].imag, 
                                                          **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpxx, 
                                                          period[nzyy], 
                                                          resp_z_obj.z[nzyy,1,1].imag, 
                                                          **kw_yy)
                                                          
                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].real,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].real,
                                             **kw_yy)
                                                         
                                rertx = mtplottools.plot_errorbar(axti, 
                                             period,
                                             resp_t_obj.tipper[ntx, 0, 0].imag,
                                             **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti, 
                                             period,
                                             resp_t_obj.tipper[nty, 0, 1].imag,
                                             **kw_yy)
                                             
                            if plot_tipper == False:
                                line_list[0] += [rerxy[0], reryx[0]]
                                line_list[1] += [rerxx[0], reryy[0]]
                                label_list[0] += ['$Z^m_{xy}$ '+
                                                   'rms={0:.2f}'.format(rms_xy),
                                                  '$Z^m_{yx}$ '+
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[1] += ['$Z^m_{xx}$ '+
                                                   'rms={0:.2f}'.format(rms_xx),
                                                  '$Z^m_{yy}$ '+
                                                  'rms={0:.2f}'.format(rms_yy)]
                            else:
                                line_list[0] += [rerxy[0], reryx[0]]
                                line_list[1] += [rerxx[0], reryy[0]]
                                line_list[2] += [rertx[0], rerty[0]]
                                label_list[0] += ['$Z^m_{xy}$ '+
                                                   'rms={0:.2f}'.format(rms_xy),
                                                  '$Z^m_{yx}$ '+
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[1] += ['$Z^m_{xx}$ '+
                                                   'rms={0:.2f}'.format(rms_xx),
                                                  '$Z^m_{yy}$ '+
                                                  'rms={0:.2f}'.format(rms_yy)]
                                label_list[2] += ['$T^m_{x}$ ', '$T^m_{y}$']
                    
            #make legends
            if self.plot_style == 1:
                legend_ax_list = ax_list[0:self.plot_component]
                if plot_tipper == True:
                    if self.plot_component == 2:
                        legend_ax_list.append(ax_list[4])
                    elif self.plot_component == 4:
                        legend_ax_list.append(ax_list[8])
                        legend_ax_list.append(ax_list[10])
                for aa, ax in enumerate(legend_ax_list):
                    ax.legend(line_list[aa],
                              label_list[aa],
                              loc=self.legend_loc,
                              markerscale=self.legend_marker_scale,
                              borderaxespad=self.legend_border_axes_pad,
                              labelspacing=self.legend_label_spacing,
                              handletextpad=self.legend_handle_text_pad,
                              borderpad=self.legend_border_pad,
                              prop={'size':max([self.font_size/(nr+1), 5])})
            if self.plot_style == 2:
                if self.plot_component == 2:
                    legend_ax_list = [ax_list[0]]
                    if plot_tipper == True:
                        legend_ax_list.append(ax_list[2])
                    for aa, ax in enumerate(legend_ax_list):
                        ax.legend(line_list[aa],
                                  label_list[aa],
                                  loc=self.legend_loc,
                                  markerscale=self.legend_marker_scale,
                                  borderaxespad=self.legend_border_axes_pad,
                                  labelspacing=self.legend_label_spacing,
                                  handletextpad=self.legend_handle_text_pad,
                                  borderpad=self.legend_border_pad,
                                  prop={'size':max([self.font_size/(nr+1), 5])})
                else:
                    legend_ax_list = ax_list[0:self.plot_component/2]
                    if plot_tipper == True:
                        if self.plot_component == 2:
                            legend_ax_list.append(ax_list[2])
                        elif self.plot_component == 4:
                            legend_ax_list.append(ax_list[4])
                    for aa, ax in enumerate(legend_ax_list):
                        ax.legend(line_list[aa],
                                  label_list[aa],
                                  loc=self.legend_loc,
                                  markerscale=self.legend_marker_scale,
                                  borderaxespad=self.legend_border_axes_pad,
                                  labelspacing=self.legend_label_spacing,
                                  handletextpad=self.legend_handle_text_pad,
                                  borderpad=self.legend_border_pad,
                                  prop={'size':max([self.font_size/(nr+1), 5])})
        
        ##--> BE SURE TO SHOW THE PLOT
        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()
        
    def save_figure(self, save_fn, file_format='pdf', orientation='portrait', 
                  fig_dpi=None, close_fig='y'):
        """
        save_plot will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')
            
        else:
            save_fn = os.path.join(save_fn, '_L2.'+
                                    file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')
        
        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)
        
        else:
            pass
        
        self.fig_fn = save_fn
        print 'Saved figure to: '+self.fig_fn
        
    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        self.fig.canvas.draw()
                          
    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """
        
        return ("Plots data vs model response computed by WS3DINV")

#==============================================================================
# plot phase tensors
#==============================================================================
class PlotPTMaps(mtplottools.MTEllipse):
    """
    Plot phase tensor maps including residual pt if response file is input.
    
    :Plot only data for one period: ::
    
        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, plot_period_list=[0])
        
    :Plot data and model response: ::
    
        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> rfn = r"/home/MT/ws3dinv/Inv1/Test_resp.00"
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, resp_fn=rfn, model_fn=mfn,
        >>> ...                 plot_period_list=[0])
        >>> # adjust colorbar
        >>> ptm.cb_res_pad = 1.25
        >>> ptm.redraw_plot()
    
    
    ========================== ================================================
    Attributes                 Description    
    ========================== ================================================
    cb_pt_pad                  percentage from top of axes to place pt 
                               color bar. *default* is .90
    cb_res_pad                 percentage from bottom of axes to place
                               resistivity color bar. *default* is 1.2
    cb_residual_tick_step      tick step for residual pt. *default* is 3
    cb_tick_step               tick step for phase tensor color bar, 
                               *default* is 45
    data                       np.ndarray(n_station, n_periods, 2, 2)
                               impedance tensors for station data                     
    data_fn                    full path to data fle               
    dscale                     scaling parameter depending on map_scale
    ellipse_cmap               color map for pt ellipses. *default* is
                               mt_bl2gr2rd
    ellipse_colorby            [ 'skew' | 'skew_seg' | 'phimin' | 'phimax'|
                                 'phidet' | 'ellipticity' ] parameter to color
                                 ellipses by. *default* is 'phimin'
    ellipse_range              (min, max, step) min and max of colormap, need
                               to input step if plotting skew_seg
    ellipse_size               relative size of ellipses in map_scale
    ew_limits                  limits of plot in e-w direction in map_scale
                               units.  *default* is None, scales to station 
                               area
    fig_aspect                 aspect of figure. *default* is 1
    fig_dpi                    resolution in dots-per-inch. *default* is 300
    fig_list                   list of matplotlib.figure instances for each
                               figure plotted.
    fig_size                   [width, height] in inches of figure window
                               *default* is [6, 6]
    font_size                  font size of ticklabels, axes labels are 
                               font_size+2. *default* is 7
    grid_east                  relative location of grid nodes in e-w direction
                               in map_scale units
    grid_north                 relative location of grid nodes in n-s direction
                               in map_scale units
    grid_z                     relative location of grid nodes in z direction
                               in map_scale units
    model_fn                 full path to initial file
    map_scale                  [ 'km' | 'm' ] distance units of map. 
                               *default* is km
    mesh_east                  np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north                 np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                   full path to model file
    nodes_east                 relative distance betwen nodes in e-w direction
                               in map_scale units
    nodes_north                relative distance betwen nodes in n-s direction
                               in map_scale units
    nodes_z                    relative distance betwen nodes in z direction
                               in map_scale units
    ns_limits                  (min, max) limits of plot in n-s direction
                               *default* is None, viewing area is station area
    pad_east                   padding from extreme stations in east direction
    pad_north                  padding from extreme stations in north direction
    period_list                list of periods from data
    plot_grid                  [ 'y' | 'n' ] 'y' to plot grid lines
                               *default* is 'n'
    plot_period_list           list of period index values to plot
                               *default* is None 
    plot_yn                    ['y' | 'n' ] 'y' to plot on instantiation
                               *default* is 'y'
    res_cmap                   colormap for resisitivity values. 
                               *default* is 'jet_r'
    res_limits                 (min, max) resistivity limits in log scale
                               *default* is (0, 4)
    res_model                  np.ndarray(n_north, n_east, n_vertical) of 
                               model resistivity values in linear scale
    residual_cmap              color map for pt residuals. 
                               *default* is 'mt_wh2or' 
    resp                       np.ndarray(n_stations, n_periods, 2, 2)
                               impedance tensors for model response  
    resp_fn                    full path to response file
    save_path                  directory to save figures to
    save_plots                 [ 'y' | 'n' ] 'y' to save plots to save_path
    station_east               location of stations in east direction in 
                               map_scale units  
    station_fn                 full path to station locations file
    station_names              station names
    station_north              location of station in north direction in 
                               map_scale units
    subplot_bottom             distance between axes and bottom of figure window
    subplot_left               distance between axes and left of figure window  
    subplot_right              distance between axes and right of figure window
    subplot_top                distance between axes and top of figure window
    title                      titiel of plot *default* is depth of slice
    xminorticks                location of xminorticks
    yminorticks                location of yminorticks
    ========================== ================================================
    """
    
    def __init__(self, data_fn=None, resp_fn=None, model_fn=None, **kwargs):
        
        self.model_fn = model_fn
        self.data_fn = data_fn
        self.resp_fn = resp_fn
        
        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
            
        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)
                
        self.save_plots = kwargs.pop('save_plots', 'y')
        self.plot_period_list = kwargs.pop('plot_period_list', None)
        self.period_dict = None
        
        self.map_scale = kwargs.pop('map_scale', 'km')
        #make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
        elif self.map_scale == 'm':
            self.dscale = 1. 
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)
        
        self.pad_east = kwargs.pop('pad_east', 2000)
        self.pad_north = kwargs.pop('pad_north', 2000)
        
        self.plot_grid = kwargs.pop('plot_grid', 'n')
        
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.fig_list = []
        
        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)
        
        self.residual_cmap = kwargs.pop('residual_cmap', 'mt_wh2or')
        self.font_size = kwargs.pop('font_size', 7)
        
        self.cb_tick_step = kwargs.pop('cb_tick_step', 45)
        self.cb_residual_tick_step = kwargs.pop('cb_residual_tick_step', 3)
        self.cb_pt_pad = kwargs.pop('cb_pt_pad', 1.2)
        self.cb_res_pad = kwargs.pop('cb_res_pad', .5)
        
        
        self.res_limits = kwargs.pop('res_limits', (0,4))
        self.res_cmap = kwargs.pop('res_cmap', 'jet_r')
        
        #--> set the ellipse properties -------------------
        self._ellipse_dict = kwargs.pop('ellipse_dict', {'size':2})
        self._read_ellipse_dict()
        
        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .2
        self.subplot_wspace = .05
        
        self.data_obj = None
        self.resp_obj = None
        self.model_obj = None
        self.period_list = None
        
        self.pt_data_arr = None
        self.pt_resp_arr = None
        self.pt_resid_arr = None
        
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()
            
    def _read_files(self):
        """
        get information from files
        """

        #--> read in data file 
        self.data_obj = Data()
        self.data_obj.read_data_file(self.data_fn)
        
        #--> read response file
        if self.resp_fn is not None:
            self.resp_obj = Data()
            self.resp_obj.read_data_file(self.resp_fn)
            
        #--> read mode file
        if self.model_fn is not None:
            self.model_obj = Model()
            self.model_obj.read_model_file(self.model_fn)
            
        self._get_plot_period_list()
        self._get_pt()
            
    def _get_plot_period_list(self):
        """
        get periods to plot from input or data file
        """
        #--> get period list to plot
        if self.plot_period_list is None:
            self.plot_period_list = self.data_obj.period_list
        else:
            if type(self.plot_period_list) is list:
                #check if entries are index values or actual periods
                if type(self.plot_period_list[0]) is int:
                    self.plot_period_list = [self.period_list[ii]
                                             for ii in self.plot_period_list]
                else:
                    pass
            elif type(self.plot_period_list) is int:
                self.plot_period_list = self.period_list[self.plot_period_list]
            elif type(self.plot_period_list) is float:
                self.plot_period_list = [self.plot_period_list]
                
        self.period_dict = dict([(key, value) for value, key in 
                                 enumerate(self.data_obj.period_list)])
                                 
    def _get_pt(self):
        """
        put pt parameters into something useful for plotting
        """
        
        ns = len(self.data_obj.mt_dict.keys())
        nf = len(self.data_obj.period_list)
        
        data_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                ('phimax', np.float),
                                                ('skew', np.float),
                                                ('azimuth', np.float),
                                                ('east', np.float),
                                                ('north', np.float)])
        if self.resp_fn is not None:
            model_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                    ('phimax', np.float),
                                                    ('skew', np.float),
                                                    ('azimuth', np.float),
                                                    ('east', np.float),
                                                    ('north', np.float)])
            
            res_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                    ('phimax', np.float),
                                                    ('skew', np.float),
                                                    ('azimuth', np.float),
                                                    ('east', np.float),
                                                    ('north', np.float),
                                                    ('geometric_mean', np.float)])
                                                
        for ii, key in enumerate(self.data_obj.mt_dict.keys()):
            east = self.data_obj.mt_dict[key].grid_east/self.dscale
            north = self.data_obj.mt_dict[key].grid_north/self.dscale
            dpt = self.data_obj.mt_dict[key].pt
            data_pt_arr[:, ii]['east'] = east
            data_pt_arr[:, ii]['north'] = north
            data_pt_arr[:, ii]['phimin'] = dpt.phimin[0]
            data_pt_arr[:, ii]['phimax'] = dpt.phimax[0]
            data_pt_arr[:, ii]['azimuth'] = dpt.azimuth[0]
            data_pt_arr[:, ii]['skew'] = dpt.beta[0]
            if self.resp_fn is not None:
                mpt = self.resp_obj.mt_dict[key].pt
                try:
                    rpt = mtpt.ResidualPhaseTensor(pt_object1=dpt, 
                                                   pt_object2=mpt)
                    rpt = rpt.residual_pt
                    res_pt_arr[:, ii]['east'] = east
                    res_pt_arr[:, ii]['north'] = north
                    res_pt_arr[:, ii]['phimin'] = rpt.phimin[0]
                    res_pt_arr[:, ii]['phimax'] = rpt.phimax[0]
                    res_pt_arr[:, ii]['azimuth'] = rpt.azimuth[0]
                    res_pt_arr[:, ii]['skew'] = rpt.beta[0]
                    res_pt_arr[:, ii]['geometric_mean'] = np.sqrt(abs(rpt.phimin[0]*\
                                                                  rpt.phimax[0]))
                except mtex.MTpyError_PT:
                    print key, dpt.pt.shape, mpt.pt.shape
                
                model_pt_arr[:, ii]['east'] = east
                model_pt_arr[:, ii]['north'] = north
                model_pt_arr[:, ii]['phimin'] = mpt.phimin[0]
                model_pt_arr[:, ii]['phimax'] = mpt.phimax[0]
                model_pt_arr[:, ii]['azimuth'] = mpt.azimuth[0]
                model_pt_arr[:, ii]['skew'] = mpt.beta[0]
                
                
                
        #make these attributes        
        self.pt_data_arr = data_pt_arr
        if self.resp_fn is not None:
            self.pt_resp_arr = model_pt_arr
            self.pt_resid_arr = res_pt_arr
        
                    

    def plot(self):
        """
        plot phase tensor maps for data and or response, each figure is of a
        different period.  If response is input a third column is added which is 
        the residual phase tensor showing where the model is not fitting the data 
        well.  The data is plotted in km.
        
        """
        #--> read in data first
        if self.data_obj is None:        
            self._read_files()
            
        # set plot properties
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        font_dict = {'size':self.font_size+2, 'weight':'bold'}
        
        # make a grid of subplots 
        gs = gridspec.GridSpec(1, 3, hspace=self.subplot_hspace,
                               wspace=self.subplot_wspace)
                               
        
        #set some parameters for the colorbar
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                ckstep = 3
        bounds = np.arange(ckmin, ckmax+ckstep, ckstep)
        
        # set plot limits to be the station area
        if self.ew_limits == None:
            east_min = self.data_obj.coord_array['rel_east'].min()-\
                                                            self.pad_east
            east_max = self.data_obj.coord_array['rel_east'].max()+\
                                                            self.pad_east
            self.ew_limits = (east_min/self.dscale, east_max/self.dscale)
            
        if self.ns_limits == None:
            north_min = self.data_obj.coord_array['rel_north'].min()-\
                                                            self.pad_north
            north_max = self.data_obj.coord_array['rel_north'].max()+\
                                                            self.pad_north
            self.ns_limits = (north_min/self.dscale, north_max/self.dscale)

        #-------------plot phase tensors------------------------------------                    
        for ff, per in enumerate(self.plot_period_list):
            data_ii = self.period_dict[per]
            
            print 'Plotting Period: {0:.5g}'.format(per)
            fig = plt.figure('{0:.5g}'.format(per), figsize=self.fig_size,
                             dpi=self.fig_dpi)
            fig.clf()
                             
            if self.resp_fn is not None:
                axd = fig.add_subplot(gs[0, 0], aspect='equal')
                axm = fig.add_subplot(gs[0, 1], aspect='equal')
                axr = fig.add_subplot(gs[0, 2], aspect='equal')
                ax_list = [axd, axm, axr]
            
            else:
                axd = fig.add_subplot(gs[0, :], aspect='equal')
                ax_list = [axd]
            
            #plot model below the phase tensors
            if self.model_fn is not None:
                approx_depth, d_index = ws.estimate_skin_depth(self.model_obj.res_model,
                                                            self.model_obj.grid_z, 
                                                            per, 
                                                            dscale=self.dscale)  
                #need to add an extra row and column to east and north to make sure 
                #all is plotted see pcolor for details.
                plot_east = np.append(self.model_obj.grid_east, 
                                      self.model_obj.grid_east[-1]*1.25)/\
                                      self.dscale
                plot_north = np.append(self.model_obj.grid_north, 
                                       self.model_obj.grid_north[-1]*1.25)/\
                                       self.dscale
                
                #make a mesh grid for plotting
                #the 'ij' makes sure the resulting grid is in east, north
                self.mesh_east, self.mesh_north = np.meshgrid(plot_east, 
                                                              plot_north,
                                                              indexing='ij')
                
                for ax in ax_list:
                    plot_res = np.log10(self.model_obj.res_model[:, :, d_index].T)
                    ax.pcolormesh(self.mesh_east,
                                   self.mesh_north, 
                                   plot_res,
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits[0],
                                   vmax=self.res_limits[1])
                    
                
            #--> plot data phase tensors
            for pt in self.pt_data_arr[data_ii]:
                eheight = pt['phimin']/\
                          self.pt_data_arr[data_ii]['phimax'].max()*\
                          self.ellipse_size
                ewidth = pt['phimax']/\
                          self.pt_data_arr[data_ii]['phimax'].max()*\
                          self.ellipse_size
                          
                ellipse = Ellipse((pt['east'],
                                   pt['north']),
                                   width=ewidth,
                                   height=eheight,
                                   angle=90-pt['azimuth'])
                
                #get ellipse color
                if self.ellipse_cmap.find('seg')>0:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax,
                                                         bounds=bounds))
                else:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax))
                
                axd.add_artist(ellipse)
                
            #-----------plot response phase tensors---------------
            if self.resp_fn is not None:
                rcmin = np.floor(self.pt_resid_arr['geometric_mean'].min())
                rcmax = np.floor(self.pt_resid_arr['geometric_mean'].max())
                for mpt, rpt in zip(self.pt_resp_arr[data_ii], 
                                    self.pt_resid_arr[data_ii]):
                    eheight = mpt['phimin']/\
                              self.pt_resp_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                    ewidth = mpt['phimax']/\
                              self.pt_resp_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                              
                    ellipsem = Ellipse((mpt['east'],
                                       mpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=90-mpt['azimuth'])
                    
                    #get ellipse color
                    if self.ellipse_cmap.find('seg')>0:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                             self.ellipse_colorby,
                                                             self.ellipse_cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                    else:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                             self.ellipse_colorby,
                                                             self.ellipse_cmap,
                                                             ckmin,
                                                             ckmax))
                
                    axm.add_artist(ellipsem)
                    
                    #-----------plot residual phase tensors---------------
                    eheight = rpt['phimin']/\
                              self.pt_resid_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                    ewidth = rpt['phimax']/\
                              self.pt_resid_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                              
                    ellipser = Ellipse((rpt['east'],
                                       rpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=rpt['azimuth'])
                    
                    #get ellipse color
                    rpt_color = np.sqrt(abs(rpt['phimin']*rpt['phimax']))
                    if self.ellipse_cmap.find('seg')>0:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                             'geometric_mean',
                                                             self.residual_cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                    else:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                             'geometric_mean',
                                                             self.residual_cmap,
                                                             ckmin,
                                                             ckmax))
                    
                    
                    axr.add_artist(ellipser)
                
            #--> set axes properties
            # data
            axd.set_xlim(self.ew_limits)
            axd.set_ylim(self.ns_limits)
            axd.set_xlabel('Easting ({0})'.format(self.map_scale), 
                           fontdict=font_dict)
            axd.set_ylabel('Northing ({0})'.format(self.map_scale),
                           fontdict=font_dict)
            #make a colorbar for phase tensors
            #bb = axd.axes.get_position().bounds
            bb = axd.get_position().bounds
            y1 = .25*(2+(self.ns_limits[1]-self.ns_limits[0])/
                     (self.ew_limits[1]-self.ew_limits[0]))
            cb_location = (3.35*bb[2]/5+bb[0], 
                            y1*self.cb_pt_pad, .295*bb[2], .02)
            cbaxd = fig.add_axes(cb_location)
            cbd = mcb.ColorbarBase(cbaxd, 
                                   cmap=mtcl.cmapdict[self.ellipse_cmap],
                                   norm=Normalize(vmin=ckmin,
                                                  vmax=ckmax),
                                   orientation='horizontal')
            cbd.ax.xaxis.set_label_position('top')
            cbd.ax.xaxis.set_label_coords(.5, 1.75)
            cbd.set_label(mtplottools.ckdict[self.ellipse_colorby])
            cbd.set_ticks(np.arange(ckmin, ckmax+self.cb_tick_step, 
                                    self.cb_tick_step))
                                    
            axd.text(self.ew_limits[0]*.95,
                     self.ns_limits[1]*.95,
                     'Data',
                     horizontalalignment='left',
                     verticalalignment='top',
                     bbox={'facecolor':'white'},
                     fontdict={'size':self.font_size+1})
                    
            #Model and residual
            if self.resp_fn is not None:
                for aa, ax in enumerate([axm, axr]):
                    ax.set_xlim(self.ew_limits)
                    ax.set_ylim(self.ns_limits)
                    ax.set_xlabel('Easting ({0})'.format(self.map_scale), 
                                   fontdict=font_dict)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    #make a colorbar ontop of axis
                    bb = ax.axes.get_position().bounds
                    y1 = .25*(2+(self.ns_limits[1]-self.ns_limits[0])/
                     (self.ew_limits[1]-self.ew_limits[0]))
                    cb_location = (3.35*bb[2]/5+bb[0], 
                                   y1*self.cb_pt_pad, .295*bb[2], .02)
                    cbax = fig.add_axes(cb_location)
                    if aa == 0:
                        cb = mcb.ColorbarBase(cbax, 
                                              cmap=mtcl.cmapdict[self.ellipse_cmap],
                                              norm=Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                               orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(mtplottools.ckdict[self.ellipse_colorby])
                        cb.set_ticks(np.arange(ckmin, ckmax+self.cb_tick_step, 
                                    self.cb_tick_step))
                        ax.text(self.ew_limits[0]*.95,
                                self.ns_limits[1]*.95,
                                'Model',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor':'white'},
                                 fontdict={'size':self.font_size+1})
                    else:
                        cb = mcb.ColorbarBase(cbax, 
                                              cmap=mtcl.cmapdict[self.residual_cmap],
                                               norm=Normalize(vmin=rcmin,
                                                              vmax=rcmax),
                                               orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(r"$\sqrt{\Phi_{min} \Phi_{max}}$")
                        cb_ticks = [rcmin, (rcmax-rcmin)/2, rcmax]
                        cb.set_ticks(cb_ticks)
                        ax.text(self.ew_limits[0]*.95,
                                self.ns_limits[1]*.95,
                                'Residual',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor':'white'},
                                fontdict={'size':self.font_size+1})
            
            if self.model_fn is not None:
                for ax in ax_list:
                    ax.tick_params(direction='out')
                    bb = ax.axes.get_position().bounds
                    y1 = .25*(2-(self.ns_limits[1]-self.ns_limits[0])/
                             (self.ew_limits[1]-self.ew_limits[0]))
                    cb_position = (3.0*bb[2]/5+bb[0], 
                                   y1*self.cb_res_pad, .35*bb[2], .02)
                    cbax = fig.add_axes(cb_position)
                    cb = mcb.ColorbarBase(cbax, 
                                          cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits[0],
                                                         vmax=self.res_limits[1]),
                                          orientation='horizontal')
                    cb.ax.xaxis.set_label_position('top')
                    cb.ax.xaxis.set_label_coords(.5, 1.5)
                    cb.set_label('Resistivity ($\Omega \cdot$m)')
                    cb_ticks = np.arange(np.floor(self.res_limits[0]), 
                                         np.ceil(self.res_limits[1]+1), 1)
                    cb.set_ticks(cb_ticks)
                    cb.set_ticklabels([mtplottools.labeldict[ctk] for ctk in cb_ticks])
                    
                            
                
            plt.show()
            self.fig_list.append(fig)
            
    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()
            
    def save_figure(self, save_path=None, fig_dpi=None, file_format='pdf', 
                    orientation='landscape', close_fig='y'):
        """
        save_figure will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_path) == False:
            try:
                os.mkdir(save_path)
            except:
                raise IOError('Need to input a correct directory path')
        
        for fig in self.fig_list:
            per = fig.canvas.get_window_title()
            save_fn = os.path.join(save_path, 'PT_DepthSlice_{0}s.{1}'.format(
                                    per, file_format))
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')
        
            if close_fig == 'y':
                plt.close(fig)
            
            else:
                pass
        
            self.fig_fn = save_fn
            print 'Saved figure to: '+self.fig_fn
            
#==============================================================================
# plot depth slices
#==============================================================================
class PlotDepthSlice(object):
    """
    Plots depth slices of resistivity model
    
    :Example: ::
    
        >>> import mtpy.modeling.ws3dinv as ws
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> sfn = r"/home/MT/ws3dinv/Inv1/WSStationLocations.txt"
        >>> # plot just first layer to check the formating        
        >>> pds = ws.PlotDepthSlice(model_fn=mfn, station_fn=sfn, 
        >>> ...                     depth_index=0, save_plots='n')
        >>> #move color bar up 
        >>> pds.cb_location
        >>> (0.64500000000000002, 0.14999999999999997, 0.3, 0.025)
        >>> pds.cb_location = (.645, .175, .3, .025)
        >>> pds.redraw_plot()
        >>> #looks good now plot all depth slices and save them to a folder
        >>> pds.save_path = r"/home/MT/ws3dinv/Inv1/DepthSlices"
        >>> pds.depth_index = None
        >>> pds.save_plots = 'y'
        >>> pds.redraw_plot()
    
    ======================= ===================================================
    Attributes              Description    
    ======================= ===================================================
    cb_location             location of color bar (x, y, width, height)
                            *default* is None, automatically locates
    cb_orientation          [ 'vertical' | 'horizontal' ] 
                            *default* is horizontal 
    cb_pad                  padding between axes and colorbar
                            *default* is None
    cb_shrink               percentage to shrink colorbar by
                            *default* is None
    climits                 (min, max) of resistivity color on log scale
                            *default* is (0, 4)
    cmap                    name of color map *default* is 'jet_r'
    data_fn                 full path to data file
    depth_index             integer value of depth slice index, shallowest
                            layer is 0
    dscale                  scaling parameter depending on map_scale 
    ew_limits               (min, max) plot limits in e-w direction in 
                            map_scale units. *default* is None, sets viewing
                            area to the station area
    fig_aspect              aspect ratio of plot. *default* is 1
    fig_dpi                 resolution of figure in dots-per-inch. *default* is
                            300
    fig_list                list of matplotlib.figure instances for each 
                            depth slice                 
    fig_size                [width, height] in inches of figure size
                            *default* is [6, 6]
    font_size               size of ticklabel font in points, labels are 
                            font_size+2. *default* is 7
    grid_east               relative location of grid nodes in e-w direction
                            in map_scale units
    grid_north              relative location of grid nodes in n-s direction
                            in map_scale units
    grid_z                  relative location of grid nodes in z direction
                            in map_scale units
    initial_fn              full path to initial file
    map_scale               [ 'km' | 'm' ] distance units of map. *default* is 
                            km
    mesh_east               np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north              np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                full path to model file
    nodes_east              relative distance betwen nodes in e-w direction
                            in map_scale units
    nodes_north             relative distance betwen nodes in n-s direction
                            in map_scale units
    nodes_z                 relative distance betwen nodes in z direction
                            in map_scale units
    ns_limits               (min, max) plot limits in n-s direction in 
                            map_scale units. *default* is None, sets viewing
                            area to the station area
    plot_grid               [ 'y' | 'n' ] 'y' to plot mesh grid lines. 
                            *default* is 'n'
    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
    res_model               np.ndarray(n_north, n_east, n_vertical) of 
                            model resistivity values in linear scale
    save_path               path to save figures to
    save_plots              [ 'y' | 'n' ] 'y' to save depth slices to save_path
    station_east            location of stations in east direction in 
                            map_scale units  
    station_fn              full path to station locations file
    station_names           station names
    station_north           location of station in north direction in 
                            map_scale units
    subplot_bottom          distance between axes and bottom of figure window
    subplot_left            distance between axes and left of figure window  
    subplot_right           distance between axes and right of figure window
    subplot_top             distance between axes and top of figure window
    title                   titiel of plot *default* is depth of slice
    xminorticks             location of xminorticks
    yminorticks             location of yminorticks
    ======================= ===================================================
    """
    
    def __init__(self, model_fn=None, data_fn=None, **kwargs):
        self.model_fn = model_fn
        self.data_fn = data_fn

        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.initial_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.initial_fn)
            
        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)
                
        self.save_plots = kwargs.pop('save_plots', 'y')
        
        self.depth_index = kwargs.pop('depth_index', None)
        self.map_scale = kwargs.pop('map_scale', 'km')
        #make map scale
        if self.map_scale=='km':
            self.dscale=1000.
        elif self.map_scale=='m':
            self.dscale=1. 
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)
        
        self.plot_grid = kwargs.pop('plot_grid', 'n')
        
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.fig_list = []
        
        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)
        
        self.climits = kwargs.pop('climits', (0,4))
        self.cmap = kwargs.pop('cmap', 'jet_r')
        self.font_size = kwargs.pop('font_size', 8)
        
        self.cb_shrink = kwargs.pop('cb_shrink', .8)
        self.cb_pad = kwargs.pop('cb_pad', .01)
        self.cb_orientation = kwargs.pop('cb_orientation', 'horizontal')
        self.cb_location = kwargs.pop('cb_location', None)
        
        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1
        
        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z  = None
        
        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None
        
        self.mesh_east = None
        self.mesh_north = None
        
        self.station_east = None
        self.station_north = None
        self.station_names = None
        
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()
            
    def read_files(self):
        """
        read in the files to get appropriate information
        """
               #--> read in model file
        if self.model_fn is not None:
            if os.path.isfile(self.model_fn) == True:
                md_model = Model()
                md_model.read_model_file(self.model_fn)
                self.res_model = md_model.res_model
                self.grid_east = md_model.grid_east/self.dscale
                self.grid_north = md_model.grid_north/self.dscale
                self.grid_z = md_model.grid_z/self.dscale
                self.nodes_east = md_model.nodes_east/self.dscale
                self.nodes_north = md_model.nodes_north/self.dscale
                self.nodes_z = md_model.nodes_z/self.dscale
            else:
                raise mtex.MTpyError_file_handling(
                        '{0} does not exist, check path'.format(self.model_fn))
        
        #--> read in data file to get station locations
        if self.data_fn is not None:
            if os.path.isfile(self.data_fn) == True:
                md_data = Data()
                md_data.read_data_file(self.data_fn)
                self.station_east = md_data.coord_array['rel_east']/self.dscale
                self.station_north = md_data.coord_array['rel_north']/self.dscale
                self.station_names = md_data.coord_array['station']
            else:
                print 'Could not find data file {0}'.format(self.data_fn)
        
    def plot(self):
        """
        plot depth slices
        """
        #--> get information from files
        self.read_files()

        fdict = {'size':self.font_size+2, 'weight':'bold'}
        
        cblabeldict={-2:'$10^{-3}$',-1:'$10^{-1}$',0:'$10^{0}$',1:'$10^{1}$',
                     2:'$10^{2}$',3:'$10^{3}$',4:'$10^{4}$',5:'$10^{5}$',
                     6:'$10^{6}$',7:'$10^{7}$',8:'$10^{8}$'}
                     
        #create an list of depth slices to plot
        if self.depth_index == None:
            zrange = range(self.grid_z.shape[0])
        elif type(self.depth_index) is int:
            zrange = [self.depth_index]
        elif type(self.depth_index) is list or \
             type(self.depth_index) is np.ndarray:
            zrange = self.depth_index
        
        #set the limits of the plot
        if self.ew_limits == None:
            if self.station_east is not None:
                xlimits = (np.floor(self.station_east.min()), 
                           np.ceil(self.station_east.max()))
            else:
                xlimits = (self.grid_east[5], self.grid_east[-5])
        else:
            xlimits = self.ew_limits
            
        if self.ns_limits == None:
            if self.station_north is not None:
                ylimits = (np.floor(self.station_north.min()), 
                           np.ceil(self.station_north.max()))
            else:
                ylimits = (self.grid_north[5], self.grid_north[-5])
        else:
            ylimits = self.ns_limits
            
            
        #make a mesh grid of north and east
        self.mesh_east, self.mesh_north = np.meshgrid(self.grid_east, 
                                                      self.grid_north,
                                                      indexing='ij')
        
        plt.rcParams['font.size'] = self.font_size
        
        #--> plot depths into individual figures
        for ii in zrange: 
            depth = '{0:.3f} ({1})'.format(self.grid_z[ii], 
                                     self.map_scale)
            fig = plt.figure(depth, figsize=self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            ax1 = fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)
            plot_res = np.log10(self.res_model[:, :, ii].T)
            mesh_plot = ax1.pcolormesh(self.mesh_east,
                                       self.mesh_north, 
                                       plot_res,
                                       cmap=self.cmap,
                                       vmin=self.climits[0],
                                       vmax=self.climits[1])
                           
            #plot the stations
            if self.station_east is not None:
                for ee, nn in zip(self.station_east, self.station_north):
                    ax1.text(ee, nn, '*', 
                             verticalalignment='center',
                             horizontalalignment='center',
                             fontdict={'size':5, 'weight':'bold'})
    
            #set axis properties
            ax1.set_xlim(xlimits)
            ax1.set_ylim(ylimits)
            ax1.xaxis.set_minor_locator(MultipleLocator(self.xminorticks/self.dscale))
            ax1.yaxis.set_minor_locator(MultipleLocator(self.yminorticks/self.dscale))
            ax1.set_ylabel('Northing ('+self.map_scale+')',fontdict=fdict)
            ax1.set_xlabel('Easting ('+self.map_scale+')',fontdict=fdict)
            ax1.set_title('Depth = {0}'.format(depth), fontdict=fdict)
                       
            #plot the grid if desired
            if self.plot_grid == 'y':
                east_line_xlist = []
                east_line_ylist = []            
                for xx in self.grid_east:
                    east_line_xlist.extend([xx, xx])
                    east_line_xlist.append(None)
                    east_line_ylist.extend([self.grid_north.min(), 
                                            self.grid_north.max()])
                    east_line_ylist.append(None)
                ax1.plot(east_line_xlist,
                              east_line_ylist,
                              lw=.25,
                              color='k')
        
                north_line_xlist = []
                north_line_ylist = [] 
                for yy in self.grid_north:
                    north_line_xlist.extend([self.grid_east.min(),
                                             self.grid_east.max()])
                    north_line_xlist.append(None)
                    north_line_ylist.extend([yy, yy])
                    north_line_ylist.append(None)
                ax1.plot(north_line_xlist,
                              north_line_ylist,
                              lw=.25,
                              color='k')
            
                
            #plot the colorbar
            if self.cb_location is None:
                if self.cb_orientation == 'horizontal':
                    self.cb_location = (ax1.axes.figbox.bounds[3]-.225,
                                        ax1.axes.figbox.bounds[1]+.05,.3,.025) 
                                            
                elif self.cb_orientation == 'vertical':
                    self.cb_location = ((ax1.axes.figbox.bounds[2]-.15,
                                        ax1.axes.figbox.bounds[3]-.21,.025,.3))
            
            ax2 = fig.add_axes(self.cb_location)
            
            cb = mcb.ColorbarBase(ax2,
                                  cmap=self.cmap,
                                  norm=Normalize(vmin=self.climits[0],
                                                 vmax=self.climits[1]),
                                  orientation=self.cb_orientation)
                                
            if self.cb_orientation == 'horizontal':
                cb.ax.xaxis.set_label_position('top')
                cb.ax.xaxis.set_label_coords(.5,1.3)
                
                
            elif self.cb_orientation == 'vertical':
                cb.ax.yaxis.set_label_position('right')
                cb.ax.yaxis.set_label_coords(1.25,.5)
                cb.ax.yaxis.tick_left()
                cb.ax.tick_params(axis='y',direction='in')
                                
            cb.set_label('Resistivity ($\Omega \cdot$m)',
                         fontdict={'size':self.font_size+1})
            cb.set_ticks(np.arange(self.climits[0],self.climits[1]+1))
            cb.set_ticklabels([cblabeldict[cc] 
                                for cc in np.arange(self.climits[0],
                                                    self.climits[1]+1)])
            
            self.fig_list.append(fig)
            
            #--> save plots to a common folder
            if self.save_plots == 'y':
                
                fig.savefig(os.path.join(self.save_path,
                            "Depth_{}_{:.4f}.png".format(ii, self.grid_z[ii])),
                            dpi=self.fig_dpi, bbox_inches='tight')
                fig.clear()
                plt.close()
    
            else:
                pass
            
    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()
        
    def update_plot(self, fig):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  
        
        Use this if you change an of the .fig or axes properties
        
        :Example: ::
            
            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()
        
        """

        fig.canvas.draw()
                          
    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """
        
        return ("Plots depth slices of model from WS3DINV")


#==============================================================================
# plot slices 
#==============================================================================
class PlotSlices(object):
    """
    plot all slices and be able to scroll through the model
    
    :Example: ::
    
        >>> import mtpy.modeling.ws3dinv as ws
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> sfn = r"/home/MT/ws3dinv/Inv1/WSStationLocations.txt"
        >>> # plot just first layer to check the formating        
        >>> pds = ws.PlotSlices(model_fn=mfn, station_fn=sfn)
        
    ======================= ===================================================
    Buttons                  Description    
    ======================= ===================================================
    'e'                     moves n-s slice east by one model block
    'w'                     moves n-s slice west by one model block
    'n'                     moves e-w slice north by one model block
    'm'                     moves e-w slice south by one model block
    'd'                     moves depth slice down by one model block
    'u'                     moves depth slice up by one model block
    ======================= ===================================================

    
    ======================= ===================================================
    Attributes              Description    
    ======================= ===================================================
    ax_en                   matplotlib.axes instance for depth slice  map view 
    ax_ez                   matplotlib.axes instance for e-w slice
    ax_map                  matplotlib.axes instance for location map
    ax_nz                   matplotlib.axes instance for n-s slice
    climits                 (min , max) color limits on resistivity in log 
                            scale. *default* is (0, 4)
    cmap                    name of color map for resisitiviy.
                            *default* is 'jet_r'
    data_fn                 full path to data file name
    dscale                  scaling parameter depending on map_scale
    east_line_xlist         list of line nodes of east grid for faster plotting
    east_line_ylist         list of line nodes of east grid for faster plotting
    ew_limits               (min, max) limits of e-w in map_scale units
                            *default* is None and scales to station area
    fig                     matplotlib.figure instance for figure
    fig_aspect              aspect ratio of plots. *default* is 1
    fig_dpi                 resolution of figure in dots-per-inch
                            *default* is 300
    fig_num                 figure instance number
    fig_size                [width, height] of figure window. 
                            *default* is [6,6]
    font_dict               dictionary of font keywords, internally created
    font_size               size of ticklables in points, axes labes are 
                            font_size+2. *default* is 7
    grid_east               relative location of grid nodes in e-w direction
                            in map_scale units
    grid_north              relative location of grid nodes in n-s direction
                            in map_scale units
    grid_z                  relative location of grid nodes in z direction
                            in map_scale units
    index_east              index value of grid_east being plotted
    index_north             index value of grid_north being plotted
    index_vertical          index value of grid_z being plotted
    initial_fn              full path to initial file
    key_press               matplotlib.canvas.connect instance
    map_scale               [ 'm' | 'km' ] scale of map. *default* is km
    mesh_east               np.meshgrid(grid_east, grid_north)[0]
    mesh_en_east            np.meshgrid(grid_east, grid_north)[0]
    mesh_en_north           np.meshgrid(grid_east, grid_north)[1]
    mesh_ez_east            np.meshgrid(grid_east, grid_z)[0]
    mesh_ez_vertical        np.meshgrid(grid_east, grid_z)[1]
    mesh_north              np.meshgrid(grid_east, grid_north)[1]
    mesh_nz_north           np.meshgrid(grid_north, grid_z)[0]
    mesh_nz_vertical        np.meshgrid(grid_north, grid_z)[1]
    model_fn                full path to model file
    ms                      size of station markers in points. *default* is 2
    nodes_east              relative distance betwen nodes in e-w direction
                            in map_scale units
    nodes_north             relative distance betwen nodes in n-s direction
                            in map_scale units
    nodes_z                 relative distance betwen nodes in z direction
                            in map_scale units
    north_line_xlist        list of line nodes north grid for faster plotting  
    north_line_ylist        list of line nodes north grid for faster plotting
    ns_limits               (min, max) limits of plots in n-s direction
                            *default* is None, set veiwing area to station area 
    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
                            *default* is 'y'
    res_model               np.ndarray(n_north, n_east, n_vertical) of 
                            model resistivity values in linear scale           
    station_color           color of station marker. *default* is black
    station_dict_east       location of stations for each east grid row
    station_dict_north      location of stations for each north grid row
    station_east            location of stations in east direction
    station_fn              full path to station file 
    station_font_color      color of station label 
    station_font_pad        padding between station marker and label
    station_font_rotation   angle of station label
    station_font_size       font size of station label
    station_font_weight     weight of font for station label
    station_id              [min, max] index values for station labels
    station_marker          station marker
    station_names           name of stations
    station_north           location of stations in north direction
    subplot_bottom          distance between axes and bottom of figure window
    subplot_hspace          distance between subplots in vertical direction
    subplot_left            distance between axes and left of figure window  
    subplot_right           distance between axes and right of figure window
    subplot_top             distance between axes and top of figure window
    subplot_wspace          distance between subplots in horizontal direction
    title                   title of plot 
    z_limits                (min, max) limits in vertical direction,
    ======================= ===================================================
    
    """
    
    def __init__(self, model_fn, data_fn=None, **kwargs):
        self.model_fn = model_fn
        self.data_fn = data_fn
        
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.font_size = kwargs.pop('font_size', 7)
        
        self.subplot_wspace = .20
        self.subplot_hspace = .30
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .97
        self.subplot_bottom = .1
        
        self.index_vertical = kwargs.pop('index_vertical', 0)
        self.index_east = kwargs.pop('index_east', 0)
        self.index_north = kwargs.pop('index_north', 0)
        
        self.cmap = kwargs.pop('cmap', 'jet_r')
        self.climits = kwargs.pop('climits', (0, 4))
        
        self.map_scale = kwargs.pop('map_scale', 'km')
        #make map scale
        if self.map_scale=='km':
            self.dscale=1000.
        elif self.map_scale=='m':
            self.dscale=1. 
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)
        self.z_limits = kwargs.pop('z_limits', None)
        
        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z  = None
        
        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None
        
        self.mesh_east = None
        self.mesh_north = None
        
        self.station_east = None
        self.station_north = None
        self.station_names = None
        
        self.station_id = kwargs.pop('station_id', None)
        self.station_font_size = kwargs.pop('station_font_size', 8)
        self.station_font_pad = kwargs.pop('station_font_pad', 1.0)
        self.station_font_weight = kwargs.pop('station_font_weight', 'bold')
        self.station_font_rotation = kwargs.pop('station_font_rotation', 60)
        self.station_font_color = kwargs.pop('station_font_color', 'k')
        self.station_marker = kwargs.pop('station_marker', 
                                         r"$\blacktriangledown$")
        self.station_color = kwargs.pop('station_color', 'k')
        self.ms = kwargs.pop('ms', 10)
        
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()
        
        
    def read_files(self):
        """
        read in the files to get appropriate information
        """
        #--> read in model file
        if self.model_fn is not None:
            if os.path.isfile(self.model_fn) == True:
                md_model = Model()
                md_model.read_model_file(self.model_fn)
                self.res_model = md_model.res_model
                self.grid_east = md_model.grid_east/self.dscale
                self.grid_north = md_model.grid_north/self.dscale
                self.grid_z = md_model.grid_z/self.dscale
                self.nodes_east = md_model.nodes_east/self.dscale
                self.nodes_north = md_model.nodes_north/self.dscale
                self.nodes_z = md_model.nodes_z/self.dscale
            else:
                raise mtex.MTpyError_file_handling(
                        '{0} does not exist, check path'.format(self.model_fn))
        
        #--> read in data file to get station locations
        if self.data_fn is not None:
            if os.path.isfile(self.data_fn) == True:
                md_data = Data()
                md_data.read_data_file(self.data_fn)
                self.station_east = md_data.coord_array['rel_east']/self.dscale
                self.station_north = md_data.coord_array['rel_north']/self.dscale
                self.station_names = md_data.coord_array['station']
            else:
                print 'Could not find data file {0}'.format(self.data_fn)
        
    def plot(self):
        """
        plot:
            east vs. vertical,
            north vs. vertical,
            east vs. north
            
        
        """
        
        self.read_files()
        
        self.get_station_grid_locations()
        
        print "=============== ==============================================="
        print "    Buttons                  Description                       "
        print "=============== ==============================================="
        print "     'e'          moves n-s slice east by one model block"
        print "     'w'          moves n-s slice west by one model block"
        print "     'n'          moves e-w slice north by one model block"
        print "     'm'          moves e-w slice south by one model block"
        print "     'd'          moves depth slice down by one model block"
        print "     'u'          moves depth slice up by one model block"
        print "=============== ==============================================="
        
        self.font_dict = {'size':self.font_size+2, 'weight':'bold'}
        
        #--> set default font size                           
        plt.rcParams['font.size'] = self.font_size
        
        #set the limits of the plot
        if self.ew_limits == None:
            if self.station_east is not None:
                self.ew_limits = (np.floor(self.station_east.min()), 
                                  np.ceil(self.station_east.max()))
            else:
                self.ew_limits = (self.grid_east[5], self.grid_east[-5])

        if self.ns_limits == None:
            if self.station_north is not None:
                self.ns_limits = (np.floor(self.station_north.min()), 
                                  np.ceil(self.station_north.max()))
            else:
                self.ns_limits = (self.grid_north[5], self.grid_north[-5])
        
        if self.z_limits == None:
                depth_limit = max([(abs(self.ew_limits[0])+abs(self.ew_limits[1])),
                                   (abs(self.ns_limits[0])+abs(self.ns_limits[1]))])
                self.z_limits = (-5000/self.dscale, depth_limit)
            
        
        self.fig = plt.figure(self.fig_num, figsize=self.fig_size,
                              dpi=self.fig_dpi)
        plt.clf()
        gs = gridspec.GridSpec(2, 2,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom, 
                               right=self.subplot_right, 
                               hspace=self.subplot_hspace)        
        
        #make subplots
        self.ax_ez = self.fig.add_subplot(gs[0, 0], aspect=self.fig_aspect)
        self.ax_nz = self.fig.add_subplot(gs[1, 1], aspect=self.fig_aspect)
        self.ax_en = self.fig.add_subplot(gs[1, 0], aspect=self.fig_aspect)
        self.ax_map = self.fig.add_subplot(gs[0, 1])
        
        #make grid meshes being sure the indexing is correct
        self.mesh_ez_east, self.mesh_ez_vertical = np.meshgrid(self.grid_east,
                                                               self.grid_z,
                                                               indexing='ij') 
        self.mesh_nz_north, self.mesh_nz_vertical = np.meshgrid(self.grid_north,
                                                                self.grid_z,
                                                                indexing='ij') 
        self.mesh_en_east, self.mesh_en_north = np.meshgrid(self.grid_east, 
                                                            self.grid_north,
                                                            indexing='ij')
                                                            
        #--> plot east vs vertical
        self._update_ax_ez()
        
        #--> plot north vs vertical
        self._update_ax_nz()
                              
        #--> plot east vs north
        self._update_ax_en()
                                 
        #--> plot the grid as a map view 
        self._update_map()
        
        #plot color bar
        cbx = mcb.make_axes(self.ax_map, fraction=.15, shrink=.75, pad = .15)
        cb = mcb.ColorbarBase(cbx[0],
                              cmap=self.cmap,
                              norm=Normalize(vmin=self.climits[0],
                                             vmax=self.climits[1]))

   
        cb.ax.yaxis.set_label_position('right')
        cb.ax.yaxis.set_label_coords(1.25,.5)
        cb.ax.yaxis.tick_left()
        cb.ax.tick_params(axis='y',direction='in')
                            
        cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size':self.font_size+1})
                     
        cb.set_ticks(np.arange(np.ceil(self.climits[0]),
                               np.floor(self.climits[1]+1)))
        cblabeldict={-2:'$10^{-3}$',-1:'$10^{-1}$',0:'$10^{0}$',1:'$10^{1}$',
                     2:'$10^{2}$',3:'$10^{3}$',4:'$10^{4}$',5:'$10^{5}$',
                     6:'$10^{6}$',7:'$10^{7}$',8:'$10^{8}$'}
        cb.set_ticklabels([cblabeldict[cc] 
                            for cc in np.arange(np.ceil(self.climits[0]),
                                                np.floor(self.climits[1]+1))])
                   
        plt.show()
        
        self.key_press = self.fig.canvas.mpl_connect('key_press_event',
                                                     self.on_key_press)


    def on_key_press(self, event):
        """
        on a key press change the slices
        
        """                                                            

        key_press = event.key
        
        if key_press == 'n':
            if self.index_north == self.grid_north.shape[0]:
                print 'Already at northern most grid cell'
            else:
                self.index_north += 1
                if self.index_north > self.grid_north.shape[0]:
                    self.index_north = self.grid_north.shape[0]
            self._update_ax_ez()
            self._update_map()
       
        if key_press == 'm':
            if self.index_north == 0:
                print 'Already at southern most grid cell'
            else:
                self.index_north -= 1 
                if self.index_north < 0:
                    self.index_north = 0
            self._update_ax_ez()
            self._update_map()
                    
        if key_press == 'e':
            if self.index_east == self.grid_east.shape[0]:
                print 'Already at eastern most grid cell'
            else:
                self.index_east += 1
                if self.index_east > self.grid_east.shape[0]:
                    self.index_east = self.grid_east.shape[0]
            self._update_ax_nz()
            self._update_map()
       
        if key_press == 'w':
            if self.index_east == 0:
                print 'Already at western most grid cell'
            else:
                self.index_east -= 1 
                if self.index_east < 0:
                    self.index_east = 0
            self._update_ax_nz()
            self._update_map()
                    
        if key_press == 'd':
            if self.index_vertical == self.grid_z.shape[0]:
                print 'Already at deepest grid cell'
            else:
                self.index_vertical += 1
                if self.index_vertical > self.grid_z.shape[0]:
                    self.index_vertical = self.grid_z.shape[0]
            self._update_ax_en()
            print 'Depth = {0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                                 self.map_scale)
       
        if key_press == 'u':
            if self.index_vertical == 0:
                print 'Already at surface grid cell'
            else:
                self.index_vertical -= 1 
                if self.index_vertical < 0:
                    self.index_vertical = 0
            self._update_ax_en()
            print 'Depth = {0:.5gf} ({1})'.format(self.grid_z[self.index_vertical],
                                                 self.map_scale)
                    
    def _update_ax_ez(self):
        """
        update east vs vertical plot
        """
        self.ax_ez.cla()
        plot_ez = np.log10(self.res_model[self.index_north, :, :]) 
        self.ax_ez.pcolormesh(self.mesh_ez_east,
                              self.mesh_ez_vertical, 
                              plot_ez,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        #plot stations
        for sx in self.station_dict_north[self.grid_north[self.index_north]]:
            self.ax_ez.text(sx,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size':self.ms,
                                      'color':self.station_color})
                                      
        self.ax_ez.set_xlim(self.ew_limits)
        self.ax_ez.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_ez.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_ez.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()
        
    def _update_ax_nz(self):
        """
        update east vs vertical plot
        """
        self.ax_nz.cla()
        plot_nz = np.log10(self.res_model[:, self.index_east, :]) 
        self.ax_nz.pcolormesh(self.mesh_nz_north,
                              self.mesh_nz_vertical, 
                              plot_nz,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        #plot stations
        for sy in self.station_dict_east[self.grid_east[self.index_east]]:
            self.ax_nz.text(sy,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size':self.ms,
                                      'color':self.station_color})
        self.ax_nz.set_xlim(self.ns_limits)
        self.ax_nz.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_nz.set_xlabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_nz.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()
        
    def _update_ax_en(self):
        """
        update east vs vertical plot
        """
        
        self.ax_en.cla()
        plot_en = np.log10(self.res_model[:, :, self.index_vertical].T) 
        self.ax_en.pcolormesh(self.mesh_en_east,
                              self.mesh_en_north, 
                              plot_en,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        self.ax_en.set_xlim(self.ew_limits)
        self.ax_en.set_ylim(self.ns_limits)
        self.ax_en.set_ylabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_en.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        #--> plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax_en.text(ee, nn, '*', 
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 fontdict={'size':5, 'weight':'bold'})

        self.fig.canvas.draw()
        self._update_map()
        
    def _update_map(self):
        self.ax_map.cla()
        self.east_line_xlist = []
        self.east_line_ylist = []            
        for xx in self.grid_east:
            self.east_line_xlist.extend([xx, xx])
            self.east_line_xlist.append(None)
            self.east_line_ylist.extend([self.grid_north.min(), 
                                         self.grid_north.max()])
            self.east_line_ylist.append(None)
        self.ax_map.plot(self.east_line_xlist,
                         self.east_line_ylist,
                         lw=.25,
                         color='k')

        self.north_line_xlist = []
        self.north_line_ylist = [] 
        for yy in self.grid_north:
            self.north_line_xlist.extend([self.grid_east.min(),
                                          self.grid_east.max()])
            self.north_line_xlist.append(None)
            self.north_line_ylist.extend([yy, yy])
            self.north_line_ylist.append(None)
        self.ax_map.plot(self.north_line_xlist,
                         self.north_line_ylist,
                         lw=.25,
                         color='k')
        #--> e-w indication line 
        self.ax_map.plot([self.grid_east.min(), 
                          self.grid_east.max()],
                         [self.grid_north[self.index_north+1], 
                          self.grid_north[self.index_north+1]],
                         lw=1,
                         color='g')
                         
        #--> e-w indication line 
        self.ax_map.plot([self.grid_east[self.index_east+1], 
                          self.grid_east[self.index_east+1]],
                         [self.grid_north.min(), 
                          self.grid_north.max()],
                         lw=1,
                         color='b')
         #--> plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax_map.text(ee, nn, '*', 
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 fontdict={'size':5, 'weight':'bold'})
                                 
        self.ax_map.set_xlim(self.ew_limits)
        self.ax_map.set_ylim(self.ns_limits)
        self.ax_map.set_ylabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_map.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        
        #plot stations                      
        self.ax_map.text(self.ew_limits[0]*.95, self.ns_limits[1]*.95,
                 '{0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                        self.map_scale),
                 horizontalalignment='left',
                 verticalalignment='top',
                 bbox={'facecolor': 'white'},
                 fontdict=self.font_dict)
        
        
        self.fig.canvas.draw()
        
    def get_station_grid_locations(self):
        """
        get the grid line on which a station resides for plotting
        
        """
        self.station_dict_east = dict([(gx, []) for gx in self.grid_east])
        self.station_dict_north = dict([(gy, []) for gy in self.grid_north])
        if self.station_east is not None:
            for ss, sx in enumerate(self.station_east):
                gx = np.where(self.grid_east <= sx)[0][-1]
                self.station_dict_east[self.grid_east[gx]].append(self.station_north[ss])
            
            for ss, sy in enumerate(self.station_north):
                gy = np.where(self.grid_north <= sy)[0][-1]
                self.station_dict_north[self.grid_north[gy]].append(self.station_east[ss])
        else:
            return 
                  
    def redraw_plot(self):
        """
        redraw plot if parameters were changed
        
        use this function if you updated some attributes and want to re-plot.
        
        :Example: ::
            
            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        
        plt.close(self.fig)
        self.plot()
            
    def save_figure(self, save_fn=None, fig_dpi=None, file_format='pdf', 
                    orientation='landscape', close_fig='y'):
        """
        save_figure will save the figure to save_fn.
        
        Arguments:
        -----------
        
            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format
                            
                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path
                            
            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 
                              
            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait
                              
            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.
                          
            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open
                          
        :Example: ::
            
            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')
            
        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')
            
        else:
            save_fn = os.path.join(save_fn, '_E{0}_N{1}_Z{2}.{3}'.format(
                                    self.index_east, self.index_north,
                                    self.index_vertical, file_format))
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')
        
        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)
        
        else:
            pass
        
        self.fig_fn = save_fn
        print 'Saved figure to: '+self.fig_fn

#==============================================================================
# Exceptions
#==============================================================================
class ModEMError(Exception):
    pass