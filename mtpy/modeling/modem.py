#!/usr/bin/env python
"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017

"""
#==============================================================================
# Imports
#==============================================================================
# general packages
import os
import numpy as np
import scipy.interpolate as spi
import scipy.stats as stats

# mtpy modules
import mtpy.core.z as mtz
import mtpy.core.mt as mt
import mtpy.utils.gis_tools as gis_tools
import mtpy.modeling.ws3dinv as ws
import mtpy.imaging.mtplottools as mtplottools
import mtpy.utils.exceptions as mtex
import mtpy.analysis.pt as mtpt
import mtpy.imaging.mtcolors as mtcl
import mtpy.utils.configfile as mtcfg

# Plotting tools
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import matplotlib.widgets as widgets
import matplotlib.colors as colors
import matplotlib.cm as cm

# vtk tools
try:
    from evtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')

#==============================================================================
class Stations(object):
    """
    station locations class
    
    ..note:: If the survey steps across multiple UTM zones, then a 
             distance will be added to the stations to place them in 
             the correct location.  This distance is 
             _utm_grid_size_north and _utm_grid_size_east.  You should 
             these parameters to place the locations in the proper spot
             as grid distances and overlaps change over the globe.
             **This is not implemented yet**
    """
    
    def __init__(self, **kwargs):
        
        self.dtype = [('station', '|S10'),
                       ('lat', np.float),
                       ('lon', np.float),
                       ('elev', np.float),
                       ('rel_east', np.float), 
                       ('rel_north', np.float), 
                       ('east', np.float),
                       ('north', np.float),
                       ('zone', 'S4')]
        self.station_locations = np.zeros(0, dtype=self.dtype)
    
    ## --> define properties that can only be returned and not set    
    @property
    def lat(self):
        return self.station_locations['lat']
    
    @property
    def lon(self):
        return self.station_locations['lon']
    
    @property
    def east(self):
        return self.station_locations['east']
    
    @property
    def north(self):
        return self.station_locations['north']

    @property
    def elev(self):
        return self.station_locations['elev']
        
    @property
    def rel_east(self):
        return self.station_locations['rel_east']
        
    @property
    def rel_north(self):
        return self.station_locations['rel_north']
        
    @property
    def utm_zone(self):
        return self.station_locations['zone']
        
    @property
    def station(self):
        return self.station_locations['station']
    
    def _get_mt_objs_from_list(self, input_list):
        """
        get mt_objects from a list of files or mt_objects
        """
        
        if type(input_list) not in [list, np.ndarray]:
            raise ValueError('Input list needs to be type list, not {0}'.format(type(input_list)))
        
        if type(input_list[0]) is mt.MT:
            return input_list
        
        if type(input_list[0]) is str:
            if input_list[0].endswith('.edi'):
                return [mt.MT(fn) for fn in input_list]
            
            else:
                raise ModEMError('file {0} not supported yet'.format(input_list[0][-4:]))
        
    def get_station_locations(self, input_list):
        """
        get station locations from a list of edi files
        
        Arguments
        -------------
            **input_list** : list
                             list of edi file names, or mt_objects
                             
        
        Returns
        ------------
            * fills station_locations array
            
        """
                       
        mt_obj_list = self._get_mt_objs_from_list(input_list)
        
        #if station locations are not input read from the edi files
        if mt_obj_list is None:
            raise AttributeError('mt_obj_list is None, need to input a list of '
                                 'mt objects to read in.')
                                 
        n_stations = len(mt_obj_list)
        
        if n_stations == 0:
            raise ModEMError('No .edi files in edi_list, please check '
                             'file locations.')
        
        #make a structured array to put station location information into
        self.station_locations = np.zeros(n_stations,
                                          dtype=self.dtype)
        #get station locations in meters
        for ii, mt_obj in enumerate(mt_obj_list):
            self.station_locations[ii]['lat'] = mt_obj.lat
            self.station_locations[ii]['lon'] = mt_obj.lon
            self.station_locations[ii]['station'] = mt_obj.station
            self.station_locations[ii]['east'] = mt_obj.east
            self.station_locations[ii]['north'] = mt_obj.north
            self.station_locations[ii]['elev'] = mt_obj.elev
            self.station_locations[ii]['zone'] = mt_obj.utm_zone
        
        # get relative station locations
        self.calculate_rel_locations()
            
    
    def calculate_rel_locations(self, shift_east=0, shift_north=0):
        """
        put station in a coordinate system relative to 
        (shift_east, shift_north)
        (+) shift right or up
        (-) shift left or down
        
        """

        #remove the average distance to get coordinates in a relative space
        self.station_locations['rel_east'] = self.east-self.east.mean()
        self.station_locations['rel_north'] = self.north-self.north.mean()
        
        #translate the stations so they are relative to 0,0
        east_center = (self.rel_east.max()-np.abs(self.rel_east.min()))/2
        north_center = (self.rel_north.max()-np.abs(self.rel_north.min()))/2

        #remove the average distance to get coordinates in a relative space
        self.station_locations['rel_east'] -= east_center+shift_east
        self.station_locations['rel_north'] -= north_center+shift_north
        
                                                   
    # make center point a get method, can't set it.
    @property
    def center_point(self):
        """
        calculate the center point from the given station locations
        
        Returns
        -------------
            **center_location** : np.ndarray
                                  structured array of length 1
                                  dtype includes (east, north, zone, lat, lon)
        """
        center_location = np.recarray(1, dtype=self.dtype)
        center_point = np.array([self.east.mean(), self.north.mean()])
        
        #translate the stations so they are relative to 0,0
        east_center = (self.rel_east.max()-np.abs(self.rel_east.min()))/2
        north_center = (self.rel_north.max()-np.abs(self.rel_north.min()))/2
        
        center_point[0] -= east_center
        center_point[1] -= north_center
        
        # calculate center point in lat, lon, easting, northing        
        center_location['east'] = center_point[0]
        center_location['north'] = center_point[1]
        center_location['zone'] = self.utm_zone[0]

        center_ll = gis_tools.project_point_utm2ll(float(center_point[0]),
                                                   float(center_point[1]),
                                                   self.utm_zone[0])
                                                   
        center_location['lat'] = center_ll[0]
        center_location['lon'] = center_ll[1]
        
        return center_location
        
    def rotate_stations(self, rotation_angle):
        """
        Rotate stations assuming N is 0
        
        Arguments
        -------------
            **rotation_angle** : float
                                 angle in degrees assuming N is 0
        
        Returns
        -------------
            * refils rel_east and rel_north in station_locations.  Does this
              because you will still need the original locations for plotting
              later.
              
        """

        cos_ang = np.cos(np.deg2rad(rotation_angle))
        sin_ang = np.sin(np.deg2rad(rotation_angle))
        rot_matrix = np.matrix(np.array([[cos_ang, sin_ang], 
                                         [-sin_ang, cos_ang]]))
                                         
        coords = np.array([self.station_locations['rel_east'],
                           self.station_locations['rel_north']])
        
        #rotate the relative station locations
        new_coords = np.array(np.dot(rot_matrix, coords))
        
        self.station_locations['rel_east'][:] = new_coords[0, :]
        self.station_locations['rel_north'][:] = new_coords[1, :]
        
        print 'Rotated stations by {0:.1f} deg clockwise from N'.format(
                                                self.mesh_rotation_angle)
        
    def check_utm_crossing(self):
        """
        If the stations cross utm zones, then estimate distance by computing
        distance on a sphere.
        """
#        
#        latMid = (Lat1+Lat2 )/2.0;  // or just use Lat1 for slightly less accurate estimate
#
#
#        m_per_deg_lat = 111132.954 - 559.822 * cos( 2.0 * latMid ) + 1.175 * cos( 4.0 * latMid);
#        m_per_deg_lon = (3.14159265359/180 ) * 6367449 * cos ( latMid );
#        
#        deltaLat = fabs(Lat1 - Lat2);
#        deltaLon = fabs(Lon1 - Lon2);
#        
#        dist_m = sqrt (  pow( deltaLat * m_per_deg_lat,2) + pow( deltaLon * m_per_deg_lon , 2) );
#        
        pass



class Data(object):
    """
    Data will read and write .dat files for ModEM and convert a WS data file 
    to ModEM format.
    
    ..note: :: the data is interpolated onto the given periods such that all
               stations invert for the same periods.  The interpolation is
               a linear interpolation of each of the real and imaginary parts 
               of the impedance tensor and induction tensor.    
               See mtpy.core.mt.MT.interpolate for more details
    
    Arguments
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
    station_locations      Stations object
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
                               * zone --> UTM zone
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
    error_type_tipper      [ 'abs' | 'floor' ] 
                           *default* is 'abs'
    error_type_z           [ 'egbert' | 'mean_od' | 'eigen' ]
                           *default* is 'egbert_floor'
                                * add '_floor' to any of the above to set the
                                  error as an error floor, otherwise all 
                                  components are give weighted the same

                                * 'egbert' sets error to  
                                           error_value_z * sqrt(abs(zxy*zyx))
                                * 'mean_od' sets error to 
                                           error_value_z * mean([Zxy, Zyx])
                                * 'eigen' sets error to
                                          error_value_z * eigenvalues(Z[ii])

                                           
    error_value_z            percentage to multiply Z by to set error
                           *default* is 5 for 5% of Z as error
    error_value_tipper     absolute error between 0 and 1.
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
    rotate_angle           Angle to rotate data to assuming 0 is N and E is 90            
    save_path              path to save data file to
    units                  [ [V/m]/[T] | [mV/km]/[nT] | Ohm ] units of Z
                           *default* is [mV/km]/[nT]
    wave_sign_impedance    [ + | - ] sign of time dependent wave.  
                           *default* is '+' as positive downwards. 
    wave_sign_tipper       [ + | - ] sign of time dependent wave.  
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
                               filling in station_locations
    read_data_file             read in a ModEM data file and fill attributes
                               data_array, station_locations, period_list, mt_dict
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
        
    :Example 5 --> rotate data: ::
        
        >>> md.rotation_angle = 60
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")
        >>> # or
        >>> md.write_data_file(save_path=r"/home/modem/Inv1", \
                               rotation_angle=60)
    
                  
    """
    
    def __init__(self, edi_list=None, **kwargs):
        self.edi_list = edi_list

        self.error_type_z = 'egbert_floor'
        self.error_value_z = 5.0
        
        self.error_value_tipper = .05
        self.error_type_tipper = 'abs'
        
        self.wave_sign_impedance = '+'
        self.wave_sign_tipper = '+'
        self.units = '[mV/km]/[nT]'
        self.inv_mode = '1'
        
        self.period_list = None
        self.period_step = 1
        self.period_min = None
        self.period_max = None
        self.period_buffer = None
        self.max_num_periods = None
        self.data_period_list = None
        
        self.data_fn = 'ModEM_Data.dat'
        self.save_path = os.getcwd()
        self.formatting = '1'
        
        self._rotation_angle = 0.0
        self._set_rotation_angle(self._rotation_angle)
        
        self.center_point = None
        
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
                       ('zone', '|S4'),
                       ('z', (np.complex, self._z_shape)),
                       ('z_err', (np.float, self._z_shape)),
                       ('z_inv_err', (np.float, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.float, self._t_shape)),
                       ('tip_inv_err', (np.float, self._t_shape))]
        
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
        
                              
        self.header_string = ' '.join(['# Period(s)',
                                       'Code',
                                       'GG_Lat',
                                       'GG_Lon',
                                       'X(m)',
                                       'Y(m)',
                                       'Z(m)',
                                       'Component',
                                       'Real',
                                       'Imag',
                                       'Error\n'])
       
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])        
        
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
                       ('zone', '|S4'),
                       ('z', (np.complex, self._z_shape)),
                       ('z_err', (np.float, self._z_shape)),
                       ('z_inv_err', (np.float, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.float, self._t_shape)),
                       ('tip_inv_err', (np.float, self._t_shape))]
                       
    def get_header_string(self, error_type, error_value, rotation_angle):
        """
        reset the header sring for file
        """
        
        h_str = ','.join(['# Created using MTpy calculated {0} error of {1:.0f}%',
                          ' data rotated {2:.1f}_deg clockwise from N\n'])
                 
        return h_str.format(error_type, error_value, rotation_angle)
            

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
        
    def get_relative_station_locations(self):
        """
        get station locations from edi files
        """
        
        stations_obj = Stations()
        mt_list = [self.mt_dict[s_key] for s_key in sorted(self.mt_dict.keys())]
        
        stations_obj.get_station_locations(mt_list)
        
        # rotate locations if needed
        if self._rotation_angle != 0:
            stations_obj.rotate_stations(self._rotation_angle)

        # fill data array
        self.data_array[:]['station'] = stations_obj.station
        self.data_array[:]['lat'] = stations_obj.lat
        self.data_array[:]['lon'] = stations_obj.lon
        self.data_array[:]['east'] = stations_obj.east
        self.data_array[:]['north'] = stations_obj.north
        self.data_array[:]['elev'] = stations_obj.elev
        self.data_array[:]['rel_east'] = stations_obj.rel_east
        self.data_array[:]['rel_north'] = stations_obj.rel_north
        
        # get center point
        self.center_point = stations_obj.center_point
        
        
    def get_period_list(self):
        """
        make a period list to invert for
        
        """
        if self.mt_dict is None:
            self.get_mt_dict()
            
        if self.period_list is not None:
            print '-'*50
            print 'Inverting for periods:'
            for per in self.period_list:
                print '     {0:<12.6f}'.format(per)
            print '-'*50
            return

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
        

    def _set_rotation_angle(self, rotation_angle):
        """
        on set rotation angle rotate mt_dict and data_array, 
        """
        if self._rotation_angle == rotation_angle:
            return
  
        print 'Changing rotation angle from {0:.1f} to {1:.1f}'.format(
                                    self._rotation_angle, rotation_angle)
        
        self._rotation_angle = -self._rotation_angle+rotation_angle
        
        if self.rotation_angle == 0:
            return
            
        print 'Changing rotation angle from {0:.1f} to {1:.1f}'.format(
                                    self._rotation_angle, rotation_angle)
        self._rotation_angle = rotation_angle
        
            
        if self.data_array is None:
            return
        if self.mt_dict is None:
            return
            
        for mt_key in sorted(self.mt_dict.keys()):
            mt_obj = self.mt_dict[mt_key]
            mt_obj.Z.rotate(self._rotation_angle)
            mt_obj.Tipper.rotate(self._rotation_angle)
            
        print 'Data rotated to align with {0:.1f} deg clockwise from N'.format(
                                                        self._rotation_angle)
                 
        print '*'*70                                       
        print '   If you want to rotate station locations as well use the'
        print '   command Data.get_relative_station_locations() '
        print '   if stations have not already been rotated in Model'
        print '*'*70
        
        self.fill_data_array()
                
    def _get_rotation_angle(self):
        return self._rotation_angle
        
    rotation_angle = property(fget=_get_rotation_angle, 
                              fset=_set_rotation_angle,
                              doc="""Rotate data assuming N=0, E=90""")
                              
    def fill_data_array(self):
        """
        fill the data array from mt_dict
        
        """
        
        if self.period_list is None:
            self.get_period_list()
            
        ns = len(self.mt_dict.keys())
        nf = len(self.period_list)

        d_array = False        
        if self.data_array is not None:
            d_arr_copy = self.data_array.copy()
            d_array = True
            
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)
        
        rel_distance = True                                    
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]
            if d_array is True:
                try:
                    d_index = np.where(d_arr_copy['station'] == s_key)[0][0]
                    self.data_array[ii]['station'] = s_key
                    self.data_array[ii]['lat'] = d_arr_copy[d_index]['lat']
                    self.data_array[ii]['lon'] = d_arr_copy[d_index]['lon']
                    self.data_array[ii]['east'] = d_arr_copy[d_index]['east']
                    self.data_array[ii]['north'] = d_arr_copy[d_index]['north']
                    self.data_array[ii]['elev'] = d_arr_copy[d_index]['elev']
                    self.data_array[ii]['rel_east'] = d_arr_copy[d_index]['rel_east']
                    self.data_array[ii]['rel_north'] = d_arr_copy[d_index]['rel_north']
                except IndexError:
                    print 'Could not find {0} in data_array'.format(s_key)
            else:
                self.data_array[ii]['station'] = mt_obj.station
                self.data_array[ii]['lat'] = mt_obj.lat
                self.data_array[ii]['lon'] = mt_obj.lon
                self.data_array[ii]['east'] = mt_obj.east
                self.data_array[ii]['north'] = mt_obj.north
                self.data_array[ii]['elev'] = mt_obj.elev
                try:
                    self.data_array[ii]['rel_east'] = mt_obj.grid_east
                    self.data_array[ii]['rel_north'] = mt_obj.grid_north
                    rel_distance = False
                except AttributeError:
                    pass
                            
            # interpolate each station onto the period list
            # check bounds of period list
            interp_periods = self.period_list[np.where(
                                (self.period_list >= 1./mt_obj.Z.freq.max()) & 
                                (self.period_list <= 1./mt_obj.Z.freq.min()))]
            
            interp_z, interp_t = mt_obj.interpolate(1./interp_periods)
            for kk, ff in enumerate(interp_periods):
                jj = np.where(self.period_list == ff)[0][0]
                self.data_array[ii]['z'][jj] = interp_z.z[kk, :, :]
                self.data_array[ii]['z_err'][jj] = interp_z.z_err[kk, :, :]

                if mt_obj.Tipper.tipper is not None:
                    self.data_array[ii]['tip'][jj] = interp_t.tipper[kk, :, :]
                    self.data_array[ii]['tip_err'][jj] = \
                                                interp_t.tipper_err[kk, :, :]
        
        if rel_distance is False:
            self.get_relative_station_locations()
            
    def _set_station_locations(self, station_obj=None):
        """
        take a station_locations array and populate data_array
        """

        if station_obj is not None:
            station_locations = station_obj.station_locations
            

        if self.data_array is None:
            self._set_dtype((len(self.period_list), 2, 2),
                            (len(self.period_list), 1, 2))
            self.data_array = np.zeros(station_locations.size, 
                                       dtype=self._dtype)
            for d_index, s_arr in enumerate(station_locations):
                self.data_array[d_index]['lat'] = s_arr['lat']
                self.data_array[d_index]['lon'] = s_arr['lon']
                self.data_array[d_index]['east'] = s_arr['east']
                self.data_array[d_index]['north'] = s_arr['north']
                self.data_array[d_index]['elev'] = s_arr['elev']
                self.data_array[d_index]['rel_east'] = s_arr['rel_east']
                self.data_array[d_index]['rel_north'] = s_arr['rel_north']
                
        else:
            for s_arr in station_locations:
                try:
                    d_index = np.where(self.data_array['station'] == 
                                        s_arr['station'])[0][0]
                except IndexError:
                    print 'Could not find {0} in data_array'.format(s_arr['station'])
                    d_index = None
                
                if d_index is not None:
                    self.data_array[d_index]['lat'] = s_arr['lat']
                    self.data_array[d_index]['lon'] = s_arr['lon']
                    self.data_array[d_index]['east'] = s_arr['east']
                    self.data_array[d_index]['north'] = s_arr['north']
                    self.data_array[d_index]['elev'] = s_arr['elev']
                    self.data_array[d_index]['rel_east'] = s_arr['rel_east']
                    self.data_array[d_index]['rel_north'] = s_arr['rel_north']
                
    def _get_station_locations(self):
        """
        extract station locations from data array
        """
        if self.data_array is None:
            return None
            
        station_locations = self.data_array[['station', 'lat', 'lon', 
                                             'north', 'east', 'elev',
                                             'rel_north', 'rel_east']]
        station_obj = Stations()
        station_obj.station_locations = station_locations
        
        return station_obj
        
    station_locations = property(_get_station_locations, 
                                  _set_station_locations,
                                  doc="""location of stations""") 
                                  
    def compute_inv_error(self):
        """
        compute the error from the given parameters
        """
        # copy values over to inversion error
        self.data_array['z_inv_err'] = self.data_array['z_err']
        self.data_array['tip_inv_err'] = self.data_array['tip_err']        
       
       #compute relative error for tipper
        if 'floor' in self.error_type_tipper:
            t_index = np.where(self.data_array['tip_err'] < self.error_value_tipper)
            self.data_array['tip_inv_err'][t_index] = self.error_value_tipper
        elif 'abs' in self.error_type_tipper:
            self.data_array['tip_inv_err'][:] = self.error_value_tipper
        
        # compute error for z
        err_value = self.error_value_z/100.
        for ss in range(self.data_array.shape[0]):
            for ff in range(max([self._t_shape[0], self._z_shape[0]])):
                d_xx = abs(self.data_array['z'][ss, ff, 0, 0])
                d_xy = abs(self.data_array['z'][ss, ff, 0, 1])
                d_yx = abs(self.data_array['z'][ss, ff, 1, 0])
                d_yy = abs(self.data_array['z'][ss, ff, 1, 1])
                d = np.array([d_xx, d_xy, d_yx, d_yy])
                nz = np.nonzero(d)
                
                if d.sum() == 0.0:
                    continue
                
                if 'egbert' in self.error_type_z:
                    if d_xy == 0.0:
                        d_xy = 1.0
                    if d_yx == 0.0:
                        d_yx = 1.0
                    err = err_value*np.sqrt(d_xy*d_yx)
                    if err == 1.0:
                        err = max([d_xx, d_xy, d_yx, d_yy])*10
        
                elif 'median' in self.error_type_z:
                    err = err_value*np.median(d[nz])
                    
                elif 'mean_od' in self.error_type_z:
                    d = np.array(d_xy, d_yx)
                    nz = np.nonzero(d)
                    err = err_value*np.mean(d[nz])
                    
                elif 'eigen' in self.error_type_z:
                    d = d.reshape((2, 2))
                    err = err_value*np.abs(np.linalg.eigvals(d)).mean()
                    if err == 0:
                        err = err_value*d.flatten()[nz].mean()
            
                else:
                    raise NameError('{0} not understood'.format(self.error_type_z))
                
                self.data_array['z_inv_err'][ss, ff, :, :] = err
                

        # if there is an error floor                
        if 'floor' in self.error_type_z:
            f_index = np.where(self.data_array['z_inv_err'] < self.data_array['z_err'])
            self.data_array['z_inv_err'][f_index] = self.data_array['z_err'][f_index]
            
                
    def write_data_file(self, save_path=None, fn_basename=None, 
                        rotation_angle=None, compute_error=True, fill=True,
                        elevation=False):
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
                              
            **rotation_angle** : float
                                angle to rotate the data by assuming N = 0, 
                                E = 90. *default* is 0.0
                              
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
        
        self.get_period_list()
        
        #rotate data if desired
        if rotation_angle is not None:
            self.rotation_angle = rotation_angle
        
        #be sure to fill in data array
        if fill is True:
            self.fill_data_array()
            # get relative station locations in grid coordinates
            self.get_relative_station_locations()
            
        if elevation is False:
            self.data_array['elev'][:] = 0.0

        dlines = []        
        for inv_mode in self.inv_mode_dict[self.inv_mode]:
            if 'impedance' in inv_mode.lower():
                dlines.append(self.get_header_string(self.error_type_z,
                                                     self.error_value_z,
                                                     self.rotation_angle))
            elif 'vertical' in inv_mode.lower():
                dlines.append(self.get_header_string(self.error_type_tipper,
                                                     self.error_value_tipper,
                                                     self.rotation_angle))
            dlines.append(self.header_string)
            dlines.append('> {0}\n'.format(inv_mode))
            
            if inv_mode.find('Impedance') > 0:
                dlines.append('> exp({0}i\omega t)\n'.format(self.wave_sign_impedance))
                dlines.append('> {0}\n'.format(self.units))
            elif inv_mode.find('Vertical') >= 0:
                dlines.append('> exp({0}i\omega t)\n'.format(self.wave_sign_tipper))
                dlines.append('> []\n')
            dlines.append('> 0\n') #oriention, need to add at some point
            dlines.append('> {0:>10.6f} {1:>10.6f}\n'.format(
                          self.center_point.lat[0], self.center_point.lon[0]))
            dlines.append('> {0} {1}\n'.format(self.data_array['z'].shape[1],
                                               self.data_array['z'].shape[0]))
                                               
            if compute_error == True:
                self.compute_inv_error()
                                               
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
                        if zz.real != 0.0 and zz.imag != 0.0 and \
                            zz.real != 1e32 and zz.imag != 1e32:
                            if self.formatting == '1':
                                per = '{0:<12.5e}'.format(self.period_list[ff])
                                sta = '{0:>7}'.format(self.data_array[ss]['station'])
                                lat = '{0:> 9.3f}'.format(self.data_array[ss]['lat'])
                                lon = '{0:> 9.3f}'.format(self.data_array[ss]['lon'])
                                eas = '{0:> 12.3f}'.format(self.data_array[ss]['rel_east'])
                                nor = '{0:> 12.3f}'.format(self.data_array[ss]['rel_north'])
                                ele = '{0:> 12.3f}'.format(self.data_array[ss]['elev'])
                                com = '{0:>4}'.format(comp.upper())
                                if self.units == 'ohm':
                                    rea = '{0:> 14.6e}'.format(zz.real/796.)
                                    ima = '{0:> 14.6e}'.format(zz.imag/796.)
                                else:
                                    rea = '{0:> 14.6e}'.format(zz.real)
                                    ima = '{0:> 14.6e}'.format(zz.imag)
                                

                            elif self.formatting == '2':
                                per = '{0:<14.6e}'.format(self.period_list[ff])
                                sta = '{0:<10}'.format(self.data_array[ss]['station'])
                                lat = '{0:> 14.6f}'.format(self.data_array[ss]['lat'])
                                lon = '{0:> 14.6f}'.format(self.data_array[ss]['lon'])
                                eas = '{0:> 12.3f}'.format(self.data_array[ss]['rel_east'])
                                nor = '{0:> 15.3f}'.format(self.data_array[ss]['rel_north'])
                                ele = '{0:> 10.3f}'.format(self.data_array[ss]['elev'])
                                com = '{0:>12}'.format(comp.upper())
                                if self.units == 'ohm':
                                    rea = '{0:> 17.6e}'.format(zz.real/796.)
                                    ima = '{0:> 17.6e}'.format(zz.imag/796.)
                                else:
                                    rea = '{0:> 17.6e}'.format(zz.real)
                                    ima = '{0:> 17.6e}'.format(zz.imag)
                            
                            # get error from inversion error
                            abs_err = self.data_array['{0}_inv_err'.format(c_key)][ss, ff, z_ii, z_jj]
                            
                            if np.isinf(abs_err) or np.isnan(abs_err):
                                abs_err = 10**(np.floor(np.log10(abs(max([float(rea), 
                                                                           float(ima)])))))
                            abs_err = '{0:> 14.6e}'.format(abs(abs_err))
                            #make sure that x==north, y==east, z==+down                            
                            dline = ''.join([per, sta, lat, lon, nor, eas, ele, 
                                             com, rea, ima, abs_err, '\n'])
                            dlines.append(dline)
        
        with open(self.data_fn, 'w') as dfid:
            dfid.writelines(dlines)
        
        print 'Wrote ModEM data file to {0}'.format(self.data_fn)
        return self.data_fn
        
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
                                   
        #--> fill data array
        for ii, d_arr in enumerate(wsd.data):
            self.data_array[ii]['station'] = d_arr['station']
            self.data_array[ii]['rel_east'] = d_arr['east']
            self.data_array[ii]['rel_north'] = d_arr['north']
            self.data_array[ii]['z'][:] = d_arr['z_data']
            self.data_array[ii]['z_err'][:] = d_arr['z_data_err'].real*\
                                                d_arr['z_err_map'].real
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
        
    def convert_modem_to_ws(self, data_fn=None, ws_data_fn=None,
                            error_map=[1, 1, 1, 1]):
        """
        convert a ModEM data file to WS format.
        
        Arguments
        -------------
            **data_fn** : string
                         full path to modem data file to convert
                         
            **ws_data_fn** : string
                             full path to write ws format data file
            
            **error_map** : [zxx, zxy, zyx, zyy] floats
                            error map that ws uses, weights for each component
                            *default* is [1, 1, 1, 1] for equal weighting
        Returns
        ------------
            **ws_data_fn** : string
                             full path of ws data file
                             
            **ws_station_fn** : string 
                                full path to ws station file
                                
        Example
        -----------
            :Convert ModEM data file to WS: ::
                >>> import mtpy.modeling.modem as modem
                >>> md = modem.Data()
                >>> md.convert_modem_to_ws(data_fn=r"/home/mt/modem/data.dat")
        """
        
        if self.data_fn is not None:
            self.read_data_file(data_fn)
            
        if ws_data_fn == None:
            save_path = os.path.dirname(self.data_fn)
            ws_data_fn = os.path.join(save_path, 'WS_Data.dat')
            
        else:
            save_path = os.path.dirname(ws_data_fn)
            
        station_info = ws.WSStation()
        station_info.east = self.data_array['rel_east']
        station_info.north = self.data_array['rel_north']
        station_info.names = self.data_array['station']
        station_info.elev = self.data_array['elev']
        station_info.save_path = save_path
        station_info.write_station_file()
        
        ws_data = ws.WSData()
        ws_data.period_list = self.period_list.copy()
        ws_data.z_err_map = error_map
        ws_data.z_err = 'data'
        z_shape = (self.period_list.size, 2, 2)
        data_dtype = [('station', '|S10'),
                      ('east', np.float),
                      ('north', np.float),
                      ('z_data', (np.complex, z_shape)),
                      ('z_data_err', (np.complex, z_shape)),
                      ('z_err_map', (np.complex, z_shape))]
        ws_data.data = np.zeros(self.data_array['station'].size, 
                                dtype=data_dtype)
        ws_data.data['station'][:] = self.data_array['station']
        ws_data.data['east'] = self.data_array['rel_east']
        ws_data.data['north'] = self.data_array['rel_north']
        ws_data.data['z_data'][:, :, :] = self.data_array['z']
        ws_data.data['z_data_err'][:, :, :] = self.data_array['z_err']*(1+1j)
        ws_data.data['z_err_map'][:, :, :] = np.array([[1, 1], [1, 1]])
        
        ws_data.write_data_file(save_path=save_path, data_fn=ws_data_fn)
        
        return ws_data.data_fn, station_info.station_fn
        
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
            self.save_path = os.path.dirname(self.data_fn)
            self.fn_basename = os.path.basename(self.data_fn)
            
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
        read_impedance = False
        read_tipper = False
        inv_list = []
        for dline in dlines:
            if dline.find('#') == 0:
                header_list.append(dline.strip())
            elif dline.find('>') == 0:
                metadata_list.append(dline[1:].strip())
                if dline.lower().find('ohm') > 0:
                    self.units = 'ohm'
                elif dline.lower().find('mv') > 0:
                    self.units =' [mV/km]/[nT]'
                elif dline.lower().find('vertical') > 0:
                    read_tipper = True
                    read_impedance = False
                    inv_list.append('Full_Vertical_Components')
                elif dline.lower().find('impedance') > 0:
                    read_impedance = True
                    read_tipper = False
                    inv_list.append('Full_Impedance')
                if dline.find('exp') > 0:
                    if read_impedance is True:
                        self.wave_sign_impedance = dline[dline.find('(')+1]
                    elif read_tipper is True:
                        self.wave_sign_tipper = dline[dline.find('(')+1]
                elif len(dline[1:].strip().split()) == 2:
                    if dline.find('.') > 0:
                        value_list = [float(value) for value in 
                                      dline[1:].strip().split()]
                            
                        self.center_point = np.recarray(1, dtype=[('station', '|S10'),
                                                                   ('lat', np.float),
                                                                   ('lon', np.float),
                                                                   ('elev', np.float),
                                                                   ('rel_east', np.float), 
                                                                   ('rel_north', np.float), 
                                                                   ('east', np.float),
                                                                   ('north', np.float),
                                                                   ('zone', 'S4')])
                        self.center_point.lat = value_list[0]
                        self.center_point.lon = value_list[1]
                        
                        ce, cn, cz = gis_tools.project_point_ll2utm(self.center_point.lat,
                                                                    self.center_point.lon)
                        
                        self.center_point.east = ce
                        self.center_point.north = cn
                        self.center_point.zone = cz

                    else:
                        pass
                        
            else:
                dline_list = dline.strip().split()
                if len(dline_list) == 11:
                    for ii, d_str in enumerate(dline_list):
                        if ii != 1:
                            try:
                                dline_list[ii] = float(d_str.strip())
                            except ValueError:
                                pass
                        # be sure the station name is a string
                        else:
                            dline_list[ii] = d_str.strip()
                    period_list.append(dline_list[0])
                    station_list.append(dline_list[1])
                    
                    data_list.append(dline_list)
        
        #try to find rotation angle
        h_list = header_list[0].split()
        for hh, h_str in enumerate(h_list):
            if h_str.find('_deg') > 0:
                try:
                    self._rotation_angle = float(h_str[0:h_str.find('_deg')])
                    print ('Set rotation angle to {0:.1f} '.format(
                             self._rotation_angle)+'deg clockwise from N')
                except ValueError:
                    pass

        # find inversion mode
        for inv_key in self.inv_mode_dict.keys():
            inv_mode_list = self.inv_mode_dict[inv_key]
            if len(inv_mode_list) != inv_list:
                continue
            else:
                tf_arr = np.zeros(len(inv_list), dtype=np.bool)
            
                for tf, data_inv in enumerate(inv_list):
                    if data_inv in self.inv_mode_dict[inv_key]:
                        tf_arr[tf] = True
                
                if np.alltrue(tf_arr) == True:
                    self.inv_mode = inv_key
                    break
                                
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
                                        z_err_array=z_dummy.copy().real, 
                                        freq=1./self.period_list)
            data_dict[station].Tipper = mtz.Tipper(tipper_array=t_dummy.copy(), 
                                                   tipper_err_array=t_dummy.copy().real, 
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
                z_err = dd[10]
                if self.wave_sign_impedance == '+':
                    z_value = dd[8]+1j*dd[9]
                elif self.wave_sign_impedance == '-':
                    z_value = dd[8]-1j*dd[9]
                    
                if self.units == 'ohm':
                    z_value *= 796.
                    z_err *= 796.
                    
                data_dict[dd[1]].Z.z[p_index, ii, jj] = z_value
                data_dict[dd[1]].Z.z_err[p_index, ii, jj] = z_err
            #fill in tipper with appropriate values
            elif dd[7].find('T') == 0:
                if self.wave_sign_tipper == '+':
                    data_dict[dd[1]].Tipper.tipper[p_index, ii, jj] = dd[8]+1j*dd[9]
                elif self.wave_sign_tipper == '-':
                    data_dict[dd[1]].Tipper.tipper[p_index, ii, jj] = dd[8]-1j*dd[9]
                data_dict[dd[1]].Tipper.tipper_err[p_index, ii, jj] = dd[10]
       
        #make mt_dict an attribute for easier manipulation later
        self.mt_dict = data_dict

        ns = len(self.mt_dict.keys())
        nf = len(self.period_list)
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)
        
        #Be sure to caclulate invariants and phase tensor for each station
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]
            
            #self.mt_dict[s_key].zinv.compute_invariants()
            self.mt_dict[s_key].pt.set_z_object(mt_obj.Z)
            self.mt_dict[s_key].Tipper._compute_amp_phase()
            self.mt_dict[s_key].Tipper._compute_mag_direction()
                                   
        
            self.data_array[ii]['station'] = mt_obj.station
            self.data_array[ii]['lat'] = mt_obj.lat
            self.data_array[ii]['lon'] = mt_obj.lon
            self.data_array[ii]['east'] = mt_obj.east
            self.data_array[ii]['north'] = mt_obj.north
            self.data_array[ii]['elev'] = mt_obj.grid_elev
            self.data_array[ii]['rel_east'] = mt_obj.grid_east
            self.data_array[ii]['rel_north'] = mt_obj.grid_north
            
            self.data_array[ii]['z'][:] = mt_obj.Z.z
            self.data_array[ii]['z_err'][:] = mt_obj.Z.z_err
            self.data_array[ii]['z_inv_err'][:] = mt_obj.Z.z_err
            
            self.data_array[ii]['tip'][:] = mt_obj.Tipper.tipper
            self.data_array[ii]['tip_err'][:] = mt_obj.Tipper.tipper_err
            self.data_array[ii]['tip_inv_err'][:] = mt_obj.Tipper.tipper_err
            
    def write_vtk_station_file(self, vtk_save_path=None, 
                               vtk_fn_basename='ModEM_stations'):
        """
        write a vtk file for station locations.  For now this in relative 
        coordinates.  
        
        Arguments:
        -------------
            **vtk_save_path** : string
                                directory to save vtk file to.  
                                *default* is Model.save_path
            **vtk_fn_basename** : string
                                  filename basename of vtk file
                                  *default* is ModEM_stations, evtk will add
                                  on the extension .vtu
        """
        
        if vtk_save_path is None:
            vtk_fn = os.path.join(self.save_path, vtk_fn_basename)
        else:
            vtk_fn = os.path.join(vtk_save_path, vtk_fn_basename)
            
        pointsToVTK(vtk_fn, 
                 self.station_locations.rel_north/1000., 
                 self.station_locations.rel_east/1000.,
                 self.station_locations.elev/1000.,
                 data={'elevation':self.station_locations.elev})
                 
        print '--> Wrote station file to {0}'.format(vtk_fn)
        print '-'*50
            
    def get_parameters(self):
        """
        get important parameters for documentation
        """
        
        parameter_list = ['error_type_z',
                          'error_value_z',
                          'error_type_tipper',
                          'error_value_tipper',
                          'wave_sign_impedance',
                          'wave_sign_tipper',
                          'rotation_angle',
                          'save_path']
                          
        parameter_dict = {}
        for parameter in parameter_list:
            key = 'data.{0}'.format(parameter)
            parameter_dict[key] = getattr(self, parameter)
            
        parameter_dict['data.period_min'] = self.period_list.min()
        parameter_dict['data.period_max'] = self.period_list.max()
        parameter_dict['data.period_num'] = self.period_list.size
        
        parameter_dict['data.inv_mode'] = self.inv_mode_dict[self.inv_mode]
        parameter_dict['data.num_stations'] = self.station_locations.station.size
        parameter_dict['data.center_point_ll'] = (self.center_point.lat[0],
                                                  self.center_point.lon[0])
        
        parameter_dict['data.center_point_utm'] = (self.center_point.north[0],
                                                   self.center_point.east[0],
                                                   self.center_point.zone[0])
        return parameter_dict
        
        
    def center_stations(self, model_fn, data_fn=None):
        """
        Center station locations to the middle of cells, might be useful for 
        topography.
    
    
        Arguments
        -----------
            **data_fn** : string
                          full path to data file
                          
            **model_fn** : string
                          full path to model file
                          
            **new_data_fn** : string
                             full path to new data file
                             *default* is None, which save as 
                             data_fn_center.dat
                             
        Returns
        -----------
            **new_data_fn** : string
                              full path to new data file 
        """
        
        if data_fn is not None:
            self.read_data_file(data_fn)
        
        m_obj = Model()
        m_obj.read_model_file(model_fn)
        
        for s_arr in self.station_locations.station_locations:
            e_index = np.where(m_obj.grid_east >= s_arr['rel_east'])[0][0]-1
            n_index = np.where(m_obj.grid_north >= s_arr['rel_north'])[0][0]-1
            
            mid_east = m_obj.grid_east[e_index:e_index+2].mean()
            mid_north = m_obj.grid_north[n_index:n_index+2].mean()
            
            s_index = np.where(self.data_array['station']==s_arr['station'])[0][0]
            
            self.data_array[s_index]['rel_east'] = mid_east
            self.data_array[s_index]['rel_north'] = mid_north
    
    def change_data_elevation(self, model_fn, data_fn=None, 
                              res_air=1e12):
        """
        At each station in the data file rewrite the elevation, so the station is
        on the surface, not floating in air.
        
        Arguments:
        ------------------
            *data_fn* : string
                        full path to a ModEM data file
                        
            *model_fn* : string
                        full path to ModEM model file that has elevation 
                        incoorporated.
                                            
            *new_data_fn* : string
                            full path to new data file name.  If None, then 
                            new file name will add _elev.dat to input filename
                            
            *res_air* : float
                        resistivity of air.  Default is 1E12 Ohm-m
        Returns:
        -------------
            *new_data_fn* : string
                            full path to new data file.
        """
        if data_fn is not None:
            self.read_data_file(data_fn)
        
        m_obj = Model()
        m_obj.read_model_file(model_fn)
        
        s_locations = self.station_locations.station_locations.copy()
        
        # need to subtract one because we are finding the cell next to it
        for s_arr in s_locations:
            e_index = np.where(m_obj.grid_east >= s_arr['rel_east'])[0][0]-1
            n_index = np.where(m_obj.grid_north >= s_arr['rel_north'])[0][0]-1
            z_index = np.where(m_obj.res_model[n_index, e_index, :] < res_air*.9)[0][0]
            s_index = np.where(self.data_array['station']==s_arr['station'])[0][0]
            self.data_array[s_index]['elev'] = m_obj.grid_z[z_index]

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
    
    The mesh is created by first making a regular grid around the station area,
    then padding cells are added that exponentially increase to the given 
    extensions.  Depth cell increase on a log10 scale to the desired depth,
    then padding cells are added that increase exponentially.  
    
    Arguments
    -------------
        **station_object** : mtpy.modeling.modem.Stations object
                            .. seealso:: mtpy.modeling.modem.Stations
                            
    Examples
    -------------
    
    :Example 1 --> create mesh first then data file: ::
    
        >>> import mtpy.modeling.modem as modem
        >>> import os
        >>> # 1) make a list of all .edi files that will be inverted for 
        >>> edi_path = r"/home/EDI_Files"
        >>> edi_list = [os.path.join(edi_path, edi) 
                        for edi in os.listdir(edi_path) 
        >>> ...         if edi.find('.edi') > 0]
        >>> # 2) Make a Stations object
        >>> stations_obj = modem.Stations()
        >>> stations_obj.get_station_locations_from_edi(edi_list)
        >>> # 3) make a grid from the stations themselves with 200m cell spacing
        >>> mmesh = modem.Model(station_obj)
        >>> # change cell sizes
        >>> mmesh.cell_size_east = 200, 
        >>> mmesh.cell_size_north = 200
        >>> mmesh.ns_ext = 300000 # north-south extension
        >>> mmesh.ew_ext = 200000 # east-west extension of model
        >>> mmesh.make_mesh()
        >>> # check to see if the mesh is what you think it should be
        >>> msmesh.plot_mesh()
        >>> # all is good write the mesh file
        >>> msmesh.write_model_file(save_path=r"/home/modem/Inv1")
        >>> # create data file
        >>> md = modem.Data(edi_list, station_locations=mmesh.station_locations)
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")
        
    :Example 2 --> Rotate Mesh: ::
    
        >>> mmesh.mesh_rotation_angle = 60
        >>> mmesh.make_mesh()
        
    .. note:: ModEM assumes all coordinates are relative to North and East, and
             does not accommodate mesh rotations, therefore, here the rotation
             is of the stations, which essentially does the same thing.  You
             will need to rotate you data to align with the 'new' coordinate
             system.
    
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
    model_fn             full path to initial file name
    n_layers             total number of vertical layers in model
    nodes_east           relative distance between nodes in east direction 
    nodes_north          relative distance between nodes in north direction 
    nodes_z              relative distance between nodes in east direction 
    pad_east             number of cells for padding on E and W sides
                         *default* is 7
    pad_north            number of cells for padding on S and N sides
                         *default* is 7
    pad_num              number of cells with cell_size with outside of 
                         station area.  *default* is 3
    pad_stretch_h        multiplicative number for padding in horizontal
                         direction.  
    pad_stretch_v        padding cells N & S will be pad_root_north**(x) 
    pad_z                number of cells for padding at bottom
                         *default* is 4
    ew_ext               E-W extension of model in meters
    ns_ext               N-S extension of model in meters
    res_list             list of resistivity values for starting model
    res_model            starting resistivity model
    mesh_rotation_angle  Angle to rotate the grid to. Angle is measured
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
    
    def __init__(self, station_object=None, **kwargs):
        
        self.station_locations = station_object
        
        # size of cells within station area in meters
        self.cell_size_east = 500
        self.cell_size_north = 500
        
        #padding cells on either side
        self.pad_east = 7
        self.pad_north = 7
        self.pad_z = 4
        
        self.pad_num = 3
        
        self.ew_ext = 100000
        self.ns_ext = 100000
        
        #root of padding cells
        self.pad_stretch_h = 1.2
        self.pad_stretch_v = 1.2
        
        self.z1_layer = 10
        self.z_target_depth = 50000
        self.z_bottom = 300000
        
        #number of vertical layers
        self.n_layers = 30
        
        #strike angle to rotate grid to
        self.mesh_rotation_angle = 0
        
        #--> attributes to be calculated
        #grid nodes
        self._nodes_east = None
        self._nodes_north = None
        self._nodes_z = None
        
        #grid locations
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None
        
        #resistivity model
        self.res_starting_value = 100.0
        self.res_model = None
        
        #inital file stuff
        self.model_fn = None
        self.save_path = os.getcwd()
        self.model_fn_basename = 'ModEM_Model_File.rho'
        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
            self.model_fn_basename = os.path.basename(self.model_fn)
        
        self.title = 'Model File written by MTpy.modeling.modem'
        self.res_scale = 'loge'
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
    
    ### --> make nodes and grid symbiotic so if you set one the other one
    ###     gets set as well 
    ## Nodes East        
    @property
    def nodes_east(self):
        if self.grid_east is not None:
            self._nodes_east = np.array([abs(self.grid_east[ii+1]-self.grid_east[ii])
                                        for ii in range(self.grid_east.size-1)])
        return self._nodes_east
                                 
    @nodes_east.setter
    def nodes_east(self, nodes):
        nodes = np.array(nodes)
        self._nodes_east = nodes
        self.grid_east = np.array([-nodes.sum()/2+nodes[0:ii].sum()
                                   for ii in range(nodes.size)]+\
                                   [nodes.sum()/2])
    
    ## Nodes North                             
    @property
    def nodes_north(self):
        if self.grid_north is not None:
            self._nodes_north = np.array([abs(self.grid_north[ii+1]-self.grid_north[ii])
                                          for ii in range(self.grid_north.size-1)])
        return self._nodes_north
                    
    @nodes_north.setter
    def nodes_north(self, nodes):
        nodes = np.array(nodes)
        self._nodes_north = nodes
        self.grid_north = np.array([-nodes.sum()/2+nodes[0:ii].sum()
                                   for ii in range(nodes.size)]+\
                                   [nodes.sum()/2])
    
    @property
    def nodes_z(self):
        if self.grid_z is not None:
            self._nodes_z = np.array([abs(self.grid_z[ii+1]-self.grid_z[ii])
                                     for ii in range(self.grid_z.size-1)])
                                 
            return self._nodes_z
            
    @nodes_z.setter
    def nodes_z(self, nodes):
        nodes = np.array(nodes)
        self._nodes_z = nodes
        self.grid_z = np.array([nodes[0:ii].sum() for ii in range(nodes.size)]+\
                                [nodes.sum()])
        
    def get_padding_cells(self, cell_width, max_distance, num_cells, stretch):
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

    def make_mesh(self):
        """ 
        create finite element mesh according to parameters set.
        
        The mesh is built by:
            1. Making a regular grid within the station area.  
            2. Adding pad_num of cell_width cells outside of station area
            3. Adding padding cells to given extension and number of padding
               cells.
            4. Making vertical cells starting with z1_layer increasing 
               logarithmically (base 10) to z_target_depth and num_layers.
            5. Add vertical padding cells to desired extension.
            6. Check to make sure none of the stations lie on a node.
               If they do then move the node by .02*cell_width
        
        """
        
        ## --> find the edges of the grid
        ## calculate the extra width of padding cells
        ## multiply by 1.5 because this is only for 1 side
        pad_width_east = self.pad_num*1.5*self.cell_size_east
        pad_width_north = self.pad_num*1.5*self.cell_size_north
        
        ## get the extremities
        west = self.station_locations.rel_east.min()-pad_width_east
        east = self.station_locations.rel_east.max()+pad_width_east
        south = self.station_locations.rel_north.min()-pad_width_north
        north = self.station_locations.rel_north.max()+pad_width_north

        # round the numbers so they are easier to read
        west = np.round(west, -2)
        east= np.round(east, -2)
        south= np.round(south, -2)
        north = np.round(north, -2)
        
        #-------make a grid around the stations from the parameters above------
        #--> make the inner grid first
        inner_east = np.arange(west,
                               east+self.cell_size_east,
                               self.cell_size_east)
        inner_north = np.arange(south,
                                north+self.cell_size_north,
                                self.cell_size_north)
        ## compute padding cells
        padding_east = self.get_padding_cells(self.cell_size_east,
                                              self.ew_ext/2-east, 
                                              self.pad_east,
                                              self.pad_stretch_h)
        padding_north = self.get_padding_cells(self.cell_size_north,
                                               self.ns_ext/2-north, 
                                               self.pad_north,
                                               self.pad_stretch_h)
   
        # make the horizontal grid                                    
        self.grid_east = np.append(np.append(-1*padding_east[::-1]+inner_east.min(), 
                                        inner_east), 
                                    padding_east+inner_east.max())
        self.grid_north = np.append(np.append(-1*padding_north[::-1]+inner_north.min(),
                                         inner_north), 
                                    padding_north+inner_north.max())    
        
        #--> need to make sure none of the stations lie on the nodes
        for s_east in sorted(self.station_locations.rel_east):
            try:
                node_index = np.where(abs(s_east-self.grid_east) < 
                                     .02*self.cell_size_east)[0][0]
                if s_east-self.grid_east[node_index] > 0:
                    self.grid_east[node_index] -= .02*self.cell_size_east
                elif s_east-self.grid_east[node_index] < 0:
                    self.grid_east[node_index] += .02*self.cell_size_east
            except IndexError:
                continue
            
          
        #--> need to make sure none of the stations lie on the nodes
        for s_north in sorted(self.station_locations.rel_north):
            try:
                node_index = np.where(abs(s_north-self.grid_north) < 
                                     .02*self.cell_size_north)[0][0]
                if s_north-self.grid_north[node_index] > 0:
                    self.grid_north[node_index] -= .02*self.cell_size_north
                elif s_north-self.grid_north[node_index] < 0:
                    self.grid_north[node_index] += .02*self.cell_size_north
            except IndexError:
                continue

        #--> make depth grid
        log_z = np.logspace(np.log10(self.z1_layer), 
                            np.log10(self.z_target_depth-np.logspace(np.log10(self.z1_layer), 
                            np.log10(self.z_target_depth), 
                            num=self.n_layers)[-2]), 
                            num=self.n_layers-self.pad_z)

        z_nodes = np.array([np.round(zz, -int(np.floor(np.log10(zz))-1)) for zz in 
                           log_z])
                           
        #padding cells in the vertical
        z_padding = self.get_padding_cells(z_nodes[-1],
                                           self.z_bottom-z_nodes.sum(),
                                           self.pad_z,
                                           self.pad_stretch_v)
        # make the blocks into nodes as oppose to total width
        z_padding = np.array([z_padding[ii+1]-z_padding[ii] 
                              for ii in range(z_padding.size-1)])

        self.nodes_z = np.append(z_nodes, z_padding)                  

        #compute grid center
        center_east = np.round(self.grid_east.min()-self.grid_east.mean(), -1)
        center_north = np.round(self.grid_north.min()-self.grid_north.mean(), -1)
        center_z = 0
        
        # this is the value to the lower left corner from the center.
        self.grid_center = np.array([center_north, center_east, center_z])
            
        #--> print out useful information                    
        self.get_mesh_params()
            
    def get_mesh_params(self):
        #--> print out useful information                    
        print '-'*15
        print '   Number of stations = {0}'.format(len(self.station_locations.station))
        print '   Dimensions: '
        print '      e-w = {0}'.format(self.grid_east.size)
        print '      n-s = {0}'.format(self.grid_north.size)
        print '       z  = {0} (without 7 air layers)'.format(self.grid_z.size)
        print '   Extensions: '
        print '      e-w = {0:.1f} (m)'.format(self.nodes_east.__abs__().sum())
        print '      n-s = {0:.1f} (m)'.format(self.nodes_north.__abs__().sum())
        print '      0-z = {0:.1f} (m)'.format(self.nodes_z.__abs__().sum())
        
        print '  Stations rotated by: {0:.1f} deg clockwise positive from N'.format(self.mesh_rotation_angle)
        print ''        
        print ' ** Note ModEM does not accommodate mesh rotations, it assumes'
        print '    all coordinates are aligned to geographic N, E'
        print '    therefore rotating the stations will have a similar effect'
        print '    as rotating the mesh.'
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
        plt.rcParams['figure.subplot.left'] = .12
        plt.rcParams['font.size'] = 7
        
        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        plt.clf()
        
        #make a rotation matrix to rotate data
        #cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
        #sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))
        
        #turns out ModEM has not accomodated rotation of the grid, so for
        #now we will not rotate anything.
        cos_ang = 1
        sin_ang = 0
        
        #--->plot map view    
        ax1 = fig.add_subplot(1, 2, 1, aspect='equal')
        
        
        #plot station locations
        plot_east = self.station_locations.rel_east
        plot_north = self.station_locations.rel_north
        
        ax1.scatter(plot_east,
                    plot_north, 
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)
                                
        
        east_line_xlist = []
        east_line_ylist = []   
        north_min = self.grid_north.min()         
        north_max = self.grid_north.max()         
        for xx in self.grid_east:
            east_line_xlist.extend([xx*cos_ang+north_min*sin_ang, 
                                    xx*cos_ang+north_max*sin_ang])
            east_line_xlist.append(None)
            east_line_ylist.extend([-xx*sin_ang+north_min*cos_ang, 
                                    -xx*sin_ang+north_max*cos_ang])
            east_line_ylist.append(None)
        ax1.plot(east_line_xlist,
                      east_line_ylist,
                      lw=line_width,
                      color=line_color)

        north_line_xlist = []
        north_line_ylist = [] 
        east_max = self.grid_east.max()
        east_min = self.grid_east.min()
        for yy in self.grid_north:
            north_line_xlist.extend([east_min*cos_ang+yy*sin_ang,
                                     east_max*cos_ang+yy*sin_ang])
            north_line_xlist.append(None)
            north_line_ylist.extend([-east_min*sin_ang+yy*cos_ang, 
                                     -east_max*sin_ang+yy*cos_ang])
            north_line_ylist.append(None)
        ax1.plot(north_line_xlist,
                      north_line_ylist,
                      lw=line_width,
                      color=line_color)
        
        if east_limits == None:
            ax1.set_xlim(plot_east.min()-10*self.cell_size_east,
                         plot_east.max()+10*self.cell_size_east)
        else:
            ax1.set_xlim(east_limits)
        
        if north_limits == None:
            ax1.set_ylim(plot_north.min()-10*self.cell_size_north,
                         plot_north.max()+ 10*self.cell_size_east)
        else:
            ax1.set_ylim(north_limits)
            
        ax1.set_ylabel('Northing (m)', fontdict={'size':9,'weight':'bold'})
        ax1.set_xlabel('Easting (m)', fontdict={'size':9,'weight':'bold'})
        
        ##----plot depth view
        ax2 = fig.add_subplot(1, 2, 2, aspect='auto', sharex=ax1)
        

        #plot the grid 
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
        ax2.scatter(plot_east,
                    [0]*self.station_locations.station.size,
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        
        if z_limits == None:
            ax2.set_ylim(self.z_target_depth, -200)
        else:
            ax2.set_ylim(z_limits)
            
        if east_limits == None:
            ax1.set_xlim(plot_east.min()-10*self.cell_size_east,
                         plot_east.max()+10*self.cell_size_east)
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
                        Prior resistivity model. 
                        
                        .. note:: again that the modeling code 
                        assumes that the first row it reads in is the southern
                        most row and the first column it reads in is the 
                        western most column.  Similarly, the first plane it 
                        reads in is the Earth's surface.
            
            **res_starting_value** : float
                                     starting model resistivity value, 
                                     assumes a half space in Ohm-m
                                     *default* is 100 Ohm-m
                        
            **res_scale** : [ 'loge' | 'log' | 'log10' | 'linear' ]
                            scale of resistivity.  In the ModEM code it 
                            converts everything to Loge, 
                            *default* is 'loge'
                            
        """
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

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
        
        # get resistivity model        
        if self.res_model is None:
            self.res_model = np.zeros((self.nodes_north.size,
                                       self.nodes_east.size,
                                      self.nodes_z.size))
            self.res_model[:, :, :] = self.res_starting_value
            
        elif type(self.res_model) in [float, int]:
            self.res_starting_value = self.res_model
            self.res_model = np.zeros((self.nodes_north.size,
                                       self.nodes_east.size,
                                      self.nodes_z.size))
            self.res_model[:, :, :] = self.res_starting_value
        

        #--> write file
        ifid = file(self.model_fn, 'w')
        ifid.write('# {0}\n'.format(self.title.upper()))
        ifid.write('{0:>5}{1:>5}{2:>5}{3:>5} {4}\n'.format(self.nodes_north.size,
                                              self.nodes_east.size,
                                              self.nodes_z.size,
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
             self.res_scale.lower() == 'log10':
            write_res_model = np.log10(self.res_model[::-1, :, :])
        elif self.res_scale.lower() == 'linear':
            write_res_model = self.res_model[::-1, :, :]
            
        #write out the layers from resmodel
        for zz in range(self.nodes_z.size):
            ifid.write('\n')
            for ee in range(self.nodes_east.size):
                for nn in range(self.nodes_north.size):
                    ifid.write('{0:>13.5E}'.format(write_res_model[nn, ee, zz]))
                ifid.write('\n')
                
                
        if self.grid_center is None:
            #compute grid center
            center_east = -self.nodes_east.__abs__().sum()/2
            center_north = -self.nodes_north.__abs__().sum()/2
            center_z = 0
            self.grid_center = np.array([center_north, center_east, center_z])
            
        ifid.write('\n{0:>16.3f}{1:>16.3f}{2:>16.3f}\n'.format(self.grid_center[0],
                   self.grid_center[1], self.grid_center[2]))
                   
        if self.mesh_rotation_angle is None:
            ifid.write('{0:>9.3f}\n'.format(0))
        else:
            ifid.write('{0:>9.3f}\n'.format(self.mesh_rotation_angle))
        ifid.close()
        
        print 'Wrote file to: {0}'.format(self.model_fn)
        
    def read_model_file(self, model_fn=None, shift_grid=False):
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
        
        # center the grids
        if self.grid_center is None:
            self.grid_center = np.array([-self.nodes_north.sum()/2, 
                                         -self.nodes_east.sum()/2,
                                         0.0])
            
        # need to shift the grid if the center is not symmetric
        shift_north = self.grid_center[0]+self.nodes_north.sum()/2
        shift_east = self.grid_center[1]+self.nodes_east.sum()/2
        
        # shift the grid.  if shift is + then that means the center is 
        self.grid_north += shift_north
        self.grid_east += shift_east
           
        # get cell size
        self.cell_size_east = stats.mode(self.nodes_east)[0][0]
        self.cell_size_north = stats.mode(self.nodes_north)[0][0]
        
        # get number of padding cells
        self.pad_east = np.where(self.nodes_east[0:int(self.nodes_east.size/2)]
                                 != self.cell_size_east)[0][-1]
        self.north_pad = np.where(self.nodes_north[0:int(self.nodes_north.size/2)]
                                 != self.cell_size_north)[0][-1]
            
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
        
            
    def write_vtk_file(self, vtk_save_path=None,
                       vtk_fn_basename='ModEM_model_res'):
        """
        write a vtk file to view in Paraview or other
        
        Arguments:
        -------------
            **vtk_save_path** : string
                                directory to save vtk file to.  
                                *default* is Model.save_path
            **vtk_fn_basename** : string
                                  filename basename of vtk file
                                  *default* is ModEM_model_res, evtk will add
                                  on the extension .vtr
        """
        
        if vtk_save_path is None:
            vtk_fn = os.path.join(self.save_path, vtk_fn_basename)
        else:
            vtk_fn = os.path.join(vtk_save_path, vtk_fn_basename)
            
        # use cellData, this makes the grid properly as grid is n+1
        gridToVTK(vtk_fn, 
                 self.grid_north/1000., 
                 self.grid_east/1000.,
                 self.grid_z/1000.,
                 cellData={'resistivity':self.res_model}) 
        
        print '-'*50
        print '--> Wrote model file to {0}\n'.format(vtk_fn)
        print '='*26
        print '  model dimensions = {0}'.format(self.res_model.shape)
        print '     * north         {0}'.format(self.nodes_north.size)
        print '     * east          {0}'.format(self.nodes_east.size)
        print '     * depth         {0}'.format(self.nodes_z.size)
        print '='*26
        
    def get_parameters(self):
        """
        get important model parameters to write to a file for documentation 
        later.  
        
        
        """

        parameter_list = ['cell_size_east',
                          'cell_size_north',
                          'ew_ext',
                          'ns_ext',
                          'pad_east',
                          'pad_north',
                          'pad_z',
                          'pad_num',
                          'z1_layer',
                          'z_target_depth',
                          'z_bottom',
                          'mesh_rotation_angle',
                          'res_starting_value',
                          'save_path']
                          
        parameter_dict = {}
        for parameter in parameter_list:
            key = 'model.{0}'.format(parameter)
            parameter_dict[key] = getattr(self, parameter)
            
        parameter_dict['model.size'] = self.res_model.shape
        
        return parameter_dict
    
    #--> read in ascii dem file
    def read_dem_ascii(self, ascii_fn, cell_size=500, model_center=(0, 0),
                       rot_90=0):
        """
        read in dem which is ascii format
        
        The ascii format is assumed to be:
        ncols         3601
        nrows         3601
        xllcorner     -119.00013888889
        yllcorner     36.999861111111
        cellsize      0.00027777777777778
        NODATA_value  -9999
        elevation data W --> E
        N
        |
        V
        S
        """
        dfid = file(ascii_fn, 'r')
        d_dict = {}
        for ii in range(6):
            dline = dfid.readline()
            dline = dline.strip().split()
            key = dline[0].strip().lower()
            value = float(dline[1].strip())
            d_dict[key] = value
            
        x0 = d_dict['xllcorner']
        y0 = d_dict['yllcorner']
        nx = int(d_dict['ncols'])
        ny = int(d_dict['nrows'])
        cs = d_dict['cellsize']
        
        # read in the elevation data
        elevation = np.zeros((nx, ny))
        
        for ii in range(1, int(ny)+2):
            dline = dfid.readline()
            if len(str(dline)) > 1:
                #needs to be backwards because first line is the furthest north row.
                elevation[:, -ii] = np.array(dline.strip().split(' '), dtype='float')
            else:
                break
    
        dfid.close()
    
        # create lat and lon arrays from the dem fle
        lon = np.arange(x0, x0+cs*(nx), cs)
        lat = np.arange(y0, y0+cs*(ny), cs)
        
        # calculate the lower left and uper right corners of the grid in meters
        ll_en = gis_tools.project_point_ll2utm(lat[0], lon[0])
        ur_en = gis_tools.project_point_ll2utm(lat[-1], lon[-1])
        
        # estimate cell sizes for each dem measurement
        d_east = abs(ll_en[0]-ur_en[0])/nx
        d_north = abs(ll_en[1]-ur_en[1])/ny
    
        # calculate the number of new cells according to the given cell size
        # if the given cell size and cs are similar int could make the value 0,
        # hence the need to make it one if it is 0.
        num_cells = max([1, int(cell_size/np.mean([d_east, d_north]))])
    
        # make easting and northing arrays in meters corresponding to lat and lon
        east = np.arange(ll_en[0], ur_en[0], d_east)
        north = np.arange(ll_en[1], ur_en[1], d_north)
        
        #resample the data accordingly
        new_east = east[np.arange(0, east.size, num_cells)]
        new_north = north[np.arange(0, north.size, num_cells)]
        new_x, new_y = np.meshgrid(np.arange(0, east.size, num_cells),
                                   np.arange(0, north.size, num_cells),
                                   indexing='ij') 
        elevation = elevation[new_x, new_y]
        # make any null values set to minimum elevation, could be dangerous
        elevation[np.where(elevation == -9999.0)] = elevation[np.where(elevation != -9999.0)].min()
    
        # estimate the shift of the DEM to relative model coordinates
        mid_east = np.where(new_east >= model_center[0])[0][0]
        mid_north = np.where(new_north >= model_center[1])[0][0]
    
        new_east -= new_east[mid_east]
        new_north -= new_north[mid_north]
     
        # need to rotate cause I think I wrote the dem backwards
        if rot_90 == 1 or rot_90 == 3:
            elevation = np.rot90(elevation, rot_90)
            return new_north, new_east, elevation
        else:
            elevation = np.rot90(elevation, rot_90)
    
            return new_east, new_north, elevation
    
    def interpolate_elevation(self, elev_east, elev_north, elevation,
                              model_east, model_north, pad=3,
                              elevation_max=None):
        """ 
        interpolate the elevation onto the model grid.
        
        Arguments:
        ---------------
        
            **elev_east** : np.ndarray(num_east_nodes)
                          easting grid for elevation model
                          
            **elev_north** : np.ndarray(num_north_nodes)
                          northing grid for elevation model 
                          
            **elevation** : np.ndarray(num_east_nodes, num_north_nodes)
                         elevation model assumes x is east, y is north
                         Units are meters
                         
            **model_east** : np.ndarray(num_east_nodes_model)
                         relative easting grid of resistivity model 
                         
            **model_north** : np.ndarray(num_north_nodes_model)
                         relative northin grid of resistivity model 
                         
            **pad** : int
                    number of cells to repeat elevation model by.  So for pad=3,
                    then the interpolated elevation model onto the resistivity
                    model grid will have the outer 3 cells will be repeats of
                    the adjacent cell.  This is to extend the elevation model
                    to the resistivity model cause most elevation models will
                    not cover the entire area.
                    
            **elevation_max** : float
                                maximum value for elevation
                                *default* is None, which will use 
                                elevation.max()
                    
        Returns:
        --------------
        
            **interp_elev** : np.ndarray(num_north_nodes_model, num_east_nodes_model)
                            the elevation model interpolated onto the resistivity 
                            model grid.
                         
        """
        # set a maximum on the elevation, used to get rid of singular high 
        # points in the model
        if type(elevation_max) in [float, int]:
            max_find = np.where(elevation > float(elevation_max))
            elevation[max_find] = elevation_max
            
        # need to line up the elevation with the model
        grid_east, grid_north = np.broadcast_arrays(elev_east[:, None],
                                                    elev_north[None, :])
        # interpolate onto the model grid
        interp_elev = spi.griddata((grid_east.ravel(), grid_north.ravel()),
                                   elevation.ravel(),
                                   (model_east[:, None], 
                                    model_north[None, :]),
                                    method='linear',
                                    fill_value=elevation.mean())
                                    
        interp_elev[0:pad, pad:-pad] = interp_elev[pad, pad:-pad]
        interp_elev[-pad:, pad:-pad] = interp_elev[-pad-1, pad:-pad]
        interp_elev[:, 0:pad] = interp_elev[:, pad].repeat(pad).reshape(
                                                    interp_elev[:, 0:pad].shape)
        interp_elev[:, -pad:] = interp_elev[:, -pad-1].repeat(pad).reshape(
                                                    interp_elev[:, -pad:].shape)
    
        # transpose the modeled elevation to align with x=N, y=E
        interp_elev = interp_elev.T
                              
        return interp_elev   
    
    def make_elevation_model(self, interp_elev, model_nodes_z, 
                             elevation_cell=30, pad=3, res_air=1e12,
                             fill_res=100, res_sea=0.3):
        """
        Take the elevation data of the interpolated elevation model and map that
        onto the resistivity model by adding elevation cells to the existing model.
        
        ..Note: that if there are large elevation gains, the elevation cell size
                might need to be increased.
                
        Arguments:
        -------------
            **interp_elev** : np.ndarray(num_nodes_north, num_nodes_east)
                            elevation model that has been interpolated onto the
                            resistivity model grid. Units are in meters.
                            
            **model_nodes_z** : np.ndarray(num_z_nodes_of_model)
                              vertical nodes of the resistivity model without
                              topography.  Note these are the nodes given in 
                              relative thickness, not the grid, which is total
                              depth.  Units are meters.
                        
            **elevation_cell** : float
                               height of elevation cells to be added on.  These
                               are assumed to be the same at all elevations. 
                               Units are in meters
                               
            **pad** : int
                    number of cells to look for maximum and minimum elevation.
                    So if you only want elevations within the survey area, 
                    set pad equal to the number of padding cells of the 
                    resistivity model grid.
                    
            **res_air** : float
                        resistivity of air.  Default is 1E12 Ohm-m
            
            **fill_res** : float
                         resistivity value of subsurface in Ohm-m.
                    
        Returns:
        -------------
            **elevation_model** : np.ndarray(num_north_nodes, num_east_nodes, 
                                           num_elev_nodes+num_z_nodes)
                             Model grid with elevation mapped onto it. 
                             Where anything above the surface will be given the
                             value of res_air, everything else will be fill_res
                             
            **new_nodes_z** : np.ndarray(num_z_nodes+num_elev_nodes)
                            a new array of vertical nodes, where any nodes smaller
                            than elevation_cell will be set to elevation_cell.
                            This can be input into a modem.Model object to
                            rewrite the model file.
                                                 
        """
    
        # calculate the max elevation within survey area
        elev_max = interp_elev[pad:-pad, pad:-pad].max()
    	
        # need to set sea level to 0 elevation
        elev_min = max([0, interp_elev[pad:-pad, pad:-pad].min()])
        
        # scale the interpolated elevations to fit within elev_max, elev_min
        interp_elev[np.where(interp_elev > elev_max)] = elev_max
        #interp_elev[np.where(interp_elev < elev_min)] = elev_min
        
        # calculate the number of elevation cells needed
        num_elev_cells = int((elev_max-elev_min)/elevation_cell)
        print 'Number of elevation cells: {0}'.format(num_elev_cells)
        
        # find sea level if it is there
        if elev_min < 0:
            sea_level_index = num_elev_cells-abs(int((elev_min)/elevation_cell))-1
        else:
            sea_level_index = num_elev_cells-1
            
        print 'Sea level index is {0}'.format(sea_level_index)
        
        
        # make an array of just the elevation for the model
        # north is first index, east is second, vertical is third
        elevation_model = np.ones((interp_elev.shape[0],
                                   interp_elev.shape[1],
                                   num_elev_cells+model_nodes_z.shape[0]))
                                   
        elevation_model[:, :, :] = fill_res
        
        
             
        # fill in elevation model with air values.  Remeber Z is positive down, so
        # the top of the model is the highest point and index 0 is highest 
        # elevation                
        for nn in range(interp_elev.shape[0]):
            for ee in range(interp_elev.shape[1]):
                # need to test for ocean
                if interp_elev[nn, ee] < 0:
                    # fill in from bottom to sea level, then rest with air
                    elevation_model[nn, ee, 0:sea_level_index] = res_air
                    dz = sea_level_index+abs(int((interp_elev[nn, ee])/elevation_cell))+1
                    elevation_model[nn, ee, sea_level_index:dz] = res_sea
                else:
                    dz = int((elev_max-interp_elev[nn, ee])/elevation_cell)
                    elevation_model[nn, ee, 0:dz] = res_air
        
        # make new z nodes array    
        new_nodes_z = np.append(np.repeat(elevation_cell, num_elev_cells), 
                                model_nodes_z) 
                                
        new_nodes_z[np.where(new_nodes_z < elevation_cell)] = elevation_cell
        
        return elevation_model, new_nodes_z    
            
    def add_topography_to_model(self, dem_ascii_fn, write_file=True, 
                                model_center=(0,0), rot_90=0, cell_size=500, 
                                elev_cell=30, pad=1, elev_max=None):
        """
        Add topography to an existing model from a dem in ascii format.      
        
        The ascii format is assumed to be:
        ncols         3601
        nrows         3601
        xllcorner     -119.00013888889
        yllcorner     36.999861111111
        cellsize      0.00027777777777778
        NODATA_value  -9999
        elevation data W --> E
        N
        |
        V
        S
        
        Arguments:
        -------------
            **dem_ascii_fn** : string
                             full path to ascii dem file
                             
            **model_fn** : string
                         full path to existing ModEM model file
             
            **model_center** : (east, north) in meters
                             Sometimes the center of the DEM and the center of the
                             model don't line up.  Use this parameter to line 
                             everything up properly.
                             
            **rot_90** : [ 0 | 1 | 2 | 3 ]
                       rotate the elevation model by rot_90*90 degrees.  Sometimes
                       the elevation model is flipped depending on your coordinate
                       system.
                       
            **cell_size** : float (meters)
                          horizontal cell size of grid to interpolate elevation
                          onto.  This should be smaller or equal to the input
                          model cell size to be sure there is not spatial aliasing
                          
            **elev_cell** : float (meters)
                          vertical size of each elevation cell.  This value should
                          be about 1/10th the smalles skin depth.
                          
        Returns:
        ---------------
            **new_model_fn** : string
                             full path to model file that contains topography
                          
        """
         ### 1.) read in the dem and center it onto the resistivity model 
        e_east, e_north, elevation = self.read_dem_ascii(dem_ascii_fn, 
                                                         cell_size=cell_size, 
                                                         model_center=model_center, 
                                                         rot_90=rot_90)
            
        ### 2.) interpolate the elevation model onto the model grid
        m_elev = self.interpolate_elevation(e_east, e_north, elevation, 
                                            self.grid_east, self.grid_north,
                                            pad=pad, elevation_max=elev_max)
        
        m_elev[np.where(m_elev == -9999.0)] = m_elev[np.where(m_elev != -9999.0)].min()    
        ### 3.) make a resistivity model that incoorporates topography
        mod_elev, elev_nodes_z = self.make_elevation_model(m_elev, 
                                                           self.nodes_z, 
                                                           elevation_cell=elev_cell) 
        
        ### 4.) write new model file  
        self.nodes_z = elev_nodes_z
        self.res_model = mod_elev
        
        if write_file == True:
            self.save_path = os.path.dirname(self.model_fn)
            self.write_model_file(model_fn_basename='{0}_topo.rho'.format(
                                   os.path.basename(self.model_fn)[0:-4]))
                                   
            return self.model_fn


#==============================================================================
# Residuals
#==============================================================================
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
        self.residual_fn = kwargs.pop('residual_fn', None)

    
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


#==============================================================================
# Control File for inversion
#==============================================================================
class Control_Inv(object):
    """
    read and write control file for how the inversion starts and how it is run
    
    """
    
    def __init__(self, **kwargs):
        
        self.output_fn = kwargs.pop('output_fn', 'MODULAR_NLCG')
        self.lambda_initial = kwargs.pop('lambda_initial', 10)
        self.lambda_step = kwargs.pop('lambda_step', 10)
        self.model_search_step = kwargs.pop('model_search_step', 1)
        self.rms_reset_search = kwargs.pop('rms_reset_search', 2.0e-3)
        self.rms_target = kwargs.pop('rms_target', 1.05)
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
# Control File for inversion
#==============================================================================
class Control_Fwd(object):
    """
    read and write control file for 
    
    This file controls how the inversion starts and how it is run
    
    """
    
    def __init__(self, **kwargs):
        
        self.num_qmr_iter = kwargs.pop('num_qmr_iter', 40)
        self.max_num_div_calls = kwargs.pop('max_num_div_calls', 20)
        self.max_num_div_iters = kwargs.pop('max_num_div_iters', 100)
        self.misfit_tol_fwd = kwargs.pop('misfit_tol_fwd', 1.0e-7)
        self.misfit_tol_adj = kwargs.pop('misfit_tol_adj', 1.0e-7)
        self.misfit_tol_div = kwargs.pop('misfit_tol_div', 1.0e-5)

        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.fn_basename = kwargs.pop('fn_basename', 'control.fwd')
        self.control_fn = kwargs.pop('control_fn', os.path.join(self.save_path,
                                                            self.fn_basename))
                                                            
        self._control_keys = ['Number of QMR iters per divergence correction',
                              'Maximum number of divergence correction calls',
                              'Maximum number of divergence correction iters',
                              'Misfit tolerance for EM forward solver',
                              'Misfit tolerance for EM adjoint solver',
                              'Misfit tolerance for divergence correction']
        
        self._control_dict = dict([(key, value) 
                                    for key, value in zip(self._control_keys,
                                    [self.num_qmr_iter, 
                                     self.max_num_div_calls,
                                     self.max_num_div_iters,
                                     self.misfit_tol_fwd,
                                     self.misfit_tol_adj,
                                     self.misfit_tol_div])])
        self._string_fmt_dict = dict([(key, value) 
                                    for key, value in zip(self._control_keys,
                                    ['<.0f', '<.0f', '<.0f', '<.1e', '<.1e', 
                                     '<.1e'])])
                                     
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
                                    [self.num_qmr_iter, 
                                     self.max_num_div_calls,
                                     self.max_num_div_iters,
                                     self.misfit_tol_fwd,
                                     self.misfit_tol_adj,
                                     self.misfit_tol_div])])
        
        clines = []
        for key in self._control_keys:
            value = self._control_dict[key]
            str_fmt = self._string_fmt_dict[key]
            clines.append('{0:<47}: {1:{2}}\n'.format(key, value, str_fmt))
            
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
        attr_list = ['num_qmr_iter','max_num_div_calls', 'max_num_div_iters', 
                     'misfit_tol_fwd', 'misfit_tol_adj', 'misfit_tol_div']
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
        self.smoothing_east =  0.3
        self.smoothing_north =  0.3
        self.smoothing_z = 0.3
        self.smoothing_num =  1
        
        self.exception_list = []
        self.mask_arr = None
        
        self.save_path = os.getcwd()
        self.cov_fn_basename = 'covariance.cov'
        
        self.cov_fn = None
                                                  
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
    
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
            
    def write_covariance_file(self, cov_fn=None, save_path=None, 
                              cov_fn_basename=None, model_fn=None,
                              sea_water=0.3, air=1e12):
        """
        write a covariance file
        """
        
        if model_fn is not None:
            mod_obj = Model()
            mod_obj.read_model_file(model_fn)
            print 'Reading {0}'.format(model_fn)
            self.grid_dimensions = mod_obj.res_model.shape
            self.mask_arr = np.ones_like(mod_obj.res_model)
            self.mask_arr[np.where(mod_obj.res_model >= air*.9)] = 0
            self.mask_arr[np.where((mod_obj.res_model <= sea_water*1.1) & 
                              (mod_obj.res_model >= sea_water*.9))] = 9
            
        
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
        for zz in range(self.grid_dimensions[0]):
            n_smooth_line += ' {0:<5.1f}'.format(self.smoothing_north)
        clines.append(n_smooth_line+'\n')

        #--> smoothing in east direction
        e_smooth_line = ''
        for zz in range(self.grid_dimensions[1]):
            e_smooth_line += ' {0:<5.1f}'.format(self.smoothing_east)
        clines.append(e_smooth_line+'\n')
        
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
          
        # need to flip north and south.                           
        write_mask_arr = self.mask_arr[::-1, :, :].copy()
        for zz in range(self.mask_arr.shape[2]):
            clines.append(' {0:<8.0f}{0:<8.0f}\n'.format(zz+1))
            for nn in range(self.mask_arr.shape[0]):
                cline = ''
                for ee in range(self.mask_arr.shape[1]):
                    cline += '{0:^3.0f}'.format(write_mask_arr[nn, ee, zz])
                clines.append(cline+'\n')
        
        cfid = file(self.cov_fn, 'w')
        cfid.writelines(clines)
        cfid.close()
        
        print 'Wrote covariance file to {0}'.format(self.cov_fn)
        
    def read_cov_file(self, cov_fn):
        """
        read a covariance file
        """
        if not os.path.isfile(cov_fn):
            raise ModEMError('{0} not found, check path'.format(cov_fn))
        
        self.cov_fn = cov_fn
        self.save_path = os.path.dirname(self.cov_fn)
        self.cov_fn_basename = os.path.basename(self.cov_fn)        
        
        
        with open(cov_fn, 'r') as fid:
            lines = fid.readlines()
            
        num_find = False
        east_find = False
        north_find = False
        count = 0
            
        for line in lines:
            if line.find('+') >= 0 or line.find('|') >= 0:
                continue
            else:
                line_list = line.strip().split()
                if len(line_list) == 0:
                    continue
                elif len(line_list) == 1 and num_find == False and \
                     line_list[0].find('.') == -1:
                    self.smoothing_num = int(line_list[0])
                    num_find = True
                elif len(line_list) == 1 and num_find == True and \
                     line_list[0].find('.') == -1:
                     self.exceptions_num = int(line_list[0])
                elif len(line_list) == 1 and line_list[0].find('.') >= 0:
                    self.smoothing_z = float(line_list[0])
                elif len(line_list) == 3:
                    nx, ny, nz = [int(ii) for ii in line_list]
                    self.grid_dimensions = (nx, ny, nz)
                    self.mask_arr = np.ones((nx, ny, nz), dtype=np.int)
                    self.smoothing_east = np.zeros(ny)
                    self.smoothing_north = np.zeros(nx)
                    
                elif len(line_list) == 2:
                    # starts at 1 but python starts at 0
                    index_00, index_01 = [int(ii)-1 for ii in line_list]
                    
                    count = 0
                elif line_list[0].find('.') >= 0 and north_find == False:
                    self.smoothing_north = np.array(line_list, dtype=np.float)
                    north_find = True
                elif line_list[0].find('.') >= 0 and north_find == True:
                    self.smoothing_east = np.array(line_list, dtype=np.float)
                    east_find = True
                elif north_find == True and east_find == True:
                    line_list = np.array(line_list, dtype=np.int)
                    line_list = line_list.reshape((ny, 1))
                    
                    self.mask_arr[count, :, index_00:index_01+1] = line_list
                    count += 1
        
    def get_parameters(self):
        
        parameter_list = ['smoothing_north',
                          'smoothing_east',
                          'smoothing_z',
                          'smoothing_num']
                          
        parameter_dict = {}
        for parameter in parameter_list:
            key = 'covariance.{0}'.format(parameter)
            parameter_dict[key] = getattr(self, parameter)
        
        return parameter_dict
        
    def write_cov_vtk_file(self, cov_vtk_fn, model_fn=None, grid_east=None,
                           grid_north=None, grid_z=None):
        """
        write a vtk file of the covariance to match things up
        """
        
        if model_fn is not None:
            m_obj = Model()
            m_obj.read_model_file(model_fn)
            grid_east = m_obj.grid_east
            grid_north = m_obj.grid_north
            grid_z = m_obj.grid_z
        
        if grid_east is not None:
            grid_east = grid_east
        if grid_north is not None:
            grid_north = grid_north
        if grid_z is not None:
            grid_z = grid_z
            
        # use cellData, this makes the grid properly as grid is n+1
        gridToVTK(cov_vtk_fn, 
                  grid_north/1000., 
                  grid_east/1000.,
                  grid_z/1000.,
                  cellData={'covariance_mask':self.mask_arr}) 
        
        print '-'*50
        print '--> Wrote covariance file to {0}\n'.format(cov_vtk_fn)
        print '='*26
##==============================================================================
## Add in elevation to the model
##==============================================================================
#        
##--> read in ascii dem file
#def read_dem_ascii(ascii_fn, cell_size=500, model_center=(0, 0), rot_90=0):
#    """
#    read in dem which is ascii format
#    
#    The ascii format is assumed to be:
#    ncols         3601
#    nrows         3601
#    xllcorner     -119.00013888889
#    yllcorner     36.999861111111
#    cellsize      0.00027777777777778
#    NODATA_value  -9999
#    elevation data W --> E
#    N
#    |
#    V
#    S
#    """
#    dfid = file(ascii_fn, 'r')
#    d_dict = {}
#    for ii in range(6):
#        dline = dfid.readline()
#        dline = dline.strip().split()
#        key = dline[0].strip().lower()
#        value = float(dline[1].strip())
#        d_dict[key] = value
#        
#    x0 = d_dict['xllcorner']
#    y0 = d_dict['yllcorner']
#    nx = int(d_dict['ncols'])
#    ny = int(d_dict['nrows'])
#    cs = d_dict['cellsize']
#    
#    # read in the elevation data
#    elevation = np.zeros((nx, ny))
#    
#    for ii in range(1, int(ny)+2):
#        dline = dfid.readline()
#        if len(str(dline)) > 1:
#            #needs to be backwards because first line is the furthest north row.
#            elevation[:, -ii] = np.array(dline.strip().split(' '), dtype='float')
#        else:
#            break
#
#    dfid.close()
#
#    # create lat and lon arrays from the dem fle
#    lon = np.arange(x0, x0+cs*(nx), cs)
#    lat = np.arange(y0, y0+cs*(ny), cs)
#    
#    # calculate the lower left and uper right corners of the grid in meters
#    ll_en = gis_tools.project_point_ll2utm(lat[0], lon[0])
#    ur_en = gis_tools.project_point_ll2utm(lat[-1], lon[-1])
#    
#    # estimate cell sizes for each dem measurement
#    d_east = abs(ll_en[0]-ur_en[0])/nx
#    d_north = abs(ll_en[1]-ur_en[1])/ny
#
#    # calculate the number of new cells according to the given cell size
#    # if the given cell size and cs are similar int could make the value 0,
#    # hence the need to make it one if it is 0.
#    num_cells = max([1, int(cell_size/np.mean([d_east, d_north]))])
#
#    # make easting and northing arrays in meters corresponding to lat and lon
#    east = np.arange(ll_en[0], ur_en[0], d_east)
#    north = np.arange(ll_en[1], ur_en[1], d_north)
#    
#    #resample the data accordingly
#    new_east = east[np.arange(0, east.size, num_cells)]
#    new_north = north[np.arange(0, north.size, num_cells)]
#    new_x, new_y = np.meshgrid(np.arange(0, east.size, num_cells),
#                               np.arange(0, north.size, num_cells),
#                               indexing='ij') 
#    elevation = elevation[new_x, new_y]
#    # make any null values set to minimum elevation, could be dangerous
#    elevation[np.where(elevation == -9999.0)] = elevation[np.where(elevation != -9999.0)].min()
#
#    # estimate the shift of the DEM to relative model coordinates
#    mid_east = np.where(new_east >= model_center[0])[0][0]
#    mid_north = np.where(new_north >= model_center[1])[0][0]
#
#    new_east -= new_east[mid_east]
#    new_north -= new_north[mid_north]
# 
#    # need to rotate cause I think I wrote the dem backwards
#    if rot_90 == 1 or rot_90 == 3:
#        elevation = np.rot90(elevation, rot_90)
#        return new_north, new_east, elevation
#    else:
#        elevation = np.rot90(elevation, rot_90)
#
#        return new_east, new_north, elevation
#
#def interpolate_elevation(elev_east, elev_north, elevation, model_east, 
#                          model_north, pad=3):
#    """ 
#    interpolate the elevation onto the model grid.
#    
#    Arguments:
#    ---------------
#    
#        *elev_east* : np.ndarray(num_east_nodes)
#                      easting grid for elevation model
#                      
#        *elev_north* : np.ndarray(num_north_nodes)
#                      northing grid for elevation model 
#                      
#        *elevation* : np.ndarray(num_east_nodes, num_north_nodes)
#                     elevation model assumes x is east, y is north
#                     Units are meters
#                     
#        *model_east* : np.ndarray(num_east_nodes_model)
#                     relative easting grid of resistivity model 
#                     
#        *model_north* : np.ndarray(num_north_nodes_model)
#                     relative northin grid of resistivity model 
#                     
#        *pad* : int
#                number of cells to repeat elevation model by.  So for pad=3,
#                then the interpolated elevation model onto the resistivity
#                model grid will have the outer 3 cells will be repeats of
#                the adjacent cell.  This is to extend the elevation model
#                to the resistivity model cause most elevation models will
#                not cover the entire area.
#                
#    Returns:
#    --------------
#    
#        *interp_elev* : np.ndarray(num_north_nodes_model, num_east_nodes_model)
#                        the elevation model interpolated onto the resistivity 
#                        model grid.
#                     
#    """
#    # need to line up the elevation with the model
#    grid_east, grid_north = np.broadcast_arrays(elev_east[:, None],
#                                                elev_north[None, :])
#    # interpolate onto the model grid
#    interp_elev = spi.griddata((grid_east.ravel(), grid_north.ravel()),
#                               elevation.ravel(),
#                               (model_east[:, None], 
#                                model_north[None, :]),
#                                method='linear',
#                                fill_value=elevation.mean())
#                                
#    interp_elev[0:pad, pad:-pad] = interp_elev[pad, pad:-pad]
#    interp_elev[-pad:, pad:-pad] = interp_elev[-pad-1, pad:-pad]
#    interp_elev[:, 0:pad] = interp_elev[:, pad].repeat(pad).reshape(
#                                                interp_elev[:, 0:pad].shape)
#    interp_elev[:, -pad:] = interp_elev[:, -pad-1].repeat(pad).reshape(
#                                                interp_elev[:, -pad:].shape)
#
#    # transpose the modeled elevation to align with x=N, y=E
#    interp_elev = interp_elev.T
#                          
#    return interp_elev   
#
#def make_elevation_model(interp_elev, model_nodes_z, elevation_cell=30, 
#                         pad=3, res_air=1e12, fill_res=100, res_sea=0.3):
#    """
#    Take the elevation data of the interpolated elevation model and map that
#    onto the resistivity model by adding elevation cells to the existing model.
#    
#    ..Note: that if there are large elevation gains, the elevation cell size
#            might need to be increased.
#            
#    Arguments:
#    -------------
#        *interp_elev* : np.ndarray(num_nodes_north, num_nodes_east)
#                        elevation model that has been interpolated onto the
#                        resistivity model grid. Units are in meters.
#                        
#        *model_nodes_z* : np.ndarray(num_z_nodes_of_model)
#                          vertical nodes of the resistivity model without
#                          topography.  Note these are the nodes given in 
#                          relative thickness, not the grid, which is total
#                          depth.  Units are meters.
#                    
#        *elevation_cell* : float
#                           height of elevation cells to be added on.  These
#                           are assumed to be the same at all elevations. 
#                           Units are in meters
#                           
#        *pad* : int
#                number of cells to look for maximum and minimum elevation.
#                So if you only want elevations within the survey area, 
#                set pad equal to the number of padding cells of the 
#                resistivity model grid.
#                
#        *res_air* : float
#                    resistivity of air.  Default is 1E12 Ohm-m
#        
#        *fill_res* : float
#                     resistivity value of subsurface in Ohm-m.
#                
#    Returns:
#    -------------
#        *elevation_model* : np.ndarray(num_north_nodes, num_east_nodes, 
#                                       num_elev_nodes+num_z_nodes)
#                         Model grid with elevation mapped onto it. 
#                         Where anything above the surface will be given the
#                         value of res_air, everything else will be fill_res
#                         
#        *new_nodes_z* : np.ndarray(num_z_nodes+num_elev_nodes)
#                        a new array of vertical nodes, where any nodes smaller
#                        than elevation_cell will be set to elevation_cell.
#                        This can be input into a modem.Model object to
#                        rewrite the model file.
#                                             
#    """
#
#    # calculate the max elevation within survey area
#    elev_max = interp_elev[pad:-pad, pad:-pad].max()
#	
#    # need to set sea level to 0 elevation
#    elev_min = max([0, interp_elev[pad:-pad, pad:-pad].min()])
#    
#    # scale the interpolated elevations to fit within elev_max, elev_min
#    interp_elev[np.where(interp_elev > elev_max)] = elev_max
#    #interp_elev[np.where(interp_elev < elev_min)] = elev_min
#    
#    # calculate the number of elevation cells needed
#    num_elev_cells = int((elev_max-elev_min)/elevation_cell)
#    print 'Number of elevation cells: {0}'.format(num_elev_cells)
#    
#    # find sea level if it is there
#    if elev_min < 0:
#        sea_level_index = num_elev_cells-abs(int((elev_min)/elevation_cell))-1
#    else:
#        sea_level_index = num_elev_cells-1
#        
#    print 'Sea level index is {0}'.format(sea_level_index)
#    
#    
#    # make an array of just the elevation for the model
#    # north is first index, east is second, vertical is third
#    elevation_model = np.ones((interp_elev.shape[0],
#                               interp_elev.shape[1],
#                               num_elev_cells+model_nodes_z.shape[0]))
#                               
#    elevation_model[:, :, :] = fill_res
#    
#    
#         
#    # fill in elevation model with air values.  Remeber Z is positive down, so
#    # the top of the model is the highest point and index 0 is highest 
#    # elevation                
#    for nn in range(interp_elev.shape[0]):
#        for ee in range(interp_elev.shape[1]):
#            # need to test for ocean
#            if interp_elev[nn, ee] < 0:
#                # fill in from bottom to sea level, then rest with air
#                elevation_model[nn, ee, 0:sea_level_index] = res_air
#                dz = sea_level_index+abs(int((interp_elev[nn, ee])/elevation_cell))+1
#                elevation_model[nn, ee, sea_level_index:dz] = res_sea
#            else:
#                dz = int((elev_max-interp_elev[nn, ee])/elevation_cell)
#                elevation_model[nn, ee, 0:dz] = res_air
#    
#    # make new z nodes array    
#    new_nodes_z = np.append(np.repeat(elevation_cell, num_elev_cells), 
#                            model_nodes_z) 
#                            
#    new_nodes_z[np.where(new_nodes_z < elevation_cell)] = elevation_cell
#    
#    return elevation_model, new_nodes_z    
#        
#def add_topography_to_model(dem_ascii_fn, model_fn, model_center=(0,0),
#                            rot_90=0, cell_size=500, elev_cell=30, pad=1):
#    """
#    Add topography to an existing model from a dem in ascii format.      
#    
#    The ascii format is assumed to be:
#    ncols         3601
#    nrows         3601
#    xllcorner     -119.00013888889
#    yllcorner     36.999861111111
#    cellsize      0.00027777777777778
#    NODATA_value  -9999
#    elevation data W --> E
#    N
#    |
#    V
#    S
#    
#    Arguments:
#    -------------
#        *dem_ascii_fn* : string
#                         full path to ascii dem file
#                         
#        *model_fn* : string
#                     full path to existing ModEM model file
#         
#        *model_center* : (east, north) in meters
#                         Sometimes the center of the DEM and the center of the
#                         model don't line up.  Use this parameter to line 
#                         everything up properly.
#                         
#        *rot_90* : [ 0 | 1 | 2 | 3 ]
#                   rotate the elevation model by rot_90*90 degrees.  Sometimes
#                   the elevation model is flipped depending on your coordinate
#                   system.
#                   
#        *cell_size* : float (meters)
#                      horizontal cell size of grid to interpolate elevation
#                      onto.  This should be smaller or equal to the input
#                      model cell size to be sure there is not spatial aliasing
#                      
#        *elev_cell* : float (meters)
#                      vertical size of each elevation cell.  This value should
#                      be about 1/10th the smalles skin depth.
#                      
#    Returns:
#    ---------------
#        *new_model_fn* : string
#                         full path to model file that contains topography
#                      
#    """
#     ### 1.) read in the dem and center it onto the resistivity model 
#    e_east, e_north, elevation = read_dem_ascii(dem_ascii_fn, 
#                                                cell_size=cell_size, 
#                                                model_center=model_center, 
#                                                rot_90=rot_90)
#    m_obj = Model()
#    m_obj.read_model_file(model_fn)
#    ### 2.) interpolate the elevation model onto the model grid
#    m_elev = interpolate_elevation(e_east, e_north, elevation, 
#                                   m_obj.grid_east, m_obj.grid_north, pad=pad)
#    
#    m_elev[np.where(m_elev == -9999.0)] = m_elev[np.where(m_elev != -9999.0)].min()    
#    ### 3.) make a resistivity model that incoorporates topography
#    mod_elev, elev_nodes_z = make_elevation_model(m_elev, m_obj.nodes_z, 
#                                                  elevation_cell=elev_cell) 
#    
#    ### 4.) write new model file  
#    m_obj.nodes_z = elev_nodes_z
#    m_obj.res_model = mod_elev
#    m_obj.model_fn = None
#    m_obj.save_path = os.path.dirname(model_fn)
#    m_obj.write_model_file(model_fn_basename='{0}_topo.rho'.format(
#                           os.path.basename(model_fn)[0:-4]))
#                           
#    return m_obj.model_fn
#
#def change_data_elevation(data_fn, model_fn, new_data_fn=None, res_air=1e12):
#    """
#    At each station in the data file rewrite the elevation, so the station is
#    on the surface, not floating in air.
#    
#    Arguments:
#    ------------------
#        *data_fn* : string
#                    full path to a ModEM data file
#                    
#        *model_fn* : string
#                    full path to ModEM model file that has elevation 
#                    incoorporated.
#                                        
#        *new_data_fn* : string
#                        full path to new data file name.  If None, then 
#                        new file name will add _elev.dat to input filename
#                        
#        *res_air* : float
#                    resistivity of air.  Default is 1E12 Ohm-m
#    Returns:
#    -------------
#        *new_data_fn* : string
#                        full path to new data file.
#    """
#    
#    d_obj = Data()
#    d_obj.read_data_file(data_fn)
#    
#    m_obj = Model()
#    m_obj.read_model_file(model_fn)
#    
#    s_locations = d_obj.station_locations.station_locations.copy()
#    
#    # need to subtract one because we are finding the cell next to it
#    for s_arr in s_locations:
#        e_index = np.where(m_obj.grid_east >= s_arr['rel_east'])[0][0]-1
#        n_index = np.where(m_obj.grid_north >= s_arr['rel_north'])[0][0]-1
#        z_index = np.where(m_obj.res_model[n_index, e_index, :] < res_air*.9)[0][0]
#        s_index = np.where(d_obj.data_array['station']==s_arr['station'])[0][0]
#        d_obj.data_array[s_index]['elev'] = m_obj.grid_z[z_index]
#        
#        print s_arr['station'], s_arr['elev'], n_index, e_index, z_index
#        print s_arr['rel_north'], s_arr['rel_east']
#        print m_obj.grid_north[n_index], m_obj.grid_east[e_index]
#        print '-'*20
##    
##    for key in d_obj.mt_dict.keys():
##        mt_obj = d_obj.mt_dict[key]
##        e_index = np.where(m_obj.grid_east > mt_obj.grid_east)[0][0]
##        n_index = np.where(m_obj.grid_north > mt_obj.grid_north)[0][0]
##        z_index = np.where(m_obj.res_model[n_index, e_index, :] < res_air*.9)[0][0]
##        s_index = np.where(d_obj.data_array['station']==key)[0][0]        
##        d_obj.data_array[s_index]['elev'] = m_obj.grid_z[z_index]
##                
##        mt_obj.grid_elev = m_obj.grid_z[z_index] 
##        
#    if new_data_fn is None:
#        new_dfn = '{0}{1}'.format(data_fn[:-4], '_elev.dat')
#    else:
#        new_dfn=new_data_fn
#        
#    d_obj.write_data_file(save_path=os.path.dirname(new_dfn), 
#                          fn_basename=os.path.basename(new_dfn),
#                          compute_error=False,
#                          fill=False, 
#                          elevation=True)
#         
#    return new_dfn
#    
#def center_stations(data_fn, model_fn, new_data_fn=None):
#    """
#    center station locations to the middle of cells, might be useful for 
#    topography.
#    """
#    
#    d_obj = Data()
#    d_obj.read_data_file(data_fn)
#    
#    m_obj = Model()
#    m_obj.read_model_file(model_fn)
#    
#    for s_arr in d_obj.station_locations.station_locations:
#        e_index = np.where(m_obj.grid_east >= s_arr['rel_east'])[0][0]-1
#        n_index = np.where(m_obj.grid_north >= s_arr['rel_north'])[0][0]-1
#        
#        mid_east = m_obj.grid_east[e_index:e_index+2].mean()
#        mid_north = m_obj.grid_north[n_index:n_index+2].mean()
#        
#        s_index = np.where(d_obj.data_array['station']==s_arr['station'])[0][0]
#        
#        d_obj.data_array[s_index]['rel_east'] = mid_east
#        d_obj.data_array[s_index]['rel_north'] = mid_north
#        
#        print s_arr['rel_east'], s_arr['rel_north']
#        print mid_east, mid_north
#        print '-'*30
#        
#    if new_data_fn is None:
#        new_dfn = '{0}{1}'.format(data_fn[:-4], '_center.dat')
#    else:
#        new_dfn=new_data_fn
#        
#    d_obj.write_data_file(save_path=os.path.dirname(new_dfn), 
#                          fn_basename=os.path.basename(new_dfn),
#                          compute_error=False,
#                          fill=False, 
#                          elevation=True)
#         
#    return new_dfn
#    
    
#==============================================================================
# Write inversion parameters to a config type file
#==============================================================================
class ModEM_Config(object):
    """
    read and write configuration files for how each inversion is run
    """
    
    def __init__(self, **kwargs):
        self.cfg_dict = {'ModEM_Inversion_Parameters':{}}
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
        
    

    def write_config_file(self, save_dir=None, 
                          config_fn_basename='ModEM_inv.cfg'):
        """
        write a config file based on provided information
        """
        
        if save_dir is None:
            save_dir = os.getcwd()
            
        cfg_fn = os.path.join(save_dir, config_fn_basename)
        
        if self.cfg_dict is not None:
            mtcfg.write_dict_to_configfile(self.cfg_dict,
                                           cfg_fn)

            
    def add_dict(self, fn=None, obj=None):
        """
        add dictionary based on file name or object
        """
        
        if fn is not None:
            if fn.endswith('.rho'):
                m_obj = Model()
                m_obj.read_model_file(fn)
            elif fn.endswith('.dat'):
                m_obj = Data()
                m_obj.read_data_file(fn)
            elif fn.endswith('.cov'):
                m_obj = Covariance()
                m_obj.read_cov_fn(fn)
        elif obj is not None:
            m_obj = obj
        
        else:
            raise ModEMError('Need to input a file name or object')
        
        add_dict = m_obj.get_parameters()
        
        for key in add_dict.keys():
            self.cfg_dict['ModEM_Inversion_Parameters'][key] = add_dict[key]
        
            
        
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
                                             'ModEM_Model_rw.ws')
        
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
                                        1000, 5000],
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
        
        self.subplot_wspace = kwargs.pop('subplot_wspace', .3)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .8)
        self.subplot_left = kwargs.pop('subplot_left', .01)
        self.subplot_top = kwargs.pop('subplot_top', .93)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)

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
            self.station_east = md_data.station_locations.rel_east
            self.station_north = md_data.station_locations.rel_north
        
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
        # set plot properties
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        font_dict = {'size':self.font_size+2, 'weight':'bold'}
        
        #make sure there is a model to plot
        if self.res_model is None:
            self.get_model()
            
        self.cmin = np.floor(np.log10(min(self.res_list)))
        self.cmax = np.ceil(np.log10(max(self.res_list)))   
        
        #-->Plot properties
        plt.rcParams['font.size'] = self.font_size
        
        #need to add an extra row and column to east and north to make sure 
        #all is plotted see pcolor for details.
        plot_east = self.grid_east/self.dscale
        plot_north = self.grid_north/self.dscale
        
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
#        self.ax2 = mcb.make_axes(self.ax1, orientation='vertical', shrink=.35)
        self.ax2 = self.fig.add_axes([.81, .45, .16, .03])
        self.ax2.xaxis.set_ticks_position('top')
        #seg_cmap = ws.cmap_discretize(self.cmap, len(self.res_list))
        self.cb = mcb.ColorbarBase(self.ax2,cmap=self.cmap,
                                   norm=colors.Normalize(vmin=self.cmin,
                                                         vmax=self.cmax),
                                    orientation='horizontal')
                                                         
                            
        self.cb.set_label('Resistivity ($\Omega \cdot$m)',
                          fontdict={'size':self.font_size})
        self.cb.set_ticks(np.arange(self.cmin, self.cmax+1))
        self.cb.set_ticklabels([mtplottools.labeldict[cc] 
                                for cc in np.arange(self.cmin, self.cmax+1)])
                            
        #make a resistivity radio button
        #resrb = self.fig.add_axes([.85,.1,.1,.2])
        #reslabels = ['{0:.4g}'.format(res) for res in self.res_list]
        #self.radio_res = widgets.RadioButtons(resrb, reslabels, 
        #                                active=self.res_dict[self.res_value])
                                        
#        slider_ax_bounds = list(self.cb.ax.get_position().bounds)
#        slider_ax_bounds[0] += .1
        slider_ax = self.fig.add_axes([.81, .5, .16, .03])
        self.slider_res = widgets.Slider(slider_ax, 'Resistivity', 
                                         self.cmin, self.cmax, 
                                         valinit=2)
        
        
        #make a rectangular selector
        self.rect_selector = widgets.RectangleSelector(self.ax1, 
                                                       self.rect_onselect,
                                                       drawtype='box',
                                                       useblit=True)

        
        plt.show()
        
        #needs to go after show()
        self.slider_res.on_changed(self.set_res_value)
        #self.radio_res.on_clicked(self.set_res_value)


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
        
#    def set_res_value(self, label):
#        self.res_value = float(label)
#        print 'set resistivity to ', label
#        print self.res_value
    def set_res_value(self, val):
        self.res_value = 10**val
        print 'set resistivity to ', self.res_value
        
        
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
        >>> rfn = r"/home/MT/ModEM/Inv1/Test_resp_000.dat"
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
        self.ms_r = kwargs.pop('ms_r', 3)
        self.lw = kwargs.pop('lw', .5)
        self.lw_r = kwargs.pop('lw_r', 1.0)
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
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')
            
            #color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')
            
        self.phase_limits_d = kwargs.pop('phase_limits_d', None)
        self.phase_limits_od = kwargs.pop('phase_limits_od', None)
        self.res_limits_d = kwargs.pop('res_limits_d', None)
        self.res_limits_od = kwargs.pop('res_limits_od', None)
        self.tipper_limits = kwargs.pop('tipper_limits', None)

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        
        self.subplot_wspace = kwargs.pop('subplot_wspace', .3)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .85)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)
        
        self.legend_loc = 'upper center'
        self.legend_pos = (.5, 1.18)
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
            h_ratio = [1, 1, .5]
        elif self.plot_z == False:
            h_ratio = [1.5, 1, .5]
        
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
            z_obj._compute_res_phase()
            
            
            
            #find locations where points have been masked
            nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
            nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
            nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
            nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
            ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
            nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]
            
            #convert to apparent resistivity and phase
            if self.plot_z == True:
                scaling = np.zeros_like(z_obj.z)
                for ii in range(2):
                    for jj in range(2):
                        scaling[:, ii, jj] = 1./np.sqrt(z_obj.freq)
                plot_res = abs(z_obj.z.real*scaling)
                plot_res_err = abs(z_obj.z_err*scaling)
                plot_phase = abs(z_obj.z.imag*scaling)
                plot_phase_err = abs(z_obj.z_err*scaling)
                h_ratio = [1, 1, .5]
                
            elif self.plot_z == False:
                plot_res = z_obj.resistivity
                plot_res_err = z_obj.resistivity_err
                plot_phase = z_obj.phase
                plot_phase_err = z_obj.phase_err
                h_ratio = [1.5, 1, .5]

                try:
                    self.res_limits_d = (10**(np.floor(np.log10(min([plot_res[nzxx, 0, 0].min(),
                                                                     plot_res[nzyy, 1, 1].min()])))),
                                         10**(np.ceil(np.log10(max([plot_res[nzxx, 0, 0].max(),
                                                                    plot_res[nzyy, 1, 1].max()])))))
                except ValueError:
                    self.res_limits_d = None
                try:
                    self.res_limits_od = (10**(np.floor(np.log10(min([plot_res[nzxy, 0, 1].min(),
                                                                      plot_res[nzyx, 1, 0].min()])))),
                                         10**(np.ceil(np.log10(max([plot_res[nzxy, 0, 1].max(),
                                                                    plot_res[nzyx, 1, 0].max()])))))
                except ValueError:
                    self.res_limits_od = None
            
            #make figure 
            fig = plt.figure(station, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            fig.suptitle(str(station), fontdict=fontdict)
            
            #set the grid of subplots
            if np.all(t_obj.tipper == 0.0) == True:
                self.plot_tipper = False
            else:
                self.plot_tipper = True
                self.tipper_limits = (np.round(min([t_obj.tipper[ntx, 0, 0].real.min(),
                                                    t_obj.tipper[nty, 0, 1].real.min(),
                                                    t_obj.tipper[ntx, 0, 0].imag.min(),
                                                    t_obj.tipper[nty, 0, 1].imag.min()]),
                                               1),
                                       np.round(max([t_obj.tipper[ntx, 0, 0].real.max(),
                                                    t_obj.tipper[nty, 0, 1].real.max(),
                                                    t_obj.tipper[ntx, 0, 0].imag.max(),
                                                    t_obj.tipper[nty, 0, 1].imag.max()]),
                                               1))


            gs = gridspec.GridSpec(3, 4,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom, 
                               right=self.subplot_right, 
                               hspace=self.subplot_hspace,
                               height_ratios=h_ratio)
                               
            axrxx = fig.add_subplot(gs[0, 0])
            axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
            axryx = fig.add_subplot(gs[0, 2], sharex=axrxx, sharey=axrxy)
            axryy = fig.add_subplot(gs[0, 3], sharex=axrxx, sharey=axrxx)
            
            axpxx = fig.add_subplot(gs[1, 0])
            axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
            axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
            axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)
            
            axtxr = fig.add_subplot(gs[2, 0], sharex=axrxx)
            axtxi = fig.add_subplot(gs[2, 1], sharex=axrxx, sharey=axtxr)
            axtyr = fig.add_subplot(gs[2, 2], sharex=axrxx)
            axtyi = fig.add_subplot(gs[2, 3], sharex=axrxx, sharey=axtyr)
            
            self.ax_list = [axrxx, axrxy, axryx, axryy,
                            axpxx, axpxy, axpyx, axpyy,
                            axtxr, axtxi, axtyr, axtyi]

            #---------plot the apparent resistivity-----------------------------------
            #plot each component in its own subplot
            # plot data response
            erxx = mtplottools.plot_errorbar(axrxx, 
                                             period[nzxx], 
                                             plot_res[nzxx, 0, 0], 
                                             plot_res_err[nzxx, 0, 0],
                                             **kw_xx)
            erxy = mtplottools.plot_errorbar(axrxy, 
                                             period[nzxy], 
                                             plot_res[nzxy, 0, 1], 
                                             plot_res_err[nzxy, 0, 1],
                                             **kw_xx)
            eryx = mtplottools.plot_errorbar(axryx, 
                                             period[nzyx], 
                                             plot_res[nzyx, 1, 0], 
                                             plot_res_err[nzyx, 1, 0],
                                             **kw_yy)
            eryy = mtplottools.plot_errorbar(axryy, 
                                             period[nzyy], 
                                             plot_res[nzyy, 1, 1], 
                                             plot_res_err[nzyy, 1, 1],
                                             **kw_yy)
            #plot phase                         
            epxx = mtplottools.plot_errorbar(axpxx, 
                                             period[nzxx], 
                                             plot_phase[nzxx, 0, 0], 
                                             plot_phase_err[nzxx, 0, 0],
                                             **kw_xx)
            epxy = mtplottools.plot_errorbar(axpxy, 
                                             period[nzxy], 
                                             plot_phase[nzxy, 0, 1], 
                                             plot_phase_err[nzxy, 0, 1],
                                             **kw_xx)
            epyx = mtplottools.plot_errorbar(axpyx, 
                                             period[nzyx], 
                                             plot_phase[nzyx, 1, 0], 
                                             plot_phase_err[nzyx, 1, 0],
                                             **kw_yy)
            epyy = mtplottools.plot_errorbar(axpyy, 
                                             period[nzyy], 
                                             plot_phase[nzyy, 1, 1], 
                                             plot_phase_err[nzyy, 1, 1],
                                             **kw_yy)
                                          
            #plot tipper
            if self.plot_tipper == True:
                ertx = mtplottools.plot_errorbar(axtxr, 
                                                 period[ntx],
                                                 t_obj.tipper[ntx, 0, 0].real,
                                                 t_obj.tipper_err[ntx, 0, 0],
                                                 **kw_xx)
                erty = mtplottools.plot_errorbar(axtyr, 
                                                 period[nty],
                                                 t_obj.tipper[nty, 0, 1].real,
                                                 t_obj.tipper_err[nty, 0, 1],
                                                 **kw_yy)
                                         
                eptx = mtplottools.plot_errorbar(axtxi, 
                                                 period[ntx],
                                                 t_obj.tipper[ntx, 0, 0].imag,
                                                 t_obj.tipper_err[ntx, 0, 0],
                                                 **kw_xx)
                epty = mtplottools.plot_errorbar(axtyi, 
                                                 period[nty],
                                                 t_obj.tipper[nty, 0, 1].imag,
                                                 t_obj.tipper_err[nty, 0, 1],
                                                 **kw_yy)
            

                                                           
            #----------------------------------------------
            # get error bar list for editing later        
            if self.plot_tipper == False: 
                try:                   
                    self._err_list = [[erxx[1][0], erxx[1][1], erxx[2][0]],
                                      [erxy[1][0], erxy[1][1], erxy[2][0]],
                                      [eryx[1][0], eryx[1][1], eryx[2][0]],
                                      [eryy[1][0], eryy[1][1], eryy[2][0]]]
                    line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
                except IndexError:
                    print 'Found no Z components for {0}'.format(self.station)
                    line_list = [[None], [None], 
                                 [None], [None]]
                                         
                    self._err_list = [[None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None]]
    
            else:
                try:                    
                    line_list = [[erxx[0]], [erxy[0]], 
                                 [eryx[0]], [eryy[0]],
                                 [ertx[0]], [erty[0]]]
                                         
                    self._err_list = [[erxx[1][0], erxx[1][1], erxx[2][0]],
                                      [erxy[1][0], erxy[1][1], erxy[2][0]],
                                      [eryx[1][0], eryx[1][1], eryx[2][0]],
                                      [eryy[1][0], eryy[1][1], eryy[2][0]],
                                      [ertx[1][0], ertx[1][1], ertx[2][0]],
                                      [erty[1][0], erty[1][1], erty[2][0]]]
                except IndexError:
                    print 'Found no Z components for {0}'.format(station)
                    line_list = [[None], [None], 
                                 [None], [None],
                                 [None], [None]]
                                         
                    self._err_list = [[None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None]]
            #------------------------------------------
            # make things look nice        
            # set titles of the Z components
            label_list = [['$Z_{xx}$'], ['$Z_{xy}$'], 
                           ['$Z_{yx}$'], ['$Z_{yy}$']] 
            for ax, label in zip(self.ax_list[0:4], label_list):
                ax.set_title(label[0],fontdict={'size':self.font_size+2, 
                                              'weight':'bold'}) 
                                              
            # set legends for tipper components
            # fake a line
            l1 = plt.Line2D([0], [0], linewidth=0, color='w', linestyle='None', 
                            marker='.')
            t_label_list = ['Re{$T_x$}', 'Im{$T_x$}', 'Re{$T_y$}', 'Im{$T_y$}']
            label_list += [['$T_{x}$'], ['$T_{y}$']]
            for ax, label in zip(self.ax_list[-4:], t_label_list):
                ax.legend([l1], [label], loc='upper left',
                          markerscale=.01,
                          borderaxespad=.05,
                          labelspacing=.01,
                          handletextpad=.05,
                          borderpad=.05,
                          prop={'size':max([self.font_size, 6])}) 
            
   

            #set axis properties
            for aa, ax in enumerate(self.ax_list):
                ax.tick_params(axis='y', pad=self.ylabel_pad)
                
                if aa < 8:
#                    ylabels[-1] = ''
#                    ylabels[0] = ''
#                    ax.set_yticklabels(ylabels)
#                    plt.setp(ax.get_xticklabels(), visible=False)
                    if self.plot_z == True:
                        ax.set_yscale('log')
    
                else:
                    ax.set_xlabel('Period (s)', fontdict=fontdict)
                    
                if aa < 4 and self.plot_z is False:
                    ax.set_yscale('log')
                    if aa == 0 or aa == 3:
                        ax.set_ylim(self.res_limits_d)
                    elif aa == 1 or aa == 2:
                        ax.set_ylim(self.res_limits_od)

                if aa > 3 and aa < 8 and self.plot_z is False:
                    ax.yaxis.set_major_formatter(MultipleLocator(10))
                    if self.phase_limits_d is not None:
                        ax.set_ylim(self.phase_limits_d)
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
                elif aa == 8:
                    ax.set_ylabel('Tipper',
                                  fontdict=fontdict)
                        
                if aa > 7:
                    ax.yaxis.set_major_locator(MultipleLocator(.1))
                    if self.tipper_limits is not None:
                        ax.set_ylim(self.tipper_limits)
                    else:
                        pass
    
                ax.set_xscale('log')
                ax.set_xlim(xmin=10**(np.floor(np.log10(period[0])))*1.01,
                         xmax=10**(np.ceil(np.log10(period[-1])))*.99)
                ax.grid(True, alpha=.25)
                
                ylabels = ax.get_yticks().tolist()
                if aa < 8:
                    ylabels[-1] = ''
                    ylabels[0] = ''
                    ax.set_yticklabels(ylabels)
                    plt.setp(ax.get_xticklabels(), visible=False)

                
             ##----------------------------------------------
            #plot model response
            if self.resp_object is not None:
                for resp_obj in self.resp_object:
                    resp_z_obj = resp_obj.mt_dict[station].Z
                    resp_z_err = np.nan_to_num((z_obj.z-resp_z_obj.z)/z_obj.z_err)
                    resp_z_obj._compute_res_phase()
                    
                    resp_t_obj = resp_obj.mt_dict[station].Tipper
                    resp_t_err = np.nan_to_num((t_obj.tipper-resp_t_obj.tipper)/t_obj.tipper_err)
                    
                    #convert to apparent resistivity and phase
                    if self.plot_z == True:
                        scaling = np.zeros_like(resp_z_obj.z)
                        for ii in range(2):
                            for jj in range(2):
                                scaling[:, ii, jj] = 1./np.sqrt(resp_z_obj.freq)
                        r_plot_res = abs(resp_z_obj.z.real*scaling)
                        r_plot_phase = abs(resp_z_obj.z.imag*scaling)
                        
                    elif self.plot_z == False:
                        r_plot_res = resp_z_obj.resistivity
                        r_plot_phase = resp_z_obj.phase
        
                    rms_xx = resp_z_err[:, 0, 0].std()
                    rms_xy = resp_z_err[:, 0, 1].std()
                    rms_yx = resp_z_err[:, 1, 0].std()
                    rms_yy = resp_z_err[:, 1, 1].std()
                    
                    #--> make key word dictionaries for plotting
                    kw_xx = {'color':self.ctem,
                             'marker':self.mtem,
                             'ms':self.ms_r,
                             'ls':':',
                             'lw':self.lw_r,
                             'e_capsize':self.e_capsize,
                             'e_capthick':self.e_capthick}        
                   
                    kw_yy = {'color':self.ctmm,
                             'marker':self.mtmm,
                             'ms':self.ms_r,
                             'ls':':',
                             'lw':self.lw_r,
                             'e_capsize':self.e_capsize,
                             'e_capthick':self.e_capthick}
                    
                    # plot data response
                    rerxx = mtplottools.plot_errorbar(axrxx, 
                                                     period[nzxx], 
                                                     r_plot_res[nzxx, 0, 0], 
                                                     None,
                                                     **kw_xx)
                    rerxy = mtplottools.plot_errorbar(axrxy, 
                                                     period[nzxy], 
                                                     r_plot_res[nzxy, 0, 1], 
                                                     None,
                                                     **kw_xx)
                    reryx = mtplottools.plot_errorbar(axryx, 
                                                     period[nzyx], 
                                                     r_plot_res[nzyx, 1, 0], 
                                                     None,
                                                     **kw_yy)
                    reryy = mtplottools.plot_errorbar(axryy, 
                                                     period[nzyy], 
                                                     r_plot_res[nzyy, 1, 1], 
                                                     None,
                                                     **kw_yy)
                    #plot phase                         
                    repxx = mtplottools.plot_errorbar(axpxx, 
                                                     period[nzxx], 
                                                     r_plot_phase[nzxx, 0, 0], 
                                                     None,
                                                     **kw_xx)
                    repxy = mtplottools.plot_errorbar(axpxy, 
                                                     period[nzxy], 
                                                     r_plot_phase[nzxy, 0, 1], 
                                                     None,
                                                     **kw_xx)
                    repyx = mtplottools.plot_errorbar(axpyx, 
                                                     period[nzyx], 
                                                     r_plot_phase[nzyx, 1, 0], 
                                                     None,
                                                     **kw_yy)
                    repyy = mtplottools.plot_errorbar(axpyy, 
                                                     period[nzyy], 
                                                     r_plot_phase[nzyy, 1, 1], 
                                                     None,
                                                     **kw_yy)
                                                  
                    #plot tipper
                    if self.plot_tipper == True:
                        rertx = mtplottools.plot_errorbar(axtxr, 
                                                         period[ntx],
                                                         resp_t_obj.tipper[ntx, 0, 0].real,
                                                         None,
                                                         **kw_xx)
                        rerty = mtplottools.plot_errorbar(axtyr, 
                                                         period[nty],
                                                         resp_t_obj.tipper[nty, 0, 1].real,
                                                         None,
                                                         **kw_yy)
                                                 
                        reptx = mtplottools.plot_errorbar(axtxi, 
                                                         period[ntx],
                                                         resp_t_obj.tipper[ntx, 0, 0].imag,
                                                         None,
                                                         **kw_xx)
                        repty = mtplottools.plot_errorbar(axtyi, 
                                                         period[nty],
                                                         resp_t_obj.tipper[nty, 0, 1].imag,
                                                         None,
                                                         **kw_yy)
                                     
                    if self.plot_tipper == False:
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
                        label_list[4] += ['$T^m_{x}$ '+
                                        'rms={0:.2f}'.format(resp_t_err[:, 0, 0].std())]
                        label_list[5] += ['$T^m_{y}$'+
                                        'rms={0:.2f}'.format(resp_t_err[:, 0, 1].std())]
                
                legend_ax_list = self.ax_list[0:4]
#                if self.plot_tipper == True:
#                    legend_ax_list += [self.ax_list[-4], self.ax_list[-2]]
                    
                for aa, ax in enumerate(legend_ax_list):
                    ax.legend(line_list[aa],
                              label_list[aa],
                              loc=self.legend_loc,
                              bbox_to_anchor=self.legend_pos,
                              markerscale=self.legend_marker_scale,
                              borderaxespad=self.legend_border_axes_pad,
                              labelspacing=self.legend_label_spacing,
                              handletextpad=self.legend_handle_text_pad,
                              borderpad=self.legend_border_pad,
                              prop={'size':max([self.font_size, 5])})
            
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

        fig = plt.gcf()
        if fig_dpi == None:
            fig_dpi = self.fig_dpi
            
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')
            
        else:
            save_fn = os.path.join(save_fn, '_L2.'+
                                    file_format)
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')
        
        if close_fig == 'y':
            plt.clf()
            plt.close(fig)
        
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
                    self.plot_period_list = [self.data_obj.period_list[ii]
                                             for ii in self.plot_period_list]
                else:
                    pass
            elif type(self.plot_period_list) is int:
                self.plot_period_list = self.data_obj.period_list[self.plot_period_list]
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
            east_min = self.data_obj.data_array['rel_east'].min()-\
                                                            self.pad_east
            east_max = self.data_obj.data_array['rel_east'].max()+\
                                                            self.pad_east
            self.ew_limits = (east_min/self.dscale, east_max/self.dscale)
            
        if self.ns_limits == None:
            north_min = self.data_obj.data_array['rel_north'].min()-\
                                                            self.pad_north
            north_max = self.data_obj.data_array['rel_north'].max()+\
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
                approx_depth, d_index = ws.estimate_skin_depth(self.model_obj.res_model.copy(),
                                                            self.model_obj.grid_z.copy()/self.dscale, 
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
                self.station_east = md_data.station_locations.rel_east/self.dscale
                self.station_north = md_data.station_locations.rel_north/self.dscale
                self.station_elev = md_data.station_locations.elev/self.dscale                
                self.station_names = md_data.station_locations.station
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
#    """
#    plot all slices and be able to scroll through the model
#    
#    :Example: ::
#    
#        >>> import mtpy.modeling.modem as modem
#        >>> mfn = r"/home/modem/Inv1/Modular_NLCG_100.rho"
#        >>> dfn = r"/home/modem/Inv1/ModEM_data.dat"       
#        >>> pds = ws.PlotSlices(model_fn=mfn, data_fn=dfn)
#        
#    ======================= ===================================================
#    Buttons                  Description    
#    ======================= ===================================================
#    'e'                     moves n-s slice east by one model block
#    'w'                     moves n-s slice west by one model block
#    'n'                     moves e-w slice north by one model block
#    'm'                     moves e-w slice south by one model block
#    'd'                     moves depth slice down by one model block
#    'u'                     moves depth slice up by one model block
#    ======================= ===================================================
#
#    
#    ======================= ===================================================
#    Attributes              Description    
#    ======================= ===================================================
#    ax_en                   matplotlib.axes instance for depth slice  map view 
#    ax_ez                   matplotlib.axes instance for e-w slice
#    ax_map                  matplotlib.axes instance for location map
#    ax_nz                   matplotlib.axes instance for n-s slice
#    climits                 (min , max) color limits on resistivity in log 
#                            scale. *default* is (0, 4)
#    cmap                    name of color map for resisitiviy.
#                            *default* is 'jet_r'
#    data_fn                 full path to data file name
#    dscale                  scaling parameter depending on map_scale
#    east_line_xlist         list of line nodes of east grid for faster plotting
#    east_line_ylist         list of line nodes of east grid for faster plotting
#    ew_limits               (min, max) limits of e-w in map_scale units
#                            *default* is None and scales to station area
#    fig                     matplotlib.figure instance for figure
#    fig_aspect              aspect ratio of plots. *default* is 1
#    fig_dpi                 resolution of figure in dots-per-inch
#                            *default* is 300
#    fig_num                 figure instance number
#    fig_size                [width, height] of figure window. 
#                            *default* is [6,6]
#    font_dict               dictionary of font keywords, internally created
#    font_size               size of ticklables in points, axes labes are 
#                            font_size+2. *default* is 7
#    grid_east               relative location of grid nodes in e-w direction
#                            in map_scale units
#    grid_north              relative location of grid nodes in n-s direction
#                            in map_scale units
#    grid_z                  relative location of grid nodes in z direction
#                            in map_scale units
#    index_east              index value of grid_east being plotted
#    index_north             index value of grid_north being plotted
#    index_vertical          index value of grid_z being plotted
#    initial_fn              full path to initial file
#    key_press               matplotlib.canvas.connect instance
#    map_scale               [ 'm' | 'km' ] scale of map. *default* is km
#    mesh_east               np.meshgrid(grid_east, grid_north)[0]
#    mesh_en_east            np.meshgrid(grid_east, grid_north)[0]
#    mesh_en_north           np.meshgrid(grid_east, grid_north)[1]
#    mesh_ez_east            np.meshgrid(grid_east, grid_z)[0]
#    mesh_ez_vertical        np.meshgrid(grid_east, grid_z)[1]
#    mesh_north              np.meshgrid(grid_east, grid_north)[1]
#    mesh_nz_north           np.meshgrid(grid_north, grid_z)[0]
#    mesh_nz_vertical        np.meshgrid(grid_north, grid_z)[1]
#    model_fn                full path to model file
#    ms                      size of station markers in points. *default* is 2
#    nodes_east              relative distance betwen nodes in e-w direction
#                            in map_scale units
#    nodes_north             relative distance betwen nodes in n-s direction
#                            in map_scale units
#    nodes_z                 relative distance betwen nodes in z direction
#                            in map_scale units
#    north_line_xlist        list of line nodes north grid for faster plotting  
#    north_line_ylist        list of line nodes north grid for faster plotting
#    ns_limits               (min, max) limits of plots in n-s direction
#                            *default* is None, set veiwing area to station area 
#    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
#                            *default* is 'y'
#    res_model               np.ndarray(n_north, n_east, n_vertical) of 
#                            model resistivity values in linear scale           
#    station_color           color of station marker. *default* is black
#    station_dict_east       location of stations for each east grid row
#    station_dict_north      location of stations for each north grid row
#    station_east            location of stations in east direction
#    station_fn              full path to station file 
#    station_font_color      color of station label 
#    station_font_pad        padding between station marker and label
#    station_font_rotation   angle of station label
#    station_font_size       font size of station label
#    station_font_weight     weight of font for station label
#    station_id              [min, max] index values for station labels
#    station_marker          station marker
#    station_names           name of stations
#    station_north           location of stations in north direction
#    subplot_bottom          distance between axes and bottom of figure window
#    subplot_hspace          distance between subplots in vertical direction
#    subplot_left            distance between axes and left of figure window  
#    subplot_right           distance between axes and right of figure window
#    subplot_top             distance between axes and top of figure window
#    subplot_wspace          distance between subplots in horizontal direction
#    title                   title of plot 
#    z_limits                (min, max) limits in vertical direction,
#    ======================= ===================================================
#    
#    """
    
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
                self.station_east = md_data.station_locations.rel_east/self.dscale
                self.station_north = md_data.station_locations.rel_north/self.dscale
                self.station_names = md_data.station_locations.station
                self.station_elev = md_data.station_locations.elev/self.dscale
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
            if self.index_north == self.grid_north.size:
                print 'Already at northern most grid cell'
            else:
                self.index_north += 1
                if self.index_north > self.grid_north.size:
                    self.index_north = self.grid_north.size
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
            if self.index_east == self.grid_east.size:
                print 'Already at eastern most grid cell'
            else:
                self.index_east += 1
                if self.index_east > self.grid_east.size:
                    self.index_east = self.grid_east.size
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
            if self.index_vertical == self.grid_z.size:
                print 'Already at deepest grid cell'
            else:
                self.index_vertical += 1
                if self.index_vertical > self.grid_z.size:
                    self.index_vertical = self.grid_z.size
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
            for ee, nn, elev, name in zip(self.station_east, 
                                          self.station_north,
                                          self.station_elev,
                                          self.station_names):
                if elev <= self.grid_z[self.index_vertical]:                
                    self.ax_en.text(ee, nn, '+', 
                                    verticalalignment='center',
                                    horizontalalignment='center',
                                    fontdict={'size':7, 'weight':'bold',
                                              'color':(.75, 0, 0)})
                    self.ax_en.text(ee, nn, name[2:], 
                                    verticalalignment='center',
                                    horizontalalignment='center',
                                    fontdict={'size':7, 'weight':'bold',
                                              'color':(.75, 0, 0)})

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
# plot rms maps
#==============================================================================
class Plot_RMS_Maps(object):
    """
    plots the RMS as (data-model)/(error) in map view for all components
    of the data file.  Gets this infomration from the .res file output
    by ModEM.
    
    Arguments:
    ------------------
    
        **residual_fn** : string
                          full path to .res file
                          
    =================== =======================================================
    Attributes                   Description    
    =================== =======================================================
    fig                 matplotlib.figure instance for a single plot                       
    fig_dpi             dots-per-inch resolution of figure *default* is 200
    fig_num             number of fig instance *default* is 1
    fig_size            size of figure in inches [width, height] 
                        *default* is [7,6]
    font_size           font size of tick labels, axis labels are +2
                        *default* is 8 
    marker              marker style for station rms, 
                        see matplotlib.line for options,
                        *default* is 's' --> square
    marker_size         size of marker in points. *default* is 10
    pad_x               padding in map units from edge of the axis to stations
                        at the extremeties in longitude. 
                        *default* is 1/2 tick_locator
    pad_y               padding in map units from edge of the axis to stations
                        at the extremeties in latitude. 
                        *default* is 1/2 tick_locator 
    period_index        index of the period you want to plot according to 
                        self.residual.period_list. *default* is 1
    plot_yn             [ 'y' | 'n' ] default is 'y' to plot on instantiation
    plot_z_list         internal variable for plotting
    residual            modem.Data instance that holds all the information 
                        from the residual_fn given
    residual_fn         full path to .res file
    rms_cmap            matplotlib.cm object for coloring the markers
    rms_cmap_dict       dictionary of color values for rms_cmap 
    rms_max             maximum rms to plot. *default* is 5.0
    rms_min             minimum rms to plot. *default* is 1.0
    save_path           path to save figures to. *default* is directory of 
                        residual_fn
    subplot_bottom      spacing from axis to bottom of figure canvas.
                        *default* is .1
    subplot_hspace      horizontal spacing between subplots.
                        *default* is .1
    subplot_left        spacing from axis to left of figure canvas.
                        *default* is .1
    subplot_right       spacing from axis to right of figure canvas.
                        *default* is .9
    subplot_top         spacing from axis to top of figure canvas.
                        *default* is .95
    subplot_vspace      vertical spacing between subplots.
                        *default* is .01
    tick_locator        increment for x and y major ticks. *default* is 
                        limits/5  
    =================== =======================================================
    
    =================== =======================================================
    Methods             Description    
    =================== =======================================================
    plot                plot rms maps for a single period 
    plot_loop           loop over all frequencies and save figures to save_path
    read_residual_fn    read in residual_fn
    redraw_plot         after updating attributes call redraw_plot to 
                        well redraw the plot
    save_figure         save the figure to a file
    =================== =======================================================
    
    
    :Example: ::
    
        >>> import mtpy.modeling.modem as modem
        >>> rms_plot = Plot_RMS_Maps(r"/home/ModEM/Inv1/mb_NLCG_030.res")
        >>> # change some attributes
        >>> rms_plot.fig_size = [6, 4]
        >>> rms_plot.rms_max = 3
        >>> rms_plot.redraw_plot()
        >>> # happy with the look now loop over all periods
        >>> rms_plot.plot_loop()
    """
    
    def __init__(self, residual_fn, **kwargs):
        self.residual_fn = residual_fn
        self.residual = None
        self.save_path = kwargs.pop('save_path', os.path.dirname(self.residual_fn))

        self.period_index = kwargs.pop('period_index', 0)        
        
        self.subplot_left = kwargs.pop('subplot_left', .1)
        self.subplot_right = kwargs.pop('subplot_right', .9)
        self.subplot_top = kwargs.pop('subplot_top', .95)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .1)
        self.subplot_vspace = kwargs.pop('subplot_vspace', .01)

        self.font_size = kwargs.pop('font_size', 8)
        
        self.fig_size = kwargs.pop('fig_size', [7.75, 6.75])
        self.fig_dpi = kwargs.pop('fig_dpi', 200)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig = None
        
        self.marker = kwargs.pop('marker', 's')
        self.marker_size = kwargs.pop('marker_size', 10)
        
        
        self.rms_max = kwargs.pop('rms_max', 5)
        self.rms_min = kwargs.pop('rms_min', 0)
        
        self.tick_locator = kwargs.pop('tick_locator', None)
        self.pad_x = kwargs.pop('pad_x', None)
        self.pad_y = kwargs.pop('pad_y', None)
        
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        
        # colormap for rms, goes white to black from 0 to rms max and 
        # red below 1 to show where the data is being over fit
        self.rms_cmap_dict = {'red':((0.0, 1.0, 1.0), 
                                    (0.2, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
                             'green':((0.0, 0.0, 0.0), 
                                      (0.2, 1.0, 1.0),
                                      (1.0, 0.0, 0.0)),
                             'blue':((0.0, 0.0, 0.0),
                                     (0.2, 1.0, 1.0),
                                     (1.0, 0.0, 0.0))}
                      
        self.rms_cmap = colors.LinearSegmentedColormap('rms_cmap', 
                                                       self.rms_cmap_dict, 
                                                       256)
                                                       
        self.plot_z_list = [{'label':r'$Z_{xx}$', 'index':(0, 0), 'plot_num':1},
                           {'label':r'$Z_{xy}$', 'index':(0, 1), 'plot_num':2},
                           {'label':r'$Z_{yx}$', 'index':(1, 0), 'plot_num':3},
                           {'label':r'$Z_{yy}$', 'index':(1, 1), 'plot_num':4},
                           {'label':r'$T_{x}$', 'index':(0, 0), 'plot_num':5},
                           {'label':r'$T_{y}$', 'index':(0, 1), 'plot_num':6}]
                           
                           
        if self.plot_yn == 'y':
            self.plot()
            
    def read_residual_fn(self):
        if self.residual is None:
            self.residual = Data()
            self.residual.read_data_file(self.residual_fn)
        else:
            pass
        
    def plot(self):
        """
        plot rms in map view
        """

        self.read_residual_fn()

        font_dict = {'size':self.font_size+2, 'weight':'bold'}
        rms_1 = 1./self.rms_max
        
        if self.tick_locator is None:
            x_locator = np.round((self.residual.data_array['lon'].max()-
                                    self.residual.data_array['lon'].min())/5, 2)
            y_locator = np.round((self.residual.data_array['lat'].max()-
                                    self.residual.data_array['lat'].min())/5, 2)
                                    
            if x_locator > y_locator:
                self.tick_locator = x_locator
            
            elif x_locator < y_locator:
                self.tick_locator = y_locator
                
            
        if self.pad_x is None:
            self.pad_x = self.tick_locator/2
        if self.pad_y is None:
            self.pad_y = self.tick_locator/2
        
        
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        plt.rcParams['figure.subplot.wspace'] = self.subplot_hspace
        plt.rcParams['figure.subplot.hspace'] = self.subplot_vspace
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        
        for p_dict in self.plot_z_list:
            ax = self.fig.add_subplot(3, 2, p_dict['plot_num'], aspect='equal')
            
            ii = p_dict['index'][0]
            jj = p_dict['index'][0]
            
            for r_arr in self.residual.data_array:
                # calulate the rms self.residual/error
                if p_dict['plot_num'] < 5:
                    rms = r_arr['z'][self.period_index, ii, jj].__abs__()/\
                                (r_arr['z_err'][self.period_index, ii, jj].real)
                    
                else: 
                    rms = r_arr['tip'][self.period_index, ii, jj].__abs__()/\
                                (r_arr['tip_err'][self.period_index, ii, jj].real)
        
                #color appropriately
                if np.nan_to_num(rms) == 0.0:
                    marker_color = (1, 1, 1)
                    marker = '.'
                    marker_size = .1
                    marker_edge_color = (1, 1, 1)
                if rms > self.rms_max:
                    marker_color = (0, 0, 0)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)
                    
                elif rms >= 1 and rms <= self.rms_max:
                    r_color = 1-rms/self.rms_max+rms_1
                    marker_color = (r_color, r_color, r_color)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)
                    
                elif rms < 1:
                    r_color = 1-rms/self.rms_max
                    marker_color = (1, r_color, r_color)
                    marker = self.marker
                    marker_size = self.marker_size
                    marker_edge_color = (0, 0, 0)
                    
                ax.plot(r_arr['lon'], r_arr['lat'], 
                        marker=marker,
                        ms=marker_size,
                        mec=marker_edge_color,
                        mfc=marker_color,
                        zorder=3)
            
            if p_dict['plot_num'] == 1 or p_dict['plot_num'] == 3:
                ax.set_ylabel('Latitude (deg)', fontdict=font_dict)
                plt.setp(ax.get_xticklabels(), visible=False)
                
            elif p_dict['plot_num'] == 2 or p_dict['plot_num'] == 4:
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                
            elif p_dict['plot_num'] == 6:
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.set_xlabel('Longitude (deg)', fontdict=font_dict)
                
            else:
                ax.set_xlabel('Longitude (deg)', fontdict=font_dict)
                ax.set_ylabel('Latitude (deg)', fontdict=font_dict)
                
            ax.text(self.residual.data_array['lon'].min()+.005-self.pad_x, 
                    self.residual.data_array['lat'].max()-.005+self.pad_y,
                    p_dict['label'],
                    verticalalignment='top',
                    horizontalalignment='left',
                    bbox={'facecolor':'white'},
                    zorder=3)
                    
            ax.tick_params(direction='out')
            ax.grid(zorder=0, color=(.75, .75, .75))
            
            #[line.set_zorder(3) for line in ax.lines]
            
            ax.set_xlim(self.residual.data_array['lon'].min()-self.pad_x, 
                        self.residual.data_array['lon'].max()+self.pad_x)
                        
            ax.set_ylim(self.residual.data_array['lat'].min()-self.pad_y, 
                        self.residual.data_array['lat'].max()+self.pad_y)
            
            ax.xaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.yaxis.set_major_locator(MultipleLocator(self.tick_locator))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%2.2f'))
            
        
        
        #cb_ax = mcb.make_axes(ax, orientation='vertical', fraction=.1)
        cb_ax = self.fig.add_axes([self.subplot_right+.02, .225, .02, .45])
        color_bar = mcb.ColorbarBase(cb_ax, 
                                     cmap=self.rms_cmap, 
                                     norm=colors.Normalize(vmin=self.rms_min, 
                                                           vmax=self.rms_max),
                                     orientation='vertical')
        
        color_bar.set_label('RMS', fontdict=font_dict)
        
        self.fig.suptitle('period = {0:.5g} (s)'.format(self.residual.period_list[self.period_index]), 
                     fontdict={'size':self.font_size+3, 'weight':'bold'})
        plt.show()
        
    def redraw_plot(self):
        plt.close('all')
        self.plot()

    def save_figure(self, save_path=None, save_fn_basename=None, 
                    save_fig_dpi=None, fig_format='.png', fig_close=True):
        """
        save figure in the desired format
        """
        if save_path is not None:
            self.save_path = save_path
        
        if save_fn_basename is not None:
            pass
        else:
            save_fn_basename = '{0:02}_RMS_{1:.5g}_s.{2}'.format(self.period_index,
                                self.residual.period_list[self.period_index],
                                fig_format)
        save_fn = os.path.join(self.save_path, save_fn_basename) 
        
        if save_fig_dpi is not None:
            self.fig_dpi = save_fig_dpi
        
        self.fig.savefig(save_fn,  dpi=self.fig_dpi)
        print 'saved file to {0}'.format(save_fn)
                    
        if fig_close == True:
            plt.close('all')
            
    def plot_loop(self, fig_format='png'):
        """
        loop over all periods and save figures accordingly
        """
        self.read_residual_fn()
        
        for f_index in range(self.residual.period_list.size):
            self.period_index = f_index
            self.plot()
            self.save_figure(fig_format=fig_format)
            


#==============================================================================
# Exceptions
#==============================================================================
class ModEMError(Exception):
    pass