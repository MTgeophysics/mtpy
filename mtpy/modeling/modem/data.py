
"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
from __future__ import print_function
import os
import sys
import csv
import numpy as np
from logging import INFO as My_Log_Level  # this module's log level

import mtpy.analysis.pt as pt
from mtpy.core import mt as mt
from mtpy.core import z as mtz
from mtpy.modeling import ws3dinv as ws
from mtpy.utils import gis_tools as gis_tools
from mtpy.utils.mtpy_decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog

from mtpy.modeling.modem.exception import ModEMError, DataError
from mtpy.modeling.modem.station import Stations
from mtpy.modeling.modem.model import Model

try:
    from pyevtk.hl import pointsToVTK
except ImportError:
    print('If you want to write a vtk file for 3d viewing, you need download '
          'and install evtk from https://bitbucket.org/pauloh/pyevtk', file=sys.stderr)

    print('Note: if you are using Windows you should build evtk first with'
          'either MinGW or cygwin using the command: \n'
          '    python setup.py build -compiler=mingw32  or \n'
          '    python setup.py build -compiler=cygwin')


# =============================================================================
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
    Attributes              Description
    ====================== ====================================================
    _dtype                 internal variable defining the data type of
                           data_array
    _logger                python logging object that put messages in logging
                           format defined in logging configure file, see MtPyLog
                           for more information
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
    error_type_z           [ 'egbert' | 'mean_od' | 'eigen' | 'median']
                           *default* is 'egbert_floor'
                                * add '_floor' to any of the above to set the
                                  error as an error floor, otherwise all
                                  components are give weighted the same

                                * 'egbert'  sets error to
                                            error_value_z * sqrt(abs(zxy*zyx))
                                * 'mean_od' sets error to
                                            error_value_z * mean([Zxy, Zyx])
                                            (non zeros)
                                * 'eigen'   sets error to
                                            error_value_z * eigenvalues(Z[ii])
                                * 'median'  sets error to
                                            error_value_z * median([Zxx, Zxy, Zyx, Zyy])
                                            (non zeros)
                           A 2x2 numpy array of error_type_z can be specified to
                           explicitly set the error_type_z for each component.

    error_value_z          percentage to multiply Z by to set error
                           *default* is 5 for 5% of Z as error
                           A 2x2 numpy array of values can be specified to
                           explicitly set the error_value_z for each component.

    error_value_tipper     absolute error between 0 and 1.
    fn_basename            basename of data file. *default* is 'ModEM_Data.dat'
    formatting             ['1' | '2'], format of the output data file, *default* is '1'
    header_strings         strings for header of data file following the format
                           outlined in the ModEM documentation
    inv_comp_dict          dictionary of inversion components
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
    model_epsg             epsg code for model projection, provide this to
                           project model to non-utm coordinates. Find the epsg
                           code for your projection on
                           http://spatialreference.org/ref/ or google search
                           epsg "your projection"
    model_utm_zone         alternative to model_epsg, choose a utm zone to
                           project all sites to (e.g. '55S')
    mt_dict                dictionary of mtpy.core.mt.MT objects with keys
                           being station names
    period_buffer          float or int
                           if specified, apply a buffer so that interpolation doesn't
                           stretch too far over periods
    period_dict            dictionary of period index for period_list
    period_list            list of periods to invert for
    period_max             maximum value of period to invert for
    period_min             minimum value of period to invert for
    period_buffer          buffer so that interpolation doesn't stretch too far
                              over periods. Provide a float or integer factor, 
                              greater than which interpolation will not stretch.
                              e.g. 1.5 means only interpolate to a maximum of
                              1.5 times each side of each frequency value
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
    center_stations            Center station locations to the middle of cells,
                               might be useful for topography.
    change_data_elevation      At each station in the data file rewrite the
                               elevation, so the station is on the surface,
                               not floating in air.
    compute_inv_error          compute the error from the given parameters
    convert_modem_to_ws        convert a ModEM data file to WS format.
    convert_ws3dinv_data_file  convert a ws3dinv file to ModEM fomrat,
                               **Note** this doesn't include tipper data and
                               you need a station location file like the one
                               output by mtpy.modeling.ws3dinv
    fill_data_array            fill the data array from mt_dict
    filter_periods             Select the periods of the mt_obj that are in
                               per_array. used to do original freq inversion.
    get_header_string          reset the header sring for file
    get_mt_dict                get mt_dict from edi file list
    get_parameters             get important parameters for documentation
    get_period_list            make a period list to invert for
    get_relative_station_locations     get station locations from edi files
    project_stations_on_topography     This method is used in add_topography().
                                       It will Re-write the data file to change
                                       the elevation column. And update
                                       covariance mask according topo elevation
                                       model.
    read_data_file             read in a ModEM data file and fill attributes
                               data_array, station_locations, period_list,
                               mt_dict
    write_data_file            write a ModEM data file
    write_vtk_station_file     write a vtk file for station locations.  For now
                               this in relative coordinates.
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
        >>> import mtpy.core.mt
        >>> import mtpy.modeling.modem as modem
        >>> edi_path = r"/home/mt/edi_files"
        >>> edi_list = [os.path.join(edi_path, edi) \
                        for edi in os.listdir(edi_path)\
                        if edi.find('.edi') > 0]
        >>> md = modem.Data(edi_list)
        >>> #get period list from an .edi file
        >>> mt_obj1 = mt.MT(edi_list[0])
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

        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)
        self._logger.setLevel(My_Log_Level)
        
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
        # self.period_step = 1
        self.period_min = None
        self.period_max = None
        self.period_buffer = None
        self.max_num_periods = None
        self.data_period_list = None
        
        self.data_fn = 'ModEM_Data.dat'
        self.save_path = os.getcwd()
        self.fn_basename = None

        self.formatting = '1'

        self._rotation_angle = 0.0
        

        self.center_point = None

        self.data_array = None
        self.mt_dict = None
        self.model_utm_zone = None
        self.model_epsg = None

        self._z_shape = (1, 2, 2)
        self._t_shape = (1, 1, 2)
        self._dtype = [('station', '|U10'),
                       ('lat', np.float),
                       ('lon', np.float),
                       ('elev', np.float),
                       ('rel_east', np.float),
                       ('rel_north', np.float),
                       ('rel_elev', np.float),
                       ('east', np.float),
                       ('north', np.float),
                       ('zone', '|S4'),
                       ('z', (np.complex, self._z_shape)),
                       ('z_err', (np.float, self._z_shape)),
                       ('z_inv_err', (np.float, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.float, self._t_shape)),
                       ('tip_inv_err', (np.float, self._t_shape))]
        
        self.inv_mode_dict = {'1': ['Full_Impedance', 'Full_Vertical_Components'],
                              '2': ['Full_Impedance'],
                              '3': ['Off_Diagonal_Impedance',
                                    'Full_Vertical_Components'],
                              '4': ['Off_Diagonal_Impedance'],
                              '5': ['Full_Vertical_Components'],
                              '6': ['Full_Interstation_TF'],
                              '7': ['Off_Diagonal_Rho_Phase']}
        self.inv_comp_dict = {'Full_Impedance': ['zxx', 'zxy', 'zyx', 'zyy'],
                              'Off_Diagonal_Impedance': ['zxy', 'zyx'],
                              'Full_Vertical_Components': ['tx', 'ty']}

        self.comp_index_dict = {'zxx': (0, 0), 'zxy': (0, 1), 'zyx': (1, 0),
                                'zyy': (1, 1), 'tx': (0, 0), 'ty': (0, 1)}

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
        
        for key in list(kwargs.keys()):
            # have to set rotation angle after period list has been set
            if key != 'rotation_angle':
                if hasattr(self, key):
                    setattr(self, key, kwargs[key])
                else:
                    self._logger.warn("Argument {}={} is not supported thus not been set.".format(key, kwargs[key]))

        # update period buffer to a default value if it is invalid
        if self.period_buffer is not None:
            try:
                self.period_buffer = float(self.period_buffer)
                if self.period_buffer < 0.:
                    self._logger.warn("Period buffer must be > 0, setting to None")
                    self.period_buffer = None
                # if provided value between 0 and 1, assume it was meant to be set to 1 + provided value
                elif self.period_buffer < 1.:
                    self._logger.warn("Period buffer must be > 1, adding 1 to provided value")
                    self.period_buffer += 1.
            except:
                self._logger.warn("Period buffer must be convertable to an integer or float, setting to None")
                
                

        if 'rotation_angle' in list(kwargs.keys()):
            setattr(self, 'rotation_angle', kwargs['rotation_angle'])
#            self._set_rotation_angle(self.rotation_angle)


    def _set_dtype(self, z_shape, t_shape):
        """
        reset dtype
        """

        self._z_shape = z_shape
        self._t_shape = t_shape

        self._dtype = [('station', '|U10'),
                       ('lat', np.float),
                       ('lon', np.float),
                       ('elev', np.float),
                       ('rel_east', np.float),
                       ('rel_north', np.float),
                       ('rel_elev', np.float),
                       ('east', np.float),
                       ('north', np.float),
                       ('zone', '|S4'),
                       ('z', (np.complex, self._z_shape)),
                       ('z_err', (np.float, self._z_shape)),
                       ('z_inv_err', (np.float, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.float, self._t_shape)),
                       ('tip_inv_err', (np.float, self._t_shape))]

    @staticmethod
    def get_header_string(error_type, error_value, rotation_angle):
        """
        reset the header sring for file
        """
        
        h_str = []
        if (np.atleast_1d(error_type).ndim == 2):
            h_str = '# Created using MTpy calculated {},{},{},{} '.format(error_type[0,0], error_type[0,1], error_type[1,0], error_type[1,1])
        else:
            h_str = '# Created using MTpy calculated {} '.format(error_type)
                          
        
        if (np.atleast_1d(error_value).ndim == 2):
            h_str += 'error floors of {0:.0f}%,{1:.0f}%,{2:.0f}%,{3:.0f}%, data rotated {4:.1f}_deg clockwise from N\n'
            return h_str.format(error_value[0,0], error_value[0,1], error_value[1,0], error_value[1,1], rotation_angle)
        else:
            h_str += 'error of {1:.0f}% data rotated {2:.1f}_deg clockwise from N\n'

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
                             '.edi files containing the full path')

        self.mt_dict = {}
        for edi in self.edi_list:
            mt_obj = mt.MT(edi)
            self.mt_dict[mt_obj.station] = mt_obj

    def get_relative_station_locations(self):
        """
        get station locations from edi files
        """
        stations_obj = Stations(model_epsg=self.model_epsg,
                                model_utm_zone=self.model_utm_zone)
        mt_list = [self.mt_dict[s_key] for s_key in sorted(self.mt_dict.keys())]
        stations_obj.get_station_locations(mt_list)

        # rotate locations if needed
        if self._rotation_angle != 0:
            # rotate stations the opposite way to the data
            stations_obj.rotate_stations(-self._rotation_angle)

        # fill data array
        self.data_array[:]['station'] = stations_obj.station
        self.data_array[:]['lat'] = stations_obj.lat
        self.data_array[:]['lon'] = stations_obj.lon
        self.data_array[:]['east'] = stations_obj.east
        self.data_array[:]['north'] = stations_obj.north
        self.data_array[:]['elev'] = stations_obj.elev
        self.data_array[:]['rel_east'] = stations_obj.rel_east
        self.data_array[:]['rel_north'] = stations_obj.rel_north
        self.data_array[:]['rel_elev'] = stations_obj.rel_elev
        self.data_array[:]['zone'] = stations_obj.utm_zone
        
        # get center point
        self.center_point = stations_obj.center_point

    def get_period_list(self):
        """
        make a period list to invert for
        """
        if self.mt_dict is None:
            self.get_mt_dict()

        if self.period_list is not None:
            print('-' * 50)
            print('Inverting for periods:')
            for per in self.period_list:
                print('     {0:<12.6f}'.format(per))
            print('-' * 50)
            return

        data_period_list = []
        for s_key in sorted(self.mt_dict.keys()):
            mt_obj = self.mt_dict[s_key]
            data_period_list.extend(list(1. / mt_obj.Z.freq))

        self.data_period_list = np.array(sorted(list(set(data_period_list)),
                                                reverse=False))

        if self.period_min is not None and self.period_max is None:
            raise DataError('Need to input period_max')
        if self.period_max is not None and self.period_min is None:
            raise DataError('Need to input period_min')
        if self.period_min is not None and self.period_max is not None and self.max_num_periods is None:
            raise DataError('Need to input number of periods to use')
        
        min_index = np.where(self.data_period_list >= self.period_min)[0][0]
        max_index = np.where(self.data_period_list <= self.period_max)[0][-1]

        pmin = np.log10(self.data_period_list[min_index])
        pmax = np.log10(self.data_period_list[max_index])
        self.period_list = np.logspace(pmin, pmax, num=self.max_num_periods)

        print('-' * 50)
        print('Inverting for periods:')
        for per in self.period_list:
            print('     {0:<12.6f}'.format(per))
        print('-' * 50)
        
        if self.period_list is None:  # YG: is this possible?
            raise ModEMError('Need to input period_min, period_max, '
                             'max_num_periods or a period_list')

    def _set_rotation_angle(self, rotation_angle):
        """
        on set rotation angle rotate mt_dict and data_array,
        """
        if self._rotation_angle == rotation_angle:
            return

#        self._logger.info('Changing rotation angle from {0:.1f} to {1:.1f}'.format(
#            self._rotation_angle, rotation_angle))
#
#        self._rotation_angle = -self._rotation_angle + rotation_angle
#
#        if self.rotation_angle == 0:
#            return

        self._logger.info('Changing rotation angle from {0:.1f} to {1:.1f}'.format(
            self._rotation_angle, rotation_angle))
        
        self._rotation_angle = rotation_angle
        

        if self.mt_dict is None:
            self.get_mt_dict()

        for mt_key in sorted(self.mt_dict.keys()):
            mt_obj = self.mt_dict[mt_key]
            # check if data already rotated
            angle_to_rotate = self._rotation_angle - mt_obj.Z.rotation_angle
            mt_obj.Z.rotate(angle_to_rotate)
            mt_obj.Tipper.rotate(angle_to_rotate)



        self._logger.info('Data rotated to align with {0:.1f} deg clockwise from N'.format(
            self._rotation_angle))
        self.fill_data_array()

    def _get_rotation_angle(self):
        return self._rotation_angle

    rotation_angle = property(fget=_get_rotation_angle,
                              fset=_set_rotation_angle,
                              doc="""Rotate data assuming N=0, E=90""")

    def _initialise_empty_data_array(self, station_locations, period_list,
                                     location_type='LL', station_names=None,
                                     epsg = None, utm_zone=None):
        """
        create an empty data array to create input files for forward modelling
        station locations is an array containing x,y coordinates of each station
        (shape = (number_of_stations,2))
        period_list = list of periods to model
        location_type = 'LL' or 'EN' - longitude/latitude or easting/northing
                        if 'EN' then utm_zone
        station_names = list or 1d array containing station names

        """
        self.period_list = period_list.copy()
        nf = len(self.period_list)
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(len(station_locations), dtype=self._dtype)
        if location_type == 'LL':
            self.data_array['lon'] = station_locations[:, 0]
            self.data_array['lat'] = station_locations[:, 1]
            for i in range(len(self.data_array['lon'])):
                lat,lon = self.data_array['lat'][i],self.data_array['lon'][i]
                east, north, zone = gis_tools.project_point_ll2utm(lat,lon,epsg=epsg,utm_zone=utm_zone)
                self.data_array['east'][i] = east
                self.data_array['north'][i] = north              
        else:
            self.data_array['east'] = station_locations[:, 0]
            self.data_array['north'] = station_locations[:, 1]
            for i in range(len(self.data_array['east'])):
                east, north = self.data_array['east'][i],self.data_array['north'][i]
                lat,lon = gis_tools.project_point_utm2ll(east,north,
                                                             utm_zone=utm_zone, 
                                                             epsg=epsg)
                self.data_array['lon'][i] = lon
                self.data_array['lat'][i] = lat

        # set non-zero values to array (as zeros will be deleted)
        # as we are setting up for forward modelling, actual values don't matter
        if self.inv_mode in '12':
            self.data_array['z'][:] = 10. + 10j
            self.data_array['z_err'][:] = 1e15
        if self.inv_mode == '1':
            self.data_array['tip'][:] = 0.1 + 0.1j
            self.data_array['tip_err'][:] = 1e15

        # set station names
        if station_names is not None:
            if len(station_names) != len(station_names):
                station_names = None

        if station_names is None:
            station_names = ['st%03i' %
                             ss for ss in range(len(station_locations))]
        self.data_array['station'] = station_names
        
        # make an mt_dict
        self.mt_dict = {}
        for i,sname in enumerate(station_names):
            mtObj = mt.MT()
            mtObj.lat = self.data_array['lat'][i]
            mtObj.lon = self.data_array['lon'][i]

            mtObj.east = self.data_array['east'][i]
            mtObj.north = self.data_array['north'][i]
            mtObj.Z = mtz.Z(z_array=self.data_array['z'][i],
                            z_err_array=self.data_array['z_err'][i],
                            freq=1./period_list)
            mtObj.Tipper = mtz.Tipper(tipper_array=self.data_array['tip'][i],
                            tipper_err_array=self.data_array['tip_err'][i],
                            freq=1./period_list)
            mtObj.station = sname
            self.mt_dict[sname] = mtObj

        self.get_relative_station_locations()

    def fill_data_array(self, new_edi_dir=None, use_original_freq=False, longitude_format='LON'):
        """
        fill the data array from mt_dict

        """

        if self.period_list is None:
            self.get_period_list()

        ns = len(list(self.mt_dict.keys()))
        nf = len(self.period_list)

        d_array = False
        if self.data_array is not None:
            d_arr_copy = self.data_array.copy()
            d_array = True

        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)

        rel_distance = False
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]
            if d_array:
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
                    self.data_array[ii]['rel_elev'] = d_arr_copy[d_index]['rel_elev']
                    self.data_array[:]['zone'] = d_arr_copy[d_index]['zone']
                except IndexError:
                    self._logger.warn('Could not find {0} in data_array'.format(s_key))
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
                    self.data_array[ii]['rel_elev'] = mt_obj.grid_elev
                    rel_distance = True
                except AttributeError:
                    self._logger.warn("Unable to set relative locations from 'mt_obj' "
                                      "- not yet implemented")
                    pass

            # interpolate each station onto the period list
            # check bounds of period list
            interp_periods = self.period_list[np.where(
                (self.period_list >= 1. / mt_obj.Z.freq.max()) &
                (self.period_list <= 1. / mt_obj.Z.freq.min()))]

            # if specified, apply a buffer so that interpolation doesn't
            # stretch too far over periods
            if type(self.period_buffer) in [float, int]:
                interp_periods_new = []
                dperiods = 1. / mt_obj.Z.freq
                for iperiod in interp_periods:
                    # find nearest data period
                    difference = np.abs(iperiod - dperiods)
                    nearestdperiod = dperiods[difference == np.amin(difference)][0]
                    if max(nearestdperiod / iperiod, iperiod / nearestdperiod) < self.period_buffer:
                        interp_periods_new.append(iperiod)

                interp_periods = np.array(interp_periods_new)

            # FZ: sort in order
            interp_periods = np.sort(interp_periods)
            self._logger.debug("station_name and its original period: %s %s %s",
                               mt_obj.station, len(mt_obj.Z.freq), 1.0 / mt_obj.Z.freq)
            self._logger.debug("station_name and interpolation period: %s %s %s",
                               mt_obj.station, len(interp_periods), interp_periods)

            # default: use_original_freq = True, each MT station edi file will use it's own frequency-filtered.
            # no new freq in the output modem.dat file. select those freq of mt_obj according to interp_periods
            if use_original_freq:
                interp_periods = self.filter_periods(mt_obj, interp_periods)
                self._logger.debug("station_name and selected/filtered periods: %s, %s, %s", mt_obj.station,
                                   len(interp_periods), interp_periods)
                # in this case the below interpolate_impedance_tensor function will degenerate into a same-freq set.
                
            if len(interp_periods) > 0:  # not empty
                interp_z, interp_t = mt_obj.interpolate(1. / interp_periods, period_buffer=self.period_buffer ,bounds_error=False)  #)
                # set rotation angle
                interp_z.rotation_angle = self.rotation_angle*np.ones(len(interp_z.z))
                interp_t.rotation_angle = self.rotation_angle*np.ones(len(interp_t.tipper))
                #                interp_z, interp_t = mt_obj.interpolate(1./interp_periods)
                for kk, ff in enumerate(interp_periods):
                    jj = np.where(self.period_list == ff)[0][0]
                    self.data_array[ii]['z'][jj] = interp_z.z[kk, :, :]
                    self.data_array[ii]['z_err'][jj] = interp_z.z_err[kk, :, :]

                    if mt_obj.Tipper.tipper is not None:
                        self.data_array[ii]['tip'][jj] = interp_t.tipper[kk, :, :]
                        self.data_array[ii]['tip_err'][jj] = \
                            interp_t.tipper_err[kk, :, :]

                # FZ: try to output a new edi files. Compare with original edi?
                if new_edi_dir is not None and os.path.isdir(new_edi_dir):
                    # new_edifile = os.path.join(new_edi_dir, mt_obj.station + '.edi')
                    mt_obj.write_mt_file(
                        save_dir=new_edi_dir,
                        fn_basename=mt_obj.station,
                        file_type='edi',
                        new_Z_obj=interp_z,
                        new_Tipper_obj=interp_t,
                        longitude_format=longitude_format)
            else:
                pass

        # BM: If we can't get relative locations from MT object, 
        #  then get them from Station object
        if not rel_distance:
            self.get_relative_station_locations()

        return

    @staticmethod
    def filter_periods(mt_obj, per_array):
        """Select the periods of the mt_obj that are in per_array.
        used to do original freq inversion.

        :param mt_obj:
        :param per_array:
        :return: array of selected periods (subset) of the mt_obj
        """

        mt_per = 1.0 / mt_obj.Z.freq

        new_per = [p for p in mt_per if any([np.isclose(p, p2, 1.e-8) 
                   for p2 in per_array])]
        # for p in mt_per:
        #     for p2 in per_array:
        #         # if abs(p - p2) < 0.00000001:  # Be aware of floating error if use ==
        #         if np.isclose(p, p2, 1.e-8):
        #             new_per.append(p)

        return np.array(new_per)

    def _set_station_locations(self, station_locations):
        """
        take a station_locations array and populate data_array
        """
        if self.data_array is None:
            self._set_dtype((len(self.period_list), 2, 2),
                            (len(self.period_list), 1, 2))
            self.data_array = np.zeros(station_locations.station_locations.size,
                                       dtype=self._dtype)
            for d_index, s_arr in enumerate(station_locations.station_locations):
                self.data_array[d_index]['lat'] = s_arr['lat']
                self.data_array[d_index]['lon'] = s_arr['lon']
                self.data_array[d_index]['east'] = s_arr['east']
                self.data_array[d_index]['north'] = s_arr['north']
                self.data_array[d_index]['elev'] = s_arr['elev']
                self.data_array[d_index]['rel_east'] = s_arr['rel_east']
                self.data_array[d_index]['rel_north'] = s_arr['rel_north']
                self.data_array[d_index]['rel_elev'] = s_arr['rel_elev']

        else:
            for s_arr in station_locations.station_locations:
                try:
                    d_index = np.where(self.data_array['station'] ==
                                       s_arr['station'])[0][0]
                except IndexError:
                    self._logger.warn('Could not find {0} in data_array'.format(s_arr['station']))
                    d_index = None

                if d_index is not None:
                    self.data_array[d_index]['lat'] = s_arr['lat']
                    self.data_array[d_index]['lon'] = s_arr['lon']
                    self.data_array[d_index]['east'] = s_arr['east']
                    self.data_array[d_index]['north'] = s_arr['north']
                    self.data_array[d_index]['elev'] = s_arr['elev']
                    self.data_array[d_index]['rel_east'] = s_arr['rel_east']
                    self.data_array[d_index]['rel_north'] = s_arr['rel_north']
                    self.data_array[d_index]['rel_elev'] = s_arr['rel_elev']

    def _get_station_locations(self):
        """
        extract station locations from data array
        """
        if self.data_array is None:
            return None

        station_locations = self.data_array[['station', 'lat', 'lon',
                                             'north', 'east', 'elev',
                                             'rel_north', 'rel_east', 
                                             'rel_elev', 'zone']]
        stations_obj = Stations(model_epsg=self.model_epsg,
                                model_utm_zone=self.model_utm_zone)
        stations_obj.station_locations = station_locations

        return stations_obj

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

        # compute relative error for tipper
        if 'floor' in self.error_type_tipper:
            t_index = np.where(self.data_array['tip_err'] < self.error_value_tipper)
            self.data_array['tip_inv_err'][t_index] = self.error_value_tipper
        elif 'abs' in self.error_type_tipper:
            self.data_array['tip_inv_err'][:] = self.error_value_tipper
        else:
            raise DataError("Unsupported error type (tipper): {}".format(self.error_type_tipper))

        # consistency checks
        error_type_z_list = np.atleast_1d(self.error_type_z)

        if (error_type_z_list.size != 1 and error_type_z_list.size != 4):
            raise DataError('Either specify a single error_type_z for all components, or ' \
                            'a 2x2 numpy array of error_type_z.')
        # end if

        # compute error for z
        err_value = self.error_value_z / 100.
        for ss in range(self.data_array.shape[0]):
            for ff in range(max([self._t_shape[0], self._z_shape[0]])):
                d_xx = abs(self.data_array['z'][ss, ff, 0, 0])
                d_xy = abs(self.data_array['z'][ss, ff, 0, 1])
                d_yx = abs(self.data_array['z'][ss, ff, 1, 0])
                d_yy = abs(self.data_array['z'][ss, ff, 1, 1])
                d = np.array([d_xx, d_xy, d_yx, d_yy])
                nz = np.nonzero(d)

                if d.sum() == 0.0:  # YG: only works if all d >= 0
                    continue

                err = np.zeros((error_type_z_list.size,
                                np.atleast_2d(err_value).shape[0],
                                np.atleast_2d(err_value).shape[1]))
                for ei, error_type_z in enumerate(error_type_z_list.flatten()):
                    if 'egbert' in error_type_z:
                        # if both components masked, then take error floor from
                        # max of z_xx or z_yy
                        if (d_xy==0.0 and d_yx==0.0):
                            err[ei] = err_value * np.max([d_xx,d_yy])
                        # else use the off diagonals depending on data availability
                        else:
                            if d_xy == 0.0:
                                d_xy = d_yx
                            if d_yx == 0.0:
                                d_yx = d_xy
                            err[ei] = err_value * np.sqrt(d_xy * d_yx)

                    elif 'median' in error_type_z:
                        err[ei] = err_value * np.median(d[nz])

                    elif 'mean_od' in error_type_z:
                        dod = np.array([d_xy, d_yx])
                        nzod = np.nonzero(dod)
                        err[ei] = err_value * np.mean(dod[nzod])

                    elif 'eigen' in error_type_z:
                        d2d = d.reshape((2, 2))
                        err[ei] = err_value * np.abs(np.linalg.eigvals(d2d)).mean()
                        if np.atleast_1d(err[ei]).sum() == 0:
                            err[ei] = err_value * d[nz].mean()

                    elif 'off_diagonals' in error_type_z:
                        # apply same error to xy and xx, and to yx and yy
                        # value is a % of xy and yx respectively
                        err[ei] = np.array([[d_xy, d_xy],[d_yx, d_yx]])*err_value

                    elif 'percent' in error_type_z:
                        # apply separate error floors to each component
                        err[ei] = err_value * np.abs(d[ei])
                    else:
                        raise DataError('error type (z) {0} not understood'.format(error_type_z))
                # end for

                if(error_type_z_list.size == 1):
                    self.data_array['z_inv_err'][ss, ff, :, :] = err[0]
                else:
                    for ei in np.arange(error_type_z_list.size):
                        ix, iy = np.divmod(ei, 2)
                        if (err.shape[1] > 1):
                            self.data_array['z_inv_err'][ss, ff, ix, iy] = err[ei, ix, iy]
                        else:
                            self.data_array['z_inv_err'][ss, ff, ix, iy] = err[ei, 0, 0]
                        # end if
                    # end for
                # end if

        # if there is an error floor
        if 'floor' in self.error_type_z:
            f_index = np.where(self.data_array['z_inv_err'] < self.data_array['z_err'])
            self.data_array['z_inv_err'][f_index] = self.data_array['z_err'][f_index]



    def write_data_file(self, save_path=None, fn_basename=None,
                        rotation_angle=None, compute_error=True, fill=True,
                        elevation=False, use_original_freq=False, longitude_format='LON'):
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
            self.data_fn = fn_basename

        print ("self.save_path, self.data_fn ==",self.save_path, self.data_fn)
        self.data_fn = os.path.join(self.save_path, self.data_fn)

        self.get_period_list()

        # rotate data if desired
        if rotation_angle is not None:
            self.rotation_angle = rotation_angle

        # be sure to fill in data array
        if fill:
            new_edi_dir = os.path.join(self.save_path, 'new_edis')  # output edi files according to selected periods
            if not os.path.exists(new_edi_dir):
                os.mkdir(new_edi_dir)
            self.fill_data_array(new_edi_dir=new_edi_dir,
                                 use_original_freq=use_original_freq,
                                 longitude_format=longitude_format)

        if not elevation:
            self.data_array['rel_elev'][:] = 0.0
            self.center_point.elev = 0.0

        d_lines = []
        for inv_mode in self.inv_mode_dict[self.inv_mode]:
            if 'impedance' in inv_mode.lower():
                d_lines.append(self.get_header_string(self.error_type_z,
                                                      self.error_value_z,
                                                      self.rotation_angle))
                d_lines.append(self.header_string)
                d_lines.append('> {0}\n'.format(inv_mode))
                d_lines.append('> exp({0}i\omega t)\n'.format(self.wave_sign_impedance))
                d_lines.append('> {0}\n'.format(self.units))

                n_sta = len(np.nonzero(np.abs(self.data_array['z']).sum(axis=(1, 2, 3)))[0])
                n_per = len(np.nonzero(np.abs(self.data_array['z']).sum(axis=(0, 2, 3)))[0])
            elif 'vertical' in inv_mode.lower():
                d_lines.append(self.get_header_string(self.error_type_tipper,
                                                      self.error_value_tipper,
                                                      self.rotation_angle))
                d_lines.append(self.header_string)
                d_lines.append('> {0}\n'.format(inv_mode))
                d_lines.append('> exp({0}i\omega t)\n'.format(self.wave_sign_tipper))
                d_lines.append('> []\n')
                n_sta = len(np.nonzero(np.abs(self.data_array['tip']).sum(axis=(1, 2, 3)))[0])
                n_per = len(np.nonzero(np.abs(self.data_array['tip']).sum(axis=(0, 2, 3)))[0])
            else:
                # maybe error here
                raise NotImplementedError("inv_mode {} is not supported yet".format(inv_mode))

            d_lines.append('> 0\n')  # orientation, need to add at some point
            if elevation:
                d_lines.append('> {0:>10.6f} {1:>10.6f} {2:>10.2f}\n'.format(
                               self.center_point.lat[0],
                               self.center_point.lon[0],
                               self.center_point.elev[0]))
            else:
                d_lines.append('> {0:>10.6f} {1:>10.6f}\n'.format(
                               self.center_point.lat[0],
                               self.center_point.lon[0]))
            d_lines.append('> {0} {1}\n'.format(n_per, n_sta))

            if compute_error:
                self.compute_inv_error()

            for ss in range(self.data_array['z'].shape[0]):
                for ff in range(self.data_array['z'].shape[1]):
                    for comp in self.inv_comp_dict[inv_mode]:
                        # index values for component with in the matrix
                        z_ii, z_jj = self.comp_index_dict[comp]

                        # get the correct key for data array according to comp
                        if comp.find('z') == 0:
                            c_key = 'z'
                        elif comp.find('t') == 0:
                            c_key = 'tip'

                        # get the value for that competent at that frequency
                        zz = self.data_array[ss][c_key][ff, z_ii, z_jj]
                        if zz.real != 0.0 and zz.imag != 0.0 and zz.real != 1e32 and zz.imag != 1e32:
                            if self.formatting == '1':
                                per = '{0:<12.5e}'.format(self.period_list[ff])
                                sta = '{0:>7}'.format(self.data_array[ss]['station'])#.decode('UTF-8'))
                                lat = '{0:> 9.3f}'.format(self.data_array[ss]['lat'])
                                lon = '{0:> 9.3f}'.format(self.data_array[ss]['lon'])
                                eas = '{0:> 12.3f}'.format(self.data_array[ss]['rel_east'])
                                nor = '{0:> 12.3f}'.format(self.data_array[ss]['rel_north'])
                                ele = '{0:> 12.3f}'.format(self.data_array[ss]['rel_elev'])
                                com = '{0:>4}'.format(comp.upper())
                                if self.units.lower() == 'ohm':
                                    rea = '{0:> 14.6e}'.format(zz.real / 796.)
                                    ima = '{0:> 14.6e}'.format(zz.imag / 796.)
                                elif self.units.lower() not in ("[v/m]/[t]", "[mv/km]/[nt]"):
                                    raise DataError("Unsupported unit \"{}\"".format(self.units))
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
                                ele = '{0:> 10.3f}'.format(self.data_array[ss]['rel_elev'])
                                com = '{0:>12}'.format(comp.upper())
                                if self.units.lower() == 'ohm':
                                    rea = '{0:> 17.6e}'.format(zz.real / 796.)
                                    ima = '{0:> 17.6e}'.format(zz.imag / 796.)
                                elif self.units.lower() not in ("[v/m]/[t]", "[mv/km]/[nt]"):
                                    raise DataError("Unsupported unit \"{}\"".format(self.units))
                                else:
                                    rea = '{0:> 17.6e}'.format(zz.real)
                                    ima = '{0:> 17.6e}'.format(zz.imag)
                            else:
                                raise NotImplementedError(
                                    "format {}({}) is not supported".format(self.formatting, type(self.formatting)))

                            # get error from inversion error
                            abs_err = self.data_array['{0}_inv_err'.format(c_key)][ss, ff, z_ii, z_jj]

                            if np.isinf(abs_err) or np.isnan(abs_err):
                                abs_err = 10 ** (np.floor(np.log10(abs(max([float(rea),
                                                                            float(ima)])))))
                            abs_err = '{0:> 14.6e}'.format(abs(abs_err))
                            # make sure that x==north, y==east, z==+down
                            dline = ''.join([per, sta, lat, lon, nor, eas, ele, com, rea, ima, abs_err, '\n'])

                            d_lines.append(dline)
        print("self.data_fn ==",  self.data_fn)
        with open(self.data_fn, 'w') as dfid:
            dfid.writelines(d_lines)

        self._logger.info('Wrote ModEM data file to {0}'.format(self.data_fn))
        return self.data_fn

    @deprecated("error type from GA implementation, not fully tested yet")
    def _impedance_components_error_meansqr(self, c_key, ss, z_ii, z_jj):
        """
        calculate the mean square of errors of a given component over all frequencies for a given station
        :param c_key:
        :param ss:
        :param z_ii:
        :param z_jj:
        :return:
        """
        abs_err = np.mean(np.square(self.data_array[ss][c_key + '_err'][:, z_ii, z_jj]))
        return abs_err

    @deprecated("error type from GA implementation, not fully tested yet")
    def _impedance_components_error_sqr(self, c_key, ff, ss, z_ii, z_jj):
        """
        use the square of the error of a given frequency and a given component at the given station
        :param c_key:
        :param ff:
        :param ss:
        :param z_ii:
        :param z_jj:
        :return:
        """
        return np.square(self.data_array[ss][c_key + '_err'][ff, z_ii, z_jj])

    @deprecated("error type from GA implementation, not fully tested yet")
    def _impedance_components_error_stddev(self, c_key, ss, z_ii, z_jj):
        """
        calculate the stddev across all frequencies on a given component
        :param c_key:
        :param ss:
        :param z_ii:
        :param z_jj:
        :return:
        """
        # errors = [self.data_array[ss][c_key + '_err'][freq, z_ii, z_jj] for freq in range(self.data_array['z'].shape[1])]
        # print errors
        # abs_err = np.std(errors)
        # print abs_err
        errors = self.data_array[ss][c_key + '_err'][:, z_ii, z_jj]
        # print errors
        abs_err = np.std(errors)
        # print abs_err
        return abs_err

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

        if not os.path.isfile(ws_data_fn):
            raise ws.WSInputError('Did not find {0}, check path'.format(ws_data_fn))

        if save_path is not None:
            self.save_path = save_path
        else:
            self.save_path = os.path.dirname(ws_data_fn)

        if fn_basename is not None:
            self.fn_basename = fn_basename

        # --> get data from data file
        wsd = ws.WSData()
        wsd.read_data_file(ws_data_fn, station_fn=station_fn)

        ns = wsd.data['station'].shape[0]
        nf = wsd.period_list.shape[0]

        self.period_list = wsd.period_list.copy()
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)

        # --> fill data array
        for ii, d_arr in enumerate(wsd.data):
            self.data_array[ii]['station'] = d_arr['station']
            self.data_array[ii]['rel_east'] = d_arr['east']
            self.data_array[ii]['rel_north'] = d_arr['north']
            self.data_array[ii]['z'][:] = d_arr['z_data']
            self.data_array[ii]['z_err'][:] = d_arr['z_data_err'].real * \
                                              d_arr['z_err_map'].real
            self.data_array[ii]['station'] = d_arr['station']
            self.data_array[ii]['lat'] = 0.0
            self.data_array[ii]['lon'] = 0.0
            self.data_array[ii]['rel_east'] = d_arr['east']
            self.data_array[ii]['rel_north'] = d_arr['north']
            self.data_array[ii]['elev'] = 0.0

        # need to change the inversion mode to be the same as the ws_data file
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

        # -->write file
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

        if ws_data_fn is None:
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
        data_dtype = [('station', '|U10'),
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
        ws_data.data['z_data_err'][:, :, :] = self.data_array['z_err'] * (1 + 1j)
        ws_data.data['z_err_map'][:, :, :] = np.array([[1, 1], [1, 1]])

        ws_data.write_data_file(save_path=save_path, data_fn=ws_data_fn)

        return ws_data.data_fn, station_info.station_fn

    def read_data_file(self, data_fn=None, center_utm=None):
        """ Read ModEM data file

       inputs:
        data_fn = full path to data file name
        center_utm = option to provide real world coordinates of the center of
                     the grid for putting the data and model back into
                     utm/grid coordinates, format [east_0, north_0, z_0]


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
            raise DataError('data_fn is None, enter a data file to read.')
        elif os.path.isfile(self.data_fn) is False:
            raise DataError(
                'Could not find {0}, check path'.format(self.data_fn))

        with open (self.data_fn, 'r') as dfid:
            dlines = dfid.readlines()

        # dfid.close()

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
                # modem outputs only 7 characters for the lat and lon
                # if there is a negative they merge together, need to split 
                # them up
                dline = dline.replace('-', ' -')
                metadata_list.append(dline[1:].strip())
                if dline.lower().find('ohm') > 0:
                    self.units = 'ohm'
                elif dline.lower().find('mv') > 0:
                    self.units = '[mV/km]/[nT]'
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
                        self.wave_sign_impedance = dline[dline.find('(') + 1]
                    elif read_tipper is True:
                        self.wave_sign_tipper = dline[dline.find('(') + 1]
                elif len(dline[1:].strip().split()) >= 2:
                    if dline.find('.') > 0:
                        value_list = [float(value) for value in
                                      dline[1:].strip().split()]

                        self.center_point = np.recarray(1, dtype=[('station', '|U10'),
                                                                  ('lat', np.float),
                                                                  ('lon', np.float),
                                                                  ('elev', np.float),
                                                                  ('rel_elev', np.float),
                                                                  ('rel_east', np.float),
                                                                  ('rel_north', np.float),
                                                                  ('east', np.float),
                                                                  ('north', np.float),
                                                                  ('zone', 'U4')])
                        self.center_point.lat = value_list[0]
                        self.center_point.lon = value_list[1]
                        try:
                            self.center_point.elev = value_list[2]
                        except IndexError:
                            self.center_point.elev = 0.0
                            print('Did not find center elevation in data file')

                        ce, cn, cz = gis_tools.project_point_ll2utm(self.center_point.lat,
                                                                    self.center_point.lon,
                                                                    epsg=self.model_epsg,
                                                                    utm_zone=self.model_utm_zone)

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

        # try to find rotation angle
        h_list = header_list[0].split()
        for hh, h_str in enumerate(h_list):
            if h_str.find('_deg') > 0:
                try:
                    self._rotation_angle = float(h_str[0:h_str.find('_deg')])
                    self._logger.info('Set rotation angle to {0:.1f} '.format(
                        self._rotation_angle) + 'deg clockwise from N')
                except ValueError:
                    pass

        # find inversion mode
        for inv_key in list(self.inv_mode_dict.keys()):
            inv_mode_list = self.inv_mode_dict[inv_key]
            if len(inv_mode_list) != inv_list:
                continue
            else:
                tf_arr = np.zeros(len(inv_list), dtype=np.bool)

                for tf, data_inv in enumerate(inv_list):
                    if data_inv in self.inv_mode_dict[inv_key]:
                        tf_arr[tf] = True

                if np.alltrue(tf_arr):
                    self.inv_mode = inv_key
                    break

        self.period_list = np.array(sorted(set(period_list)))
        station_list = sorted(set(station_list))

        # make a period dictionary to with key as period and value as index
        period_dict = dict([(per, ii) for ii, per in enumerate(self.period_list)])

        # --> need to sort the data into a useful fashion such that each station
        #    is an mt object

        data_dict = {}
        z_dummy = np.zeros((len(self.period_list), 2, 2), dtype='complex')
        t_dummy = np.zeros((len(self.period_list), 1, 2), dtype='complex')

        index_dict = {'zxx': (0, 0), 'zxy': (0, 1), 'zyx': (1, 0), 'zyy': (1, 1),
                      'tx': (0, 0), 'ty': (0, 1)}

        # dictionary for true false if station data (lat, lon, elev, etc)
        # has been filled already so we don't rewrite it each time
        tf_dict = {}
        for station in station_list:
            data_dict[station] = mt.MT()
            data_dict[station].Z = mtz.Z(z_array=z_dummy.copy(),
                                         z_err_array=z_dummy.copy().real,
                                         freq=1. / self.period_list)
            data_dict[station].Tipper = mtz.Tipper(tipper_array=t_dummy.copy(),
                                                   tipper_err_array=t_dummy.copy().real,
                                                   freq=1. / self.period_list)
            # make sure that the station data starts out with false to fill
            # the data later
            tf_dict[station] = False

        # fill in the data for each station
        for dd in data_list:
            # get the period index from the data line
            p_index = period_dict[dd[0]]
            # get the component index from the data line
            ii, jj = index_dict[dd[7].lower()]

            # if the station data has not been filled yet, fill it
            if not tf_dict[dd[1]]:
                data_dict[dd[1]].lat = dd[2]
                data_dict[dd[1]].lon = dd[3]
                data_dict[dd[1]].grid_north = dd[4]
                data_dict[dd[1]].grid_east = dd[5]
                data_dict[dd[1]].grid_elev = dd[6]
                data_dict[dd[1]].elev = dd[6]
                data_dict[dd[1]].station = dd[1]
                tf_dict[dd[1]] = True
            # fill in the impedance tensor with appropriate values
            if dd[7].find('Z') == 0:
                z_err = dd[10]
                if self.wave_sign_impedance == '+':
                    z_value = dd[8] + 1j * dd[9]
                elif self.wave_sign_impedance == '-':
                    z_value = dd[8] - 1j * dd[9]
                else:
                    raise DataError("Incorrect wave sign \"{}\" (impedance)".format(self.wave_sign_impedance))

                if self.units.lower() == 'ohm':
                    z_value *= 796.
                    z_err *= 796.
                elif self.units.lower() not in ("[v/m]/[t]", "[mv/km]/[nt]"):
                    raise DataError("Unsupported unit \"{}\"".format(self.units))

                data_dict[dd[1]].Z.z[p_index, ii, jj] = z_value
                data_dict[dd[1]].Z.z_err[p_index, ii, jj] = z_err
            # fill in tipper with appropriate values
            elif dd[7].find('T') == 0:
                if self.wave_sign_tipper == '+':
                    data_dict[dd[1]].Tipper.tipper[p_index, ii, jj] = dd[8] + 1j * dd[9]
                elif self.wave_sign_tipper == '-':
                    data_dict[dd[1]].Tipper.tipper[p_index, ii, jj] = dd[8] - 1j * dd[9]
                else:
                    raise DataError("Incorrect wave sign \"{}\" (tipper)".format(self.wave_sign_tipper))
                data_dict[dd[1]].Tipper.tipper_err[p_index, ii, jj] = dd[10]

        # make mt_dict an attribute for easier manipulation later
        self.mt_dict = data_dict

        ns = len(list(self.mt_dict.keys()))
        nf = len(self.period_list)
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)

        # Be sure to caclulate invariants and phase tensor for each station
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]

            # self.mt_dict[s_key].zinv.compute_invariants()
            self.mt_dict[s_key].pt.set_z_object(mt_obj.Z)
            self.mt_dict[s_key].Tipper.compute_amp_phase()
            self.mt_dict[s_key].Tipper.compute_mag_direction()
            
            self.data_array[ii]['station'] = mt_obj.station
            self.data_array[ii]['lat'] = mt_obj.lat
            self.data_array[ii]['lon'] = mt_obj.lon
            #east,north,zone = gis_tools.project_point_ll2utm(mt_obj.lat,mt_obj.lon,epsg=self.model_epsg)
            self.data_array[ii]['east'] = mt_obj.east
            self.data_array[ii]['north'] = mt_obj.north
            self.data_array[ii]['elev'] = mt_obj.elev
            self.data_array[ii]['rel_elev'] = mt_obj.grid_elev
            self.data_array[ii]['rel_east'] = mt_obj.grid_east
            self.data_array[ii]['rel_north'] = mt_obj.grid_north

            self.data_array[ii]['z'][:] = mt_obj.Z.z
            self.data_array[ii]['z_err'][:] = mt_obj.Z.z_err
            self.data_array[ii]['z_inv_err'][:] = mt_obj.Z.z_err

            self.data_array[ii]['tip'][:] = mt_obj.Tipper.tipper
            self.data_array[ii]['tip_err'][:] = mt_obj.Tipper.tipper_err
            self.data_array[ii]['tip_inv_err'][:] = mt_obj.Tipper.tipper_err

        
        # option to provide real world coordinates in eastings/northings
        # (ModEM data file contains real world center in lat/lon but projection
        # is not provided so utm is assumed, causing errors when points cross
        # utm zones. And lat/lon cut off to 3 d.p. causing errors in smaller
        # areas)
        if center_utm is not None:
            self.data_array['east'] = self.data_array[
                                          'rel_east'] + center_utm[0]
            self.data_array['north'] = self.data_array[
                                           'rel_north'] + center_utm[1]

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
                    self.station_locations.rel_north / 1000.,
                    self.station_locations.rel_east / 1000.,
                    self.station_locations.elev / 1000.,
                    data={'elevation': self.station_locations.elev / 1000.})

        self._logger.info('Wrote station file to {0}'.format(vtk_fn))

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
            e_index = np.where(m_obj.grid_east >= s_arr['rel_east'])[0][0] - 1
            n_index = np.where(m_obj.grid_north >= s_arr['rel_north'])[0][0] - 1

            mid_east = m_obj.grid_east[e_index:e_index + 2].mean()
            mid_north = m_obj.grid_north[n_index:n_index + 2].mean()

            s_index = np.where(self.data_array['station'] == s_arr['station'])[0][0]

            self.data_array[s_index]['rel_east'] = mid_east
            self.data_array[s_index]['rel_north'] = mid_north

    def change_data_elevation(self, model_obj, data_fn=None, res_air=1e12):
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

        s_locations = self.station_locations.station_locations.copy()

        # need to subtract one because we are finding the cell next to it
        for s_arr in s_locations:
            e_index = np.where(model_obj.grid_east >= s_arr['rel_east'])[0][0] - 1
            n_index = np.where(model_obj.grid_north >= s_arr['rel_north'])[0][0] - 1
            z_index = np.where(model_obj.res_model[n_index, e_index, :] < res_air * .9)[0][0]
            s_index = np.where(self.data_array['station'] == s_arr['station'])[0][0]
            
            # elevation needs to be relative to the point (0, 0, 0), where
            # 0 is the top of the model, the highest point, so the station 
            # elevation is relative to that.
            self.data_array[s_index]['rel_elev'] = model_obj.nodes_z[0:z_index].sum()
            
            # need to add in the max elevation to the center point so that
            # when we read the file later we can adjust the stations
            self.center_point.elev = model_obj.grid_z[0]
            
    def project_stations_on_topography(self, model_object, air_resistivity=1e12):
        """
        Re-write the data file to change the elevation column.
        And update covariance mask according topo elevation model.
        :param model_object:
        :param air_resistivity:
        :return:
        """

        # sx = self.station_locations.station_locations['rel_east']
        # sy = self.station_locations.station_locations['rel_north']

        # find index of each station on grid
        station_index_x = []
        station_index_y = []
        for sname in self.station_locations.station_locations['station']:
            ss = np.where(self.station_locations.station_locations['station'] == sname)[0][0]
            # relative locations of stations
            sx, sy = self.station_locations.station_locations['rel_east'][ss], \
                     self.station_locations.station_locations['rel_north'][ss]
            # indices of stations on model grid
            sxi = np.where((sx <= model_object.grid_east[1:]) & (
                sx > model_object.grid_east[:-1]))[0][0]
            syi = np.where((sy <= model_object.grid_north[1:]) & (
                sy > model_object.grid_north[:-1]))[0][0]

            # first, check if there are any air cells
            if np.any(model_object.res_model[syi, sxi] > 0.95 * air_resistivity):
                szi = np.amin(
                    np.where((model_object.res_model[syi, sxi] < 0.95 * air_resistivity))[0])
            # otherwise place station at the top of the model
            else:
                szi = 0
            
            # get relevant grid point elevation
            topoval = model_object.grid_z[szi]
            
            station_index_x.append(sxi)
            station_index_y.append(syi)

            # update elevation in station locations and data array, +1 m as
            # data elevation needs to be below the topography (as advised by Naser)
            # ====================== ====================================================
            # The following line have been commented as
            # self.db_array and elf.station_locations.station_locations['elev'][ss]
            # point to same location
            # ====================== ====================================================
            # self.station_locations.station_locations['elev'][ss] = topoval + 0.1
            self.data_array['rel_elev'][ss] = topoval + 0.001
            
            print('{0} at E={1}, N={2}, z={3}, model_z={4}'.format(sname,
                                                                   sxi,
                                                                   syi,
                                                                   topoval,
                                                                   self.data_array['rel_elev'][ss]))

        # BM: After applying topography, center point of grid becomes
        #  highest point of surface model.
        self.center_point.elev = model_object.grid_z[0]

        # logger.debug("Re-write data file after adding topo")
        self.write_data_file(fn_basename=os.path.basename(self.data_fn)[:-4]+'_topo.dat',
                             fill=False, elevation=True)  # (Xi, Yi, Zi) of each station-i may be shifted

        # debug self.Data.write_data_file(save_path='/e/tmp', fill=False)

        return station_index_x, station_index_y

    # FZ: moved from the modem_data_to_phase_tensor.py ref: AUSLAMP-112
    def compute_phase_tensor(self, datfile, outdir):
        """
        Compute the phase tensors from a ModEM dat file
        :param datfile: path2/file.dat
        :return: path2csv  created by this method
        """

        # Data file
        data_file = datfile
        dest_dir = outdir

        # Create a new ModEM data instance
        md = Data()
        # Read the datafile
        md.read_data_file(data_fn=data_file)

        num_sites = md.data_array.shape[0]
        print ("ModEM data file number of sites:", num_sites)

        first_site_periods = md.data_array[0][9]  # (23L, 2L, 2L)
        print ("first_site_periods = %s" % str(first_site_periods.shape[0]))

        period_list = md.period_list
        freq_list = 1.0 / period_list
        num_periods = len(period_list)
        print ("ModEM data file number of periods:", num_periods)

        csv_basename ="modem_data_to_phase_tensor"
        csvfname = os.path.join(dest_dir, "%s.csv" % csv_basename)

        csv_header = [
            'Freq', 'Station', 'Lat', 'Long', 'Phimin', 'Phimax', 'Ellipticity', 'Azimuth']

        with open(csvfname, "wb") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

        for period_num in range(num_periods):
            per= period_list[period_num]
            freq = freq_list[period_num]
            self._logger.info("Working on period %s; frequency: %s", per, freq )
            csvrows = []
            for num_site in range(num_sites):
                # Obtain the site for this number
                this_site = md.data_array[num_site]

                # Station Label is first record in data array for a site
                site_label = this_site[0]
                # Latitude is the second record in data array for a site
                site_lat = this_site[1]
                # Longitude is the third record in data array for a site
                site_long = this_site[2]

                # Actual data is a list of 2x2 matrices (ordered by period), 10th record in data array for a site
                this_site_data = this_site[9]

                # Now get the data for this period only, going off the period number we're looping through currently
                this_period_data = this_site_data[period_num]
                # Create a Z object based on this 2x2 array of data
                this_period_data = mtz.Z(z_array=this_period_data, freq=freq_list)

                # Given the Z object we just created, give us a PhaseTensor object
                this_phase_tensor = pt.PhaseTensor(z_object=this_period_data)

                # Get the four parameters we care about
                this_phimin = this_phase_tensor.phimin[0]
                this_phimax = this_phase_tensor.phimax[0]
                this_ellipticity = this_phase_tensor.ellipticity[0]
                this_azimuth = this_phase_tensor.azimuth[0]

                # Print out comma delimited version of the parameters: label, lat, long, phimin, phimax, ellipticity, azimuth
                arow = [freq, site_label, site_lat, site_long, this_phimin, this_phimax, this_ellipticity, this_azimuth]
                # Done for this site

                csvrows.append(arow)

            with open(csvfname, "ab") as csvf:  # append to this summary csv file for all freqs
                writer = csv.writer(csvf)
                writer.writerows(csvrows)

            csv_basename2 = "%s_%sHz.csv" % (csv_basename, str(freq))
            csvfile2 = os.path.join(dest_dir, csv_basename2)

            with open(csvfile2, "wb") as csvf:  # csvfile  for eachindividual freq
                writer = csv.writer(csvf)
                writer.writerow(csv_header)
                writer.writerows(csvrows)

        # Done with all sites and periods
        self._logger.info("CSV files created in %s", outdir)

        return csvfname
