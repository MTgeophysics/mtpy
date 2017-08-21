#!/usr/bin/env python
"""
Description:
   This module is refactored from modem.py which is too big to manage edit
    Define the Data class

Author: fei.zhang@ga.gov.au

Date: 2017-06-05
"""

__author__ = 'fei.zhang@ga.gov.au'

import glob
# from __future__ import print_function
import os
import os.path as op
import sys

import numpy as np

import mtpy.core.mt as mt
import mtpy.core.z as mtz
import mtpy.modeling.ws3dinv as ws
import mtpy.utils.latlon_utm_conversion as utm2ll
from mtpy import constants
from mtpy.core.edi_collection import EdiCollection
from mtpy.utils.mtpylog import MtPyLog

try:
    from evtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')

    print ('Note: if you are using Windows you should build evtk first with'
           'either MinGW or cygwin using the command: \n'
           '    python setup.py build -compiler=mingw32  or \n'
           '    python setup.py build -compiler=cygwin')

logger = MtPyLog().get_mtpy_logger(__name__)


class DataError(Exception):
    """Raise for ModEM Data class specific exception"""
    pass


# =============================================================================
class Data(object):
    """ Read and write .dat files for ModEM and convert a WS data file
    to ModEM format.
    Output effective new edi files (selected inversion periods) into the output dir

    ..note: :: the data is interpolated onto the given periods such that all
               stations invert for the same periods.  The interpolation is
               a linear interpolation of each of the real and imaginary parts
               of the impedance tensor and induction tensor.
               See mtpy.core.mt.MT.interpolate for more details

    Arguments
    ------------
        **edi_list** : list of full paths to .edi files you want to invert for

    ====================== ====================================================
    Attributes/Key Words   Description
    ====================== ====================================================
    _dtype                 internal variable defining the data type of
                           data_array
    _t_shape               internal variable defining shape of tipper array in
                           _dtype
    _z_shape               internal variable defining shape of Z array in
                           _dtype
    center_position_EN        (east, north, evel) for center point of station
                           array.  All stations are relative to this location
                           for plotting purposes.
    comp_index_dict        dictionary for index values of component of Z and T
    station_locations      numpy.ndarray structured to store station
                           location values.  Keys are:
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
    error_egbert           percentage to multiply sqrt(Z_xy*Zyx) by.
                           *default* is 3 as prescribed by Egbert & Kelbert
    error_floor            percentage to set the error floor at, anything below
                           this number will be set to error_floor.
                           *default* is 10
    error_tipper           absolute tipper error, all tipper error will be
                           set to this value unless you specify error_type as
                           'floor' or 'floor_egbert'.
                           *default* is .05 for 5%
    error_type             [ 'floor' | 'value' | 'egbert' | 'floor_egbert' |'stddev' | 'sqr' | 'meansqr' ]
                           *default* is 'egbert'
                                * 'floor' sets the error floor to error_floor
                                * 'value' sets error to error_value
                                * 'egbert' sets error to
                                           error_egbert * sqrt(abs(zxy*zyx))
                                * 'floor_egbert' sets error floor to
                                           error_egbert * sqrt(abs(zxy*zyx))
                                * 'stddev' use the stddev of the errors of a given component of a station
                                    across all frequencies
                                * 'sqr' use square error of a frequency of a component of a station
                                * 'meansqr' mean sqr of the errors for a given component of a station across
                                    all crequenction
    comp_error_type        dictionary contains some of the following keys
                            { 'zxx', 'zxy', 'zyx', 'zyy'}
                            where each key is associated with an error type from error_type
                            the error of each component will be calculated using its specified error
                            type. If the error type is not specified for a component, the value of error_type
                            will be used
                            *default* is None

    error_value            percentage to multiply Z by to set error
                           *default* is 5 for 5% of Z as error
    fn_basename            basename of data file. *default* is 'ModEM_Data.dat'
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
                               filling in station_locations
    read_data_file             read in a ModEM data file and fill attributes
                               data_array, station_locations, period_list, mt_dict
    write_data_file            write a ModEM data file
    ========================== ================================================


    :Example 1 --> create inversion period list: ::

        >>> import os
        >>> from mtpy.modeling.modem_data import Data
        >>> edi_path = r"/home/mt/edi_files"
        >>> edi_list = [os.path.join(edi_path, edi) \
                        for edi in os.listdir(edi_path)\
                        if edi.find('.edi') > 0]
        >>> md = Data(edi_list, period_min=.1, period_max=300,\
                            max_num_periods=12)
        >>> md.write_data_file(save_path=r"/home/modem/inv1")

    :Example 2 --> set inverions period list from data: ::

        >>> import os
        >>> from mtpy.modeling.modem_data import Data
        >>> edi_path = r"/home/mt/edi_files"
        >>> edi_list = [os.path.join(edi_path, edi) \
                        for edi in os.listdir(edi_path)\
                        if edi.find('.edi') > 0]
        >>> md = Data(edi_list)
        >>> #get period list from an .edi file
        >>> mt_obj1 = mt.MT(edi_list[0])
        >>> inv_period_list = 1./mt_obj1.Z.freq
        >>> #invert for every third period in inv_period_list
        >>> inv_period_list = inv_period_list[np.arange(0, len(inv_period_list, 3))]
        >>> md.period_list = inv_period_list
        >>> md.write_data_file(save_path=r"/home/modem/inv1")

    :Example 3 --> change error values: ::

        >>> from mtpy.modeling.modem_data import Data
        >>> mdr = Data()
        >>> mdr.read_data_file(r"/home/modem/inv1/ModEM_Data.dat")
        >>> mdr.error_type = 'floor'
        >>> mdr.error_floor = 10
        >>> mdr.error_tipper = .03
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")

    :Example 4 --> change inversion type: ::

        >>> from mtpy.modeling.modem_data import Data
        >>> mdr = Data()
        >>> mdr.read_data_file(r"/home/modem/inv1/ModEM_Data.dat")
        >>> mdr.inv_mode = '3'
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")

    """

    def __init__(self, edi_list=None, **kwargs):
        self.edi_list = edi_list

        self.error_type = kwargs.pop('error_type', 'egbert')
        self.comp_error_type = kwargs.pop('comp_error_type', None)
        self.error_floor = kwargs.pop('error_floor', 10.0)
        self.error_value = kwargs.pop('error_value', 5.0)
        self.error_egbert = kwargs.pop('error_egbert', 3.0)
        self.error_tipper = kwargs.pop('error_tipper', .05)

        self.wave_sign_impedance = kwargs.pop('wave_sign_impedance', '+')
        self.wave_sign_tipper = kwargs.pop('wave_sign_tipper', '+')
        self.units = kwargs.pop('units', '[mV/km]/[nT]')
        self.inv_mode = kwargs.pop('inv_mode', '1')
        self.period_list = kwargs.pop('period_list', None)
        self.period_step = kwargs.pop('period_step', 1)
        self.period_min = kwargs.pop('period_min', None)
        self.period_max = kwargs.pop('period_max', None)
        self.period_buffer = kwargs.pop('period_buffer', None)
        self.max_num_periods = kwargs.pop('max_num_periods', None)
        self.data_period_list = None

        self.fn_basename = kwargs.pop('fn_basename', 'ModEM_Data.dat')
        self.save_path = kwargs.pop('save_path', os.getcwd())
        self.formatting = kwargs.pop('format', '1')

        self._rotation_angle = kwargs.pop('rotation_angle', 0.0)
        self._set_rotation_angle(self._rotation_angle)

        self._station_locations = None
        self.center_position = np.array([0.0, 0.0])
        self.epsg = kwargs.pop('epsg', None)

        self.data_array = None
        self.mt_dict = None
        self.data_fn = kwargs.pop('data_fn', 'ModEM_Data.dat')

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
                       ('z_err', (np.complex, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.complex, self._t_shape))]

        self.inv_mode_dict = {'1': ['Full_Impedance', 'Full_Vertical_Components'],
                              '2': ['Full_Impedance'],
                              '3': ['Off_Diagonal_Impedance', 'Full_Vertical_Components'],
                              '4': ['Off_Diagonal_Impedance'],
                              '5': ['Full_Vertical_Components'],
                              '6': ['Full_Interstation_TF'],
                              '7': ['Off_Diagonal_Rho_Phase']}
        self.inv_comp_dict = {'Full_Impedance': ['zxx', 'zxy', 'zyx', 'zyy'],
                              'Off_Diagonal_Impedance': ['zxy', 'zyx'],
                              'Full_Vertical_Components': ['tx', 'ty']}

        self.comp_index_dict = {'zxx': (0, 0), 'zxy': (0, 1), 'zyx': (1, 0),
                                'zyy': (1, 1), 'tx': (0, 0), 'ty': (0, 1)}

        self.header_strings = \
            ['# Created using MTpy error type {0} of {1:.0f}%, data rotated {2:.1f} deg clockwise from N\n'.format(
                self.error_type, self.error_floor, self._rotation_angle),
                '# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n'
            ]
        if self.comp_error_type is not None:
            self.header_strings += [
                "# component \'{}\' used error type: {}\n".format(comp, err_type) for comp, err_type in
                self.comp_error_type.iteritems()
            ]

        # size of a utm grid
        self._utm_grid_size_north = 888960.0
        self._utm_grid_size_east = 640000.0
        self._utm_cross = False
        self._utm_ellipsoid = 23

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
                       ('z_err', (np.complex, self._z_shape)),
                       ('tip', (np.complex, self._t_shape)),
                       ('tip_err', (np.complex, self._t_shape))]

    def _set_header_string(self):
        """
        reset the header string for file
        """

        h_str = '# Created using MTpy error type {0} of {1:.0f}%, data rotated {2:.1f} deg clockwise from N\n'
        if self.error_type == 'egbert':
            self.header_strings[0] = h_str.format(self.error_type,
                                                  self.error_egbert,
                                                  self._rotation_angle)
        elif self.error_type == 'floor':
            self.header_strings[0] = h_str.format(self.error_type,
                                                  self.error_floor,
                                                  self._rotation_angle)
        elif self.error_type == 'value':
            self.header_strings[0] = h_str.format(self.error_type,
                                                  self.error_value,
                                                  self._rotation_angle)

    def get_mt_dict(self):
        """
        get mt_dict from edi file list
        """

        if self.edi_list is None:
            raise DataError('edi_list is None, please input a list of '
                            '.edi files containing the full path')

        if len(self.edi_list) == 0:
            raise DataError('edi_list is empty, please input a list of '
                            '.edi files containing the full path')

        self.mt_dict = {}
        for edi in self.edi_list:
            mt_obj = mt.MT(edi)
            self.mt_dict[mt_obj.station] = mt_obj

    def project_sites(self):
        """
        function to project sites from lat/long to eastings/northing.
        no dependency on external projection modules (e.g. pyproj) but
        limited flexibility for projection.
        """

        utm_zones_dict = {'M': 9, 'L': 8, 'K': 7, 'J': 6, 'H': 5, 'G': 4, 'F': 3,
                          'E': 2, 'D': 1, 'C': 0, 'N': 10, 'P': 11, 'Q': 12, 'R': 13,
                          'S': 14, 'T': 15, 'U': 16, 'V': 17, 'W': 18, 'X': 19}

        # --> need to convert lat and lon to east and north
        for c_arr in self.data_array:
            if c_arr['lat'] != 0.0 and c_arr['lon'] != 0.0:
                c_arr['zone'], c_arr['east'], c_arr['north'] = \
                    utm2ll.LLtoUTM(self._utm_ellipsoid,
                                   c_arr['lat'],
                                   c_arr['lon'])

        # --> need to check to see if all stations are in the same zone
        utm_zone_list = list(set(self.data_array['zone']))

        # if there are more than one zone, figure out which zone is the odd
        # ball
        utm_zone_dict = dict([(utmzone, 0) for utmzone in utm_zone_list])

        if len(utm_zone_list) != 1:
            self._utm_cross = True
            for c_arr in self.data_array:
                utm_zone_dict[c_arr['zone']] += 1

            # flip keys and values so the key is the number of zones and
            # the value is the utm zone
            utm_zone_dict = dict([(utm_zone_dict[key], key)
                                  for key in utm_zone_dict.keys()])

            # get the main utm zone as the one with the most stations in it
            main_utm_zone = utm_zone_dict[max(utm_zone_dict.keys())]

            # Get a list of index values where utm zones are not the
            # same as the main zone
            diff_zones = np.where(self.data_array['zone'] != main_utm_zone)[0]
            for c_index in diff_zones:
                c_arr = self.data_array[c_index]
                c_utm_zone = c_arr['zone']

                print '{0} utm_zone is {1} and does not match {2}'.format(
                    c_arr['station'], c_arr['zone'], main_utm_zone)

                zone_shift = 1 - abs(utm_zones_dict[c_utm_zone[-1]] -
                                     utm_zones_dict[main_utm_zone[-1]])

                # --> check to see if the zone is in the same latitude
                # if odd ball zone is north of main zone, add 888960 m
                if zone_shift > 1:
                    north_shift = self._utm_grid_size_north * zone_shift
                    print ('--> adding {0:.2f}'.format(north_shift) +
                           ' meters N to place station in ' +
                           'proper coordinates relative to all other ' +
                           'staions.')
                    c_arr['north'] += north_shift

                # if odd ball zone is south of main zone, subtract 88960 m
                elif zone_shift < -1:
                    north_shift = self._utm_grid_size_north * zone_shift
                    print ('--> subtracting {0:.2f}'.format(north_shift) +
                           ' meters N to place station in ' +
                           'proper coordinates relative to all other ' +
                           'staions.')
                    c_arr['north'] -= north_shift

                # --> if zone is shifted east or west
                if int(c_utm_zone[0:-1]) > int(main_utm_zone[0:-1]):
                    east_shift = self._utm_grid_size_east * \
                                 abs(int(c_utm_zone[0:-1]) - int(main_utm_zone[0:-1]))
                    print ('--> adding {0:.2f}'.format(east_shift) +
                           ' meters E to place station in ' +
                           'proper coordinates relative to all other ' +
                           'staions.')
                    c_arr['east'] += east_shift
                elif int(c_utm_zone[0:-1]) < int(main_utm_zone[0:-1]):
                    east_shift = self._utm_grid_size_east * \
                                 abs(int(c_utm_zone[0:-1]) - int(main_utm_zone[0:-1]))
                    print ('--> subtracting {0:.2f}'.format(east_shift) +
                           ' meters E to place station in ' +
                           'proper coordinates relative to all other ' +
                           'staions.')
                    c_arr['east'] -= east_shift

    def project_sites_pyproj(self):
        """
        project site locations from lat/long to eastings/northings (defined by
        epsg in data object). Uses pyproj so this needs to be installed before
        this module will work.


        """

        if self.epsg not in constants.epsg_dict.keys():
            self.epsg = None

        if self.epsg is None:
            logger.debug("can't project, invalid (or no) epsg provided")
            return

        for c_arr in self.data_array:
            if c_arr['lat'] != 0.0 and c_arr['lon'] != 0.0:
                c_arr['zone'] = constants.epsg_dict[self.epsg][1]
                c_arr['east'], c_arr['north'] = \
                    utm2ll.project(c_arr['lon'], c_arr['lat'],
                                   4326, self.epsg)

    def project_xy(self, x, y, epsg_from=None, epsg_to=4326):
        """
        project some xy points
        """
        if epsg_from is None:
            epsg_from = self.epsg

        return np.array(utm2ll.project(x, y, epsg_from, epsg_to))

    def get_relative_station_locations(self):
        """
        get station locations from edi files and project to local coordinates

        ..note:: There are two options for projection method. If pyproj is
                 installed, you can use the method that uses pyproj. In this
                 case, specify the epsg number as an attribute to the model
                 object or when setting it up. The epsg can generally be found
                 through a google search. If epsg is specified then **all**
                 sites are projected to that epsg. It is up to the user to
                 make sure all sites are in the bounds of projection.
                 **note** epsg 3112 (Geoscience Australia Lambert) covers all
                 of Australia but may cause signficiant rotation at some
                 locations.

                ***If pyproj is not used:***
                If the survey steps across multiple UTM zones, then a
                 distance will be added to the stations to place them in
                 the correct location.  This distance is
                 _utm_grid_size_north and _utm_grid_size_east. You should
                 these parameters to place the locations in the proper spot
                 as grid distances and overlaps change over the globe.

        """

        if self.epsg is not None:
            use_pyproj = True
        else:
            use_pyproj = False

        if use_pyproj:
            self.project_sites_pyproj()
        else:
            self.project_sites()

        # define center of the grid in east/north coordinates
        self.center_position_EN = 0.5 * np.array([self.data_array['east'].min() + self.data_array['east'].max(),
                                                  self.data_array['north'].min() + self.data_array['north'].max()])

        # update center_position by projecting center xy
        # try to use pyproj if desired, if not then have to use inbuilt
        # projection module but may give bad results if crossing more than one zone
        self.center_position = self.project_xy(*self.center_position_EN)

        # get center position of the stations in lat and lon
        self.center_position2 = 0.5 * np.array([self.data_array['lon'].min() + self.data_array['lon'].max(),
                                                self.data_array['lat'].min() + self.data_array['lat'].max()])

        logger.debug("compare center positions: %s, %s", self.center_position, self.center_position2)

        # remove the average distance to get coordinates in a relative space
        self.data_array['rel_east'] = self.data_array[
                                          'east'] - self.center_position_EN[0]
        self.data_array['rel_north'] = self.data_array[
                                           'north'] - self.center_position_EN[1]

        # --> rotate grid if necessary
        # to do this rotate the station locations because ModEM assumes the
        # input mesh is a lateral grid.
        # needs to be 90 - because North is assumed to be 0 but the rotation
        # matrix assumes that E is 0.
        if self.rotation_angle != 0:
            cos_ang = np.cos(np.deg2rad(self.rotation_angle))
            sin_ang = np.sin(np.deg2rad(self.rotation_angle))
            rot_matrix = np.matrix(np.array([[cos_ang, sin_ang],
                                             [-sin_ang, cos_ang]]))

            coords = np.array([self.data_array['rel_east'],
                               self.data_array['rel_north']])

            # rotate the relative station locations
            new_coords = np.array(np.dot(rot_matrix, coords))

            self.data_array['rel_east'][:] = new_coords[0, :]
            self.data_array['rel_north'][:] = new_coords[1, :]

            print 'Rotated stations by {0:.1f} deg clockwise from N'.format(
                self.rotation_angle)

    def get_period_list(self):
        """
        make a period list to invert for
        """
        if self.mt_dict is None:
            self.get_mt_dict()

        if self.period_list is None:
            raise DataError('Need to input period_min, period_max, '
                            'max_num_periods or a period_list')
        else:
            print '-' * 50
            print ('Inverting for these periods:', len(self.period_list))
            for per in self.period_list:
                print '     {0:<12.6f}'.format(per)
            print '-' * 50

            return  # finished

        # FZ: why here ? log space interpolation???
        # YG: NOTE the code below never reached
        data_period_list = []
        for s_key in sorted(self.mt_dict.keys()):
            mt_obj = self.mt_dict[s_key]
            data_period_list.extend(list(1. / mt_obj.Z.freq))

        self.data_period_list = np.array(sorted(list(set(data_period_list)),
                                                reverse=False))

        if self.period_min is not None:
            if self.period_max is None:
                raise DataError('Need to input period_max')
        if self.period_max is not None:
            if self.period_min is None:
                raise DataError('Need to input period_min')
        if self.period_min is not None and self.period_max is not None:
            if self.max_num_periods is None:
                raise DataError('Need to input number of periods to use')

            min_index = np.where(self.data_period_list >=
                                 self.period_min)[0][0]
            max_index = np.where(self.data_period_list <=
                                 self.period_max)[0][-1]

            pmin = np.log10(self.data_period_list[min_index])
            pmax = np.log10(self.data_period_list[max_index])
            self.period_list = np.logspace(
                pmin, pmax, num=self.max_num_periods)

            print '-' * 50
            print ('Inverting for periods:', len(self.period_list))
            for per in self.period_list:
                print '     {0:<12.6f}'.format(per)
            print '-' * 50

    def _set_rotation_angle(self, rotation_angle):
        """
        on set rotation angle rotate mt_dict and data_array,
        """
        if self._rotation_angle == rotation_angle:
            return

        new_rotation_angle = -self._rotation_angle + rotation_angle

        if new_rotation_angle == 0:
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
            mt_obj.Z.rotate(new_rotation_angle)
            mt_obj.Tipper.rotate(new_rotation_angle)

        print 'Data rotated to align with {0:.1f} deg clockwise from N'.format(
            self._rotation_angle)

        print '*' * 70
        print '   If you want to rotate station locations as well use the'
        print '   command Data.get_relative_station_locations() '
        print '   if stations have not already been rotated in Model'
        print '*' * 70

        self._fill_data_array()

    def _get_rotation_angle(self):
        return self._rotation_angle

    rotation_angle = property(fget=_get_rotation_angle,
                              fset=_set_rotation_angle,
                              doc="""Rotate data assuming N=0, E=90""")

    def _initialise_empty_data_array(self, stationlocations, period_list,
                                     location_type='LL', stationnames=None):
        """
        create an empty data array to create input files for forward modelling
        station locations is an array containing x,y coordinates of each station
        (shape = (number_of_stations,2))
        period_list = list of periods to model
        location_type = 'LL' or 'EN' - longitude/latitude or easting/northing

        """
        self.period_list = period_list.copy()
        nf = len(self.period_list)
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(len(stationlocations), dtype=self._dtype)
        if location_type == 'LL':
            self.data_array['lon'] = stationlocations[:, 0]
            self.data_array['lat'] = stationlocations[:, 1]
        else:
            self.data_array['east'] = stationlocations[:, 0]
            self.data_array['north'] = stationlocations[:, 1]

            # set non-zero values to array (as zeros will be deleted)
        if self.inv_mode in '12':
            self.data_array['z'][:] = 100. + 100j
            self.data_array['z_err'][:] = 1e15
        if self.inv_mode == '1':
            self.data_array['tip'][:] = 0.1 + 0.1j
            self.data_array['tip_err'][:] = 1e15

        # set station names
        if stationnames is not None:
            if len(stationnames) != len(stationnames):
                stationnames = None

        if stationnames is None:
            stationnames = ['st%03i' %
                            ss for ss in range(len(stationlocations))]
        self.data_array['station'] = stationnames

        self.get_relative_station_locations()

    def _fill_data_array(self, new_edi_dir=None, use_original_freq=True):
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
            logger.debug("mt_dict key: %s and ii= %s", s_key, ii)  # s_key is station name

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
                logger.debug("populate d_array from edi file mt_obj of station %s !!!", s_key)
                self.data_array[ii]['station'] = mt_obj.station
                self.data_array[ii]['lat'] = mt_obj.lat
                self.data_array[ii]['lon'] = mt_obj.lon
                self.data_array[ii]['east'] = mt_obj.east
                self.data_array[ii]['north'] = mt_obj.north
                self.data_array[ii]['elev'] = mt_obj.elev
                try:  # this block will raise exception. see get_relative_station_locations()
                    self.data_array[ii]['rel_east'] = mt_obj.grid_east  # does not exist attribute grid.east
                    self.data_array[ii]['rel_north'] = mt_obj.grid_north
                    rel_distance = False
                except AttributeError:
                    logger.debug("skipping - self.data_array[ii]['rel_east'] was not assigned here !!!")
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
            logger.debug("station_name and its original period: %s %s %s", mt_obj.station, len(mt_obj.Z.freq),
                         1.0 / mt_obj.Z.freq)
            logger.debug("station_name and interpolation period: %s %s %s", mt_obj.station, len(interp_periods),
                         interp_periods)

            # default: use_original_freq = True, each MT station edi file will use it's own frequency-filtered.
            # no new freq in the output modem.dat file. select those freq of mt_obj according to interp_periods
            if use_original_freq:
                interp_periods = self.filter_periods(mt_obj, interp_periods)
                logger.debug("station_name and selected/filtered periods: %s, %s, %s", mt_obj.station,
                             len(interp_periods), interp_periods)
                # in this case the below interpolate_impedance_tensor function will degenerate into a same-freq set.

            if len(interp_periods) > 0:  # not empty
                interp_z, interp_t = mt_obj.interpolate_impedance_tensor(1. / interp_periods)  # ,bounds_error=False)
                for kk, ff in enumerate(interp_periods):
                    jj = np.where(self.period_list == ff)[0][0]
                    self.data_array[ii]['z'][jj] = interp_z.z[kk, :, :]
                    self.data_array[ii]['z_err'][jj] = interp_z.z_err[kk, :, :]

                    if mt_obj.Tipper.tipper is not None:
                        self.data_array[ii]['tip'][jj] = interp_t.tipper[kk, :, :]
                        self.data_array[ii]['tip_err'][jj] = \
                            interp_t.tipper_err[kk, :, :]

                # FZ: try to output a new edi files. Compare with original edi?
                if new_edi_dir is not None:
                    new_edifile = os.path.join(new_edi_dir, mt_obj.station + '.edi')
                    mt_obj.write_edi_file(new_fn=new_edifile, new_Z=interp_z, new_Tipper=interp_t)
            else:
                pass

        if rel_distance is False:
            self.get_relative_station_locations()

        return

    def filter_periods(self, mt_obj, per_array):
        """Select the periods of the mt_obj that are in per_array.
        used to do original freq inversion.

        :param mt_obj:
        :param per_array:
        :return: array of selected periods (subset) of the mt_obj
        """

        new_per = []

        mt_per = 1.0 / mt_obj.Z.freq
        for p in mt_per:
            for p2 in per_array:
                if abs(p - p2) < 0.00000001:  # Be aware of floating error if use ==
                    new_per.append(p)

        return np.array(new_per)

    def _set_station_locations(self, station_locations):
        """
        take a station_locations array and populate data_array
        """

        if self.data_array is None:
            self.get_mt_dict()
            self.get_period_list()
            self._fill_data_array()

        for s_arr in station_locations:
            try:
                d_index = np.where(self.data_array['station'] ==
                                   s_arr['station'])[0][0]
            except IndexError:
                logger.debug('Could not find {0} in data_array'.format(s_arr['station']))
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
                                             'north', 'east', 'elev', 'zone',
                                             'rel_north', 'rel_east']]
        return station_locations

    station_locations = property(_get_station_locations,
                                 _set_station_locations,
                                 doc="""location of stations""")

    def write_data_file(self, save_path=None, fn_basename=None,
                        rotation_angle=None, compute_error=True,
                        fill=True):
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
            **data_fn** : string full path to created data file
        """

        if save_path is not None:
            self.save_path = save_path
        if fn_basename is not None:
            self.fn_basename = fn_basename

        self.data_fn = os.path.join(self.save_path, self.fn_basename)

        if fill:
            self.get_period_list()

        # rotate data if required
        if rotation_angle is not None:
            self.rotation_angle = rotation_angle

        # be sure to fill in data array
        if fill is True:
            new_edi_dir = os.path.join(self.save_path, 'new_edis')  # output edi files according to selected periods
            if not os.path.exists(new_edi_dir):
                os.mkdir(new_edi_dir)
            self._fill_data_array(new_edi_dir=new_edi_dir)
            # get relative station locations in grid coordinates
            self.get_relative_station_locations()

        # reset the header string to be informational
        self._set_header_string()

        # number of periods - subtract periods with all zero components
        nper = len(np.where(np.mean(
            np.mean(np.mean(np.abs(self.data_array['z']), axis=0), axis=1), axis=1) > 0)[0])

        dlines = []
        for inv_mode in self.inv_mode_dict[self.inv_mode]:
            # dlines.append(self.header_strings[0])
            # dlines.append(self.header_strings[1])
            dlines += self.header_strings
            dlines.append('> {0}\n'.format(inv_mode))

            if inv_mode.find('Impedance') > 0:
                dlines.append('> exp({0}i\omega t)\n'.format(
                    self.wave_sign_impedance))
                dlines.append('> {0}\n'.format(self.units))
            elif inv_mode.find('Vertical') >= 0:
                dlines.append('> exp({0}i\omega t)\n'.format(
                    self.wave_sign_tipper))
                dlines.append('> []\n')
            dlines.append('> 0.00\n')  # oriention, need to add at some point
            dlines.append('> {0: >10.6f} {1:>10.6f}\n'.format(
                self.center_position[1], self.center_position[0]))  # (lat,long) correct order
            dlines.append('> {0} {1}\n'.format(nper,
                                               self.data_array['z'].shape[0]))

            # YG: create new list for sorting data
            data_lines = []

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

                        # get the value for that compenent at that frequency
                        zz = self.data_array[ss][c_key][ff, z_ii, z_jj]
                        if zz.real != 0.0 and zz.imag != 0.0 and \
                                        zz.real != 1e32 and zz.imag != 1e32:
                            if self.formatting == '1':
                                per = '{0:<12.5e}'.format(self.period_list[ff])
                                sta = '{0:>7}'.format(
                                    self.data_array[ss]['station'])
                                lat = '{0:> 9.3f}'.format(
                                    self.data_array[ss]['lat'])
                                lon = '{0:> 9.3f}'.format(
                                    self.data_array[ss]['lon'])
                                eas = '{0:> 12.3f}'.format(
                                    self.data_array[ss]['rel_east'])
                                nor = '{0:> 12.3f}'.format(
                                    self.data_array[ss]['rel_north'])
                                ele = '{0:> 12.3f}'.format(
                                    -self.data_array[ss]['elev'])
                                com = '{0:>4}'.format(comp.upper())
                                if self.units == 'ohm':
                                    rea = '{0:> 14.6e}'.format(zz.real / 796.)
                                    ima = '{0:> 14.6e}'.format(zz.imag / 796.)
                                else:
                                    rea = '{0:> 14.6e}'.format(zz.real)
                                    ima = '{0:> 14.6e}'.format(zz.imag)

                            elif self.formatting == '2':
                                per = '{0:<14.6e}'.format(self.period_list[ff])
                                sta = '{0:<10}'.format(
                                    self.data_array[ss]['station'])
                                lat = '{0:> 14.6f}'.format(
                                    self.data_array[ss]['lat'])
                                lon = '{0:> 14.6f}'.format(
                                    self.data_array[ss]['lon'])
                                eas = '{0:> 12.3f}'.format(
                                    self.data_array[ss]['rel_east'])
                                nor = '{0:> 15.3f}'.format(
                                    self.data_array[ss]['rel_north'])
                                ele = '{0:> 10.3f}'.format(
                                    -self.data_array[ss]['elev'])
                                com = '{0:>12}'.format(comp.upper())
                                if self.units == 'ohm':
                                    rea = '{0:> 17.6e}'.format(zz.real / 796.)
                                    ima = '{0:> 17.6e}'.format(zz.imag / 796.)
                                else:
                                    rea = '{0:> 17.6e}'.format(zz.real)
                                    ima = '{0:> 17.6e}'.format(zz.imag)
                            if compute_error:
                                # compute relative error
                                if comp.find('t') == 0:
                                    abs_err = self._vertical_components_error_floor(ff, c_key, ss, z_ii, z_jj)
                                elif comp.find('z') == 0:
                                    comp_error_type = self.comp_error_type[comp] if self.comp_error_type is not None \
                                                                                    and comp in self.comp_error_type \
                                        else self.error_type
                                    # use the type from comp_error_type if the type of comp is specified otherwise
                                    # use self.error_type

                                    if comp_error_type == 'floor':
                                        abs_err = self._impedance_components_error_floor(c_key, ff, ss, z_ii, z_jj, zz)
                                    elif comp_error_type == 'value':
                                        abs_err = self._impedance_components_error_value(zz)
                                    elif comp_error_type == 'egbert':
                                        abs_err = self._impedance_components_error_egbert(ff, ss)
                                    elif comp_error_type == 'floor_egbert':
                                        abs_err = self._impedance_components_error_floor_egbert(c_key, ff, ss, z_ii,
                                                                                                z_jj)
                                    elif comp_error_type == 'stddev':
                                        abs_err = self._impedance_components_error_stddev(c_key, ss, z_ii, z_jj)
                                    elif comp_error_type == 'sqr':
                                        abs_err = self._impedance_components_error_sqr(c_key, ff, ss, z_ii, z_jj)
                                    elif comp_error_type == 'meansqr':
                                        abs_err = self._impedance_components_error_meansqr(c_key, ss, z_ii, z_jj)

                                if abs_err == 0.0:
                                    abs_err = 1e3
                                    print ('error at {0} is 0 for period {1}'.format(
                                        sta, per) + 'set to 1e3')
                                    if self.units == 'ohm':
                                        abs_err /= 796.

                            else:
                                abs_err = self.data_array[ss][
                                    c_key + '_err'][ff, z_ii, z_jj].real
                                if ((c_key.find('z') >= 0) and (self.units == 'ohm')):
                                    abs_err /= 796.

                            abs_err = '{0:> 14.6e}'.format(abs(abs_err))
                            # make sure that x==north, y==east, z==+down
                            # YG: populate data that needs to be sorted
                            data_lines.append([per, sta, lat, lon, nor, eas, ele,
                                               com, rea, ima, abs_err, '\n'])
                            # dline = ''.join([per, sta, lat, lon, nor, eas, ele,
                            #                  com, rea, ima, abs_err, '\n'])
                            # dlines.append(dline)

            # YG: sort by station then by period
            data_lines = sorted(data_lines, key=lambda line: (line[1], float(line[0])))
            # add data table to output lines
            dlines += [''.join(row) for row in data_lines]

        # self.data_fn = MTfh.make_unique_filename(self.data_fn) # make a unique file in debug stage
        # if os.path.exists(self.data_fn):
        #     data_fn1= self.data_fn[:-3]+"OLD"
        #     os.rename(self.data_fn, data_fn1)

        # dfid = file(self.data_fn, 'w')
        dfid = file(self.data_fn, 'w')
        dfid.writelines(dlines)
        dfid.close()

        # write epsg and center position to a file, if they exist
        np.savetxt(op.join(self.save_path, 'center_position.txt'),
                   self.center_position_EN, fmt='%.1f')
        np.savetxt(op.join(self.save_path, 'epsg.txt'),
                   np.array([self.epsg]), fmt='%1i')

        logger.debug('Wrote ModEM data file to %s', self.data_fn)

        return self.data_fn

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

    def _impedance_components_error_egbert(self, ff, ss):
        d_zxy = self.data_array[
            ss]['z'][ff, 0, 1]
        d_zyx = self.data_array[
            ss]['z'][ff, 1, 0]
        abs_err = np.sqrt(abs(d_zxy * d_zyx)) * \
                  self.error_egbert / 100.
        return abs_err

    def _impedance_components_error_value(self, zz):
        abs_err = abs(zz) * \
                  self.error_value / 100.
        return abs_err

    def _impedance_components_error_floor_egbert(self, c_key, ff, ss, z_ii, z_jj):
        # abs_err = self.data_array[ss][
        #     c_key + '_err'][ff, z_ii, z_jj]
        d_zxy = self.data_array[ss]['z'][ff, 0, 1]
        d_zyx = self.data_array[ss]['z'][ff, 1, 0]
        # if abs_err < np.sqrt(abs(d_zxy * d_zyx)) * self.error_egbert / 100.:
        #     abs_err = np.sqrt(
        #         abs(d_zxy * d_zyx)) * self.error_egbert / 100.
        abs_err = max(
            np.sqrt(abs(d_zxy * d_zyx)) * self.error_egbert / 100.,
            self.data_array[ss][c_key + '_err'][ff, z_ii, z_jj]
        )
        return abs_err

    def _impedance_components_error_floor(self, c_key, ff, ss, z_ii, z_jj, zz):
        rel_err = max(
            self.error_floor / 100.,
            self.data_array[ss][c_key + '_err'][ff, z_ii, z_jj] / abs(zz)
        )
        # if rel_err < self.error_floor / 100.:
        #     rel_err = self.error_floor / 100.
        abs_err = rel_err * abs(zz)
        return abs_err

    def _vertical_components_error_floor(self, ff, c_key, ss, z_ii, z_jj):
        if 'floor' in self.error_type:
            # abs_err = max(self.error_tipper,
            #               self.data_array[ss][c_key + '_err'][ff, 0, z_ii])  # this may be wrong as z_ii is always 0
            abs_err = max(self.error_tipper, self.data_array[ss][c_key + '_err'][ff, z_ii, z_jj])
        else:
            abs_err = self.error_tipper
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

            >>> from mtpy.modeling.modem_data import Data as modem
            >>> mdr = modem.Data()
            >>> mdr.convert_ws3dinv_data_file(r"/home/ws3dinv/inv1/WSData.dat",
                    station_fn=r"/home/ws3dinv/inv1/WS_Station_Locations.txt")
        """

        if os.path.isfile(ws_data_fn) == False:
            raise ws.WSInputError(
                'Did not find {0}, check path'.format(ws_data_fn))

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
        linecount = 0
        print "reading data lines"
        for dline in dlines:
            linecount += 1
            if dline.find('#') == 0:
                header_list.append(dline.strip())
            elif dline.find('>') == 0:
                metadata_list.append(dline[1:].strip())
                if dline.lower().find('ohm') > 0:
                    self.units = 'ohm'
                if dline.lower().find('mv') > 0:
                    self.units = ' [mV/km]/[nT]'
                if dline.lower().find('vertical') > 0:
                    read_tipper = True
                    read_impedance = False
                elif dline.lower().find('impedance') > 0:
                    read_impedance = True
                    read_tipper = False
                if linecount == 7:
                    print "getting center position", dline
                    self.center_position = [
                        float(val) for val in dline.strip().replace('>', '').split()]
                    print self.center_position
                if dline.find('exp') > 0:
                    if read_impedance is True:
                        self.wave_sign_impedance = dline[dline.find('(') + 1]
                    elif read_tipper is True:
                        self.wave_sign_tipper = dline[dline.find('(') + 1]
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
                    print ('Set rotation angle to {0:.1f} '.format(
                        self._rotation_angle) + 'deg clockwise from N')
                except ValueError:
                    pass

        self.period_list = np.array(sorted(set(period_list)))
        station_list = sorted(set(station_list))

        # make a period dictionary to with key as period and value as index
        period_dict = dict([(per, ii)
                            for ii, per in enumerate(self.period_list)])

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
            if tf_dict[dd[1]] == False:
                data_dict[dd[1]].lat = dd[2]
                data_dict[dd[1]].lon = dd[3]
                data_dict[dd[1]].grid_north = dd[4]
                data_dict[dd[1]].grid_east = dd[5]
                data_dict[dd[1]].grid_elev = dd[6]
                data_dict[dd[1]].station = dd[1]
                tf_dict[dd[1]] = True
            # fill in the impedance tensor with appropriate values
            if dd[7].find('Z') == 0:
                z_err = dd[10]
                if self.wave_sign_impedance == '+':
                    z_value = dd[8] + 1j * dd[9]
                elif self.wave_sign_impedance == '-':
                    z_value = dd[8] - 1j * dd[9]

                if self.units == 'ohm':
                    z_value *= 796.
                    z_err *= 796.

                data_dict[dd[1]].Z.z[p_index, ii, jj] = z_value
                data_dict[dd[1]].Z.z_err[p_index, ii, jj] = z_err
            # fill in tipper with appropriate values
            elif dd[7].find('T') == 0:
                if self.wave_sign_tipper == '+':
                    data_dict[dd[1]].Tipper.tipper[
                        p_index, ii, jj] = dd[8] + 1j * dd[9]
                elif self.wave_sign_tipper == '-':
                    data_dict[dd[1]].Tipper.tipper[
                        p_index, ii, jj] = dd[8] - 1j * dd[9]
                data_dict[dd[1]].Tipper.tipper_err[p_index, ii, jj] = dd[10]

        # make mt_dict an attribute for easier manipulation later
        self.mt_dict = data_dict

        ns = len(self.mt_dict.keys())
        nf = len(self.period_list)
        self._set_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(ns, dtype=self._dtype)

        # Be sure to caclulate invariants and phase tensor for each station
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]

            self.mt_dict[s_key].zinv.compute_invariants()
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

            self.data_array[ii]['tip'][:] = mt_obj.Tipper.tipper
            self.data_array[ii]['tip_err'][:] = mt_obj.Tipper.tipper_err

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

        if vtk_save_path is not None:
            vtk_fn = os.path.join(self.save_path, vtk_fn_basename)
        else:
            vtk_fn = os.path.join(vtk_save_path, vtk_fn_basename)

        pointsToVTK(vtk_fn,
                    self.station_locations['rel_north'],
                    self.station_locations['rel_east'],
                    -self.station_locations['elev'],
                    pointData={'elevation': self.station_locations['elev']})

        logger.debug('Wrote file to %s', vtk_fn)


# ====================================================================

def select_periods(edifiles_list, percent=10.0):
    """
    FZ: Use edi_collection to analyse the whole set of EDI files
    :param edifiles:
    :return:
    """
    import matplotlib.pyplot as plt

    edis_obj = EdiCollection(edifiles_list)

    uniq_period_list = edis_obj.all_unique_periods  # filtered list of periods ?
    logger.info("Number of Unique Periods= %s", len(uniq_period_list))

    plt.hist(edis_obj.mt_periods, bins=uniq_period_list)
    # plt.hist(edis_obj.mt_periods, bins=1000)
    plt.title("Histogram with uniq_periods bins")
    plt.xlabel("Periods")
    plt.ylabel("Occurance in number of MT stations")
    plt.show()

    # 1 ASK user to input a Pmin and Pmax

    # 2 percetage stats
    # select commonly occured frequencies from all stations.
    # This could miss some slightly varied frquencies in the middle range.
    select_period_list = np.array(edis_obj.get_periods_by_stats(percentage=percent))
    logger.info("Number of Selected Periods= %s", len(select_period_list))

    return select_period_list


####################################################################################
if __name__ == '__main__':
    """ Quick test of this script for selection of inversion periods, writing a data file,
    and output new effective edi files in the output directory. No topo file is used.
    USAGE examples:
    python mtpy/modeling/modem_data.py tests/data/edifiles/  /e/tmp/modem_data_test
    python mtpy/modeling/modem_data.py /e/Data/MT_Datasets/WenPingJiang_EDI /e/tmp/WenPingTest
    """
    # epsg to project to. Google epsg 'your projection'
    # epsg_code = 28354  # UTM 54
    epsg_code = 3112  # GA-LCC

    if len(sys.argv) < 3:
        print("USAGE: %s  path2edifiles path2outdir" % sys.argv[0])
        sys.exit(1)
    else:
        edipath = sys.argv[1]  # edi files to be inversioned
        outputdir = sys.argv[2]  # path to save to

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    edi_list = glob.glob(edipath + '/*.edi')

    if edi_list is None or (edi_list) < 1:
        print("Error: No edi files found in the dir %s" % edipath)
        sys.exit(2)

    # period list (can take periods from one of the edi files, or just specify
    # periods directly using the logspace function (commented out))

    # eo = mtedi.Edi(edi_list[0])  # this may miss some periods?
    # period_list = 1. / eo.Z.freq # period_list = np.logspace(-3,3)

    period_list = select_periods(edi_list, percent=30.0)

    datob = Data(edi_list=edi_list,
                 inv_mode='1',
                 period_list=period_list,
                 epsg=epsg_code,
                 error_type='floor',
                 error_floor=10)
    # period_buffer=0.000001)

    datob.write_data_file(save_path=outputdir)
