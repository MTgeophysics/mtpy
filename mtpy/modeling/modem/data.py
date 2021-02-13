"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch
# revised by JP 2021 adding functionality and updating.

"""
from __future__ import print_function
import csv
import numpy as np
from pathlib import Path
from copy import deepcopy

import matplotlib.pyplot as plt

import mtpy.analysis.pt as pt
from mtpy.core import mt as mt
from mtpy.core import z as mtz
from mtpy.modeling import ws3dinv as ws
from mtpy.utils import gis_tools as gis_tools
from mtpy.utils.mtpy_logger import get_mtpy_logger

from mtpy.modeling.modem.exception import ModEMError, DataError
from mtpy.modeling.modem.station import Stations


try:
    from pyevtk.hl import pointsToVTK
except ImportError:
    print("If you want to write a vtk file for 3d viewing, you need to install pyevtk")


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
    
    :param edi_list: list of edi files to read

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


    :Example 1 --> create inversion period list: ::

        >>> from pathlib import Path
        >>> import mtpy.modeling.modem as modem
        >>> edi_path = Path(r"/home/mt/edi_files")
        >>> edi_list = list(edi_path.glob("*.edi"))
        >>> md = modem.Data(edi_list, period_min=.1, period_max=300,\
        >>> ...             max_num_periods=12)
        >>> md.write_data_file(save_path=r"/home/modem/inv1")
        >>> md
        

    :Example 2 --> set inverions period list from data: ::

        >>> md = modem.Data(edi_list)
        >>> #get period list from an .edi file
        >>> inv_period_list = 1./md.mt_dict["mt01"].Z.freq
        >>> #invert for every third period in inv_period_list
        >>> inv_period_list = inv_period_list[np.arange(0, len(inv_period_list, 3))]
        >>> md.period_list = inv_period_list
        >>> md.write_data_file(save_path=r"/home/modem/inv1")

    :Example 3 --> change error values: ::

        >>> mdr.error_type = 'floor'
        >>> mdr.error_floor = 10
        >>> mdr.error_tipper = .03
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")

    :Example 4 --> change inversion type: ::

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

        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self.mt_dict = None
        self.edi_list = None

        self.error_type_z = "egbert_floor"
        self.error_value_z = 5.0

        self.error_value_tipper = 0.05
        self.error_type_tipper = "abs"

        self.wave_sign_impedance = "+"
        self.wave_sign_tipper = "+"
        self.units = "[mV/km]/[nT]"
        self.inv_mode = "1"

        self.period_list = None
        # self.period_step = 1
        self.period_min = None
        self.period_max = None
        self.period_buffer = None
        self.max_num_periods = None
        self.data_period_list = None

        self.data_fn = "ModEM_Data.dat"
        self.save_path = Path.cwd()

        self.formatting = "1"

        self._rotation_angle = 0.0

        self.data_array = None
        self.model_utm_zone = None
        self.model_epsg = None

        self._center_lat = None
        self._center_lon = None
        self._center_elev = None

        self.inv_mode_dict = {
            "1": ["Full_Impedance", "Full_Vertical_Components"],
            "2": ["Full_Impedance"],
            "3": ["Off_Diagonal_Impedance", "Full_Vertical_Components"],
            "4": ["Off_Diagonal_Impedance"],
            "5": ["Full_Vertical_Components"],
            "6": ["Full_Interstation_TF"],
            "7": ["Off_Diagonal_Rho_Phase"],
        }
        self.inv_comp_dict = {
            "Full_Impedance": ["zxx", "zxy", "zyx", "zyy"],
            "Off_Diagonal_Impedance": ["zxy", "zyx"],
            "Full_Vertical_Components": ["tx", "ty"],
        }

        self.comp_index_dict = {
            "zxx": (0, 0),
            "zxy": (0, 1),
            "zyx": (1, 0),
            "zyy": (1, 1),
            "tx": (0, 0),
            "ty": (0, 1),
        }

        self.header_string = " ".join(
            [
                "# Period(s)",
                "Code",
                "GG_Lat",
                "GG_Lon",
                "X(m)",
                "Y(m)",
                "Z(m)",
                "Component",
                "Real",
                "Imag",
                "Error\n",
            ]
        )

        # fill mt dict if edi list is not None
        if edi_list is not None:
            self.mt_dict = self.make_mt_dict(edi_list)

        for key, value in kwargs.items():
            # have to set rotation angle after period list has been set
            if key != "rotation_angle":
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    self.logger.debug(f"Argument {key}={value} will not be set")

        # update period buffer to a default value if it is invalid
        if self.period_buffer is not None:
            try:
                self.period_buffer = float(self.period_buffer)
                if self.period_buffer < 0.0:
                    self.logger.warning("Period buffer must be > 0, setting to None")
                    self.period_buffer = None
                # if provided value between 0 and 1, assume it was meant to be set to 1 + provided value
                elif self.period_buffer < 1.0:
                    self.logger.warning(
                        "Period buffer must be > 1, adding 1 to provided value"
                    )
                    self.period_buffer += 1.0
            except:
                self.logger.warning(
                    "Period buffer must be convertable to an integer or float, setting to None"
                )

        if "rotation_angle" in list(kwargs.keys()):
            setattr(self, "rotation_angle", kwargs["rotation_angle"])

    #            self._set_rotation_angle(self.rotation_angle)

    def __str__(self):
        lines = ["ModEM Data Object:"]
        if self.data_array is not None:
            lines += [f"\tNumber of stations: {self.data_array.shape[0]}"]
            lines += [f"\tNumber of periods:  {self.period_list.shape[0]}"]
            lines += ["\tPeriod range:  "]
            lines += [f"\t\tMin: {self.period_list.min()} s"]
            lines += [f"\t\tMax: {self.period_list.max()} s"]
            lines += [f"\tRotation angle:     {self.rotation_angle}"]
            lines += ["\tData center:        "]
            lines += [f"\t\t latitude:  {self.center_point.lat[0]:.4f} deg"]
            lines += [f"\t\t longitude: {self.center_point.lon[0]:.4f} deg"]
            lines += [f"\t\t Elevation: {self.center_point.elev[0]:.1f} m"]
            lines += [f"\t\t Easting:   {self.center_point.east[0]:.4f} m"]
            lines += [f"\t\t Northing:  {self.center_point.north[0]:.4f} m"]
            lines += [f"\t\t UTM zone:  {self.center_point.zone[0]}"]
            lines += [f"\tModel EPSG:         {self.model_epsg}"]
            lines += [f"\tModel UTM zone:     {self.model_utm_zone}"]
            lines += [f"\tImpedance data:     {self.data_array['z'].mean() != 0.0}"]
            lines += [f"\tTipper data:        {self.data_array['tip'].mean() != 0.0}"]
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @staticmethod
    def make_dtype(z_shape, t_shape):
        """
        Create data type given shapes of the impedance and tipper arrays
        
        :param z_shape: (number of periods, 2, 2)
        :type z_shape: tuple
        
        :param t_shape: (number of periods, 2, 2)
        :type t_shape: tuple
        
        """

        dtype = [
            ("station", "|U50"),
            ("lat", np.float),
            ("lon", np.float),
            ("elev", np.float),
            ("rel_east", np.float),
            ("rel_north", np.float),
            ("rel_elev", np.float),
            ("east", np.float),
            ("north", np.float),
            ("zone", "|S4"),
            ("z", (np.complex, z_shape)),
            ("z_err", (np.float, z_shape)),
            ("z_inv_err", (np.float, z_shape)),
            ("tip", (np.complex, t_shape)),
            ("tip_err", (np.float, t_shape)),
            ("tip_inv_err", (np.float, t_shape)),
        ]

        return dtype

    @staticmethod
    def get_header_string(error_type, error_value, rotation_angle):
        """
        Create the header strings
        
        # Created using MTpy calculated egbert_floor error of 5% data rotated 0.0_deg
        clockwise from N
        
        :param error_type: The method to calculate the errors
        :type error_type: string
        :param error_value: value of error or error floor
        :type error_value: float
        :param rotation_angle: angle data have been rotated by
        :type rotation_angle: float

        """

        h_str = []
        if np.atleast_1d(error_type).ndim == 2:
            h_str = "# Created using MTpy calculated {},{},{},{} ".format(
                error_type[0, 0], error_type[0, 1], error_type[1, 0], error_type[1, 1]
            )
        else:
            h_str = "# Created using MTpy calculated {} ".format(error_type)

        if np.atleast_1d(error_value).ndim == 2:
            h_str += "error floors of {0:.0f}%,{1:.0f}%,{2:.0f}%,{3:.0f}%, data rotated {4:.1f}_deg clockwise from N\n"
            return h_str.format(
                error_value[0, 0],
                error_value[0, 1],
                error_value[1, 0],
                error_value[1, 1],
                rotation_angle,
            )
        else:
            h_str += "error of {1:.0f}% data rotated {2:.1f}_deg clockwise from N\n"

            return h_str.format(error_type, error_value, rotation_angle)

    def make_mt_dict(self, edi_list=None):
        """
        Create a dictionary of :class:`mtpy.core.mt.MT` objects to pull data from
        
        :param edi_list: list of edi files to read
        :type edi_list: list of full paths to files
        
        """

        if edi_list is not None:
            self.edi_list = edi_list
        if self.edi_list is None:
            return

        mt_dict = {}
        if self.edi_list is None:
            raise ModEMError(
                "mt_list is None, please input a list of "
                "mt files containing the full path"
            )

        if len(self.edi_list) == 0:
            raise ModEMError(
                "edi_list is empty, please input a list of "
                ".edi files containing the full path"
            )

        for edi in self.edi_list:
            mt_obj = mt.MT(edi)
            if mt_obj.station is None:
                continue
            # should be more efficient for appending to a dictionary
            mt_dict.update({mt_obj.station: mt_obj})

        return mt_dict

    def get_relative_station_locations(self, mt_dict, data_array):
        """
        Compute the relative station locations on a grid where the center is (0, 0) 
        
        Computes from station locations in mt_dict. 
        Calls modem.Station().get_station_locations()  
        
        If Data._center_lat, Data._center_lon are assigned the center will be 
        relative to that point.  
        
        :param mt_dict: dictionary of :class:`mtpy.core.mt.MT` objects, keys are 
        station names.
        :type mt_dict: dictionary
        :param data_array: data array 
        :type data_array: np.ndarray 
        :return: data_array with relative locations in keys labels 
        rel_east, rel_north, rel_elev
        :rtype: np.ndarray 
        
        .. seealso:: `mtpy.modeling.modem.station.Stations`
        
        """
        stations_obj = Stations(
            model_epsg=self.model_epsg, model_utm_zone=self.model_utm_zone
        )
        mt_list = [mt_dict[s_key] for s_key in sorted(mt_dict.keys())]
        stations_obj.get_station_locations(mt_list)

        # rotate locations if needed
        if self._rotation_angle != 0:
            # rotate stations the opposite way to the data
            stations_obj.rotate_stations(-self._rotation_angle)

        # fill data array
        data_array[:]["station"] = stations_obj.station
        data_array[:]["lat"] = stations_obj.lat
        data_array[:]["lon"] = stations_obj.lon
        data_array[:]["east"] = stations_obj.east
        data_array[:]["north"] = stations_obj.north
        data_array[:]["elev"] = stations_obj.elev
        data_array[:]["rel_east"] = stations_obj.rel_east
        data_array[:]["rel_north"] = stations_obj.rel_north
        data_array[:]["rel_elev"] = stations_obj.rel_elev
        data_array[:]["zone"] = stations_obj.utm_zone

        return data_array

    def get_data_periods(self, mt_dict):
        """
         Get an array of unique periods from the data
         
        :param mt_dict: dictionary of :class:`mtpy.core.mt.MT` objects, keys are 
        station names.
        :type mt_dict: dictionary
        :return: array of unique periods from all stations provided
        :rtype: np.ndarray

        """
        data_period_list = []
        for s_key, mt_obj in mt_dict.items():
            data_period_list.extend(list(1.0 / mt_obj.Z.freq))

        return np.array(sorted(list(set(data_period_list)), reverse=False))

    def make_period_list(self, mt_dict):
        """
        Create an array of periods to invert for.
        
        If these parameters are not None, uses them to compute the period array
            - Data.period_min
            - Data.period_max
            - Data.max_num_periods
            
        otherwise the period max and period min is estimated from the data.
        
        :param mt_dict: dictionary of :class:`mtpy.core.mt.MT` objects, keys are 
        station names.
        :type mt_dict: dictionary
        :raises: :class:`mtpy.utils.exceptions.DataError` if a parameter is missing
        
        .. code-block::
            :linenos:
                
            >>> md = Data()
            >>> md.period_min = 0.01
            >>> md.period_max = 1000
            >>> md.max_num_periods = 23
            >>> inversion_periods = md.make_period_list(mt_dict)
            
        """
        
        if self.period_list is not None:
            self.logger.debug(
                "Inverting periods "
                + ", ".join([f"{pp:.5E}" for pp in self.period_list])
            )
            return self.period_list

        data_period_list = self.get_data_periods(mt_dict)

        if self.period_min is not None and self.period_max is None:
            raise DataError("Need to input period_max")
        if self.period_max is not None and self.period_min is None:
            raise DataError("Need to input period_min")
        if (
            self.period_min is not None
            and self.period_max is not None
            and self.max_num_periods is None
        ):
            raise DataError("Need to input number of periods to use")

        min_index = np.where(data_period_list >= self.period_min)[0][0]
        max_index = np.where(data_period_list <= self.period_max)[0][-1]

        pmin = np.log10(data_period_list[min_index])
        pmax = np.log10(data_period_list[max_index])
        period_list = np.logspace(pmin, pmax, num=self.max_num_periods)

        self.logger.debug(
            "Inverting periods " + ", ".join([f"{pp:.5E}" for pp in self.period_list])
        )

        if period_list is None:  # YG: is this possible?
            raise ModEMError(
                "Need to input period_min, period_max, "
                "max_num_periods or a period_list"
            )
        return period_list

    @property
    def rotation_angle(self):
        """ angle to rotated the data by """
        return self._rotation_angle

    @rotation_angle.setter
    def rotation_angle(self, rotation_angle):
        """
        When the rotation angle is set rotate MT objects and fill data array with
        rotated values. 
        :param rotation_angle: angle 0 is N, 90 E, positive clockwise (degrees)
        :type rotation_angle: float
        
        """

        if self._rotation_angle == rotation_angle:
            return

        self.logger.debug(
            "Changing rotation angle from {0:.1f} to {1:.1f}".format(
                self._rotation_angle, rotation_angle
            )
        )

        self._rotation_angle = rotation_angle

        if self.mt_dict is None:
            self.logger.warning("mt_dict is None, rotation will not be applied to data")
            return

        for mt_key in sorted(self.mt_dict.keys()):
            mt_obj = self.mt_dict[mt_key]
            # check if data already rotated
            angle_to_rotate = self._rotation_angle - mt_obj.Z.rotation_angle
            mt_obj.Z.rotate(angle_to_rotate)
            mt_obj.Tipper.rotate(angle_to_rotate)

        self.logger.debug(
            "Data rotated to align with {0:.1f} deg clockwise from N".format(
                self._rotation_angle
            )
        )
        self.data_array = self.fill_data_array(self.mt_dict)

    def _initialise_empty_data_array(
        self,
        station_locations,
        period_list,
        location_type="LL",
        station_names=None,
        epsg=None,
        utm_zone=None,
    ):
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
        dtype = self.make_dtype((nf, 2, 2), (nf, 1, 2))
        self.data_array = np.zeros(len(station_locations), dtype=dtype)
        if location_type == "LL":
            self.data_array["lon"] = station_locations[:, 0]
            self.data_array["lat"] = station_locations[:, 1]
            for i in range(len(self.data_array["lon"])):
                lat, lon = self.data_array["lat"][i], self.data_array["lon"][i]
                east, north, zone = gis_tools.project_point_ll2utm(
                    lat, lon, epsg=epsg, utm_zone=utm_zone
                )
                self.data_array["east"][i] = east
                self.data_array["north"][i] = north
        else:
            self.data_array["east"] = station_locations[:, 0]
            self.data_array["north"] = station_locations[:, 1]
            for i in range(len(self.data_array["east"])):
                east, north = self.data_array["east"][i], self.data_array["north"][i]
                lat, lon = gis_tools.project_point_utm2ll(
                    east, north, utm_zone=utm_zone, epsg=epsg
                )
                self.data_array["lon"][i] = lon
                self.data_array["lat"][i] = lat

        # set non-zero values to array (as zeros will be deleted)
        # as we are setting up for forward modelling, actual values don't matter
        if self.inv_mode in "12":
            self.data_array["z"][:] = 10.0 + 10j
            self.data_array["z_err"][:] = 1e15
        if self.inv_mode == "1":
            self.data_array["tip"][:] = 0.1 + 0.1j
            self.data_array["tip_err"][:] = 1e15

        # set station names
        if station_names is not None:
            if len(station_names) != len(station_names):
                station_names = None

        if station_names is None:
            station_names = ["st%03i" % ss for ss in range(len(station_locations))]
        self.data_array["station"] = station_names

        # make an mt_dict
        self.mt_dict = {}
        for i, sname in enumerate(station_names):
            mtObj = mt.MT()
            mtObj.lat = self.data_array["lat"][i]
            mtObj.lon = self.data_array["lon"][i]

            mtObj.east = self.data_array["east"][i]
            mtObj.north = self.data_array["north"][i]
            mtObj.Z = mtz.Z(
                z_array=self.data_array["z"][i],
                z_err_array=self.data_array["z_err"][i],
                freq=1.0 / period_list,
            )
            mtObj.Tipper = mtz.Tipper(
                tipper_array=self.data_array["tip"][i],
                tipper_err_array=self.data_array["tip_err"][i],
                freq=1.0 / period_list,
            )
            mtObj.station = sname
            self.mt_dict[sname] = mtObj

        self.data_array = self.get_relative_station_locations(
            self.mt_dict, self.data_array
        )

    def fill_data_array(
        self, mt_dict, new_edi_dir=None, use_original_freq=False, longitude_format="LON"
    ):
        """
        Populate the data array from values in the :class:`mtpy.core.mt.MT` objects
        
        :param mt_dict: dictionary of :class:`mtpy.core.mt.MT` objects, keys are 
        station names.
        :type mt_dict: dictionary
        :param new_edi_dir: full path to a new folder to write EDI files with the
        inversion periods, defaults to None
        :type new_edi_dir: string or Path, optional
        :param use_original_freq: If True uses original frequencies in the data,
        defaults to False
        :type use_original_freq: Boolean, optional
        :param longitude_format: How to write the EDI file longitude, defaults to "LON"
        useful if you want to read into Winglink.
        :type longitude_format: string, optional
        :raises ValueError: If cannot compute locations
        :return: data array
        :rtype: np.ndarray
        
        .. code-block::
            :linenos:
            
            >>> from pathlib import Path
            >>> from mtpy.modeling.modem import Data
            >>> md = Data()
            >>> md.period_list = [.01, .1, 1, 10, 100]
            >>> edi_path = Path(r"/home/mt")
            >>> edi_list = list(edi_path.glob("*.edi"))
            >>> mt_dict = md.make_mt_dict(edi_list)
            >>> md.data_array = md.fill_data_array(mt_dict)

        """

        if self.period_list is None:
            self.period_list = self.make_period_list(mt_dict)

        ns = len(list(mt_dict.keys()))
        nf = len(self.period_list)

        # d_array = False
        # if self.data_array is not None:
        #     d_arr_copy = self.data_array.copy()
        #     d_array = True

        dtype = self.make_dtype((nf, 2, 2), (nf, 1, 2))
        data_array = np.zeros(ns, dtype=dtype)

        rel_distance = False
        for ii, s_key in enumerate(sorted(mt_dict.keys())):
            mt_obj = mt_dict[s_key]
            data_array[ii]["station"] = mt_obj.station
            data_array[ii]["lat"] = mt_obj.latitude
            data_array[ii]["lon"] = mt_obj.longitude
            data_array[ii]["east"] = mt_obj.east
            data_array[ii]["north"] = mt_obj.north
            data_array[ii]["elev"] = mt_obj.elevation
            data_array[ii]["zone"] = mt_obj.utm_zone
            try:
                data_array[ii]["rel_east"] = mt_obj.grid_east
                data_array[ii]["rel_north"] = mt_obj.grid_north
                data_array[ii]["rel_elev"] = mt_obj.grid_elev
                rel_distance = True
            except AttributeError:
                self.logger.debug(
                    "Unable to set relative locations from 'mt_obj' "
                    "- not filled yet."
                )
                pass

            # interpolate each station onto the period list
            # check bounds of period list
            interp_periods = self.period_list[
                np.where(
                    (self.period_list >= 1.0 / mt_obj.Z.freq.max())
                    & (self.period_list <= 1.0 / mt_obj.Z.freq.min())
                )
            ]

            # if specified, apply a buffer so that interpolation doesn't
            # stretch too far over periods
            if type(self.period_buffer) in [float, int]:
                interp_periods_new = []
                dperiods = 1.0 / mt_obj.Z.freq
                for iperiod in interp_periods:
                    # find nearest data period
                    difference = np.abs(iperiod - dperiods)
                    nearestdperiod = dperiods[difference == np.amin(difference)][0]
                    if (
                        max(nearestdperiod / iperiod, iperiod / nearestdperiod)
                        < self.period_buffer
                    ):
                        interp_periods_new.append(iperiod)

                interp_periods = np.array(interp_periods_new)

            # FZ: sort in order
            interp_periods = np.sort(interp_periods)
            self.logger.debug(
                "station_name and its original period: %s %s %s",
                mt_obj.station,
                len(mt_obj.Z.freq),
                1.0 / mt_obj.Z.freq,
            )
            self.logger.debug(
                "station_name and interpolation period: %s %s %s",
                mt_obj.station,
                len(interp_periods),
                interp_periods,
            )

            # default: use_original_freq = True, each MT station edi file will use it's own frequency-filtered.
            # no new freq in the output modem.dat file. select those freq of mt_obj according to interp_periods
            if use_original_freq:
                interp_periods = self.filter_periods(mt_obj, interp_periods)
                self.logger.debug(
                    "station_name and selected/filtered periods: %s, %s, %s",
                    mt_obj.station,
                    len(interp_periods),
                    interp_periods,
                )
                # in this case the below interpolate_impedance_tensor function will degenerate into a same-freq set.

            if len(interp_periods) > 0:  # not empty
                interp_z, interp_t = mt_obj.interpolate(
                    1.0 / interp_periods,
                    period_buffer=self.period_buffer,
                    bounds_error=False,
                )  # )
                # set rotation angle
                interp_z.rotation_angle = self.rotation_angle * np.ones(len(interp_z.z))
                interp_t.rotation_angle = self.rotation_angle * np.ones(
                    len(interp_t.tipper)
                )
                #                interp_z, interp_t = mt_obj.interpolate(1./interp_periods)
                for kk, ff in enumerate(interp_periods):
                    jj = np.where(self.period_list == ff)[0][0]
                    data_array[ii]["z"][jj] = interp_z.z[kk, :, :]
                    data_array[ii]["z_err"][jj] = interp_z.z_err[kk, :, :]
                    data_array[ii]["z_inv_err"][jj] = interp_z.z_err[kk, :, :]

                    if mt_obj.Tipper.tipper is not None:
                        data_array[ii]["tip"][jj] = interp_t.tipper[kk, :, :]
                        data_array[ii]["tip_err"][jj] = interp_t.tipper_err[kk, :, :]
                        data_array[ii]["tip_inv_err"][jj] = interp_t.tipper_err[
                            kk, :, :
                        ]

                # need to set the mt_object to have Z and T with same periods
                # as the data file, otherwise adding a station will not work.
                mt_obj.Z = interp_z
                mt_obj.Tipper = interp_t
                # FZ: try to output a new edi files. Compare with original edi?
                if new_edi_dir is not None and Path(new_edi_dir).is_dir():
                    # new_edifile = os.path.join(new_edi_dir, mt_obj.station + '.edi')

                    mt_obj.write_mt_file(
                        save_dir=new_edi_dir,
                        file_type="edi",
                        longitude_format=longitude_format,
                    )
            else:
                pass

        # BM: If we can't get relative locations from MT object,
        #  then get them from Station object
        if not rel_distance:
            try:
                data_array = self.get_relative_station_locations(mt_dict, data_array)
            except ValueError as error:
                if self.model_epsg is None and self.model_utm_zone is None:
                    msg = (
                        "Cannot compute relative locations without a "
                        + "model_epsg or model_utm_zone set."
                    )
                    self.logger.error(msg)
                    raise ValueError(msg)
                else:
                    self.logger.error(error)
                    raise ValueError(error)

        return data_array

    @staticmethod
    def filter_periods(mt_obj, per_array):
        """
        Select the periods of the mt_obj that are in per_array.
        used to do original freq inversion.

        :param mt_obj: MT object for single station
        :type mt_obj: :class:`mtpy.core.mt.MT`
        :param per_array: array of periods to map to
        :return: array of selected periods (subset) of the mt_obj
        """

        mt_per = 1.0 / mt_obj.Z.freq

        new_per = [
            p for p in mt_per if any([np.isclose(p, p2, 1.0e-8) for p2 in per_array])
        ]

        return np.array(new_per)

    @property
    def station_locations(self):
        """
        extract station locations from data array 
        
        :returns: :class:`mtpy.modeling.modem.station.Stations`
        """
        if self.data_array is None:
            return None

        station_locations = self.data_array[
            [
                "station",
                "lat",
                "lon",
                "north",
                "east",
                "elev",
                "rel_north",
                "rel_east",
                "rel_elev",
                "zone",
            ]
        ]
        input_dict = {
            "model_epsg": self.model_epsg,
            "model_utm_zone": self.model_utm_zone,
            "_center_lat": self._center_lat,
            "_center_lon": self._center_lon,
            "_center_elev": self._center_elev,
        }
        stations_obj = Stations(**input_dict)
        stations_obj.station_locations = station_locations

        return stations_obj

    @station_locations.setter
    def station_locations(self, station_locations):
        """
        take a station_locations array and populate data_array
        
        :param station_locations: array of station locations
        :type station_locations: :class:`mtpy.modeling.modem.Station`
        """
        if self.data_array is None:
            dtype = self.make_dtype(
                (len(self.period_list), 2, 2), (len(self.period_list), 1, 2)
            )
            self.data_array = np.zeros(
                station_locations.station_locations.size, dtype=dtype
            )
            for d_index, s_arr in enumerate(station_locations.station_locations):
                self.data_array[d_index]["lat"] = s_arr["lat"]
                self.data_array[d_index]["lon"] = s_arr["lon"]
                self.data_array[d_index]["east"] = s_arr["east"]
                self.data_array[d_index]["north"] = s_arr["north"]
                self.data_array[d_index]["elev"] = s_arr["elev"]
                self.data_array[d_index]["rel_east"] = s_arr["rel_east"]
                self.data_array[d_index]["rel_north"] = s_arr["rel_north"]
                self.data_array[d_index]["rel_elev"] = s_arr["rel_elev"]

        else:
            for s_arr in station_locations.station_locations:
                try:
                    d_index = np.where(self.data_array["station"] == s_arr["station"])[
                        0
                    ][0]
                except IndexError:
                    self.logger.warning(
                        "Could not find {0} in data_array".format(s_arr["station"])
                    )
                    d_index = None

                if d_index is not None:
                    self.data_array[d_index]["lat"] = s_arr["lat"]
                    self.data_array[d_index]["lon"] = s_arr["lon"]
                    self.data_array[d_index]["east"] = s_arr["east"]
                    self.data_array[d_index]["north"] = s_arr["north"]
                    self.data_array[d_index]["elev"] = s_arr["elev"]
                    self.data_array[d_index]["rel_east"] = s_arr["rel_east"]
                    self.data_array[d_index]["rel_north"] = s_arr["rel_north"]
                    self.data_array[d_index]["rel_elev"] = s_arr["rel_elev"]

    @property
    def center_point(self):
        """ center point derived from the data unless otherwise specified """

        if self.data_array is not None:
            return self.station_locations.center_point
        return None

    def compute_inv_error(self, data_array):
        """
        compute the error from the given parameters for a given data array
        
        :param data_array: data array to invert
        :type data_array: np.ndarray
        
        Uses parameters:
            - Data.error_type_z
            - Data.error_value_z
            - Data.error_type_tipper
            - Data.error_value_tiper
         
        **Impedance Error Types**
        
        =========== ==================================================================
        Error Type  Calculation                          
        =========== ================================================================== 
        egbert      error_value_z $\cdot \sqrt(|(Z_{xy}\cdot Z_{yx}|$)) 
        mean_od     error_value_z $\cdot (Z_{xy} + Z_{yx})/2$ 
        eigen       error_value_z $\cdot$ eigenvalues($Z(\omega)$) 
        median      error_value_z $\cdot$ median($Z(\omega)$)
        =========== ==================================================================
        
        **Tipper Error Types**
        
        =========== ==================================================================
        Error Type  Description                          
        =========== ================================================================== 
        abs         A value given to all tipper data 
        =========== ==================================================================

        .. note:: If floor is added to an error type then any value below that value 
        will be set to the floor and anything above the floor will remain above.  For
        example if the error floor is 5 but the measurement error is 7, the error 
        will be left at 7.  If the measurement error is 3 then it will be changed
        to 5.
        
        """
        # copy values over to inversion error
        data_array["z_inv_err"] = data_array["z_err"]
        data_array["tip_inv_err"] = data_array["tip_err"]

        # compute relative error for tipper
        if "floor" in self.error_type_tipper:
            t_index = np.where(data_array["tip_err"] < self.error_value_tipper)
            data_array["tip_inv_err"][t_index] = self.error_value_tipper
        elif "abs" in self.error_type_tipper:
            data_array["tip_inv_err"][:] = self.error_value_tipper
        else:
            raise DataError(
                "Unsupported error type (tipper): {}".format(self.error_type_tipper)
            )

        # consistency checks
        error_type_z_list = np.atleast_1d(self.error_type_z)

        if error_type_z_list.size != 1 and error_type_z_list.size != 4:
            raise DataError(
                "Either specify a single error_type_z for all components, or "
                "a 2x2 numpy array of error_type_z."
            )
        # end if

        # compute error for z
        err_value = self.error_value_z / 100.0
        for ss in range(data_array.shape[0]):
            for ff in range(
                max([data_array["z"].shape[1], data_array["tip"].shape[1]])
            ):
                d_xx = abs(data_array["z"][ss, ff, 0, 0])
                d_xy = abs(data_array["z"][ss, ff, 0, 1])
                d_yx = abs(data_array["z"][ss, ff, 1, 0])
                d_yy = abs(data_array["z"][ss, ff, 1, 1])
                d = np.array([d_xx, d_xy, d_yx, d_yy])
                nz = np.nonzero(d)

                if d.sum() == 0.0:  # YG: only works if all d >= 0
                    continue

                err = np.zeros(
                    (
                        error_type_z_list.size,
                        np.atleast_2d(err_value).shape[0],
                        np.atleast_2d(err_value).shape[1],
                    )
                )
                for ei, error_type_z in enumerate(error_type_z_list.flatten()):
                    if "egbert" in error_type_z:
                        # if both components masked, then take error floor from
                        # max of z_xx or z_yy
                        if d_xy == 0.0 and d_yx == 0.0:
                            err[ei] = err_value * np.max([d_xx, d_yy])
                        # else use the off diagonals depending on data availability
                        else:
                            if d_xy == 0.0:
                                d_xy = d_yx
                            if d_yx == 0.0:
                                d_yx = d_xy
                            err[ei] = err_value * np.sqrt(d_xy * d_yx)

                    elif "median" in error_type_z:
                        err[ei] = err_value * np.median(d[nz])

                    elif "mean_od" in error_type_z:
                        dod = np.array([d_xy, d_yx])
                        nzod = np.nonzero(dod)
                        err[ei] = err_value * np.mean(dod[nzod])

                    elif "eigen" in error_type_z:
                        d2d = d.reshape((2, 2))
                        err[ei] = err_value * np.abs(np.linalg.eigvals(d2d)).mean()
                        if np.atleast_1d(err[ei]).sum() == 0:
                            err[ei] = err_value * d[nz].mean()

                    elif "off_diagonals" in error_type_z:
                        # apply same error to xy and xx, and to yx and yy
                        # value is a % of xy and yx respectively
                        err[ei] = np.array([[d_xy, d_xy], [d_yx, d_yx]]) * err_value

                    elif "percent" in error_type_z:
                        # apply separate error floors to each component
                        err[ei] = err_value * np.abs(d[ei])
                    else:
                        raise DataError(
                            "error type (z) {0} not understood".format(error_type_z)
                        )
                # end for

                if error_type_z_list.size == 1:
                    data_array["z_inv_err"][ss, ff, :, :] = err[0]
                else:
                    for ei in np.arange(error_type_z_list.size):
                        ix, iy = np.divmod(ei, 2)
                        if err.shape[1] > 1:
                            data_array["z_inv_err"][ss, ff, ix, iy] = err[ei, ix, iy]
                        else:
                            data_array["z_inv_err"][ss, ff, ix, iy] = err[ei, 0, 0]

        # if there is an error floor
        if "floor" in self.error_type_z:
            f_index = np.where(data_array["z_inv_err"] < data_array["z_err"])
            data_array["z_inv_err"][f_index] = data_array["z_err"][f_index]

        return data_array

    def write_data_file(
        self,
        save_path=None,
        fn_basename=None,
        rotation_angle=None,
        compute_error=True,
        fill=True,
        elevation=False,
        use_original_freq=False,
        longitude_format="LON",
        new_edis=False,
    ):
        """
        
        :param save_path: full directory to save file to, defaults to None
        :type save_path: string or Path, optional
        :param fn_basename: Basename of the saved file, defaults to None
        :type fn_basename: string, optional
        :param rotation_angle: Angle to rotate the data to (positive clockwise, N=0),
        defaults to None
        :type rotation_angle: float, optional
        :param compute_error: If True recomputes error give parameters, defaults to True
        :type compute_error: Boolean, optional
        :param fill: If True recomputes the data array from given dictionary of 
        :class:`mtpy.core.mt.MT` objects, defaults to True
        :type fill: boolean, optional
        :param elevation: If True adds in elevation from 'rel_elev' column in data
        array, defaults to False
        :type elevation: boolean, optional
        :param use_original_freq: If True use original periods in 
        :class:`mtpy.core.mt.MT` objects, defaults to False
        :type use_original_freq: boolean, optional
        :param longitude_format: If new edis is True uses this format for the 
        longitude key, defaults to "LON"
        :type longitude_format: string, optional
        :param new_edis: if True writes new EDI files with the inversion data
        to a folder called save_path/new_edis, defaults to False
        :type new_edis: boolean, optional

        :raises NotImplementedError: If the inversion mode is not supported
        :raises DataError: :class:`mtpy.utils.exceptions.DataError` if a parameter
        is missing
        :return: full path to data file
        :rtype: Path

        .. code-block::
            :linenos:

            >>> from pathlib import Path
            >>> import mtpy.modeling.modem as modem
            >>> edi_path = Path(r"/home/mt/edi_files")
            >>> edi_list = list(edi_path.glob("*.ed"))
            >>> md = modem.Data(edi_list, period_min=.1, period_max=300,\
            >>> ...             max_num_periods=12)
            >>> md.write_data_file(save_path=r"/home/modem/inv1")
            /home/modem/inv1/ModemDataFile.dat
            
        """

        if save_path is not None:
            self.save_path = Path(save_path)
        if fn_basename is not None:
            self.data_fn = fn_basename

        self.data_fn = Path(self.save_path, self.data_fn)
        if self.mt_dict is None:
            self.mt_dict = self.make_mt_dict()

        self.period_list = self.make_period_list(self.mt_dict)

        # rotate data if desired
        if rotation_angle is not None:
            self.rotation_angle = rotation_angle

        # make edis along prescribe periods
        if new_edis:
            new_edi_dir = self.save_path.joinpath("new_edis")
            if not new_edi_dir.exists():
                new_edi_dir.mkdir()
        else:
            new_edi_dir = None

        # be sure to fill in data array
        if fill:
            self.data_array = self.fill_data_array(
                self.mt_dict,
                new_edi_dir=new_edi_dir,
                use_original_freq=use_original_freq,
                longitude_format=longitude_format,
            )

        if not elevation:
            self.data_array["rel_elev"][:] = 0.0

        d_lines = []
        for inv_mode in self.inv_mode_dict[self.inv_mode]:
            if "impedance" in inv_mode.lower():
                d_lines.append(
                    self.get_header_string(
                        self.error_type_z, self.error_value_z, self.rotation_angle
                    )
                )
                d_lines.append(self.header_string)
                d_lines.append("> {0}\n".format(inv_mode))
                d_lines.append("> exp({0}i\omega t)\n".format(self.wave_sign_impedance))
                d_lines.append("> {0}\n".format(self.units))

                n_sta = len(
                    np.nonzero(np.abs(self.data_array["z"]).sum(axis=(1, 2, 3)))[0]
                )
                n_per = len(
                    np.nonzero(np.abs(self.data_array["z"]).sum(axis=(0, 2, 3)))[0]
                )
            elif "vertical" in inv_mode.lower():
                d_lines.append(
                    self.get_header_string(
                        self.error_type_tipper,
                        self.error_value_tipper,
                        self.rotation_angle,
                    )
                )
                d_lines.append(self.header_string)
                d_lines.append("> {0}\n".format(inv_mode))
                d_lines.append("> exp({0}i\omega t)\n".format(self.wave_sign_tipper))
                d_lines.append("> []\n")
                n_sta = len(
                    np.nonzero(np.abs(self.data_array["tip"]).sum(axis=(1, 2, 3)))[0]
                )
                n_per = len(
                    np.nonzero(np.abs(self.data_array["tip"]).sum(axis=(0, 2, 3)))[0]
                )
            else:
                # maybe error here
                raise NotImplementedError(
                    "inv_mode {} is not supported yet".format(inv_mode)
                )

            d_lines.append(
                f"> {self.rotation_angle:.3g}\n"
            )  # orientation, need to add at some point
            if elevation:
                d_lines.append(
                    "> {0:>10.6f} {1:>10.6f} {2:>10.2f}\n".format(
                        self.center_point.lat[0],
                        self.center_point.lon[0],
                        self.center_point.elev[0],
                    )
                )
            else:
                d_lines.append(
                    "> {0:>10.6f} {1:>10.6f}\n".format(
                        self.center_point.lat[0], self.center_point.lon[0]
                    )
                )
            d_lines.append("> {0} {1}\n".format(n_per, n_sta))

            if compute_error:
                self.data_array = self.compute_inv_error(self.data_array)

            for ss in range(self.data_array["z"].shape[0]):
                for ff in range(self.data_array["z"].shape[1]):
                    for comp in self.inv_comp_dict[inv_mode]:
                        # index values for component with in the matrix
                        z_ii, z_jj = self.comp_index_dict[comp]

                        # get the correct key for data array according to comp
                        if comp.find("z") == 0:
                            c_key = "z"
                        elif comp.find("t") == 0:
                            c_key = "tip"

                        # get the value for that competent at that frequency
                        zz = self.data_array[ss][c_key][ff, z_ii, z_jj]
                        if (
                            zz.real != 0.0
                            and zz.imag != 0.0
                            and zz.real != 1e32
                            and zz.imag != 1e32
                        ):
                            if self.formatting == "1":
                                per = "{0:<12.5e}".format(self.period_list[ff])
                                sta = "{0:>7}".format(
                                    self.data_array[ss]["station"]
                                )  # .decode('UTF-8'))
                                lat = "{0:> 9.3f}".format(self.data_array[ss]["lat"])
                                lon = "{0:> 9.3f}".format(self.data_array[ss]["lon"])
                                eas = "{0:> 12.3f}".format(
                                    self.data_array[ss]["rel_east"]
                                )
                                nor = "{0:> 12.3f}".format(
                                    self.data_array[ss]["rel_north"]
                                )
                                ele = "{0:> 12.3f}".format(
                                    self.data_array[ss]["rel_elev"]
                                )
                                com = "{0:>4}".format(comp.upper())
                                if self.units.lower() == "ohm":
                                    rea = "{0:> 14.6e}".format(zz.real / 796.0)
                                    ima = "{0:> 14.6e}".format(zz.imag / 796.0)
                                elif self.units.lower() not in (
                                    "[v/m]/[t]",
                                    "[mv/km]/[nt]",
                                ):
                                    raise DataError(
                                        'Unsupported unit "{}"'.format(self.units)
                                    )
                                else:
                                    rea = "{0:> 14.6e}".format(zz.real)
                                    ima = "{0:> 14.6e}".format(zz.imag)

                            elif self.formatting == "2":
                                per = "{0:<14.6e}".format(self.period_list[ff])
                                sta = "{0:<10}".format(self.data_array[ss]["station"])
                                lat = "{0:> 14.6f}".format(self.data_array[ss]["lat"])
                                lon = "{0:> 14.6f}".format(self.data_array[ss]["lon"])
                                eas = "{0:> 12.3f}".format(
                                    self.data_array[ss]["rel_east"]
                                )
                                nor = "{0:> 15.3f}".format(
                                    self.data_array[ss]["rel_north"]
                                )
                                ele = "{0:> 10.3f}".format(
                                    self.data_array[ss]["rel_elev"]
                                )
                                com = "{0:>12}".format(comp.upper())
                                if self.units.lower() == "ohm":
                                    rea = "{0:> 17.6e}".format(zz.real / 796.0)
                                    ima = "{0:> 17.6e}".format(zz.imag / 796.0)
                                elif self.units.lower() not in (
                                    "[v/m]/[t]",
                                    "[mv/km]/[nt]",
                                ):
                                    raise DataError(
                                        'Unsupported unit "{}"'.format(self.units)
                                    )
                                else:
                                    rea = "{0:> 17.6e}".format(zz.real)
                                    ima = "{0:> 17.6e}".format(zz.imag)
                            else:
                                raise NotImplementedError(
                                    "format {}({}) is not supported".format(
                                        self.formatting, type(self.formatting)
                                    )
                                )

                            # get error from inversion error
                            abs_err = self.data_array["{0}_inv_err".format(c_key)][
                                ss, ff, z_ii, z_jj
                            ]

                            if np.isinf(abs_err) or np.isnan(abs_err):
                                abs_err = 10 ** (
                                    np.floor(
                                        np.log10(abs(max([float(rea), float(ima)])))
                                    )
                                )
                            abs_err = "{0:> 14.6e}".format(abs(abs_err))
                            # make sure that x==north, y==east, z==+down
                            dline = "".join(
                                [
                                    per,
                                    sta,
                                    lat,
                                    lon,
                                    nor,
                                    eas,
                                    ele,
                                    com,
                                    rea,
                                    ima,
                                    abs_err,
                                    "\n",
                                ]
                            )

                            d_lines.append(dline)
        with open(self.data_fn, "w") as dfid:
            dfid.writelines(d_lines)

        self.logger.info("Wrote ModEM data file to {0}".format(self.data_fn))
        return self.data_fn

    def convert_ws3dinv_data_file(
        self,
        ws_data_fn,
        station_fn=None,
        save_path=None,
        fn_basename="ws_data_file.dat",
    ):
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

        if not Path(ws_data_fn).is_file():
            raise ws.WSInputError("Did not find {0}, check path".format(ws_data_fn))

        if save_path is not None:
            save_path = Path(save_path)
        else:
            save_path = self.save_path

        # --> get data from data file
        wsd = ws.WSData()
        wsd.read_data_file(ws_data_fn, station_fn=station_fn)

        ns = wsd.data["station"].shape[0]
        nf = wsd.period_list.shape[0]

        self.period_list = wsd.period_list.copy()
        self.data_array = np.zeros(ns, dtype=self.make_dtype((nf, 2, 2), (nf, 1, 2)))

        # --> fill data array
        for ii, d_arr in enumerate(wsd.data):
            self.data_array[ii]["station"] = d_arr["station"]
            self.data_array[ii]["rel_east"] = d_arr["east"]
            self.data_array[ii]["rel_north"] = d_arr["north"]
            self.data_array[ii]["z"][:] = d_arr["z_data"]
            self.data_array[ii]["z_err"][:] = (
                d_arr["z_data_err"].real * d_arr["z_err_map"].real
            )
            self.data_array[ii]["station"] = d_arr["station"]
            self.data_array[ii]["lat"] = 0.0
            self.data_array[ii]["lon"] = 0.0
            self.data_array[ii]["rel_east"] = d_arr["east"]
            self.data_array[ii]["rel_north"] = d_arr["north"]
            self.data_array[ii]["elev"] = 0.0

        # need to change the inversion mode to be the same as the ws_data file
        if self.data_array["z"].all() == 0.0:
            if self.data_array["tip"].all() == 0.0:
                self.inv_mode = "4"
            else:
                self.inv_mode = "3"
        else:
            if self.data_array["tip"].all() == 0.0:
                self.inv_mode = "2"
            else:
                self.inv_mode = "1"

        # -->write file
        self.write_data_file()

    def convert_modem_to_ws(
        self, data_fn=None, ws_data_fn=None, error_map=[1, 1, 1, 1]
    ):
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
            save_path = self.save_path
            ws_data_fn = Path(save_path, "WS_Data.dat")

        else:
            save_path = Path(ws_data_fn).parent

        station_info = ws.WSStation()
        station_info.east = self.data_array["rel_east"]
        station_info.north = self.data_array["rel_north"]
        station_info.names = self.data_array["station"]
        station_info.elev = self.data_array["elev"]
        station_info.save_path = save_path
        station_info.write_station_file()

        ws_data = ws.WSData()
        ws_data.period_list = self.period_list.copy()
        ws_data.z_err_map = error_map
        ws_data.z_err = "data"
        z_shape = (self.period_list.size, 2, 2)
        data_dtype = [
            ("station", "|U50"),
            ("east", np.float),
            ("north", np.float),
            ("z_data", (np.complex, z_shape)),
            ("z_data_err", (np.complex, z_shape)),
            ("z_err_map", (np.complex, z_shape)),
        ]
        ws_data.data = np.zeros(self.data_array["station"].size, dtype=data_dtype)
        ws_data.data["station"][:] = self.data_array["station"]
        ws_data.data["east"] = self.data_array["rel_east"]
        ws_data.data["north"] = self.data_array["rel_north"]
        ws_data.data["z_data"][:, :, :] = self.data_array["z"]
        ws_data.data["z_data_err"][:, :, :] = self.data_array["z_err"] * (1 + 1j)
        ws_data.data["z_err_map"][:, :, :] = np.array([[1, 1], [1, 1]])

        ws_data.write_data_file(save_path=save_path, data_fn=ws_data_fn)

        return ws_data.data_fn, station_info.station_fn

    def read_data_file(self, data_fn, center_utm=None):
        """
        
        :param data_fn: full path to data file name
        :type data_fn: string or Path
        :param center_utm: option to provide real world coordinates of the center of
        the grid for putting the data and model back into utm/grid coordinates,
        format [east_0, north_0, z_0], defaults to None
        :type center_utm: list or tuple, optional
        :raises DataError: If cannot compute component

        Fills attributes:
            * data_array
            * period_list
            * mt_dict
            
        .. code-block::
                
            >>> md = Data()
            >>> md.read_data_file(r"/home/modem_data.dat")
            >>> md
            ModEM Data Object:
                Number of stations: 169
                Number of periods:  22
               	Period range:  
               		Min: 0.01 s
               		Max: 15230.2 s
               	Rotation angle:     0.0
               	Data center:        
               		 latitude:  39.6351 deg
               		 longitude: -119.8039 deg
               		 Elevation: 0.0 m
               		 Easting:   259368.9746 m
               		 Northing:  4391021.1981 m
               		 UTM zone:  11S
               	Model EPSG:         None
               	Model UTM zone:     None
               	Impedance data:     True
               	Tipper data:        True
            

        """

        self.data_fn = Path(data_fn)
        self.save_path = self.data_fn.parent
        self.fn_basename = self.data_fn.name

        if self.data_fn is None:
            raise DataError("data_fn is None, enter a data file to read.")
        elif not self.data_fn.is_file():
            raise DataError("Could not find {0}, check path".format(self.data_fn))

        with open(self.data_fn, "r") as dfid:
            dlines = dfid.readlines()

        header_list = []
        metadata_list = []
        data_list = []
        period_list = []
        station_list = []
        read_impedance = False
        read_tipper = False
        inv_list = []
        for dline in dlines:
            if dline.find("#") == 0:
                header_list.append(dline.strip())
            elif dline.find(">") == 0:
                # modem outputs only 7 characters for the lat and lon
                # if there is a negative they merge together, need to split
                # them up
                dline = dline.replace("-", " -")
                metadata_list.append(dline[1:].strip())
                if dline.lower().find("ohm") > 0:
                    self.units = "ohm"
                    continue
                elif dline.lower().find("mv") > 0:
                    self.units = "[mV/km]/[nT]"
                    continue
                elif dline.lower().find("vertical") > 0:
                    read_tipper = True
                    read_impedance = False
                    inv_list.append("Full_Vertical_Components")
                    continue
                elif dline.lower().find("impedance") > 0:
                    read_impedance = True
                    read_tipper = False
                    inv_list.append("Full_Impedance")
                    continue
                if dline.find("exp") > 0:
                    if read_impedance is True:
                        self.wave_sign_impedance = dline[dline.find("(") + 1]
                    elif read_tipper is True:
                        self.wave_sign_tipper = dline[dline.find("(") + 1]
                elif len(dline[1:].strip().split()) >= 2:
                    if dline.find(".") > 0:
                        value_list = [
                            float(value) for value in dline[1:].strip().split()
                        ]
                        if value_list[0] != 0.0:
                            self._center_lat = value_list[0]
                        if value_list[1] != 0.0:
                            self._center_lon = value_list[1]
                        try:
                            self._center_elev = value_list[2]
                        except IndexError:
                            self._center_elev = 0.0
                            self.logger.debug(
                                "Did not find center elevation in data file"
                            )
                elif len(dline[1:].strip().split()) < 2:
                    try:
                        self.rotation_angle = float(dline[1:].strip())
                    except ValueError:
                        continue

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
            if h_str.find("_deg") > 0:
                try:
                    self._rotation_angle = float(h_str[0 : h_str.find("_deg")])
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

        data_dict = {}
        z_dummy = np.zeros((len(self.period_list), 2, 2), dtype="complex")
        t_dummy = np.zeros((len(self.period_list), 1, 2), dtype="complex")

        index_dict = {
            "zxx": (0, 0),
            "zxy": (0, 1),
            "zyx": (1, 0),
            "zyy": (1, 1),
            "tx": (0, 0),
            "ty": (0, 1),
        }

        # dictionary for true false if station data (lat, lon, elev, etc)
        # has been filled already so we don't rewrite it each time
        tf_dict = {}
        for station in station_list:
            data_dict[station] = mt.MT()
            data_dict[station].Z = mtz.Z(
                z_array=z_dummy.copy(),
                z_err_array=z_dummy.copy().real,
                freq=1.0 / self.period_list,
            )
            data_dict[station].Tipper = mtz.Tipper(
                tipper_array=t_dummy.copy(),
                tipper_err_array=t_dummy.copy().real,
                freq=1.0 / self.period_list,
            )
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
                data_dict[dd[1]].latitude = dd[2]
                data_dict[dd[1]].longitude = dd[3]
                data_dict[dd[1]].grid_north = dd[4]
                data_dict[dd[1]].grid_east = dd[5]
                data_dict[dd[1]].grid_elev = dd[6]
                data_dict[dd[1]].elevation = dd[6]
                data_dict[dd[1]].station = dd[1]
                tf_dict[dd[1]] = True
            # fill in the impedance tensor with appropriate values
            if dd[7].find("Z") == 0:
                z_err = dd[10]
                if self.wave_sign_impedance == "+":
                    z_value = dd[8] + 1j * dd[9]
                elif self.wave_sign_impedance == "-":
                    z_value = dd[8] - 1j * dd[9]
                else:
                    raise DataError(
                        'Incorrect wave sign "{}" (impedance)'.format(
                            self.wave_sign_impedance
                        )
                    )

                if self.units.lower() == "ohm":
                    z_value *= 796.0
                    z_err *= 796.0
                elif self.units.lower() not in ("[v/m]/[t]", "[mv/km]/[nt]"):
                    raise DataError('Unsupported unit "{}"'.format(self.units))

                data_dict[dd[1]].Z.z[p_index, ii, jj] = z_value
                data_dict[dd[1]].Z.z_err[p_index, ii, jj] = z_err
            # fill in tipper with appropriate values
            elif dd[7].find("T") == 0:
                if self.wave_sign_tipper == "+":
                    data_dict[dd[1]].Tipper.tipper[p_index, ii, jj] = dd[8] + 1j * dd[9]
                elif self.wave_sign_tipper == "-":
                    data_dict[dd[1]].Tipper.tipper[p_index, ii, jj] = dd[8] - 1j * dd[9]
                else:
                    raise DataError(
                        'Incorrect wave sign "{}" (tipper)'.format(
                            self.wave_sign_tipper
                        )
                    )
                data_dict[dd[1]].Tipper.tipper_err[p_index, ii, jj] = dd[10]

        # make mt_dict an attribute for easier manipulation later
        self.mt_dict = data_dict

        ns = len(list(self.mt_dict.keys()))
        nf = len(self.period_list)
        self.data_array = np.zeros(ns, dtype=self.make_dtype((nf, 2, 2), (nf, 1, 2)))

        # Be sure to caclulate invariants and phase tensor for each station
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]

            # self.mt_dict[s_key].zinv.compute_invariants()
            self.mt_dict[s_key].pt.set_z_object(mt_obj.Z)
            self.mt_dict[s_key].Tipper.compute_amp_phase()
            self.mt_dict[s_key].Tipper.compute_mag_direction()

            self.data_array[ii]["station"] = mt_obj.station
            self.data_array[ii]["lat"] = mt_obj.latitude
            self.data_array[ii]["lon"] = mt_obj.longitude
            self.data_array[ii]["east"] = mt_obj.east
            self.data_array[ii]["north"] = mt_obj.north
            self.data_array[ii]["zone"] = mt_obj.utm_zone
            self.data_array[ii]["elev"] = mt_obj.elev
            self.data_array[ii]["rel_elev"] = mt_obj.grid_elev
            self.data_array[ii]["rel_east"] = mt_obj.grid_east
            self.data_array[ii]["rel_north"] = mt_obj.grid_north

            self.data_array[ii]["z"][:] = mt_obj.Z.z
            self.data_array[ii]["z_err"][:] = mt_obj.Z.z_err
            self.data_array[ii]["z_inv_err"][:] = mt_obj.Z.z_err

            self.data_array[ii]["tip"][:] = mt_obj.Tipper.tipper
            self.data_array[ii]["tip_err"][:] = mt_obj.Tipper.tipper_err
            self.data_array[ii]["tip_inv_err"][:] = mt_obj.Tipper.tipper_err

        # option to provide real world coordinates in eastings/northings
        # (ModEM data file contains real world center in lat/lon but projection
        # is not provided so utm is assumed, causing errors when points cross
        # utm zones. And lat/lon cut off to 3 d.p. causing errors in smaller
        # areas)
        if center_utm is not None:
            self.data_array["east"] = self.data_array["rel_east"] + center_utm[0]
            self.data_array["north"] = self.data_array["rel_north"] + center_utm[1]


    def write_vtk_station_file(
        self,
        vtk_save_path=None,
        vtk_fn_basename="ModEM_stations",
        geographic=False,
        shift_east=0,
        shift_north=0,
        shift_elev=0,
        units="km",
        coordinate_system="nez+",
    ):
        """
        
        :param vtk_save_path: directory to save vtk file to, defaults to None
        :type vtk_save_path: string or Path, optional
        :param vtk_fn_basename: filename basename of vtk file, note that .vtr 
        extension is automatically added, defaults to "ModEM_stations"
        :type vtk_fn_basename: string, optional
        :param geographic: If true puts the grid on geographic coordinates based 
        on the model_utm_zone, defaults to False
        :type geographic: boolean, optional
        :param shift_east: shift in east directions in meters, defaults to 0
        :type shift_east: float, optional
        :param shift_north: shift in north direction in meters, defaults to 0
        :type shift_north: float, optional
        :param shift_elev: shift in elevation + down in meters, defaults to 0
        :type shift_elev: float, optional
        :param units: Units of the spatial grid [ km | m | ft ], defaults to "km"
        :type units: string, optional
        :type : string
        :param coordinate_system: coordinate system for the station, either the
        normal MT right-hand coordinate system with z+ down or the sinister 
        z- down [ nez+ | enz- ], defaults to nez+
        :return: full path to VTK file
        :rtype: Path
        
        Write VTK file   
        >>> md.write_vtk_station_file(vtk_fn_basename="modem_stations")
        
        Write VTK file in geographic coordinates
        >>> md.write_vtk_station_file(vtk_fn_basename="modem_stations",
        >>> ...                       geographic=True)
        
        Write VTK file in geographic coordinates with z+ up
        >>> md.write_vtk_station_file(vtk_fn_basename="modem_stations",
        >>> ...                       geographic=True,
        >>> ...                       coordinate_system='enz-')

        """
        if isinstance(units, str):
            if units.lower() == "km":
                scale = 1.0 / 1000.00
            elif units.lower() == "m":
                scale = 1.0
            elif units.lower() == "ft":
                scale = 3.2808
        elif isinstance(units, (int, float)):
            scale = units

        if vtk_save_path is None:
            vtk_fn = self.save_path.joinpath(vtk_fn_basename)
        else:
            vtk_fn = Path(vtk_save_path, vtk_fn_basename)

        if not geographic:
            if coordinate_system == 'nez+':
                vtk_x = (self.station_locations.rel_north + shift_north) * scale
                vtk_y = (self.station_locations.rel_east + shift_east) * scale
                vtk_z = (self.station_locations.rel_elev + shift_elev) * scale
                extra = (self.station_locations.rel_elev + shift_elev) * scale
            elif coordinate_system == 'enz-':
                vtk_x = (self.station_locations.rel_north + shift_north) * scale
                vtk_y = (self.station_locations.rel_east + shift_east) * scale
                vtk_z = (self.station_locations.rel_elev + shift_elev) * scale
                extra = (self.station_locations.rel_elev + shift_elev) * scale
                
        else:
            self.station_locations.model_utm_zone = self.center_point.zone[0]
            if coordinate_system == 'nez+':
                vtk_y = (self.station_locations.north + shift_north) * scale
                vtk_x = (self.station_locations.east + shift_east) * scale
                vtk_z = -1 * (self.station_locations.elev + shift_elev) * scale
                extra = -1 * (self.station_locations.elev + shift_elev)
            elif coordinate_system == 'enz-':
                vtk_y = (self.station_locations.north + shift_north) * scale
                vtk_x = (self.station_locations.east + shift_east) * scale
                vtk_z = -1 * (self.station_locations.elev + shift_elev) * scale
                extra = -1 * (self.station_locations.elev + shift_elev)
        
        # write file
        pointsToVTK(vtk_fn.as_posix(), vtk_x, vtk_y, vtk_z, data={"elevation": extra})

        self.logger.info("Wrote station VTK file to {0}".format(vtk_fn))
        return vtk_fn

    def get_parameters(self):
        """
        get important parameters for documentation
        """

        parameter_list = [
            "error_type_z",
            "error_value_z",
            "error_type_tipper",
            "error_value_tipper",
            "wave_sign_impedance",
            "wave_sign_tipper",
            "rotation_angle",
            "save_path",
        ]

        parameter_dict = {}
        for parameter in parameter_list:
            key = "data.{0}".format(parameter)
            parameter_dict[key] = getattr(self, parameter)

        parameter_dict["data.period_min"] = self.period_list.min()
        parameter_dict["data.period_max"] = self.period_list.max()
        parameter_dict["data.period_num"] = self.period_list.size

        parameter_dict["data.inv_mode"] = self.inv_mode_dict[self.inv_mode]
        parameter_dict["data.num_stations"] = self.station_locations.station.size
        parameter_dict["data.center_point_ll"] = (
            self.center_point.lat[0],
            self.center_point.lon[0],
        )

        parameter_dict["data.center_point_utm"] = (
            self.center_point.north[0],
            self.center_point.east[0],
            self.center_point.zone[0],
        )
        return parameter_dict

    def center_stations(self, model_obj):
        """
        Center station locations to the middle of cells, is useful for
        topography cause it reduces edge effects of stations close to cell edges.
        Recalculates rel_east, rel_north to center of model cell.
        
        :param model_obj: :class:`mtpy.modeling.modem.Model` object of the model
        :type model_obj: :class:`mtpy.modeling.modem.Model`

        
        """
        
        for s_arr in self.station_locations.station_locations:
            e_index = np.where(model_obj.grid_east >= s_arr["rel_east"])[0][0] - 1
            n_index = np.where(model_obj.grid_north >= s_arr["rel_north"])[0][0] - 1

            mid_east = model_obj.grid_east[e_index : e_index + 2].mean()
            mid_north = model_obj.grid_north[n_index : n_index + 2].mean()

            s_index = np.where(self.data_array["station"] == s_arr["station"])[0][0]

            self.data_array[s_index]["rel_east"] = mid_east
            self.data_array[s_index]["rel_north"] = mid_north

    def project_stations_on_topography(self, model_object, air_resistivity=1e12,
                                       sea_resistivity=0.3, ocean_bottom=False):
        """
        Project stations on topography of a given model

        :param model_obj: :class:`mtpy.modeling.modem.Model` object of the model
        :type model_obj: :class:`mtpy.modeling.modem.Model`
        :param air_resistivity: resistivity value of air cells in the model
        :type air_resistivity:  float
        :param sea_resistivity: resistivity of sea
        :type sea_resistivity: float
        :param ocean_bottom: If True places stations at bottom of sea cells
        :type ocean_bottom: boolean
        
        Recaluclates rel_elev
        """

        # find index of each station on grid
        station_index_x = []
        station_index_y = []
        for sname in self.station_locations.station_locations["station"]:
            ss = np.where(self.station_locations.station_locations["station"] == sname)[
                0
            ][0]
            # relative locations of stations
            sx, sy = (
                self.station_locations.station_locations["rel_east"][ss],
                self.station_locations.station_locations["rel_north"][ss],
            )
            # indices of stations on model grid
            sxi = np.where(
                (sx <= model_object.grid_east[1:]) & (sx > model_object.grid_east[:-1])
            )[0][0]
            syi = np.where(
                (sy <= model_object.grid_north[1:])
                & (sy > model_object.grid_north[:-1])
            )[0][0]

            # first, check if there are any air cells
            if np.any(model_object.res_model[syi, sxi] > 0.95 * air_resistivity):
                szi = np.amin(
                    np.where(
                        (model_object.res_model[syi, sxi] < 0.95 * air_resistivity)
                    )[0]
                )
            # otherwise place station at the top of the model
            else:
                szi = 0
            
            # JP: estimate ocean bottom stations if requested
            if ocean_bottom:
                if np.any(model_object.res_model[syi, sxi] <= sea_resistivity):
                    szi = np.amax(
                        np.where(
                            (model_object.res_model[syi, sxi] <= sea_resistivity)
                        )[0]
                    )
                # if the stations are not in the ocean let the previous szi estimation
                # be used
            
            # get relevant grid point elevation
            topoval = model_object.grid_z[szi]

            station_index_x.append(sxi)
            station_index_y.append(syi)

            # update elevation in station locations and data array, +1 m as
            # data elevation needs to be below the topography (as advised by Naser)
            self.data_array["rel_elev"][ss] = topoval + 0.001

        # BM: After applying topography, center point of grid becomes
        #  highest point of surface model.
        self._center_elev = model_object.grid_z[0]

        self.logger.debug("Re-writing data file after adding topo to "
                          + self.data_fn.stem + "_topo.dat")
        self.write_data_file(
            fn_basename=self.data_fn.stem + "_topo.dat", fill=False, elevation=True,)

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
        self.logger.debug("ModEM data file number of sites:", num_sites)

        first_site_periods = md.data_array[0][9]  # (23L, 2L, 2L)
        self.logger.debug("first_site_periods = %s" % str(first_site_periods.shape[0]))

        period_list = md.period_list
        freq_list = 1.0 / period_list
        num_periods = len(period_list)
        self.logger.debug("ModEM data file number of periods:", num_periods)

        csv_basename = "modem_data_to_phase_tensor"
        csvfname = Path(dest_dir, f"{csv_basename}.csv")

        csv_header = [
            "Freq",
            "Station",
            "Lat",
            "Long",
            "Phimin",
            "Phimax",
            "Ellipticity",
            "Azimuth",
        ]

        with open(csvfname, "wb") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

        for period_num in range(num_periods):
            per = period_list[period_num]
            freq = freq_list[period_num]
            self.logger.debug("Working on period %s; frequency: %s", per, freq)
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
                arow = [
                    freq,
                    site_label,
                    site_lat,
                    site_long,
                    this_phimin,
                    this_phimax,
                    this_ellipticity,
                    this_azimuth,
                ]
                # Done for this site

                csvrows.append(arow)

            with open(
                csvfname, "ab"
            ) as csvf:  # append to this summary csv file for all freqs
                writer = csv.writer(csvf)
                writer.writerows(csvrows)

            csv_basename2 = f"{csv_basename}_{freq:.0f}Hz.csv"
            csvfile2 = Path(dest_dir, csv_basename2)

            with open(csvfile2, "wb") as csvf:  # csvfile  for eachindividual freq
                writer = csv.writer(csvf)
                writer.writerow(csv_header)
                writer.writerows(csvrows)

        # Done with all sites and periods
        self.logger.info("CSV files created in %s", outdir)

        return csvfname

    def add_station(self, fn=None, mt_object=None, new_edi_dir=None):
        """
        Add a station to an existing data object.
        
        :param fn: file name or list of files names to add, defaults to None
        :type fn: str, Path, list of str or Path, optional
        :param mt_object: MT objects or list of MT objects, defaults to None
        :type mt_object: :class:`mtpy.core.mt.MT` or list of :class:`mtpy.core.mt.MT, 
        optional
        :return: new data array 
        :rtype: np.ndarray
        :return: new_mt_dict
        :type: dictionary of :class:`mtpy.core.mt.MT` objects
        
        .. note:: As long as the added station(s) are within the existing station
        area the center point will remain the same, if they are not, then the center
        will change and you will need to adjust your mesh.
        
        Example
        
        >>> dfn = "/home/mt/test.dat"
        >>> d = modem.Data()
        >>> d.read_data_file(dfn)
        >>> d.data_array, d.mt_dict = d.add_station(fn=[r"/01.edi", r"/02.edi"])
        >>> d.write_data_file(fn_basename="test_added.dat", 
        >>> ...               fill=False,
        >>> ...               compute_error=False,
        >>> ...               elevation=True)

        """

        if fn is not None:
            if isinstance(fn, (list, tuple)):
                mt_object = [mt.MT(f) for f in fn]
            elif isinstance(fn, (str, Path)):
                mt_object = [mt.MT(fn)]
            else:
                msg = f"Do not understand input type {type(fn)}"
                self.logger.error(msg)
                raise ValueError(msg)

        if mt_object is None:
            return None
        elif isinstance(mt_object, mt.MT):
            mt_object = [mt_object]
        elif not isinstance(mt_object, (list, tuple)):
            msg = f"Do not understand input type {type(fn)}"
            self.logger.error(msg)
            raise ValueError(msg)

        add_mt_dict = {}
        new_mt_dict = deepcopy(self.mt_dict)
        for mt_obj in mt_object:
            if mt_obj.station in self.mt_dict.keys():
                self.logger.warning(
                    f"Station {mt_obj.station} already exists, skipping"
                )
                continue
            add_mt_dict[mt_obj.station] = mt_obj
            new_mt_dict[mt_obj.station] = mt_obj

        add_data_array = self.fill_data_array(add_mt_dict, new_edi_dir=new_edi_dir)
        add_data_array = self.compute_inv_error(add_data_array)

        new_data_array = np.append(self.data_array, add_data_array)
        # need to sort stations because that is how things are calculated
        # for station location
        new_data_array.sort(kind="station")
        new_data_array = self.get_relative_station_locations(
            new_mt_dict, new_data_array
        )
        return new_data_array, new_mt_dict

    def remove_station(self, name):
        """
        Remove a station or stations
        
        :param name: name(s) of station(s) to remove
        :type name: string or list of strings
        :return: new data array with stations removed
        :rtype: np.ndarray
        :return: new dictionary of :class:`mtpy.core.mt` objects with stations removed
        :rtype: dictionary
        
        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> d.data_array, d.mt_dict = d.remove_station(["mt666", "mt013"]

        """
        if isinstance(name, str):
            name = [name]
        new_data_array = self.data_array.copy()
        new_mt_dict = dict(self.mt_dict)

        for b_station in name:
            try:
                s_find = np.where(self.data_array["station"] == b_station)[0][0]
            except IndexError:
                msg = f"Could not find {b_station} in data file"
                self.logger.warn(msg)
                continue

            new_data_array[s_find] = 0
            new_mt_dict.pop(b_station)

        return new_data_array[np.where(new_data_array["station"] != "0")], new_mt_dict

    def add_error(self, station, comp=[], z_value=5, t_value=0.05, periods=None):
        """
        
        Add error to a station's components for given period range
        
        :param station: name of station(s) to add error to
        :type station: string or list of strings
        :param comp: list of components to add data to, valid components are 
        zxx, zxy, zyx, zyy, tx, ty
        :type comp: string or list of strings
        :param periods: the period range to add to, if None all periods, otherwise
        enter as a tuple as (minimum, maximum) period in seconds
        :type periods: tuple (minimum, maxmum)
        :return: data array with added errors
        :rtype: np.ndarray
        
        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> d.data_array = d.add_error("mt01", comp=["zxx", "zxy", "tx"], z_value=7, t_value=.05)

        """
        c_dict = {
            "zxx": (0, 0),
            "zxy": (0, 1),
            "zyx": (1, 0),
            "zyy": (1, 1),
            "tx": (0, 0),
            "ty": (0, 1),
        }
        if isinstance(station, str):
            station = [station]
        if isinstance(comp, str):
            comp = [comp]
        if periods is not None:
            if len(periods) != 2:
                msg = "Must enter a minimum and maximum period value"
                self.logger.error(msg)
                raise ValueError(msg)
            p_min = np.where(self.period_list >= min(periods))[0][0]
            p_max = np.where(self.period_list <= max(periods))[0][-1]
        else:
            p_min = 0
            p_max = len(self.period_list) - 1

        new_data_array = self.data_array.copy()
        new_mt_dict = deepcopy(self.mt_dict)
        for ss in station:
            try:
                s_find = np.where(self.data_array["station"] == ss)[0][0]
            except IndexError:
                msg = f"Could not find {ss} in data file"
                self.logger.warn(msg)
                continue

            for cc in comp:
                try:
                    ii, jj = c_dict[cc]
                except KeyError:
                    msg = f"Component {cc} is not a valid component, skipping"
                    self.logger.warning(msg)
                    continue
                if "z" in cc:
                    new_data_array[s_find]["z_err"][p_min:p_max, ii, jj] *= z_value
                    new_mt_dict[ss].Z.z_err[p_min:p_max, ii, jj] *= z_value
                elif "t" in cc:
                    new_data_array[s_find]["tip_err"][p_min:p_max, ii, jj] += t_value
                    new_mt_dict[ss].Tipper.tipper_err[p_min:p_max, ii, jj] += t_value

        return new_data_array, new_mt_dict

    def flip_phase(self, station, comp=[]):
        """
        Flip the phase of a station in case its plotting in the wrong quadrant
        
        :param station: name(s) of station to flip phase
        :type station: string or list of strings
        :param comp: components to flip, valid inputs are zx, zy, tx, ty, defaults to []
        :type comp: list, optional
        :return: new_data_array
        :rtype: np.ndarray
        
        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> d.data_array, d.mt_dict = d.flip_phase("mt01", comp=["zx", "tx"])

        """
        c_dict = {"zx": 0, "zy": 1, "tx": 0, "ty": 1}
        if isinstance(station, str):
            station = [station]
        if isinstance(comp, str):
            comp = [comp]
        new_data_array = self.data_array.copy()
        new_mt_dict = deepcopy(self.mt_dict)
        for ss in station:
            try:
                s_find = np.where(self.data_array["station"] == ss)[0][0]
            except IndexError:
                msg = f"Could not find {ss} in data file"
                self.logger.warn(msg)
                continue
            for cc in comp:
                try:
                    index = c_dict[cc]
                except KeyError:
                    msg = f"Component {cc} is not a valid component, skipping"
                    self.logger.warning(msg)
                    continue
                if "z" in cc:
                    new_data_array[s_find]["z"][:, index, :] *= -1
                    new_mt_obj = new_mt_dict[ss].copy()
                    new_mt_obj.Z.z[:, index, :] *= -1
                    new_mt_dict[ss] = new_mt_obj

                elif "t" in cc:
                    new_data_array[s_find]["tip"][:, 0, index] *= -1
                    new_mt_dict[ss].Tipper.tipper[:, 0, index] *= -1

        return new_data_array, new_mt_dict

    def remove_static_shift(self, station, ss_x=1, ss_y=1):
        """
        Remove static shift from impedance components
        
        The application is  Z  / ss_i 
        
        .. note:: The factors are in resistivity scale, so the
                  entries of  the matrix "S" need to be given by their
                  square-roots!
        
        :param name: name(s) of station(s) to apply static shift
        :type name: string or list of strings
        :param ss_x: static correction in electric x components, defaults to 1
        :type ss_x: float, optional
        :param ss_y: static correction in electric y components, defaults to 1
        :type ss_y: float, optional
        :return: new data array
        :rtype: np.ndarray
        :return: new dictionary of :class:`mtpy.core.mt` objects
        rtype: dictionary
        
        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> d.data_array, d.mt_dict = d.remove_static_shift("mt01", ss_x=1.5, ss_y=.4)

        """

        if isinstance(station, str):
            station = [station]

        new_data_array = self.data_array.copy()
        new_mt_dict = deepcopy(self.mt_dict)
        for ss in station:
            try:
                s_find = np.where(self.data_array["station"] == ss)[0][0]
            except IndexError:
                msg = f"Could not find {ss} in data file"
                self.logger.warn(msg)
                continue
            new_data_array[s_find]["z"][:, 0, :] *= 1.0 / np.sqrt(ss_x)
            new_data_array[s_find]["z"][:, 1, :] *= 1.0 / np.sqrt(ss_y)

            new_mt_obj = new_mt_dict[ss].copy()
            new_mt_obj.Z = new_mt_obj.remove_static_shift(ss_x=ss_x, ss_y=ss_y)
            new_mt_dict[ss] = new_mt_obj

        return new_data_array, new_mt_dict

    def remove_component(
        self, station, zxx=False, zxy=False, zyy=False, zyx=False, tx=False, ty=False
    ):
        """
        Remove a component for a given station(s)
        
        :param station: station name or list of station names
        :type station: string or list
        :param zxx: Z_xx, defaults to False
        :type zxx: TYPE, optional
        :param zxy: Z_xy, defaults to False
        :type zxy: TYPE, optional
        :param zyy: Z_yx, defaults to False
        :type zyy: TYPE, optional
        :param zyx: Z_yy, defaults to False
        :type zyx: TYPE, optional
        :param tx: T_zx, defaults to False
        :type tx: TYPE, optional
        :param ty: T_zy, defaults to False
        :type ty: TYPE, optional
        :return: new data array with components removed
        :rtype: np.ndarray 
        :return: new mt_dict with components removed
        :rtype: dictionary
        
        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> d.data_array, d.mt_dict = d.remove_component("mt01", zxx=True, tx=True)

        """
        c_dict = {
            "zxx": {"index": (0, 0), "bool": zxx},
            "zxy": {"index": (0, 1), "bool": zxy},
            "zyx": {"index": (1, 0), "bool": zyx},
            "zyy": {"index": (1, 1), "bool": zyy},
            "tx": {"index": (0, 0), "bool": tx},
            "ty": {"index": (0, 1), "bool": ty},
        }
        if not isinstance(station, (list, tuple)):
            station = [station]

        new_data_array = self.data_array.copy()
        new_mt_dict = deepcopy(self.mt_dict)

        for ss in station:
            try:
                s_find = np.where(self.data_array["station"] == ss)[0][0]
            except IndexError:
                msg = f"Could not find {ss} in data file"
                self.logger.warn(msg)
                continue
            for ckey, dd in c_dict.items():
                if dd["bool"]:
                    if "z" in ckey:
                        new_data_array[s_find]["z"][
                            :, dd["index"][0], dd["index"][1]
                        ] = 0
                        new_data_array[s_find]["z_err"][
                            :, dd["index"][0], dd["index"][1]
                        ] = 0
                        new_mt_dict[ss].Z.z[:, dd["index"][0], dd["index"][1]] = 0
                        new_mt_dict[ss].Z.z_err[:, dd["index"][0], dd["index"][1]] = 0
                    elif "t" in ckey:
                        new_data_array[s_find]["tip"][
                            :, dd["index"][0], dd["index"][1]
                        ] = 0
                        new_data_array[s_find]["tip_err"][
                            :, dd["index"][0], dd["index"][1]
                        ] = 0
                        new_mt_dict[ss].Tipper.tipper[
                            :, dd["index"][0], dd["index"][1]
                        ] = 0
                        new_mt_dict[ss].Tipper.tipper_err[
                            :, dd["index"][0], dd["index"][1]
                        ] = 0

        return new_data_array, new_mt_dict

    def estimate_starting_rho(self):
        """
        Estimate starting resistivity from the data.
        Creates a plot of the mean and median apparent resistivity values.
        
        :return: array of the median rho per period
        :rtype: np.ndarray(n_periods)
        :return: array of the mean rho per period
        :rtype: np.ndarray(n_periods)
        
        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> rho_median, rho_mean = d.estimate_starting_rho()
        
        """
        rho = np.zeros((self.data_array.shape[0], self.period_list.shape[0]))

        for ii, d_arr in enumerate(self.data_array):
            z_obj = mtz.Z(d_arr["z"], freq=1.0 / self.period_list)
            rho[ii, :] = z_obj.res_det

        mean_rho = np.apply_along_axis(lambda x: x[np.nonzero(x)].mean(), 0, rho)
        median_rho = np.apply_along_axis(lambda x: np.median(x[np.nonzero(x)]), 0, rho)

        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        (l1,) = ax.loglog(self.period_list, mean_rho, lw=2, color=(0.75, 0.25, 0))
        (l2,) = ax.loglog(self.period_list, median_rho, lw=2, color=(0, 0.25, 0.75))

        ax.loglog(
            self.period_list,
            np.repeat(mean_rho.mean(), self.period_list.size),
            ls="--",
            lw=2,
            color=(0.75, 0.25, 0),
        )
        ax.loglog(
            self.period_list,
            np.repeat(np.median(median_rho), self.period_list.size),
            ls="--",
            lw=2,
            color=(0, 0.25, 0.75),
        )

        ax.set_xlabel("Period (s)", fontdict={"size": 12, "weight": "bold"})
        ax.set_ylabel("Resistivity (Ohm-m)", fontdict={"size": 12, "weight": "bold"})

        ax.legend(
            [l1, l2],
            [
                "Mean = {0:.1f}".format(mean_rho.mean()),
                "Median = {0:.1f}".format(np.median(median_rho)),
            ],
            loc="upper left",
        )
        ax.grid(which="both", ls="--", color=(0.75, 0.75, 0.75))

        plt.show()
        
        return median_rho, mean_rho
