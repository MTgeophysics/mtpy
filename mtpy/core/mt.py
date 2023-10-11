# -*- coding: utf-8 -*-
"""
.. module:: MT
   :synopsis: The main container for MT response functions.

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
"""

# ==============================================================================
from pathlib import Path
from copy import deepcopy

import numpy as np

from mt_metadata.transfer_functions.core import TF

from mtpy.core import Z, Tipper
from mtpy.core.mt_location import MTLocation
from mtpy.core.mt_dataframe import MTDataFrame

from mtpy.imaging import (
    PlotMTResponse,
    PlotPhaseTensor,
    PlotPenetrationDepth1D,
)
from mtpy.modeling.errors import ModelErrors


# =============================================================================
class MT(TF, MTLocation):
    """
    Basic MT container to hold all information necessary for a MT station
    including the following parameters.


    """

    def __init__(self, fn=None, **kwargs):
        TF.__init__(self)
        MTLocation.__init__(self, survey_metadata=self._survey_metadata)

        # MTLocation.__init__(self)
        # TF.__init__(self)

        self.fn = fn

        self._Z = Z()
        self._Tipper = Tipper()
        self._rotation_angle = 0

        self.save_dir = Path.cwd()

        for key, value in kwargs.items():
            setattr(self, key, value)

    def clone_empty(self):
        """
        copy metadata but not the transfer function estimates
        """

        new_mt_obj = MT()
        new_mt_obj.survey_metadata.update(self.survey_metadata)
        new_mt_obj.station_metadata.update(self.station_metadata)
        new_mt_obj.station_metadata.runs = self.station_metadata.runs
        new_mt_obj._datum_crs = self._datum_crs
        new_mt_obj._utm_crs = self._utm_crs
        new_mt_obj._east = self._east
        new_mt_obj._north = self._north
        new_mt_obj.model_east = self.model_east
        new_mt_obj.model_north = self.model_north
        new_mt_obj.model_elevation = self.model_elevation
        new_mt_obj._rotation_angle = self._rotation_angle

        return new_mt_obj

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            if k in ["logger"]:
                continue

            setattr(result, k, deepcopy(v, memo))
        return result

    def copy(self):
        return deepcopy(self)

    @property
    def rotation_angle(self):
        """rotation angle in degrees from north"""
        return self._rotation_angle

    @rotation_angle.setter
    def rotation_angle(self, theta_r):
        """
        set rotation angle in degrees assuming North is 0 measuring clockwise
        positive to East as 90.

        upon setting rotates Z and Tipper

        TODO figure this out with xarray
        """

        self.rotate(theta_r)
        self._rotation_angle += theta_r

    def rotate(self, theta_r, inplace=True):
        """
        Rotate the data in degrees assuming North is 0 measuring clockwise
        positive to East as 90.

        :param theta_r: DESCRIPTION
        :type theta_r: TYPE
        :param inplace: DESCRIPTION, defaults to True
        :type inplace: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if inplace:
            if self.has_impedance():
                self.Z = self.Z.rotate(theta_r)
            if self.has_tipper():
                self.Tipper = self.Tipper.rotate(theta_r)

            self._rotation_angle = theta_r

            self.logger.info(
                f"Rotated transfer function by: {self._rotation_angle:.3f} degrees clockwise"
            )
        else:
            new_m = self.clone_empty()
            if self.has_impedance():
                new_m.Z = self.Z.rotate(theta_r)
            if self.has_tipper():
                new_m.Tipper = self.Tipper.rotate(theta_r)
            new_m._rotation_angle = self._rotation_angle
            return new_m

    @property
    def Z(self):
        """mtpy.core.z.Z object to hold impedance tensor"""

        if self.has_impedance():
            return Z(
                z=self.impedance.to_numpy(),
                z_error=self.impedance_error.to_numpy(),
                frequency=self.frequency,
                z_model_error=self.impedance_model_error.to_numpy(),
            )
        return Z()

    @Z.setter
    def Z(self, z_object):
        """
        set z_object

        recalculate phase tensor and invariants, which shouldn't change except
        for strike angle
        """
        if not isinstance(z_object.frequency, type(None)):
            if self.frequency.size != z_object.frequency.shape:

                self.frequency = z_object.frequency

            elif not (self.frequency == z_object.frequency).all():
                self.frequency = z_object.frequency
        self.impedance = z_object.z
        self.impedance_error = z_object.z_error
        self.impedance_model_error = z_object.z_model_error

    @property
    def Tipper(self):
        """mtpy.core.z.Tipper object to hold tipper information"""

        if self.has_tipper():
            return Tipper(
                tipper=self.tipper.to_numpy(),
                tipper_error=self.tipper_error.to_numpy(),
                frequency=self.frequency,
                tipper_model_error=self.tipper_model_error.to_numpy(),
            )

    @Tipper.setter
    def Tipper(self, t_object):
        """
        set tipper object

        recalculate tipper angle and magnitude
        """

        if t_object is None:
            return

        if not isinstance(t_object.frequency, type(None)):
            if not (self.frequency == t_object.frequency).all():
                self.frequency = t_object.frequency
        self.tipper = t_object.tipper
        self.tipper_error = t_object.tipper_error
        self.tipper_model_error = t_object.tipper_model_error

    @property
    def pt(self):
        """mtpy.analysis.pt.PhaseTensor object to hold phase tensor"""
        return self.Z.phase_tensor

    @property
    def ex_metadata(self):
        """EX metadata"""
        return self.station_metadata.runs[0].ex

    @ex_metadata.setter
    def ex_metadata(self, value):
        """set EX metadata"""
        self.station_metadata.runs[0].ex = value

    @property
    def ey_metadata(self):
        """EY metadata"""
        return self.station_metadata.runs[0].ey

    @ey_metadata.setter
    def ey_metadata(self, value):
        """set EY metadata"""
        self.station_metadata.runs[0].ey = value

    @property
    def hx_metadata(self):
        """HX metadata"""
        return self.station_metadata.runs[0].hx

    @hx_metadata.setter
    def hx_metadata(self, value):
        """set hx metadata"""
        self.station_metadata.runs[0].hx = value

    @property
    def hy_metadata(self):
        """HY metadata"""
        return self.station_metadata.runs[0].hy

    @hy_metadata.setter
    def hy_metadata(self, value):
        """set hy metadata"""
        self.station_metadata.runs[0].hy = value

    @property
    def hz_metadata(self):
        """HZ metadata"""
        return self.station_metadata.runs[0].hz

    @hz_metadata.setter
    def hz_metadata(self, value):
        """set hz metadata"""
        self.station_metadata.runs[0].hz = value

    @property
    def rrhx_metadata(self):
        """RRHX metadata"""
        return self.station_metadata.runs[0].rrhx

    @property
    def rrhy_metadata(self):
        """RRHY metadata"""
        return self.station_metadata.runs[0].rrhy

    def remove_distortion(
        self, n_frequencies=None, comp="det", only_2d=False, inplace=False
    ):
        """
        remove distortion following Bibby et al. [2005].

        :param n_frequencies: number of frequencies to look for distortion from the
                         highest frequency
        :type n_frequencies: int

        :returns: Distortion matrix
        :rtype: np.ndarray(2, 2, dtype=real)

        :returns: Z with distortion removed
        :rtype: mtpy.core.z.Z

        :Remove distortion and write new .edi file: ::

            >>> import mtpy.core.mt as mt
            >>> mt1 = mt.MT(fn=r"/home/mt/edi_files/mt01.edi")
            >>> D, new_z = mt1.remove_distortion()
            >>> mt1.write_mt_file(new_fn=r"/home/mt/edi_files/mt01_dr.edi",\
            >>>                    new_Z=new_z)

        """
        if inplace:
            self.Z = self.Z.remove_distortion(
                n_frequencies=n_frequencies,
                comp=comp,
                only_2d=only_2d,
                inplace=False,
            )
        else:
            new_mt = self.clone_empty()
            new_mt.Z = self.Z.remove_distortion(
                n_frequencies=n_frequencies,
                comp=comp,
                only_2d=only_2d,
                inplace=False,
            )
            new_mt.Tipper = self.Tipper
            return new_mt

    def remove_static_shift(self, ss_x=1.0, ss_y=1.0, inplace=False):
        """
        Remove static shift from the apparent resistivity

        Assume the original observed tensor Z is built by a static shift S
        and an unperturbated "correct" Z0 :

             * Z = S * Z0

        therefore the correct Z will be :
            * Z0 = S^(-1) * Z


        :param ss_x: correction factor for x component
        :type ss_x: float

        :param ss_y: correction factor for y component
        :type ss_y: float

        :returns: new Z object with static shift removed
        :rtype: mtpy.core.z.Z

        .. note:: The factors are in resistivity scale, so the
                  entries of  the matrix "S" need to be given by their
                  square-roots!

        :Remove Static Shift: ::

            >>> import mtpy.core.mt as mt
            >>> mt_obj = mt.MT(r"/home/mt/mt01.edi")
            >>> new_z_obj = mt.remove_static_shift(ss_x=.5, ss_y=1.2)
            >>> mt_obj.write_mt_file(new_fn=r"/home/mt/mt01_ss.edi",
            >>> ...                   new_Z_obj=new_z_obj)
        """

        if inplace:
            self.Z = self.Z.remove_ss(
                reduce_res_factor_x=ss_x,
                reduce_res_factor_y=ss_y,
                inplace=inplace,
            )

        else:
            new_mt = self.clone_empty()
            new_mt.Z = self.Z.remove_ss(
                reduce_res_factor_x=ss_x,
                reduce_res_factor_y=ss_y,
                inplace=inplace,
            )
            new_mt.Tipper = self.Tipper
            return new_mt

    def interpolate(
        self,
        new_period,
        method="slinear",
        bounds_error=True,
        f_type="period",
        **kwargs,
    ):
        """
        Interpolate the impedance tensor onto different frequencies

        :param new_period: a 1-d array of frequencies to interpolate on
         to.  Must be with in the bounds of the existing frequency range,
         anything outside and an error will occur.
        :type new_period: np.ndarray
        :param method: method to interpolate by, defaults to "cubic"
        :type method: string, optional
        :param bounds_error: check for if input frequencies are within the
         original frequencies, defaults to True
        :type bounds_error: boolean, optional
        :param f_type: frequency type can be [ 'frequency' | 'period' ]
        :type f_type: string, defaults to 'period'
        :param **kwargs: key word arguments for `interp`
        :type **kwargs: dictionary
        :raises ValueError: If input frequencies are out of bounds
        :return: New MT object with interpolated values.
        :rtype: :class:`mtpy.core.MT`


        .. note:: 'cubic' seems to work the best, the 'slinear' seems to do
         the same as 'linear' when using the `interp` in xarray.

        :Interpolate over frequency: ::

            >>> mt_obj = MT()
            >>> new_frequency = np.logspace(-3, 3, 20)
            >>> new_mt_obj = mt_obj.interpolate(new_frequency, f_type="frequency")

        """

        if f_type not in ["frequency", "freq", "period", "per"]:
            raise ValueError(
                "f_type must be either 'frequency' or 'period' not {f_type}"
            )

        # make sure the input is a numpy array
        if not isinstance(new_period, np.ndarray):
            new_period = np.array(new_period)

        if f_type in ["frequency", "freq"]:
            new_period = 1.0 / new_period

        # check the bounds of the new frequency array
        if bounds_error:

            if self.period.min() > new_period.min():
                raise ValueError(
                    f"New period minimum of {new_period.min():.5g} "
                    "is smaller than old period minimum of "
                    f"{self.period.min():.5g}.  The new period range "
                    "needs to be within the bounds of the old one."
                )
            if self.period.max() < new_period.max():
                raise ValueError(
                    f"New period maximum of {new_period.max():.5g} "
                    "is smaller than old frequency maximum of "
                    f"{self.period.max():.5g}.  The new period range "
                    "needs to be within the bounds of the old one."
                )

        new_m = self.clone_empty()
        if self.has_impedance():
            new_m.Z = self.Z.interpolate(new_period, method=method, **kwargs)
            if new_m.has_impedance():
                if np.all(np.isnan(new_m.Z.z)):

                    self.logger.warning(
                        f"Station {self.station}: Interpolated Z values are all NaN, "
                        "consider an alternative interpolation method. "
                        "See scipy.interpolate.interp1d for more information."
                    )
        if self.has_tipper():
            new_m.Tipper = self.Tipper.interpolate(
                new_period, method=method, **kwargs
            )
            if new_m.has_tipper():
                if np.all(np.isnan(new_m.Tipper.tipper)):
                    self.logger.warning(
                        f"Station {self.station}: Interpolated T values are all NaN, "
                        "consider an alternative interpolation method. "
                        "See scipy.interpolate.interp1d for more information."
                    )

        return new_m

    def plot_mt_response(self, **kwargs):
        """
        Returns a mtpy.imaging.plotresponse.PlotResponse object

        :Plot Response: ::

            >>> mt_obj = mt.MT(edi_file)
            >>> pr = mt.plot_mt_response()
            >>> # if you need more info on plot_mt_response
            >>> help(pr)

        """

        plot_obj = PlotMTResponse(
            z_object=self.Z,
            t_object=self.Tipper,
            pt_obj=self.pt,
            station=self.station,
            **kwargs,
        )

        return plot_obj

    def plot_phase_tensor(self, **kwargs):
        """

        :return: DESCRIPTION
        :rtype: TYPE

        """
        kwargs["ellipse_size"] = 0.5
        return PlotPhaseTensor(self.pt, station=self.station, **kwargs)

    def plot_depth_of_penetration(self, **kwargs):
        """
        Plot Depth of Penetration estimated from Niblett-Bostick estimation

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return PlotPenetrationDepth1D(self, **kwargs)

    def to_dataframe(self, utm_crs=None, cols=None):
        """
        Create a dataframe from the transfer function for use with plotting
        and modeling.

        :parameter utm_crs: the utm zone to project station to, could be a
         name, pyproj.CRS, EPSG number, or anything that pyproj.CRS can intake.
        :type utm_crs: string, int, :class:`pyproj.CRS`

        """
        if utm_crs is not None:
            self.utm_crs = utm_crs

        n_entries = self.period.size
        mt_df = MTDataFrame(n_entries=n_entries)

        mt_df.survey = self.survey
        mt_df.station = self.station
        mt_df.latitude = self.latitude
        mt_df.longitude = self.longitude
        mt_df.elevation = self.elevation
        mt_df.datum_epsg = self.datum_epsg
        mt_df.east = self.east
        mt_df.north = self.north
        mt_df.utm_epsg = self.utm_epsg
        mt_df.model_east = self.model_east
        mt_df.model_north = self.model_north
        mt_df.model_elevation = self.model_elevation
        mt_df.profile_offset = self.profile_offset

        mt_df.dataframe.loc[:, "period"] = self.period
        if self.has_impedance():
            mt_df.from_z_object(self.Z)
        if self.has_tipper():
            mt_df.from_t_object(self.Tipper)

        return mt_df

    def from_dataframe(self, mt_df):
        """
        fill transfer function attributes from a dataframe for a single station

        :param df: DESCRIPTION
        :type df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not isinstance(mt_df, MTDataFrame):
            try:
                mt_df = MTDataFrame(mt_df)
            except TypeError:
                raise TypeError(
                    f"Input dataframe must be an MTDataFrame not {type(mt_df)}"
                )
            except ValueError as error:
                raise ValueError(error)

        for key in [
            "survey",
            "station",
            "latitude",
            "longitude",
            "elevation",
            "east",
            "north",
            "utm_epsg",
            "model_north",
            "model_east",
            "model_elevation",
            "profile_offset",
        ]:
            try:
                setattr(self, key, getattr(mt_df, key))
            except KeyError:
                continue

        self.tf_id = self.station

        self.Z = mt_df.to_z_object()
        self.Tipper = mt_df.to_t_object()

    def compute_model_z_errors(
        self, error_value=5, error_type="geometric_mean", floor=True
    ):
        """
        Compute mode errors based on the error type

        ========================== ===========================================
        key                        definition
        ========================== ===========================================
        egbert                     error_value * sqrt(Zxy * Zyx)
        geometric_mean             error_value * sqrt(Zxy * Zyx)
        arithmetic_mean            error_value * (Zxy + Zyx) / 2
        mean_od                    error_value * (Zxy + Zyx) / 2
        off_diagonals              zxx_error == zxy_error, zyx_error == zyy_error
        median                     error_value * median(z)
        eigen                      error_value * mean(eigen(z))
        percent                    error_value * z
        absolute                   error_value
        ========================== ===========================================

        :param error_value: DESCRIPTION, defaults to 5
        :type error_value: TYPE, optional
        :param error_type: DESCRIPTION, defaults to "geometric_mean"
        :type error_type: TYPE, optional
        :param floor: DESCRIPTION, defaults to True
        :type floor: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not self.has_impedance():
            self.logger.warning(
                "MT Object contains no impedance data, cannot comput errors"
            )
            return

        z_model_error = ModelErrors(
            data=self.impedance,
            measurement_error=self.impedance_error,
            error_value=error_value,
            error_type=error_type,
            floor=floor,
            mode="impedance",
        )

        err = z_model_error.compute_error()

        if len(err.shape) == 1:
            z_error = np.zeros_like(self.impedance, dtype=float)
            z_error[:, 0, 0] = err
            z_error[:, 0, 1] = err
            z_error[:, 1, 0] = err
            z_error[:, 1, 1] = err

        else:
            z_error = err

        self.impedance_model_error = z_error

    def compute_model_t_errors(
        self, error_value=0.02, error_type="absolute", floor=False
    ):
        """
        Compute mode errors based on the error type

        ========================== ===========================================
        key                        definition
        ========================== ===========================================
        percent                    error_value * t
        absolute                   error_value
        ========================== ===========================================

        :param error_value: DESCRIPTION, defaults to .02
        :type error_value: TYPE, optional
        :param error_type: DESCRIPTION, defaults to "absolute"
        :type error_type: TYPE, optional
        :param floor: DESCRIPTION, defaults to True
        :type floor: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if not self.has_tipper():
            self.logger.warning(
                f"MT object for {self.station} contains no Tipper, cannot "
                "compute model errors"
            )
            return

        t_model_error = ModelErrors(
            data=self.tipper,
            measurement_error=self.tipper_error,
            error_value=error_value,
            error_type=error_type,
            floor=floor,
            mode="tipper",
        )

        err = t_model_error.compute_error()

        if len(err.shape) == 1:
            t_error = np.zeros_like(self.tipper, dtype=float)
            t_error[:, 0, 0] = err
            t_error[:, 0, 1] = err

        else:
            t_error = err

        self.tipper_model_error = t_error

    def add_model_error(self, comp=[], z_value=5, t_value=0.05, periods=None):
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
        >>> d.data = d.add_error("mt01", comp=["zxx", "zxy", "tx"], z_value=7, t_value=.05)

        """
        c_dict = {
            "zxx": (0, 0),
            "zxy": (0, 1),
            "zyx": (1, 0),
            "zyy": (1, 1),
            "tx": (0, 0),
            "ty": (0, 1),
        }

        if isinstance(comp, str):
            comp = [comp]
        if periods is not None:
            if len(periods) != 2:
                msg = "Must enter a minimum and maximum period value"
                self.logger.error(msg)
                raise ValueError(msg)
            p_min = np.where(self.period >= min(periods))[0][0]
            p_max = np.where(self.period <= max(periods))[0][-1]
        else:
            p_min = 0
            p_max = len(self.period) - 1

        z_model_error = self.impedance_model_error.copy().data
        t_model_error = self.tipper_model_error.copy().data
        for cc in comp:
            try:
                ii, jj = c_dict[cc]
            except KeyError:
                msg = f"Component {cc} is not a valid component, skipping"
                self.logger.warning(msg)
                continue
            if "z" in cc:
                z_model_error[p_min:p_max, ii, jj] *= z_value

            elif "t" in cc:
                t_model_error[p_min:p_max, ii, jj] += t_value

        self.impedance_model_error = z_model_error
        self.tipper_model_error = t_model_error

    def flip_phase(
        self, zxx=False, zxy=False, zyx=False, zyy=False, tzx=False, tzy=False
    ):
        """
        Flip the phase of a station in case its plotting in the wrong quadrant

        :param station: name(s) of station to flip phase
        :type station: string or list of strings
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
        :return: new_data
        :rtype: np.ndarray
        :return: new mt_dict with components removed
        :rtype: dictionary

        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> d.data, d.mt_dict = d.flip_phase("mt01", comp=["zx", "tx"])

        """
        c_dict = {
            "zxx": {"index": (0, 0), "bool": zxx},
            "zxy": {"index": (0, 1), "bool": zxy},
            "zyx": {"index": (1, 0), "bool": zyx},
            "zyy": {"index": (1, 1), "bool": zyy},
            "tzx": {"index": (0, 0), "bool": tzx},
            "tzy": {"index": (0, 1), "bool": tzy},
        }

        z_obj = self.Z.copy()
        t_obj = self.Tipper.copy()

        z_change = False
        t_change = False
        for ckey, dd in c_dict.items():
            if dd["bool"]:
                ii, jj = dd["index"]
                if "z" in ckey:
                    z_change = True
                    z_obj.z[:, ii, jj] *= -1
                    z_obj.z_error[:, ii, jj] *= -1
                    z_obj.z_model_error[:, ii, jj] *= -1

                elif "t" in ckey:
                    t_change = True
                    t_obj.tipper[:, ii, jj] *= -1
                    t_obj.tipper_error[:, ii, jj] *= -1
                    t_obj.tipper_model_error[:, ii, jj] *= -1

        if z_change:
            self.Z = z_obj
        if t_change:
            self.Tipper = t_obj

    def remove_component(
        self, zxx=False, zxy=False, zyy=False, zyx=False, tzx=False, tzy=False
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
        >>> d.data, d.mt_dict = d.remove_component("mt01", zxx=True, tx=True)

        """
        c_dict = {
            "zxx": {"index": (0, 0), "bool": zxx},
            "zxy": {"index": (0, 1), "bool": zxy},
            "zyx": {"index": (1, 0), "bool": zyx},
            "zyy": {"index": (1, 1), "bool": zyy},
            "tzx": {"index": (0, 0), "bool": tzx},
            "tzy": {"index": (0, 1), "bool": tzy},
        }

        z_obj = self.Z.copy()
        t_obj = self.Tipper.copy()

        z_change = False
        t_change = False
        for ckey, dd in c_dict.items():
            if dd["bool"]:
                ii, jj = dd["index"]
                if "z" in ckey:
                    z_change = True
                    z_obj.z[:, ii, jj] = 0
                    z_obj.z_error[:, ii, jj] = 0
                    z_obj.z_model_error[:, ii, jj] = 0

                elif "t" in ckey:
                    t_change = True
                    t_obj.tipper[:, ii, jj] = 0
                    t_obj.tipper_error[:, ii, jj] = 0
                    t_obj.tipper_model_error[:, ii, jj] = 0

        if z_change:
            self.Z = z_obj
        if t_change:
            self.Tipper = t_obj

    def add_white_noise(self, value, inplace=True):
        """
        Add white noise to the data, useful for synthetic tests.

        :param value: DESCRIPTION
        :type value: TYPE
        :param inplace: DESCRIPTION, defaults to True
        :type inplace: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if value > 1:
            value = value / 100.0

        if not inplace:
            new_mt_obj = self.clone_empty()

        tf_shape = self._transfer_function.transfer_function.shape
        noise_real = 1 + np.random.random(tf_shape) * value * (-1) ** (
            np.random.randint(0, 3, tf_shape)
        )
        noise_imag = 1 + np.random.random(tf_shape) * value * (-1) ** (
            np.random.randint(0, 3, tf_shape)
        )

        if inplace:
            self._transfer_function[
                "transfer_function"
            ] = self._transfer_function.transfer_function.real * (
                noise_real
            ) + (
                1j
                * self._transfer_function.transfer_function.imag
                * noise_imag
            )

            self._transfer_function["transfer_function_error"] = (
                self._transfer_function.transfer_function_error + value
            )

        else:
            new_mt_obj._transfer_function = self._transfer_function.copy()
            new_mt_obj._transfer_function[
                "transfer_function"
            ] = self._transfer_function.transfer_function.real * (
                noise_real
            ) + (
                1j
                * self._transfer_function.transfer_function.imag
                * noise_imag
            )

            self._transfer_function["transfer_function_error"] = (
                self._transfer_function.transfer_function_error + value
            )
            return new_mt_obj
