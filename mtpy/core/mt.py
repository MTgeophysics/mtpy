# -*- coding: utf-8 -*-
"""
.. module:: MT
   :synopsis: The main container for MT response functions.

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
"""

# ==============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate as spi

from mt_metadata.transfer_functions.core import TF

from mtpy.core import Z, Tipper
from mtpy.core.mt_location import MTLocation
from mtpy.core.mt_dataframe import MTDataFrame

import mtpy.analysis.pt as MTpt
import mtpy.analysis.distortion as MTdistortion
from mtpy.imaging import PlotMTResponse, PlotPhaseTensor
from mtpy.modeling.errors import ModelErrors


# =============================================================================
class MT(TF, MTLocation):
    """
    Basic MT container to hold all information necessary for a MT station
    including the following parameters.


    """

    def __init__(self, fn=None, **kwargs):
        TF.__init__(self, fn=fn, **kwargs)
        MTLocation.__init__(self, **kwargs)

        self._Z = Z()
        self._Tipper = Tipper()
        self._rotation_angle = 0

        self.save_dir = Path.cwd()

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
        self._rotation_angle = theta_r
        z_copy = self.Z.copy()
        t_copy = self.Tipper.copy()

        z_copy.rotate(theta_r)
        self.Z = z_copy

        t_copy.rotate(theta_r)
        self.Tipper = t_copy

        self.logger.info(
            f"Rotated transfer function by: {self._rotation_angle:.3f} degrees clockwise"
        )

    @property
    def Z(self):
        """mtpy.core.z.Z object to hold impedance tensor"""

        if self.has_impedance():
            return Z(
                z_array=self.impedance.to_numpy(),
                z_err_array=self.impedance_error.to_numpy(),
                frequency=self.frequency,
                z_model_err_array=self.impedance_model_error.to_numpy(),
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
        self.impedance_error = z_object.z_err
        self.impedance_model_error = z_object.z_model_err

    @property
    def Tipper(self):
        """mtpy.core.z.Tipper object to hold tipper information"""

        if self.has_tipper():
            return Tipper(
                tipper_array=self.tipper.to_numpy(),
                tipper_err_array=self.tipper_error.to_numpy(),
                frequency=self.frequency,
                tipper_model_err_array=self.tipper_model_error.to_numpy(),
            )

    @Tipper.setter
    def Tipper(self, t_object):
        """
        set tipper object

        recalculate tipper angle and magnitude
        """

        if not isinstance(t_object.frequency, type(None)):
            if not (self.frequency == t_object.frequency).all():
                self.frequency = t_object.frequency
        self.tipper = t_object.tipper
        self.tipper_error = t_object.tipper_err
        self.tipper_model_error = t_object.tipper_model_err

    @property
    def pt(self):
        """mtpy.analysis.pt.PhaseTensor object to hold phase tensor"""
        return MTpt.PhaseTensor(z_object=self.Z)

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

    def remove_distortion(self, num_freq=None):
        """
        remove distortion following Bibby et al. [2005].

        :param num_freq: number of frequencies to look for distortion from the
                         highest frequency
        :type num_freq: int

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
        dummy_z_obj = Z.copy.deepcopy(self.Z)
        D, new_z_object = MTdistortion.remove_distortion(
            z_object=dummy_z_obj, num_freq=num_freq
        )

        return D, new_z_object

    def remove_static_shift(self, ss_x=1.0, ss_y=1.0):
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

        s_array, new_z = self.Z.remove_ss(
            reduce_res_factor_x=ss_x, reduce_res_factor_y=ss_y
        )

        new_z_obj = Z(
            z_array=new_z,
            z_err_array=self.Z.z_err.copy(),
            frequency=self.Z.frequency.copy(),
        )

        return new_z_obj

    def get_interp1d_functions_z(
        self, interp_type="slinear", bounds_error=False, fill_value=np.nan
    ):
        """

        :param interp_type: DESCRIPTION, defaults to "slinear"
        :type interp_type: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE

        """
        if self.Z is None:
            return None

        # interpolate the impedance tensor
        zmap = {0: "x", 1: "y"}
        interp_dict = {}
        for ii in range(2):
            for jj in range(2):
                comp = f"z{zmap[ii]}{zmap[jj]}"
                interp_dict[comp] = {}
                # need to look out for zeros in the impedance
                # get the indicies of non-zero components
                nz_index = np.nonzero(self.Z.z[:, ii, jj])

                if len(nz_index[0]) == 0:
                    continue
                # get the non-zero components
                z_real = self.Z.z[nz_index, ii, jj].real
                z_imag = self.Z.z[nz_index, ii, jj].imag
                z_err = self.Z.z_err[nz_index, ii, jj]

                # get the frequencies of non-zero components
                f = self.Z.frequency[nz_index]

                # create a function that does 1d interpolation
                interp_dict[comp]["real"] = spi.interp1d(
                    f,
                    z_real,
                    kind=interp_type,
                    bounds_error=bounds_error,
                )
                interp_dict[comp]["imag"] = spi.interp1d(
                    f,
                    z_imag,
                    kind=interp_type,
                    bounds_error=bounds_error,
                    fill_value=fill_value,
                )
                interp_dict[comp]["err"] = spi.interp1d(
                    f,
                    z_err,
                    kind=interp_type,
                    bounds_error=bounds_error,
                    fill_value=fill_value,
                )

        return interp_dict

    def get_interp1d_functions_t(
        self, interp_type="slinear", bounds_error=False, fill_value=np.nan
    ):
        """

        :param interp_type: DESCRIPTION, defaults to "slinear"
        :type interp_type: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE

        """
        if self.Tipper is None:
            return None

        # interpolate the impedance tensor
        zmap = {0: "x", 1: "y"}
        interp_dict = {}
        for jj in range(2):
            comp = f"tz{zmap[jj]}"
            interp_dict[comp] = {}
            # need to look out for zeros in the impedance
            # get the indicies of non-zero components
            nz_index = np.nonzero(self.Tipper.tipper[:, 0, jj])

            if len(nz_index[0]) == 0:
                continue
            # get the non-zero components
            t_real = self.Tipper.tipper[nz_index, 0, jj].real
            t_imag = self.Tipper.tipper[nz_index, 0, jj].imag
            t_err = self.Tipper.tipper_err[nz_index, 0, jj]

            # get the frequencies of non-zero components
            f = self.Tipper.frequency[nz_index]

            # create a function that does 1d interpolation
            interp_dict[comp]["real"] = spi.interp1d(
                f,
                t_real,
                kind=interp_type,
                bounds_error=bounds_error,
                fill_value=fill_value,
            )
            interp_dict[comp]["imag"] = spi.interp1d(
                f,
                t_imag,
                kind=interp_type,
                bounds_error=bounds_error,
                fill_value=fill_value,
            )
            interp_dict[comp]["err"] = spi.interp1d(
                f,
                t_err,
                kind=interp_type,
                bounds_error=bounds_error,
                fill_value=fill_value,
            )

        return interp_dict

    def interpolate(
        self,
        new_freq_array,
        interp_type="slinear",
        bounds_error=True,
    ):
        """
        Interpolate the impedance tensor onto different frequencies

        :param new_freq_array: a 1-d array of frequencies to interpolate on
                               to.  Must be with in the bounds of the existing
                               frequency range, anything outside and an error
                               will occur.
        :type new_freq_array: np.ndarray
        :param period_buffer: maximum ratio of a data period and the closest
                              interpolation period. Any points outside this
                              ratio will be excluded from the interpolated
                              impedance array.

        :returns: a new impedance object with the corresponding
                               frequencies and components.
        :rtype: mtpy.core.z.Z

        :returns: a new tipper object with the corresponding
                                    frequencies and components.
        :rtype: mtpy.core.z.Tipper

        :Interpolate: ::

            >>> import mtpy.core.mt as mt
            >>> edi_fn = r"/home/edi_files/mt_01.edi"
            >>> mt_obj = mt.MT(edi_fn)
            >>> # create a new frequency range to interpolate onto
            >>> new_freq = np.logspace(-3, 3, 24)
            >>> new_z_object, new_tipper_obj = mt_obj.interpolate(new_freq)
            >>> mt_obj.write_mt_file(new_fn=r"/home/edi_files/mt_01_interp.edi",
            >>> ...                   new_Z_obj=new_z_object,
            >>> ...                   new_Tipper_obj=new_tipper_object)

        """

        # make sure the input is a numpy array
        if not isinstance(new_freq_array, np.ndarray):
            new_freq_array = np.array(new_freq_array)

        # check the bounds of the new frequency array
        if bounds_error:

            if self.frequency.min() > new_freq_array.min():
                raise ValueError(
                    f"New frequency minimum of {new_freq_array.min():.5g} "
                    "is smaller than old frequency minimum of "
                    f"{self.frequency.min():.5g}.  The new frequency range "
                    "needs to be within the bounds of the old one."
                )
            if self.frequency.max() < new_freq_array.max():
                raise ValueError(
                    f"New frequency maximum of {new_freq_array.max():.5g} "
                    "is smaller than old frequency maximum of "
                    f"{self.frequency.max():.5g}.  The new frequency range "
                    "needs to be within the bounds of the old one."
                )
        # make a new Z object
        new_Z = Z(
            z_array=np.zeros((new_freq_array.shape[0], 2, 2), dtype="complex"),
            z_err_array=np.zeros((new_freq_array.shape[0], 2, 2)),
            frequency=new_freq_array,
        )

        new_Tipper = Tipper(
            tipper_array=np.zeros(
                (new_freq_array.shape[0], 1, 2), dtype="complex"
            ),
            tipper_err_array=np.zeros((new_freq_array.shape[0], 1, 2)),
            frequency=new_freq_array,
        )

        z_dict = self.get_interp1d_functions_z(interp_type=interp_type)
        t_dict = self.get_interp1d_functions_t(interp_type=interp_type)

        # interpolate the impedance tensor
        zmap = {0: "x", 1: "y"}
        for ii in range(2):
            for jj in range(2):
                comp = f"z{zmap[ii]}{zmap[jj]}"
                # get frequencies to interpolate on to, making sure the
                # bounds are with in non-zero components
                new_nz_index = np.where(
                    (new_freq_array >= self.Z.frequency.min())
                    & (new_freq_array <= self.Z.frequency.max())
                )[0]
                new_f = new_freq_array[new_nz_index]

                # interpolate onto new frequency range
                new_Z.z[new_nz_index, ii, jj] = z_dict[comp]["real"](
                    new_f
                ) + 1j * z_dict[comp]["imag"](new_f)
                new_Z.z_err[new_nz_index, ii, jj] = z_dict[comp]["err"](new_f)

        # if there is not tipper than skip
        if self.Tipper is None:
            return new_Z, new_Tipper

        # interpolate the Tipper
        for jj in range(2):
            comp = f"tz{zmap[jj]}"

            # get new frequency to interpolate over, making sure bounds are
            # for non-zero components
            new_nz_index = np.where(
                (new_freq_array >= self.Tipper.frequency.min())
                & (new_freq_array <= self.Tipper.frequency.max())
            )
            new_f = new_freq_array[new_nz_index]

            # interpolate onto new frequency range
            new_Tipper.tipper[new_nz_index, 0, jj] = t_dict[comp]["real"](
                new_f
            ) + 1j * t_dict[comp]["imag"](new_f)

            new_Tipper.tipper_err[new_nz_index, 0, jj] = t_dict[comp]["err"](
                new_f
            )
        new_Tipper.compute_mag_direction()

        return new_Z, new_Tipper

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
        mt_df = MTDataFrame()
        if cols is not None:
            mt_df.df_dtypes = mt_df._get_dtypes(cols)

        entry = mt_df.make_empty_entry(n_entries)

        entry["station"][:] = self.station
        entry["latitude"][:] = self.latitude
        entry["longitude"][:] = self.longitude
        entry["elevation"][:] = self.elevation
        entry["datum_epsg"][:] = self.datum_epsg
        entry["east"][:] = self.east
        entry["north"][:] = self.north
        entry["utm_epsg"][:] = self.utm_epsg
        entry["model_east"][:] = self.model_east
        entry["model_north"][:] = self.model_north
        entry["model_elevation"][:] = self.model_elevation

        entry["period"][:] = self.period
        if self.has_impedance():
            entry = mt_df.from_z_object(self.Z, entry)
        if self.has_tipper():
            entry = mt_df.from_t_object(self.Tipper, entry)

        return pd.DataFrame(entry)

    def from_dataframe(self, df):
        """
        fill transfer function attributes from a dataframe for a single station

        :param df: DESCRIPTION
        :type df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        mt_df = MTDataFrame()
        for key in [
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
        ]:
            try:
                setattr(self, key, df[key].unique()[0])
            except KeyError:
                continue

        self.Z = mt_df.to_z_object(df)
        self.Tipper = mt_df.to_t_object(df)

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
        off_diagonals              zxx_err == zxy_err, zyx_err == zyy_err
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
            array=self.impedance_error,
            error_value=error_value,
            error_type=error_type,
            floor=floor,
        )

        err = z_model_error.compute_error()

        if len(err.shape) == 1:
            z_err = np.zeros_like(self.impedance, dtype=float)
            z_err[:, 0, 0] = err
            z_err[:, 0, 1] = err
            z_err[:, 1, 0] = err
            z_err[:, 1, 1] = err

        else:
            z_err = err

        self.impedance_model_error = z_err

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
                "MT object containse no Tipper, cannot compute model errors"
            )
            return

        t_model_error = ModelErrors(
            array=self.tipper_error,
            error_value=error_value,
            error_type=error_type,
            floor=floor,
        )

        err = t_model_error.compute_error()

        if len(err.shape) == 1:
            t_err = np.zeros_like(self.tipper, dtype=float)
            t_err[:, 0, 0] = err
            t_err[:, 0, 1] = err

        else:
            t_err = err

        self.tipper_model_error = t_err

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

        z_model_err = self.impedance_model_error.copy().data
        t_model_err = self.tipper_model_error.copy().data
        for cc in comp:
            try:
                ii, jj = c_dict[cc]
            except KeyError:
                msg = f"Component {cc} is not a valid component, skipping"
                self.logger.warning(msg)
                continue
            if "z" in cc:
                z_model_err[p_min:p_max, ii, jj] *= z_value

            elif "t" in cc:
                t_model_err[p_min:p_max, ii, jj] += t_value

        self.impedance_model_error = z_model_err
        self.tipper_model_error = t_model_err

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
        :return: new_data_array
        :rtype: np.ndarray
        :return: new mt_dict with components removed
        :rtype: dictionary

        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> d.data_array, d.mt_dict = d.flip_phase("mt01", comp=["zx", "tx"])

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
                    z_obj.z_err[:, ii, jj] *= -1
                    z_obj.z_model_err[:, ii, jj] *= -1

                elif "t" in ckey:
                    t_change = True
                    t_obj.tipper[:, ii, jj] *= -1
                    t_obj.tipper_err[:, ii, jj] *= -1
                    t_obj.tipper_model_err[:, ii, jj] *= -1

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
        >>> d.data_array, d.mt_dict = d.remove_component("mt01", zxx=True, tx=True)

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
                    z_obj.z_err[:, ii, jj] = 0
                    z_obj.z_model_err[:, ii, jj] = 0

                elif "t" in ckey:
                    t_change = True
                    t_obj.tipper[:, ii, jj] = 0
                    t_obj.tipper_err[:, ii, jj] = 0
                    t_obj.tipper_model_err[:, ii, jj] = 0

        if z_change:
            self.Z = z_obj
        if t_change:
            self.Tipper = t_obj
