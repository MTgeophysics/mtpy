# -*- coding: utf-8 -*-
"""
.. module:: MT
   :synopsis: The main container for MT response functions.

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
"""

# ==============================================================================
from pathlib import Path

import numpy as np
from scipy import interpolate as spi

from mt_metadata.transfer_functions.core import TF

from mtpy.core.z import Z, Tipper
import mtpy.analysis.pt as MTpt
import mtpy.analysis.distortion as MTdistortion
from mtpy.utils import gis_tools
from mtpy.imaging import PlotMTResponse, PlotPhaseTensor


# =============================================================================
class MT(TF):
    """
    Basic MT container to hold all information necessary for a MT station
    including the following parameters.

    
    """

    def __init__(self, fn=None, **kwargs):
        super().__init__(fn=fn, **kwargs)

        self._Z = Z()
        self._Tipper = Tipper()
        self._rotation_angle = 0
        self._utm_location = {"east": 0, "north": 0, "zone": None}

        self.save_dir = Path.cwd()

        self.project_to_utm()

    def project_to_utm(self):
        """
        project point to utm
        
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.latitude and self.longitude:
            east, north, zone = gis_tools.project_point_ll2utm(
                self.latitude, self.longitude
            )
            self._utm_location = {"east": east, "north": north, "zone": zone}

    def project_to_ll(self):
        """
        project point to utm
        
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.east != 0 and self.north != 0 and self.utm_zone != None:
            self.latitude, self.longitude = gis_tools.project_point_utm2ll(
                self.east, self.north, self.utm_zone
            )

    @property
    def east(self):
        """ easting """
        return self._utm_location["east"]

    @east.setter
    def east(self, value):
        """ set east """
        self._utm_location["east"] = value

    @property
    def north(self):
        """ northing """
        return self._utm_location["north"]

    @north.setter
    def north(self, value):
        """ set north"""
        self._utm_location["north"] = value

    @property
    def utm_zone(self):
        """ utm zone """
        return self._utm_location["zone"]

    @utm_zone.setter
    def utm_zone(self, value):
        """ set utm_zone """
        self._utm_location["utm_zone"] = value

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
                freq=self.frequency,
            )
        return Z()

    @Z.setter
    def Z(self, z_object):
        """
        set z_object

        recalculate phase tensor and invariants, which shouldn't change except
        for strike angle
        """
        if not isinstance(z_object.freq, type(None)):
            if not (self.frequency == z_object.freq).all():
                self.frequency = z_object.freq
        self.impedance = z_object.z
        self.impedance_error = z_object.z_err

    @property
    def Tipper(self):
        """mtpy.core.z.Tipper object to hold tipper information"""

        if self.has_tipper():
            return Tipper(
                tipper_array=self.tipper.to_numpy(),
                tipper_err_array=self.tipper_error.to_numpy(),
                freq=self.frequency,
            )

    @Tipper.setter
    def Tipper(self, t_object):
        """
        set tipper object

        recalculate tipper angle and magnitude
        """

        if not isinstance(t_object.freq, type(None)):
            if not (self.frequency == t_object.freq).all():
                self.frequency = t_object.freq
        self.tipper = t_object.tipper
        self.tipper_error = t_object.tipper_err

    @property
    def pt(self):
        """mtpy.analysis.pt.PhaseTensor object to hold phase tensor"""
        return MTpt.PhaseTensor(z_object=self.Z)

    @property
    def ex_metadata(self):
        """ EX metadata """
        return self.station_metadata.runs[0].ex

    @ex_metadata.setter
    def ex_metadata(self, value):
        """ set EX metadata """
        self.station_metadata.runs[0].ex = value

    @property
    def ey_metadata(self):
        """ EY metadata """
        return self.station_metadata.runs[0].ey

    @ey_metadata.setter
    def ey_metadata(self, value):
        """ set EY metadata """
        self.station_metadata.runs[0].ey = value

    @property
    def hx_metadata(self):
        """ HX metadata """
        return self.station_metadata.runs[0].hx

    @hx_metadata.setter
    def hx_metadata(self, value):
        """ set hx metadata """
        self.station_metadata.runs[0].hx = value

    @property
    def hy_metadata(self):
        """ HY metadata """
        return self.station_metadata.runs[0].hy

    @hy_metadata.setter
    def hy_metadata(self, value):
        """ set hy metadata """
        self.station_metadata.runs[0].hy = value

    @property
    def hz_metadata(self):
        """ HZ metadata """
        return self.station_metadata.runs[0].hz

    @hz_metadata.setter
    def hz_metadata(self, value):
        """ set hz metadata """
        self.station_metadata.runs[0].hz = value

    @property
    def rrhx_metadata(self):
        """ RRHX metadata """
        return self.station_metadata.runs[0].rrhx

    @property
    def rrhy_metadata(self):
        """ RRHY metadata """
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
            z_array=new_z, z_err_array=self.Z.z_err.copy(), freq=self.Z.freq.copy()
        )

        return new_z_obj

    def interpolate(
        self,
        new_freq_array,
        interp_type="slinear",
        bounds_error=True,
        period_buffer=None,
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
        if period_buffer is not None:
            if 0.0 < period_buffer < 1.0:
                period_buffer += 1.0
                self.logger.warning(
                    "Period buffer must be > 1. Updating to", period_buffer
                )
        # check the bounds of the new frequency array
        if bounds_error:

            # logger.debug("new freq array %s", new_freq_array)
            if self.frequency.min() > new_freq_array.min():
                raise ValueError(
                    "New frequency minimum of {0:.5g}".format(new_freq_array.min())
                    + " is smaller than old frequency minimum of {0:.5g}".format(
                        self.frequency.min()
                    )
                    + ".  The new frequency range needs to be within the "
                    + "bounds of the old one."
                )
            if self.frequency.max() < new_freq_array.max():
                raise ValueError(
                    "New frequency maximum of {0:.5g}".format(new_freq_array.max())
                    + "is smaller than old frequency maximum of {0:.5g}".format(
                        self.frequency.max()
                    )
                    + ".  The new frequency range needs to be within the "
                    + "bounds of the old one."
                )
        # make a new Z object
        new_Z = Z(
            z_array=np.zeros((new_freq_array.shape[0], 2, 2), dtype="complex"),
            z_err_array=np.zeros((new_freq_array.shape[0], 2, 2)),
            freq=new_freq_array,
        )

        new_Tipper = Tipper(
            tipper_array=np.zeros((new_freq_array.shape[0], 1, 2), dtype="complex"),
            tipper_err_array=np.zeros((new_freq_array.shape[0], 1, 2)),
            freq=new_freq_array,
        )

        # interpolate the impedance tensor
        for ii in range(2):
            for jj in range(2):
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
                f = self.Z.freq[nz_index]

                # get frequencies to interpolate on to, making sure the
                # bounds are with in non-zero components
                new_nz_index = np.where(
                    (new_freq_array >= f.min()) & (new_freq_array <= f.max())
                )[0]
                new_f = new_freq_array[new_nz_index]

                # apply period buffer
                if type(period_buffer) in [float, int]:
                    new_f_update = []
                    new_nz_index_update = []
                    for ifidx, ifreq in enumerate(new_f):
                        # find nearest data period
                        difference = np.abs(np.log10(ifreq) - np.log10(f))
                        fidx = np.where(difference == np.amin(difference))[0][0]
                        if max(f[fidx] / ifreq, ifreq / f[fidx]) < period_buffer:
                            new_f_update.append(ifreq)
                            new_nz_index_update.append(new_nz_index[ifidx])
                    new_f = np.array(new_f_update)
                    new_nz_index = np.array(new_nz_index_update)
                # create a function that does 1d interpolation
                z_func_real = spi.interp1d(f, z_real, kind=interp_type)
                z_func_imag = spi.interp1d(f, z_imag, kind=interp_type)
                z_func_err = spi.interp1d(f, z_err, kind=interp_type)

                # interpolate onto new frequency range
                new_Z.z[new_nz_index, ii, jj] = z_func_real(new_f) + 1j * z_func_imag(
                    new_f
                )
                new_Z.z_err[new_nz_index, ii, jj] = z_func_err(new_f)
        # compute resistivity and phase for new Z object
        new_Z.compute_resistivity_phase()

        # if there is not tipper than skip
        if self.Tipper.tipper is None:
            return new_Z, new_Tipper
        # interpolate the Tipper
        for jj in range(2):
            # get indicies of non-zero components
            nz_index = np.nonzero(self.Tipper.tipper[:, 0, jj])

            if len(nz_index[0]) < 2:
                continue
            # get non-zero components
            t_real = self.Tipper.tipper[nz_index, 0, jj].real
            t_imag = self.Tipper.tipper[nz_index, 0, jj].imag
            t_err = self.Tipper.tipper_err[nz_index, 0, jj]

            # get frequencies for non-zero components
            f = self.Tipper.freq[nz_index]

            # create interpolation functions
            t_func_real = spi.interp1d(f, t_real, kind=interp_type)
            t_func_imag = spi.interp1d(f, t_imag, kind=interp_type)
            t_func_err = spi.interp1d(f, t_err, kind=interp_type)

            # get new frequency to interpolate over, making sure bounds are
            # for non-zero components
            new_nz_index = np.where(
                (new_freq_array >= f.min()) & (new_freq_array <= f.max())
            )
            new_f = new_freq_array[new_nz_index]

            # interpolate onto new frequency range
            new_Tipper.tipper[new_nz_index, 0, jj] = t_func_real(
                new_f
            ) + 1j * t_func_imag(new_f)

            new_Tipper.tipper_err[new_nz_index, 0, jj] = t_func_err(new_f)
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


# ==============================================================================
#             Error
# ==============================================================================


class MTError(Exception):
    pass
