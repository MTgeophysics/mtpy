# -*- coding: utf-8 -*-
"""
===============
MT
===============

    * Basic MT-object containing the very basic and necessary information for a 
      single MT station

Created on Tue Jan 07 12:42:34 2014

@author: jpeacock-pr
"""

import os

import mtpy.utils.format as MTformat
import numpy as np

import mtpy.analysis.distortion as MTdistortion
import mtpy.analysis.pt as MTpt
import mtpy.analysis.zinvariants as MTinv

# ==============================================================================
import mtpy.core.edi as MTedi
import mtpy.core.z as MTz
import mtpy.imaging.plotresponse as plotresponse
import mtpy.utils.exceptions as MTex
import mtpy.utils.gis_tools

try:
    import scipy

    scipy_version = int(scipy.__version__.replace(".", ""))

    if scipy_version < 140:
        print (
            "Note: need scipy version 0.14.0 or higher or interpolation "
            + "might not work."
        )
    import scipy.interpolate as spi

    interp_import = True

except ImportError:
    print (
        "Could not find scipy.interpolate, cannot use method interpolate"
        + "check installation you can get scipy from scipy.org."
    )
    interp_import = False

# ==============================================================================


class MT(object):
    """
    Basic object containing all information necessary for a single MT station
    including
    
    Attributes
    ------------
    ===================== =====================================================
    Attribute             Description
    ===================== =====================================================
    name                  station name
    lat                   station latitude in decimal degrees
    lon                   station longitude in decimal degrees
    elev                  station elevation in meters
    Z                     mtpy.core.z.Z object for impedance tensor
    Tipper                mtpy.core.z.Tipper object for tipper
    date                  date collected
    east                  station location in UTM coordinates assuming WGS-84
    north                 station location in UTM coordinates assuming WGS-84 
    utm_zone              zone of UTM coordinates assuming WGS-84
    ===================== =====================================================
        
    .. note:: 
        * can change the utm grid by changing _utm_ellipsoid.  See 
          mtpy.utils.latlongutmconversion for details on reference 
          ellipsoids
        
        * currently the information is assumed to come from an edi file
          but can be extended later to .j files or something else or can
          be input by hand
              
        * if you set the coordinates east, north or utm_zone, be sure to 
          run _get_ll() to recalculate the latitude and longitude.
          
        * can input the following key words to fill values in Z and Tipper:
            - z_object        --> mtpy.core.z.Z object
            - z_array         --> np.ndarray(n_freq, 2, 2, dtype='complex')
            - z_err_array      --> np.ndarray(n_freq, 2, 2)
            - freq            --> np.ndarray(n_freq)
            - resistivity     --> np.ndarray(n_freq, 2, 2) (linear scale)
            - resistivity_err --> np.ndarray(n_freq, 2, 2) 
            - phase           --> np.ndarray(n_freq, 2, 2)           
            - phase_err       --> np.ndarray(n_freq, 2, 2) 
            - tipper_object   --> mtpy.core.z.Tipper object
            - tipper          --> np.ndarray(n_freq, 1, 2, dtype='complex') 
            - tipper_err       --> np.ndarray(n_freq, 1, 2)
        
    Methods
    ------------
    ===================== =====================================================
    Methods               Description
    ===================== =====================================================
    write_edi_file        write an edi_file from the MT data
    remove_distortion     remove distortion from the data following 
                          Bibby et al. [2005]
    remove_static_shift   Shifts apparent resistivity curves up or down
    interpolate           interpolates the impedance tensor and induction
                          vectors onto a specified frequency array.
    plot_mt_response      plots the MT response using mtpy.imaging.plotresponse
    ===================== =====================================================
    

    Examples
    -------------------
    :Read from an .edi File: ::
        
        >>> import mtpy.core.mt as mt
        >>> mt_obj = mt.MT(r"/home/edi_files/s01.edi")
    
    :Plot MT response: ::

        >>> # plot all components of mt response and phase tensor
        >>> plot_obj = mt_obj.plot_mt_response(plot_num=2, plot_pt='y')        
        >>> # plot the tipper as well
        >>> plot_obj.plot_tipper = 'yri'
        >>> plot_obj.redraw_plot()
      
    :Remove Distortion: ::
         
        >>> D, new_z = mt_obj.remove_distortion()
        >>> print D
        >>> np.array([[0.1, .9],
        >>> ...       [0.98, .43]])
        >>> # write a new edi file
        >>> mt_obj.write_edi_file(new_Z=new_z)
        >>> wrote file to: /home/edi_files/s01_RW.edi
    
    :Remove Static Shift: ::
         
        >>> new_z_obj = mt_obj.remove_static_shift(ss_x=.78, ss_y=1.1)
        >>> # write a new edi file
        >>> mt_obj.write_edi_file(new_fn=r"/home/edi_files/s01_ss.edi",
        >>>                       new_Z=new_z)
        >>> wrote file to: /home/edi_files/s01_ss.edi
        
    :Interpolate: ::
    
        >>> new_freq = np.logspace(-3, 3, num=24)
        >>> new_z_obj, new_tipper_obj = mt_obj.interpolate(new_freq)
        >>> mt_obj.write_edi_file(new_Z=new_z_obj, new_Tipper=new_tipper_obj)
        >>> wrote file to: /home/edi_files/s01_RW.edi
    """

    def __init__(self, fn=None, **kwargs):

        self._fn = fn
        self.station = kwargs.pop("station", None)
        self._lat = kwargs.pop("lat", None)
        self._lon = kwargs.pop("lon", None)
        self._elev = kwargs.pop("elev", None)
        self._Z = kwargs.pop("Z", MTz.Z())
        self._Tipper = kwargs.pop("Tipper", MTz.Tipper())
        self._utm_zone = kwargs.pop("utm_zone", None)
        self._east = kwargs.pop("east", None)
        self._north = kwargs.pop("north", None)
        self._rotation_angle = kwargs.pop("rotation_angle", 0)

        self.edi_object = MTedi.Edi()
        self.pt = None
        self.zinv = None
        self._utm_ellipsoid = 23

        # provide key words to fill values if an edi file does not exist
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

        # --> read in the file name given
        if self._fn is not None:
            self._set_fn(fn)

    # ==========================================================================
    # set functions
    # ==========================================================================
    def _set_lat(self, latitude):
        """
        set latitude making sure the input is in decimal degrees
        
        upon setting utm coordinates are recalculated
        """

        self._lat = MTformat._assert_position_format("lat", latitude)

        if self._lon is not None and self._lat is not None:
            self._get_utm()

    def _set_lon(self, longitude):
        """
        set longitude making sure the input is in decimal degrees
        
        upon setting utm coordinates are recalculated
        """

        self._lon = MTformat._assert_position_format("lon", longitude)

        if self._lon is not None and self._lat is not None:
            self._get_utm()

    def _set_elev(self, elevation):
        """
        set elevation, should be input as meters
        """

        self._elev = elevation

    def _set_east(self, easting):
        """
        set easting in meters
        
        upon setting lat and lon are recalculated
        """

        self._east = easting

    def _set_north(self, northing):
        """
        set northing in meters
        
        upon setting lat and lon are recalculated
        """

        self._north = northing

    def _set_utm_zone(self, utm_zone):
        """
        set UTM zone
        
        upon setting lat and lon are recalculated
        """

        self._utm_zone = utm_zone

    def _set_fn(self, filename):
        """
        set filename, currently only support .edi files
        """

        self._fn = filename
        if self._fn.lower().endswith(".edi"):
            self._read_edi_file()
        else:
            not_fn = self._fn[os.path.basename(self._fn).find["."] :]
            raise MTex.MTpyError_file_handling(
                "File " + "type {0} not supported yet.".format(not_fn)
            )

    def _set_rotation_angle(self, theta_r):
        """
        set rotation angle in degrees assuming North is 0 measuring clockwise
        positive to East as 90.
        
        upon setting rotates Z and Tipper
        """

        self._rotation_angle = theta_r
        self._Z.rotate(theta_r)
        self._Tipper.rotate(theta_r)
        self.pt.rotate(theta_r)
        self.zinv.rotate(theta_r)

        print (
            "Rotated Z, Tipper, Phase Tensor and Zinvariants by"
            "{0:.3f} degrees".format(self._rotation_angle)
        )

    def _set_Z(self, z_object):
        """
        set z_object
        
        recalculate phase tensor and invariants, which shouldn't change except
        for strike angle
        """

        self._Z = z_object
        self._Z._compute_res_phase()

        # --> compute phase tensor
        self.pt = MTpt.PhaseTensor(z_object=self._Z, freq=self._Z.freq)

        # --> compute invariants
        self.zinv = MTinv.Zinvariants(z_object=self._Z)

    def _set_Tipper(self, t_object):
        """
        set tipper object
        
        recalculate tipper angle and magnitude
        """

        self._Tipper = t_object
        if self._Tipper is not None:
            self._Tipper.compute_amp_phase()
            self._Tipper.compute_mag_direction()

    # ==========================================================================
    # get functions
    # ==========================================================================
    def _get_lat(self):
        return self._lat

    def _get_lon(self):
        return self._lon

    def _get_elev(self):
        return self._elev

    def _get_east(self):
        return self._east

    def _get_north(self):
        return self._north

    def _get_utm_zone(self):
        return self._utm_zone

    def _get_fn(self):
        return self._fn

    def _get_rotation_angle(self):
        return self._rotation_angle

    def _get_Z(self):
        return self._Z

    def _get_Tipper(self):
        return self._Tipper

    # ==========================================================================
    # set properties
    # ==========================================================================
    lat = property(_get_lat, _set_lat, doc="latitude of station in decimal degrees")

    lon = property(_get_lon, _set_lon, doc="longitude of station in decimal degrees")

    elev = property(_get_elev, _set_elev, doc="elevation in meters")

    east = property(
        _get_east, _set_east, doc="easting in meters of station location on UTM grid"
    )

    north = property(
        _get_north, _set_north, doc="northing in meters of station location on UTM grid"
    )

    utm_zone = property(_get_utm_zone, _set_utm_zone, doc="UTM zone")

    fn = property(_get_fn, _set_fn, doc="name of file containing MT info")

    rotation_angle = property(
        _get_rotation_angle, _set_rotation_angle, doc="rotation angle of Z and Tipper"
    )

    Z = property(_get_Z, _set_Z, doc="impedence tensor object")

    Tipper = property(_get_Tipper, _set_Tipper, doc="Tipper object")

    # --> conversion between utm and ll
    def _get_utm(self):
        """
        get utm coordinates from lat and lon
        """

        self.utm_zone, self.east, self.north = mtpy.utils.gis_tools.ll_to_utm(
            self._utm_ellipsoid, self.lat, self.lon
        )

    def _get_ll(self):
        """
        get lat and long from utm
        """

        self.lat, self.lon = mtpy.utils.gis_tools.utm_to_ll(
            self._utm_ellipsoid, self.north, self.east, self.utm_zone
        )

    # --> read in edi file
    def _read_edi_file(self):
        """
        read in edi file and set attributes accordingly
        
        """

        self.edi_object = MTedi.Edi(edi_fn=self.fn)
        self._lat = self.edi_object.lat
        self._lon = self.edi_object.lon
        self._elev = self.edi_object.elev
        self._Z = self.edi_object.Z
        self._Tipper = self.edi_object.Tipper
        self.station = self.edi_object.station

        # --> get utm coordinates from lat and lon
        self._get_utm()

        # --> make sure things are ordered from high frequency to low
        self._check_freq_order()

        # --> compute phase tensor
        self.pt = MTpt.PhaseTensor(z_object=self.Z, freq=self.Z.freq)

        # --> compute invariants
        self.zinv = MTinv.Zinvariants(z_object=self.Z)

    # --> write edi file
    def write_edi_file(self, new_fn=None, new_Z=None, new_Tipper=None):
        """
        write a new edi file if things have changed.  Note if new_Z or
        new_Tipper are not None, they are not changed in MT object, you
        need to change them manually if you want them to be changed.
        Similarly, the new function name does not change the MT objecte fn
        attribute but does change MT.edi_object.fn attribute.
        
        Arguments
        --------------
            
            *new_fn* : string
                       full path to new file name
                       
            *new_Z* : mtpy.core.Z object
                      a new impedance tensor object to be written
                      
            *new_Tipper* : mtpy.core.Z.Tipper object
                           a new Tipper object to be written
        """

        if new_Z is not None:
            self.edi_object.Z = new_Z
        else:
            self.edi_object.Z = self._Z

        if new_Tipper is not None:
            self.edi_object.Tipper = new_Tipper
        else:
            self.edi_object.Tipper = self._Tipper

        self.edi_object.lat = self._lat
        self.edi_object.lon = self._lon
        self.edi_object.station = self.station
        self.edi_object.zrot = self.rotation_angle

        if new_fn is None:
            new_fn = self.fn[:-4] + "_RW" + ".edi"

        self.edi_object.write_edi_file(new_edi_fn=new_fn)

    # --> check the order of frequencies
    def _check_freq_order(self):
        """
        check to make sure the Z and Tipper arrays are ordered such that
        the first index corresponds to the highest frequency and the last
        index corresponds to the lowest freqeuncy.
        
        """

        if self.Z.freq[0] < self.Z.freq[1]:
            print "Flipping arrays to be ordered from short period to long"
            self.Z.z = self.Z.z.copy()[::-1]
            self.Z.z_err = self.Z.z_err.copy()[::-1]
            self.Z.freq = self.Z.freq.copy()[::-1]

        if self.Tipper.tipper is not None:
            if self.Tipper.freq[0] < self.Tipper.freq[1]:
                self.Tipper.tipper = self.Tipper.tipper.copy()[::-1]
                self.Tipper.tipper_err = self.Tipper.tipper_err.copy()[::-1]
                self.Tipper.freq = self.Tipper.freq.copy()[::-1]

    def remove_distortion(self, num_freq=None):
        """
        remove distortion following Bibby et al. [2005].
        
        Example
        ----------
        :Remove Distortion and Write New .edi: ::
        
            >>> import mtpy.core.mt as mt
            >>> mt1 = mt.MT(fn=r"/home/mt/edi_files/mt01.edi")
            >>> D, new_z = mt1.remove_distortion()
            >>> mt1.write_edi_file(new_fn=r"/home/mt/edi_files/mt01_dr.edi",\
            >>>                    new_Z=new_z)
        """
        dummy_z_obj = MTz.copy.deepcopy(self.Z)
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
            
        
        Arguments
        ------------
        
            *ss_x* : float
                    correction factor for x component
            
            *ss_y* : float
                   correction factor for y component
                   
        .. note:: The factors are in resistivity scale, so the
                  entries of  the matrix "S" need to be given by their
                  square-roots!
                   
        Returns
        ------------
           
           *new_z* : new z array
           

                  
        Examples
        ----------
        :Remove Static Shift: ::
        
            >>> import mtpy.core.mt as mt
            >>> mt_obj = mt.MT(r"/home/mt/mt01.edi")
            >>> new_z_obj = mt.remove_static_shift(ss_x=.5, ss_y=1.2)
            >>> mt_obj.write_edi_file(new_fn=r"/home/mt/mt01_ss.edi",
            >>> ...                   new_Z=new_z_obj)
        """

        s_array, new_z = self.Z.remove_ss(
            reduce_res_factor_x=ss_x, reduce_res_factor_y=ss_y
        )

        new_z_obj = MTz.copy.deepcopy(self.Z)
        new_z_obj.z = new_z

        return new_z_obj

    def interpolate(self, new_freq_array, interp_type="slinear"):
        """
        Interpolate the impedance tensor onto different frequencies
        
        
        Arguments
        ------------
        
            *new_freq_array* : np.ndarray 
                               a 1-d array of frequencies to interpolate on
                               to.  Must be with in the bounds of the existing
                               frequency range, anything outside and an error
                               will occur.
                               
        Returns
        -----------
        
            *new_z_object* : mtpy.core.z.Z object
                             a new impedance object with the corresponding
                             frequencies and components.
                             
            *new_tipper_object* : mtpy.core.z.Tipper object
                             a new tipper object with the corresponding
                             frequencies and components.
                             
       
        Examples
        ----------
       
        :Interpolate: ::
         
            >>> import mtpy.core.mt as mt
            >>> edi_fn = r"/home/edi_files/mt_01.edi"
            >>> mt_obj = mt.MT(edi_fn)
            >>> # create a new frequency range to interpolate onto
            >>> new_freq = np.logspace(-3, 3, 24)
            >>> new_z_object, new_tipper_obj = mt_obj.interpolate(new_freq)
            >>> mt_obj.write_edi_file(new_fn=r"/home/edi_files/mt_01_interp.edi",
            >>> ...                   new_Z=new_z_object,
            >>> ...                   new_Tipper=new_tipper_object)
            
        """
        # if the interpolation module has not been loaded return
        if interp_import is False:
            print ("could not interpolate, need to install scipy")
            return

        # make sure the input is a numpy array
        if type(new_freq_array) != np.ndarray:
            new_freq_array = np.array(new_freq_array)

        # check the bounds of the new frequency array
        if self.Z.freq.min() > new_freq_array.min():
            raise ValueError(
                "New frequency minimum of {0:.5g}".format(new_freq_array.min())
                + " is smaller than old frequency minimum of {0:.5g}".format(
                    self.Z.freq.min()
                )
                + ".  The new frequency range needs to be within the "
                + "bounds of the old one."
            )
        if self.Z.freq.max() < new_freq_array.max():
            raise ValueError(
                "New frequency maximum of {0:.5g}".format(new_freq_array.max())
                + "is smaller than old frequency maximum of {0:.5g}".format(
                    self.Z.freq.max()
                )
                + ".  The new frequency range needs to be within the "
                + "bounds of the old one."
            )

        # make a new Z object
        new_Z = MTz.Z(
            z_array=np.zeros((new_freq_array.shape[0], 2, 2), dtype="complex"),
            z_err_array=np.zeros((new_freq_array.shape[0], 2, 2)),
            freq=new_freq_array,
        )

        new_Tipper = MTz.Tipper(
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
                )
                new_f = new_freq_array[new_nz_index]

                # create a function that does 1d interpolation
                z_func_real = spi.interp1d(f, z_real, kind=interp_type)
                z_func_imag = spi.interp1d(f, z_imag, kind=interp_type)
                z_func_err = spi.interp1d(f, z_err, kind=interp_type)

                # interpolate onto new frequency range
                new_Z.z[new_nz_index, ii, jj] = z_func_real(new_f) + 1j * z_func_imag(
                    new_f
                )
                new_Z.z_err[new_nz_index, ii, jj] = z_func_err(new_f)

        # if there is not tipper than skip
        if self.Tipper.tipper is None:
            return new_Z, new_Tipper

        # interpolate the Tipper
        for jj in range(2):
            # get indicies of non-zero components
            nz_index = np.nonzero(self.Tipper.tipper[:, 0, jj])

            if len(nz_index[0]) == 0:
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

        return new_Z, new_Tipper

    def plot_mt_response(self, **kwargs):
        """ 
        Returns a mtpy.imaging.plotresponse.PlotResponse object
        
        Examples
        ------------
        :Plot Response: ::
        
            >>> mt_obj = mt.MT(edi_file)
            >>> pr = mt.plot_mt_response()
            >>> # if you need more infor on plot_mt_response 
            >>> help(pr)
            
        """

        plot_obj = plotresponse.PlotResponse(fn=self.fn, **kwargs)

        return plot_obj
