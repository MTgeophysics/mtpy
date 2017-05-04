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

#==============================================================================
import mtpy.core.edi as MTedi
import mtpy.core.z as MTz
import mtpy.utils.latlongutmconversion as MTutm
import mtpy.utils.exceptions as MTex
import mtpy.utils.gis_tools as gis_tools
import mtpy.analysis.pt as MTpt
import mtpy.analysis.zinvariants as MTinv
import mtpy.analysis.distortion as MTdistortion
import os
import numpy as np
import mtpy.imaging.plotresponse as plotresponse

try:
    import scipy
    scipy_version = int(scipy.__version__.replace('.', ''))
    
    if scipy_version < 140:
        print ('Note: need scipy version 0.14.0 or higher or interpolation '+\
               'might not work.')
    import scipy.interpolate as spi
    interp_import = True
    
except ImportError:
    print('Could not find scipy.interpolate, cannot use method interpolate'+\
          'check installation you can get scipy from scipy.org.' )
    interp_import = False

#==============================================================================

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
        
        # important information held in objects
        self.Site = Site()
        self.FieldNotes = FieldNotes()
        self.Provenance = Provenance()
        self.Notes = MTedi.Information()
        self.Processing = Processing()

        self._Z = kwargs.pop('Z', MTz.Z())
        self._Tipper = kwargs.pop('Tipper', MTz.Tipper())
        self._rotation_angle = 0
        self._fn = None
        self._edi_obj = MTedi.Edi()
        
        self.fn = fn
        
        #provide key words to fill values if an edi file does not exist
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

    #==========================================================================
    # set functions                        
    #==========================================================================
    def _set_lat(self, latitude):
        """
        set latitude making sure the input is in decimal degrees
        
        upon setting utm coordinates are recalculated
        """
        self.Site.Location.latitude = latitude
        self.Site.Location.project_location2utm()
        
    def _set_lon(self, longitude):
        """
        set longitude making sure the input is in decimal degrees
        
        upon setting utm coordinates are recalculated
        """
        self.Site.Location.longitude = longitude
        self.Site.Location.project_location2utm()
        
        
    def _set_elev(self, elevation):
        """
        set elevation, should be input as meters
        """
        self._elev = self.Site.Location.elevation = elevation

        
    def _set_east(self, easting):
        """
        set easting in meters
        
        upon setting lat and lon are recalculated
        """
        self.Site.Location.easting = easting
        self.Site.Location.project_location2ll()
        
    def _set_north(self, northing):
        """
        set northing in meters
        
        upon setting lat and lon are recalculated
        """
        self.Site.Location.northing = northing
        self.Site.Location.project_location2ll()
    
        
    def _set_utm_zone(self, utm_zone):
        """
        set UTM zone
        
        upon setting lat and lon are recalculated
        """
        self.Site.Location.utm_zone = utm_zone
        self.Site.Location.project_location2ll()
        
    def _set_fn(self, filename):
        """
        set filename, currently only support .edi files
        """
        
        self._fn = filename
        if self._fn is None:
            return 
            
        if self._fn.lower().endswith('.edi'):
            self.read_edi_file()
        else:
            not_fn = self._fn[os.path.basename(self._fn).find['.']:]
            raise MTex.MTpyError_file_handling('File '+\
                              'type {0} not supported yet.'.format(not_fn))
            
    
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
        
        
        print ("Rotated Z, Tipper, Phase Tensor and Zinvariants by"
               "{0:.3f} degrees".format(self._rotation_angle))
               
    def _set_Z(self, z_object):
        """
        set z_object
        
        recalculate phase tensor and invariants, which shouldn't change except
        for strike angle
        """
        
        self._Z = z_object
        self._Z._compute_res_phase()
        
    def _set_Tipper(self, t_object):
        """
        set tipper object
        
        recalculate tipper angle and magnitude
        """
        
        self._Tipper = t_object
        if self._Tipper is not None:
            self._Tipper._compute_amp_phase()
            self._Tipper._compute_mag_direction()
            
    def _set_station(self, station_name):
        """
        set station name
        """
        self.Site.id = station_name
        
    #==========================================================================
    # get functions                         
    #==========================================================================    
    def _get_lat(self):
        return self.Site.Location.latitude

    def _get_lon(self):
        return self.Site.Location.longitude
        
    def _get_elev(self):
        return self.Site.Location.elevation
        
    def _get_east(self):
        return self.Site.Location.easting
        
    def _get_north(self):
        return self.Site.Location.northing
    
    def _get_utm_zone(self):
        return self.Site.Location.utm_zone
    
    def _get_fn(self):
        return self._fn
    
    def _get_rotation_angle(self):
        return self._rotation_angle
    
    def _get_Z(self):
        return self._Z
        
    def _get_Tipper(self):
        return self._Tipper
        
    def _get_station(self):
        return self.Site.id
    
    def _get_pt(self):
        return MTpt.PhaseTensor(z_object=self.Z)
    #==========================================================================
    # set properties                          
    #==========================================================================
    lat = property(_get_lat, _set_lat, 
                   doc="latitude of station in decimal degrees")
    
    lon = property(_get_lon, _set_lon, 
                   doc="longitude of station in decimal degrees")
                   
    elev = property(_get_elev, _set_elev, 
                    doc="elevation in meters")
                   
    east = property(_get_east, _set_east,
                    doc="easting in meters of station location on UTM grid")
    
    north = property(_get_north, _set_north,
                    doc="northing in meters of station location on UTM grid")
                    
    utm_zone = property(_get_utm_zone, _set_utm_zone,
                        doc="UTM zone")
                        
    fn = property(_get_fn, _set_fn, doc="name of file containing MT info")
    
    rotation_angle = property(_get_rotation_angle, _set_rotation_angle,
                              doc="rotation angle of Z and Tipper")
                              
    Z = property(_get_Z, _set_Z, doc="impedence tensor object")
    
    Tipper = property(_get_Tipper, _set_Tipper, doc="Tipper object")
    
    station = property(_get_station, _set_station, doc="Station name")
    
    pt = property(_get_pt, doc="""Phase Tensor, can only get, not set""")
                                       
    #==========================================================================
    #  read in files   
    #==========================================================================
    #--> read in edi file                                                    
    def read_edi_file(self):
        """
        read in edi file and set attributes accordingly
        
        """        
        edi_obj = MTedi.Edi(edi_fn=self.fn)

        # read in site information from the header
        self.lat = edi_obj.lat
        self.lon = edi_obj.lon
        self.elev = edi_obj.elev
        self.Site.acquired_by = edi_obj.Header.acqby
        self.Site.survey = edi_obj.Header.loc
        self.Site.start_date = edi_obj.Header.acqdate
        self.Site.project = edi_obj.Header.project
        self.Site.Location.datum = edi_obj.Header.datum
        self.Site.Location.elev_units = edi_obj.Define_measurement.units
        self.Site.Location.coordinate_system = edi_obj.Header.coordinate_system
        
        # get information about different sensors
        try:
            for key in edi_obj.Define_measurement.meas_hx.__dict__.keys():
                setattr(self.FieldNotes.magnetometer_hx,
                        key, 
                        edi_obj.Define_measurement.meas_hx.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in edi_obj.Define_measurement.meas_hy.__dict__.keys():
                setattr(self.FieldNotes.magnetometer_hy,
                        key, 
                        edi_obj.Define_measurement.meas_hy.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in edi_obj.Define_measurement.meas_hz.__dict__.keys():
                setattr(self.FieldNotes.magnetometer_hz,
                        key, 
                        edi_obj.Define_measurement.meas_hz.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in edi_obj.Define_measurement.meas_ex.__dict__.keys():
                setattr(self.FieldNotes.electrode_ex,
                        key, 
                        edi_obj.Define_measurement.meas_ex.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in edi_obj.Define_measurement.meas_ey.__dict__.keys():
                setattr(self.FieldNotes.electrode_ey,
                        key, 
                        edi_obj.Define_measurement.meas_ey.__dict__[key])
        except AttributeError:
            pass
        
        # get info 
        self.Notes = edi_obj.Info
        try: 
            self.FieldNotes.magnetometer_hx.id = self.Notes.hx
        except AttributeError:
            pass
        try: 
            self.FieldNotes.magnetometer_hy.id = self.Notes.hy
        except AttributeError:
            pass
        try: 
            self.FieldNotes.magnetometer_hz.id = self.Notes.hz
        except AttributeError:
            pass
        
        try:
            self.FieldNotes.magnetometer_hx.type = self.Notes.b_instrument_type
            self.FieldNotes.magnetometer_hy.type = self.Notes.b_instrument_type
            self.FieldNotes.magnetometer_hz.type = self.Notes.b_instrument_type
        except AttributeError:
            pass
        
        try:
            self.FieldNotes.electrode_ex.type = self.Notes.e_instrument_type
            self.FieldNotes.electrode_ey.type = self.Notes.e_instrument_type
        except AttributeError:
            pass
        
        try:
            self.FieldNotes.data_logger = self.Notes.b_logger_type
        except AttributeError:
            pass
                
        
        self.Z = edi_obj.Z
        self.Tipper = edi_obj.Tipper
        self.station = edi_obj.station
        
        #--> make sure things are ordered from high frequency to low
        self._check_freq_order()
        
        self._edi_obj = edi_obj
        
    #--> write edi file 
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
        
        # get header information, mostly from site
        self._edi_obj = MTedi.Edi()
        self._edi_obj.Header = self._get_edi_header()
        
        # get information 
        self._edi_obj.Info = MTedi.Information()
        self._edi_obj.Info.info_list = self._get_edi_info_list()
        
        # get define measurement
        self._edi_obj.Define_measurement = self._get_edi_define_measurement()
        
        # get mtsec
        self._edi_obj.Data_sect = self._get_edi_data_sect()
        
        if new_Z is not None:
            self._edi_obj.Z = new_Z
        else:
            self._edi_obj.Z = self._Z
       
        if new_Tipper is not None:
            self._edi_obj.Tipper = new_Tipper
        else:
            self._edi_obj.Tipper = self._Tipper
            
        self._edi_obj.zrot = self.rotation_angle
        
        if new_fn is None:
            new_fn = self.fn[:-4]+'_RW'+'.edi'
            
        self._edi_obj.write_edi_file(new_edi_fn=new_fn)
        
    def _get_edi_header(self):
        """
        make an edi header class
        """
        
        header = MTedi.Header()
        header.acqby = self.Site.acquired_by
        header.acqdate = self.Site.start_date
        header.coordinate_system = self.Site.Location.coordinate_system
        header.dataid = self.Site.id
        header.datum = self.Site.Location.datum
        header.elev = self.elev
        header.fileby = self.Site.acquired_by
        header.lat = self.lat
        header.loc = self.Site.project
        header.lon = self.lon
        header.project = self.Site.project
        header.survey = self.Site.survey
        header.units = self.Site.Location.elev_units
        
        return header
        
    #--> get information list for edi
    def _get_edi_info_list(self):
        """
        get the information for an edi file
        """
        
        info_list = []
        # write previous information first
        for key in self.Notes.__dict__.keys():
            if key.lower() not in ['edi_lines', 'info_list', 'edi_fn']:
                l_key = key.lower()
                l_value = getattr(self.Notes, key)
                info_list.append('{0} = {1}'.format(l_key, l_value))
        
        # get instrument information
        for key in ['ex', 'ey', 'hx', 'hy', 'hz']:
            if 'e' in key:
                key = 'electrode_{0}'.format(key)
            elif 'h' in key:
                key = 'magnetometer_{0}'.format(key)
                
            instrument_obj = getattr(self.FieldNotes, key) 
            for mkey in ['manufacturer', 'id', 'type']:
                l_key = '{0}.{1}'.format(key, mkey)
                line = '{0} = {1}'.format(l_key.lower(), 
                                          getattr(instrument_obj, mkey))
                info_list.append(line)
                
        # get processing information  
        for p_key in self.Processing.__dict__.keys():
            if p_key.lower() == 'software':
                for s_key in self.Processing.Software.__dict__.keys():
                    l_key = 'processing.software.{0}'.format(s_key)
                    l_value = getattr(self.Processing.Software, s_key)
                    info_list.append('{0} = {1}'.format(l_key.lower(), l_value))
            else:
                l_key = 'processing.{0}'.format(p_key)
                l_value = getattr(self.Processing, p_key)
                info_list.append('{0} = {1}'.format(l_key, l_value))
            
        return info_list
    
    # get edi define measurement
    def _get_edi_define_measurement(self):
        """
        get define measurement block for an edi file
        """       
        
        define_meas = MTedi.DefineMeasurement()
        define_meas.refelev = self.elev
        define_meas.reflat = self.lat
        define_meas.reflon = self.lon
        define_meas.reftype = self.Site.Location.coordinate_system
        define_meas.units = self.Site.Location.elev_units
        
        define_meas.meas_ex = MTedi.EMeasurement()
        for key in define_meas.meas_ex._kw_list:
            setattr(define_meas.meas_ex,
                    key, 
                    getattr(self.FieldNotes.electrode_ex, key))
        
        define_meas.meas_ey = MTedi.EMeasurement()
        for key in define_meas.meas_ey._kw_list:
            setattr(define_meas.meas_ey,
                    key, 
                    getattr(self.FieldNotes.electrode_ey, key))
            
        define_meas.meas_hx = MTedi.HMeasurement()
        for key in define_meas.meas_hx._kw_list:
            setattr(define_meas.meas_hx,
                    key, 
                    getattr(self.FieldNotes.magnetometer_hx, key))
 
        define_meas.meas_hy = MTedi.HMeasurement()
        for key in define_meas.meas_hy._kw_list:
            setattr(define_meas.meas_hy,
                    key, 
                    getattr(self.FieldNotes.magnetometer_hy, key)) 
            
        if np.all(self.Tipper.tipper == 0) == False:    
            define_meas.meas_hz = MTedi.HMeasurement()
            for key in define_meas.meas_hz._kw_list:
                setattr(define_meas.meas_hz,
                        key, 
                        getattr(self.FieldNotes.magnetometer_hz, key))
            
        return define_meas
    
    def _get_edi_data_sect(self):
        """
        get mt data section for edi file
        """
        
        sect = MTedi.DataSection()
        sect.data_type = 'MT'
        sect.nfreq = self.Z.z.shape[0]
        sect.sectid = self.station
        nchan = 5
        if np.all(self.Tipper.tipper == 0) == True:
            nchan = 4
        sect.nchan = nchan
        sect.maxblks = 999
        sect.ex = self.FieldNotes.electrode_ex.id
        sect.ey = self.FieldNotes.electrode_ey.id
        sect.hx = self.FieldNotes.magnetometer_hx.id
        sect.hy = self.FieldNotes.magnetometer_hy.id
        if np.all(self.Tipper.tipper == 0) == False:
            sect.hz = self.FieldNotes.magnetometer_hz.id
            
        return sect
            
        
    #--> check the order of frequencies
    def _check_freq_order(self):
        """
        check to make sure the Z and Tipper arrays are ordered such that
        the first index corresponds to the highest frequency and the last
        index corresponds to the lowest freqeuncy.
        
        """

        if self.Z.freq[0] < self.Z.freq[1]:
            print 'Flipping arrays to be ordered from short period to long'
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
        D, new_z_object = MTdistortion.remove_distortion(z_object=dummy_z_obj,
                                                         num_freq=num_freq)
        
        return D, new_z_object
        
    def remove_static_shift(self, ss_x=1.0, ss_y =1.0):
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
        
        s_array, new_z = self.Z.remove_ss(reduce_res_factor_x=ss_x,
                                          reduce_res_factor_y=ss_y)
                                      
        new_z_obj = MTz.copy.deepcopy(self.Z)
        new_z_obj.z = new_z
        
        return new_z_obj
        
        
    def interpolate(self, new_freq_array, interp_type='slinear'):
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
            print('could not interpolate, need to install scipy')            
            return
        
        #make sure the input is a numpy array
        if type(new_freq_array) != np.ndarray:
            new_freq_array = np.array(new_freq_array)

        # check the bounds of the new frequency array
        if self.Z.freq.min() > new_freq_array.min():
            raise ValueError('New frequency minimum of {0:.5g}'.format(new_freq_array.min())+\
                             ' is smaller than old frequency minimum of {0:.5g}'.format(self.Z.freq.min())+\
                             '.  The new frequency range needs to be within the '+\
                             'bounds of the old one.')
        if self.Z.freq.max() < new_freq_array.max():
            raise ValueError('New frequency maximum of {0:.5g}'.format(new_freq_array.max())+\
                             'is smaller than old frequency maximum of {0:.5g}'.format(self.Z.freq.max())+\
                             '.  The new frequency range needs to be within the '+\
                             'bounds of the old one.')

        # make a new Z object
        new_Z = MTz.Z(z_array=np.zeros((new_freq_array.shape[0], 2, 2), 
                                       dtype='complex'),
                      z_err_array=np.zeros((new_freq_array.shape[0], 2, 2)), 
                      freq=new_freq_array)
                      
        new_Tipper = MTz.Tipper(tipper_array=np.zeros((new_freq_array.shape[0], 1, 2), 
                                             dtype='complex'),
                      tipper_err_array=np.zeros((new_freq_array.shape[0], 1, 2)), 
                      freq=new_freq_array)
        
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
                new_nz_index = np.where((new_freq_array >= f.min()) & 
                                        (new_freq_array <= f.max()))
                new_f = new_freq_array[new_nz_index]
                
                # create a function that does 1d interpolation
                z_func_real = spi.interp1d(f, z_real, kind=interp_type)
                z_func_imag = spi.interp1d(f, z_imag, kind=interp_type)
                z_func_err = spi.interp1d(f, z_err, kind=interp_type)
                
                # interpolate onto new frequency range
                new_Z.z[new_nz_index, ii, jj] = z_func_real(new_f)+1j*z_func_imag(new_f)
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
            new_nz_index = np.where((new_freq_array >= f.min()) & 
                                    (new_freq_array <= f.max()))
            new_f = new_freq_array[new_nz_index] 
                                                 
            # interpolate onto new frequency range
            new_Tipper.tipper[new_nz_index, 0, jj] = t_func_real(new_f)+\
                                                  1j*t_func_imag(new_f)

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
        
#==============================================================================
# Site details    
#==============================================================================
class Site(object):
    """
    Information on the site, including location, id, etc.
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    aqcuired_by       string       name of company or person whom aqcuired the
                                   data.
    id                string       station name
    Location          object       Holds location information, lat, lon, elev
                      Location     datum, easting, northing see Location class  
    start_date        string       YYYY-MM-DD start date of measurement
    end_date          string       YYYY-MM-DD end date of measurement
    year_collected    string       year data collected
    survey            string       survey name
    project           string       project name
    run_list          string       list of measurment runs ex. [mt01a, mt01b]
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>> Site(**{'state':'Nevada', 'Operator':'MTExperts'})

    """
    
    def __init__(self, **kwargs):
        
        self.acquired_by = None
        self.end_date = None
        self.id = None
        self.Location = Location()
        self.project = None
        self.run_list = None
        self.start_date = None
        self.survey = None
        self.year_collected = None
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
#==============================================================================
# Location class, be sure to put locations in decimal degrees, and note datum            
#==============================================================================
class Location(object):
    """
    location details
    """
    
    def __init__(self, **kwargs):
        self.datum = 'WGS84'
        self.declination = None
        self.declination_epoch = None
        
        self._elevation = None
        self._latitude = None
        self._longitude = None
        
        self.northing = None
        self.easting = None
        self.utm_zone = None
        self.elev_units = 'm'
        self.coordinate_system = 'Geographic North'
        
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    def _get_latitude(self):
        return self._latitude
    def _set_latitude(self, lat):
        self._latitude = gis_tools.assert_lat_value(lat)
        
    latitude = property(fget=_get_latitude, 
                        fset=_set_latitude,
                        doc="""Latitude in decimal degrees""") 
                        
    def _get_longitude(self):
        return self._longitude
    def _set_longitude(self, lon):
        self._longitude = gis_tools.assert_lon_value( lon)
        
    longitude = property(fget=_get_longitude, 
                        fset=_set_longitude,
                        doc="""Longitude in decimal degrees""") 
                        
    def _get_elevation(self):
        return self._elevation
    def _set_elevation(self, elev):
        self._elevation = gis_tools.assert_elevation_value(elev)
        
    elevation = property(fget=_get_elevation, 
                        fset=_set_elevation,
                        doc="""Elevation in floating point""") 
    
    def project_location2utm(self):
        """
        project location coordinates into meters given the reference ellipsoid,
        for now that is constrained to WGS84
        
        Returns East, North, Zone
        """
        utm_point = gis_tools.project_point_ll2utm(self._latitude, 
                                                   self._longitude,
                                                   datum=self.datum)
        
        self.easting = utm_point[0]
        self.northing = utm_point[1]
        self.utm_zone = utm_point[2]
        
    def project_location2ll(self):
        """
        project location coordinates into meters given the reference ellipsoid,
        for now that is constrained to WGS84
        
        Returns East, North, Zone
        """
        ll_point = gis_tools.project_point_utm2ll(self.easting, 
                                                   self.northing,
                                                   self.utm_zone,
                                                   datum=self.datum)
        
        self.latitude = ll_point[0]
        self.longitude = ll_point[1]
 
#==============================================================================
# Field Notes    
#==============================================================================
class FieldNotes(object):
    """
    Field note information.
    
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    data_quality      DataQuality notes on data quality
    electrode         Instrument      type of electrode used
    data_logger       Instrument      type of data logger
    magnetometer      Instrument      type of magnetotmeter
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>> FieldNotes(**{'electrode_ex':'Ag-AgCl 213', 'magnetometer_hx':'102'})
    """
    
    def __init__(self, **kwargs):
        null_emeas = MTedi.EMeasurement()
        null_hmeas = MTedi.HMeasurement()
        
        self.data_quality = DataQuality()
        self.data_logger = Instrument()
        self.electrode_ex = Instrument(**null_emeas.__dict__)
        self.electrode_ey = Instrument(**null_emeas.__dict__)
        self.magnetometer_hx = Instrument(**null_hmeas.__dict__)
        self.magnetometer_hy = Instrument(**null_hmeas.__dict__)
        self.magnetometer_hz = Instrument(**null_hmeas.__dict__)

        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            

#==============================================================================
# Instrument            
#==============================================================================
class Instrument(object):
    """
    Information on an instrument that was used.
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    id                string      serial number or id number of data logger
    manufacturer      string      company whom makes the instrument
    type              string      Broadband, long period, something else 
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>> Instrument(**{'ports':'5', 'gps':'time_stamped'})
    """
    
    def __init__(self, **kwargs):
        self.id = None
        self.manufacturer = None
        self.type = None
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

#==============================================================================
# Data Quality            
#==============================================================================
class DataQuality(object):
    """
    Information on data quality.
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    comments          string      comments on data quality
    good_from_period  float       minimum period data are good 
    good_to_period    float       maximum period data are good
    rating            int         [1-5]; 1 = poor, 5 = excellent
    warrning_comments string      any comments on warnings in the data
    warnings_flag     int         [0-#of warnings]       
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>>DataQuality(**{'time_series_comments':'Periodic Noise'})
    """
    
    def __init__(self, **kwargs):
        self.comments = None
        self.good_from_period = None
        self.good_to_period = None
        self.rating = None
        self.warnings_comments = None
        self.warnings_flag = 0
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            

#==============================================================================
# Citation        
#==============================================================================
class Citation(object):
    """
    Information for a citation.
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    author            string      Author names
    title             string      Title of article, or publication
    journal           string      Name of journal
    doi               string      DOI number (doi:10.110/sf454)
    year              int         year published
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>> Citation(**{'volume':56, 'pages':'234--214'})
    """
    
    def __init__(self, **kwargs):
        self.author = None
        self.title = None
        self.journal = None
        self.volume = None
        self.doi = None
        self.year = None
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            

#==============================================================================
# Copyright            
#==============================================================================
class Copyright(object):
    """
    Information of copyright, mainly about how someone else can use these 
    data. Be sure to read over the conditions_of_use.
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    citation          Citation    citation of published work using these data
    conditions_of_use string      conditions of use of these data
    release_status    string      release status [ open | public | proprietary]
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>> Copyright(**{'owner':'University of MT', 'contact':'Cagniard'})
    """
    
    def __init__(self, **kwargs):
        self.citation = Citation()
        self.conditions_of_use = ''.join(['All data and metadata for this survey are ',
                                          'available free of charge and may be copied ',
                                          'freely, duplicated and further distributed ',
                                          'provided this data set is cited as the ',
                                          'reference. While the author(s) strive to ',
                                          'provide data and metadata of best possible ',
                                          'quality, neither the author(s) of this data ',
                                          'set, not IRIS make any claims, promises, or ', 
                                          'guarantees about the accuracy, completeness, ',
                                          'or adequacy of this information, and expressly ',
                                          'disclaim liability for errors and omissions in ',
                                          'the contents of this file. Guidelines about ',
                                          'the quality or limitations of the data and ',
                                          'metadata, as obtained from the author(s), are ',
                                          'included for informational purposes only.'])
        self.release_status = None
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
#==============================================================================
# Provenance        
#==============================================================================
class Provenance(object):
    """
    Information of the file history, how it was made
    
    Holds the following information:
    
    ====================== =========== ========================================
    Attributes             Type        Explanation    
    ====================== =========== ========================================
    creation_time          string      creation time of file YYYY-MM-DD,hh:mm:ss
    creating_application   string      name of program creating the file
    creator                Person      person whom created the file
    submitter              Person      person whom is submitting file for 
                                       archiving 
    ====================== =========== ========================================

    More attributes can be added by inputing a key word dictionary
    
    >>> Provenance(**{'archive':'IRIS', 'reprocessed_by':'grad_student'}) 
    """
    
    def __init__(self, **kwargs):
        
        self.creation_time = None
        self.creating_application = None
        self.creator = Person()
        self.submitter = Person()
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
 
#==============================================================================
# Person
#==============================================================================
class Person(object):
    """
    Information for a person
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    email             string      email of person 
    name              string      name of person
    organization      string      name of person's organization
    organization_url  string      organizations web address
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>> Person(**{'phone':'650-888-6666'})
    """    
    
    def __init__(self, **kwargs):
        self.email = None
        self.name = None
        self.organization = None
        self.organization_url = None
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
#==============================================================================
# Processing
#==============================================================================
class Processing(object):
    """
    Information for a processing
    
    Holds the following information:
    
    ================= =========== =============================================
    Attributes         Type        Explanation    
    ================= =========== =============================================
    email             string      email of person 
    name              string      name of person
    organization      string      name of person's organization
    organization_url  string      organizations web address
    ================= =========== =============================================

    More attributes can be added by inputing a key word dictionary
    
    >>> Person(**{'phone':'650-888-6666'})
    """    
    
    def __init__(self, **kwargs):
        self.Software = Software()
        self.author = None
        self.organization = None
        self.organization_url = None
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
        
class Software(object):
    """
    software
    """        
    
    def __init__(self, **kwargs):
        self.name = None
        self.version = None
        
        for key in kwargs:
            setattr(self, key, kwargs[key])
        