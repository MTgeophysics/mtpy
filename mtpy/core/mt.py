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
import mtpy.utils.format as MTformat
import mtpy.analysis.pt as MTpt
import mtpy.analysis.zinvariants as MTinv
import mtpy.analysis.distortion as MTdistortion
import os

#==============================================================================

class MT(object):
    """
    Basic object containing all information necessary for a single MT station
    including:
    
    ===================== =====================================================
    **Attribute**         Description
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
            - zerr_array      --> np.ndarray(n_freq, 2, 2)
            - freq            --> np.ndarray(n_freq)
            - resistivity     --> np.ndarray(n_freq, 2, 2) (linear scale)
            - resistivity_err --> np.ndarray(n_freq, 2, 2) 
            - phase           --> np.ndarray(n_freq, 2, 2)           
            - phase_err       --> np.ndarray(n_freq, 2, 2) 
            - tipper_object   --> mtpy.core.z.Tipper object
            - tipper          --> np.ndarray(n_freq, 1, 2, dtype='complex') 
            - tippererr       --> np.ndarray(n_freq, 1, 2)
        
    ===================== =====================================================
    **methods**           Description
    ===================== =====================================================
    write_edi_file        write an edi_file from the MT data
    remove_distortion     remove distortion from the data following 
                          Bibby et al. [2005]
    ===================== =====================================================
    
        
    
    """
    
    def __init__(self, fn=None, **kwargs):
        
        self._fn = fn
        self.station = kwargs.pop('station', None)
        self._lat = kwargs.pop('lat', None)
        self._lon = kwargs.pop('lon', None)
        self.elev = kwargs.pop('elev', None)
        self._Z = kwargs.pop('Z', MTz.Z())
        self.Tipper = kwargs.pop('Tipper', MTz.Tipper())
        self._utm_zone = kwargs.pop('utm_zone', None)
        self._east = kwargs.pop('east', None)
        self._north = kwargs.pop('north', None)
        self._rotation_angle = kwargs.pop('rotation_angle', 0)
        
        #provide key words to fill values if an edi file does not exist
        if 'z_object' in kwargs:
            self._Z = kwargs['z_object']
            
        if 'z_array' in kwargs:
            self._Z.z = kwargs['z_array']
        
        if 'zerr_array' in kwargs:
            self._Z.zerr = kwargs['zerr_array']
        
        if 'freq' in kwargs:
            self._Z.freq = kwargs['freq']
            self.Tipper.freq = kwargs['freq']
            
        if 'tipper_object' in kwargs:
            self.Tipper = kwargs['tipper_object']
            
        if 'tipper' in kwargs:
            self.Tipper.tipper = kwargs['tipper']
        
        if 'tippererr' in kwargs:
            self.Tipper.tippererr = kwargs['tippererr']
            
        if 'resisitivity' in kwargs:
            self._Z.resistivity = kwargs['resistivity']
        
        if 'resisitivity_err' in kwargs:
            self._Z.resistivity_err = kwargs['resistivity_err']
        
        if 'phase' in kwargs:
            self._Z.phase = kwargs['phase']
            
        if 'phase_err' in kwargs:
            self._Z.phase = kwargs['phase_err']
        
        
        self.edi_object = MTedi.Edi()
        self.pt = None
        self.zinv = None
        self._utm_ellipsoid = 23

        #--> read in the edi file if its given
        if self._fn is not None:
            if self._fn[-3:] == 'edi':
                self._read_edi_file()
            else:
                not_fn = self._fn[os.path.basename(self._fn).find['.']:]
                raise MTex.MTpyError_file_handling('File '+\
                          'type {0} not supported yet.'.format(not_fn))
    
    #==========================================================================
    # set functions                        
    #==========================================================================
    def _set_lat(self, latitude):
        """
        set latitude making sure the input is in decimal degrees
        
        upon setting utm coordinates are recalculated
        """
        
        self._lat = MTformat._assert_position_format('lat', latitude)
        
        if self._lon is not None and self._lat is not None:
            self._get_utm()
        
    def _set_lon(self, longitude):
        """
        set longitude making sure the input is in decimal degrees
        
        upon setting utm coordinates are recalculated
        """
        
        self._lon = MTformat._assert_position_format('lon', longitude)
        
        if self._lon is not None and self._lat is not None:
            self._get_utm()
        
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
        if self._fn[-3:] == 'edi':
            self._read_edi_file()
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
        self.Z.rotate(theta_r)
        self.Tipper.rotate(theta_r)
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
        
        #--> compute phase tensor
        self.pt = MTpt.PhaseTensor(z_object=self._Z, freq=self._Z.freq)
        
        #--> compute invariants 
        self.zinv = MTinv.Zinvariants(z_object=self._Z) 
    #==========================================================================
    # get functions                         
    #==========================================================================    
    def _get_lat(self):
        return self._lat

    def _get_lon(self):
        return self._lon
        
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
        
    #==========================================================================
    # set properties                          
    #==========================================================================
    lat = property(_get_lat, _set_lat, 
                   doc="latitude of station in decimal degrees")
    
    lon = property(_get_lon, _set_lon, 
                   doc="longitude of station in decimal degrees")
                   
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
    
    #--> conversion between utm and ll
    def _get_utm(self):
        """
        get utm coordinates from lat and lon
        """
        
        self.utm_zone, self.east, self.north = MTutm.LLtoUTM(self._utm_ellipsoid,
                                                             self.lat, self.lon)
                                                         
    def _get_ll(self):
        """
        get lat and long from utm
        """
        
        self.lat, self.lon = MTutm.UTMtoLL(self._utm_ellipsoid, 
                                           self.north, 
                                           self.east, 
                                           self.utm_zone)
                                           
    #--> read in edi file                                                    
    def _read_edi_file(self):
        """
        read in edi file and set attributes accordingly
        
        """
        
        self.edi_object = MTedi.Edi(self.fn)
        self.lat = self.edi_object.lat
        self.lon = self.edi_object.lon
        self.elev = self.edi_object.elev
        self.Z = self.edi_object.Z
        self.Tipper = self.edi_object.Tipper
        self.station = self.edi_object.station
        
        #--> get utm coordinates from lat and lon        
        self._get_utm()
        
        #--> make sure things are ordered from high frequency to low
        self._check_freq_order()
        
        #--> compute phase tensor
        self.pt = MTpt.PhaseTensor(z_object=self.Z, freq=self.Z.freq)
        
        #--> compute invariants 
        self.zinv = MTinv.Zinvariants(z_object=self.Z)
        
    #--> write edi file 
    def write_edi_file(self, new_fn=None, new_Z=None):
        """
        write a new edi file if things have changed.
        """
        
        if new_Z is not None:
            self.edi_object.Z = new_Z
        else:
            self.edi_object.Z = self._Z
        self.edi_object.Tipper = self.Tipper
        self.edi_object.lat = self._lat
        self.edi_object.lon = self._lon
        self.edi_object.station = self.station
        self.edi_object.zrot = self.rotation_angle
        
        if new_fn is None:
            new_fn = self.fn[:-4]+'_RW'+'.edi'
            
        self.edi_object.writefile(new_fn)
        
        
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
            self.Z.zerr = self.Z.zerr.copy()[::-1]
            self.Z.freq = self.Z.freq.copy()[::-1]
            
        if self.Tipper.tipper is not None:
            if self.Tipper.freq[0] < self.Tipper.freq[1]:
                self.Tipper.tipper = self.Tipper.tipper.copy()[::-1]
                self.Tipper.tippererr = self.Tipper.tippererr.copy()[::-1]
                self.Tipper.freq = self.Tipper.freq.copy()[::-1]
                
    def remove_distortion(self):
        """
        remove distortion following Bibby et al. [2005].
        
        if you want to write a new edi file with distortion removed you can 
        do this by:
        
            >>> import mtpy.core.mt as mt
            >>> mt1 = mt.MT(fn=r"/home/mt/edi_files/mt01.edi")
            >>> D, new_z = mt1.remove_distortion()
            >>> mt1.write_edi_file(new_fn=r"/home/mt/edi_files/mt01_dr.edi",\
                                   new_Z=new_z)
        """
        
        D, new_z_object = MTdistortion.remove_distortion(z_object=self.Z)
        
        return D, new_z_object
        
    
        
        
        
    