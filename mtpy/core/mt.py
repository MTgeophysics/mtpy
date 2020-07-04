# -*- coding: utf-8 -*-
"""
.. module:: MT
   :synopsis: The main container for MT response functions.

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
"""

# ==============================================================================
import numpy as np
import os
import time
import warnings
from dateutil import parser as dt_parser
from pathlib import Path

import mtpy.core.edi as MTedi
import mtpy.core.z as MTz
import mtpy.utils.gis_tools as gis_tools
import mtpy.analysis.pt as MTpt
import mtpy.analysis.distortion as MTdistortion
import mtpy.core.jfile as MTj
import mtpy.core.mt_xml as MTxml
import mtpy.core.zmm as MTzmm

from mtpy.utils.mtpylog import MtPyLog

_logger = MtPyLog.get_mtpy_logger(__name__)
# _logger.setLevel(logging.DEBUG)

try:
    import scipy

    scipy_version = [int(ss) for ss in scipy.__version__.split('.')]
    if scipy_version[0] == 0:
        if scipy_version[1] < 14:
            warnings.warn('Note: need scipy version 0.14.0 or higher or interpolation '
                          'might not work.', ImportWarning)
            _logger.warning('Note: need scipy version 0.14.0 or higher or interpolation '
                            'might not work.')
    import scipy.interpolate as spi

    interp_import = True

except ImportError:  # pragma: no cover
    warnings.warn('Could not find scipy.interpolate, cannot use method interpolate'
                  'check installation you can get scipy from scipy.org.')
    _logger.warning('Could not find scipy.interpolate, cannot use method interpolate'
                    'check installation you can get scipy from scipy.org.')
    interp_import = False


# =============================================================================
class MT(object):
    """
    Basic MT container to hold all information necessary for a MT station
    including the following parameters.

        * Site --> information on site details (lat, lon, name, etc)
        * FieldNotes --> information on instruments, setup, etc.
        * Copyright --> information on how the data can be used and citations
        * Provenance --> where the data come from and how they are stored
        * Processing --> how the data were processed.

    The most used attributes are made available from MT, namely the following.

    ===================== =====================================================
    Attribute             Description
    ===================== =====================================================
    station               station name
    lat                   station latitude in decimal degrees
    lon                   station longitude in decimal degrees
    elev                  station elevation in meters
    Z                     mtpy.core.z.Z object for impedance tensor
    Tipper                mtpy.core.z.Tipper object for tipper
    pt                    mtpy.analysis.pt.PhaseTensor for phase tensor
    east                  station location in UTM coordinates assuming WGS-84
    north                 station location in UTM coordinates assuming WGS-84
    utm_zone              zone of UTM coordinates assuming WGS-84
    rotation_angle        rotation angle of the data
    fn                    absolute path to the data file
    ===================== =====================================================

    Other information is contained with in the different class attributes. For
    instance survey name is in MT.Site.survey

    .. note::

        * The best way to see what all the information is and where it is
          contained would be to write out a configuration file ::

              >>> import mtpy.core.mt as mt
              >>> mt_obj = mt.MT()
              >>> mt_obj.write_cfg_file(r"/home/mt/generic.cfg")

        * Currently EDI, XML, and j file are supported to read in information,
          and can write out EDI and XML formats.  Will be extending to j and
          Egberts Z format.

    ===================== =====================================================
    Methods               Description
    ===================== =====================================================
    read_mt_file          read in a MT file [ EDI | XML | j ]
    write_mt_file         write a MT file [ EDI | XML ]
    read_cfg_file         read a configuration file
    write_cfg_file        write a configuration file
    remove_distortion     remove distortion  following Bibby et al. [2005]
    remove_static_shift   Shifts apparent resistivity curves up or down
    interpolate           interpolates Z and T onto specified frequency array.
    ===================== =====================================================


    Examples
    -------------------
    :Read from an .edi File: ::

        >>> import mtpy.core.mt as mt
        >>> mt_obj = mt.MT(r"/home/edi_files/s01.edi")

    :Remove Distortion: ::

        >>> import mtpy.core.mt as mt
        >>> mt1 = mt.MT(fn=r"/home/mt/edi_files/mt01.edi")
        >>> D, new_z = mt1.remove_distortion()
        >>> mt1.write_mt_file(new_fn=r"/home/mt/edi_files/mt01_dr.edi",\
        >>>                    new_Z=new_z)

    :Remove Static Shift: ::

        >>> new_z_obj = mt_obj.remove_static_shift(ss_x=.78, ss_y=1.1)
        >>> # write a new edi file
        >>> mt_obj.write_mt_file(new_fn=r"/home/edi_files/s01_ss.edi",
        >>>                       new_Z=new_z)
        >>> wrote file to: /home/edi_files/s01_ss.edi

    :Interpolate: ::

        >>> new_freq = np.logspace(-3, 3, num=24)
        >>> new_z_obj, new_tipper_obj = mt_obj.interpolate(new_freq)
        >>> mt_obj.write_mt_file(new_Z=new_z_obj, new_Tipper=new_tipper_obj)
        >>> wrote file to: /home/edi_files/s01_RW.edi
    """

    def __init__(self, fn=None, **kwargs):
        self._logging = MtPyLog.get_mtpy_logger(self.__class__.__name__)
        # important information held in objects
        self.Site = Site()
        self.FieldNotes = FieldNotes()
        self.Provenance = Provenance()
        self.Notes = MTedi.Information()
        self.Processing = Processing()
        self.Copyright = Copyright()

        self._Z = MTz.Z()
        self._Tipper = MTz.Tipper()
        self._rotation_angle = 0

        self.save_dir = os.getcwd()
        self.original_file_type = None
        if fn is not None:
            self.read_mt_file(fn)
            self._fn = os.path.normpath(os.path.abspath(fn))  # store file reference

        # provide key words to fill values if an edi file does not exist
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    # ==========================================================================
    # get functions
    # ==========================================================================
    @property
    def fn(self):
        """ reference to original data file"""
        return self._fn

    @property
    def lat(self):
        """Latitude"""
        return self.Site.Location.latitude

    @property
    def lon(self):
        """Longitude"""
        return self.Site.Location.longitude

    @property
    def elev(self):
        """Elevation"""
        return self.Site.Location.elevation

    @property
    def east(self):
        """easting (m)"""
        return self.Site.Location.easting

    @property
    def north(self):
        """northing (m)"""
        return self.Site.Location.northing

    @property
    def utm_zone(self):
        """utm zone"""
        return self.Site.Location.utm_zone

    @property
    def rotation_angle(self):
        """rotation angle in degrees from north"""
        return self._rotation_angle

    @property
    def Z(self):
        """mtpy.core.z.Z object to hole impedance tensor"""
        return self._Z

    @property
    def Tipper(self):
        """mtpy.core.z.Tipper object to hold tipper information"""
        return self._Tipper

    @property
    def station(self):
        """station name"""
        return self.Site.id

    @property
    def pt(self):
        """mtpy.analysis.pt.PhaseTensor object to hold phase tensor"""
        return MTpt.PhaseTensor(z_object=self.Z)

    # ==========================================================================
    # set functions
    # ==========================================================================
    @lat.setter
    def lat(self, latitude):
        """
        set latitude making sure the input is in decimal degrees

        upon setting utm coordinates are recalculated
        """
        self.Site.Location.latitude = latitude
        self.Site.Location.project_location2utm()

    @lon.setter
    def lon(self, longitude):
        """
        set longitude making sure the input is in decimal degrees

        upon setting utm coordinates are recalculated
        """
        self.Site.Location.longitude = longitude
        self.Site.Location.project_location2utm()

    @elev.setter
    def elev(self, elevation):
        """
        set elevation, should be input as meters
        """
        self.Site.Location.elevation = elevation

    @east.setter
    def east(self, easting):
        """
        set easting in meters

        upon setting lat and lon are recalculated
        """
        self.Site.Location.easting = easting
        self.Site.Location.project_location2ll()

    @north.setter
    def north(self, northing):
        """
        set northing in meters

        upon setting lat and lon are recalculated
        """
        self.Site.Location.northing = northing
        self.Site.Location.project_location2ll()

    @utm_zone.setter
    def utm_zone(self, utm_zone):
        """
        set UTM zone

        upon setting lat and lon are recalculated
        """
        self.Site.Location.utm_zone = utm_zone
        self.Site.Location.project_location2ll()

    @rotation_angle.setter
    def rotation_angle(self, theta_r):
        """
        set rotation angle in degrees assuming North is 0 measuring clockwise
        positive to East as 90.

        upon setting rotates Z and Tipper
        """

        self._rotation_angle = theta_r
        self._Z.rotate(theta_r)
        self._Tipper.rotate(theta_r)
        self.pt.rotate(theta_r)

        print(("Rotated Z, Tipper, Phase Tensor and Zinvariants by"
               "{0:.3f} degrees".format(self._rotation_angle)))

    @Z.setter
    def Z(self, z_object):
        """
        set z_object

        recalculate phase tensor and invariants, which shouldn't change except
        for strike angle
        """

        self._Z = z_object
        self._Z.compute_resistivity_phase()

    @Tipper.setter
    def Tipper(self, t_object):
        """
        set tipper object

        recalculate tipper angle and magnitude
        """

        self._Tipper = t_object
        if self._Tipper is not None:
            self._Tipper.compute_amp_phase()
            self._Tipper.compute_mag_direction()

    @station.setter
    def station(self, station_name):
        """
        set station name
        """
        self.Site.id = station_name

    # ==========================================================================
    #  read in files
    # ==========================================================================
    def read_mt_file(self, fn, file_type=None):
        """
        Read an MT response file.

        .. note:: Currently only .edi, .xml, and .j files are supported

        :param fn: full path to input file
        :type fn: string

        :param file_type: ['edi' | 'j' | 'xml' | ... ]
                          if None, automatically detects file type by
                          the extension.
        :type file_type: string

        :Example: ::

            >>> import mtpy.core.mt as mt
            >>> mt_obj = mt.MT()
            >>> mt_obj.read_mt_file(r"/home/mt/mt01.xml")

        """

        if file_type is None:
            file_type = os.path.splitext(fn)[1][1:].lower()

        if file_type.lower() == 'edi':
            self._read_edi_file(fn)
        elif file_type.lower() == 'j':
            self._read_j_file(fn)
        elif file_type.lower() == 'xml':
            self._read_xml_file(fn)
        elif file_type.lower() == 'zmm':
            self._read_zmm_file(fn)
        else:
            raise MTError('File type not supported yet')

    def write_mt_file(self, save_dir=None, fn_basename=None, file_type='edi',
                      new_Z_obj=None, new_Tipper_obj=None, longitude_format='LON',
                      latlon_format='dms'
                      ):
        """
        Write an mt file, the supported file types are EDI and XML.

        .. todo:: jtype and Gary Egberts z format

        :param save_dir: full path save directory
        :type save_dir: string

        :param fn_basename: name of file with or without extension
        :type fn_basename: string

        :param file_type: [ 'edi' | 'xml' ]
        :type file_type: string

        :param new_Z_obj: new Z object
        :type new_Z_obj: mtpy.core.z.Z

        :param new_Tipper_obj: new Tipper object
        :type new_Tipper_obj: mtpy.core.z.Tipper

        :param longitude_format:  whether to write longitude as LON or LONG. 
                                  options are 'LON' or 'LONG', default 'LON'
        :type longitude_format:  string
        :param latlon_format:  format of latitude and longitude in output edi,
                               degrees minutes seconds ('dms') or decimal 
                               degrees ('dd')
        :type latlon_format:  string
        
        :returns: full path to file
        :rtype: string

        :Example: ::

            >>> mt_obj.write_mt_file(file_type='xml')

        """

        if save_dir is not None:
            self.save_dir = save_dir

        if fn_basename is not None:
            ext = os.path.splitext(fn_basename)[1][1:].lower()
            fn_basename = os.path.splitext(fn_basename)[0]
            if ext == '':
                fn_basename = '{0}.{1}'.format(fn_basename, file_type.lower())
            elif ext in ['xml', 'edi']:
                fn_basename = '{0}.{1}'.format(fn_basename, ext)
                file_type = ext
            else:
                raise MTError('File type {0} not supported yet.'.format(ext))
        else:
            fn_basename = '{0}.{1}'.format(self.station, file_type)

        fn = os.path.join(self.save_dir, fn_basename)

        if file_type == 'edi':
            fn = self._write_edi_file(fn,
                                      new_Z=new_Z_obj,
                                      new_Tipper=new_Tipper_obj,
                                      longitude_format=longitude_format,
                                      latlon_format=latlon_format)
        elif file_type == 'xml':
            fn = self._write_xml_file(fn,
                                      new_Z=new_Z_obj,
                                      new_Tipper=new_Tipper_obj)

        return fn

    # --> read in edi file
    def _read_edi_file(self, edi_fn):
        """
        read in edi file and set attributes accordingly

        """
        if not os.path.isfile(edi_fn):
            raise MTError('Could not find {0}, check path.'.format(edi_fn))

        self.save_dir = os.path.dirname(edi_fn)

        edi_obj = MTedi.Edi(edi_fn=edi_fn)

        self._edi_get_site(edi_obj)
        
        # get info
        self.Notes = edi_obj.Info
        self._parse_notes()
        
        # get field notes
        self._edi_get_field_notes(edi_obj)
        
        self.Z = edi_obj.Z
        self.Tipper = edi_obj.Tipper
        self.station = edi_obj.station

        # --> make sure things are ordered from high frequency to low
        self._check_freq_order()

    def _edi_get_site(self, edi_obj):
        """
        Get Site information
        """

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
        if hasattr(edi_obj.Header,'enddate'):
            self.Site.end_date = edi_obj.Header.enddate

        self.Site.Location.declination = edi_obj.Header.declination

    def _edi_get_field_notes(self, edi_obj):
        """
        get FieldNotes attributes from edi
        """

        # get information about different sensors
        try:
            for key in list(edi_obj.Define_measurement.meas_hx.__dict__.keys()):
                setattr(self.FieldNotes.Magnetometer_hx,
                        key,
                        edi_obj.Define_measurement.meas_hx.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in list(edi_obj.Define_measurement.meas_hy.__dict__.keys()):
                setattr(self.FieldNotes.Magnetometer_hy,
                        key,
                        edi_obj.Define_measurement.meas_hy.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in list(edi_obj.Define_measurement.meas_hz.__dict__.keys()):
                setattr(self.FieldNotes.Magnetometer_hz,
                        key,
                        edi_obj.Define_measurement.meas_hz.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in list(edi_obj.Define_measurement.meas_ex.__dict__.keys()):
                setattr(self.FieldNotes.Electrode_ex,
                        key,
                        edi_obj.Define_measurement.meas_ex.__dict__[key])
        except AttributeError:
            pass
        try:
            for key in list(edi_obj.Define_measurement.meas_ey.__dict__.keys()):
                setattr(self.FieldNotes.Electrode_ey,
                        key,
                        edi_obj.Define_measurement.meas_ey.__dict__[key])
        except AttributeError:
            pass

        try:
            self.FieldNotes.Magnetometer_hx.id = self.Notes.hx
        except AttributeError:
            pass
        try:
            self.FieldNotes.Magnetometer_hy.id = self.Notes.hy
        except AttributeError:
            pass
        try:
            self.FieldNotes.Magnetometer_hz.id = self.Notes.hz
        except AttributeError:
            pass

        try:
            self.FieldNotes.Magnetometer_hx.type = self.Notes.b_instrument_type
            self.FieldNotes.Magnetometer_hy.type = self.Notes.b_instrument_type
            self.FieldNotes.Magnetometer_hz.type = self.Notes.b_instrument_type
        except AttributeError:
            pass

        try:
            self.FieldNotes.Electrode_ex.type = self.Notes.e_instrument_type
            self.FieldNotes.Electrode_ey.type = self.Notes.e_instrument_type
        except AttributeError:
            pass

        try:
            self.FieldNotes.DataLogger = self.Notes.b_logger_type
        except AttributeError:
            pass
        # keep the edi object around, should be able to deprecate this later
        self._edi_obj = edi_obj

    def _parse_notes(self):
        """
        parse the notes section if there is any information that is useful
        """

        for a_key in list(self.Notes.info_dict.keys()):
            a_value = self.Notes.info_dict[a_key]
            try:
                a_value = float(a_value)
            except (ValueError, TypeError):
                pass
            a_list = a_key.strip().lower().split('.')
            if a_key.find('mtft') == 0 or a_key.find('birrp') == 0:
                a_list = ['processing'] + a_list

            a1 = a_list[0]
            obj_attr = a_list[-1]
            if a1 in ['processing', 'fieldnotes', 'copyright', 'provenance']:
                cl = a_list[0].capitalize()
                if a1 == 'fieldnotes':
                    cl = 'FieldNotes'

                obj = getattr(self, cl)
                count = 1
                while count < len(a_list) - 1:
                    cl_attr = a_list[count]
                    if cl_attr == 'dataquality':
                        cl_attr = 'DataQuality'
                    elif cl_attr == 'datalogger':
                        cl_attr = 'DataLogger'
                    try:
                        obj = getattr(obj, cl_attr)
                    except AttributeError:
                        try:
                            obj = getattr(obj, cl_attr.capitalize())
                        except AttributeError:
                            # print 'Could not get {0}'.format(cl_attr)
                            pass

                    count += 1

                setattr(obj, obj_attr, a_value)
                self.Notes.info_dict.pop(a_key)

    # --> write edi file
    def _write_edi_file(self, new_edi_fn, new_Z=None, new_Tipper=None, 
                        longitude_format='LON', latlon_format='dms'):
        """
        write a new edi file if things have changed.  Note if new_Z or
        new_Tipper are not None, they are not changed in MT object, you
        need to change them manually if you want them to be changed.
        Similarly, the new function name does not change the MT objecte fn
        attribute but does change MT.edi_object.fn attribute.

        :param new_edi_fn: full path to new edi file
        :type new_edi_fn: string

        :param new_Z: new Z object
        :type new_Z: mtpy.core.z.Z

        :param new_Tipper: new Tipper object
        :type new_Tipper: mtpy.core.z.Tipper

        :returns edi_fn: full path to edi file written
        :rtype edi_fn: string
        """

        # get header information, mostly from site
        edi_obj = MTedi.Edi()
        edi_obj.Header = self._edi_set_header()

        # get information
        edi_obj.Info = MTedi.Information()
        edi_obj.Info.info_list = self._edi_set_info_list()

        # get define measurement
        edi_obj.Define_measurement = self._edi_set_define_measurement()

        # get mtsec
        edi_obj.Data_sect = self._edi_set_data_sect()

        # set Z and T data
        if new_Z is not None:
            edi_obj.Z = new_Z
        else:
            edi_obj.Z = self._Z

        if new_Tipper is not None:
            edi_obj.Tipper = new_Tipper
        else:
            edi_obj.Tipper = self._Tipper

        # set rotation angle
        # edi_obj.zrot = self.rotation_angle

        # --> write edi file
        edi_fn = edi_obj.write_edi_file(new_edi_fn=new_edi_fn, 
                                        longitude_format=longitude_format,
                                        latlon_format=latlon_format)

        return edi_fn

    def _edi_set_header(self):
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
        if type(self.Site.survey) is list:
            header.survey = ','.join(self.Site.survey)
        else:
            header.survey = self.Site.survey
        header.units = '[mV/km]/[nT]' 
        header.declination = self.Site.Location.declination
        header.progvers = 'MTpy'
        header.progdate = time.strftime('%Y-%m-%d', time.gmtime())

        return header

    # --> get information list for edi
    def _edi_set_info_list(self):
        """
        get the information for an edi file
        """

        info_list = []
        # write previous information first
        for key in sorted(self.Notes.info_dict.keys()):
            l_key = key.lower()
            l_value = self.Notes.info_dict[key]
            info_list.append(
                '{0} = {1}'.format(
                    l_key.strip(),
                    str(l_value).strip()))

        # get field notes information includes data quality
        for f_key in sorted(self.FieldNotes.__dict__.keys()):
            obj = getattr(self.FieldNotes, f_key)
            for t_key in sorted(obj.__dict__.keys()):
                if t_key in ['_kw_list', '_fmt_list', '_logger']:
                    continue
                l_key = 'fieldnotes.{0}.{1}'.format(f_key.lower(),
                                                    t_key.lower())
                l_value = getattr(obj, t_key)
                if l_value in [None, 'None', 'none']:
                    continue
                info_list.append('{0} = {1}'.format(l_key, l_value))

        # get processing information
        for p_key in sorted(self.Processing.__dict__.keys()):
            if p_key.lower() == 'software':
                for s_key in sorted(self.Processing.Software.__dict__.keys()):
                    if s_key.lower() == 'author':
                        for a_key in sorted(
                                self.Processing.Software.Author.__dict__.keys()):
                            l_key = 'processing.software.author.{0}'.format(
                                a_key)
                            l_value = getattr(self.Processing.Software.Author,
                                              a_key)
                            if l_value in [None, 'None', 'none']:
                                continue
                            info_list.append('{0} = {1}'.format(l_key,
                                                                l_value))
                    else:
                        l_key = 'processing.software.{0}'.format(s_key)
                        l_value = getattr(self.Processing.Software, s_key)
                        if l_value in [None, 'None', 'none']:
                            continue
                        info_list.append('{0} = {1}'.format(l_key,
                                                            l_value))
            elif p_key.lower() == 'remotesite':
                if self.Processing.RemoteSite.id is None:
                    continue
                else:
                    for s_key in sorted(
                            self.Processing.RemoteSite.__dict__.keys()):
                        if s_key == 'Location':
                            for a_key in sorted(
                                    self.Processing.RemoteSite.Location.__dict__.keys()):
                                l_key = 'processing.remote_site.location.{0}'.format(
                                    a_key)
                                l_value = getattr(self.Processing.RemoteSite.Location,
                                                  a_key)
                                if l_value in [None, 'None', 'none']:
                                    continue
                                info_list.append('{0} = {1}'.format(l_key,
                                                                    l_value))
                        else:
                            l_key = 'processing.remote_site.{0}'.format(s_key)
                            l_value = getattr(self.Processing.RemoteSite, s_key)
                            if l_value in [None, 'None', 'none']:
                                continue
                            info_list.append('{0} = {1}'.format(l_key,
                                                                l_value))
            elif p_key.lower() in ['datum', 'coordinate_system',
                                   'sign_convention', 'remote_reference',
                                   'processed_by']:
                l_key = 'processing.{0}'.format(p_key)
                l_value = getattr(self.Processing, p_key)
                if l_value in [None, 'None', 'none']:
                    continue
                info_list.append('{0} = {1}'.format(l_key, l_value))

        # get copyright information
        for c_key in sorted(self.Copyright.__dict__.keys()):
            if c_key.lower() == 'citation':
                if self.Copyright.Citation.author is not None:
                    for p_key in sorted(self.Copyright.Citation.__dict__.keys()):
                        l_key = 'copyright.citation.{0}'.format(p_key.lower())
                        l_value = getattr(self.Copyright.Citation, p_key)
                        if l_value in [None, 'None', 'none']:
                            continue
                        info_list.append('{0} = {1}'.format(l_key, l_value))
                else:
                    continue
            else:
                l_key = 'copyright.{0}'.format(c_key.lower())
                l_value = getattr(self.Copyright, c_key)
                if type(l_value) is list:
                    l_value = ''.join(l_value)
                if l_value in [None, 'None', 'none']:
                    continue
                info_list.append('{0} = {1}'.format(l_key, l_value))

        # get provenance
        for p_key in sorted(self.Provenance.__dict__.keys()):
            if p_key.lower() == 'creator':
                for s_key in list(self.Provenance.Creator.__dict__.keys()):
                    l_key = 'provenance.creator.{0}'.format(s_key)
                    l_value = getattr(self.Provenance.Creator, s_key)
                    if l_value in [None, 'None', 'none']:
                        continue
                    info_list.append('{0} = {1}'.format(l_key, l_value))
            elif p_key.lower() == 'submitter':
                for s_key in list(self.Provenance.Submitter.__dict__.keys()):
                    l_key = 'provenance.submitter.{0}'.format(s_key)
                    l_value = getattr(self.Provenance.Submitter, s_key)
                    if l_value in [None, 'None', 'none']:
                        continue
                    info_list.append('{0} = {1}'.format(l_key, l_value))
            else:
                l_key = 'provenance.{0}'.format(p_key)
                l_value = getattr(self.Provenance, p_key)
                if l_value in [None, 'None', 'none']:
                    continue
                info_list.append('{0} = {1}'.format(l_key, l_value))

        return info_list

    # get edi define measurement
    def _edi_set_define_measurement(self):
        """
        get define measurement block for an edi file
        """

        define_meas = MTedi.DefineMeasurement()
        define_meas.maxchan = 7
        define_meas.refelev = self.elev
        define_meas.reflat = self.lat
        define_meas.reflon = self.lon
        define_meas.reftype = self.Site.Location.coordinate_system
        define_meas.units = self.Site.Location.elev_units

        define_meas.meas_ex = MTedi.EMeasurement()
        for key in define_meas.meas_ex._kw_list:
            setattr(define_meas.meas_ex,
                    key,
                    getattr(self.FieldNotes.Electrode_ex, key))

        define_meas.meas_ey = MTedi.EMeasurement()
        for key in define_meas.meas_ey._kw_list:
            setattr(define_meas.meas_ey,
                    key,
                    getattr(self.FieldNotes.Electrode_ey, key))

        define_meas.meas_hx = MTedi.HMeasurement()
        for key in define_meas.meas_hx._kw_list:
            setattr(define_meas.meas_hx,
                    key,
                    getattr(self.FieldNotes.Magnetometer_hx, key))

        define_meas.meas_hy = MTedi.HMeasurement()
        for key in define_meas.meas_hy._kw_list:
            setattr(define_meas.meas_hy,
                    key,
                    getattr(self.FieldNotes.Magnetometer_hy, key))

        if self.Tipper is not None:
            if np.all(self.Tipper.tipper == 0) == False:
                define_meas.meas_hz = MTedi.HMeasurement()
                for key in define_meas.meas_hz._kw_list:
                    setattr(define_meas.meas_hz,
                            key,
                            getattr(self.FieldNotes.Magnetometer_hz, key))

        return define_meas

    def _edi_set_data_sect(self):
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
        sect.ex = self.FieldNotes.Electrode_ex.acqchan
        sect.ey = self.FieldNotes.Electrode_ey.acqchan
        sect.hx = self.FieldNotes.Magnetometer_hx.acqchan
        sect.hy = self.FieldNotes.Magnetometer_hy.acqchan
        if np.all(self.Tipper.tipper == 0) == False:
            sect.hz = self.FieldNotes.Magnetometer_hz.acqchan

        return sect

    # --> check the order of frequencies
    def _check_freq_order(self):
        """
        check to make sure the Z and Tipper arrays are ordered such that
        the first index corresponds to the highest frequency and the last
        index corresponds to the lowest freqeuncy.

        """

        if self.Z.freq[0] < self.Z.freq[1]:
            print('Flipping arrays to be ordered from short period to long')
            self.Z.z = self.Z.z.copy()[::-1]
            self.Z.z_err = self.Z.z_err.copy()[::-1]
            self.Z.freq = self.Z.freq.copy()[::-1]

        if self.Tipper.tipper is not None:
            if self.Tipper.freq[0] < self.Tipper.freq[1]:
                self.Tipper.tipper = self.Tipper.tipper.copy()[::-1]
                self.Tipper.tipper_err = self.Tipper.tipper_err.copy()[::-1]
                self.Tipper.freq = self.Tipper.freq.copy()[::-1]

    def _read_j_file(self, j_fn):
        """
        read j file
        """

        j_obj = MTj.JFile(j_fn)

        self.save_dir = os.path.dirname(j_fn)
        self.station = os.path.splitext(os.path.basename(j_fn))[0]

        self.Z = j_obj.Z
        self.Tipper = j_obj.Tipper

        self._check_freq_order()

        self.Site.Location.latitude = j_obj.metadata_dict['latitude']
        self.Site.Location.longitude = j_obj.metadata_dict['longitude']
        self.Site.Location.elevation = j_obj.metadata_dict['elevation']

        self.Notes.info_dict = j_obj.header_dict
        
    def _read_zmm_file(self, zmm_fn):
        """
        read zmm file
        """
        if not isinstance(zmm_fn, Path):
            zmm_fn = Path(zmm_fn)
            
        zmm_obj = MTzmm.ZMM(zmm_fn)
        zmm_obj.read_zmm_file()
        
        self.save_dir = zmm_fn.parent
        self.Z = zmm_obj.Z
        self.Tipper = zmm_obj.Tipper
        
        # set location
        self.Site.Location.latitude = zmm_obj.lat
        self.Site.Location.longitude = zmm_obj.lon
        self.Site.Location.declination = zmm_obj.declination
        
        # set station name
        self.Site.id = zmm_obj.station
        

    def _read_xml_file(self, xml_fn):
        """
        read xml file
        """

        if not os.path.isfile(xml_fn):
            raise MTError('Could not find {0}, check path.'.format(xml_fn))

        self.save_dir = os.path.dirname(xml_fn)

        xml_obj = MTxml.MT_XML()
        xml_obj.read_xml_file(xml_fn)

        self.Z = xml_obj.Z
        self.Tipper = xml_obj.Tipper

        # check order
        self._check_freq_order()

        # get information and fill attributes
        self._xml_get_site(xml_obj)
        self._xml_get_notes(xml_obj)
        self._xml_get_field_notes(xml_obj)
        self._xml_get_copyright(xml_obj)
        self._xml_get_provenance(xml_obj)
        self._xml_get_processing(xml_obj)

    def _xml_get_site(self, xml_obj):
        """
        get Site information from xml Site
        """
        # get information
        for s_attr in list(xml_obj.Site.__dict__.keys()):
            if s_attr in ['_name', '_attr', '_value']:
                continue
            x_obj = getattr(xml_obj.Site, s_attr)
            name = x_obj.name.lower()
            if name == 'acquiredby':
                name = 'acquired_by'
            elif name == 'end':
                name = 'end_date'
            elif name == 'start':
                name = 'start_date'
            elif name == 'runlist':
                name = 'run_list'
            elif name == 'yearcollected':
                name = 'year_collected'
            elif name == 'datecollected':
                name = 'date_collected'

            value = x_obj.value
            if name == 'location':
                for l_attr in list(xml_obj.Site.Location.__dict__.keys()):
                    if l_attr in ['_name', '_attr', '_value']:
                        continue
                    l_obj = getattr(xml_obj.Site.Location, l_attr)
                    name = l_obj.name.lower()
                    value = l_obj.value
                    if name == 'elevation':
                        units = l_obj.attr['units']
                        self.Site.Location.elev_units = units

                    elif name == 'declination':
                        units = l_obj.attr['epoch']
                        self.Site.Location.declination_epoch = units
                    setattr(self.Site.Location, name, value)
            else:
                setattr(self.Site, name, value)

    def _xml_get_notes(self, xml_obj):
        """
        get notes and set in a dictionary that is similar to edi notes
        """

        self.Notes = MTedi.Information()
        self.Notes.info_dict = {}
        self.Notes.info_list = []

        self.Notes.info_dict['notes'] = str(xml_obj.Notes.value)
        self.Notes.info_list = ['notes = {0}'.format(str(xml_obj.Notes.value))]

    def _xml_get_field_notes(self, xml_obj):
        """
        get field notes information
        """

        for f_attr in list(xml_obj.FieldNotes.__dict__.keys()):
            if f_attr.lower() == 'instrument':
                for i_attr in list(xml_obj.FieldNotes.Instrument.__dict__.keys()):
                    if i_attr in ['_name', '_attr', '_value']:
                        continue
                    i_obj = getattr(xml_obj.FieldNotes.Instrument, i_attr)
                    name = i_obj.name.lower()
                    value = i_obj.value
                    setattr(self.FieldNotes.DataLogger, name, value)

            elif 'dipole' in f_attr.lower():
                xml_d_obj = getattr(xml_obj.FieldNotes, f_attr)
                azm = 0.0
                length = 0.0
                try:
                    comp = xml_d_obj.attr['name'].lower()
                except KeyError:
                    comp = 'ex'
                try:
                    t = xml_d_obj.attr['type'].lower()
                    if comp == 'ex':
                        setattr(self.FieldNotes.Electrode_ex, 'type', t)
                    elif comp == 'ey':
                        setattr(self.FieldNotes.Electrode_ey, 'type', t)
                except KeyError:
                    pass

                for e_attr in list(xml_d_obj.__dict__.keys()):
                    if e_attr in ['_name', '_attr', '_value']:
                        continue
                    e_obj = getattr(xml_d_obj, e_attr)
                    name = e_obj.name.lower()
                    value = e_obj.value

                    if name == 'azimuth':
                        azm = float(value)

                    if name == 'length':
                        length = float(value)

                    if comp == 'ex':
                        setattr(self.FieldNotes.Electrode_ex, name, value)
                    elif comp == 'ey':
                        setattr(self.FieldNotes.Electrode_ey, name, value)

                # need to set x, y, x2, y2
                x = 0
                y = 0
                x2 = length * np.cos(np.deg2rad(azm))
                y2 = length * np.sin(np.deg2rad(azm))

                for name, value in zip(['x', 'y', 'x2', 'y2'], [x, y, x2, y2]):
                    if comp == 'ex':
                        setattr(self.FieldNotes.Electrode_ex, name, value)
                    elif comp == 'ey':
                        setattr(self.FieldNotes.Electrode_ey, name, value)

            elif 'magnetometer' in f_attr.lower():
                xml_d_obj = getattr(xml_obj.FieldNotes, f_attr)
                try:
                    comp = xml_d_obj.attr['name'].lower()
                except KeyError:
                    try:
                        comp = xml_d_obj.attr['type'].lower()
                    except KeyError:
                        comp = 'hx'

                try:
                    t = xml_d_obj.attr['type'].lower()
                    if comp == 'hx':
                        setattr(self.FieldNotes.Magnetometer_hx, 'type', t)
                    elif comp == 'hy':
                        setattr(self.FieldNotes.Magnetometer_hy, 'type', t)
                    elif comp == 'hz':
                        setattr(self.FieldNotes.Magnetometer_hz, 'type', t)
                    elif comp == 'fluxgate':
                        setattr(self.FieldNotes.Magnetometer_hx, 'type', t)
                        setattr(self.FieldNotes.Magnetometer_hy, 'type', t)
                        setattr(self.FieldNotes.Magnetometer_hz, 'type', t)
                    else:
                        pass
                except KeyError:
                    pass

                for m_attr in list(xml_d_obj.__dict__.keys()):
                    if m_attr in ['_name', '_attr', '_value']:
                        continue
                    m_obj = getattr(xml_obj.FieldNotes.Magnetometer, m_attr)
                    name = m_obj.name.lower()
                    value = m_obj.value

                    if comp == 'hx':
                        setattr(self.FieldNotes.Magnetometer_hx, name, value)
                    elif comp == 'hy':
                        setattr(self.FieldNotes.Magnetometer_hy, name, value)
                    elif comp == 'hz':
                        setattr(self.FieldNotes.Magnetometer_hz, name, value)
                    elif comp == 'fluxgate':
                        setattr(self.FieldNotes.Magnetometer_hx, name, value)
                        setattr(self.FieldNotes.Magnetometer_hy, name, value)
                        setattr(self.FieldNotes.Magnetometer_hz, name, value)

            elif 'dataquality' in f_attr.lower():
                obj = getattr(xml_obj.FieldNotes, f_attr)
                for d_attr in list(obj.__dict__.keys()):
                    if d_attr in ['_name', '_attr', '_value']:
                        continue
                    d_obj = getattr(obj, d_attr)
                    name = d_obj.name.lower()
                    if name == 'goodfromperiod':
                        name = 'good_from_period'
                    elif name == 'goodtoperiod':
                        name = 'good_to_period'
                    elif name == 'comments':
                        setattr(self.FieldNotes.DataQuality,
                                'author',
                                d_obj.attr['author'])
                    if name == 'comments' and \
                            f_attr.lower() == 'dataqualitywarnings':
                        name = 'warnings_' + name
                    value = d_obj.value

                    setattr(self.FieldNotes.DataQuality, name, value)

    def _xml_get_copyright(self, xml_obj):
        """
        get copyright information
        """

        for f_attr in list(xml_obj.Copyright.__dict__.keys()):
            if f_attr in ['_name', '_attr', '_value']:
                continue
            if f_attr.lower() == 'citation':
                for i_attr in list(xml_obj.Copyright.Citation.__dict__.keys()):
                    if i_attr in ['_name', '_attr', '_value']:
                        continue
                    i_obj = getattr(xml_obj.Copyright.Citation, i_attr)
                    name = i_obj.name.lower()
                    value = i_obj.value
                    setattr(self.Copyright.Citation, name, value)
            else:
                obj = getattr(xml_obj.Copyright, f_attr)
                name = obj.name.lower()
                if name == 'releasestatus':
                    name = 'release_status'
                elif name == 'conditionsofuse':
                    name = 'conditions_of_use'
                elif name == 'additionalinfo':
                    name = 'additional_info'
                value = obj.value.replace('\n', '')

                setattr(self.Copyright, name, value)

    def _xml_get_provenance(self, xml_obj):
        """
        get provenance infor
        """
        for f_attr in list(xml_obj.Provenance.__dict__.keys()):
            if f_attr in ['_name', '_attr', '_value']:
                continue
            if f_attr.lower() in ['creator', 'submitter']:
                obj = getattr(xml_obj.Provenance, f_attr)
                s_obj = getattr(self.Provenance, f_attr)
                for i_attr in list(obj.__dict__.keys()):
                    if i_attr in ['_name', '_attr', '_value']:
                        continue
                    i_obj = getattr(obj, i_attr)
                    name = i_obj.name.lower()
                    value = i_obj.value
                    setattr(s_obj, name, value)
            else:
                obj = getattr(xml_obj.Provenance, f_attr)
                name = obj.name.lower()
                if name == 'creationtime':
                    name = 'creation_time'
                elif name == 'creatingapplication':
                    name = 'creating_application'
                value = obj.value

                setattr(self.Provenance, name, value)

    def _xml_get_processing(self, xml_obj):
        """
        get processing info
        """

        for f_attr in list(xml_obj.ProcessingInfo.__dict__.keys()):
            if f_attr in ['_name', '_attr', '_value']:
                continue
            if 'software' in f_attr.lower():
                obj = getattr(xml_obj.ProcessingInfo, f_attr)
                for i_attr in list(obj.__dict__.keys()):
                    if i_attr in ['_name', '_attr', '_value']:
                        continue
                    i_obj = getattr(obj, i_attr)
                    name = i_obj.name.lower()
                    if name == 'lastmod':
                        name = 'last_modified'
                    value = i_obj.value
                    if name == 'author':
                        value = Person()
                        value.name = i_obj.value
                    setattr(self.Processing.Software, name, value)
            elif 'remoteinfo' in f_attr.lower():
                obj = getattr(xml_obj.ProcessingInfo, f_attr)
                for i_attr in list(obj.__dict__.keys()):
                    if i_attr in ['_name', '_attr', '_value']:
                        continue
                    if i_attr.lower() == 'location':
                        loc_obj = getattr(obj, i_attr)
                        
                        for l_attr in list(loc_obj.__dict__.keys()):
                            if l_attr in ['_name', '_attr', '_value']:
                                continue
                            l_obj = getattr(loc_obj, l_attr)
                            name = l_obj.name.lower()
                            value = l_obj.value
                            setattr(self.Processing.RemoteSite.Location,
                                    name, value)

                    else:
                        i_obj = getattr(obj, i_attr)
                        name = i_obj.name.lower()
                        if name == 'yearcollected':
                            name = 'year_collected'
                            value = i_obj.value
                        setattr(self.Processing.RemoteSite, name, value)
            else:
                obj = getattr(xml_obj.ProcessingInfo, f_attr)
                name = obj.name.lower()
                if name == 'signconvention':
                    name = 'sign_convention'
                elif name == 'remoteref':
                    name = 'remote_ref'
                elif name == 'processedby':
                    name = 'processed_by'
                elif name == 'processingtag':
                    name = 'processing_tag'
                value = obj.value

                setattr(self.Processing, name, value)

    def _write_xml_file(self, xml_fn, new_Z=None, new_Tipper=None):
        """
        Write a xml file.
        """

        if new_Z is not None:
            self.Z = new_Z
        if new_Tipper is not None:
            self.Tipper = new_Tipper

        xml_obj = MTxml.MT_XML()
        xml_obj.Attachment.Filename.value = os.path.basename(self.fn)
        xml_obj.PrimaryData.Filename.value = os.path.basename(self.fn)[:-4]+'.png'

        xml_obj.ProductId.value = '{0}.{1}'.format(self.station.upper(), 
                                                   self.Site.year_collected)
        
        xml_obj.Z = self.Z
        xml_obj.Tipper = self.Tipper

        xml_obj = self._xml_set_provenance(xml_obj)
        xml_obj = self._xml_set_copyright(xml_obj)
        xml_obj = self._xml_set_site(xml_obj)
        xml_obj = self._xml_set_field_notes(xml_obj)
        xml_obj = self._xml_set_processing(xml_obj)
        xml_obj = self._xml_set_site_layout(xml_obj)
        
        xml_obj.write_xml_file(xml_fn)

    def _xml_set_site(self, xml_obj):
        """
        set the Site attributes in the xml object
        """
        xml_obj.Site.Project.value = self.Site.project
        if type(self.Site.survey) is list:
            xml_obj.Site.Survey.value = ','.join(self.Site.survey)
        else:
            xml_obj.Site.Survey.value = self.Site.survey
        xml_obj.Site.Id.value = self.Site.id
        xml_obj.Site.AcquiredBy.value = self.Site.acquired_by
        xml_obj.Site.Start.value = self.Site.start_date
        xml_obj.Site.End.value = self.Site.end_date
        xml_obj.Site.RunList.value = self.Site.run_list
        xml_obj.Site.Orientation.value = 'geomagnetic'
        try:
            xml_obj.Site.Orientation.attr = {'angle_to_geographic_north':\
                                             '{0:.2f}'.format(self.Site.Location.declination)}
        except ValueError:
            xml_obj.Site.Orientation.attr = {'angle_to_geographic_north': '0.00'}
        
        xml_obj.Site.Location.Latitude.value = self.lat
        xml_obj.Site.Location.Longitude.value = self.lon
        xml_obj.Site.Location.Elevation.value = self.elev
        xml_obj.Site.Location.Elevation.attr = {
            'units': self.Site.Location.elev_units}
        xml_obj.Site.Location.Declination.value = self.Site.Location.declination
        xml_obj.Site.Location.Declination.attr = {
            'epoch': self.Site.Location.declination_epoch}

        return xml_obj

    def _xml_set_field_notes(self, xml_obj):
        """
        Set the FieldNotes attributes of the xml object
        """

        xml_obj.FieldNotes.Instrument.Type.value = self.FieldNotes.DataLogger.type
        xml_obj.FieldNotes.Instrument.Id.value = self.FieldNotes.DataLogger.id
        xml_obj.FieldNotes.Instrument.Manufacturer.value = self.FieldNotes.DataLogger.manufacturer

        # EX
        xml_obj.FieldNotes.Dipole.Type.value = self.FieldNotes.Electrode_ex.type
        xml_obj.FieldNotes.Dipole.Id.value = self.FieldNotes.Electrode_ex.id
        xml_obj.FieldNotes.Dipole.Manufacturer.value = self.FieldNotes.Electrode_ex.manufacturer
        xml_obj.FieldNotes.Dipole.attr = {'name': 'EX'}
        
        length = np.sqrt((self.FieldNotes.Electrode_ex.x2 - self.FieldNotes.Electrode_ex.x) ** 2 +
                         (self.FieldNotes.Electrode_ex.y2 - self.FieldNotes.Electrode_ex.y) ** 2)
        xml_obj.FieldNotes.Dipole.Length.value = length
        try:
            azm = np.arctan((self.FieldNotes.Electrode_ex.y2 - self.FieldNotes.Electrode_ex.y) /
                            (self.FieldNotes.Electrode_ex.x2 - self.FieldNotes.Electrode_ex.x))
        except ZeroDivisionError:
            azm = 0.0
        xml_obj.FieldNotes.Dipole.Azimuth.value = np.degrees(azm)
        xml_obj.FieldNotes.Dipole.Channel.value = self.FieldNotes.Electrode_ex.acqchan

        # EY
        xml_obj.FieldNotes.Dipole_00.Type.value = self.FieldNotes.Electrode_ey.type
        xml_obj.FieldNotes.Dipole_00.Id.value = self.FieldNotes.Electrode_ey.id
        xml_obj.FieldNotes.Dipole_00.Manufacturer.value = self.FieldNotes.Electrode_ey.manufacturer
        xml_obj.FieldNotes.Dipole_00.attr = {'name': 'EY'}
        length = np.sqrt((self.FieldNotes.Electrode_ey.x2 - self.FieldNotes.Electrode_ey.x) ** 2 +
                         (self.FieldNotes.Electrode_ey.y2 - self.FieldNotes.Electrode_ey.y) ** 2)
        xml_obj.FieldNotes.Dipole_00.Length.value = length
        try:
            azm = np.arctan((self.FieldNotes.Electrode_ey.y2 - self.FieldNotes.Electrode_ey.y) /
                            (self.FieldNotes.Electrode_ey.x2 - self.FieldNotes.Electrode_ey.x))
        except ZeroDivisionError:
            azm = 90.0
        xml_obj.FieldNotes.Dipole_00.Azimuth.value = np.degrees(min([np.pi/2, azm]))
        xml_obj.FieldNotes.Dipole_00.Channel.value = self.FieldNotes.Electrode_ey.acqchan

        # HX
        xml_obj.FieldNotes.Magnetometer.Type.value = self.FieldNotes.Magnetometer_hx.type
        xml_obj.FieldNotes.Magnetometer.Id.value = self.FieldNotes.Magnetometer_hx.id
        xml_obj.FieldNotes.Magnetometer.Manufacturer.value = self.FieldNotes.Magnetometer_hx.manufacturer
        xml_obj.FieldNotes.Magnetometer.attr = {'name': 'HX'}
        xml_obj.FieldNotes.Magnetometer.Azimuth.value = self.FieldNotes.Magnetometer_hx.azm
        xml_obj.FieldNotes.Magnetometer.Channel.value = self.FieldNotes.Magnetometer_hx.acqchan

        # HY
        xml_obj.FieldNotes.Magnetometer_00.Type.value = self.FieldNotes.Magnetometer_hy.type
        xml_obj.FieldNotes.Magnetometer_00.Id.value = self.FieldNotes.Magnetometer_hy.id
        xml_obj.FieldNotes.Magnetometer_00.Manufacturer.value = self.FieldNotes.Magnetometer_hy.manufacturer
        xml_obj.FieldNotes.Magnetometer_00.attr = {'name': 'HY'}
        xml_obj.FieldNotes.Magnetometer_00.Azimuth.value = self.FieldNotes.Magnetometer_hy.azm
        xml_obj.FieldNotes.Magnetometer_00.Channel.value = self.FieldNotes.Magnetometer_hy.acqchan

        # HZ
        xml_obj.FieldNotes.Magnetometer_01.Type.value = self.FieldNotes.Magnetometer_hz.type
        xml_obj.FieldNotes.Magnetometer_01.Id.value = self.FieldNotes.Magnetometer_hz.id
        xml_obj.FieldNotes.Magnetometer_01.Manufacturer.value = self.FieldNotes.Magnetometer_hz.manufacturer
        xml_obj.FieldNotes.Magnetometer_01.attr = {'name': 'HZ'}
        xml_obj.FieldNotes.Magnetometer_01.Azimuth.value = self.FieldNotes.Magnetometer_hz.azm
        xml_obj.FieldNotes.Magnetometer_01.Channel.value = self.FieldNotes.Magnetometer_hz.acqchan

        # Data Quality Notes
        xml_obj.FieldNotes.DataQualityNotes.Rating.value = self.FieldNotes.DataQuality.rating
        xml_obj.FieldNotes.DataQualityNotes.GoodFromPeriod.value = self.FieldNotes.DataQuality.good_from_period
        xml_obj.FieldNotes.DataQualityNotes.GoodToPeriod.value = self.FieldNotes.DataQuality.good_to_period
        xml_obj.FieldNotes.DataQualityNotes.Comments.value = self.FieldNotes.DataQuality.comments
        xml_obj.FieldNotes.DataQualityNotes.Comments.attr = {
            'author': self.FieldNotes.DataQuality.author}
        # Data Quality Warnings
        xml_obj.FieldNotes.DataQualityWarnings.Flag.value = self.FieldNotes.DataQuality.warnings_flag
        xml_obj.FieldNotes.DataQualityWarnings.Comments.value = self.FieldNotes.DataQuality.warnings_comments
        xml_obj.FieldNotes.DataQualityWarnings.Comments.attr = {
            'author': self.FieldNotes.DataQuality.author}

        return xml_obj

    def _xml_set_processing(self, xml_obj):
        """
        Set the Processing attributes of the xml object
        """

        xml_obj.ProcessingInfo.ProcessedBy.value = self.Processing.processed_by
        xml_obj.ProcessingInfo.ProcessingSoftware.Name.value = self.Processing.Software.name
        xml_obj.ProcessingInfo.ProcessingSoftware.Author.value = self.Processing.Software.Author.name
        xml_obj.ProcessingInfo.ProcessingSoftware.Version.value = self.Processing.Software.version

        # TODO: Need to find a way to put in processing parameters.

        xml_obj.ProcessingInfo.SignConvention.value = self.Processing.sign_convention

        xml_obj.ProcessingInfo.RemoteRef.value = self.Processing.remote_reference
        xml_obj.ProcessingInfo.RemoteInfo.Project.value = self.Processing.RemoteSite.project
        xml_obj.ProcessingInfo.RemoteInfo.Survey.value = self.Processing.RemoteSite.survey
        xml_obj.ProcessingInfo.RemoteInfo.ID.value = self.Processing.RemoteSite.id
        xml_obj.ProcessingInfo.RemoteInfo.YearCollected.value = self.Processing.RemoteSite.year_collected
        xml_obj.ProcessingInfo.RemoteInfo.AcquiredBy.value = self.Processing.RemoteSite.acquired_by
        xml_obj.ProcessingInfo.RemoteInfo.Location.Latitude.value = self.Processing.RemoteSite.Location.latitude
        xml_obj.ProcessingInfo.RemoteInfo.Location.Longitude.value = self.Processing.RemoteSite.Location.longitude
        xml_obj.ProcessingInfo.RemoteInfo.Location.Elevation.value = self.Processing.RemoteSite.Location.elevation
        xml_obj.ProcessingInfo.RemoteInfo.Location.Elevation.attr = {
            'units': self.Processing.RemoteSite.Location.elev_units}
        xml_obj.ProcessingInfo.RemoteInfo.Location.attr = {
            'datum': self.Processing.RemoteSite.Location.datum}

        return xml_obj

    def _xml_set_provenance(self, xml_obj):
        """
        Set the Provenance attributes of the xml object
        """

        xml_obj.Provenance.CreatingApplication.value = 'MTpy 0.1.0'

        xml_obj.Provenance.Submitter.Name.value = self.Provenance.Submitter.name
        xml_obj.Provenance.Submitter.Email.value = self.Provenance.Submitter.email
        xml_obj.Provenance.Submitter.Org.value = self.Provenance.Submitter.organization
        xml_obj.Provenance.Submitter.OrgURL.value = self.Provenance.Submitter.organization_url

        xml_obj.Provenance.Creator.Name.value = self.Provenance.Creator.name
        xml_obj.Provenance.Creator.Email.value = self.Provenance.Creator.email
        xml_obj.Provenance.Creator.Org.value = self.Provenance.Creator.organization
        xml_obj.Provenance.Creator.OrgURL.value = self.Provenance.Creator.organization_url

        return xml_obj

    def _xml_set_copyright(self, xml_obj):
        """
        Set the Copyright attributes of the xml object
        """

        xml_obj.Copyright.Citation.Authors.value = self.Copyright.Citation.author
        xml_obj.Copyright.Citation.Title.value = self.Copyright.Citation.title
        xml_obj.Copyright.Citation.Year.value = self.Copyright.Citation.year
        xml_obj.Copyright.Citation.Journal.value = self.Copyright.Citation.journal
        xml_obj.Copyright.Citation.Volume.value = self.Copyright.Citation.volume
        xml_obj.Copyright.Citation.DOI.value = self.Copyright.Citation.doi
        if type(self.Copyright.conditions_of_use) is list:
            xml_obj.Copyright.ConditionsOfUse.value = ''.join(self.Copyright.conditions_of_use)
        else:
            xml_obj.Copyright.ConditionsOfUse.value = self.Copyright.conditions_of_use
        xml_obj.Copyright.ReleaseStatus.value = self.Copyright.release_status
        xml_obj.Copyright.AdditionalInfo.value = self.Copyright.additional_info

        return xml_obj

    def _xml_set_site_layout(self, xml_obj):
        """
        set the site layout from define measurement
        """
        
        xml_obj.SiteLayout.InputChannels.Magnetic_hx.attr = {'name':"Hx",
                                                             'orientation':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hx.azm), 
                                                             'x':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hx.x),
                                                             'y':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hx.y),
                                                             'z':'0.00'}
        xml_obj.SiteLayout.InputChannels.Magnetic_hy.attr = {'name':"Hy",
                                                             'orientation':'{0:.2f}'.format(max([90, self.FieldNotes.Magnetometer_hy.azm])), 
                                                             'x':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hy.x),
                                                             'y':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hy.y),
                                                             'z':'0.00'}
        xml_obj.SiteLayout.OutputChannels.Magnetic_hz.attr = {'name':"Hz",
                                                             'orientation':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hz.azm), 
                                                             'x':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hz.x),
                                                             'y':'{0:.2f}'.format(self.FieldNotes.Magnetometer_hz.y),
                                                             'z':'0.00'}
        try:
            azm = np.arctan((self.FieldNotes.Electrode_ex.y2 - self.FieldNotes.Electrode_ex.y) /
                            (self.FieldNotes.Electrode_ex.x2 - self.FieldNotes.Electrode_ex.x))
        except ZeroDivisionError:
            azm = 90.0
        xml_obj.SiteLayout.OutputChannels.Electric_ex.attr = {'name':"Ex",
                                                             'orientation':'{0:.2f}'.format(azm), 
                                                             'x':'{0:.2f}'.format(self.FieldNotes.Electrode_ex.x),
                                                             'y':'{0:.2f}'.format(self.FieldNotes.Electrode_ex.y),
                                                             'z':'0.00',
                                                             'x2':'{0:.2f}'.format(self.FieldNotes.Electrode_ex.x2),
                                                             'y2':'{0:.2f}'.format(self.FieldNotes.Electrode_ex.y2),
                                                             'z2':'0.00'}
        try:
            azm = np.arctan((self.FieldNotes.Electrode_ey.y2 - self.FieldNotes.Electrode_ey.y) /
                            (self.FieldNotes.Electrode_ey.x2 - self.FieldNotes.Electrode_ey.x))
        except ZeroDivisionError:
            azm = 90.0
        xml_obj.SiteLayout.OutputChannels.Electric_ey.attr = {'name':"Ey",
                                                             'orientation':'{0:.2f}'.format(azm), 
                                                             'x':'{0:.2f}'.format(self.FieldNotes.Electrode_ey.x),
                                                             'y':'{0:.2f}'.format(self.FieldNotes.Electrode_ey.y),
                                                             'z':'0.00',
                                                             'x2':'{0:.2f}'.format(self.FieldNotes.Electrode_ey.x2),
                                                             'y2':'{0:.2f}'.format(self.FieldNotes.Electrode_ey.y2),
                                                             'z2':'0.00'}
        return xml_obj
    
    def read_cfg_file(self, cfg_fn):
        """
        Read in a configuration file and populate attributes accordingly.

        The configuration file should be in the form:
            | Site.Location.latitude = 46.5
            | Site.Location.longitude = 122.7
            | Site.Location.datum = 'WGS84'
            |
            | Processing.Software.name = BIRRP
            | Processing.Software.version = 5.2.1
            |
            | Provenance.Creator.name = L. Cagniard
            | Provenance.Submitter.name = I. Larionov


        :param cfg_fn: full path to configuration file
        :type cfg_fn: string

        .. note:: The best way to make a configuration file would be to save
                  a configuration file first from MT, then filling in the
                  fields.

        :Make configuration file: ::

            >>> import mtpy.core.mt as mt
            >>> mt_obj = mt.MT()
            >>> mt_obj.write_cfg_file(r"/mt/generic_config.cfg")

        :Read in configuration file: ::

            >>> import mtpy.core.mt as mt
            >>> mt_obj = mt.MT()
            >>> mt_obj.read_cfg_file(r"/home/mt/survey_config.cfg")
        """

        with open(cfg_fn, 'r') as fid:
            cfg_lines = fid.readlines()

        for line in cfg_lines:
            if line[0] == '#':
                continue
            line_list = line.strip().split('=')
            if len(line) < 2:
                continue
            else:
                name_list = line_list[0].strip().split('.')
                cl_name = name_list[0]
                if cl_name.lower() == 'fieldnotes':
                    cl_name = 'FieldNotes'
                else:
                    cl_name = cl_name.capitalize()

                cl_obj = getattr(self, cl_name)
                cl_attr = name_list[-1]
                cl_value = line_list[1].strip()
                if cl_value in ['None', 'none']:
                    cl_value = None
                try:
                    cl_value = float(cl_value)
                except ValueError:
                    if cl_value.find(',') >= 0:
                        cl_value = cl_value.replace('[', '').replace(']', '')
                        cl_value = cl_value.strip().split(',')
                        try:
                            cl_value = [float(vv) for vv in cl_value]
                        except ValueError:
                            pass
                except TypeError:
                    cl_value = None

                count = 1
                while count < len(name_list) - 1:
                    cl_name = name_list[count].lower()
                    if cl_name == 'dataquality':
                        cl_name = 'DataQuality'
                    elif cl_name == 'datalogger':
                        cl_name = 'DataLogger'
                    elif cl_name == 'remotesite':
                        cl_name = 'RemoteSite'
                    try:
                        cl_obj = getattr(cl_obj, cl_name)
                    except AttributeError:
                        try:
                            cl_obj = getattr(cl_obj, cl_name.capitalize())
                            cl_name = cl_name.capitalize()
                        except AttributeError:
                            print('Could not get {0}'.format(cl_name))

                    count += 1
                setattr(cl_obj, cl_attr, cl_value)

    def write_cfg_file(self, cfg_fn):
        """
        Write a configuration file for the MT sections

        :param cfg_fn: full path to configuration file to write to
        :type cfg_fn: string

        :return cfg_fn: full path to configuration file
        :rtype cfg_fn: string

        :Write configuration file: ::

            >>> import mtpy.core.mt as mt
            >>> mt_obj = mt.MT()
            >>> mt_obj.read_mt_file(r"/home/mt/edi_files/mt01.edi")
            >>> mt_obj.write_cfg_file(r"/home/mt/survey_config.cfg")

        """

        cfg_lines = []
        for obj_name in sorted(['Site', 'FieldNotes', 'Provenance',
                                'Processing', 'Copyright']):
            obj = getattr(self, obj_name)
            l_key = obj_name
            for obj_key in sorted(obj.__dict__.keys()):
                obj_attr = getattr(obj, obj_key)
                l_key = '{0}.{1}'.format(obj_name, obj_key)

                if not isinstance(obj_attr, (str, float, int, list)) and \
                        obj_attr is not None:
                    for a_key in sorted(obj_attr.__dict__.keys()):
                        if a_key in ['_kw_list', '_fmt_list', '_logger']:
                            continue
                        obj_attr_01 = getattr(obj_attr, a_key)
                        l_key = '{0}.{1}.{2}'.format(obj_name, obj_key, a_key)
                        if not isinstance(obj_attr_01, (str, float, int,
                                                        list, np.float64)) and \
                                obj_attr_01 is not None:
                            for b_key in sorted(obj_attr_01.__dict__.keys()):
                                obj_attr_02 = getattr(obj_attr_01, b_key)
                                l_key = '{0}.{1}.{2}.{3}'.format(obj_name,
                                                                 obj_key,
                                                                 a_key,
                                                                 b_key)

                                cfg_lines.append('{0} = {1}'.format(l_key,
                                                                    obj_attr_02))
                        else:
                            cfg_lines.append('{0} = {1}'.format(l_key,
                                                                obj_attr_01))
                else:
                    cfg_lines.append('{0} = {1}'.format(l_key, obj_attr))

            cfg_lines.append('')

        with open(cfg_fn, 'w') as fid:
            fid.write('\n'.join(cfg_lines))

        print('--> Wrote MT configuration file to {0}'.format(cfg_fn))

        return cfg_fn

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
        dummy_z_obj = MTz.copy.deepcopy(self.Z)
        D, new_z_object = MTdistortion.remove_distortion(z_object=dummy_z_obj,
                                                         num_freq=num_freq)

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

        s_array, new_z = self.Z.remove_ss(reduce_res_factor_x=ss_x,
                                          reduce_res_factor_y=ss_y)

        new_z_obj = MTz.Z(z_array=new_z,
                          z_err_array=self.Z.z_err.copy(),
                          freq=self.Z.freq.copy())

        return new_z_obj

    def interpolate(self, new_freq_array, interp_type='slinear', bounds_error=True, period_buffer=None):
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
        # if the interpolation module has not been loaded return
        if interp_import is False:
            raise ImportError('could not interpolate, need to install scipy')

        # make sure the input is a numpy array
        if not isinstance(new_freq_array, np.ndarray):
            new_freq_array = np.array(new_freq_array)
            
        if period_buffer is not None:
            if 0. < period_buffer < 1.:
                period_buffer += 1.
                print("Warning: period buffer must be > 1. Updating to",period_buffer)

        # check the bounds of the new frequency array
        if bounds_error:
            # YG: the commented block below seems no longer necessary.
            # floater = 1.e-8  # FZ: a small offset to avoid out-of-bound error in spi interpolation module.
            # self._logger.info("massage the new_freq_array's min and max to avoid out-of-bound interp")
            # minindex = np.argmin(new_freq_array)
            # maxindex = np.argmax(new_freq_array)
            # new_freq_array[minindex] += floater
            # new_freq_array[maxindex] -= floater

            # logger.debug("new freq array %s", new_freq_array)
            if self.Z.freq.min() > new_freq_array.min():
                raise ValueError('New frequency minimum of {0:.5g}'.format(new_freq_array.min()) + \
                                 ' is smaller than old frequency minimum of {0:.5g}'.format(self.Z.freq.min()) + \
                                 '.  The new frequency range needs to be within the ' +
                                 'bounds of the old one.')
            if self.Z.freq.max() < new_freq_array.max():
                raise ValueError('New frequency maximum of {0:.5g}'.format(new_freq_array.max()) + \
                                 'is smaller than old frequency maximum of {0:.5g}'.format(self.Z.freq.max()) + \
                                 '.  The new frequency range needs to be within the ' +
                                 'bounds of the old one.')

        # make a new Z object
        new_Z = MTz.Z(z_array=np.zeros((new_freq_array.shape[0], 2, 2),
                                       dtype='complex'),
                      z_err_array=np.zeros((new_freq_array.shape[0], 2, 2)),
                      freq=new_freq_array)

        new_Tipper = MTz.Tipper(tipper_array=np.zeros((new_freq_array.shape[0], 1, 2),
                                                      dtype='complex'),
                                tipper_err_array=np.zeros(
                                    (new_freq_array.shape[0], 1, 2)),
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
                                        (new_freq_array <= f.max()))[0]
                new_f = new_freq_array[new_nz_index]
                
                # apply period buffer
                if type(period_buffer) in [float, int]:
                    new_f_update = []
                    new_nz_index_update = []
                    for ifidx,ifreq in enumerate(new_f):
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
                new_Z.z[new_nz_index, ii, jj] = z_func_real(
                    new_f) + 1j * z_func_imag(new_f)
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
            new_Tipper.tipper[new_nz_index, 0, jj] = t_func_real(new_f) + \
                                                     1j * t_func_imag(new_f)

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

        from mtpy.imaging import plot_mt_response
        # todo change this to the format of the new imaging API
        plot_obj = plot_mt_response.PlotMTResponse(z_object=self.Z,
                                                   t_object=self.Tipper,
                                                   pt_obj=self.pt,
                                                   station=self.station,
                                                   **kwargs)

        return plot_obj
        # raise NotImplementedError


# ==============================================================================
# Site details
# ==============================================================================


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
        self._end_dt = None
        self.id = None
        self.Location = Location()
        self.project = None
        self.run_list = None
        self._start_dt = None
        self.survey = None
        self._date_fmt = '%Y-%m-%d'

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
    @property
    def year_collected(self):
        try:
            return self._start_dt.year
        except TypeError:
            return None
    @year_collected.setter
    def year_collected(self, value):
        pass
    
    @property
    def start_date(self):
        try:
            return self._start_dt.strftime(self._date_fmt)
        except AttributeError:
            return None
        
    @start_date.setter
    def start_date(self, start_str):
        self._start_dt = self._read_date(start_str)
            
    @property
    def end_date(self):
        try:
            return self._end_dt.strftime(self._date_fmt)
        except AttributeError:
            return None
        
    @end_date.setter
    def end_date(self, end_str):
        self._end_dt = self._read_date(end_str)
            
    def _read_date(self, date_str):
        """
        read a date string
        """
        if date_str in [None, 'None', 'none', 'NONE']:
            return None
        try:
            return dt_parser.parse(date_str)
        except dt_parser.ParserError:
            try:
                return dt_parser.parse(date_str, dayfirst=True)
            except dt_parser.ParserError as error:
                raise ValueError(error)
        
# ==============================================================================
# Location class, be sure to put locations in decimal degrees, and note datum
# ==============================================================================


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

        self._northing = None
        self._easting = None
        self.utm_zone = None
        self.elev_units = 'm'
        self.coordinate_system = 'Geographic North'

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, lat):
        self._latitude = gis_tools.assert_lat_value(lat)

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, lon):
        self._longitude = gis_tools.assert_lon_value(lon)

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, elev):
        self._elevation = gis_tools.assert_elevation_value(elev)

    @property
    def easting(self):
        return self._easting

    @easting.setter
    def easting(self, easting):
        try:
            self._easting = float(easting)
        except TypeError:
            self._easting = None

    @property
    def northing(self):
        return self._northing

    @northing.setter
    def northing(self, northing):
        try:
            self._northing = float(northing)
        except TypeError:
            self._northing = None

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


# ==============================================================================
# Field Notes
# ==============================================================================


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

        self.DataQuality = DataQuality()
        self.DataLogger = Instrument()
        self.Electrode_ex = Instrument(**null_emeas.__dict__)

        self.Electrode_ey = Instrument(**null_emeas.__dict__)

        self.Magnetometer_hx = Instrument(**null_hmeas.__dict__)
        self.Magnetometer_hy = Instrument(**null_hmeas.__dict__)
        self.Magnetometer_hz = Instrument(**null_hmeas.__dict__)

        self.Electrode_ex.chtype = 'ex'
        self.Electrode_ey.chtype = 'ey'
        self.Magnetometer_hx.chtype = 'hx'
        self.Magnetometer_hy.chtype = 'hy'
        self.Magnetometer_hz.chtype = 'hz'

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


# ==============================================================================
# Instrument
# ==============================================================================
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

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


# ==============================================================================
# Data Quality
# ==============================================================================


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
        self.author = None

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


# ==============================================================================
# Citation
# ==============================================================================
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

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


# ==============================================================================
# Copyright
# ==============================================================================
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
        self.Citation = Citation()
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
        self.additional_info = None
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


# ==============================================================================
# Provenance
# ==============================================================================


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
        self.creation_time = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
        self.creating_application = 'MTpy'
        self.Creator = Person()
        self.Submitter = Person()

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


# ==============================================================================
# Person
# ==============================================================================


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

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


# ==============================================================================
# Processing
# ==============================================================================


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

    >>> Person(**{'phone':'888-867-5309'})
    """

    def __init__(self, **kwargs):
        self.Software = Software()
        self.notes = None
        self.processed_by = None
        self.sign_convention = 'exp(+i \omega t)'
        self.remote_reference = None
        self.RemoteSite = Site()

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])


class Software(object):
    """
    software
    """

    def __init__(self, **kwargs):
        self.name = None
        self.version = None
        self.Author = Person()

        for key in kwargs:
            setattr(self, key, kwargs[key])


# ==============================================================================
#             Error
# ==============================================================================


class MTError(Exception):
    pass
