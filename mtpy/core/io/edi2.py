# -*- coding: utf-8 -*-
"""
EDI
=======

Read and write EDI file from an MT object


Created on Fri Nov 13 10:58:20 2020

:author: Jared Peacock

:license: MIT

"""


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
    selt.station_metadata.locationdatum = edi_obj.Header.datum
    selt.station_metadata.locationelev_units = edi_obj.Define_measurement.units
    selt.station_metadata.locationcoordinate_system = edi_obj.Header.coordinate_system
    if hasattr(edi_obj.Header,'enddate'):
        self.Site.end_date = edi_obj.Header.enddate

    selt.station_metadata.locationdeclination = edi_obj.Header.declination

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
    header.coordinate_system = selt.station_metadata.locationcoordinate_system
    header.dataid = self.Site.id
    header.datum = selt.station_metadata.locationdatum
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
    header.declination = selt.station_metadata.locationdeclination
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
    define_meas.reftype = selt.station_metadata.locationcoordinate_system
    define_meas.units = selt.station_metadata.locationelev_units

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

