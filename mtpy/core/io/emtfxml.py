# -*- coding: utf-8 -*-
"""
MTXML
====================

Read and write EMTF XML files


Created on Fri Nov 13 11:00:04 2020

:author: Jared Peacock

:license: MIT

"""


def _read_xml_file(self, xml_fn):
    """
    read xml file
    """

    if not os.path.isfile(xml_fn):
        raise MTError("Could not find {0}, check path.".format(xml_fn))

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
        if s_attr in ["_name", "_attr", "_value"]:
            continue
        x_obj = getattr(xml_obj.Site, s_attr)
        name = x_obj.name.lower()
        if name == "acquiredby":
            name = "acquired_by"
        elif name == "end":
            name = "end_date"
        elif name == "start":
            name = "start_date"
        elif name == "runlist":
            name = "run_list"
        elif name == "yearcollected":
            name = "year_collected"
        elif name == "datecollected":
            name = "date_collected"

        value = x_obj.value
        if name == "location":
            for l_attr in list(xml_obj.Site.Location.__dict__.keys()):
                if l_attr in ["_name", "_attr", "_value"]:
                    continue
                l_obj = getattr(xml_obj.Site.Location, l_attr)
                name = l_obj.name.lower()
                value = l_obj.value
                if name == "elevation":
                    units = l_obj.attr["units"]
                    selt.station_metadata.locationelev_units = units

                elif name == "declination":
                    units = l_obj.attr["epoch"]
                    selt.station_metadata.locationdeclination_epoch = units
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

    self.Notes.info_dict["notes"] = str(xml_obj.Notes.value)
    self.Notes.info_list = ["notes = {0}".format(str(xml_obj.Notes.value))]


def _xml_get_field_notes(self, xml_obj):
    """
    get field notes information
    """

    for f_attr in list(xml_obj.FieldNotes.__dict__.keys()):
        if f_attr.lower() == "instrument":
            for i_attr in list(xml_obj.FieldNotes.Instrument.__dict__.keys()):
                if i_attr in ["_name", "_attr", "_value"]:
                    continue
                i_obj = getattr(xml_obj.FieldNotes.Instrument, i_attr)
                name = i_obj.name.lower()
                value = i_obj.value
                setattr(self.FieldNotes.DataLogger, name, value)

        elif "dipole" in f_attr.lower():
            xml_d_obj = getattr(xml_obj.FieldNotes, f_attr)
            azm = 0.0
            length = 0.0
            try:
                comp = xml_d_obj.attr["name"].lower()
            except KeyError:
                comp = "ex"
            try:
                t = xml_d_obj.attr["type"].lower()
                if comp == "ex":
                    setattr(self.FieldNotes.Electrode_ex, "type", t)
                elif comp == "ey":
                    setattr(self.FieldNotes.Electrode_ey, "type", t)
            except KeyError:
                pass

            for e_attr in list(xml_d_obj.__dict__.keys()):
                if e_attr in ["_name", "_attr", "_value"]:
                    continue
                e_obj = getattr(xml_d_obj, e_attr)
                name = e_obj.name.lower()
                value = e_obj.value

                if name == "azimuth":
                    azm = float(value)

                if name == "length":
                    length = float(value)

                if comp == "ex":
                    setattr(self.FieldNotes.Electrode_ex, name, value)
                elif comp == "ey":
                    setattr(self.FieldNotes.Electrode_ey, name, value)

            # need to set x, y, x2, y2
            x = 0
            y = 0
            x2 = length * np.cos(np.deg2rad(azm))
            y2 = length * np.sin(np.deg2rad(azm))

            for name, value in zip(["x", "y", "x2", "y2"], [x, y, x2, y2]):
                if comp == "ex":
                    setattr(self.FieldNotes.Electrode_ex, name, value)
                elif comp == "ey":
                    setattr(self.FieldNotes.Electrode_ey, name, value)

        elif "magnetometer" in f_attr.lower():
            xml_d_obj = getattr(xml_obj.FieldNotes, f_attr)
            try:
                comp = xml_d_obj.attr["name"].lower()
            except KeyError:
                try:
                    comp = xml_d_obj.attr["type"].lower()
                except KeyError:
                    comp = "hx"

            try:
                t = xml_d_obj.attr["type"].lower()
                if comp == "hx":
                    setattr(self.FieldNotes.Magnetometer_hx, "type", t)
                elif comp == "hy":
                    setattr(self.FieldNotes.Magnetometer_hy, "type", t)
                elif comp == "hz":
                    setattr(self.FieldNotes.Magnetometer_hz, "type", t)
                elif comp == "fluxgate":
                    setattr(self.FieldNotes.Magnetometer_hx, "type", t)
                    setattr(self.FieldNotes.Magnetometer_hy, "type", t)
                    setattr(self.FieldNotes.Magnetometer_hz, "type", t)
                else:
                    pass
            except KeyError:
                pass

            for m_attr in list(xml_d_obj.__dict__.keys()):
                if m_attr in ["_name", "_attr", "_value"]:
                    continue
                m_obj = getattr(xml_obj.FieldNotes.Magnetometer, m_attr)
                name = m_obj.name.lower()
                value = m_obj.value

                if comp == "hx":
                    setattr(self.FieldNotes.Magnetometer_hx, name, value)
                elif comp == "hy":
                    setattr(self.FieldNotes.Magnetometer_hy, name, value)
                elif comp == "hz":
                    setattr(self.FieldNotes.Magnetometer_hz, name, value)
                elif comp == "fluxgate":
                    setattr(self.FieldNotes.Magnetometer_hx, name, value)
                    setattr(self.FieldNotes.Magnetometer_hy, name, value)
                    setattr(self.FieldNotes.Magnetometer_hz, name, value)

        elif "dataquality" in f_attr.lower():
            obj = getattr(xml_obj.FieldNotes, f_attr)
            for d_attr in list(obj.__dict__.keys()):
                if d_attr in ["_name", "_attr", "_value"]:
                    continue
                d_obj = getattr(obj, d_attr)
                name = d_obj.name.lower()
                if name == "goodfromperiod":
                    name = "good_from_period"
                elif name == "goodtoperiod":
                    name = "good_to_period"
                elif name == "comments":
                    setattr(self.FieldNotes.DataQuality, "author", d_obj.attr["author"])
                if name == "comments" and f_attr.lower() == "dataqualitywarnings":
                    name = "warnings_" + name
                value = d_obj.value

                setattr(self.FieldNotes.DataQuality, name, value)


def _xml_get_copyright(self, xml_obj):
    """
    get copyright information
    """

    for f_attr in list(xml_obj.Copyright.__dict__.keys()):
        if f_attr in ["_name", "_attr", "_value"]:
            continue
        if f_attr.lower() == "citation":
            for i_attr in list(xml_obj.Copyright.Citation.__dict__.keys()):
                if i_attr in ["_name", "_attr", "_value"]:
                    continue
                i_obj = getattr(xml_obj.Copyright.Citation, i_attr)
                name = i_obj.name.lower()
                value = i_obj.value
                setattr(self.Copyright.Citation, name, value)
        else:
            obj = getattr(xml_obj.Copyright, f_attr)
            name = obj.name.lower()
            if name == "releasestatus":
                name = "release_status"
            elif name == "conditionsofuse":
                name = "conditions_of_use"
            elif name == "additionalinfo":
                name = "additional_info"
            value = obj.value.replace("\n", "")

            setattr(self.Copyright, name, value)


def _xml_get_provenance(self, xml_obj):
    """
    get provenance infor
    """
    for f_attr in list(xml_obj.Provenance.__dict__.keys()):
        if f_attr in ["_name", "_attr", "_value"]:
            continue
        if f_attr.lower() in ["creator", "submitter"]:
            obj = getattr(xml_obj.Provenance, f_attr)
            s_obj = getattr(self.Provenance, f_attr)
            for i_attr in list(obj.__dict__.keys()):
                if i_attr in ["_name", "_attr", "_value"]:
                    continue
                i_obj = getattr(obj, i_attr)
                name = i_obj.name.lower()
                value = i_obj.value
                setattr(s_obj, name, value)
        else:
            obj = getattr(xml_obj.Provenance, f_attr)
            name = obj.name.lower()
            if name == "creationtime":
                name = "creation_time"
            elif name == "creatingapplication":
                name = "creating_application"
            value = obj.value

            setattr(self.Provenance, name, value)


def _xml_get_processing(self, xml_obj):
    """
    get processing info
    """

    for f_attr in list(xml_obj.ProcessingInfo.__dict__.keys()):
        if f_attr in ["_name", "_attr", "_value"]:
            continue
        if "software" in f_attr.lower():
            obj = getattr(xml_obj.ProcessingInfo, f_attr)
            for i_attr in list(obj.__dict__.keys()):
                if i_attr in ["_name", "_attr", "_value"]:
                    continue
                i_obj = getattr(obj, i_attr)
                name = i_obj.name.lower()
                if name == "lastmod":
                    name = "last_modified"
                value = i_obj.value
                if name == "author":
                    value = Person()
                    value.name = i_obj.value
                setattr(self.Processing.Software, name, value)
        elif "remoteinfo" in f_attr.lower():
            obj = getattr(xml_obj.ProcessingInfo, f_attr)
            for i_attr in list(obj.__dict__.keys()):
                if i_attr in ["_name", "_attr", "_value"]:
                    continue
                if i_attr.lower() == "location":
                    loc_obj = getattr(obj, i_attr)

                    for l_attr in list(loc_obj.__dict__.keys()):
                        if l_attr in ["_name", "_attr", "_value"]:
                            continue
                        l_obj = getattr(loc_obj, l_attr)
                        name = l_obj.name.lower()
                        value = l_obj.value
                        setattr(self.Processing.RemoteSite.Location, name, value)

                else:
                    i_obj = getattr(obj, i_attr)
                    name = i_obj.name.lower()
                    if name == "yearcollected":
                        name = "year_collected"
                        value = i_obj.value
                    setattr(self.Processing.RemoteSite, name, value)
        else:
            obj = getattr(xml_obj.ProcessingInfo, f_attr)
            name = obj.name.lower()
            if name == "signconvention":
                name = "sign_convention"
            elif name == "remoteref":
                name = "remote_ref"
            elif name == "processedby":
                name = "processed_by"
            elif name == "processingtag":
                name = "processing_tag"
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
    xml_obj.PrimaryData.Filename.value = os.path.basename(self.fn)[:-4] + ".png"

    xml_obj.ProductId.value = "{0}.{1}".format(
        self.station.upper(), self.Site.year_collected
    )

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
        xml_obj.Site.Survey.value = ",".join(self.Site.survey)
    else:
        xml_obj.Site.Survey.value = self.Site.survey
    xml_obj.Site.Id.value = self.Site.id
    xml_obj.Site.AcquiredBy.value = self.Site.acquired_by
    xml_obj.Site.Start.value = self.Site.start_date
    xml_obj.Site.End.value = self.Site.end_date
    xml_obj.Site.RunList.value = self.Site.run_list
    xml_obj.Site.Orientation.value = "geomagnetic"
    try:
        xml_obj.Site.Orientation.attr = {
            "angle_to_geographic_north": "{0:.2f}".format(
                selt.station_metadata.locationdeclination
            )
        }
    except ValueError:
        xml_obj.Site.Orientation.attr = {"angle_to_geographic_north": "0.00"}

    xml_obj.Site.Location.Latitude.value = self.lat
    xml_obj.Site.Location.Longitude.value = self.lon
    xml_obj.Site.Location.Elevation.value = self.elev
    xml_obj.Site.Location.Elevation.attr = {
        "units": selt.station_metadata.locationelev_units
    }
    xml_obj.Site.Location.Declination.value = selt.station_metadata.locationdeclination
    xml_obj.Site.Location.Declination.attr = {
        "epoch": selt.station_metadata.locationdeclination_epoch
    }

    return xml_obj


def _xml_set_field_notes(self, xml_obj):
    """
    Set the FieldNotes attributes of the xml object
    """

    xml_obj.FieldNotes.Instrument.Type.value = self.FieldNotes.DataLogger.type
    xml_obj.FieldNotes.Instrument.Id.value = self.FieldNotes.DataLogger.id
    xml_obj.FieldNotes.Instrument.Manufacturer.value = (
        self.FieldNotes.DataLogger.manufacturer
    )

    # EX
    xml_obj.FieldNotes.Dipole.Type.value = self.FieldNotes.Electrode_ex.type
    xml_obj.FieldNotes.Dipole.Id.value = self.FieldNotes.Electrode_ex.id
    xml_obj.FieldNotes.Dipole.Manufacturer.value = (
        self.FieldNotes.Electrode_ex.manufacturer
    )
    xml_obj.FieldNotes.Dipole.attr = {"name": "EX"}

    length = np.sqrt(
        (self.FieldNotes.Electrode_ex.x2 - self.FieldNotes.Electrode_ex.x) ** 2
        + (self.FieldNotes.Electrode_ex.y2 - self.FieldNotes.Electrode_ex.y) ** 2
    )
    xml_obj.FieldNotes.Dipole.Length.value = length
    try:
        azm = np.arctan(
            (self.FieldNotes.Electrode_ex.y2 - self.FieldNotes.Electrode_ex.y)
            / (self.FieldNotes.Electrode_ex.x2 - self.FieldNotes.Electrode_ex.x)
        )
    except ZeroDivisionError:
        azm = 0.0
    xml_obj.FieldNotes.Dipole.Azimuth.value = np.degrees(azm)
    xml_obj.FieldNotes.Dipole.Channel.value = self.FieldNotes.Electrode_ex.acqchan

    # EY
    xml_obj.FieldNotes.Dipole_00.Type.value = self.FieldNotes.Electrode_ey.type
    xml_obj.FieldNotes.Dipole_00.Id.value = self.FieldNotes.Electrode_ey.id
    xml_obj.FieldNotes.Dipole_00.Manufacturer.value = (
        self.FieldNotes.Electrode_ey.manufacturer
    )
    xml_obj.FieldNotes.Dipole_00.attr = {"name": "EY"}
    length = np.sqrt(
        (self.FieldNotes.Electrode_ey.x2 - self.FieldNotes.Electrode_ey.x) ** 2
        + (self.FieldNotes.Electrode_ey.y2 - self.FieldNotes.Electrode_ey.y) ** 2
    )
    xml_obj.FieldNotes.Dipole_00.Length.value = length
    try:
        azm = np.arctan(
            (self.FieldNotes.Electrode_ey.y2 - self.FieldNotes.Electrode_ey.y)
            / (self.FieldNotes.Electrode_ey.x2 - self.FieldNotes.Electrode_ey.x)
        )
    except ZeroDivisionError:
        azm = 90.0
    xml_obj.FieldNotes.Dipole_00.Azimuth.value = np.degrees(min([np.pi / 2, azm]))
    xml_obj.FieldNotes.Dipole_00.Channel.value = self.FieldNotes.Electrode_ey.acqchan

    # HX
    xml_obj.FieldNotes.Magnetometer.Type.value = self.FieldNotes.Magnetometer_hx.type
    xml_obj.FieldNotes.Magnetometer.Id.value = self.FieldNotes.Magnetometer_hx.id
    xml_obj.FieldNotes.Magnetometer.Manufacturer.value = (
        self.FieldNotes.Magnetometer_hx.manufacturer
    )
    xml_obj.FieldNotes.Magnetometer.attr = {"name": "HX"}
    xml_obj.FieldNotes.Magnetometer.Azimuth.value = self.FieldNotes.Magnetometer_hx.azm
    xml_obj.FieldNotes.Magnetometer.Channel.value = (
        self.FieldNotes.Magnetometer_hx.acqchan
    )

    # HY
    xml_obj.FieldNotes.Magnetometer_00.Type.value = self.FieldNotes.Magnetometer_hy.type
    xml_obj.FieldNotes.Magnetometer_00.Id.value = self.FieldNotes.Magnetometer_hy.id
    xml_obj.FieldNotes.Magnetometer_00.Manufacturer.value = (
        self.FieldNotes.Magnetometer_hy.manufacturer
    )
    xml_obj.FieldNotes.Magnetometer_00.attr = {"name": "HY"}
    xml_obj.FieldNotes.Magnetometer_00.Azimuth.value = (
        self.FieldNotes.Magnetometer_hy.azm
    )
    xml_obj.FieldNotes.Magnetometer_00.Channel.value = (
        self.FieldNotes.Magnetometer_hy.acqchan
    )

    # HZ
    xml_obj.FieldNotes.Magnetometer_01.Type.value = self.FieldNotes.Magnetometer_hz.type
    xml_obj.FieldNotes.Magnetometer_01.Id.value = self.FieldNotes.Magnetometer_hz.id
    xml_obj.FieldNotes.Magnetometer_01.Manufacturer.value = (
        self.FieldNotes.Magnetometer_hz.manufacturer
    )
    xml_obj.FieldNotes.Magnetometer_01.attr = {"name": "HZ"}
    xml_obj.FieldNotes.Magnetometer_01.Azimuth.value = (
        self.FieldNotes.Magnetometer_hz.azm
    )
    xml_obj.FieldNotes.Magnetometer_01.Channel.value = (
        self.FieldNotes.Magnetometer_hz.acqchan
    )

    # Data Quality Notes
    xml_obj.FieldNotes.DataQualityNotes.Rating.value = (
        self.FieldNotes.DataQuality.rating
    )
    xml_obj.FieldNotes.DataQualityNotes.GoodFromPeriod.value = (
        self.FieldNotes.DataQuality.good_from_period
    )
    xml_obj.FieldNotes.DataQualityNotes.GoodToPeriod.value = (
        self.FieldNotes.DataQuality.good_to_period
    )
    xml_obj.FieldNotes.DataQualityNotes.Comments.value = (
        self.FieldNotes.DataQuality.comments
    )
    xml_obj.FieldNotes.DataQualityNotes.Comments.attr = {
        "author": self.FieldNotes.DataQuality.author
    }
    # Data Quality Warnings
    xml_obj.FieldNotes.DataQualityWarnings.Flag.value = (
        self.FieldNotes.DataQuality.warnings_flag
    )
    xml_obj.FieldNotes.DataQualityWarnings.Comments.value = (
        self.FieldNotes.DataQuality.warnings_comments
    )
    xml_obj.FieldNotes.DataQualityWarnings.Comments.attr = {
        "author": self.FieldNotes.DataQuality.author
    }

    return xml_obj


def _xml_set_processing(self, xml_obj):
    """
    Set the Processing attributes of the xml object
    """

    xml_obj.ProcessingInfo.ProcessedBy.value = self.Processing.processed_by
    xml_obj.ProcessingInfo.ProcessingSoftware.Name.value = self.Processing.Software.name
    xml_obj.ProcessingInfo.ProcessingSoftware.Author.value = (
        self.Processing.Software.Author.name
    )
    xml_obj.ProcessingInfo.ProcessingSoftware.Version.value = (
        self.Processing.Software.version
    )

    # TODO: Need to find a way to put in processing parameters.

    xml_obj.ProcessingInfo.SignConvention.value = self.Processing.sign_convention

    xml_obj.ProcessingInfo.RemoteRef.value = self.Processing.remote_reference
    xml_obj.ProcessingInfo.RemoteInfo.Project.value = self.Processing.RemoteSite.project
    xml_obj.ProcessingInfo.RemoteInfo.Survey.value = self.Processing.RemoteSite.survey
    xml_obj.ProcessingInfo.RemoteInfo.ID.value = self.Processing.RemoteSite.id
    xml_obj.ProcessingInfo.RemoteInfo.YearCollected.value = (
        self.Processing.RemoteSite.year_collected
    )
    xml_obj.ProcessingInfo.RemoteInfo.AcquiredBy.value = (
        self.Processing.RemoteSite.acquired_by
    )
    xml_obj.ProcessingInfo.RemoteInfo.Location.Latitude.value = (
        self.Processing.RemoteSite.Location.latitude
    )
    xml_obj.ProcessingInfo.RemoteInfo.Location.Longitude.value = (
        self.Processing.RemoteSite.Location.longitude
    )
    xml_obj.ProcessingInfo.RemoteInfo.Location.Elevation.value = (
        self.Processing.RemoteSite.Location.elevation
    )
    xml_obj.ProcessingInfo.RemoteInfo.Location.Elevation.attr = {
        "units": self.Processing.RemoteSite.Location.elev_units
    }
    xml_obj.ProcessingInfo.RemoteInfo.Location.attr = {
        "datum": self.Processing.RemoteSite.Location.datum
    }

    return xml_obj


def _xml_set_provenance(self, xml_obj):
    """
    Set the Provenance attributes of the xml object
    """

    xml_obj.Provenance.CreatingApplication.value = "MTpy 0.1.0"

    xml_obj.Provenance.Submitter.Name.value = self.Provenance.Submitter.name
    xml_obj.Provenance.Submitter.Email.value = self.Provenance.Submitter.email
    xml_obj.Provenance.Submitter.Org.value = self.Provenance.Submitter.organization
    xml_obj.Provenance.Submitter.OrgURL.value = (
        self.Provenance.Submitter.organization_url
    )

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
        xml_obj.Copyright.ConditionsOfUse.value = "".join(
            self.Copyright.conditions_of_use
        )
    else:
        xml_obj.Copyright.ConditionsOfUse.value = self.Copyright.conditions_of_use
    xml_obj.Copyright.ReleaseStatus.value = self.Copyright.release_status
    xml_obj.Copyright.AdditionalInfo.value = self.Copyright.additional_info

    return xml_obj


def _xml_set_site_layout(self, xml_obj):
    """
    set the site layout from define measurement
    """

    xml_obj.SiteLayout.InputChannels.Magnetic_hx.attr = {
        "name": "Hx",
        "orientation": "{0:.2f}".format(self.FieldNotes.Magnetometer_hx.azm),
        "x": "{0:.2f}".format(self.FieldNotes.Magnetometer_hx.x),
        "y": "{0:.2f}".format(self.FieldNotes.Magnetometer_hx.y),
        "z": "0.00",
    }
    xml_obj.SiteLayout.InputChannels.Magnetic_hy.attr = {
        "name": "Hy",
        "orientation": "{0:.2f}".format(max([90, self.FieldNotes.Magnetometer_hy.azm])),
        "x": "{0:.2f}".format(self.FieldNotes.Magnetometer_hy.x),
        "y": "{0:.2f}".format(self.FieldNotes.Magnetometer_hy.y),
        "z": "0.00",
    }
    xml_obj.SiteLayout.OutputChannels.Magnetic_hz.attr = {
        "name": "Hz",
        "orientation": "{0:.2f}".format(self.FieldNotes.Magnetometer_hz.azm),
        "x": "{0:.2f}".format(self.FieldNotes.Magnetometer_hz.x),
        "y": "{0:.2f}".format(self.FieldNotes.Magnetometer_hz.y),
        "z": "0.00",
    }
    try:
        azm = np.arctan(
            (self.FieldNotes.Electrode_ex.y2 - self.FieldNotes.Electrode_ex.y)
            / (self.FieldNotes.Electrode_ex.x2 - self.FieldNotes.Electrode_ex.x)
        )
    except ZeroDivisionError:
        azm = 90.0
    xml_obj.SiteLayout.OutputChannels.Electric_ex.attr = {
        "name": "Ex",
        "orientation": "{0:.2f}".format(azm),
        "x": "{0:.2f}".format(self.FieldNotes.Electrode_ex.x),
        "y": "{0:.2f}".format(self.FieldNotes.Electrode_ex.y),
        "z": "0.00",
        "x2": "{0:.2f}".format(self.FieldNotes.Electrode_ex.x2),
        "y2": "{0:.2f}".format(self.FieldNotes.Electrode_ex.y2),
        "z2": "0.00",
    }
    try:
        azm = np.arctan(
            (self.FieldNotes.Electrode_ey.y2 - self.FieldNotes.Electrode_ey.y)
            / (self.FieldNotes.Electrode_ey.x2 - self.FieldNotes.Electrode_ey.x)
        )
    except ZeroDivisionError:
        azm = 90.0
    xml_obj.SiteLayout.OutputChannels.Electric_ey.attr = {
        "name": "Ey",
        "orientation": "{0:.2f}".format(azm),
        "x": "{0:.2f}".format(self.FieldNotes.Electrode_ey.x),
        "y": "{0:.2f}".format(self.FieldNotes.Electrode_ey.y),
        "z": "0.00",
        "x2": "{0:.2f}".format(self.FieldNotes.Electrode_ey.x2),
        "y2": "{0:.2f}".format(self.FieldNotes.Electrode_ey.y2),
        "z2": "0.00",
    }
    return xml_obj
