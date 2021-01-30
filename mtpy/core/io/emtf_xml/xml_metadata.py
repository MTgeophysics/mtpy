# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 20:27:29 2020

@author: jrpeacock
"""

from mtpy.core.metadata.metadata import (
    Base,
    Citation,
    Comment,
    Location,
    Orientation,
    Instrument,
    write_lines,
)

from mtpy.core.io.emtf_xml.xml_schema import XMLStandards
from mtpy.utils.mttime import MTime

ATTR_DICT = XMLStandards().ATTR_DICT

# =============================================================================
# Metadata objects
# =============================================================================
class ExternalUrl(Base):
    __doc__ = write_lines(ATTR_DICT["xml_external_url"])

    def __init__(self, **kwargs):

        self.description = None
        self.url = None
        super().__init__(attr_dict=ATTR_DICT["xml_external_url"], **kwargs)


class PrimaryData(Base):
    __doc__ = write_lines(ATTR_DICT["xml_primary_data"])

    def __init__(self, **kwargs):

        self.filename = None
        super().__init__(attr_dict=ATTR_DICT["xml_primary_data"], **kwargs)


class Attachment(Base):
    __doc__ = write_lines(ATTR_DICT["xml_attachment"])

    def __init__(self, **kwargs):

        self.filename = None
        self.description = None
        super().__init__(attr_dict=ATTR_DICT["xml_attachment"], **kwargs)


class Person(Base):
    __doc__ = write_lines(ATTR_DICT["xml_person"])

    def __init__(self, **kwargs):

        self.name = None
        self.org = None
        self.org_url = None
        self.email = None
        super().__init__(attr_dict=ATTR_DICT["xml_person"], **kwargs)


class Provenance(Base):
    __doc__ = write_lines(ATTR_DICT["xml_provenance"])

    def __init__(self, **kwargs):

        self._creation_dt = MTime()
        self.creating_application = None
        self.submitter = Person()
        self.creator = Person()
        super().__init__(attr_dict=ATTR_DICT["xml_provenance"], **kwargs)

    @property
    def create_time(self):
        return self._creation_dt.iso_str

    @create_time.setter
    def create_time(self, dt_str):
        self._creation_dt.from_str(dt_str)


class Copyright(Base):
    __doc__ = write_lines(ATTR_DICT["xml_copyright"])

    def __init__(self, **kwargs):

        self.release_status = None
        self.conditions_of_use = None
        self.creating_application = None
        self.citation = Citation()
        super().__init__(attr_dict=ATTR_DICT["xml_copyright"], **kwargs)


class DataQualityNotes(Base):
    __doc__ = write_lines(ATTR_DICT["xml_data_quality_notes"])

    def __init__(self, **kwargs):

        self.good_from_period = None
        self.good_to_period = None
        self.rating = 0
        self.comments = Comment()
        super().__init__(attr_dict=ATTR_DICT["xml_data_quality_notes"], **kwargs)


class DataQualityWarnings(Base):
    __doc__ = write_lines(ATTR_DICT["xml_data_quality_warnings"])

    def __init__(self, **kwargs):

        self.flag = 0
        self.comments = Comment()
        super().__init__(attr_dict=ATTR_DICT["xml_data_quality_notes"], **kwargs)


class Site(Base):
    __doc__ = write_lines(ATTR_DICT["xml_site"])

    def __init__(self, **kwargs):
        self.project = None
        self.survey = None
        self.year_collected = None
        self.country = None
        self.id = None
        self.name = None
        self.acquired_by = None
        self.location = Location()
        self.orientation = Orientation()
        self.run_list = None
        self._start_dt = MTime()
        self._end_dt = MTime()

        super().__init__(attr_dict=ATTR_DICT["xml_site"], **kwargs)

    @property
    def start(self):
        return self._start_dt.iso_str

    @start.setter
    def start(self, value):
        self._start_dt.from_str(value)

    @property
    def end(self):
        return self._end_dt.iso_str

    @end.setter
    def end(self, value):
        self._end_dt.from_str(value)


class Electrode(Base):
    __doc__ = write_lines(ATTR_DICT["xml_electrode"])

    def __init__(self, **kwargs):
        self.location = None
        self.number = 0

        super().__init__(attr_dict=ATTR_DICT["xml_electrode"], **kwargs)


class Dipole(Base):
    __doc__ = write_lines(ATTR_DICT["xml_dipole"])

    def __init__(self, **kwargs):
        self.manufacturer = None
        self.length = None
        self.azimuth = None
        self.electrode_01 = Electrode()
        self.electrode_02 = Electrode()

        super().__init__(attr_dict=ATTR_DICT["xml_dipole"], **kwargs)


class FieldNotes(Base):
    __doc__ = write_lines(ATTR_DICT["xml_field_notes"])

    def __init__(self, **kwargs):
        self.errors = None
        self.run = None
        self._start_dt = MTime()
        self._end_dt = MTime()
        self.instrument = Instrument()
        self.hx = Instrument()
        self.hy = Instrument()
        self.hz = Instrument()
        self.dipole_01 = Dipole()
        self.dipole_02 = Dipole()

        super().__init__(attr_dict=ATTR_DICT["xml_field_notes"], **kwargs)

    @property
    def start(self):
        return self._start_dt.iso_str

    @start.setter
    def start(self, value):
        self._start_dt.from_str(value)

    @property
    def end(self):
        return self._end_dt.iso_str

    @end.setter
    def end(self, value):
        self._end_dt.from_str(value)


class Software(Base):
    __doc__ = write_lines(ATTR_DICT["xml_software"])

    def __init__(self, **kwargs):
        self.author = None
        self._last_mod_dt = MTime()
        self.remote_ref = None

        super().__init__(attr_dict=ATTR_DICT["xml_software"], **kwargs)

    @property
    def last_mod(self):
        return self._last_mod_dt.iso_str

    @last_mod.setter
    def last_mod(self, value):
        self._last_mod_dt.from_str(value)


class ProcessingInfo(Base):
    __doc__ = write_lines(ATTR_DICT["xml_processing_info"])

    def __init__(self, **kwargs):
        self.sign_convention = None
        self.processed_by = None
        self.remote_ref = None
        self.processing_software = Software()
        self.processing_tag = None

        super().__init__(attr_dict=ATTR_DICT["xml_processing_info"], **kwargs)


class Estimate(Base):
    __doc__ = write_lines(ATTR_DICT["xml_estimate"])

    def __init__(self, **kwargs):
        self.name = None
        self.type = None
        self.description = None
        self.tag = None
        self.external_url = None
        self.intention = None

        super().__init__(attr_dict=ATTR_DICT["xml_estimate"], **kwargs)


class StatisticalEstimates(Base):
    __doc__ = write_lines(ATTR_DICT["xml_statistical_estimates"])

    def __init__(self, **kwargs):
        self.estimates_list = []

        super().__init__(attr_dict=ATTR_DICT["xml_statistical_estimates"], **kwargs)
