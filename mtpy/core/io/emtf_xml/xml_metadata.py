# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 20:27:29 2020

@author: jrpeacock
"""

from mtpy.core.metadata.metadata import (Base, Citation, Comment, Location, 
                                         Orientation, write_lines)
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
        super().__init__(attr_dict=ATTR_DICT["xml_external_url"],
                         **kwargs)

class PrimaryData(Base):
    __doc__ = write_lines(ATTR_DICT["xml_primary_data"])

    def __init__(self, **kwargs):

        self.filename = None
        super().__init__(attr_dict=ATTR_DICT["xml_primary_data"],
                         **kwargs)

class Attachment(Base):
    __doc__ = write_lines(ATTR_DICT["xml_attachment"])

    def __init__(self, **kwargs):

        self.filename = None
        self.description = None
        super().__init__(attr_dict=ATTR_DICT["xml_attachment"],
                         **kwargs)

class Person(Base):
    __doc__ = write_lines(ATTR_DICT["xml_person"])

    def __init__(self, **kwargs):

        self.name = None
        self.org = None
        self.org_url = None
        self.email = None
        super().__init__(attr_dict=ATTR_DICT["xml_person"],
                         **kwargs)
        
class Provenance(Base):
    __doc__ = write_lines(ATTR_DICT["xml_provenance"])

    def __init__(self, **kwargs):

        self._creation_dt = MTime()
        self.creating_application = None
        self.submitter = Person()
        self.creator = Person()
        super().__init__(attr_dict=ATTR_DICT["xml_provenance"],
                         **kwargs)
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
        super().__init__(attr_dict=ATTR_DICT["xml_copyright"],
                         **kwargs)
    
class DataQualityNotes(Base):
    __doc__ = write_lines(ATTR_DICT["xml_data_quality_notes"])

    def __init__(self, **kwargs):

        self.good_from_period = None
        self.good_to_period = None
        self.rating = 0
        self.comments = Comment()
        super().__init__(attr_dict=ATTR_DICT["xml_data_quality_notes"],
                         **kwargs)
        
class DataQualityWarnings(Base):
    __doc__ = write_lines(ATTR_DICT["xml_data_quality_warnings"])

    def __init__(self, **kwargs):

        self.flag = 0
        self.comments = Comment()
        super().__init__(attr_dict=ATTR_DICT["xml_data_quality_notes"],
                         **kwargs)
        
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
        
        super().__init__(attr_dict=ATTR_DICT["xml_site"],
                         **kwargs)
        
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
        
        
        