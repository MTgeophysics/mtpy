# -*- coding: utf-8 -*-
"""
=======================
schema
=======================

Convenience Classes and Functions to deal with the base metadata standards
described by the csv file.

The hope is that only the csv files will need to be changed as the standards
are modified.  The attribute dictionaries are stored in ATTRICT

Created on Wed Apr 29 11:11:31 2020

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import logging

from copy import deepcopy

from mtpy.core.io.emtf_xml import XML_CSV_FN_PATHS
from mtpy.core.metadata.standards import schema

# =============================================================================
#
# =============================================================================


class XMLStandards:
    """
    Helper container to read in csv files and make the appropriate
    dictionaries used in metadata.

    The thought is that only the csv files need to be changed if there is
    a change in standards.

    """

    def __init__(self):
        self.standards_dict = {}
        self.required_keys = schema.REQUIRED_KEYS
        self.accepted_styles = schema.ACCEPTED_STYLES

        self.logger = logging.getLogger("{__name__}.{self.__class__.__name__}")
        self.logger.debug("Initiating Standards")
        self.mt_standards = schema.Standards()

    @property
    def xml_external_url_dict(self):
        return schema.from_csv(
            schema.get_level_fn("xml_external_url", XML_CSV_FN_PATHS)
        )

    @property
    def xml_primary_data_dict(self):
        return schema.from_csv(
            schema.get_level_fn("xml_primary_data", XML_CSV_FN_PATHS)
        )

    @property
    def xml_attachment_dict(self):
        return schema.from_csv(schema.get_level_fn("xml_attachment", XML_CSV_FN_PATHS))

    @property
    def xml_person_dict(self):
        return schema.from_csv(schema.get_level_fn("xml_person", XML_CSV_FN_PATHS))

    @property
    def xml_provenance_dict(self):
        provenance_dict = schema.from_csv(
            schema.get_level_fn("xml_provenance", XML_CSV_FN_PATHS)
        )
        provenance_dict.add_dict(self.xml_person_dict, "creator")
        provenance_dict.add_dict(self.xml_person_dict, "submitter")
        return provenance_dict

    @property
    def xml_copyright_dict(self):
        copyright_dict = schema.from_csv(
            schema.get_level_fn("xml_copyright", XML_CSV_FN_PATHS)
        )
        copyright_dict.add_dict(self.mt_standards.citation_dict, "citation")
        return copyright_dict

    @property
    def xml_emtf_dict(self):
        emtf_dict = schema.from_csv(schema.get_level_fn("xml_emtf", XML_CSV_FN_PATHS))
        emtf_dict.add_dict(self.xml_external_url_dict, "external_url")
        emtf_dict.add_dict(self.xml_primary_data_dict, "primary_data")
        emtf_dict.add_dict(self.xml_attachment_dict, "attachment")
        emtf_dict.add_dict(self.xml_provenance_dict, "provenance")
        emtf_dict.add_dict(self.xml_copyright_dict, "copyright")

        return emtf_dict

    @property
    def xml_data_quality_notes_dict(self):
        dq_notes = schema.from_csv(
            schema.get_level_fn("xml_data_quality_notes", XML_CSV_FN_PATHS)
        )
        dq_notes.add_dict(self.mt_standards.comment_dict, "comments")
        return dq_notes

    @property
    def xml_data_quality_warnings_dict(self):
        dq_warnings = schema.from_csv(
            schema.get_level_fn("xml_data_quality_warnings", XML_CSV_FN_PATHS)
        )
        dq_warnings.add_dict(self.mt_standards.comment_dict, "comments")
        return dq_warnings

    @property
    def xml_site_dict(self):
        site_dict = schema.from_csv(schema.get_level_fn("xml_site", XML_CSV_FN_PATHS))
        site_dict.add_dict(self.mt_standards.location_dict.copy(), "location")
        site_dict.add_dict(
            self.xml_data_quality_notes_dict.copy(), "data_quality_notes"
        )
        site_dict.add_dict(
            self.xml_data_quality_warnings_dict.copy(), "data_quality_warnings"
        )
        return site_dict
    
    @property
    def xml_electrode_dict(self):
        return schema.from_csv(schema.get_level_fn("xml_electrode",
                                                   XML_CSV_FN_PATHS))  
    
    @property
    def xml_dipole_dict(self):
        dipole_dict = schema.from_csv(schema.get_level_fn("xml_dipole",
                                                   XML_CSV_FN_PATHS)) 
        dipole_dict.add_dict(self.xml_electrode_dict, "electrode_01")
        dipole_dict.add_dict(self.xml_electrode_dict, "electrode_02")
        
        return dipole_dict
    
    @property
    def xml_field_notes_dict(self):
        fn_dict = schema.from_csv(schema.get_level_fn("xml_field_notes",
                                                      XML_CSV_FN_PATHS))
        fn_dict.add_dict(self.mt_standards.instrument_dict, "instrument")
        fn_dict.add_dict(self.xml_dipole_dict, "dipole_01")
        fn_dict.add_dict(self.xml_dipole_dict, "dioole_02")
        fn_dict.add_dict(self.mt_standards.instrument_dict, "hx")
        fn_dict.add_dict(self.mt_standards.instrument_dict, "hy")
        fn_dict.add_dict(self.mt_standards.instrument_dict, "hz")
        fn_dict.add_dict(self.mt_standards.time_period_dict)
        
        return fn_dict

    @property
    def ATTR_DICT(self):
        keys = [fn.stem for fn in XML_CSV_FN_PATHS]
        return dict(
            [(key, deepcopy(getattr(self, "{0}_dict".format(key)))) for key in keys]
        )
        self.logger.debug("Successfully made ATTR_DICT")
        

    def summarize_standards(
        self, levels=["survey", "station", "run", "auxiliary", "electric", "magnetic"]
    ):
        """
        Summarize the metadata definitions
        
        :return: DESCRIPTION
        :rtype: TYPE

        """
        summary_dict = schema.BaseDict()
        for name in levels:
            summary_dict.add_dict(getattr(self, "{0}_dict".format(name)), name)

        return summary_dict
