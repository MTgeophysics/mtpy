# package file
import os
from pathlib import Path

CSV_PATH = Path(__file__).parent

CSV_LIST = [
    "xml_external_url.csv",
    "xml_primary_data.csv",
    "xml_attachment.csv",
    "xml_person.csv",
    "xml_provenance.csv",
    "xml_copyright.csv",
    "xml_data_quality_notes.csv",
    "xml_data_quality_warnings.csv",
    "xml_emtf.csv",
    "xml_site.csv",
    "xml_electrode.csv",
    "xml_dipole.csv",
    "xml_field_notes.csv",
    "xml_software.csv",
    "xml_processing_info.csv",
    "xml_statistical_estimates.csv",
    "xml_estimate.csv",
]

XML_CSV_FN_PATHS = [CSV_PATH.joinpath(fn) for fn in CSV_LIST]
