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
    "xml_emtf.csv",
	"xml_site.csv",
]

XML_CSV_FN_PATHS = [CSV_PATH.joinpath(fn) for fn in CSV_LIST]
