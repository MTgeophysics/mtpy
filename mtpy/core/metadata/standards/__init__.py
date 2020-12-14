# package file
import os
from pathlib import Path

CSV_PATH = Path(os.path.dirname(os.path.abspath(__file__)))

CSV_LIST = [
    "auxiliary.csv",
    "battery.csv",
    "channel.csv",
    "citation.csv",
	"comment.csv",
    "copyright.csv",
    "data_quality.csv",
    "datalogger.csv",
    "fdsn.csv",
    "rating.csv",
    "declination.csv",
    "diagnostic.csv",
    "electric.csv",
    "filter.csv",
    "filtered.csv",
    "instrument.csv",
    "location.csv",
    "magnetic.csv",
    "person.csv",
    "provenance.csv",
    "run.csv",
    "software.csv",
    "station.csv",
    "survey.csv",
    "timing_system.csv",
    "time_period.csv",
    "orientation.csv",
    "transfer_function.csv",
]

CSV_FN_PATHS = [CSV_PATH.joinpath(fn) for fn in CSV_LIST]
