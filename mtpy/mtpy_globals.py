"""
Description:
    keep all the mtpy global params constants in this module.

Author: fei.zhang@ga.gov.au

FZ Last Updated: 2017-12-04
JP (2021-01-18) updated to use Path and get relative path locations.

"""
from pathlib import Path
import tempfile

epsg_dict = {
    28350: [
        "+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        50,
    ],
    28351: [
        "+proj=utm +zone=51 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        51,
    ],
    28352: [
        "+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        52,
    ],
    28353: [
        "+proj=utm +zone=53 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        53,
    ],
    28354: [
        "+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        54,
    ],
    28355: [
        "+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        55,
    ],
    28356: [
        "+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        56,
    ],
    3112: [
        "+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
        0,
    ],
    4326: ["+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 0],
}


MTPY_ROOT = Path(__file__).parent.parent
EDI_DATA_DIR = Path(MTPY_ROOT, "examples/data/edi_files")
EDI_DATA_DIR2 = Path(MTPY_ROOT, "examples/data/edi_files_2")
AUS_TOPO_FILE = Path(MTPY_ROOT, "examples/data/AussieContinent_etopo1.asc")
SAMPLE_DIR = Path(MTPY_ROOT, "examples/model_files")

SYSTEM_TEMP_DIR = tempfile.gettempdir()

NEW_TEMP_DIR = tempfile.mkdtemp(prefix="mtpy_tmpdir_")
