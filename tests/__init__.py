"""
Set paths for testing
"""
from pathlib import Path

# assume tests is on the root level of mtpy
TEST_MTPY_ROOT = Path(__file__).parent.parent

EDI_DATA_DIR = Path(TEST_MTPY_ROOT, "examples/data/edi_files")
EDI_DATA_DIR2 = Path(TEST_MTPY_ROOT, "examples/data/edi_files_2")
EDI_DATA_DIR_BB = Path(TEST_MTPY_ROOT, "data/BBMT")
AUS_TOPO_FILE = Path(TEST_MTPY_ROOT, "examples/data/AussieContinent_etopo1.asc")
SAMPLE_DIR = Path(TEST_MTPY_ROOT, "examples/model_files")
M2D_DIR = Path(TEST_MTPY_ROOT, "examples/data/mare2dem")

# set temporary directory for tests
TEST_DIR = Path(__file__).parent
TEST_TEMP_DIR = Path(TEST_DIR, "temp")
if not TEST_TEMP_DIR.is_dir():
    TEST_TEMP_DIR.mkdir()


def make_temp_dir(dir_name, base_dir=TEST_TEMP_DIR):
    _temp_dir = Path(base_dir, dir_name)
    if _temp_dir.is_dir():
        _temp_dir.rmdir()
    _temp_dir.mkdir()
    return _temp_dir
