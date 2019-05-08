

import os
import shutil

from mtpy.utils.mtpylog import MtPyLog

TEST_MTPY_ROOT = os.path.normpath(
    os.path.abspath(
        os.path.dirname(
            os.path.dirname(__file__)
        )
    )
)  # assume tests is on the root level of mtpy

TEST_DIR = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
TEST_TEMP_DIR = os.path.normpath(os.path.join(TEST_DIR, "temp"))

if not os.path.isdir(TEST_TEMP_DIR):
    os.mkdir(TEST_TEMP_DIR)


def make_temp_dir(dir_name, base_dir=TEST_TEMP_DIR):
    _temp_dir = os.path.normpath(os.path.join(base_dir, dir_name))
    if os.path.isdir(_temp_dir):
        shutil.rmtree(_temp_dir)
    os.mkdir(_temp_dir)
    return _temp_dir


EDI_DATA_DIR = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files'))
EDI_DATA_DIR2 = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files_2'))
AUS_TOPO_FILE = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/AussieContinent_etopo1.asc'))
SAMPLE_DIR = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/model_files'))  # r'E:\Githubz\mtpy\examples\model_files'


# set test logging configure
MtPyLog.load_configure(os.path.join(TEST_DIR, "logging.yml"))
