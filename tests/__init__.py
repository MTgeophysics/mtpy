from __future__ import print_function

import os
import shutil
import sys
import matplotlib

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive svg backend", file=sys.stderr)
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt

    plt.ion()

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


def plt_wait(seconds):
    if plt.isinteractive() and plt.get_fignums():
        try:
            plt.pause(seconds)
        except:
            pass


def plt_close(to_close=None):
    if plt.get_fignums():
        if to_close is None:
            plt.close()
        else:
            plt.close(to_close)


EDI_DATA_DIR = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files'))
EDI_DATA_DIR2 = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/edi_files_2'))
AUS_TOPO_FILE = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/data/AussieContinent_etopo1.asc'))
SAMPLE_DIR = os.path.normpath(
    os.path.join(TEST_MTPY_ROOT, 'examples/model_files'))  # r'E:\Githubz\mtpy\examples\model_files'
