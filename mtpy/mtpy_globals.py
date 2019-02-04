"""
Description:
    keep all the mtpy global params constants in this module.

Author: fei.zhang@ga.gov.au

FZ Last Updated: 2017-12-04
"""

import os
import tempfile

epsg_dict = {28350: ['+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 50],
             28351: ['+proj=utm +zone=51 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 51],
             28352: ['+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 52],
             28353: ['+proj=utm +zone=53 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 53],
             28354: ['+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 54],
             28355: ['+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 55],
             28356: ['+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 56],
             3112: [
                 '+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',
                 0],
             4326: ['+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', 0]}


# hack to get the path2/repositoryDir like: C:/Github/mtpy
MTPY_ROOT = os.path.normpath(
    os.path.abspath(
        os.path.dirname(
            os.path.dirname(__file__)
        )
    )
)

# print("MTPY_ROOT = ", MTPY_ROOT)

EDI_DATA_DIR = os.path.normpath(
    os.path.join(MTPY_ROOT, 'examples/data/edi_files'))
EDI_DATA_DIR2 = os.path.normpath(
    os.path.join(MTPY_ROOT, 'examples/data/edi_files_2'))

AUS_TOPO_FILE = os.path.normpath(
    os.path.join(MTPY_ROOT, 'examples/data/AussieContinent_etopo1.asc'))
SAMPLE_DIR = os.path.normpath(
    os.path.join(MTPY_ROOT, 'examples/model_files'))  # r'E:\Githubz\mtpy\examples\model_files'


SYSTEM_TEMP_DIR = tempfile.gettempdir()

NEW_TEMP_DIR=tempfile.mkdtemp(prefix="mtpy_tmpdir_")

# print("SYSTEM_TEMP_DIR = ", SYSTEM_TEMP_DIR) # in my Windows =  c:\users\u25656\appdata\local\temp\1
#
# print ("NEW_TEMP_DIR = ", NEW_TEMP_DIR)

