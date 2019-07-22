#! /usr/bin/env python
"""
Description:
    Example python script
    plot 3D penetration depth for a folder of EDI files

    input = path2edifolder

CreationDate:   23/03/2018
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     23/03/2018   FZ

"""

import os
import sys
import tempfile
from mtpy.imaging import penetration_depth2d as pen2d
import matplotlib.pyplot as plt


try:
    # PACK_ROOT = os.environ['PACK_ROOT']
    # mtpy_path = os.path.join(PACK_ROOT, 'mtpy')
    mtpy_path = os.environ['MTPY_ROOT']
except:
    print("Warn: The environment variable MTPY_ROOT is not defined. We will guess")
    mtpy_path = os.path.abspath('../..')

if not os.path.isdir(mtpy_path):
    raise Exception("the guessed mtpy dir %s is not a folder!"% mtpy_path)

# change the variable below according to your edi files folder !!!
# edidir = r'C:/mtpywin/mtpy/data/edifiles'  # / is Unix and Win-Dos compatible
# or get this variable from the cmdline:  edidir = sys.argv[1]

edidir = os.path.join(mtpy_path,'data','edifiles2')

# savepath = r'C:\tmp'
temp_dir = tempfile.gettempdir()
print('Using temporary directory ' + temp_dir)
savepath = temp_dir


if not os.path.isdir(edidir):
    print ("Error: please provide the path to edi folder")
    sys.exit(1)

period_index_list = [0, 1, 10, 20, 30, 40, 50, 59]  # user to customise

# show three different kind of calculated pen-depth
pen2d.plot2Dprofile(edidir, period_index_list, 'det')

pen2d.plot2Dprofile(edidir, period_index_list, 'zxy')

pen2d.plot2Dprofile(edidir, period_index_list, 'zyx')


plt.savefig(os.path.join(savepath,'penetration_depth_profile.png'), # change to your preffered filename
            dpi=400) # change to your preferred file resolution

