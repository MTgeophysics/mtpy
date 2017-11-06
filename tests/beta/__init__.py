
# package tests.beta scope global params
import os
import matplotlib.pyplot as plt
plt.ion() # all beta/test_ plot disappear automatically

# MTPY_ROOT='/Softlab/Githubz/mtpy'    # source code root dir
MTPY_ROOT='E:/Githubz/mtpy'    # source code root dir

EDI_DATA_DIR = os.path.join(MTPY_ROOT,'examples/data/edi_files')
EDI_DATA_DIR2 = os.path.join(MTPY_ROOT,'examples/data/edi_files_2')

AUS_TOPO_FILE = os.path.join(MTPY_ROOT,'examples/data/AussieContinent_etopo1.asc')

# path to directory containing model input files - samples reference for compare
SAMPLE_DIR = os.path.join(MTPY_ROOT,'examples/model_files') # r'E:\Githubz\mtpy\examples\model_files'

# test runs output directory where to save plots/files to
TEMP_OUT_DIR = os.path.join(MTPY_ROOT,'temp/beta_out_dir')  # r'E:\Githubz\mtpy\temp'
