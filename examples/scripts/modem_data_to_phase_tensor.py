#!/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Compute Phase Tensors from ModEM Dat File and output to CSV files

Usage Examples:
    python scripts/modem_data_to_phase_tensor.py examples/data/ModEM_files/Modular_MPI_NLCG_028.dat [OutDir]
    python scripts/modem_data_to_phase_tensor.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat [OutDir]


Developer:      fei.zhang@ga.gov.au
LastUpdate:     08/09/2017
LastUpdate:     05/12/2017 FZ moved the function into the module mtpy.modeling.modem.Data
LastUpdate:     21/02/2018 Added command for running the script
"""

import sys, os
from mtpy.modeling.modem import Data
from mtpy.mtpy_globals import NEW_TEMP_DIR
import click


if __name__ == "__main__old":

    file_dat = sys.argv[1]
    if len(sys.argv)>2:
        outdir = sys.argv[2]
    else:
        outdir=NEW_TEMP_DIR

    obj = Data()
    
    obj.compute_phase_tensor(file_dat, outdir)

# =============================================================================================
# Command line wrapper for processing phase tensors and output to csv file
# =============================================================================================

@click.command()
@click.option('-i','--dat_file',type=str,
              default='examples/data/ModEM_files/Modular_MPI_NLCG_028.dat', \
              help='input path/datafile')
@click.option('-o','--output_dir',type=str,default="temp",help='Output directory')
def process_phase_tensors(dat_file,output_dir):
    print ("Input path/datfile   --------->     {}".format(dat_file))
    print ("Output directory     --------->     {}".format(output_dir))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if os.path.isfile(dat_file):
        obj = Data()
        obj.compute_phase_tensor(dat_file, output_dir)
    else:
        print("Please provide an input dat file !")

if __name__ == '__main__':
    process_phase_tensors()
