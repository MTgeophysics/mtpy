# -*- coding: utf-8 -*-
"""
Created on Mon Feb 02 11:38:58 2015

@author: Alison Kirkby

create modem input files

"""

import os.path as op
import os
os.chdir(r'C:/mtpywin/mtpy') # change to the directory where your mtpy is installed to
                             # ensure you are using the correct mtpy


from mtpy.modeling.modem import Model
import numpy as np

# path to save to
wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'

# provide centre position of model in real world coordinates (eastings/northings)
centre = np.array([0., 0., 0.])

# modem .rho file
iterfn = 'Modular_MPI_NLCG_004.rho'

moo = Model(model_fn=op.join(wd,iterfn))
moo.read_model_file()
moo.write_gocad_sgrid_file(origin=centre,
                           fn = r'C:\test\ModEM\ModEM_output.sg' # filename to save to, optional
                                                                 # saves in model directory if not provided
                           )