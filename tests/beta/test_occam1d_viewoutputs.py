"""
Script to visualise the 1D model output.

"""
import os
import os.path as op
import numpy as np

import mtpy.modeling.occam1d as mtoc1d
from tests.beta import *


# import matplotlib.pyplot as plt
# plt.ion() # make figure disappear automatically
# plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc



def test_fun():
    """
    test function
    :return: T/F
    """
    # AK's working directory
    # idir = r'C:\Git\mtpy\examples\model_files\Occam1d'
    # savepath = r'C:\Git\mtpy\examples\plots\Occam1d'

    #FZ's workdir
    idir = os.path.join(SAMPLE_DIR,'Occam1d')
    savepath = TEST_TEMP_DIR

    # model and data file names
    modelfn=op.join(idir,'Model1D')
    datafn=op.join(idir,'Occam1d_DataFile_DET.dat')

    # list of all the 'smooth' files, exclude the log file
    fn_list = np.array([ff for ff in os.listdir(idir)])

    #go to model results directory and find the latest iteration file
    #iterfile = 'ITER_135.iter'
    #respfile = 'ITER_135.resp'
    iterfile = 'ITER_97.iter'
    respfile = 'ITER_97.resp'


    # get maximum iteration
    iterfn  = op.join(idir,iterfile)
    respfn  = op.join(idir,respfile)

    # read in the model, don't need these lines to view the output but useful if you want to analyse the data
    oc1m = mtoc1d.Model(model_fn=modelfn)
    oc1m.read_iter_file(iterfn)

    # read in the data file
    oc1d = mtoc1d.Data(data_fn=datafn)
    oc1d.read_data_file(data_fn=datafn)
    oc1d.read_resp_file(resp_fn=respfn,data_fn=datafn)

    # plot the model output
    pr = mtoc1d.Plot1DResponse(data_te_fn = datafn,
                               data_tm_fn = datafn,
                               model_fn = modelfn,
                               resp_te_fn = respfn,
                               iter_te_fn = iterfn,
                               resp_tm_fn = respfn,
                               iter_tm_fn = iterfn,
                               depth_limits = (0,1),

                               )
    pr.axm.set_xlim(1e-1,1e3)
    pr.axr.set_ylim(1,100)

    p2file = op.join(savepath, 'occam1dplot.png')
    pr.save_figure(p2file)

    assert (os.path.exists(p2file))

############################################
if __name__ == "__main__":
    test_fun()
