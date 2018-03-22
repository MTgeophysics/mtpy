
import os
import os.path as op
os.chdir(r'C:\mtpywin\mtpy')

from mtpy.modeling.modem import PlotSlices
from mtpy.modeling.modem import Data

wd = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'

model_fn = op.join(wd,'Modular_MPI_NLCG_004.rho')
data_fn = op.join(wd,'ModEM_Data.dat')



pObj = PlotSlices(model_fn=model_fn,
                  data_fn=data_fn,
                  save_path=wd,
                  cmap='seismic_r',
                  plot_stations=True,
                  font_size=6,
                  fig_size=(6,3)
                  )
figs,fpaths = pObj.export_slices(plane='N-E',
                                 indexlist=[30],
                                 station_buffer=20e3,
                                 save=False
                                 )#range(20,40))