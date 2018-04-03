"""
# ==============================================================================
Plots data vs model response computed by an Inversion Method

Examples



# ==============================================================================
"""

import click
from mtpy.modeling.modem.plot_response import PlotResponse
from mtpy.mtpy_globals import *



stem_data_file = 'Modular_MPI_NLCG_004.dat'
input_data_file = 'ModEM_Data.dat'
directory = r'C:\mtpywin\mtpy\examples\model_files\ModEM_2'
collection_station = 'Synth00'
font_size = 6
    
ro = PlotResponse(data_fn=os.path.join(directory, input_data_file),
                  resp_fn=os.path.join(directory, stem_data_file),
                  plot_type=[collection_station],
                  plot_style=2,
                  plot_z=True,
                  font_size=font_size)
ro.plot()

