"""
# ==============================================================================
Plots data vs model response computed by an Inversion Method

Examples

python examples/cmdline/plot_response.py
python examples/cmdline/plot_response.py -d /e/Githubz/example_plot_response -s Modular_MPI_NLCG_094.dat -i  ModEM_Data.dat -c GB08 -p False
python examples/cmdline/plot_response.py -d /e/Githubz/example_plot_response -s Modular_MPI_NLCG_094.dat -i  ModEM_Data.dat -c GB08 -p True
python examples/cmdline/plot_response.py -d /e/Githubz/example_plot_response -s Modular_MPI_NLCG_094.dat -i  ModEM_Data.dat -c GB09 -p True
python examples/cmdline/plot_response.py -d /e/Githubz/example_plot_response -s Modular_MPI_NLCG_094.dat -i  ModEM_Data.dat -c GB09 -p False

See also
python mtpy/modeling/modem/plot_response.py 2col
python mtpy/modeling/modem/plot_response.py

# ==============================================================================
"""

import click
from mtpy.modeling.modem.plot_response import PlotResponse
from mtpy.mtpy_globals import *

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-d','--directory',type=str,default=r'examples/model_files/ModEM_2',help='directory for data files')
@click.option('-s','--stem_data_file',type=str,default='Modular_MPI_NLCG_004.dat', help='file stem')
@click.option('-i','--input_data_file',type=str,default='ModEM_Data.dat', help='Data File')
@click.option('-c','--collection_station',type=str,default='Synth02', help='Data Collection station')
@click.option('-p','--plot_z',type=bool,default=False, help=
                            '[True | False ] Plot True for Impedence, False for Resistivity and Phsse')
@click.option('-f','--font_size',type=int,default=3, help='Plot Text Fond Size ')
def merge_plotting(directory, stem_data_file, input_data_file, collection_station,plot_z, font_size):

    print("============================================================================")
    print("")
    print("Following are the examples for running plot_response : ")
    print("")
    print("python examples/cmdline/plot_response.py  [--help | -h ]")
    print("python examples/cmdline/plot_response.py")
    print("python examples/cmdline/plot_response.py -d examples/model_files/ModEM_2 " +
          "-s Modular_MPI_NLCG_004.dat -i ModEM_Data.dat -c Synth02 -p False -f 3")
    print("python examples/cmdline/plot_response.py -d examples/data/ModEM_files -p True ")
    print("")
    print("============================================================================")


    ro = PlotResponse(data_fn=os.path.join(directory, input_data_file),
                      resp_fn=os.path.join(directory, stem_data_file),
                      plot_type=[collection_station],
                      plot_style=2,
                      plot_z=plot_z,
                      font_size=font_size)
    ro.plot_2col()

if __name__ == "__main__":

    merge_plotting()
