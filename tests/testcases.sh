#! /bin/bash

# env setup --- please modify accordingly

export MTPYPATH=/Softlab/Githubz/mtpy2
export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH
#windows git-bash: export  PYTHONPATH=/e/Githubz/mtpy2


echo "------------Be mindful which python interpreter? ---------"
which python
python -V

cd $MTPYPATH

echo "-----------Test running Alison's scripts---------"
echo "-----------in mtpy2 develop branch---------"

python examples/modem_plotptresidual.py

python examples/plot_edis.py examples/data/edi_files/
python examples/plot_edis.py examples/data/edi_files/georgina

python examples/plot_phase_tensor_map.py examples/data/edi_files
python examples/plot_phase_tensor_map.py examples/data/edi_files/georgina

python examples/plot_phase_tensor_section.py examples/data/edi_files
python examples/plot_phase_tensor_section.py examples/data/edi_files/georgina


####################################
# Change log
# 2016-12-01

# PTMap works by a single index instead of loop-over all [0,1,2,,...18]
# could not save figure, but now can save figure.
# python examples/modem_plotmodel2.py examples/data/ModEM_files/VicSynthetic07 PTMap 8

# changed ticks default, colorbar better positioned, depth_index => depth in KM? how to get values?
# one plot each time for a given depth tested OK.  All plots (0-80) not sure.
# python examples/modem_plotmodel2.py ./examples/data/ModEM_files/VicSynthetic07 DepthSlice

# save figure not working, Now OK.
# python examples/modem_plotmodel2.py ./examples/data/ModEM_files/VicSynthetic07 Response

# will loop over all periods,
# save to a png file for each period.
# python examples/modem_plotmodel2.py ./examples/data/ModEM_files/VicSynthetic07 RMSMap

# fix github issue #12: plt show in commandline run.
# python examples/modem_plotmodel_vertical.py

python mtpy/imaging/modem_plot_vertical_slice.py examples/data/ModEM_files/VicSynthetic07/

echo "Finished running the script $0"

