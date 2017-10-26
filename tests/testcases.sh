#! /bin/bash

# env setup --- please modify accordingly

# export MTPYPATH=/Softlab/Githubz/mtpy2
# export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH
# OR cd /e/Githubz/mtpy2 (develop)
# export PYTHONPATH=.
#windows git-bash: export  PYTHONPATH=/e/Githubz/mtpy2


echo "------------Be mindful which python interpreter? ---------"
which python
python -V

cd $MTPYPATH

echo "-----------in mtpy2 develop branch---------"

python examples/modem_plotptresidual.py

# plot the response from edi files
python examples/plot_edis.py examples/data/edi_files/
python examples/plot_edis.py examples/data/edi_files/georgina

# phase tensor map
python examples/plot_phase_tensor_map.py tests/data/edifiles/ 1
python examples/plot_phase_tensor_map.py examples/data/edi_files/georgina 1

#phase tensor psuedo section
python examples/plot_phase_tensor_section.py examples/data/edi_files
python examples/plot_phase_tensor_section.py examples/data/edi_files/georgina
python examples/plot_phase_tensor_section.py tests/data/edifiles/


echo "Finished running the script $0"

