#! /bin/bash

# env set up --- please modify accordingly

export MTPYPATH=/Softlab/Githubz/mtpy2
export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH


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

echo "Finished running the script $0"






