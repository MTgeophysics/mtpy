#! /bin/bash

export MTPYPATH=/Softlab/Githubz/mtpy2
export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH


python -V

cd $MTPYPATH

python examples/modem_plotptresidual.py

python examples/plot_edis.py examples/data/edi_files/georgina
python examples/plot_edis.py examples/data/edi_files/

python examples/plot_phase_tensor_map.py examples/data/edi_files


python mtpy/imaging/plot_phase_tensor_maps_with_tipper.py

python mtpy/imaging/mtplottools.py

python mtpy/core/mt.py




