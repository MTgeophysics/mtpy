#! /bin/bash

export MTPYPATH=/Softlab/Githubz/mtpy2
export PYTHONPATH=/Softlab/Githubz/mtpy2:$PYTHONPATH

echo "------------ which python interpreter? ---------"
which python
python -V

cd $MTPYPATH

echo "----------- Running Alison's scripts---------"

python examples/modem_plotptresidual.py

python examples/plot_edis.py examples/data/edi_files/
python examples/plot_edis.py examples/data/edi_files/georgina

python examples/plot_phase_tensor_map.py examples/data/edi_files
python examples/plot_phase_tensor_map.py examples/data/edi_files/georgina

python examples/plot_phase_tensor_section.py examples/data/edi_files
python examples/plot_phase_tensor_section.py examples/data/edi_files/georgina

#------------------------------
echo "---------- Running mtpy1 scripts---------"

python mtpy/imaging/plot_phase_tensor_maps_with_tipper.py

python mtpy/core/mt.py
# replaced by
python testz/test_mt.py examples/data/edi_files/





