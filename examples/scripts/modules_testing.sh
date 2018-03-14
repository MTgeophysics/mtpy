#!/bin/bash

cd $MTPYPATH
echo `pwd`
export MTPYPYTHON="${HOME}/mtpy"
cd ${MTPYPYTHON}

if [[ -d "${MTPYPYTHON}/temp" ]]; then
    rm -f "${MTPYPYTHON}/temp/*"
    echo "Cleaned the tmp directory .... " 
fi

#----------------------------------------------
# Testing plot_edsi
#----------------------------------------------
echo "edi files plotting   ........ "
python examples/scripts/plot_edis.py

#----------------------------------------------
# Testing edi_collection
#----------------------------------------------
echo "edi_collection       ........ "
python mtpy/core/edi_collection.py

#----------------------------------------------
# Testing modem_data_to_phase_tensor
#----------------------------------------------
echo "modem_data_to_phase_tensor Testig ........ "
python examples/scripts/modem_data_to_phase_tensor.py


#----------------------------------------------
# Testing penetration depth 1d
#----------------------------------------------
echo "penetration depth testing 1d ........ "
python mtpy/imaging/penetration_depth1d.py


#----------------------------------------------
# Testing penetration depth 2d
#----------------------------------------------
echo "penetration depth testing 1d ........ "
python mtpy/imaging/penetration_depth2d.py


#----------------------------------------------
# Testing plot response
#----------------------------------------------
echo "plot response ........ "
python mtpy/imaging/plot_response.py


#----------------------------------------------
# Testing PlotMcmcResults 
#----------------------------------------------
echo "PlotMcmcResults ............ "
python mtpy/imaging2/PlotMcmcResults.py


#----------------------------------------------
# Testing concatenate input
#----------------------------------------------
echo "concatenate input  ........ "
python mtpy/utils/concatenate_input.py


#----------------------------------------------
# Testing shape files generation
#----------------------------------------------
echo "Shape files genetation  ........ "
python mtpy/utils/shapefiles.py


#----------------------------------------------
# Testing shape files creator
#----------------------------------------------
echo "shaope files creator  ........ "
python mtpy/utils/shapefiles_creator.py

#----------------------------------------------
# Testing Creating Modem Input
#----------------------------------------------
echo "Modem Input           ........ "
python examples/scripts/create_modem_input.py

