#! /usr/bin/env bash
#######################################
# test-running all python scripts.
#######################################


cd mtpy
export PYTHONPATH=.
export GDAL_DATA=/c/Anaconda2/Library/share/gdal  # a directory containing many CSV files used by GDAL
# export PYTHONPATH=/e/Githubz/mtpy

python examples/plot_edis.py data/edifiles/15125A.edi

# New shape file creator script: create csv and shape files
python mtpy/utils/shapefiles_creator.py /e/Data/MT_Datasets/WenPingJiang_EDI /e/Data/MT_Datasets/WenPingJiang_SHP2 # long time
python mtpy/utils/shapefiles_creator.py /e/Data/MT_Datasets/WenPingJiang_EDI /e/tmp_20/
python mtpy/utils/shapefiles_creator.py /e/Data//MT_Datasets/GA_UA_edited_10s-10000s/ /e/Data//MT_Datasets/GA_UA_edited_10s-10000s_SHP2

# **** compare with the Old shape files generation: user please tune the __main__ section to provide suitable params
python mtpy/utils/shapefiles.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ /e/Data/MT_Datasets/GA_UA_edited_10s-10000s_SHP/ # not working after merge
python mtpy/utils/shapefiles.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM_SHP/

# EDI collection (surveys) properties: plot stations and create CSV files
python mtpy/core/edi_collection.py data/edifiles/ /e/tmp0
python mtpy/core/edi_collection.py examples/data/edi2/ /e/tmp0
python mtpy/core/edi_collection.py examples/data/edi_files /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/WenPingJiang_EDI/ /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/728889/EDI_files/ /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/75098/EDI_files/ /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/75099_Youanmi/EDI_Files_edited/BB_edi_edited/ /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/75099_Youanmi/EDI_Files_edited/BB_edi_edited/YM1/ /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/75099_Youanmi/EDI_Files_edited/BB_edi_edited/YM2 /e/tmp0
python mtpy/core/edi_collection.py /e/Data/MT_Datasets/75099_Youanmi/EDI_Files_edited/LP_edi_edited/YML123/ /e/tmp0


python mtpy/imaging/penetration_depth3d.py data/edifiles/ 10
python mtpy/imaging/penetration_depth3d.py data/edifiles/ 2.8571s
python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10

python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 0.0625s
python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 10s
python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 16s
python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 40s
python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 341.0s
python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 4096s
python mtpy/imaging/penetration_depth3d.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 2048s

# plot_phase_tensor_map.py has been moved to examples/scripts/ modified ( 21/02/2018 )

python examples/scripts/plot_phase_tensor_map.py  /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10 /e/MTPY2_Outputs

python examples/scripts/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.025  /e/MTPY2_Outputs
python examples/scripts/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.01  /e/MTPY2_Outputs
python examples/scripts/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0625 /e/MTPY2_Outputs/
python examples/scripts/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0005 /e/MTPY2_Outputs/
#python examples/scripts/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0004883 /e/MTPY2_Outputs/

#python examples/plot_phase_tensor_map.py  examples/data/edi_files 10 /e/MTPY2_Outputs/
#python examples/plot_phase_tensor_map.py  examples/data/edi2 10
#python examples/plot_phase_tensor_map.py  data/edifiles/ 10
#python examples/plot_phase_tensor_map.py  data/edifiles/ 10 /e/MTPY2_Outputs/
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10 /e/MTPY2_Outputs/
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10 /e/MTPY2_Outputs/ptmap3.jpg
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10 /e/MTPY2_Outputs/ptmap3deg.jpg
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0001
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.003
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0625
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0625 /e/MTPY2_Outputs/
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.1
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.1
#python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.1 /e/MTPY2_Outputs/
#python examples/plot_phase_tensor_map.py examples/data/edi2 10 /e/MTPY2_Outputs/
#python examples/plot_phase_tensor_map.py examples/data/edi2 10 /e/MTPY2_Outputs/abc.jpg
#python examples/plot_phase_tensor_map.py examples/data/edi2 10 /e/MTPY2_Outputs/ptmap2.jpg
#python examples/plot_phase_tensor_map.py data/edifiles 20
#python examples/plot_phase_tensor_map.py data/edifiles/ 10 /e/MTPY2_Outputs/ptmap3deg.png


# phase tensor pseudo sections
# plot_phase_tensor_section.py has been moved to examples/scripts/ modified ( 21/02/2018 )

python examples/scripts/plot_phase_tensor_section.py data/edifiles/
python examples/scripts/plot_phase_tensor_section.py examples/data/edi2
python examples/scripts/plot_phase_tensor_section.py examples/data/edi_files/
python examples/scripts/plot_phase_tensor_section.py examples/data/edi_files_2/

python examples/plot_phase_tensor_section.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/
python examples/plot_phase_tensor_section.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s  #change to stretch=(20, 40), y more

# generate inputs for MODEM, and output effective new_edi/ files inside the result folder.
# prepare input for ModEM input files: examples/create_modem_input.py
# create_modem_input.py has been moved to examples/scripts/ modified ( 21/02/2018 )

python examples/scripts/create_modem_input.py data/edifiles/ examples/etopo1.asc /e/tmp/modem_inputs/

python examples/scripts/create_modem_input.py  /e/Data/MT_Datasets/concurry_EDI_files/ /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc /e/tmp/
python examples/scripts/create_modem_input.py  /e/Data/MT_Datasets/GA_UA_edited_10s-10000s /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc /e/tmp/GA_UA_edited_10s-10000s_B/
python examples/scripts/create_modem_input.py  /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc /e/tmp/3D_MT_data_modemin/
python examples/scripts/create_modem_input.py  /e/Data/MT_Datasets/Isa_EDI_edited_10Hz_1000s/ /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc /e/tmp/Isa_EDI_edited_10Hz_1000s
python examples/scripts/create_modem_input.py  /e/Data/MT_Datasets/WenPingJiang_EDI /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc /e/tmp/WenPingJiang_EDI_modemin

python examples/create_modem_input.py  /e/Data/MT_Datasets/GA_UA_edited_10s-10000s /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc /e/tmp/GA_UA_edited_10s-10000s_16
python examples/create_modem_input.py /e/Data/MT_Datasets/concurry_EDI_files/ /e/Data/MT_Datasets/concurry_topo/AussieContinent_etopo1.asc /e/tmp/concurry_modem

# FZ: Todo 2017-09 Visualize ModEM outputs and create CSV files for the data
# modem_plot_models.py has been moved to examples/scripts/ modified ( 21/02/2018 )
#
python examples/scripts/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07
python examples/scripts/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07 RMSMap
python examples/scripts/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07 PTMap
python examples/scripts/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07 Response
python examples/scripts/modem_plot_models.py ./examples/data/ModEM_files/VicSynthetic07 DepthSlice

python examples/scripts/modem_plot_models.py /e/Data/Modeling/Isa/100hs_flat_BB/
python examples/scripts/modem_plot_models.py /e/Data/Modeling/Isa/100hs_flat_BB/ Response
python examples/scripts/modem_plot_models.py /e/Data/Modeling/Isa/100hs_flat_BB/ DepthSlice

# view horizontal slices of a rho file
# ( 21/02/2018 )
# plot_depth_slice.py exist in two directories (/mtpy/imaging/plot_depth_slice.py/mtpy/modeling/modem/plot_depth_slice.py)
#
python mtpy/imaging/plot_depth_slice.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho 20
python mtpy/imaging/plot_depth_slice.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Model.ws


# View all or selected multiple horizontal slices of an inversion model (output of ModEM)
python mtpy/imaging/modem_plot_slices.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho z 300
python mtpy/imaging/modem_plot_slices.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho ns 300
python mtpy/imaging/modem_plot_slices.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho ew 300

# View all or selected horizontal slices of an initial model
python mtpy/imaging/modem_plot_slices.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat  /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Model.ws ew
python mtpy/imaging/modem_plot_slices.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat  /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Model.ws z -1000 0 100
python mtpy/imaging/modem_plot_slices.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Model.ws z -1000 1000


# create CSV files
# Following scripts were not found in the repository ( 20/02/2018 )
#python mtpy/modeling/modem_output_to_views.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho 20
#python mtpy/modeling/modem_output_to_views.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat  /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Model.ws -1000 1000

# create csv files from modem.dat
# New loactions for the following scripts ( 20/02/2018 )
#python mtpy/modeling/modem_data_to_phase_tensor.py examples/data/ModEM_files/Modular_MPI_NLCG_028.dat [OutDir]
#python mtpy/modeling/modem_data_to_phase_tensor.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat [OutDir]
python examples/ecripts/modem_data_to_phase_tensor.py examples/data/ModEM_files/Modular_MPI_NLCG_028.dat [OutDir]
python examples/ecripts/modem_data_to_phase_tensor.py /e/MTPY2_Outputs/GA_UA_edited_10s-10000s_modem_inputs/ModEM_Data.dat [OutDir]

