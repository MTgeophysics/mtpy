#! /usr/bin/env bash
#######################################
# commandline for test-runing all plots
#######################################
#WinPC: python examples/plot_phase_tensor_map.py  E:/Datasets/MT_Datasets/GA_UA_edited_10s-10000s 0.0625 E:/MTPY2_Outputs/

python examples/Modem_buildinputfiles.py

python examples/plot_edis.py tests/data/edifiles/15125A.edi


python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10
python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 0.0625s
python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 10s
python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 16s
python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 40s
python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 341.0s
python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 4096s
python mtpy/imaging/penetration_depth_3d_profile.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 2048s
python mtpy/imaging/penetration_depth_3d_profile.py tests/data/edifiles/ 10
python mtpy/imaging/penetration_depth_3d_profile.py tests/data/edifiles/ 2.8571s

python examples/plot_phase_tensor_map.py  /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/ 10 /e/MTPY2_Outputs

python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 0.025  /e/MTPY2_Outputs
python examples/plot_phase_tensor_map.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s 0.01  /e/MTPY2_Outputs
python examples/plot_phase_tensor_map.py  /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0625 /e/MTPY2_Outputs/
python examples/plot_phase_tensor_map.py  /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ 0.0004883 /e/MTPY2_Outputs/

#python examples/plot_phase_tensor_map.py  examples/data/edi_files 10 /e/MTPY2_Outputs/
#python examples/plot_phase_tensor_map.py  examples/data/edi2 10
#python examples/plot_phase_tensor_map.py  tests/data/edifiles/ 10
#python examples/plot_phase_tensor_map.py  tests/data/edifiles/ 10 /e/MTPY2_Outputs/
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
#python examples/plot_phase_tensor_map.py tests/data/edifiles 20
#python examples/plot_phase_tensor_map.py tests/data/edifiles/ 10 /e/MTPY2_Outputs/ptmap3deg.png


# phase tensor pseudo sections
python examples/plot_phase_tensor_section.py tests/data/edifiles/
python examples/plot_phase_tensor_section.py examples/data/edi2
python examples/plot_phase_tensor_section.py examples/data/edi_files/
python examples/plot_phase_tensor_section.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM/
python examples/plot_phase_tensor_section.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s  #change to stretch=(20, 40), y more


# shape files generation: fine-tune the __main__ section about calling params
python mtpy/utils/shapefiles.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s/ /e/Data/MT_Datasets/GA_UA_edited_10s-10000s_SHP/
python mtpy/utils/shapefiles.py /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM /e/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM_SHP/

