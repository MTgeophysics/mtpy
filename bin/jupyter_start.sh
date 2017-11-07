#! /usr/bin/env bash

#  script to start a localhost jupyter notebook server, using anaconda2 python:

# If want to use C:/Anaconda2:
#    export PATH=/c/Anaconda2:/c/Anaconda2/Scripts:/c/Anaconda2/Library/bin:$PATH

# provide python lib search path for import the mtpy modules in scripts:
export  PYTHONPATH=/e/Githubz/mtpy2  #  path to mtpy2 (github repo) dir  
#export  PYTHONPATH=/Softlab/Githubz/mtpy2  # for import mtpy.module 



# start a python notebook server in background

jupyter notebook &> ./_jupyter.log &
