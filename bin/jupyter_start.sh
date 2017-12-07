#! /usr/bin/env bash

#  script to start a localhost jupyter notebook server, using anaconda2 python:

# provide python lib search path for import the mtpy modules in scripts:
export  PYTHONPATH=/e/Githubz/mtpy  #  path to mtpy (github repo) dir


# start a python notebook server in background

jupyter notebook &> ./_jupyter.log &
