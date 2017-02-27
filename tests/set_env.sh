# Setup script to use anaconda2 in Windows PC bash shell:
#    source set_env.sh

# use C:/Anaconda2:
    export PATH=/c/Anaconda2:/c/Anaconda2/Scripts:/c/Anaconda2/Library/bin:$PATH
    export  PYTHONPATH=/e/Githubz/mtpy2  # use mtpy module in export

# Test:
#   cd /e/Githubz/mtpy2
#   python mtpy/imaging/penetration_depth_3d_profile.py tests/data/edifiles/ 50

# start a python notebook server in background
#	jupyter notebook &> localdir/jupyter.log &
