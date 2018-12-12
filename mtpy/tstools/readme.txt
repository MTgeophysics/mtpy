Windows:

>conda update conda
>conda create -c conda-forge --name tstools python=3.6 pyqt pyasdf matplotlib
>source activate tstools
>python example.py
>conda deactivate

Linux:

$conda update conda
$conda env create -f environment.yml
$conda activate tstools
$python example.py
$conda deactivate
