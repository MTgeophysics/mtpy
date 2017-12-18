@ECHO off
setlocal

pushd "%~dp0"
cd ..
set mtpy_root="%CD%"
popd
pushd %mtpy_root%

:Home
CLS
ECHO 1.Setup new conda environment for MTpy
ECHO 2.Install conda
ECHO 3.Quit
ECHO.

CHOICE /C 123 /M "Enter your choice:"

:: Note - list ERRORLEVELS in decreasing order
IF ERRORLEVEL 3 GOTO Quit
IF ERRORLEVEL 2 GOTO InstallConda
IF ERRORLEVEL 1 GOTO Setup

:Setup
ECHO Setup new conda environment for MTpy
:: check conda
where conda exe >nul 2>nul
IF %ERRORLEVEL%==1 (
    ECHO conda not found in path. Please check your conda installation.
    ECHO If conda is freshly installed and added to PATH, you may need to log-off and log back in.
    GOTO End
)
SET mtpyenv=mtpy
SET INPUTSTRING=
SET /p INPUTSTRING="Enter conda environment name for MTpy (Default: %mtpyenv%):"
IF NOT "%INPUTSTRING%" == "" SET mtpyenv=%INPUTSTRING%

:: add channel: conda-forge
conda config --add channels conda-forge
conda create -y --name %mtpyenv% python=2 --file "%mtpy_root%"\requirements.txt
conda install -y pyqt=5 -n %mtpyenv%

GOTO End

:InstallConda
ECHO Install conda
ECHO We will redirect you to https://conda.io/miniconda.html and please download the installer for your system from the webpage and install conda.
PAUSE
START https://conda.io/miniconda.html
::SET /p condaInstaller="Please enter the path to the installer or drag and drop the installer to this window: "
::ECHO Installing conda using "%condaInstaller%"
::"%condaInstaller%" /InstallationType=JustMe /AddToPath=1 /RegisterPython=0
GOTO End

:Quit
GOTO End

:End
popd
PAUSE