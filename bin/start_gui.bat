@ECHO off
setlocal

SET mtpyenv=mtpy
IF NOT "%1"=="" SET mtpyenv="%1"
ECHO calling activate %mtpyenv%
call activate %mtpyenv%
pushd "%~dp0"
cd ..
set mtpy_root=%CD%
popd
pushd %mtpy_root%
set PYTHONPATH=%mtpy_root%;

ECHO mtpy_root is %mtpy_root%
ECHO PYTHONPATH is %PYTHONPATH%

python -OO %mtpy_root%\mtpy\gui\SmartMT\start.py
popd
call deactivate
PAUSE
