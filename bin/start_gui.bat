setlocal

call activate mtpy
pushd "%~dp0"
cd ..
set mtpy_root="%CD%"
popd
pushd %mtpy_root%
set PYTHONPATH=%mtpy_root%;%PYTHONPATH%
python -OO %mtpy_root%\mtpy\gui\SmartMT\start.py
popd
call deactivate
PAUSE