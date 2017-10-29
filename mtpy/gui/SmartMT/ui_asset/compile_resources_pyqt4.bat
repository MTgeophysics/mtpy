@echo Compiling .qrc files
@for %%f in (*.qrc) do @echo %%f -^> %%~nf_py2_qt4_rc.py & @call pyrcc4 -py2 "%%f" -o "%%~nf_py2_qt4_rc.py"
@for %%f in (*.qrc) do @echo %%f -^> %%~nf_py3_qt5_rc.py & @call pyrcc4 -py3 "%%f" -o "%%~nf_py3_qt4_rc.py"
