@echo Converting .ui files
@for %%f in (*.ui) do @echo %%f -^> %%~nf.py & @call pyuic4.bat "%%f" -o "%%~nf.py"

compile_resources_pyqt4.bat
