@echo Converting .ui files
@for %%f in (*.ui) do @echo %%f -^> %%~nf.py & @call pyuic4.bat "%%f" -o "%%~nf.py"

@echo Converting .qrc files
@for %%f in (*.qrc) do @echo %%f -^> %%~nf_rc.py & @call pyrcc4 -py2 "%%f" -o "%%~nf_rc.py"
