@echo Compiling .qrc files
@for %%f in (*.qrc) do @echo %%f -^> %%~nf_qt5_rc.py & @call pyrcc5 "%%f" -o "%%~nf_qt5_rc.py"
