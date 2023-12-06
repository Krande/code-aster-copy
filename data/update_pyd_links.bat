@echo off
title Create testing link
mode con: cols=90 lines=30
Color 0F
cls
set aster_lib=%~dp0lib\aster
for %%A IN (aster,aster_core,aster_fonctions,libaster,med_aster) do (
    if exist "%aster_lib%\%%A.pyd" (
        del %aster_lib%\%%A.pyd
    )
mklink /h %aster_lib%\%%A.pyd %aster_lib%\aster.dll
)
exit
