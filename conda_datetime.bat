@echo off

:: Get the full year, month, day, hour, and minute
for /f "tokens=2 delims==" %%a in ('wmic os get localdatetime /value') do set datetime=%%a

:: Extract parts of the date and time
set year=%datetime:~2,2%
set month=%datetime:~4,2%
set day=%datetime:~6,2%
set hour=%datetime:~8,2%
set minute=%datetime:~10,2%

:: Combine them into the desired format
set datetimeString=%year%%month%%day%%hour%%minute%
