@echo off
setlocal
cd /d "%~dp0"
if not "%~1"=="" goto explicit
where py >nul 2>nul
if errorlevel 1 goto missing
py -3 install_molsimplify.py
goto done
:explicit
"%~1" install_molsimplify.py
goto done
:missing
echo Install Python 3.10-3.12 x64, or run install_molsimplify.bat "C:\path\python.exe"
pause
exit /b 1
:done
if errorlevel 1 pause
