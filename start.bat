@echo off
setlocal
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" goto check_tk
echo Creating a local Python environment...
where py >nul 2>nul
if errorlevel 1 goto use_python
py -3 -m venv .venv
if errorlevel 1 goto failed
goto check_tk
:use_python
where python >nul 2>nul
if errorlevel 1 goto no_python
python -m venv .venv
if errorlevel 1 goto failed
:check_tk
".venv\Scripts\python.exe" -c "import tkinter" >nul 2>nul
if errorlevel 1 goto no_tk
".venv\Scripts\python.exe" -c "import rdkit, numpy, matplotlib, scipy, ase, gemmi" >nul 2>nul
if not errorlevel 1 goto run
echo Installing application dependencies. Internet access is required on first launch...
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto failed
:run
".venv\Scripts\python.exe" gaussian_input.py %*
if errorlevel 1 goto failed
exit /b 0
:no_python
echo Python was not found. Install Python 3.10-3.13 from https://www.python.org/downloads/
echo Enable the Python launcher and Tcl/Tk during installation.
goto failed
:no_tk
echo Tcl/Tk is missing. Install a standard Python distribution with Tcl/Tk enabled.
:failed
echo The program could not finish. See the message above and README.md.
pause
exit /b 1
