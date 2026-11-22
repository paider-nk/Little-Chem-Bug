@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo First start start.bat once to install the application dependencies.
  pause
  exit /b 1
)
if "%~1"=="" goto gui
".venv\Scripts\python.exe" gaussian_workflow.py %*
if errorlevel 1 pause
exit /b %errorlevel%

:gui
".venv\Scripts\python.exe" workflow_gui.py
if errorlevel 1 pause
exit /b %errorlevel%
