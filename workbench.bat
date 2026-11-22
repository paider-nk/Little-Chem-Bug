@echo off
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
  ".venv\Scripts\python.exe" chem_workbench.py
) else (
  python chem_workbench.py
)
if errorlevel 1 pause
