@echo off
title CLAIREscope Server Controller
cd /d "%~dp0"

echo ===================================================
echo     Launching CLAIREscope Server Controller...
echo ===================================================

:: Try local virtual environment first
if exist ".venv\Scripts\python.exe" (
    start "" ".venv\Scripts\python.exe" launcher.py
    exit /b 0
)

if exist "..\.venv\Scripts\python.exe" (
    start "" "..\.venv\Scripts\python.exe" launcher.py
    exit /b 0
)

:: Fallback to system Python
python launcher.py
