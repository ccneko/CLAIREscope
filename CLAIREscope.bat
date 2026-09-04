@echo off
title CLAIREscope Server Controller
cd /d "%~dp0"

echo ===================================================
echo     Launching CLAIREscope Server Controller...
echo ===================================================

:: 1. Try local Windows virtual environment if present
if exist ".venv\Scripts\pythonw.exe" (
    start "" ".venv\Scripts\pythonw.exe" launcher.py
    exit /b 0
)
if exist ".venv\Scripts\python.exe" (
    start "" ".venv\Scripts\python.exe" launcher.py
    exit /b 0
)

:: 2. Try WSL managed bioinfo environment
where wsl >nul 2>nul
if %ERRORLEVEL% equ 0 (
    wsl /home/claire/Software/pyenvs/bioinfo/.venv/bin/python /home/claire/dev/CLAIREscope/launcher.py
    exit /b 0
)

:: 3. Fallback to Windows system Python
python launcher.py
