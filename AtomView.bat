@echo off
REM Double-click to start AtomView on Windows.
cd /d "%~dp0"
".venv\Scripts\python.exe" -m atomview %*
