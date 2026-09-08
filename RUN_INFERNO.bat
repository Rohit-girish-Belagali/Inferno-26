@echo off
REM Double-clickable launcher. A .ps1 cannot be run by double-clicking --
REM Windows opens it in Notepad instead -- but a .bat can, and this one
REM just starts PowerShell with the right flags.
REM
REM The OpenRouter key is read from openrouter_key.txt sitting next to
REM this file, if it exists. That keeps the key out of the repository
REM while still letting the whole thing run from one double-click. With
REM no key file the application still runs: the AI panel reports itself
REM unavailable and the eight predefined models solve with no API calls.
setlocal
cd /d "%~dp0"

set "KEY="
if exist "%~dp0openrouter_key.txt" (
  set /p KEY=<"%~dp0openrouter_key.txt"
)

if defined KEY (
  echo Starting Inferno with the AI copilot enabled.
  powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_windows.ps1" -ApiKey "%KEY%"
) else (
  echo No openrouter_key.txt found - starting Inferno without the AI copilot.
  echo The solver and all 8 example models still work.
  powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_windows.ps1"
)

echo.
echo Inferno has stopped. Press any key to close this window.
pause >nul
