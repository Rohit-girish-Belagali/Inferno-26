<#
.SYNOPSIS
    One command: install prerequisites, build Inferno, and run the copilot.

.DESCRIPTION
    Installs Python, CMake and the MSVC C++ toolchain via winget (skipping
    anything already present), builds the solver, then starts the server
    and opens a browser.

    The OpenRouter key is a PARAMETER, never a value stored in this file.
    A key committed to a repository is a key that has to be rotated, and
    this script is committed. Without a key the application still runs --
    the AI panel reports itself unavailable and the eight predefined
    models solve with zero API calls.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\run_windows.ps1 -ApiKey "sk-or-v1-..."

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File .\run_windows.ps1
#>
param(
    [string]$ApiKey = "",
    [string]$TextModel = "nvidia/nemotron-3-super-120b-a12b:free",
    [int]$Port = 8420,
    [switch]$SkipInstall
)

$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

function Say($msg, $color = "Cyan") { Write-Host "==> $msg" -ForegroundColor $color }
function Have($cmd) { $null -ne (Get-Command $cmd -ErrorAction SilentlyContinue) }

function Refresh-Path {
    $env:Path = [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" +
                [Environment]::GetEnvironmentVariable("Path", "User")
}

# ---------------------------------------------------------------- prereqs
if (-not $SkipInstall) {
    if (-not (Have winget)) {
        Write-Host @"
winget is not available on this machine.

Install "App Installer" from the Microsoft Store, then re-run this script.
Or install these three by hand and re-run with -SkipInstall:
  * Python 3      https://www.python.org/downloads/  (tick "Add to PATH")
  * CMake         https://cmake.org/download/
  * VS Build Tools with the "Desktop development with C++" workload
                  https://visualstudio.microsoft.com/downloads/
"@ -ForegroundColor Yellow
        exit 1
    }

    if (Have python) { Say "Python already installed - skipping" Green }
    else {
        Say "Installing Python 3.12 (~30 MB)"
        winget install --id Python.Python.3.12 -e --accept-package-agreements `
                       --accept-source-agreements --silent
        Refresh-Path
    }

    if (Have cmake) { Say "CMake already installed - skipping" Green }
    else {
        Say "Installing CMake (~30 MB)"
        winget install --id Kitware.CMake -e --accept-package-agreements `
                       --accept-source-agreements --silent
        Refresh-Path
    }

    # The C++ toolchain is the big one. Several GB, and there is no way
    # around it: this project compiles a real solver, not a script.
    $vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
    $haveCxx = (Test-Path $vswhere) -and
               ((& $vswhere -latest -products * `
                   -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
                   -property installationPath) -ne $null)
    if ($haveCxx) { Say "MSVC C++ toolchain already installed - skipping" Green }
    else {
        Say "Installing Visual Studio Build Tools with the C++ workload."
        Say "This is 3-6 GB and can take 10-20 minutes. Leave it running." Yellow
        winget install --id Microsoft.VisualStudio.2022.BuildTools -e `
            --accept-package-agreements --accept-source-agreements `
            --override "--quiet --wait --norestart --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
        Refresh-Path
    }
}

Refresh-Path
foreach ($tool in @("python", "cmake")) {
    if (-not (Have $tool)) {
        Write-Host "'$tool' is still not on PATH. Close this window, open a NEW PowerShell, and re-run." -ForegroundColor Red
        exit 1
    }
}

# ------------------------------------------------------------------ build
# A build/ directory copied from another machine is worse than no build
# directory: CMakeCache.txt records the absolute path it was generated in,
# so CMake stops with an error about the cache belonging elsewhere. Copying
# the project folder off a Mac or another PC is the normal way to get it
# here, which makes this the normal first failure. Clear it and move on.
$cache = Join-Path $PSScriptRoot "build\CMakeCache.txt"
if (Test-Path $cache) {
    $recorded = (Select-String -Path $cache -Pattern "^CMAKE_CACHEFILE_DIR:INTERNAL=(.*)$" |
                 Select-Object -First 1).Matches.Groups[1].Value
    $expected = (Join-Path $PSScriptRoot "build") -replace "\\", "/"
    if ($recorded -and ($recorded.TrimEnd('/') -ne $expected.TrimEnd('/'))) {
        Say "build\ was generated on another machine ($recorded) - removing it" Yellow
        Remove-Item -Path (Join-Path $PSScriptRoot "build") -Recurse -Force
    }
}

Say "Configuring"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
if ($LASTEXITCODE -ne 0) { Write-Host "CMake configure failed." -ForegroundColor Red; exit 1 }

Say "Building the solver"
cmake --build build --config Release
if ($LASTEXITCODE -ne 0) { Write-Host "Build failed." -ForegroundColor Red; exit 1 }

$dll = Get-ChildItem -Path build -Filter inferno.dll -Recurse -ErrorAction SilentlyContinue |
       Select-Object -First 1
if ($null -eq $dll) {
    Write-Host "Built, but inferno.dll was not produced. Cannot continue." -ForegroundColor Red
    exit 1
}
$env:INFERNO_LIBRARY = $dll.FullName
Say "Solver library: $($dll.FullName)" Green

Say "Running the test suite"
ctest --test-dir build -C Release --output-on-failure
if ($LASTEXITCODE -ne 0) {
    Write-Host "Some tests failed. Starting anyway - note what failed." -ForegroundColor Yellow
}

# ------------------------------------------------------------------- run
if ($ApiKey -ne "") {
    $env:OPENROUTER_API_KEY = $ApiKey
    $env:OPENROUTER_TEXT_MODEL = $TextModel
    $env:OPENROUTER_TIMEOUT_SECONDS = "180"
    $env:OPENROUTER_MAX_ATTEMPTS = "2"
    Say "AI copilot enabled - model $TextModel" Green
} else {
    Say "No API key given. AI is off; the solver and all 8 examples still work." Yellow
}

Say "Starting on http://127.0.0.1:$Port/  (Ctrl+C to stop)" Green
Start-Job -ScriptBlock {
    param($p) Start-Sleep -Seconds 3; Start-Process "http://127.0.0.1:$p/"
} -ArgumentList $Port | Out-Null

python server\app.py --port $Port
