param(
    [string]$PythonPath = "",
    [ValidateRange(1024, 65535)][int]$Port = 8501,
    [switch]$Live,
    [switch]$SkipInstall
)

$ErrorActionPreference = "Stop"
$taskRepoRoot = $PSScriptRoot
$taskVenvPython = Join-Path $taskRepoRoot ".venv\Scripts\python.exe"
$taskPreparedPackages = Join-Path $taskRepoRoot ".runtime\python"
$taskBundledPython = Join-Path ([Environment]::GetFolderPath("UserProfile")) ".cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"

# Reuse the dependencies prepared for this local demonstration when available.
# The portable release excludes .runtime and follows the normal venv path below.
if ([string]::IsNullOrWhiteSpace($PythonPath) -and
    (Test-Path -LiteralPath (Join-Path $taskPreparedPackages "streamlit\__init__.py")) -and
    (Test-Path -LiteralPath $taskBundledPython)) {
    $taskPreviousPythonPath = $env:PYTHONPATH
    Push-Location -LiteralPath $taskRepoRoot
    try {
        $env:PYTHONPATH = $taskPreparedPackages
        $env:RE4_ENABLE_LIVE_EXECUTION = if ($Live) { "1" } else { "0" }
        Write-Host "Re4 workspace: http://localhost:$Port (prepared local dependencies)" -ForegroundColor Cyan
        & $taskBundledPython -m streamlit run streamlit_app.py --server.port $Port --server.address 127.0.0.1
    } finally {
        $env:PYTHONPATH = $taskPreviousPythonPath
        Pop-Location
    }
    return
}

if (-not (Test-Path -LiteralPath $taskVenvPython)) {
    if ([string]::IsNullOrWhiteSpace($PythonPath)) {
        if (Get-Command py -ErrorAction SilentlyContinue) {
            $PythonPath = (& py -3.12 -c "import sys; print(sys.executable)").Trim()
        } elseif (Get-Command python -ErrorAction SilentlyContinue) {
            $PythonPath = (Get-Command python).Source
        } else {
            throw "Install Python 3.12, or pass -PythonPath C:\path\to\python.exe."
        }
    }
    & $PythonPath -m venv (Join-Path $taskRepoRoot ".venv")
    if ($LASTEXITCODE -ne 0) { throw "Could not create the Python virtual environment." }
}

Push-Location -LiteralPath $taskRepoRoot
try {
    if (-not $SkipInstall) {
        $taskRequirements = if ($Live) { "requirements-live.txt" } else { "requirements.txt" }
        & $taskVenvPython -m pip install -r $taskRequirements
        if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed." }
    }
    $env:RE4_ENABLE_LIVE_EXECUTION = if ($Live) { "1" } else { "0" }
    Write-Host "Re4 workspace: http://localhost:$Port" -ForegroundColor Cyan
    & $taskVenvPython -m streamlit run streamlit_app.py --server.port $Port --server.address 127.0.0.1
} finally {
    Pop-Location
}
