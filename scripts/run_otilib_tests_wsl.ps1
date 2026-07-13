<#
.SYNOPSIS
    Run the genuine OTILib tests inside WSL, from Windows. One copy-paste command.

.DESCRIPTION
    OTILib (GPLv3, external) only builds on Linux, so on Windows the tests must run
    inside WSL. This launcher translates the repo path to its WSL form and hands off
    to scripts/run_otilib_tests_wsl.sh, which does the deterministic activation
    (conda env + OTILIB_ROOT + RUN_OTILIB_TESTS=1).

    It either RUNS the real OTILib tests or exits non-zero explaining exactly why
    it cannot. It never lets them silently skip.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File scripts\run_otilib_tests_wsl.ps1

.NOTES
    Pass-through env overrides (optional):
      $env:OTILIB_ROOT, $env:OTILIB_CONDA_ENV, $env:CONDA_SH
#>

$ErrorActionPreference = "Stop"

# repo root = parent of this script's directory
$repo = Split-Path -Parent $PSScriptRoot

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
    Write-Host ""
    Write-Host "CANNOT RUN THE OTILIB TESTS" -ForegroundColor Red
    Write-Host "--------------------------"
    Write-Host "WSL is not available on this machine. OTILib is GPLv3 and only builds"
    Write-Host "on Linux; on Windows it must run inside WSL."
    Write-Host "Install WSL (wsl --install), then: bash scripts/setup_otilib.sh"
    exit 1
}

# translate C:\... -> /mnt/c/...
$wslRepo = (wsl.exe wslpath -a ("{0}" -f $repo.Replace('\','/'))) 2>$null
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($wslRepo)) {
    Write-Host "ERROR: could not translate '$repo' to a WSL path (wslpath failed)." -ForegroundColor Red
    exit 1
}
$wslRepo = $wslRepo.Trim()

# forward optional overrides into the WSL shell
$fwd = @()
foreach ($v in @("OTILIB_ROOT", "OTILIB_CONDA_ENV", "CONDA_SH")) {
    $val = [Environment]::GetEnvironmentVariable($v)
    if ($val) { $fwd += ("{0}='{1}'" -f $v, $val) }
}
$prefix = if ($fwd.Count) { ($fwd -join " ") + " " } else { "" }

Write-Host "== launching OTILib tests inside WSL ==" -ForegroundColor Cyan
Write-Host "windows repo : $repo"
Write-Host "wsl repo     : $wslRepo"
Write-Host ""

$cmd = "cd '$wslRepo' && ${prefix}bash scripts/run_otilib_tests_wsl.sh"
wsl.exe -e bash -lc $cmd
exit $LASTEXITCODE
