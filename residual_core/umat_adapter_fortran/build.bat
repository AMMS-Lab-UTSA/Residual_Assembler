@echo off
rem =====================================================================
rem  build.bat  -- Windows (cmd.exe) gfortran build for the standalone
rem                UMAT replay adapter.  Mirror of build.sh.
rem
rem  Produces (in .\build):
rem    umat_driver_mock.exe  -- driver + mock UMAT.  Always builds.
rem    umat_driver.exe       -- driver + REAL Grilli UMAT + stubs.  Only
rem                             links if umat.for compiles (needs ifort +
rem                             Abaqus / MKL, or a patched tree; see README).
rem
rem  Usage:  build.bat  [path-to-ngrilli_Oxford_Crystal_Plasticity]
rem =====================================================================
setlocal
set HERE=%~dp0
if "%~1"=="" (
  set "SRC=%HERE%..\..\sources\permissive\ngrilli_Oxford_Crystal_Plasticity"
) else (
  set "SRC=%~1"
)
set "STUBS=%HERE%aba_stubs"
set "BUILD=%HERE%build"
if not defined FC set FC=gfortran

if not exist "%BUILD%" mkdir "%BUILD%"
cd /d "%BUILD%"

echo == toolchain ==
%FC% --version
echo SRC   = %SRC%
echo STUBS = %STUBS%
echo.

set "FREE=-O2 -fbacktrace"
set "FIXED=-O2 -ffixed-line-length-none"
set "UMATFLAGS=-c -cpp -fcray-pointer -ffixed-line-length-none -fno-range-check -std=legacy -fallow-argument-mismatch -fallow-invalid-boz -w -I"%STUBS%" -I"%SRC%""

echo == compile driver + mock + stubs (our code) ==
%FC% %FREE%  -c "%HERE%umat_driver.f90" -o umat_driver.o   || goto :err
%FC% %FREE%  -c "%STUBS%\umat_mock.f90"  -o umat_mock.o     || goto :err
%FC% %FIXED% -c "%STUBS%\aba_stubs.f"    -o aba_stubs.o     || goto :err
%FC% %FIXED% -c "%STUBS%\lapack_stub.f"  -o lapack_stub.o   || goto :err

echo.
echo == link umat_driver_mock (plumbing test binary) ==
%FC% %FREE% umat_driver.o umat_mock.o -o umat_driver_mock.exe
if exist umat_driver_mock.exe ( echo OK  -^> build\umat_driver_mock.exe ) else ( echo FAILED to link mock. )

echo.
echo == attempt to compile the REAL Grilli UMAT (umat.for) ==
echo    log -^> build\umat_compile.log
if exist "%SRC%\umat.for" (
  %FC% %UMATFLAGS% "%SRC%\umat.for" -o umat.o 2> umat_compile.log
  if exist umat.o (
    echo umat.for COMPILED. Linking real umat_driver ...
    %FC% %FREE% umat_driver.o umat.o aba_stubs.o lapack_stub.o -o umat_driver.exe
    if exist umat_driver.exe ( echo OK  -^> build\umat_driver.exe ) else ( echo umat.o built but link failed. )
  ) else (
    echo umat.for did NOT compile with gfortran ^(expected^).
    echo See build\umat_compile.log.  Build the real UMAT with Intel ifort + Abaqus:
    echo    abaqus make library=umat.for      ^(or ifort via the Abaqus env^)
  )
) else (
  echo SKIP: umat.for not present at "%SRC%".
)

echo.
echo == done. binaries in %BUILD% ==
goto :eof

:err
echo BUILD ERROR compiling our own sources ^(unexpected^). See messages above.
exit /b 1
