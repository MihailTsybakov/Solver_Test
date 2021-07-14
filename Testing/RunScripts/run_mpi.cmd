@ECHO off

CALL config.cmd

SET NumProcs=2
SET MPIExecExe=%I_MPI_ROOT%\intel64\bin\mpiexec.exe
SET MPIPPrefix="%MPIExecExe%" -n %NumProcs% -wdir %MyWorkDir% -envall -print-all-exitcodes 

pushd %MyWorkDir%

%MPIPPrefix% %MyProjectDir%\%TypeBuild%\%MyExeName%.exe

popd