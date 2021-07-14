SET PATH=C:\Windows\system32;C:\Windows

SET PARCH=64
SET IntelCS=ICS2017UP2
SET versionVS=2015
SET MyRepository=%cd%\Repository
SET MyProjectDir=%cd%\build
SET MyWorkDir=%cd%\workdir
SET MyExeName=Solver_Test

SET TypeBuild=RelWithDebInfo
SET TypeBuild=Release
SET TypeBuild=Debug

rem save system MPI root, ipsxe-comp-vars.bat changes to internal ICS mpi root
SET I_MPI_ROOT_SYSTEM=%I_MPI_ROOT%
SET I_MPI_ONEAPI_ROOT_SYSTEM=%I_MPI_ONEAPI_ROOT%

IF %IntelCS%==ICS2017UP2 (
 set ICPP_COMPILER13=
 set IFORT_COMPILER13=
 set ICPP_COMPILER14=
 set IFORT_COMPILER14=
 set ICPP_COMPILER15=
 set IFORT_COMPILER15=
 set ICPP_COMPILER16=
 set IFORT_COMPILER16=

 call "%ICPP_COMPILER17%bin\ipsxe-comp-vars.bat" intel64 vs%versionVS%
 SET IntelConvertTool="C:\Program Files (x86)\Common Files\Intel\shared files\ia32\Bin\ICProjConvert170.exe"
)

rem restore system MPI root
SET I_MPI_ROOT=%I_MPI_ROOT_SYSTEM%
SET I_MPI_ONEAPI_ROOT=%I_MPI_ONEAPI_ROOT_SYSTEM%

IF %versionVS%==2012 (
  SET CMAKE_GENERATOR_NAME="Visual Studio 11 2012 Win64"
)
IF %versionVS%==2013 (
  SET CMAKE_GENERATOR_NAME="Visual Studio 12 2013 Win64"
)
IF %versionVS%==2015 (
  SET CMAKE_GENERATOR_NAME="Visual Studio 14 2015 Win64"
)
IF %versionVS%==2019 (
  SET CMAKE_GENERATOR_NAME="Visual Studio 16 2019 Win64"
)

set CMAKE_EXE="C:\Program Files\CMake\bin\cmake.exe"
set CMAKE_EXE_GUI="C:\Program Files\CMake\bin\cmake-gui.exe"