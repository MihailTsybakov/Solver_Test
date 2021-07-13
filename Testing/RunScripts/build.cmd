@ECHO off

CALL config.cmd

IF EXIST "%MyProjectDir%" RD /s /q "%MyProjectDir%"
MKDIR "%MyProjectDir%"
PUSHD "%MyProjectDir%"

%CMAKE_EXE% -G%CMAKE_GENERATOR_NAME% %MyRepository% ^
-DCMAKE_C_COMPILER=cl ^
-DCMAKE_CXX_COMPILER=cl 


POPD