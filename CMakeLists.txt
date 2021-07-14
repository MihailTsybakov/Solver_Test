cmake_minimum_required(VERSION 3.1)

project(Solver_Test)

set(SOLVTEST_SRC main.cpp)

include_directories(<path to mpi>\include)
include_directories(<path to mkl>\include)

link_directories(<path to mpi>\lib\debug)
link_directories(<path to mkl>\lib\intel64)

target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES} impi.lib mkl_core.lib mkl_intel_thread.lib mkl_intel_ilp64.lib mkl_blacs_intelmpi_ilp64.lib)

add_executable(${PROJECT_NAME} ${SOLVTEST_SRC})
