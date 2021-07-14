cmake_minimum_required(VERSION 3.1)

project(Solver_Test)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

set(SOLVTEST_SRC main.cpp)

find_package(MPILibraries REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDES})
link_directories(${MKL_LIBRARYDIR})
link_directories(${MKL_BINARYDIR})
link_directories(${OMP_LIBRARYDIR})

add_executable(${PROJECT_NAME} ${SOLVTEST_SRC})
target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES} ${MKL_LIBRARIES_ILP})


#target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES} impi.lib mkl_core.lib mkl_intel_thread.lib mkl_intel_ilp64.lib mkl_blacs_intelmpi_ilp64.lib)
