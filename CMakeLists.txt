cmake_minimum_required(VERSION 3.1)

project(Solver_Test)

include_directories(${PROJECT_SOURCE_DIR})

set(SOLVTEST_SRC main.cpp)

add_executable(${PROJECT_NAME} ${SOLVTEST_SRC})