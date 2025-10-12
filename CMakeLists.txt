cmake_minimum_required(VERSION 3.16)

project(solver VERSION 1.2.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS LinguistTools Core Widgets Charts PrintSupport)
include(CPack)
qt_standard_project_setup()

add_subdirectory(solver)
add_subdirectory(advection)
