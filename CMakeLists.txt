cmake_minimum_required(VERSION 3.16)

# use ifx not ifort
set (CMAKE_Fortran_COMPILER "ifx")

project(solver VERSION 1.2.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_PREFIX_PATH "C:/Qt/6.8.0/msvc2022_64")

find_package(Qt6 REQUIRED COMPONENTS Core Widgets Charts PrintSupport)
include(CPack)
qt_standard_project_setup()

add_subdirectory(solver)
add_subdirectory(advection)
