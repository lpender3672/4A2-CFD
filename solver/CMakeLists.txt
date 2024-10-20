

file(GLOB SOURCES *.cpp gui/*.cpp fortran/*.f90)

qt_add_executable(solverApp ${SOURCES})

target_link_libraries(solverApp PRIVATE Qt6::Core Qt6::Widgets)

include_directories(${CMAKE_Fortran_INCLUDE_DIRECTORIES})

if (WIN32)
    set_target_properties(solverApp PROPERTIES WIN32_EXECUTABLE OFF)
endif()

include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-fast" HAS_FAST_MATH_FLAG)

if (HAS_FAST_MATH_FLAG)
    # If the flag is supported, apply it to the target
    target_compile_options(solverApp PRIVATE "-fast")
endif()
