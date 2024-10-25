

file(GLOB SOURCES *.cpp gui/*.cpp fortran/*.f90)

qt_add_resources(RESOURCE_FILES gui/resources.qrc)



include_directories(${CMAKE_Fortran_INCLUDE_DIRECTORIES})

if (WIN32)
    qt_add_executable(solverApp ${SOURCES} ${RESOURCE_FILES} gui/app.rc)
    set_target_properties(solverApp PROPERTIES WIN32_EXECUTABLE OFF)
else()
    qt_add_executable(solverApp ${SOURCES} ${RESOURCE_FILES})
endif()

target_link_libraries(solverApp PRIVATE Qt6::Core Qt6::Widgets Qt6::Charts Qt6::PrintSupport)

include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-fast" HAS_FAST_MATH_FLAG)

if (HAS_FAST_MATH_FLAG)
    # If the flag is supported, apply it to the target
    target_compile_options(solverApp PRIVATE "-fast")
endif()
