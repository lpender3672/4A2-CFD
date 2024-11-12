

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

set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} /O2 /fp:fast")
