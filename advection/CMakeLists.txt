

qt_add_executable(advectionApp
    advection.cpp
    ../solver/block_mesh_solver/routines.f90
    advection.f90
    example.f90
)

target_link_libraries(advectionApp PRIVATE Qt6::Core Qt6::Widgets)

include_directories(${CMAKE_Fortran_INCLUDE_DIRECTORIES})

if (WIN32)
    set_target_properties(advectionApp PROPERTIES WIN32_EXECUTABLE OFF)
endif()
