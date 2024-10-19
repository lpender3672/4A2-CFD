

qt_add_executable(solverApp
    main.cpp
    Code/types.f90
    Code/stencils.f90
    Code/routines.f90
    Code/qtio.f90

    Code/solver.f90
    Code/read_settings.f90
    Code/allocate_arrays.f90
    Code/read_geom.f90
    Code/read_mesh.f90
    Code/generate_mesh.f90
    Code/calc_areas.f90
    Code/check_mesh.f90
    Code/flow_guess.f90
    Code/set_timestep.f90
    Code/set_secondary.f90
    Code/apply_bconds.f90
    Code/euler_iteration.f90
    Code/check_conv.f90
    Code/check_stop.f90
    Code/write_output.f90
)

target_link_libraries(solverApp PRIVATE Qt6::Core Qt6::Widgets)

include_directories(${CMAKE_Fortran_INCLUDE_DIRECTORIES})

if (WIN32)
    set_target_properties(solverApp PROPERTIES WIN32_EXECUTABLE OFF)
endif()
