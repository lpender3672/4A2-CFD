

file(GLOB SOURCES *.cpp **.f90 gui/*.cpp block_mesh_solver/*.f90 curve_fill_solver/*.f90)

qt_add_resources(RESOURCE_FILES gui/resources.qrc)

include_directories(${CMAKE_Fortran_INCLUDE_DIRECTORIES})

if (WIN32)
    qt_add_executable(solverApp ${SOURCES} ${RESOURCE_FILES} gui/app.rc)
    set_target_properties(solverApp PROPERTIES WIN32_EXECUTABLE OFF)
else()
    qt_add_executable(solverApp ${SOURCES} ${RESOURCE_FILES})
endif()

target_link_libraries(solverApp PRIVATE Qt6::Core Qt6::Widgets Qt6::Charts Qt6::PrintSupport)

if(WIN32)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} /O2")
    # Only add the /STACK linker flag for MSVC targets. MinGW's ld doesn't accept /STACK.
    if (MSVC)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:100485760") # 10 MB stack
    else()
        message(VERBOSE "Not MSVC: skipping /STACK linker flag")
    endif()

    find_program(WINDEPLOYQT_EXECUTABLE windeployqt
        HINTS "${QT_INSTALL_DIR}/bin"
    )
    if(WINDEPLOYQT_EXECUTABLE)
        add_custom_command(TARGET solverApp POST_BUILD
            COMMAND ${WINDEPLOYQT_EXECUTABLE} --dir $<TARGET_FILE_DIR:solverApp>/bin $<TARGET_FILE:solverApp>
            COMMENT "Deploying Qt dependencies to bin directory using windeployqt"
        )
    else()
        message(WARNING "windeployqt not found! Ensure Qt DLLs are included manually.")
    endif()

endif()

install(TARGETS solverApp
    RUNTIME DESTINATION bin              # For executables
    LIBRARY DESTINATION lib              # For shared libraries (.so/.dll)
    ARCHIVE DESTINATION lib/static       # For static libraries (.a/.lib)
)

install(FILES gui/resources.qrc
    DESTINATION share/solverApp/resources # Adjust destination as needed
)

install(TARGETS solverApp
    RUNTIME DESTINATION bin
)
install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/
    DESTINATION bin
    FILES_MATCHING PATTERN "*.dll"
)
