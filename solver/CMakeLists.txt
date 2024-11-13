

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
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:10485760") # 10 MB stack

if(WIN32)
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

    set(INTEL_DLLS
        "C:/Program Files (x86)/Intel/oneAPI/compiler/2024.2/bin/libifcoremd.dll"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/2024.2/bin/libmmd.dll"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/2024.2/bin/libiomp5md.dll"
    )
    
    foreach(DLL ${INTEL_DLLS})
        if(EXISTS ${DLL})
            file(COPY ${DLL} DESTINATION ${CMAKE_BINARY_DIR}/bin)
        else()
            message(WARNING "Intel DLL not found: ${DLL}")
        endif()
    endforeach()
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
