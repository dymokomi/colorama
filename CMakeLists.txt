cmake_minimum_required(VERSION 3.12)
project(colorama VERSION 0.1.0 LANGUAGES CXX DESCRIPTION "Colorama is a library for color manipulation in C++.")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include(GNUInstallDirs)

# Generate export header
include(GenerateExportHeader)

# Add colorama subdirectory
add_subdirectory(external/spacely)

set(SPDI_SOURCES
    src/color.cpp
    src/colorspace.cpp
    src/observer.cpp
    src/spectrum.cpp
    src/whitepoint.cpp
)
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()
# Create library target
add_library(colorama SHARED ${SPDI_SOURCES})

generate_export_header(colorama
    BASE_NAME COLORAMA
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/colorama_export.h
)

set_target_properties(colorama PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER include/color.h
    WINDOWS_EXPORT_ALL_SYMBOLS ON
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/bin/Debug"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/bin/Release"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/lib/Debug"
    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/lib/Release"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/lib/Debug"
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/lib/Release"
)

configure_file(colorama.pc.in colorama.pc @ONLY)

target_include_directories(colorama 
    PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        ${CMAKE_CURRENT_SOURCE_DIR}/external/colorama/include
)

target_link_libraries(colorama 
    PUBLIC 
        spacely
)

add_dependencies(colorama spacely)

option(BUILD_TESTS "Build tests" ON)

if(BUILD_TESTS)
  include(CTest)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tests)
endif()

install(TARGETS colorama
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${CMAKE_BINARY_DIR}/colorama.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)