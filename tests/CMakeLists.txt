enable_testing()

# Fetch doctest
include(FetchContent)
FetchContent_Declare(
    doctest
    GIT_REPOSITORY "https://github.com/onqtam/doctest"
    GIT_TAG "ae7a13539fb71f270b87eb2e874fbac80bc8dda2"
)
FetchContent_MakeAvailable(doctest)

# Create test executable with a unique name
add_executable(colorama_tests 
    ${PROJECT_SOURCE_DIR}/tests/tests.cpp
    color.cpp
)

# Include doctest cmake utilities
include(${doctest_SOURCE_DIR}/scripts/cmake/doctest.cmake)

# Link against both doctest and colorama
target_link_libraries(colorama_tests 
    PRIVATE 
        doctest
        colorama
)

# Make sure tests can find colorama headers
target_include_directories(colorama_tests
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

# Add post-build commands to copy DLLs next to the test executable
if(WIN32)
    add_custom_command(TARGET colorama_tests POST_BUILD
        # Copy colorama DLL
        COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:colorama_tests>"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:colorama>
            $<TARGET_FILE:spacely>
            $<TARGET_FILE_DIR:colorama_tests>
    )
endif()

# For discovery of tests
set_property(TARGET colorama_tests PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:colorama_tests>")

# Discover and add tests
doctest_discover_tests(colorama_tests
    WORKING_DIRECTORY "$<TARGET_FILE_DIR:colorama_tests>"
)

# Register test with correct working directory
add_test(NAME colorama_tests 
    COMMAND colorama_tests
    WORKING_DIRECTORY "$<TARGET_FILE_DIR:colorama_tests>"
)