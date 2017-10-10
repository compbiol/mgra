# -*- cmake -*-

# Add basic selfmade io library
add_subdirectory(${MGRA_MAIN_LIB_DIR}/io)

# Add jsoncpp library
add_subdirectory(${MGRA_MAIN_LIB_DIR}/JsonCpp)

set(MGRA_COMMON_LIBRARIES json io_lib)

# Compile blossom V by Vladimir Kolmogorov
add_subdirectory(${MGRA_MAIN_LIB_DIR}/blossom5)
set_property(TARGET BLOSSOM5 APPEND_STRING PROPERTY COMPILE_FLAGS " -w")

# Google Test & Mock
if (MGRA_TESTS_ON)
    #message("Start build test library")
    #set(MGRA_TESTS_PROJECT_NAME mgra_tests)
    #set(MGRA_TESTS_DIRECTORY tests)

    #FILE(GLOB_RECURSE MGRA_TESTS_SOURCE ${MGRA_TESTS_DIRECTORY}/*.cpp)

    enable_testing()

    # We need thread support
    find_package(Threads REQUIRED)

    # Enable ExternalProject CMake module
    include(ExternalProject)

    # Download and install GoogleTest
    ExternalProject_Add(
        gtest
        URL https://github.com/google/googletest/archive/master.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
        # Disable install step
        INSTALL_COMMAND ""
    )

    # Get GTest source and binary directories from CMake project
    ExternalProject_Get_Property(gtest source_dir binary_dir)

    # Create a libgtest target to be used as a dependency by test programs
    add_library(libgtest IMPORTED STATIC GLOBAL)
    add_dependencies(libgtest gtest)

    # Set libgtest properties
    set_target_properties(libgtest PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/googlemock/gtest/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    )

    # Create a libgmock target to be used as a dependency by test programs
    add_library(libgmock IMPORTED STATIC GLOBAL)
    add_dependencies(libgmock gtest)

    # Set libgmock properties
    set_target_properties(libgmock PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/googlemock/libgmock.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    )

    include_directories("${source_dir}/googletest/include"
                        "${source_dir}/googlemock/include")

    set(TEST_PROJECT_NAME ${PROJECT_NAME}_test)
    set(TESTS_DIRECTORY tests)
    FILE(GLOB_RECURSE TEST_SOURCE ${TESTS_DIRECTORY}/*.cpp)

    include_directories(${TESTS_DIRECTORY}/include ${MGRA_MAIN_INCLUDE_DIR} ${MGRA_MAIN_LIB_DIR}
    )

    add_executable(${TEST_PROJECT_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_PROJECT_NAME}
            libgtest
            libgmock
        )

    add_test(lol ${TEST_PROJECT_NAME})

    
endif()

