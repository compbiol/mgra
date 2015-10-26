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
    message("Start build test library")
    set(MGRA_TESTS_PROJECT_NAME mgra_tests)
    set(MGRA_TESTS_DIRECTORY tests)

    FILE(GLOB_RECURSE MGRA_TESTS_SOURCE ${MGRA_TESTS_DIRECTORY}/*.cpp)

    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    add_subdirectory(server_test_library/gtest-1.7.0)

    enable_testing()

    include_directories(${MGRA_TESTS_DIRECTORY}/include ${MGRA_MAIN_INCLUDE_DIR} ${MGRA_MAIN_LIB_DIR})
    add_executable(${MGRA_TESTS_PROJECT_NAME} ${MGRA_TESTS_SOURCE})
    target_link_libraries(${MGRA_TESTS_PROJECT_NAME} gtest gtest_main)

    # We need thread support
    #	find_package(Threads REQUIRED)

    # Enable ExternalProject CMake module
    #    include(ExternalProject)

    # Download and install GoogleTest
    #    ExternalProject_Add(
    #        gtest
    #        URL https://googletest.googlecode.com/files/gtest-1.7.0.zip
    #        DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest
    # Disable install step
    #	UPDATE_COMMAND ""
    #	INSTALL_COMMAND ""
    #	PATCH_COMMAND ""
    #    )

    # Create a libgtest target to be used as a dependency by test programs
    #	add_library(libgtest IMPORTED STATIC GLOBAL)
    #	add_dependencies(libgtest gtest)

    # Set gtest properties
    #    ExternalProject_Get_Property(gtest source_dir binary_dir)
    #    set_target_properties(libgtest PROPERTIES
    #        "IMPORTED_LOCATION" "${binary_dir}/libgtest.a"
    #        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    #    )
    #    include_directories("${source_dir}/include")

    # Download and install GoogleMock
    #    ExternalProject_Add(
    #        gmock
    #        URL https://github.com/google/googlemock/archive/release-1.7.0.zip
    #	DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}/gmock
    # Disable install step
    #	UPDATE_COMMAND ""
    #	INSTALL_COMMAND ""
    #	PATCH_COMMAND ""
    #)

    # Create a libgmock target to be used as a dependency by test programs
    #    add_library(libgmock IMPORTED STATIC GLOBAL)
    #    add_dependencies(libgmock gmock)

    #    Set gmock properties
    #    ExternalProject_Get_Property(gmock source_dir binary_dir)
    #    set_target_properties(libgmock PROPERTIES
    #        "IMPORTED_LOCATION" "${binary_dir}/libgmock.a"
    #        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    #    )
    #    include_directories("${source_dir}/include")

    #    set(TEST_PROJECT_NAME ${PROJECT_NAME}_test)
    #    set(TESTS_DIRECTORY tests)
    #    FILE(GLOB_RECURSE TEST_SOURCE ${TESTS_DIRECTORY}/*.cpp)

    #    include_directories(${TESTS_DIRECTORY}/include)

    #    add_executable(${TEST_PROJECT_NAME} ${TEST_SOURCE})
    #    target_link_libraries(${TEST_PROJECT_NAME} libgtest)

    #    add_test(lol ${TEST_PROJECT_NAME})
endif()

