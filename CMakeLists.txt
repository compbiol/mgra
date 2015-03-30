cmake_minimum_required(VERSION 2.8.10)

project (mgra)

set(PROJECT_NAME mgra)

set (MGRA_VERSION_MAJOR 2)
set (MGRA_VERSION_MINOR 2)
set (MGRA_VERSION_PATCH 1)

# So we can check build type throughout the script
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Google Test

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    # We need thread support
    find_package(Threads REQUIRED)

    # Enable ExternalProject CMake module
    include(ExternalProject)

    # Download and install GoogleTest
    ExternalProject_Add(
        gtest
        URL https://googletest.googlecode.com/files/gtest-1.7.0.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
        # Disable install step
        INSTALL_COMMAND ""
    )

    # Create a libgtest target to be used as a dependency by test programs
    add_library(libgtest IMPORTED STATIC GLOBAL)
    add_dependencies(libgtest gtest)

    # Set gtest properties
    ExternalProject_Get_Property(gtest source_dir binary_dir)
    set_target_properties(libgtest PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    #    "INTERFACE_INCLUDE_DIRECTORIES" "${source_dir}/include"
    )
    # I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
    include_directories("${source_dir}/include")

    # Download and install GoogleMock
    ExternalProject_Add(
        gmock
        URL https://googlemock.googlecode.com/files/gmock-1.7.0.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gmock
        # Disable install step
        INSTALL_COMMAND ""
    )

    # Create a libgmock target to be used as a dependency by test programs
    add_library(libgmock IMPORTED STATIC GLOBAL)
    add_dependencies(libgmock gmock)

    # Set gmock properties
    ExternalProject_Get_Property(gmock source_dir binary_dir)
    set_target_properties(libgmock PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/libgmock.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
    #    "INTERFACE_INCLUDE_DIRECTORIES" "${source_dir}/include"
    )
    # I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
    include_directories("${source_dir}/include")
endif()

# The 4 following lines are for CLion and do nothing (as far as I, nkartashov know) for compilation
FILE(GLOB_RECURSE ALL_HEADERS ${CMAKE_SOURCE_DIR}/*.h)
ADD_CUSTOM_TARGET(ADD_HEADERS SOURCES ${ALL_HEADERS})
FILE(GLOB_RECURSE ALL_NEW_HEADERS ${CMAKE_SOURCE_DIR}/*.hpp)
ADD_CUSTOM_TARGET(ADD_NEW_HEADERS SOURCES ${ALL_NEW_HEADERS})

set(MGRA_SOURCE src/mgra.cpp src/reader.cpp)

#	Wstats.cpp

configure_file (
  src/include/defined.h.in
  ${PROJECT_BINARY_DIR}/defined.h
  [@ONLY]
  )

set(ROOT_DIR "${CMAKE_SOURCE_DIR}/..")
set(COMPILE_OPTION_FLAG "-std=c++11" )
set(DEBUG_FLAG "-g3")
set(OPTIMIZE_FLAG "-O2 -pipe -finline-functions -finline-limit=1000") # -DPSEUDO_EDGE
set(WARNS "-Wall -Werror -Wextra -pedantic-errors")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILE_OPTION_FLAG} ${OPTIMIZE_FLAG} ${WARNS}")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILE_OPTION_FLAG} ${DEBUG_FLAG} ${WARNS}")
endif()

include_directories(${CMAKE_BINARY_DIR}
			src/include
			src/utils
			src/libs
		  )

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7")
        message(SEND_ERROR "GCC version must be at least 4.7!")
        return()
    endif()
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.2")
#        message(SEND_ERROR "Clang version must be at least 3.2!")
#        return()
#    endif()
else()
    message(SEND_ERROR "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
    return()
endif()

add_subdirectory(src/utils/io)
add_subdirectory(src/libs/blossom5)
set_property(TARGET BLOSSOM5 APPEND_STRING PROPERTY COMPILE_FLAGS " -w")

add_executable(${PROJECT_NAME} ${MGRA_SOURCE})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

add_dependencies(${PROJECT_NAME} BLOSSOM5)
target_link_libraries(${PROJECT_NAME} BLOSSOM5)

add_dependencies(${PROJECT_NAME} input)
target_link_libraries(${PROJECT_NAME} input)


if(WIN32)
	message(SEND_ERROR "Windows is unsupported system, please use UNIX or OS X systems")
else()
	set(LIB_DIR "lib/MGRA")
	set(SHARE_DIR "share/MGRA")
	set(DOC_DIR "share/MGRA/doc")
	install(TARGETS mgra RUNTIME DESTINATION bin)
	install(FILES ${ROOT_DIR}/NEWS.md ${ROOT_DIR}/README.md ${ROOT_DIR}/USAGE.md ${ROOT_DIR}/INSTALL.md ${ROOT_DIR}/COPYING DESTINATION ${DOC_DIR})
endif()

set(CPACK_NSIS_MODIFY_PATH "ON")
set(CPACK_SOURCE_IGNORE_FILES "\\\\.git.*;/build/.+;${CPACK_SOURCE_IGNORE_FILES};")
#set(CPACK_RESOURCE_FILE_LICENSE "${ROOT_DIR}/COPYING")
set(CPACK_PACKAGE_NAME "MGRA")
set(CPACK_PACKAGE_VENDOR "George Washington University")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MGRA is a tool used for reconstruction genomes and transformation history")
set(CPACK_PACKAGE_VERSION_MAJOR ${MGRA_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${MGRA_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${MGRA_VERSION_PATCH})
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "MGRA")
set(CMAKE_SOURCE_DIR ${ROOT_DIR})
include(CPack)