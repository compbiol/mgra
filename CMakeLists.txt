############################################################################
# Copyright (c) 2009-2015 by Pavel Avdeyev, Nikita Kartashov, Max Alekseyev <maxal@gwu.edu>
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# Hard prerequisites
cmake_minimum_required(VERSION 2.8.12)

# Global name of project
project(MGRA CXX)

# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# The 4 following lines are for CLion and do nothing for compilation
FILE(GLOB_RECURSE ALL_HEADERS ${CMAKE_SOURCE_DIR}/*.h)
ADD_CUSTOM_TARGET(ADD_HEADERS SOURCES ${ALL_HEADERS})
FILE(GLOB_RECURSE ALL_NEW_HEADERS ${CMAKE_SOURCE_DIR}/*.hpp)
ADD_CUSTOM_TARGET(ADD_NEW_HEADERS SOURCES ${ALL_NEW_HEADERS})

# Modules to find additional libraries and flags (like blossom V, boost and etc)
# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/src/cmake")

# Define various dirs
set(MGRA_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(MGRA_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(MGRA_MAIN_LIB_DIR ${MGRA_MAIN_SRC_DIR}/libs)
set(MGRA_MAIN_INCLUDE_DIR ${MGRA_MAIN_SRC_DIR}/include)

# Define some parameters
set(MGRA_VERSION_MAJOR 2)
set(MGRA_VERSION_MINOR 4)
set(MGRA_VERSION_PATCH 0)

# Sanity check our source directory to make sure that we are not trying to
# generate an in-tree build
if ((CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR) OR (CMAKE_SOURCE_DIR STREQUAL "${CMAKE_BINARY_DIR}/src"))
	message(FATAL_ERROR "In-source builds are not allowed.
			Please create a directory and run cmake from there, passing the path
			to this source directory as the last argument. This process created the file `CMakeCache.txt'
			and the directory `CMakeFiles'. Please delete them.")
endif()

# Default configuration
set(MGRA_DEFAULT_BUILD_TYPE "RelWithAsserts" CACHE STRING "MGRA default build type")
if (NOT CMAKE_BUILD_TYPE)
	message("Setting default build configuration: ${MGRA_DEFAULT_BUILD_TYPE}")
	set(CMAKE_BUILD_TYPE "${MGRA_DEFAULT_BUILD_TYPE}" CACHE STRING
			"Choose the type of build, options are: None Debug Release RelWithAsserts RelWithDebInfo."
			FORCE)
endif()

# Add flags if necessary
include(flags)

# Add common libs
include(libs)

configure_file (${MGRA_MAIN_SRC_DIR}/include/version.hpp.in ${PROJECT_BINARY_DIR}/version.hpp)

# Build recover tree module
add_subdirectory(${MGRA_MAIN_SRC_DIR}/recover_tree)

# Build main module
add_subdirectory(${MGRA_MAIN_SRC_DIR}/mgra)

# Build linearization module
#add_subdirectory(${MGRA_MAIN_SRC_DIR}/linearization)

# Add pack option
include(pack)


