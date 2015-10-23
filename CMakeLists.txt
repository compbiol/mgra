############################################################################
# Copyright (c) 2009-2015 by Pavel Avdeyev, Nikita Kartashov, Max Alekseyev <maxal@gwu.edu>
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# Hard prerequisites
cmake_minimum_required(VERSION 2.8.12)

# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# The 4 following lines are for CLion and do nothing for compilation
FILE(GLOB_RECURSE ALL_HEADERS ${CMAKE_SOURCE_DIR}/*.h)
ADD_CUSTOM_TARGET(ADD_HEADERS SOURCES ${ALL_HEADERS})
FILE(GLOB_RECURSE ALL_NEW_HEADERS ${CMAKE_SOURCE_DIR}/*.hpp)
ADD_CUSTOM_TARGET(ADD_NEW_HEADERS SOURCES ${ALL_NEW_HEADERS})

# modules to find additional libraries (like zlib, boost and etc)
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
set(MGRA_DEFAULT_BUILD_TYPE "RelWithAsserts" CACHE STRING "SPAdes default build type")
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

#add_subdirectory(${MGRA_MAIN_SRC_DIR}/linearization)


#include_directories(${CMAKE_BINARY_DIR}
#			src/libs/JsonCpp/
#			src/include
#			src/recover_tree
#			src/mgra
#			src/linearization
#			src/libs
#			)



#if(WIN32)
#	message(SEND_ERROR "Windows is unsupported system, please use UNIX or OS X systems")
#else()
#	set(LIB_DIR "lib/MGRA")
#	set(SHARE_DIR "share/MGRA")
#	set(DOC_DIR "share/MGRA/doc")
#	install(TARGETS mgra RUNTIME DESTINATION bin)
#	install(FILES ${ROOT_DIR}/NEWS.md ${ROOT_DIR}/README.md ${ROOT_DIR}/USAGE.md ${ROOT_DIR}/INSTALL.md ${ROOT_DIR}/COPYING DESTINATION ${DOC_DIR})
#endif()


