cmake_minimum_required(VERSION 2.8)

project (MGRA)

set(PROJECT mgra)

set(MGRA_SOURCE mgra.cpp
	reader.cpp
	pconf.cpp
	mpbgraph.cpp
	genome.cpp
	genome_match.cpp
	gen_dist.cpp
	Wstats.cpp
	Wdots.cpp
	mcolor.cpp
	2break.cpp
   )

set (ERR_NO_UNIX "Cannot build on non Unix systems")

set (DEBUG -DVERSION2)
set (WARNS -Wall -Wextra -pedantic-errors)
set (OPTIMIZE_FLAG " -std=c++0x -O2 -pipe -finline-functions -finline-limit=1000") #-m64 -march=core2

include_directories(${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include/writer ${PROJECT_SOURCE_DIR}/include/algo) 

if (CMAKE_COMPILER_IS_GNUCXX)	
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG}")# ${DEBUG}")
else (CMAKE_COMPILER_IS_GNUCXX)
	message (FATAL_ERROR ${ERR_NO_UNIX})
endif (CMAKE_COMPILER_IS_GNUCXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
	message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}") 

add_executable(${PROJECT} ${MGRA_SOURCE})
