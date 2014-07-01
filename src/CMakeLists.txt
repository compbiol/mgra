cmake_minimum_required(VERSION 2.8.6)

project (MGRA)

set(PROJECT mgra)

set(MGRA_SOURCE mgra.cpp
	reader.cpp
	Wstats.cpp
   )

set (ERR_NO_UNIX "Cannot build on non Unix systems")
set (DEBUG -DVERSION2)
set (OPTIMIZE_FLAG " -std=c++0x -O2 -pipe -finline-functions -finline-limit=1000")
set (WARNS "-Wall -Werror -Wextra -pedantic-errors")


find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

include_directories(${CMAKE_BINARY_DIR} 
			${PROJECT_SOURCE_DIR}/include 
			${PROJECT_SOURCE_DIR}/include/writer 
			${PROJECT_SOURCE_DIR}/include/algo 
			${PROJECT_SOURCE_DIR}/include/graph
			${PROJECT_SOURCE_DIR}/include/structures
			${PROJECT_SOURCE_DIR}/include/event
		   ) 


if (CMAKE_COMPILER_IS_GNUCXX)	
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG} ${DEBUG}") #${WARNS} 
else (CMAKE_COMPILER_IS_GNUCXX)
	message (FATAL_ERROR ${ERR_NO_UNIX})
endif (CMAKE_COMPILER_IS_GNUCXX)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
	message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}") 

add_executable(${PROJECT} ${MGRA_SOURCE})
target_link_libraries(${PROJECT} ${Boost_LIBRARIES})