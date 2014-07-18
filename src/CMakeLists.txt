cmake_minimum_required(VERSION 2.8.10)

project (mgra)

set(PROJECT_NAME mgra)

#set(CMAKE_BUILD_TYPE Release)

set (MGRA_VERSION_MAJOR 2)
set (MGRA_VERSION_MINOR 2)
set (MGRA_VERSION_PATCH 0)

set(MGRA_SOURCE mgra.cpp
	reader.cpp
	Wstats.cpp
   )

configure_file (
  ${PROJECT_SOURCE_DIR}/include/defined.h.in
  ${PROJECT_BINARY_DIR}/defined.h
  [@ONLY]
  )

set(ROOT_DIR "${CMAKE_SOURCE_DIR}/../")
set (OPTIMIZE_FLAG " -std=c++11 -O2 -pipe -finline-functions -finline-limit=1000")
set (WARNS "-Wall -Werror -Wextra -pedantic-errors")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG}") #" ${WARNS}") 

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG} ${DEBUG}")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG} ${DEBUG}")
# -DBOOST_LOG_DYN_LINK")


include_directories(${CMAKE_BINARY_DIR} 
			${PROJECT_SOURCE_DIR}/include 
			${PROJECT_SOURCE_DIR}/include/writer 
			${PROJECT_SOURCE_DIR}/include/algo 
			${PROJECT_SOURCE_DIR}/include/graph
			${PROJECT_SOURCE_DIR}/include/structures
			${PROJECT_SOURCE_DIR}/include/event
			${PROJECT_SOURCE_DIR}/libs 
		   ) 

find_package(Boost COMPONENTS program_options filesystem system REQUIRED) #log log_setup thread date_time 
if(NOT ${Boost_FOUND})
    message(SEND_ERROR "Failed to find boost::program_options or boost::filesystem or boost::system.")
    return()
else()
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIR})
endif()

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

add_subdirectory(libs/blossom5)
set_property(TARGET BLOSSOM5 APPEND_STRING PROPERTY COMPILE_FLAGS " -w")

add_executable(${PROJECT_NAME} ${MGRA_SOURCE})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
add_dependencies(${PROJECT_NAME} BLOSSOM5)
target_link_libraries(${PROJECT_NAME} BLOSSOM5)

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
set(CPACK_RESOURCE_FILE_LICENSE "${ROOT_DIR}/COPYING")
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

