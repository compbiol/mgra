cmake_minimum_required(VERSION 2.8.10)

project (mgra)

set(PROJECT_NAME mgra)

#set(CMAKE_BUILD_TYPE Release)

set (MGRA_VERSION_MAJOR 2)
set (MGRA_VERSION_MINOR 2)
set (MGRA_VERSION_PATCH 1)

# The 4 following lines are for CLion and do nothing (as far as I, nkartashov know) for compilation
FILE(GLOB_RECURSE ALL_HEADERS ${CMAKE_SOURCE_DIR}/*.h)
ADD_CUSTOM_TARGET(ADD_HEADERS SOURCES ${ALL_HEADERS})
FILE(GLOB_RECURSE ALL_NEW_HEADERS ${CMAKE_SOURCE_DIR}/*.hpp)
ADD_CUSTOM_TARGET(ADD_NEW_HEADERS SOURCES ${ALL_NEW_HEADERS})

set(MGRA_SOURCE mgra.cpp reader.cpp)

#	Wstats.cpp


configure_file (
  ${PROJECT_SOURCE_DIR}/include/defined.h.in
  ${PROJECT_BINARY_DIR}/defined.h
  [@ONLY]
  )

set(ROOT_DIR "${CMAKE_SOURCE_DIR}/..")
set (OPTIMIZE_FLAG " -std=c++11 -O2 -pipe -finline-functions -finline-limit=1000") # -DPSEUDO_EDGE
set (WARNS "-Wall -Werror -Wextra -pedantic-errors")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG} ${WARNS}")

include_directories(${CMAKE_BINARY_DIR}
			${PROJECT_SOURCE_DIR}/include
			${PROJECT_SOURCE_DIR}/utils
			${PROJECT_SOURCE_DIR}/libs
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

add_subdirectory(utils/io)
add_subdirectory(libs/blossom5)
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