# -*- cmake -*-

project(recover_tree CXX)

set(PROJECT_NAME recover_tree)
include_directories(${MGRA_SOURCE_DIR}/include)
include_directories(${MGRA_SOURCE_DIR}/libs)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(recover_tree ${COMMON_LIBRARIES})

add_executable(recover_tree main.cpp
                command_line_parsing.cpp
                ${MGRA_SOURCE_DIR}/include/reader/reader.cpp
                )

