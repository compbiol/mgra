# -*- cmake -*-

project(recover_tree CXX)

set(PROJECT_NAME recover_tree)
include_directories(${MGRA_MAIN_INCLUDE_DIR})
include_directories(${MGRA_MAIN_LIB_DIR})
include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(recover_tree recover_tree.cpp command_line_parsing.cpp) #${MGRA_MAIN_SRC_DIR}/include/reader/reader.cpp)
target_link_libraries(recover_tree ${MGRA_COMMON_LIBRARIES} ${Boost_LIBRARIES})


