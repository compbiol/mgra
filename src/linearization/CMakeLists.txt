# -*- cmake -*-

project(linearization CXX)

set(PROJECT_NAME linearization)
include_directories(${MGRA_MAIN_INCLUDE_DIR})
include_directories(${MGRA_MAIN_LIB_DIR})
include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(linearization linearization.cpp) # command_line_parsing.cpp ${MGRA_MAIN_SRC_DIR}/include/reader/reader.cpp)
target_link_libraries(linearization ${MGRA_COMMON_LIBRARIES} ${Boost_LIBRARIES})


