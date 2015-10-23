############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(io_lib CXX)

include_directories(${MGRA_MAIN_INCLUDE_DIR}
                    ${MGRA_MAIN_LIB_DIR})

add_library(io_lib STATIC path_helper.cpp copy_file.cpp)


