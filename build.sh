#!/bin/bash

HOME_DIR=$PWD

#############################
## MAIN FUNCTION ############
#############################

# Make Build Directory
mkdir build
BUILD_DIR=${HOME_DIR}/build

# Copy CMakeLists.txt and build MGRA
cp ${HOME_DIR}/CMakeLists.txt ${BUILD_DIR}
cd ${BUILD_DIR}
cmake -DCMAKE_BUILD_TYPE=Debug ${HOME_DIR}
make

# Run tests
# ./MGRA_test

#############################
## END ######################
#############################

