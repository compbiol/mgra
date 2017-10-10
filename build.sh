#!/bin/bash

HOME_DIR=$PWD

#############################
## MAIN FUNCTION ############
#############################

if [ "$#" -ne 1 ]; then
	echo "./build.sh {option}
Where option is a build type for MGRA.
Pick one: None, Debug, Release, RelWithAsserts, RelWithDebInfo
	"
	exit
fi
BUILD_TYPE=$1

# Make Build Directory
rm -rf build
mkdir build
BUILD_DIR=${HOME_DIR}/build

# Copy CMakeLists.txt and build MGRA
cp ${HOME_DIR}/CMakeLists.txt ${BUILD_DIR}
cd ${BUILD_DIR}
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${HOME_DIR}
make

#############################
## END ######################
#############################

