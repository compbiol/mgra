#!/bin/bash 

cd ..

HOME_DIR=$PWD
TEST_DIRECTORY=${1}

if [ ! -d src ]; then
	echo ERROR
	exit 
fi

if [ ! -d build ]; then
   mkdir build
fi

if [ ! -d current ]; then
   mkdir current
fi

echo "MAKE PROJECT"
cd build 
rm -rf *
cmake ../src/ >/dev/null
make >make.log
MGRA_TOOL=$PWD/mgra
cd ..

echo "COPY FILES"
cd current
rm -rf *

cd ..

cp -rf $HOME_DIR/once_test/$TEST_DIRECTORY $HOME_DIR/current/
cp $MGRA_TOOL $HOME_DIR/current/$TEST_DIRECTORY

echo "RUN MGRA"
cd $HOME_DIR/current/$TEST_DIRECTORY

NAME_CFG=$(find -name *.cfg)
./mgra $NAME_CFG 


