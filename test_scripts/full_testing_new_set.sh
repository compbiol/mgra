#!/bin/bash 

cd ..

HOME_DIR=$PWD

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

for x in $HOME_DIR/once_test/*
do
	cp -rf $x $HOME_DIR/current/
done

for x in $HOME_DIR/current/*
do 
	cp $MGRA_TOOL $x/
done

echo "RUN MGRA AND COMPARE"
for x in $HOME_DIR/current/*
do 
	cd $x 
	NAME_DIR=${PWD##*/}

	echo "WORK in $NAME_DIR directory"

	NAME_CFG=$(find -name *.cfg)
	./mgra $NAME_CFG >/dev/null 2>logerror.txt
done 
