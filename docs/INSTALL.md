Installation instructions for MGRA
==================================
Supported operating systems
------------
MGRA was tested under Linux.

Binary distribution
-------------------
While we recommend to build MGRA from the source code, it is also possible to 
use pre-compiled binaries, which are available for Linux and Mac OS from 
at Github: https://github.com/ablab/mgra/releases

Build requirements
------------------
* CMake 2.8.10+
* Make
* GCC C++ compiler with C++11 support (version 4.8.3+ works fine)

Building from the source code
-----------------------------
Building and installing MGRA from the source code is done with the following commands: 

	mkdir build
	cd build
	cmake ../
	make
	sudo make install 

If user does not have the root access or wants to install MGRA to a
location other than /usr/local/bin, one needs to run CMake with -DCMAKE_INSTALL_PREFIX
option:

	cd build
	cmake ../ -DCMAKE_INSTALL_PREFIX="<install destination>"
	make
	make install

