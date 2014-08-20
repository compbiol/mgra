Installation instructions for MGRA
==================================
Availability
------------
MGRA was tested under Linux and Mac OS.

Build requirements
------------------
* CMake 2.8.10+
* Make
* GCC C++ compiler with C++11 support (version 4.6.3+ works fine)
* BOOST library 1.54.0+

Binary distribution
-------------------
While we recommend to build MGRA from source on each machine, you also can
use pre-compiled binaries which are available for Linux and Mac OS from 
Releases page on Github: https://github.com/ablab/mgra/releases

In this case, you do not need any installation procedures.

Building from the Source Code
-----------------------------
To build and install the program under Linux, type following:

	mkdir build
	cd build
	cmake ../src
	make
	sudo make install 

If you do not have the root access or you want to install MGRA to a
location other than /usr/local/bin, run CMake with -DCMAKE_INSTALL_PREFIX
option set:

	cd build
	cmake ../src -DCMAKE_INSTALL_PREFIX="<install destination>"
	make
	make install

