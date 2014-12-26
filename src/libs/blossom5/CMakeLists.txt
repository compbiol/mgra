project (BLOSSOM5)

set(PROJECT_NAME BLOSSOM5)

set(BLOSSOM5_SOURCE 
		PMduals.cpp
		PMexpand.cpp
		PMinit.cpp
		PMinterface.cpp
		PMmain.cpp
		PMrepair.cpp
		PMshrink.cpp
		misc.cpp
		MinCost/MinCost.cpp
	 ) 

set (OPTIMIZE_FLAG " -O2 -pipe -finline-functions -finline-limit=1000")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZE_FLAG}")

include_directories(${CMAKE_BINARY_DIR} 
					${PROJECT_SOURCE_DIR}/MinCost
					)

add_library(${PROJECT_NAME} STATIC ${BLOSSOM5_SOURCE})

