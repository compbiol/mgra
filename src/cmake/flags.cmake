# -*- cmake -*-

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7")
        message(SEND_ERROR "GCC version must be at least 4.7!")
        return()
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(SEND_ERROR "Support for Clang will be in next versions!")
    return()
    #if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.2")
        #message(SEND_ERROR "Clang version must be at least 3.2!")
        #return()
    #endif()
else()
    message(SEND_ERROR "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
    return()
endif()


# So we can check build type throughout the script
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(TESTS_ON true)
endif()
