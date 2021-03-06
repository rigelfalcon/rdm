# src-level CMakeLists.txt
cmake_minimum_required(VERSION 3.10)


# Specify the project, its version and the used language
project(rdm
        VERSION 0.0.1
        LANGUAGES CXX)

# To enable the C++11 standard, we have to add flags like -std=gnu+11 or so. CMake can take care of this:
set(CMAKE_CXX_STANDARD 11)

# Set some compiler optimizations (https://stackoverflow.com/a/41361741/7930415)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Set CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/SetCMakeEnvironment.cmake)

# Find all required packages
include(${CMAKE_SOURCE_DIR}/cmake/FindPackages.cmake)

# Add the sources for the library
add_subdirectory(${PROJECT_SOURCE_FOLDER})



# Set up testing
enable_testing()
add_subdirectory(${PROJECT_TESTS_FOLDER})
