# In this CMake file, we will find all required packages


# Find the Boost package - needed for unittests
find_package(Boost REQUIRED)

# Find Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Find bmqc for bitset manipulations
find_package(bmqc 0.1.0 REQUIRED)
