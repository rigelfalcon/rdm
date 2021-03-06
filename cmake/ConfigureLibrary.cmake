# In this CMake file, we will include the headers and link to the necessary libraries


# Include this project's headers
target_include_directories(${LIBRARY_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})

# Include the boost headers (dynamic bitset)
target_include_directories(${LIBRARY_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# Include Eigen
target_link_libraries(${LIBRARY_NAME} PUBLIC Eigen3::Eigen)

# Include bmqc
target_include_directories(${LIBRARY_NAME} PUBLIC ${bmqc_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC bmqc)

# Use Eigen with MKL
# target_include_directories(${LIBRARY_NAME} PRIVATE $ENV{MKLROOT}/include)
# target_link_libraries(${LIBRARY_NAME} PRIVATE $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_sequential.a $ENV{MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)
