cmake_minimum_required(VERSION 3.15)
project(SING CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Boost and Eigen
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph)

# Collect sources
set(SOURCES
    src/main.cpp
    src/io.cpp
    src/distances.cpp
    src/utils.cpp
)

# Add executable
add_executable(sing_3D ${SOURCES})

# Add include directories (for io.h and for gudhi headers)
target_include_directories(sing_3D PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    $ENV{CONDA_PREFIX}/include
)


# Link libraries
target_link_libraries(sing_3D Eigen3::Eigen Boost::boost Boost::graph)
