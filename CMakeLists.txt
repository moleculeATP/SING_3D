cmake_minimum_required(VERSION 3.15)
project(SING CXX)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(BOOST_ROOT $ENV{CONDA_PREFIX})
set(BOOST_INCLUDEDIR $ENV{CONDA_PREFIX}/include)
set(BOOST_LIBRARYDIR $ENV{CONDA_PREFIX}/lib)

# Find packages
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph)
find_package(OpenMP REQUIRED)

set(SOURCES
    src/main.cpp
    src/io_sing.cpp
    src/distances.cpp
    src/utils.cpp
)

add_executable(sing_3D ${SOURCES})

target_include_directories(sing_3D PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    $ENV{CONDA_PREFIX}/include
)

target_link_libraries(sing_3D
    PUBLIC
    Eigen3::Eigen
    Boost::boost
    Boost::graph
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(sing_3D PUBLIC OpenMP::OpenMP_CXX)
endif()
