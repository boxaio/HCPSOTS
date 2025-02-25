project(SOTS)
cmake_minimum_required(VERSION 3.12)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

set(CMAKE_BUILD_TYPE "Debug")


# add eigen
include_directories(/usr/local/include/eigen3)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP)

include(polyscope)
include(geometry-central)
include(eigen)
include(openmp)

add_subdirectory(deps/bvh)

file(GLOB SRC
     "src/**"
     "deps/CLI11.hpp"
)

message(STATUS ${SRC})


# add_executable(spherical_bluenoise experiments/spherical_bluenoise.cpp ${SRC})
# target_link_libraries(spherical_bluenoise polyscope geometry-central  bvh OpenMP::OpenMP_CXX)

# add_executable(spherical_mesh_sampling experiments/spherical_mesh_sampling.cpp ${SRC})
# target_link_libraries(spherical_mesh_sampling polyscope geometry-central  bvh OpenMP::OpenMP_CXX)


# add_executable(vizu_pc experiments/vizu_pc.cpp)
# target_link_libraries(vizu_pc polyscope geometry-central)

add_executable(tests 
               tests/HCPSOTS_bluenoise.cpp 
               # tests/HCPSOTS_mesh.cpp 
               # tests/test_NESOTS_bluenoise.cpp 
               # tests/test_NESOTS_mesh.cpp 
               ${SRC})
target_link_libraries(tests polyscope geometry-central bvh OpenMP::OpenMP_CXX)