cmake_minimum_required ( VERSION 2.8.12 )
project ( shock_tube )

set ( CMAKE_CXX_STANDARD 17 )
set ( CMAKE_CXX_STANDARD_REQUIRED True )
set( CMAKE_BUILD_TYPE Debug)

# list(PREPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/external/yaml-cpp/build")
# find_package ( yaml-cpp REQUIRED )
# find_package(yaml-cpp REQUIRED
#     PATHS "${CMAKE_SOURCE_DIR}/external/yaml-cpp/install"
#     NO_DEFAULT_PATH
# )
find_package(OpenMP REQUIRED)
add_subdirectory(external/yaml-cpp)
include_directories (
    include
)

add_compile_options ( -std=c++17 )

add_executable(1d_shock_tube src/euler_1d.cpp src/param_loader.cpp) 
target_link_libraries(1d_shock_tube yaml-cpp)

add_executable(2d_shock_tube src/euler_2d.cpp src/param_loader.cpp) 
target_link_libraries(2d_shock_tube yaml-cpp)

add_executable(parallel_2d_shock_tube src/euler_2d_parallel.cpp src/param_loader.cpp) 
target_link_libraries(parallel_2d_shock_tube yaml-cpp OpenMP::OpenMP_CXX)

