
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(Final_Project_GPU LANGUAGES CXX CUDA)

add_library(drivers STATIC
    canopy_air_space.cu
    canopy_air_space.h
    numerics.cu
    numerics.h
    photosynthesis.cu
    photosynthesis.h
    radiation_solver.cu
    radiation_solver.h
    parameters.h
)

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(drivers PUBLIC cxx_std_11)

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties( drivers
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON
                     )

add_executable(gpu cohort_serial.cu)
set_target_properties(gpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(gpu PRIVATE drivers)

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60 --expt-extended-lambda")
