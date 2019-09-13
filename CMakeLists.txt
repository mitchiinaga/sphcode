cmake_minimum_required(VERSION 3.13)
project(sphcode CXX)

find_package(OpenMP REQUIRED)

add_executable(sph)
target_compile_features(sph PUBLIC cxx_std_14)
target_compile_options(sph
  PUBLIC
    -Wall
    -Wno-sign-compare
    -Wno-maybe-uninitialized
    -funroll-loops
    -ffast-math
  )
target_include_directories(sph PUBLIC include)
target_link_libraries(sph PUBLIC OpenMP::OpenMP_CXX)

add_subdirectory(include)
add_subdirectory(src)
