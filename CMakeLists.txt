cmake_minimum_required(VERSION 3.16)
project(pz_7)

set(CMAKE_CXX_STANDARD 17)

add_executable(pz_7
        src/Fourier.cpp
        src/Operators.cpp
        src/Wavelet_Analysis.cpp
        src/Wavelets.cpp
        headers/consts.h
)

target_include_directories(pz_7 PRIVATE ${CMAKE_SOURCE_DIR}/headers)
