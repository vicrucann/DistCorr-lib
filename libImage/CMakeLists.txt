cmake_minimum_required(VERSION 2.6)
include_directories(${CMAKE_SOURCE_DIR}/distortion)
add_library(image correction.h image.cpp image.h pgm_io.cpp pgm_io.h spline.cpp spline.h)