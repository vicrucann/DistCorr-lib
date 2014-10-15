cmake_minimum_required(VERSION 2.6)
include_directories(${CMAKE_SOURCE_DIR}/distortion ${CMAKE_SOURCE_DIR}/libNumerics)
add_library(distortion distortion.cpp distortion.h LM.h)