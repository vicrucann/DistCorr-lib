cmake_minimum_required(VERSION 2.6)
include_directories(${CMAKE_SOURCE_DIR}/distortion ${CMAKE_SOURCE_DIR}/libImage)
add_library(linedetect devernay.cpp devernay.h gauss.cpp gauss.h gaussian_convol_on_curve.cpp gaussian_convol_on_curve.h lsd.cpp lsd.h ntuple.cpp ntuple.h ntuple_ll.cpp ntuple_ll.h straight_edge_points.cpp straight_edge_points.h)