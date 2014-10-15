cmake_minimum_required(VERSION 2.6)

project(Distortion_Correction)

include_directories(${CMAKE_SOURCE_DIR}/libDistortion ${CMAKE_SOURCE_DIR}/libLineDetection ${CMAKE_SOURCE_DIR}/libImage ${CMAKE_SOURCE_DIR}/libNumerics)

add_executable(polyestim polyEstimation.cpp misc.cpp misc.h)
add_executable(distcorrect distCorrection.cpp misc.cpp misc.h)
 
target_link_libraries(polyestim image numerics linedetect distortion)
target_link_libraries(distcorrect image numerics linedetect distortion)