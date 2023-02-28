cmake_minimum_required(VERSION 3.13)

project(project_template CXX C)
#############################################
# load analysis_program
#############################################
find_package(analysis_program)
include_directories(${ANALYSIS_PROGRAM_INCLUDE_DIR})
MESSAGE("ANALYSIS_PROGRAM_INCLUDE_DIR   ${ANALYSIS_PROGRAM_INCLUDE_DIR}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARIES: ${ANALYSIS_PROGRAM_LIBRARIES}")
MESSAGE("ANALYSIS_PROGRAM_LIBRARY_DIR: ${ANALYSIS_PROGRAM_LIBRARY_DIR}")
IF(NOT analysis_program_FOUND)
    MESSAGE("find_package faild to find yaml-cpp: you need to do it manually")
    MESSAGE("hint:  analysis_program_FOUND_DIR or CMAKE_PREFIX_PATH need to be global path or relative to the source dir")
    MESSAGE(FATAL_ERROR)
ENDIF()
find_package(analysis_program)

find_package(Eigen3 REQUIRED) # we need eigen

#############################################
# add a local library for this project
#############################################
add_library(
    functions_project_template STATIC
    functions_project_template.hpp
    functions_project_template.cpp
)
target_link_libraries(functions_project_template PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})
target_include_directories(functions_project_template PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(project_template project_template.cpp)
target_link_libraries(project_template PUBLIC functions_project_template)
target_link_libraries(project_template PUBLIC ${ANALYSIS_PROGRAM_LIBRARIES})
