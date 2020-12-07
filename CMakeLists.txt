cmake_minimum_required(VERSION 3.13)
project(analysis_program CXX)
#add_subdirectory(../modules  ../modules)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "-fopenmp   -pedantic  -g -O0 -pg -no-pie  -lm -lgmp -lfftw3 ")
#-g # then compile , then open gdb, write the commands : file $namefile; run  $options ; bt
# valgrind --leak-check=yes --track-origins=yes  ./program   #compile with -g to have line number


# no need to install Eigen with cmake and use find package
#SET( EIGEN3_INCLUDE_DIR /home/marco/library/eigen )
#IF( NOT EIGEN3_INCLUDE_DIR )
#    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
#ENDIF()
#INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

#using the globally installed version of Eigen: sudo apt-get libeigen3-dev
# Eigen
IF(NOT DEFINED EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 REQUIRED)
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )





add_subdirectory(include  include)

# this has been taken straight from the build system of Balint's Dslash (github.com/bjoo/KokkosDslash)
macro(add_target_with_lib TARGETNAME)
  # Make the executable
  add_executable(${TARGETNAME} ${ARGN})
  target_link_libraries(${TARGETNAME} PUBLIC global)
  target_link_libraries(${TARGETNAME} PUBLIC rand)
  target_link_libraries(${TARGETNAME} PUBLIC linear_fit)
  target_link_libraries(${TARGETNAME} PUBLIC resampling)
  target_link_libraries(${TARGETNAME} PUBLIC mutils)
 # target_link_libraries (${TARGETNAME} PUBLIC Eigen3::Eigen)
  target_link_libraries (${TARGETNAME} PUBLIC indices)
  target_link_libraries (${TARGETNAME} PUBLIC m_eff)
  target_link_libraries(${TARGETNAME} PUBLIC fftw3)
  target_link_libraries(${TARGETNAME} PUBLIC read)
  target_link_libraries(${TARGETNAME} PUBLIC continuum_reph)
  target_link_libraries(${TARGETNAME} PUBLIC gnuplot)
  target_link_libraries(${TARGETNAME} PUBLIC various_fits)
  target_link_libraries(${TARGETNAME} PUBLIC bessel)
  target_link_libraries(${TARGETNAME} PUBLIC eigensystem)
  target_link_libraries(${TARGETNAME} PUBLIC fve)
  target_link_libraries(${TARGETNAME} PUBLIC global_barion_QED)  
 # target_link_libraries(${TARGETNAME} PUBLIC global_fit_KandD)
  target_link_libraries(${TARGETNAME} PUBLIC global_reph)
  target_link_libraries(${TARGETNAME} PUBLIC KandD)
  target_link_libraries(${TARGETNAME} PUBLIC non_linear_fit)
  target_link_libraries(${TARGETNAME} PUBLIC pion)
  target_link_libraries(${TARGETNAME} PUBLIC m_eff1)
  target_link_libraries(${TARGETNAME} PUBLIC virtualphoton_routines)
  target_link_libraries(${TARGETNAME} PUBLIC gamma_analysis)
#target_link_libraries(${TARGETNAME} PUBLIC Kokkos::kokkos)
endmacro()


add_subdirectory(projects/realphoton  realphoton)
add_subdirectory(projects/clover  clover)
add_subdirectory(projects/virtualphoton  virtualphoton)
add_subdirectory(tools tools)
add_subdirectory(projects/phi4 phi4)



add_subdirectory(tests  tests)

