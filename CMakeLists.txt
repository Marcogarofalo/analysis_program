cmake_minimum_required(VERSION 3.13)
#add_subdirectory(../modules  ../modules)

#set(CMAKE_CXX_COMPILER /usr/bin/g++)
#set(CMAKE_C_COMPILER /usr/bin/gcc)
#set(CMAKE_CXX_FLAGS "-fopenmp   -pedantic  -g -O3   -lm -lgmp  ")
MESSAGE("CMAKE_CXX_FLAGS:   ${CMAKE_CXX_FLAGS} ")
set(CMAKE_C_FLAGS "-fopenmp   -pedantic  -g -O3  -lm -lgmp  ")

project(analysis_program CXX C)
#-pg -no-pie , for gproof, then run , gproof exe  gmout.
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


IF(WITH_PYTHON)
    SET(python_dir "/usr/include/python3.9/")  #sudo apt-get install python3.9-dbg
    #find_package( PythonInterp REQUIRED )
    #find_package(PythonLibs REQUIRED)
    find_package (Python3 COMPONENTS Interpreter  Development )
    
    include_directories(${Python3_INCLUDE_DIRS})
   # include_directories(${PYTHON_INCLUDE_DIRS})
    SET(Python3_LINK_OPTIONS "/usr/lib/python3.9/config-3.9-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu  -lcrypt -lpthread -ldl  -lutil -lm -lm")
    set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -DPYTHON   ")
    MESSAGE("link options: " ${Python3_LINK_OPTIONS})
ENDIF()


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
  target_link_libraries(${TARGETNAME} PUBLIC correlators_analysis)
  target_link_libraries(${TARGETNAME} PUBLIC tower)
  target_link_libraries(${TARGETNAME} PUBLIC init_zeta_function)
  target_link_libraries(${TARGETNAME} PUBLIC dzeta_function)
  #target_link_libraries(${TARGETNAME} PUBLIC header_phi4)
  target_link_libraries(${TARGETNAME} PUBLIC zeta_interpolation)
  target_include_directories(${TARGETNAME} PUBLIC "${CMAKE_SOURCE_DIR}/include")
  target_include_directories(${TARGETNAME} PUBLIC "${CMAKE_BINARY_DIR}/include")
  IF(WITH_PYTHON)
        #target_include_directories(${TARGETNAME} PUBLIC pyhelper)
        #target_include_directories(${TARGETNAME} PUBLIC python3.9)
        target_link_libraries(${TARGETNAME} PUBLIC ${PYTHON_LIBRARIES})
        target_link_libraries(${TARGETNAME} PUBLIC QC3_interface)
        target_link_libraries(${TARGETNAME} PUBLIC ${Python3_LIBRARIES})
        
  ENDIF()
#target_link_libraries(${TARGETNAME} PUBLIC Kokkos::kokkos)
endmacro()


add_subdirectory(projects/realphoton  realphoton)
add_subdirectory(projects/clover  clover)
add_subdirectory(projects/virtualphoton  virtualphoton)
add_subdirectory(tools tools)
add_subdirectory(projects/phi4 phi4)
add_subdirectory(projects/BSM BSM)
add_subdirectory(projects/g-2 g-2)

if(WITH_PYTHON)
    add_subdirectory(projects/python_wrapper python_wrapper)
ENDIF()



add_subdirectory(tests  tests)

