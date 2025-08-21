cmake_minimum_required(VERSION 3.13)
#add_subdirectory(../modules  ../modules)


MESSAGE("CMAKE_CXX_FLAGS:   ${CMAKE_CXX_FLAGS} ")
#set(CMAKE_C_FLAGS "-fopenmp   -pedantic  -g -O3  -lm -lgmp  ")

project(analysis_program  VERSION 0.1.0 LANGUAGES CXX C)
#-pg -no-pie , for gprof, then run , gprof exe  gmout.out > output.txt.
#-g # then compile , then open gdb, write the commands : file $namefile; run  $options ; bt
# valgrind --leak-check=yes --track-origins=yes  ./program   #compile with -g to have line number

# Enable C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # optional: disables compiler-specific extensio
MESSAGE("using CXX:" ${CMAKE_CXX_STANDARD} )


IF(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment")
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -static-libasan -g   ")#-fsanitize=float-divide-by-zero
  set(LDFLAGS "${LDFLAGS}  -fsanitize=address")
  set(CMAKE_LDFLAGS_FLAGS  " ${CMAKE_LDFLAGS_FLAGS} ${LDFLAGS} ")
  #
ENDIF()

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
  MESSAGE("Eigen: " ${EIGEN3_INCLUDE_DIR})
ENDIF()


#using the globally installed version of Eigen: sudo apt-get libeigen3-dev
# Eigen
IF(NOT DEFINED GSL_INCLUDE_DIR)
  find_package(GSL)
  IF(NOT GSL_FOUND) 
    MESSAGE("\nGSL NOT FOUND TRY with PkgConfig otherwise set manually:")
    MESSAGE("        -DCMAKE_EXE_LINKER_FLAGS=\"-L $GSL_LIB \" " )
    MESSAGE("        -DCMAKE_C_FLAGS=\"$GSL_INCLUDE \"\n " )
    find_package(PkgConfig REQUIRED)
    pkg_check_modules( GSL  gsl )
  ENDIF()
ENDIF()
# if (PkgConfig_FOUND)
#     pkg_check_modules(GSL gsl)
# endif ()


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



IF(WITH_ARB)
  if(NOT DEFINED ARB_ROOT)
        message(FATAL_ERROR "WITH_ARB=ON but ARB_ROOT is not specified. Use -DARB_ROOT=/path/to/flint")
  endif()
  set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${CXXFLAGS} -DWITH_ARB  ")#-larb -lflint -lgmp -lmpfr 
  # Find GMP and MPFR (required by Flint/Arb)
  # find_package(GMP REQUIRED)
  # find_package(MPFR REQUIRED)
  
  # Include directories
  include_directories(
      ${ARB_ROOT}/include
      ${ARB_ROOT}/include/flint
  )

  # Library directories
  link_directories(${ARB_ROOT}/lib)
  MESSAGE(${ARB_ROOT}/lib)
  set(ARB_LIBS flint)

  # include_directories(${CMAKE_INCLUDE_PATH})
ENDIF()


add_subdirectory(include  include)
#add_subdirectory(external/rzeta/src  zeta)
#############################################
# install package
#############################################
include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)
include(GNUInstallDirs)


target_include_directories(analysis_program 
PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include;${PROJECT_SOURCE_DIR}/external/rzeta/src>"
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
PRIVATE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
) 
target_include_directories(zeta_function 
PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external/rzeta/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
PRIVATE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
) 



set(CONFIG_EXPORT_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/analysis_program")
set(EXPORT_TARGETS analysis_program zeta_function)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/analysis_program-config.cmake.in"
  "${PROJECT_BINARY_DIR}/analysis_program-config.cmake"
  INSTALL_DESTINATION "${CONFIG_EXPORT_DIR}"
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR CONFIG_EXPORT_DIR)
unset(EXPORT_TARGETS)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/analysis_program-config-version.cmake"
  COMPATIBILITY AnyNewerVersion)

configure_file("${PROJECT_SOURCE_DIR}/analysis_program.pc.in" "${PROJECT_BINARY_DIR}/analysis_program.pc" @ONLY)

MESSAGE("PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR}")
MESSAGE("CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}")
MESSAGE("CONFIG_EXPORT_DIR ${CONFIG_EXPORT_DIR}")
MESSAGE("CMAKE_INSTALL_LIBDIR  ${CMAKE_INSTALL_LIBDIR}")
MESSAGE("EXPORT_TARGETS  ${EXPORT_TARGETS}")

install(TARGETS analysis_program zeta_function
  EXPORT analysis_program-targets 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
                  ${PROJECT_SOURCE_DIR}/external/rzeta/src/
                  ${PROJECT_SOURCE_DIR}/projects/g-2/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING 
              PATTERN "*.hpp"
              PATTERN "*.h")
install(EXPORT analysis_program-targets
  DESTINATION "${CONFIG_EXPORT_DIR}")
install(FILES
    "${PROJECT_BINARY_DIR}/analysis_program-config.cmake"
    "${PROJECT_BINARY_DIR}/analysis_program-config-version.cmake"
  DESTINATION "${CONFIG_EXPORT_DIR}")
install(FILES "${PROJECT_BINARY_DIR}/analysis_program.pc"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

unset(CONFIG_EXPORT_DIR)


#############################################

# this has been taken straight from the build system of Balint's Dslash (github.com/bjoo/KokkosDslash)
macro(add_target_with_lib TARGETNAME)
  # Make the executable
  add_executable(${TARGETNAME} ${ARGN})
  # links
  target_link_libraries(${TARGETNAME}  PUBLIC analysis_program)
endmacro()

if(realphoton_PROJECT OR ALL_PROJECT)
  add_subdirectory(projects/realphoton  realphoton)
ENDIF()
if(clover_PROJECT OR ALL_PROJECT)
  add_subdirectory(projects/clover  clover)
ENDIF()
if(virtualphoton_PROJECT OR ALL_PROJECT)
  add_subdirectory(projects/virtualphoton  virtualphoton)
ENDIF()
add_subdirectory(tools tools)
add_subdirectory(projects/phi4 phi4)
if(BSM_PROJECT OR ALL_PROJECT)
  add_subdirectory(projects/BSM BSM)
ENDIF()
if(GM2_PROJECT OR ALL_PROJECT)
  add_subdirectory(projects/g-2 g-2)
ENDIF()
if(NISSA_MPCAC_PROJECT OR ALL_PROJECT)
  add_subdirectory(projects/nissa_mpcac nissa_mpcac)
ENDIF()

if(WITH_PYTHON)
    add_subdirectory(projects/python_wrapper python_wrapper)
ENDIF()



add_subdirectory(tests  tests)

IF(renormalization_D_PROJECT OR ALL_PROJECT)
    add_subdirectory(projects/renormalization_D  renormalization_D)
ENDIF()
