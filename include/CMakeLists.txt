
add_library(
    global STATIC
    global.hpp
) 
set_target_properties(global PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(global PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(
    tower STATIC
    tower.hpp
    ../modules/tower.cpp
) 

target_include_directories(tower PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")




add_library(
    resampling STATIC
    resampling.hpp
    ../modules/resampling.cpp
) 

target_link_libraries(resampling PUBLIC tower)
target_link_libraries(resampling PUBLIC rand)
target_include_directories(resampling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(
    mutils STATIC
    mutils.hpp
    ../modules/mutils.cpp
) 

target_include_directories(mutils PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    rand STATIC
    rand.hpp
    ../modules/rand.cpp
) 

target_include_directories(rand PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    linear_fit STATIC
    linear_fit.hpp
    ../modules/linear_fit.cpp
) 
target_link_libraries(linear_fit PUBLIC resampling)
target_include_directories(linear_fit PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    non_linear_fit STATIC
    non_linear_fit.hpp
    ../modules/non_linear_fit.cpp
    ../modules/non_linear_fit_sigmax.cpp
) 
target_link_libraries(non_linear_fit PUBLIC linear_fit)
target_link_libraries(non_linear_fit PUBLIC mutils)
target_include_directories(non_linear_fit PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    read STATIC
    read.hpp
    ../modules/read.cpp
) 

target_include_directories(read PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    m_eff STATIC
    m_eff.hpp
    ../modules/m_eff.cpp
) 

target_include_directories(m_eff PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    gnuplot STATIC
    gnuplot.hpp
    ../modules/gnuplot.cpp
) 

target_include_directories(gnuplot PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    eigensystem STATIC
    eigensystem.hpp
    ../modules/eigensystem.cpp
) 

target_include_directories(eigensystem PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    various_fits STATIC
    various_fits.hpp
    ../modules/various_fits.cpp
) 
target_link_libraries(various_fits PUBLIC indices)
target_link_libraries(various_fits PUBLIC linear_fit)
target_link_libraries(various_fits PUBLIC gnuplot)
target_include_directories(various_fits PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    indices STATIC
    indices.hpp
    ../modules/indices.cpp
) 

target_include_directories(indices PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    continuum_reph STATIC
    continuum_reph.hpp
    ../modules/reph/FAandV.cpp
    ../modules/reph/automatic_plateau_R.cpp
    ../modules/reph/plot_reph.cpp 
    ../modules/reph/FAV_interpolation_msc.cpp
    ../modules/reph/FAV_treshold.cpp
) 

target_include_directories(continuum_reph PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    bessel STATIC
    bessel.hpp
    ../modules/bessel/dbesk0.cpp 
    ../modules/bessel/dbesk1.cpp
) 

target_include_directories(bessel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    fve STATIC
    fve.hpp
    ../modules/FVE/fve.cpp
    ../modules/FVE/fve_K.cpp
#    ../modules/FVE/fve_Mpi.cpp
    
) 
target_link_libraries(fve PUBLIC bessel)
target_include_directories(fve PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(
    global_barion_QED STATIC
    global_barion_QED.hpp
) 
set_target_properties(global_barion_QED PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(global_barion_QED PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")





add_library(
    global_reph STATIC
    global_reph.hpp
) 
set_target_properties(global_reph PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(global_reph PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    KandD STATIC
    KandD.hpp
    ../modules/KandD/fKoverfpi_clover.cpp
    ../modules/KandD/global_fit_KandD.cpp
    ../modules/KandD/KandD_clover.cpp
    ../modules/KandD/KandD.cpp
    ../modules/KandD/KandD_from_M.cpp
    ../modules/KandD/K_chiral_fit.cpp
    ../modules/KandD/Omega_jacob.cpp
    ../modules/KandD/KandD.cpp
) 
target_link_libraries(KandD PUBLIC fve)
target_include_directories(KandD PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    pion STATIC
    pion.hpp
    ../modules/pion/Pion_clover.cpp
    ../modules/pion/Pion_clover_treshold.cpp
    ../modules/pion/pion_FVE.cpp
    ../modules/pion/pion.cpp
) 

target_include_directories(pion PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    m_eff1 STATIC
    m_eff1.hpp
    ../modules/m_eff1.cpp
) 

target_include_directories(m_eff1 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    virtualphoton_routines STATIC
    virtualphoton_routines.hpp
    ../modules/virtualphoton/routines_H.cpp
) 

target_include_directories(virtualphoton_routines PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_library(
    gamma_analysis STATIC
    gamma_analysis.hpp
    ../modules/Gamma.cpp
) 

target_include_directories(gamma_analysis PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



add_library(
    correlators_analysis STATIC
    correlators_analysis.hpp
    ../modules/correlators_analysis.cpp
) 
target_link_libraries(correlators_analysis PUBLIC non_linear_fit)
target_include_directories(correlators_analysis PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")




