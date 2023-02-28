
macro(include_and_install TARGETNAME)
    target_include_directories(${TARGETNAME} PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>  
        $<INSTALL_INTERFACE:include>  # <prefix>/include/mylib
    )
    install( 
        TARGETS ${TARGETNAME}
        EXPORT analysis_program-targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endmacro()

# add_library(
#     global STATIC
#     global.hpp
# ) 
# set_target_properties(global PROPERTIES LINKER_LANGUAGE CXX)
# include_and_install(global)

# add_library(
#     tower STATIC
#     tower.hpp
#     ../modules/tower.cpp
# ) 

# include_and_install(tower)



# add_library(
#     resampling STATIC
#     resampling.hpp
#     ../modules/resampling.cpp
# ) 

# target_link_libraries(resampling PUBLIC tower)
# target_link_libraries(resampling PUBLIC rand)
# target_link_libraries(resampling PUBLIC resampling_new)
# include_and_install(resampling)


# add_library(
#     resampling_new STATIC
#     resampling_new.hpp
#     ../modules/resampling_new.cpp
# ) 

# target_link_libraries(resampling_new PUBLIC tower)
# target_link_libraries(resampling_new PUBLIC rand)
# include_and_install(resampling_new)


# add_library(
#     mutils STATIC
#     mutils.hpp
#     ../modules/mutils.cpp
# ) 

# include_and_install(mutils)

# add_library(
#     rand STATIC
#     rand.hpp
#     ../modules/rand.cpp
# ) 
# target_link_libraries(rand PUBLIC linear_fit)
# include_and_install(rand)



# add_library(
#     linear_fit STATIC
#     linear_fit.hpp
#     ../modules/linear_fit.cpp
# ) 
# target_link_libraries(linear_fit PUBLIC resampling)
# target_link_libraries(linear_fit PUBLIC non_linear_fit)
# include_and_install(linear_fit)


# add_library(
#     non_linear_fit STATIC
#     non_linear_fit.hpp
#     ../modules/non_linear_fit.cpp
#     ../modules/non_linear_fit_sigmax.cpp
# )
# target_link_libraries(non_linear_fit PUBLIC fit_all) 
# target_link_libraries(non_linear_fit PUBLIC linear_fit)
# target_link_libraries(non_linear_fit PUBLIC mutils)
# target_link_libraries(non_linear_fit PUBLIC sorting)
# include_and_install(non_linear_fit)


# add_library(
#     minimizer STATIC
#     minimizer.hpp
#     ../modules/minimizer.cpp
# )
# target_link_libraries(minimizer PUBLIC non_linear_fit)
# include_and_install(minimizer)


# add_library(
#     read STATIC
#     read.hpp
#     ../modules/read.cpp
# ) 

# include_and_install(read)



# add_library(
#     m_eff STATIC
#     m_eff.hpp
#     ../modules/m_eff.cpp
# ) 

# include_and_install(m_eff)



# add_library(
#     gnuplot STATIC
#     gnuplot.hpp
#     ../modules/gnuplot.cpp
# ) 

# include_and_install(gnuplot)


# add_library(
#     eigensystem STATIC
#     eigensystem.hpp
#     ../modules/eigensystem.cpp
# ) 

# target_link_libraries(eigensystem PUBLIC sorting)
# include_and_install(eigensystem)


# add_library(
#     various_fits STATIC
#     various_fits.hpp
#     ../modules/various_fits.cpp
# ) 
# target_link_libraries(various_fits PUBLIC indices)
# target_link_libraries(various_fits PUBLIC linear_fit)
# target_link_libraries(various_fits PUBLIC gnuplot)
# include_and_install(various_fits)


# add_library(
#     indices STATIC
#     indices.hpp
#     ../modules/indices.cpp
# ) 

# include_and_install(indices)


# add_library(
#     continuum_reph STATIC
#     continuum_reph.hpp
#     ../modules/reph/FAandV.cpp
#     ../modules/reph/automatic_plateau_R.cpp
#     ../modules/reph/plot_reph.cpp 
#     ../modules/reph/FAV_interpolation_msc.cpp
#     ../modules/reph/FAV_treshold.cpp
# ) 

# include_and_install(continuum_reph)



# add_library(
#     bessel STATIC
#     bessel.hpp
#     ../modules/bessel/dbesk0.cpp 
#     ../modules/bessel/dbesk1.cpp
# ) 

# include_and_install(bessel)



# add_library(
#     fve STATIC
#     fve.hpp
#     ../modules/FVE/fve.cpp
#     ../modules/FVE/fve_K.cpp
# #    ../modules/FVE/fve_Mpi.cpp
    
# ) 
# target_link_libraries(fve PUBLIC bessel)
# include_and_install(fve)

# add_library(
#     global_barion_QED STATIC
#     global_barion_QED.hpp
# ) 
# set_target_properties(global_barion_QED PROPERTIES LINKER_LANGUAGE CXX)
# include_and_install(global_barion_QED)





# add_library(
#     global_reph STATIC
#     global_reph.hpp
# ) 
# set_target_properties(global_reph PROPERTIES LINKER_LANGUAGE CXX)
# target_lin(global_reph)


# add_library(
#     KandD STATIC
#     KandD.hpp
#     ../modules/KandD/fKoverfpi_clover.cpp
#     ../modules/KandD/global_fit_KandD.cpp
#     ../modules/KandD/KandD_clover.cpp
#     ../modules/KandD/KandD.cpp
#     ../modules/KandD/KandD_from_M.cpp
#     ../modules/KandD/K_chiral_fit.cpp
#     ../modules/KandD/Omega_jacob.cpp
#     ../modules/KandD/KandD.cpp
# ) 
# target_link_libraries(KandD PUBLIC fve)
# target_include_directories(KandD PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")



# add_library(
#     pion STATIC
#     pion.hpp
#     ../modules/pion/Pion_clover.cpp
#     ../modules/pion/Pion_clover_treshold.cpp
#     ../modules/pion/pion_FVE.cpp
#     ../modules/pion/pion.cpp
# ) 

# target_include_directories(pion PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


# add_library(
#     m_eff1 STATIC
#     m_eff1.hpp
#     ../modules/m_eff1.cpp
# ) 

# include_and_install(m_eff1)



# add_library(
#     virtualphoton_routines STATIC
#     virtualphoton_routines.hpp
#     ../modules/virtualphoton/routines_H.cpp
# ) 

# target_include_directories(virtualphoton_routines PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


# add_library(
#     gamma_analysis STATIC
#     gamma_analysis.hpp
#     ../modules/Gamma.cpp
# ) 

# include_and_install(gamma_analysis)



# add_library(
#     correlators_analysis STATIC
#     correlators_analysis.hpp
#     ../modules/correlators_analysis.cpp
# ) 
# target_link_libraries(correlators_analysis PUBLIC read)
# target_link_libraries(correlators_analysis PUBLIC non_linear_fit)
# target_link_libraries(correlators_analysis PUBLIC resampling_new)
# set_target_properties(correlators_analysis PROPERTIES LINKER_LANGUAGE CXX)
# include_and_install(correlators_analysis)



# add_library(
#     init_zeta_function STATIC
#     ../external/rzeta/src/init_zeta_function.h
#     ../external/rzeta/src/init_zeta_function.c
    
# ) 
# set_target_properties(init_zeta_function PROPERTIES LINKER_LANGUAGE C)
# include_and_install(init_zeta_function)
 
# add_library(
#     dzeta_function STATIC
#     ../external/rzeta/src/dzeta_function.h
#     ../external/rzeta/src/dzeta_function.c
#     ../external/rzeta/src/qzeta_function.h
#     ../external/rzeta/src/qzeta_function.c
#     #qzeta_function.h
#     #qzeta_function.c
# ) 
# set_target_properties(dzeta_function PROPERTIES LINKER_LANGUAGE C)
# target_link_libraries(dzeta_function PUBLIC quadmath)
# target_link_libraries(dzeta_function PUBLIC init_zeta_function)
# target_link_libraries(dzeta_function PUBLIC gsl gslcblas)
# include_and_install(dzeta_function)
 

# add_library(
#     header_phi4 STATIC
#     header_phi4.hpp
#     ../modules/header_phi4.cpp
# )
# target_link_libraries(header_phi4 PUBLIC fit_all)
# set_target_properties(header_phi4 PROPERTIES LINKER_LANGUAGE CXX)
# include_and_install(header_phi4)
 


# add_library(
#     fit_all STATIC
#     fit_all.hpp
#     ../modules/fit_all.cpp
# )
# target_link_libraries(fit_all PUBLIC non_linear_fit)
# target_link_libraries(fit_all PUBLIC global)
# #set_target_properties(fit_all PROPERTIES LINKER_LANGUAGE CXX)
# # target_include_directories(fit_all PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
# include_and_install(fit_all)

 
# add_library(
#     zeta_interpolation STATIC
#     zeta_interpolation.hpp
#     ../modules/Zeta/zeta_interpolation.cpp
# )
# target_link_libraries(zeta_interpolation PUBLIC dzeta_function)
# target_link_libraries(zeta_interpolation PUBLIC header_phi4)
# include_and_install(zeta_interpolation)

# if(WITH_PYTHON)
#     add_library(
#         QC3_interface STATIC
#         QC3_interface.hpp
#         ../modules/QC3_interface.cpp
#     )
#     set_target_properties(QC3_interface PROPERTIES LINKER_LANGUAGE CXX)
#     target_link_libraries(QC3_interface PUBLIC  ${PYTHON_LIBRARIES})
#     include_and_install(QC3_interface)

# ENDIF()



# add_library(
#     integral_eq_QC3 STATIC
#     integral_eq_QC3.hpp
#     ../modules/three_particle_QC/integral_eq_QC3.cpp
# )
# target_link_libraries(integral_eq_QC3 PUBLIC tower)
# include_and_install(integral_eq_QC3)

# add_library(
#     sorting STATIC
#     sorting.hpp
#     ../modules/sorting.cpp
# ) 
# include_and_install(sorting)
#########################################################

set(analysis_header
    global.hpp 
    tower.hpp
    resampling.hpp
    resampling_new.hpp
    mutils.hpp
    rand.hpp
    linear_fit.hpp
    non_linear_fit.hpp
    minimizer.hpp
    read.hpp
    m_eff.hpp
    gnuplot.hpp
    eigensystem.hpp
    various_fits.hpp
    indices.hpp
    bessel.hpp
    fve.hpp
    global_reph.hpp
    KandD.hpp
    pion.hpp
    m_eff1.hpp
    virtualphoton_routines.hpp
    gamma_analysis.hpp
    correlators_analysis.hpp
    header_phi4.hpp
    fit_all.hpp
    zeta_interpolation.hpp
    sorting.hpp
    ../external/rzeta/src/init_zeta_function.h
    ../external/rzeta/src/dzeta_function.h
    ../external/rzeta/src/qzeta_function.h
    global_reph.hpp
    continuum_reph.hpp
    integral_eq_QC3.hpp


)
set(analysis_scr
    ../modules/tower.cpp
    ../modules/resampling.cpp
    ../modules/resampling_new.cpp
    ../modules/mutils.cpp
    ../modules/rand.cpp
    ../modules/linear_fit.cpp
    ../modules/non_linear_fit.cpp
    ../modules/non_linear_fit_sigmax.cpp
    ../modules/minimizer.cpp
    ../modules/read.cpp
    ../modules/m_eff.cpp
    ../modules/gnuplot.cpp
    ../modules/eigensystem.cpp
    ../modules/various_fits.cpp
    ../modules/indices.cpp
    ../modules/bessel/dbesk0.cpp 
    ../modules/bessel/dbesk1.cpp
    ../modules/FVE/fve.cpp
    ../modules/FVE/fve_K.cpp
    ../modules/KandD/fKoverfpi_clover.cpp
    ../modules/KandD/global_fit_KandD.cpp
    ../modules/KandD/KandD_clover.cpp
    ../modules/KandD/KandD.cpp
    ../modules/KandD/KandD_from_M.cpp
    ../modules/KandD/K_chiral_fit.cpp
    ../modules/KandD/Omega_jacob.cpp
    ../modules/KandD/KandD.cpp
    ../modules/pion/Pion_clover.cpp
    ../modules/pion/Pion_clover_treshold.cpp
    ../modules/pion/pion_FVE.cpp
    ../modules/pion/pion.cpp
    ../modules/m_eff1.cpp
    ../modules/virtualphoton/routines_H.cpp
    ../modules/Gamma.cpp
    ../modules/correlators_analysis.cpp
    ../external/rzeta/src/init_zeta_function.c
    ../external/rzeta/src/dzeta_function.c
    ../external/rzeta/src/qzeta_function.c
    ../modules/header_phi4.cpp
    ../modules/fit_all.cpp
    ../modules/Zeta/zeta_interpolation.cpp
    ../modules/sorting.cpp
    ../modules/reph/FAandV.cpp
    ../modules/reph/automatic_plateau_R.cpp
    ../modules/reph/plot_reph.cpp 
    ../modules/reph/FAV_interpolation_msc.cpp
    ../modules/reph/FAV_treshold.cpp
    ../modules/three_particle_QC/integral_eq_QC3.cpp
)
if(WITH_PYTHON)
  set(analysis_header 
      ${analysis_header}
      QC3_interface.hpp
      )
  set(analysis_scr
      ${analysis_scr}
      ../modules/QC3_interface.cpp
      )
ENDIF()

add_library(analysis_program STATIC ${analysis_header} ${analysis_scr})
target_link_libraries(analysis_program PUBLIC gsl gslcblas)
target_link_libraries(analysis_program PRIVATE Eigen3::Eigen)

if(WITH_PYTHON)
  target_link_libraries(analysis_program PUBLIC ${Python3_LIBRARIES})
ENDIF()
set_target_properties(analysis_program PROPERTIES LINKER_LANGUAGE CXX)


# install( 
#     TARGETS analysis_program
#     EXPORT analysis_program-targets
#     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
# )
