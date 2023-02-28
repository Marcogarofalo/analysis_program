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


