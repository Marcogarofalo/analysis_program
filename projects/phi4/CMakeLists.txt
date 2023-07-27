# add_target_with_lib(phi4 phi4.cpp)
add_target_with_lib(test test.cpp)
add_target_with_lib(phi4_new phi4_new.cpp)
add_target_with_lib(analysis_SDE analysis_SDE.cpp)

IF(PHI4_PROJECT OR ALL_PROJECT)
    IF(WITH_PYTHON)
        add_target_with_lib(fit_all_g0.25 fit_all_g0.25.cpp)
        add_target_with_lib(fit_all_g1.5 fit_all_g1.5.cpp)
        add_target_with_lib(fit_all_g0.25_t1 fit_all_g0.25_t1.cpp)
        add_target_with_lib(fit_all_ising_g0.005 fit_all_ising_g0.005.cpp)

        add_library(
            extra_func_phi4 STATIC
            extra_func_phi4.hpp
            extra_func_phi4.cpp
        )
        target_link_libraries(extra_func_phi4 PUBLIC analysis_program)
        target_include_directories(extra_func_phi4 PUBLIC ${CMAKE_SOURCE_DIR}/include)
        target_include_directories(extra_func_phi4 PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

        add_target_with_lib(fit_all_complex_g0 fit_all_complex_g0.cpp)
        target_link_libraries(fit_all_complex_g0 PUBLIC extra_func_phi4)

        add_target_with_lib(fit_all_complex_g5 fit_all_complex_g5.cpp)
        target_link_libraries(fit_all_complex_g5 PUBLIC extra_func_phi4)

        add_target_with_lib(fit_all_complex_g5_largeL fit_all_complex_g5_largeL.cpp)
        target_link_libraries(fit_all_complex_g5_largeL PUBLIC extra_func_phi4)

        add_target_with_lib(fit_all_complex_g10 fit_all_complex_g10.cpp)
        target_link_libraries(fit_all_complex_g10 PUBLIC extra_func_phi4)

        add_target_with_lib(fit_all_complex_g20 fit_all_complex_g20.cpp)
        target_link_libraries(fit_all_complex_g20 PUBLIC extra_func_phi4)

        add_library(
            common_integral_eq STATIC
            common_integral_eq.hpp
            common_integral_eq.cpp
        )
        target_link_libraries(common_integral_eq PUBLIC analysis_program)
        target_include_directories(common_integral_eq PUBLIC ${CMAKE_SOURCE_DIR}/include)

        macro(add_target_with_int_eq_lib TARGETNAME)
            # Make the executable
            add_target_with_lib(${TARGETNAME} ${ARGN})

            # links
            target_link_libraries(${TARGETNAME} PUBLIC common_integral_eq)
        endmacro()

        add_target_with_int_eq_lib(integral_eq_g5_largeL integral_eq_g5_largeL.cpp)
        add_target_with_int_eq_lib(integral_eq_g5_2par_largeL integral_eq_g5_2par_largeL.cpp)
        add_target_with_int_eq_lib(integral_eq_g10 integral_eq_g10.cpp)
        add_target_with_int_eq_lib(integral_eq_g10_2par integral_eq_g10_2par.cpp)
        add_target_with_int_eq_lib(integral_eq_g20 integral_eq_g20.cpp)
        add_target_with_int_eq_lib(integral_eq_g20_2par integral_eq_g20_2par.cpp)
        add_target_with_int_eq_lib(integral_eq_scratchpad integral_eq_scratchpad.cpp)
        add_target_with_int_eq_lib(checks_phi4 checks_phi4.cpp)
    ENDIF()

    
ENDIF()