

macro(add_target_with_lib_clover TARGETNAME)
    # Make the executable
    add_target_with_lib(${TARGETNAME} ${ARGN})
    # target_link_libraries(${TARGETNAME} PUBLIC bessel)
    # target_link_libraries(${TARGETNAME} PUBLIC pion)
    # target_link_libraries(${TARGETNAME} PUBLIC KandD)
    # target_link_libraries(${TARGETNAME} PUBLIC fve)
endmacro()




add_target_with_lib_clover(binning binning.cpp)
add_target_with_lib_clover(convert convert.cpp)
add_target_with_lib_clover(fit_all_beta fit_all_beta.cpp)
add_target_with_lib_clover(fit_all_clover fit_all_clover.cpp)
add_target_with_lib_clover(form_factors form_factors.cpp)
add_target_with_lib_clover(form_factors_out_max_twist form_factors_out_max_twist.cpp)
add_target_with_lib_clover(petros_correlators petros_correlators.cpp)
add_target_with_lib_clover(prace2gm2 prace2gm2.cpp)




