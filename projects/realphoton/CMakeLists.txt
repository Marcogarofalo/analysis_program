



macro(add_target_realphoton TARGETNAME)
    # Make the executable
    add_target_with_lib(${TARGETNAME} ${ARGN})
    target_link_libraries(${TARGETNAME} PUBLIC global_reph)
    target_link_libraries(${TARGETNAME} PUBLIC continuum_reph)
endmacro()

add_target_realphoton(reph reph.cpp)
add_target_realphoton(fit_all_reph fit_all_reph.cpp)
add_target_realphoton(binning_reph binning_reph.cpp)

