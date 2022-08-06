



macro(add_target_gm2 TARGETNAME)
    # Make the executable
    add_target_with_lib(${TARGETNAME} ${ARGN})
    target_link_libraries(${TARGETNAME} PUBLIC dzeta_function)

endmacro()

add_target_gm2(g-2 g-2.cpp)
add_target_gm2(fit_all_g-2 fit_all_g-2.cpp)


