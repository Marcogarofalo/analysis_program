



macro(add_target_gm2 TARGETNAME)
    # Make the executable
    add_target_with_lib(${TARGETNAME} ${ARGN})
    target_link_libraries(${TARGETNAME} PUBLIC dzeta_function)

endmacro()

add_target_gm2(g-2 g-2.cpp)
add_target_gm2(fit_all_g-2 fit_all_g-2.cpp)

add_library(
    systematics_charm STATIC
    do_analysis_charm.hpp
)
target_link_libraries(systematics_charm PUBLIC mutils)
target_include_directories(systematics_charm PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(fit_all_g-2 PUBLIC systematics_charm)

