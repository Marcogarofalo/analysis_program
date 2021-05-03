

add_library(
    init_zeta_function SHARED
    init_zeta_function.h
    init_zeta_function.c
    
) 
set_target_properties(init_zeta_function PROPERTIES LINKER_LANGUAGE C)
target_include_directories(init_zeta_function PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
 


add_library(
    dzeta_function SHARED
    dzeta_function.h
    dzeta_function.c
    #qzeta_function.h
    #qzeta_function.c
) 
set_target_properties(dzeta_function PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(dzeta_function PUBLIC init_zeta_function)
target_link_libraries(dzeta_function PUBLIC gsl)
target_include_directories(dzeta_function PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
 
