


#add_target_with_lib(phi4 phi4.cpp)

add_target_with_lib(test test.cpp)
add_target_with_lib(phi4_new phi4_new.cpp)
IF(WITH_PYTHON)
    add_target_with_lib(fit_all_quantization_condition fit_all_quantization_condition.cpp)
    add_target_with_lib(fit_all_phi4 fit_all_phi4.cpp)
ENDIF()
add_target_with_lib(checks_phi4 checks_phi4.cpp)




