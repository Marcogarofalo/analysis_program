



add_target_with_lib(test_covariance_generator test_covariance_generator.cpp)
add_target_with_lib(test_linear_fit test_linear_fit.cpp)
add_target_with_lib(test_non_linear_fit_sigmax test_non_linear_fit_sigmax.cpp)
add_target_with_lib(test_non_linear_fit test_non_linear_fit.cpp)
add_target_with_lib(fit_as_minimiser fit_as_minimiser.cpp)
add_target_with_lib(test_fake_resempling test_fake_resempling.cpp)
if(WITH_ARB)
    add_target_with_lib(test_HLT test_HLT.cpp)
ENDIF()

add_target_with_lib(benchmark benchmark.cpp)

add_target_with_lib(test_GEVP test_GEVP.cpp)



