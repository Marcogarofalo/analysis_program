add_library(
    functions_renormalization_D STATIC
    functions_renormalization_D.hpp
    functions_renormalization_D.cpp
)
# target_link_libraries(functions_renormalization_D PUBLIC non_linear_fit)
# target_include_directories(functions_renormalization_D PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(functions_renormalization_D  PUBLIC analysis_program)
target_include_directories(common_integral_eq PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_target_with_lib(renormalization_D renormalization_D.cpp)
target_link_libraries(renormalization_D PUBLIC functions_renormalization_D)
