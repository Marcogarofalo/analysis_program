add_library(
    functions_project_template STATIC
    functions_project_template.hpp
    functions_project_template.cpp
)
target_link_libraries(functions_project_template PUBLIC non_linear_fit)
target_include_directories(functions_project_template PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_target_with_lib(project_template project_template.cpp)
target_link_libraries(project_template PUBLIC functions_project_template)
