


add_library(
    header_file_virph STATIC
    header_file_virph.hpp
) 
set_target_properties(header_file_virph PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(header_file_virph PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


add_target_with_lib(virph virph.cpp)
target_link_libraries(virph PUBLIC header_file_virph)
