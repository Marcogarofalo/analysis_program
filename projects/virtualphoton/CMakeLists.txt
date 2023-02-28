


add_library(
    header_file_virph STATIC
    header_file_virph.hpp
) 
set_target_properties(header_file_virph PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(header_file_virph PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")


macro(add_target_virtualphoton TARGETNAME)
    # Make the executable
    add_target_with_lib(${TARGETNAME} ${ARGN})
    # target_link_libraries(${TARGETNAME} PUBLIC global_reph)
    # target_link_libraries(${TARGETNAME} PUBLIC continuum_reph)
    # target_link_libraries(${TARGETNAME} PUBLIC header_file_virph)
endmacro()


add_target_virtualphoton(virph virph.cpp)
