#!/bin/bash

if [ -z "$1" ]
then
    echo "usage: create_project.sh \$name"
    exit 1
fi
name=$1

grep -q "$name" CMakeLists.txt && echo "project already exist" && exit 1

echo "IF($name)" >> CMakeLists.txt
echo "    add_subdirectory(projects/$name  $name)" >> CMakeLists.txt
echo "ENDIF()" >> CMakeLists.txt

cp -r  projects/project_template  projects/$name
sed -i "s/project_template/$name/g"  projects/$name/CMakeLists.txt

mv -v projects/$name/functions_project_template.hpp     projects/$name/functions_$name.hpp
mv -v projects/$name/functions_project_template.cpp     projects/$name/functions_$name.cpp
mv -v projects/$name/project_template.cpp     projects/$name/$name.cpp

sed -i "s/project_template/$name/g"  projects/$name/functions_$name.hpp 
sed -i "s/project_template/$name/g"  projects/$name/functions_$name.cpp
sed -i "s/project_template/$name/g"  projects/$name/$name.cpp
# cat >projects/$name/CMakeLists.txt << EOF
# add_library(
#     functions_$name STATIC
#     functions_$name.hpp
#     functions_$name.cpp
# )
# add_target_with_lib(${name} ${name}.cpp)
# target_link_libraries(${name} PUBLIC functions_$name)
# EOF

# cat >projects/$name/functions_$name.hpp << EOF
# #ifndef functions_${name}_H
# #define functions_${name}_H

# #endif
# EOF

# cat >projects/$name/functions_$name.cpp << EOF
# #define functions_${name}_C
# #include "functions_$name.hpp"
# EOF

# cat > projects/$name/$name.cpp<< EOF
# #define CONTROL
# #include <stdlib.h>
# #include <stdio.h>
# #include <math.h>
# #include <time.h>
# #include <string.h>
# #include <complex.h>
# #include <iostream>

# #include "functions_$name.hpp"
# EOF
