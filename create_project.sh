#!/bin/bash

if [ -z "$1" ]
then
    echo "usage: create_project.sh \$name"
    exit 1
fi
name=${1/_PROJECT/}

grep -q "${name}_PROJECT" CMakeLists.txt && echo "project already exist" && exit 1

echo "IF(${name}_PROJECT OR ALL_PROJECT)" >> CMakeLists.txt
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
