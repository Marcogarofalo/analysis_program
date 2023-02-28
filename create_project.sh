#!/bin/bash

if [ -z "$1" ]; then
    echo "usage: create_project.sh \$name"
    exit 1
fi
name=${1/_PROJECT/}

ls .. | grep -q "${name}" && echo "project already exist" && exit 1

cp -r projects/project_template ../$name
sed -i "s/project_template/$name/g" ../$name/CMakeLists.txt

mv -v ../$name/functions_project_template.hpp ../$name/functions_$name.hpp
mv -v ../$name/functions_project_template.cpp ../$name/functions_$name.cpp
mv -v ../$name/project_template.cpp ../$name/$name.cpp

sed -i "s/project_template/$name/g" ../$name/functions_$name.hpp
sed -i "s/project_template/$name/g" ../$name/functions_$name.cpp
sed -i "s/project_template/$name/g" ../$name/$name.cpp

mkdir ../$name/build

here=$(pwd)
cat > ../$name/build/do_cmake.sh << EOF
cmake -DCMAKE_PREFIX_PATH=${here}/build/install_dir \\
      -DCMAKE_CXX_FLAGS="-g -O2"  \\
      $2 \\
      ..

EOF

cp .gitignore ../$name/
