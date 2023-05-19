rm -rv CMakeCache.txt CMakeFiles

source load_modules_marconi.sh

GSL="${GSL_HOME}/lib/pkgconfig/"
src_dir=..
cmake   ${src_dir}  -DCMAKE_MODULE_PATH=$GSL \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_EXE_LINKER_FLAGS="-L $GSL_LIB " \
    -DCMAKE_CXX_FLAGS="-fopenmp   -pedantic  -g -O2    " \
    -DCMAKE_C_FLAGS="-fopenmp   -pedantic  -g -O2   -I  $GSL_INCLUDE " \

 
 
