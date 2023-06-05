rm -rv CMakeCache.txt CMakeFiles

source load_modules_marconi.sh

GSL="${GSL_HOME}/lib/pkgconfig/"
src_dir=..
cmake   ${src_dir}  -DCMAKE_MODULE_PATH=../eigen \
    -DCMAKE_CXX_COMPILER=g++ \
    -DEigen3_DIR="/m100_work/INF23_lqcd123_0/mgarofal/eigen/install_dir"\
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_EXE_LINKER_FLAGS="-L $GSL_LIB " \
    -DCMAKE_CXX_FLAGS="-fopenmp   -pedantic  -g -O2    " \
    -DCMAKE_C_FLAGS="-fopenmp   -pedantic  -g -O2   -I  $GSL_INCLUDE " \

 
 
