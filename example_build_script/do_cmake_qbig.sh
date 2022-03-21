rm -rv CMakeCache.txt CMakeFiles

#source load_modules.sh
module purge
module load fosscuda/2020b CMake/3.18.4-GCCcore-10.2.0
module load libarchive/3.4.3-GCCcore-10.2.0
module load GSL/2.7-GCC-11.2.0 
module load Eigen/3.3.9-GCCcore-11.2.0
module load Python/3.9.6-GCCcore-11.2.0

src_dir=..
#CXXFLAGS="-O3 -mtune=native -march=native -g" \
cmake   ${src_dir}   \
  -DWITH_PYTHON=ON  \
  -DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_PREFIX_PATH=/opt/eb/software/GSL/2.7-GCC-11.2.0 \
  -DCMAKE_CXX_FLAGS="-fopenmp   -pedantic  -g -O3   -lm -lgmp "

make -j8

  
 # -D CMAKE_FIND_DEBUG_MODE=ON



