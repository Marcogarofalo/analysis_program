rm -rv CMakeCache.txt
rm -rv CMakeFiles

#source load_modules.sh


src_dir=..
#CXXFLAGS="-O3 -DNDEBUG" \

cmake   ${src_dir} \
	-DCMAKE_BUILD_TYPE=RELEASE \
	-DCMAKE_INSTALL_PREFIX=$(pwd)/install_dir \
	-DBUILD_SHARED_LIBS=ON \
	-DWITH_PYTHON=ON \
	-DWITH_ARB=ON \
	-DARB_ROOT=/home/garofalo/programs/flint/install_dir/ \
	-DPHI4_PROJECT=ON \
	-DGM2_PROJECT=ON \
	-DBSM_PROJECT=OFF \
	-Dclover_PROJECT=ON \
	-DALL_PROJECT=ON \
        -DCMAKE_CXX_FLAGS="  -g -O2   -lm" \
	-DCMAKE_VERBOSE_MAKEFILE=OFF \
	-DCMAKE_CXX_COMPILER=/usr/bin/g++ \
	
	#-GNinja \
	#-G "Unix Makefiles"
 # -D CMAKE_FIND_DEBUG_MODE=ON
 # -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
#-DCMAKE_PREFIX_PATH="/usr/lib/x86_64-linux-gnu" \
# -lgmp 
# -DCMAKE_CXX_FLAGS="-fopenmp   -g -O2   -lm  -lmpir -I/home/garofalo/programs/mpir/install/include -L/home/garofalo/programs/mpir/install/lib" \


