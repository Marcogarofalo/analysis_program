#  Eigen is required

````
git clone https://gitlab.com/libeigen/eigen.git
````

we also need the rzeta package https://github.com/HISKP-LQCD/rzeta as extra 
dependency
```
git submodule update --init --recursive
```

create a build directory and compile

   ````
  mkdir build
  cp example_build_script/do_cmake.sh  build
  cd build
  bash do_cmake.sh
   ````
