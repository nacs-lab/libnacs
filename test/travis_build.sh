#!/bin/bash

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DENABLE_LLVM=On -DENABLE_SIMD=On
make -j $(nproc)
