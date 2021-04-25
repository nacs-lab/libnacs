#!/bin/bash

mkdir build
cd build
if [ -n "${LLVM_VER}" ]; then
    cmake .. -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DENABLE_LLVM=On -DENABLE_SIMD=On
else
    cmake .. -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DENABLE_SIMD=On
fi
make -j $(nproc)
