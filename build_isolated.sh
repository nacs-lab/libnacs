#!/bin/bash
# Isolated build: own build dir, install prefix is a local folder; never run `make install`.
set -e
cd /c/projects/libnacs-coilcheck
mkdir -p build
cd build
cmake .. -G "MSYS Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/c/projects/libnacs-coilcheck/install-local-only \
  -DINSTALL_PKGCONFIG_DIR=/c/projects/libnacs-coilcheck/install-local-only/lib/pkgconfig \
  -DENABLE_TESTING=OFF -DENABLE_KNACS=OFF -DENABLE_LLVM=On -DENABLE_SIMD=On \
  -DSTATIC_LLVM=On -DSEQ_IN_UTILS=OFF -DENABLE_INSTALL_NONLIB=OFF \
  -DLLVM_DIR=C:/msys64/mingw64/lib/cmake/llvm
make -j8 "$@"
