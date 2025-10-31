#!/bin/bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

source ci/env_cmake.sh
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON .. || fatal
make -j5 || fatal
