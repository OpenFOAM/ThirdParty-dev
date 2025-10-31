#!/bin/bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

source ci/env_cmake.sh
cd build
ctest -T test -O ctest.log --output-junit ../report.xml --timeout 20 || fatal
cd ..
