#!/bin/bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

source ci/env_makefile.sh
cd src
make check${JOBCHECK} || fatal
make ptcheck${JOBCHECK} || fatal
make escheck || fatal
find . -print0 | xargs -0 touch
cp libscotch/parser_ll.c libscotch/lex.yy.c
cp libscotch/parser_yy.c libscotch/y.tab.c
cat check/test_scotch_graph_dump.c /tmp/m16x16_b1.c > check/test_scotch_graph_dump2.c
cd ..
gcovr --xml-pretty --exclude-unreachable-branches --print-summary -o ${JOBNAME}.cov --root . || fatal
#if [[ $CI_PIPELINE_SOURCE == "schedule" ]]
#then
#  lcov --capture --directory . --output-file scotch-${JOBNAME}.lcov || fatal
#fi
