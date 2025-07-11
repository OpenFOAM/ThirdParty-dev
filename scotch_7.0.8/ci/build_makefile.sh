#!/bin/bash

fatal() {
    echo "$0: error occurred, exit"
    exit 1
}

source ci/env_makefile.sh
cp ci/Makefile.inc src/Makefile.inc
cd src
sed -i -e "s#^\(CFLAGS.*\)#\1 ${JOBBUILDFLAGS}#1" Makefile.inc
find . -name Makefile | xargs sed -i -e "s#-c \$(<)#-c \$\(shell pwd\)/\$\(<\)#g"

if [[ $CI_PIPELINE_SOURCE == "schedule" ]]
then
  SCAN="scan-build -plist --intercept-first --analyze-headers -o ../analyzer_reports "
else
  SCAN=""
fi
eval '${SCAN}make scotch ptscotch esmumps ptesmumps 2>&1 | tee ../scotch-build-${JOBNAME}.log' || fatal
