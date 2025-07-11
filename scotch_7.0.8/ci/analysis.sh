#!/bin/sh

#export LCOVFILES=""
#for filename in $(ls -1 scotch-*.lcov); do export LCOVFILES="$LCOVFILES -a $filename"; done
#lcov $LCOVFILES -o scotch.lcov
#lcov_cobertura scotch.lcov --output scotch-coverage.xml
gcovr --xml-pretty --exclude-unreachable-branches --print-summary -o coverage.xml --root .

export CPPCHECK_DEFINITIONS="$(grep SCOTCH_GITLAB_SEPARATOR < src/Makefile.inc | sed -e 's#^CFLAGS.*SCOTCH_GITLAB_SEPARATOR##1' | sed -e 's#[ ][^-][^ ]*##g' -e 's#[ ][-][^D][^ ]*##g')"
export CPPCHECK_INCLUDES="-Isrc/scotch -Isrc/misc -Isrc/libscotch -Isrc/esmumps -Isrc/libscotchmetis"

./ci/gitlab-ci-filelist.sh

cppcheck -v --max-configs=1 --language=c ${CPPCHECK_DEFINITIONS_VM:---platform=native} --enable=all --xml --xml-version=2 --suppress=missingIncludeSystem --suppress=varFuncNullUB --suppress=invalidPrintfArgType_sint ${CPPCHECK_DEFINITIONS} ${CPPCHECK_INCLUDES} --file-list=./filelist.txt 2> scotch-cppcheck.xml

rats -w 3 --xml `cat filelist.txt` > scotch-rats.xml

# sonar.projectKey=tadaam:scotch:gitlab:$CI_PROJECT_NAMESPACE:$CI_COMMIT_REF_NAME
# sonar.projectKey=tadaam:scotch:gitlab:$CI_PROJECT_NAMESPACE
# sonar.branch.name=$CI_COMMIT_REF_NAME
# sonar.exclusions=build/**,cmake_modules/**,doc/**,tools/**,**/*.xml
# "sonar.lang.patterns.c++=**/*.cpp,**/*.hpp" : needed so that C files are not considered as C++ files
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=https://gitlab.inria.fr/$CI_PROJECT_NAMESPACE/scotch/pipelines
sonar.links.issue=https://gitlab.inria.fr/$CI_PROJECT_NAMESPACE/scotch/issues

sonar.projectKey=tadaam:scotch:gitlab:$CI_PROJECT_NAMESPACE
sonar.projectDescription=Package for graph and mesh/hypergraph partitioning, graph clustering, and sparse matrix ordering.
sonar.projectVersion=7.1
sonar.branch.name=$CI_COMMIT_REF_NAME

sonar.sourceEncoding=UTF-8
sonar.sources=include,src/check,src/esmumps,src/libscotch,src/libscotchmetis,src/misc,src/scotch
sonar.inclusions=`cat filelist.txt | sed ':a;N;$!ba;s/\n/, /g'`
sonar.exclusions=src/libscotch/dummysizes.c,src/libscotch/parser_ll.c,src/libscotch/parser_yy.c,src/libscotch/last_resort/**
sonar.cxx.file.suffixes=.h,.c
sonar.cxx.errorRecoveryEnabled=true
sonar.cxx.gcc.encoding=UTF-8
sonar.cxx.gcc.regex=(?<file>.*):(?<line>[0-9]+):[0-9]+:\\\x20warning:\\\x20(?<message>.*)\\\x20\\\[(?<id>.*)\\\]
sonar.cxx.gcc.reportPaths=scotch-build*.log
sonar.cxx.clangsa.reportPaths=analyzer_reports/*/*.plist
sonar.cxx.cppcheck.reportPaths=scotch-cppcheck.xml
sonar.cxx.rats.reportPaths=scotch-rats.xml
sonar.cxx.cobertura.reportPaths=*.cov
EOF
