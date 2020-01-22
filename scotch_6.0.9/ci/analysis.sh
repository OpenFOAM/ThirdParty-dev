#!/bin/sh

export CPPCHECK_DEFINITIONS="$(grep SCOTCH_GITLAB_SEPARATOR < src/Makefile.inc | sed -e 's#^CFLAGS.*SCOTCH_GITLAB_SEPARATOR##1' | sed -e 's#[ ][^-][^ ]*##g' -e 's#[ ][-][^D][^ ]*##g')"
export CPPCHECK_INCLUDES="-Isrc/scotch -Isrc/misc -Isrc/libscotch -Isrc/esmumps -Isrc/libscotchmetis"

./ci/gitlab-ci-filelist.sh

cppcheck -v --max-configs=1 --language=c ${CPPCHECK_DEFINITIONS_VM:---platform=native} --enable=all --xml --xml-version=2 --suppress=missingIncludeSystem --suppress=varFuncNullUB --suppress=invalidPrintfArgType_sint ${CPPCHECK_DEFINITIONS} ${CPPCHECK_INCLUDES} --file-list=./filelist.txt 2> scotch-cppcheck.xml

rats -w 3 --xml `cat filelist.txt` > scotch-rats.xml

cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.bordeaux.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN
sonar.links.homepage=$CI_PROJECT_URL
sonar.links.scm=$CI_REPOSITORY_URL
sonar.links.ci=https://gitlab.inria.fr/$CI_PROJECT_NAMESPACE/scotch/pipelines
sonar.links.issue=https://gitlab.inria.fr/$CI_PROJECT_NAMESPACE/scotch/issues
sonar.projectKey=tadaam:scotch:gitlab:$CI_PROJECT_NAMESPACE:$CI_COMMIT_REF_NAME
sonar.projectDescription=Package for graph and mesh/hypergraph partitioning, graph clustering, and sparse matrix ordering.
sonar.projectVersion=6.0
sonar.language=c
sonar.sourceEncoding=UTF-8
sonar.sources=include,src/check,src/esmumps,src/libscotch,src/libscotchmetis,src/misc,src/scotch
sonar.c.includeDirectories=$(echo | gcc -E -Wp,-v - 2>&1 | grep "^ " | tr '\n' ',')include,src/scotch,src/misc,src/libscotch,src/esmumps,src/libscotchmetis
sonar.c.errorRecoveryEnabled=true
sonar.c.compiler.charset=UTF-8
sonar.c.compiler.parser=GCC
sonar.c.compiler.regex=^(.*):(\\\d+):\\\d+: warning: (.*)\\\[(.*)\\\]$
sonar.c.compiler.reportPath=scotch-build.log
sonar.c.coverage.reportPath=scotch-coverage.xml
sonar.c.cppcheck.reportPath=scotch-cppcheck.xml
sonar.c.rats.reportPath=scotch-rats.xml
sonar.c.clangsa.reportPath=analyzer_reports/*/*.plist
EOF
