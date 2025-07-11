## Copyright 2014-2016,2021-2023 IPB, Universite de Bordeaux, INRIA & CNRS
##
## This file is part of the Scotch software package for static mapping,
## graph partitioning and sparse matrix ordering.
##
## This software is governed by the CeCILL-C license under French law
## and abiding by the rules of distribution of free software. You can
## use, modify and/or redistribute the software under the terms of the
## CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
## URL: "http://www.cecill.info".
##
## As a counterpart to the access to the source code and rights to copy,
## modify and redistribute granted by the license, users are provided
## only with a limited warranty and the software's author, the holder of
## the economic rights, and the successive licensors have only limited
## liability.
##
## In this respect, the user's attention is drawn to the risks associated
## with loading, using, modifying and/or developing or reproducing the
## software by the user in light of its specific status of free software,
## that may mean that it is complicated to manipulate, and that also
## therefore means that it is reserved for developers and experienced
## professionals having in-depth computer knowledge. Users are therefore
## encouraged to load and test the software's suitability as regards
## their requirements in conditions enabling the security of their
## systems and/or data to be ensured and, more generally, to use and
## operate it in the same conditions as regards security.
##
## The fact that you are presently reading this means that you have had
## knowledge of the CeCILL-C license and that you accept its terms.
##
############################################################
##                                                        ##
##   AUTHORS    : Cedric LACHAT                           ##
##                Amaury JACQUES                          ##
##                Florent PRUVOST                         ##
##                Marc FUENTES                            ##
##                Francois PELLEGRINI                     ##
##                                                        ##
##   FUNCTION   : Secondary configuration file for CMake  ##
##                                                        ##
##   DATES      : # Version 6.0  : from : 01 sep 2014     ##
##                                 to     01 sep 2021     ##
##                # Version 7.0  : from : 01 sep 2021     ##
##                                 to     21 apr 2022     ##
##                # Version 7.1  : from : 11 dec 2022     ##
##                                 to     01 nov 2023     ##
##                                                        ##
############################################################

set_source_files_properties(${GENERATED_INCLUDE_DIR}/scotch.h PROPERTIES GENERATED 1)

################
#  libEsMUMPS  #
################

# esmumps.h
add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/esmumps.h
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/library.h ${GENERATED_INCLUDE_DIR}/esmumps.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/library.h
  COMMENT "Copy esmumps.h to include directory")

add_custom_target(esmumps_h
  DEPENDS "${GENERATED_INCLUDE_DIR}/esmumps.h")

# libEsMUMPS
set(esmumps_src
  dof.c
  esmumps.c
  graph_graph.c
  library_esmumps.c
  library_esmumps_f.c
  library_esmumps_strats.c
  order.c
  order_check.c
  order_io.c
  order_scotch_graph.c
  symbol.c
  symbol_check.c
  symbol_cost.c
  symbol_fax_graph.c
  symbol_fax.c)

set(esmumps_headers
  module.h
  library.h
  ${GENERATED_INCLUDE_DIR}/scotch.h
  ${LIBSCOTCH_INCLUDE_DIR}/common.h
  graph.h
  dof.h
  symbol.h
  order.h
  fax.h
  ${GENERATED_INCLUDE_DIR}/esmumps.h)

add_library(esmumps ${esmumps_src} ${esmumps_headers})
add_dependencies(esmumps esmumps_h)
set_target_properties(esmumps PROPERTIES PUBLIC_HEADER "${GENERATED_INCLUDE_DIR}/esmumps.h")
set_target_properties(esmumps PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
  SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})

if(APPLE)
  set_target_properties(esmumps PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

target_include_directories(esmumps PRIVATE
  $<BUILD_INTERFACE:${LIBSCOTCH_INCLUDE_DIR}>
  $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(esmumps PRIVATE scotch)
if(NOT "${LIBSCOTCHERR}" STREQUAL "")
  target_link_libraries(esmumps PRIVATE ${LIBSCOTCHERR})
endif()

add_dependencies(esmumps scotch_h)

add_executable(test_esmumps test_esmumps.c)
set_target_properties(test_esmumps PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(test_esmumps PRIVATE ${GENERATED_INCLUDE_DIR}/)
target_link_libraries(test_esmumps PRIVATE esmumps scotch scotcherrexit)
add_dependencies(test_esmumps scotch_h)

add_executable(test_fax test_fax.c)
set_target_properties(test_fax PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(test_fax PRIVATE ${GENERATED_INCLUDE_DIR}/)
target_link_libraries(test_fax PRIVATE esmumps scotch scotcherrexit)
add_dependencies(test_fax scotch_h)

# libEsMUMPS targets install
install(EXPORT esmumpsTargets
  COMPONENT libscotch
  FILE esmumpsTargets.cmake
  NAMESPACE SCOTCH::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scotch)

install(TARGETS esmumps
  COMPONENT libscotch
  EXPORT esmumpsTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Dummy PTEsMUMPS library file
if(BUILD_PTSCOTCH)
  add_library(ptesmumps ${esmumps_src} ${esmumps_headers})

  if(APPLE)
    set_target_properties(ptesmumps PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif()

  set_target_properties(ptesmumps PROPERTIES PUBLIC_HEADER "${GENERATED_INCLUDE_DIR}/esmumps.h")

  set_target_properties(ptesmumps PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
    SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})

  target_include_directories(ptesmumps PRIVATE
    $<BUILD_INTERFACE:${LIBSCOTCH_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  target_link_libraries(ptesmumps PRIVATE scotch)
  if(NOT "${LIBSCOTCHERR}" STREQUAL "")
    target_link_libraries(ptesmumps PRIVATE ${LIBSCOTCHERR})
  endif()
  add_dependencies(ptesmumps scotch_h)

  # libptesmumps targets install
  install(EXPORT ptesmumpsTargets
    COMPONENT libptscotch
    FILE ptesmumpsTargets.cmake
    NAMESPACE SCOTCH::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scotch)

  install(TARGETS ptesmumps
    COMPONENT libptscotch
    EXPORT ptesmumpsTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()
