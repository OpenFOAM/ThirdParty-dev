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
##                                                        ##
##   FUNCTION   : Secondary configuration file for CMake  ##
##                                                        ##
##   DATES      : # Version 6.0  : from : 01 sep 2014     ##
##                                 to     01 sep 2021     ##
##                # Version 7.0  : from : 01 sep 2021     ##
##                                 to     01 nov 2023     ##
##                                                        ##
############################################################

####################
#  libScotchMeTiS  #
####################

set_source_files_properties(${GENERATED_INCLUDE_DIR}/scotch.h PROPERTIES GENERATED 1)
set_source_files_properties(${GENERATED_INCLUDE_DIR}/metis.h PROPERTIES GENERATED 1)

# metis.h
add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/metis.h
  COMMAND $<TARGET_FILE:dummysizes> ${CMAKE_CURRENT_SOURCE_DIR}/library_metis.h ${GENERATED_INCLUDE_DIR}/metis.h DEPENDS dummysizes
  DEPENDS library_metis.h
  COMMENT "Generating metis.h")

add_custom_target(metis_h
  DEPENDS "${GENERATED_INCLUDE_DIR}/metis.h")

# metisf.h
add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/metisf.h
  COMMAND $<TARGET_FILE:dummysizes> ${CMAKE_CURRENT_SOURCE_DIR}/library_metis_f.h ${GENERATED_INCLUDE_DIR}/metisf.h DEPENDS dummysizes
  DEPENDS library_metis_f.h
  COMMENT "Generating metisf.h")

add_custom_target(metisf_h
  DEPENDS "${GENERATED_INCLUDE_DIR}/metisf.h")

# libScotchMeTiS
foreach(version 3 5)
  add_library(scotchmetisv${version}
    metis_graph_order.c
    metis_graph_order_f.c
    metis_graph_part.c
    metis_graph_part_f.c
    metis_graph_dual.c
    metis_graph_dual_f.c
    metis_graph_part_dual.c
    metis_graph_part_dual_f.c
    metis_options.c
    metis_options_f.c
    ${LIBSCOTCH_INCLUDE_DIR}/module.h
    ${LIBSCOTCH_INCLUDE_DIR}/common.h
    ${GENERATED_INCLUDE_DIR}/scotch.h
    ${GENERATED_INCLUDE_DIR}/metis.h
    ${GENERATED_INCLUDE_DIR}/metisf.h)

  set_target_properties(scotchmetisv${version} PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
    SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})
  target_compile_definitions(scotchmetisv${version} PUBLIC SCOTCH_METIS_VERSION=${version})
  if(SCOTCH_METIS_PREFIX)
    target_compile_definitions(scotchmetisv${version} PUBLIC SCOTCH_METIS_PREFIX)
  endif()

  target_include_directories(scotchmetisv${version} PRIVATE
    $<BUILD_INTERFACE:${LIBSCOTCH_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  target_link_libraries(scotchmetisv${version} PRIVATE scotch)
  if(NOT "${LIBSCOTCHERR}" STREQUAL "")
    target_link_libraries(scotchmetisv${version} PRIVATE ${LIBSCOTCHERR})
  endif()

  add_dependencies(scotchmetisv${version} scotch_h metisf_h)
  if(INSTALL_METIS_HEADERS)
    set_target_properties(scotchmetisv${version} PROPERTIES PUBLIC_HEADER "${GENERATED_INCLUDE_DIR}/metis.h;${GENERATED_INCLUDE_DIR}/metisf.h")
  endif()

  if(APPLE)
    set_target_properties(scotchmetisv${version} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif()

  # libScotchMeTiS targets install
  install(EXPORT scotchmetisTargets
    COMPONENT libscotch
    FILE scotchmetisTargets.cmake
    NAMESPACE SCOTCH::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scotch)
  if(INSTALL_METIS_HEADERS)
    install(TARGETS scotchmetisv${version}
      COMPONENT libscotch
      EXPORT scotchmetisTargets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  else()
    install(TARGETS scotchmetisv${version}
      COMPONENT libscotch
      EXPORT scotchmetisTargets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endforeach(version 3 5)

#########################
#  libPTScotchParMeTiS  #
#########################

if(BUILD_PTSCOTCH)
  set_source_files_properties(${GENERATED_INCLUDE_DIR}/ptscotch.h PROPERTIES GENERATED 1)

  # parmetis.h
  add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/parmetis.h
    COMMAND $<TARGET_FILE:dummysizes> ${CMAKE_CURRENT_SOURCE_DIR}/library_parmetis.h ${GENERATED_INCLUDE_DIR}/parmetis.h DEPENDS dummysizes
    DEPENDS library_parmetis.h
    COMMENT "Generating parmetis.h")

  add_custom_target(parmetis_h
    DEPENDS "${GENERATED_INCLUDE_DIR}/parmetis.h")

  # libPTScotchParMeTiS
  foreach(version 3)
    add_library(ptscotchparmetisv${version}
      parmetis_dgraph_order.c
      parmetis_dgraph_order_f.c
      parmetis_dgraph_part.c
      parmetis_dgraph_part_f.c
      ${LIBSCOTCH_INCLUDE_DIR}/module.h
      ${LIBSCOTCH_INCLUDE_DIR}/common.h
      ${GENERATED_INCLUDE_DIR}/ptscotch.h
      ${GENERATED_INCLUDE_DIR}/parmetis.h)

    set_target_properties(ptscotchparmetisv${version} PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
      SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})
    target_compile_definitions(ptscotchparmetisv${version} PUBLIC SCOTCH_PTSCOTCH SCOTCH_METIS_VERSION=${version})
    if(SCOTCH_METIS_PREFIX)
      target_compile_definitions(ptscotchparmetisv${version} PUBLIC SCOTCH_METIS_PREFIX)
    endif()
    target_include_directories(ptscotchparmetisv${version} PRIVATE
      $<BUILD_INTERFACE:${LIBSCOTCH_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

    target_link_libraries(ptscotchparmetisv${version} PRIVATE scotch)
    target_link_libraries(ptscotchparmetisv${version} PRIVATE ptscotch)
    target_link_libraries(ptscotchparmetisv${version} PUBLIC MPI::MPI_C)
    if(NOT "${LIBPTSCOTCHERR}" STREQUAL "")
      target_link_libraries(ptscotchparmetisv${version} PRIVATE ${LIBPTSCOTCHERR})
    endif()

    add_dependencies(ptscotchparmetisv${version} ptscotch_h)
    if(INSTALL_METIS_HEADERS)
      set_target_properties(ptscotchparmetisv${version} PROPERTIES PUBLIC_HEADER "${GENERATED_INCLUDE_DIR}/parmetis.h")
    endif()

    if(APPLE)
      set_target_properties(ptscotchparmetisv${version} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()

    # libPTScotchParMeTiS targets install
    install(EXPORT ptscotchparmetisTargets
      COMPONENT libptscotch
      FILE ptscotchparmetisTargets.cmake
      NAMESPACE SCOTCH::
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scotch)

    if(INSTALL_METIS_HEADERS)
      install(TARGETS ptscotchparmetisv${version}
        COMPONENT libptscotch
        EXPORT ptscotchparmetisTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    else()
      install(TARGETS ptscotchparmetisv${version}
        COMPONENT libptscotch
        EXPORT ptscotchparmetisTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
  endforeach(version 3)
endif()
