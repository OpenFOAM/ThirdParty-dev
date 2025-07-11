## Copyright 2014-2016,2021,2022,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
##   AUTHORS    : Marc FUENTES                            ##
##                Florent PRUVOST                         ##
##                Cedric LACHAT                           ##
##                Amaury JACQUES                          ##
##                                                        ##
##   FUNCTION   : Secondary configuration file for CMake  ##
##                                                        ##
##   DATES      : # Version 6.0  : from : 01 sep 2014     ##
##                                 to     01 sep 2021     ##
##                # Version 7.0  : from : 01 sep 2021     ##
##                                 to     26 jan 2023     ##
##                                                        ##
############################################################

# Set directory aliases
set(dat ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(tgt ${PROJECT_SOURCE_DIR}/tgt)

# Macro to add quickly a Scotch test executable
function(add_test_scotch)
  add_executable(${ARGV0} ${ARGV0}.c)
  set_target_properties(${ARGV0} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(${ARGV0} PRIVATE scotch scotcherr)
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    target_compile_definitions(${ARGV0} PRIVATE SCOTCH_DEBUG_ALL)
  endif()
  if(${ARGC} GREATER_EQUAL 2)
    add_test(NAME ${ARGV0} COMMAND ${ARGV0})
  endif()
endfunction(add_test_scotch)

# Macro to suffix a fortran test file
function(suffix_fortran_file src)
  get_filename_component(file_we ${src} NAME_WLE)
  add_custom_command(OUTPUT
    ${CMAKE_CURRENT_BINARY_DIR}/${file_we}.tmp
    COMMAND
    $<TARGET_FILE:dummysizes> "-s${SCOTCH_NAME_SUFFIX}"
    ${CMAKE_CURRENT_SOURCE_DIR}/${file_we}.in
    ${CMAKE_CURRENT_BINARY_DIR}/${file_we}.tmp
    DEPENDS
    dummysizes
    DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/${file}
    COMMENT "Generate ${file}.f90")

  add_custom_target(${file_we}_tmp ALL
    DEPENDS
    "${CMAKE_CURRENT_BINARY_DIR}/${file_we}.tmp")

  file(WRITE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/cmake_sed_${file_we}.cmake
"
file(READ ${file_we}.tmp FILE_CONTENTS)
string(REGEX REPLACE \"SCOTCHMETISNAMEFU *\\\\(([^()]*)\\\\) *\" "${SCOTCH_METIS_PREFIXFU}\\\\1" FILE_CONTENTS \${FILE_CONTENTS})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${file_we} \${FILE_CONTENTS})
"
)

  add_custom_command(
    OUTPUT
      ${CMAKE_CURRENT_BINARY_DIR}/${file_we}
    COMMAND
      ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/cmake_sed_${file_we}.cmake
    DEPENDS
      ${file_we}_tmp
  )
  add_custom_target(${file_we}_src ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${file_we})
endfunction(suffix_fortran_file)


# Bash needed for running some tests
include(FindUnixCommands)

#######################################
#  Tests of the sequential libraries  #
#######################################

# test_common_file_compress
add_test_scotch(test_common_file_compress)
set(_map_bz2 "BZIP2")
set(_map_gz "ZLIB")
set(_map_lzma "LIBLZMA")
foreach(ext "bz2" "gz" "lzma")
  if(${_map_${ext}}_FOUND)
    set(cur ${CMAKE_CURRENT_BINARY_DIR})
    set(src_cur ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME test_common_file_compress_${ext} COMMAND ${BASH} -c
      "'${CMAKE_COMMAND}' -E copy ${src_cur}/data/bump_b1.grf ${cur}/bump_b1.grf && \
      '$<TARGET_FILE:test_common_file_compress>' ${cur}/bump_b1.grf ${cur}/bump_b1.grf.${ext} && \
      '$<TARGET_FILE:test_common_file_compress>' ${cur}/bump_b1.grf.${ext} ${cur}/bump_b1.grf && \
      '${CMAKE_COMMAND}' -E compare_files ${cur}/bump_b1.grf ${src_cur}/data/bump_b1.grf")
  endif(${_map_${ext}}_FOUND)
endforeach()

# test_common_random
add_test_scotch(test_common_random)
add_test(NAME test_common_random_0 COMMAND test_common_random ${CMAKE_CURRENT_BINARY_DIR}/rand.dat 0)
add_test(NAME test_common_random_1 COMMAND test_common_random ${CMAKE_CURRENT_BINARY_DIR}/rand.dat 1)

# test_common_thread
add_test_scotch(test_common_thread 1)

# test_fibo
add_test_scotch(test_fibo 1)

# Test EsMumps
if(BUILD_LIBESMUMPS AND NOT USE_SUFFIX)
  add_executable(test_libesmumps test_libesmumps.c)
  add_dependencies(test_libesmumps esmumps_h)
  set_target_properties(test_libesmumps PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  target_include_directories(test_libesmumps PRIVATE ${GENERATED_INCLUDE_DIR})
  target_link_libraries(test_libesmumps scotch esmumps scotcherr)
  add_test(NAME test_libesmumps_1 COMMAND $<TARGET_FILE:test_libesmumps> ${dat}/bump.grf)
  add_test(NAME test_libesmumps_2 COMMAND $<TARGET_FILE:test_libesmumps> ${dat}/bump_imbal_32.grf)
endif(BUILD_LIBESMUMPS AND NOT USE_SUFFIX)

# Test libScotchMeTiS
if(BUILD_LIBSCOTCHMETIS)
  # manage prefix
  if(NOT SCOTCH_METIS_PREFIXFU)
    if(NOT SCOTCH_METIS_PREFIX)
      set(SCOTCH_METIS_PREFIXFU "")
    else()
      set(SCOTCH_METIS_PREFIXFU "SCOTCHF")
    endif(NOT SCOTCH_METIS_PREFIX)
  endif(NOT SCOTCH_METIS_PREFIXFU)

  suffix_fortran_file(test_libmetis_dual_f.f90.in)
  foreach(version 3 5)
    add_executable(test_libmetisv${version} test_libmetis.c)
    set_target_properties(test_libmetisv${version} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    target_compile_definitions(test_libmetisv${version} PRIVATE SCOTCH_METIS_VERSION=${version})
    target_include_directories(test_libmetisv${version} PRIVATE ${GENERATED_INCLUDE_DIR})
    target_link_libraries(test_libmetisv${version} scotchmetisv${version} scotch scotcherr)

    add_test(NAME test_libmetisv${version}_1 COMMAND $<TARGET_FILE:test_libmetisv${version}> ${dat}/bump.grf)
    add_test(NAME test_libmetisv${version}_2 COMMAND $<TARGET_FILE:test_libmetisv${version}> ${dat}/bump_b1.grf)
    add_test(NAME test_libmetisv${version}_3 COMMAND $<TARGET_FILE:test_libmetisv${version}> ${dat}/bump_imbal_32.grf)

    add_executable(test_libmetis_dualv${version} test_libmetis_dual.c)
    target_compile_definitions(test_libmetis_dualv${version} PRIVATE SCOTCH_METIS_VERSION=${version})
    target_include_directories(test_libmetis_dualv${version} PRIVATE ${GENERATED_INCLUDE_DIR})
    target_link_libraries(test_libmetis_dualv${version} scotchmetisv${version} scotch scotcherr)
    add_test(NAME test_libmetis_dualv${version} COMMAND test_libmetis_dualv${version})
    if(CMAKE_Fortran_COMPILER)
      add_executable(test_libmetis_dual_f_v${version} ${CMAKE_CURRENT_BINARY_DIR}/test_libmetis_dual_f.f90)
      target_compile_definitions(test_libmetis_dual_f_v${version} PRIVATE SCOTCH_METIS_VERSION=${version})
      target_include_directories(test_libmetis_dual_f_v${version} PRIVATE ${GENERATED_INCLUDE_DIR})
      target_link_libraries(test_libmetis_dual_f_v${version} scotch scotchmetisv${version} scotcherr)
      add_test(NAME test_libmetis_dual_f_v${version} COMMAND test_libmetis_dual_f_v${version})
    endif()
  endforeach(version 3 5)
endif(BUILD_LIBSCOTCHMETIS)

# test_scotch_arch
add_test_scotch(test_scotch_arch)
add_test(NAME test_scotch_arch COMMAND test_scotch_arch  arch.tgt)

# test_scotch_arch_deco
add_test_scotch(test_scotch_arch_deco)
set(cur_src ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME test_scotch_arch_deco COMMAND test_scotch_arch_deco ${cur_src}/data/m4x4.grf arch2.tgt)

# test_scotch_context
add_test_scotch(test_scotch_context)
add_test(NAME test_scotch_context COMMAND test_scotch_context ${cur_src}/data/bump_b1.grf)
if(Threads_FOUND)
  target_link_libraries(test_scotch_context PRIVATE Threads::Threads ${THREADS_PTHREADS_WIN32_LIBRARY})
endif(Threads_FOUND)

# test_graph_coarsen
add_test_scotch(test_scotch_graph_coarsen)
add_test(NAME test_scotch_graph_coarsen_1 COMMAND test_scotch_graph_coarsen ${cur_src}/data/bump.grf)
add_test(NAME test_scotch_graph_coarsen_2 COMMAND test_scotch_graph_coarsen ${cur_src}/data/bump_b100000.grf)

# test_graph_color
add_test_scotch(test_scotch_graph_color)
add_test(NAME test_scotch_graph_color_1 COMMAND test_scotch_graph_color ${cur_src}/data/bump.grf)
add_test(NAME test_scotch_graph_color_2 COMMAND test_scotch_graph_color ${cur_src}/data/bump_b100000.grf)

# test_graph_diam
add_test_scotch(test_scotch_graph_diam)
add_test(NAME test_scotch_graph_diam_1 COMMAND test_scotch_graph_diam ${cur_src}/data/bump.grf)
add_test(NAME test_scotch_graph_diam_2 COMMAND test_scotch_graph_diam ${cur_src}/data/m16x16_b1.grf)
add_test(NAME test_scotch_graph_diam_3 COMMAND test_scotch_graph_diam ${cur_src}/data/p2-p2.grf)
add_test(NAME test_scotch_graph_diam_4 COMMAND test_scotch_graph_diam ${cur_src}/data/bump_b100000.grf)

# test_graph_dump
#TODO

# test_graph_induce
add_test_scotch(test_scotch_graph_induce)
add_test(NAME test_scotch_graph_induce_1 COMMAND test_scotch_graph_induce ${cur_src}/data/bump.grf)
add_test(NAME test_scotch_graph_induce_2 COMMAND test_scotch_graph_induce ${cur_src}/data/bump_b100000.grf)

# test_scotch_graph_map
add_test_scotch(test_scotch_graph_map)
add_test(NAME test_scotch_graph_1 COMMAND test_scotch_graph_map ${cur_src}/data/m4x4.grf)
add_test(NAME test_scotch_graph_2 COMMAND test_scotch_graph_map ${cur_src}/data/m4x4_b1.grf)
add_test(NAME test_scotch_graph_3 COMMAND test_scotch_graph_map ${cur_src}/data/m16x16.grf)
add_test(NAME test_scotch_graph_4 COMMAND test_scotch_graph_map ${cur_src}/data/m16x16_b1.grf)

# test_scotch_graph_map_copy
add_test_scotch(test_scotch_graph_map_copy)
add_test(NAME test_scotch_graph_map_copy_1 COMMAND test_scotch_graph_map_copy ${cur_src}/data/bump.grf)
add_test(NAME test_scotch_graph_map_copy_2 COMMAND test_scotch_graph_map_copy ${cur_src}/data/bump_b100000.grf)

# test_scotch_graph_order
add_test_scotch(test_scotch_graph_order)
add_test(NAME test_scotch_graph_order_1 COMMAND test_scotch_graph_order ${cur_src}/data/bump.grf)
add_test(NAME test_scotch_graph_order_2 COMMAND test_scotch_graph_order ${cur_src}/data/bump_b100000.grf)

# test_scotch_graph_part_ovl
add_test_scotch(test_scotch_graph_part_ovl)
add_test(NAME test_scotch_graph_part_ovl_1 COMMAND test_scotch_graph_part_ovl 4 ${cur_src}/data/m16x16.grf ${dev_null})
add_test(NAME test_scotch_graph_part_ovl_2 COMMAND test_scotch_graph_part_ovl 4 ${cur_src}/data/m16x16_b1.grf ${dev_null})

# test_scotch_mesh_graph
add_test_scotch(test_scotch_mesh_graph)
add_test(NAME test_scotch_mesh_graph COMMAND test_scotch_mesh_graph  ${cur_src}/data/ship001.msh)

# test_strat_seq
add_test_scotch(test_strat_seq 1)

#####################################
#  Tests of the parallel libraries  #
#####################################

if(BUILD_PTSCOTCH)
  # Macro to add quickly a PT-Scotch test executable
  function(add_exec_ptscotch)
    add_executable(${ARGV0} ${ARGV0}.c)
    set_target_properties(${ARGV0} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(${ARGV0} ptscotch ptscotcherr)
  endfunction(add_exec_ptscotch)

  function(add_test_ptscotch)
    set(args ${ARGV})
    list(REMOVE_AT args 0) # in place of POP_FRONT which is available only for cmake >3.15
    string(REGEX REPLACE "^(.+)/data/(.+)\\.grf" "\\2" arg_ext ${args})
    add_test(NAME ${ARGV0}_${arg_ext} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:${ARGV0}> ${args})
  endfunction(add_test_ptscotch)

  # test_strat_par
  add_exec_ptscotch(test_strat_par)
  add_test(NAME test_strat_par COMMAND test_strat_par)

  # test_scotch_dgraph_band
  add_exec_ptscotch(test_scotch_dgraph_band)
  add_test_ptscotch(test_scotch_dgraph_band ${CMAKE_CURRENT_SOURCE_DIR}/data/bump.grf test_scotch_dgraph_band.map)
  add_test_ptscotch(test_scotch_dgraph_band ${CMAKE_CURRENT_SOURCE_DIR}/data/bump_b100000.grf test_scotch_dgraph_band.map)

  # test_scotch_dgraph_coarsen
  add_exec_ptscotch(test_scotch_dgraph_coarsen)
  add_test_ptscotch(test_scotch_dgraph_coarsen ${CMAKE_CURRENT_SOURCE_DIR}/data/bump.grf)
  add_test_ptscotch(test_scotch_dgraph_coarsen ${CMAKE_CURRENT_SOURCE_DIR}/data/bump_b100000.grf)
  add_test_ptscotch(test_scotch_dgraph_coarsen ${CMAKE_CURRENT_SOURCE_DIR}/data/m4x4_b1.grf)

  # test_scotch_dgraph_check
  add_exec_ptscotch(test_scotch_dgraph_check)
  add_test_ptscotch(test_scotch_dgraph_check ${CMAKE_CURRENT_SOURCE_DIR}/data/bump.grf)
  add_test_ptscotch(test_scotch_dgraph_check ${CMAKE_CURRENT_SOURCE_DIR}/data/bump_b100000.grf)

  # test_scotch_dgraph_grow
  add_exec_ptscotch(test_scotch_dgraph_grow)
  add_test_ptscotch(test_scotch_dgraph_grow ${CMAKE_CURRENT_SOURCE_DIR}/data/bump.grf test_scotch_dgraph_grow.map)
  add_test_ptscotch(test_scotch_dgraph_grow ${CMAKE_CURRENT_SOURCE_DIR}/data/bump_b100000.grf test_scotch_dgraph_grow.map)

  # test_scotch_dgraph_induce
  add_exec_ptscotch(test_scotch_dgraph_induce)
  add_test_ptscotch(test_scotch_dgraph_induce ${CMAKE_CURRENT_SOURCE_DIR}/data/bump.grf)
  add_test_ptscotch(test_scotch_dgraph_induce ${CMAKE_CURRENT_SOURCE_DIR}/data/bump_b100000.grf)

  # test_scotch_dgraph_redist
  add_exec_ptscotch(test_scotch_dgraph_redist)
  add_test_ptscotch(test_scotch_dgraph_redist ${CMAKE_CURRENT_SOURCE_DIR}/data/bump.grf)
  add_test_ptscotch(test_scotch_dgraph_redist ${CMAKE_CURRENT_SOURCE_DIR}/data/bump_b100000.grf)
endif(BUILD_PTSCOTCH)

######################################
#  Tests of the sequential programs  #
######################################

# check_prog_acpl
add_test(NAME amk_acpl COMMAND ${BASH} -c "'$<TARGET_FILE:amk_m2>' 5 12 | '$<TARGET_FILE:acpl>' | '$<TARGET_FILE:atst>'")

# check_prog_amk
add_test(NAME amk_ccc COMMAND $<TARGET_FILE:amk_ccc> 5 ${dev_null})
add_test(NAME amk_fft2 COMMAND $<TARGET_FILE:amk_fft2> 5 ${dev_null})
add_test(NAME amk_hy COMMAND $<TARGET_FILE:amk_hy> 5 ${dev_null})
add_test(NAME amk_m2 COMMAND $<TARGET_FILE:amk_m2> 5 12 ${dev_null})
add_test(NAME amk_m2_mo COMMAND $<TARGET_FILE:amk_m2> -mo 5 12 ${dev_null})
add_test(NAME amk_p2 COMMAND $<TARGET_FILE:amk_p2> 5 7 ${dev_null})

# check_prog_amk_grf
add_test(NAME amk_grf_m16x16 COMMAND $<TARGET_FILE:amk_grf> ${dat}/m16x16.grf ${dev_null})
add_test(NAME amk_grf_m16x16_2 COMMAND $<TARGET_FILE:amk_grf> -2 ${dat}/m16x16.grf ${dev_null})
add_test(NAME amk_grf_m4x4 COMMAND $<TARGET_FILE:amk_grf> -l${dat}/m4x4_vertlist.txt ${dat}/m4x4.grf ${dev_null})
add_test(NAME amk_grf_m4x4_2 COMMAND $<TARGET_FILE:amk_grf> -l${dat}/m4x4_vertlist.txt -2 ${dat}/m4x4.grf ${dev_null})

# check_prog_atst
add_test(NAME atst_4x4x4 COMMAND $<TARGET_FILE:atst> ${tgt}/t4x4x4.tgt)

# check_prog_gbase
add_test(NAME test_gbase_1 COMMAND ${BASH} -c
  "'$<TARGET_FILE:gbase>' 1 '${dat}/m4x4.grf' m4x4_b1.grf && \
  '${CMAKE_COMMAND}' -E compare_files  '${dat}/m4x4_b1.grf' m4x4_b1.grf")
add_test(NAME test_gbase_2 COMMAND ${BASH} -c
  "'$<TARGET_FILE:gbase>' 1 '${dat}/m4x4_b1.grf' m4x4_b1.grf && \
  '${CMAKE_COMMAND}' -E compare_files  '${dat}/m4x4_b1.grf' m4x4_b1.grf")
add_test(NAME test_gbase_3 COMMAND ${BASH} -c
  "'$<TARGET_FILE:gbase>' 0 '${dat}/m4x4_b1.grf' m4x4.grf && \
  '${CMAKE_COMMAND}' -E compare_files  '${dat}/m4x4.grf' m4x4.grf")

# check_prog_gmk
add_test(NAME gmk_hy COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_hy>' 5 | '$<TARGET_FILE:gtst>'")
add_test(NAME gmk_m2 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_m2>' 5 1 | '$<TARGET_FILE:gtst>'")
add_test(NAME gmk_m2_b1 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_m2>' -b1 5 4 -g${dev_null} | '$<TARGET_FILE:gtst>'")
add_test(NAME gmk_m2_t COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_m2>' -t 5 4 | '$<TARGET_FILE:gtst>'")
add_test(NAME gmk_m2_b1_t COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_m2>' -b1 -t 5 1 | '$<TARGET_FILE:gtst>'")
add_test(NAME gmk_m3 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_m3>' -b1 5 4 3 -g${dev_null} | '$<TARGET_FILE:gtst>'")
add_test(NAME gmk_ub2 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_ub2>' 5 | '$<TARGET_FILE:gtst>'")

# check_prog_gmk_gmsh
foreach(msh_arg "5_1" "5_4")
  string(REPLACE "_" " " args ${msh_arg})
  add_test(NAME gmk_msh_${msh_arg} COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m2>' ${args} ${msh_arg}.msh && \
    '$<TARGET_FILE:gmk_msh>' ${msh_arg}.msh ${msh_arg}.grf && '$<TARGET_FILE:gtst>' ${msh_arg}.grf")
endforeach()
add_test(NAME gmk_msh_3_5_4_3 COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m3>' 5 4 3 3_5_4_3.msh && \
  '$<TARGET_FILE:gmk_msh>' 3_5_4_3.msh 3_5_4_3.grf && '$<TARGET_FILE:gtst>' 3_5_4_3.grf")
add_test(NAME gmk_msh_ship001 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_msh>' '${dat}/ship001.msh' ship001.grf && \
  '$<TARGET_FILE:gtst>' ship001.grf")

foreach(i RANGE 1 2)
  add_test(NAME gmk_msh_dual_${i} COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_msh>' '${dat}/small_55.msh' small_55_d${i}.grf -d${i} && \
    '${CMAKE_COMMAND}' -E compare_files  '${dat}/small_55_d${i}.grf' small_55_d${i}.grf")
endforeach()
add_test(NAME gmk_msh_dual_3 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_msh>' '${dat}/small_55.msh' small_55_d3.grf -d3 && \
  '${CMAKE_COMMAND}' -E compare_files  '${dat}/small_55_d2.grf' small_55_d3.grf")

# check_prog_mmk
add_test(NAME mmk_m2_5_1 COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m2>' 5 1 m2_5_1.msh && \
  '$<TARGET_FILE:mtst>' m2_5_1.msh")
add_test(NAME mmk_m2_5_4 COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m2>' 5 4 m2_5_4.msh && \
  '$<TARGET_FILE:mtst>' m2_5_4.msh")
add_test(NAME mmk_m3_5_4_3 COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m3>' 5 4 3 m3_5_4_3.msh && \
  '$<TARGET_FILE:mtst>' m3_5_4_3.msh")

# check_prog_mcv
foreach(msh_arg "5_1" "5_4")
  string(REPLACE "_" " " args ${msh_arg})
  add_test(NAME mcv_msh_${msh_arg} COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m2>' ${args} m2_${msh_arg}.msh && \
    '$<TARGET_FILE:mcv>' m2_${msh_arg}.msh output_m2_${msh_arg}.msh -is -os && \
    '${CMAKE_COMMAND}' -E compare_files m2_${msh_arg}.msh output_m2_${msh_arg}.msh")
endforeach()
add_test(NAME mcv_2_5_4_3 COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m3>' 5 4 3 m2_5_4_3.msh && \
  '$<TARGET_FILE:mcv>' m2_5_4_3.msh output_m2_5_4_3.msh -is -os && \
  '${CMAKE_COMMAND}' -E compare_files m2_5_4_3.msh output_m2_5_4_3.msh")

# check_prog_mord
foreach(msh_arg "5_1" "5_4")
  string(REPLACE "_" " " args ${msh_arg})
  add_test(NAME mord_${msh_arg} COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m2>' ${args} m2_${msh_arg}.msh && \
    '$<TARGET_FILE:gmk_msh>' m2_${msh_arg}.msh m2_${msh_arg}.grf && \
    '$<TARGET_FILE:mord>' m2_${msh_arg}.msh | '$<TARGET_FILE:gotst>' m2_${msh_arg}.grf")
endforeach()
add_test(NAME mord_5_4_3 COMMAND ${BASH} -c "'$<TARGET_FILE:mmk_m3>' 5 4 3 m2_5_4_3.msh && \
  '$<TARGET_FILE:gmk_msh>' m2_5_4_3.msh m2_5_4_3.grf && \
  '$<TARGET_FILE:mord>' m2_5_4_3.msh | '$<TARGET_FILE:gotst>' m2_5_4_3.grf")

add_test(NAME mord_ship001 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_msh>' '${dat}/ship001.msh' ship001.grf && \
  '$<TARGET_FILE:mord>' '${dat}/ship001.msh' | '$<TARGET_FILE:gotst>' ship001.grf")

# check_prog_gmap
add_test(NAME gmap_bump COMMAND $<TARGET_FILE:gmap> ${dat}/bump.grf ${tgt}/h3.tgt bump_h3.map -vmt)
add_test(NAME gmap_small COMMAND $<TARGET_FILE:gmap> ${dat}/small.grf ${tgt}/m11x13.tgt small_m11x13.map -vmt)

# check_prog_gord
add_test(NAME gord_bump COMMAND $<TARGET_FILE:gord> ${dat}/bump.grf ${dev_null} -vt)
add_test(NAME gord_bump_b1 COMMAND $<TARGET_FILE:gord> ${dat}/bump_b1.grf ${dev_null} -vt)
add_test(NAME gord_cmplx COMMAND ${BASH} -c "'$<TARGET_FILE:gord>' '${dat}/bump_b1.grf' '-On{sep=(/((vert)>(240))?(e{strat=g}|g|h);),ole=d{cmin=15,cmax=100000,frat=0},ose=g{pass=3}}' bump_b1_2.ord -vt")
add_test(NAME gord_bump_imbal COMMAND $<TARGET_FILE:gord> ${dat}/bump_imbal_32.grf bump_imbal_32.ord -vt)
add_test(NAME gord_ship001 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_msh>' '${dat}/ship001.msh' ship001.grf && \
  '$<TARGET_FILE:gord>' ship001.grf ship001.ord -vt")

# check_prog_gotst
foreach(grf bump bump_b1 bump_imbal_32)
  add_test(NAME gotst_${grf} COMMAND ${BASH} -c "'$<TARGET_FILE:gord>' '${dat}/${grf}.grf' ${grf}.ord -vt && \
    '$<TARGET_FILE:gotst>' '${dat}/${grf}.grf' ${grf}.ord")
endforeach()
add_test(NAME gotst_ship001 COMMAND ${BASH} -c "'$<TARGET_FILE:gmk_msh>' '${dat}/ship001.msh' ship001.grf && \
  '$<TARGET_FILE:gord>' ship001.grf ship001.ord -vt && \
  '$<TARGET_FILE:gotst>' ship001.grf ship001.ord ")

# check_prog_gord_block
add_test(NAME gord_block_1 COMMAND $<TARGET_FILE:gord> ${dat}/m4x4_b1_ev.grf m4x4_b1.ord '-ob{cmin=3,strat=s}')
add_test(NAME gord_block_2 COMMAND $<TARGET_FILE:gord> ${dat}/m4x4_b1_ev.grf m4x4_b1_ev.ord '-ob{cmin=2,strat=s}')
add_test(NAME gord_block_3 COMMAND $<TARGET_FILE:gord> ${dat}/m4x4_b1_ev.grf m4x4_b1_ev.ord '-ob{cmin=1,strat=s}')

# check_prog_gout
add_test(NAME gout_1 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' 9 '${dat}/bump.grf' bump_k9.map -vmt && \
  '$<TARGET_FILE:gord>' '${dat}/bump.grf' bump.ord -vt && \
  '$<TARGET_FILE:gout>' '-Om{e}' -gn '${dat}/bump.grf' ${dev_null} bump.ord bump_ord.eps")
add_test(NAME gout_2 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' 9 '${dat}/bump.grf' bump_k9.map -vmt && \
  '$<TARGET_FILE:gout>' '-Op{d,e,g,r,s}' '${dat}/bump.grf' '${dat}/bump.xyz' bump_k9.map bump_k9.eps")
add_test(NAME gout_3 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' 9 '${dat}/bump.grf' bump_k9.map -vmt && \
  '$<TARGET_FILE:gout>' '-Oi{c,r}' '${dat}/bump.grf' '${dat}/bump.xyz' bump_k9.map bump_k9.iv")

# check_prog_gpart
add_test(NAME gpart_1 COMMAND $<TARGET_FILE:gpart> 9 ${dat}/bump.grf bump_k9.map -vmt)
add_test(NAME gpart_2 COMMAND $<TARGET_FILE:gpart> 9 ${dat}/bump_b100000.grf bump_k9.map -vmt)

# check_prog_gpart_clustering
add_test(NAME gpart_cluster_1 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' -q 1 '${dat}/bump.grf' bump_part_cls_9.map -vmt && \
  '$<TARGET_FILE:gpart>' -q 1 '${dat}/bump_b100000.grf' bump_b100000_part_cls_9.map -vmt && \
  '$<TARGET_FILE:gmtst>' '${dat}/bump_b100000.grf' '${tgt}/vcmplt.tgt' bump_b100000_part_cls_9.map")
add_test(NAME gpart_cluster_2 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' -q 1 '${dat}/bump_imbal_32.grf' bump_imbal_32_part_cls_32.map -vmt && \
  '$<TARGET_FILE:gmtst>' '${dat}/bump_imbal_32.grf' '${tgt}/vhcub.tgt' bump_imbal_32_part_cls_32.map")

# check_prog_gpart_overlap
add_test(NAME gpart_overlap_1 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' -o 9 '${dat}/bump.grf' bump_part_ovl_9.map -vmt && \
  '$<TARGET_FILE:gpart>' -o 9 '${dat}/bump_b100000.grf' bump_b100000_part_ovl_9.map -vmt && \
  '$<TARGET_FILE:gmtst>' -o '${dat}/bump_b100000.grf' '${tgt}/k9.tgt' bump_b100000_part_ovl_9.map")
add_test(NAME gpart_overlap_2 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' -o 32 '${dat}/bump_imbal_32.grf' bump_imbal_32_part_ovl_32.map -vmt && \
  '$<TARGET_FILE:gmtst>' '${dat}/bump_imbal_32.grf' '${tgt}/k32.tgt' bump_imbal_32_part_ovl_32.map")
add_test(NAME gpart_overlap_3 COMMAND ${BASH} -c "'$<TARGET_FILE:gpart>' -o 6 '${dat}/bump_imbal_32.grf' bump_imbal_32_part_ovl_32_6.map '-me{strat=m{asc=b{width=3,bnd=d{pass=40,dif=1,rem=0}f{move=80,pass=-1,bal=0.01},org=f{move=80,pass=-1,bal=0.01}},low=r{job=t,bal=0.01,map=t,poli=S,sep=(m{asc=b{bnd=f{move=120,pass=-1,bal=0.01,type=b},org=f{move=120,pass=-1,bal=0.01,type=b},width=3},low=h{pass=10}f{move=120,pass=-1,bal=0.01,type=b},vert=120,rat=0.8}|m{asc=b{bnd=f{move=120,pass=-1,bal=0.01,type=b},org=f{move=120,pass=-1,bal=0.01,type=b},width=3},low=h{pass=10}f{move=120,pass=-1,bal=0.01,type=b},vert=120,rat=0.8})},vert=10000,rat=0.8,type=h}}f' -vmt")

# check_prog_gpart_remap
add_test(NAME gpart_remap COMMAND $<TARGET_FILE:gpart> 32 ${dat}/bump_imbal_32.grf ${dev_null} -ro${dat}/bump_old.map -vmt)

# check_prog_gscat
add_test(NAME gscat_1 COMMAND $<TARGET_FILE:gscat> 4 ${dat}/small.grf small_%r.grf)
add_test(NAME gscat_2 COMMAND $<TARGET_FILE:gscat> 4 -i1 ${dat}/bump.grf bump4_%r.grf)

# check_prog_gtst
add_test(NAME gtst COMMAND ${BASH} -c "'$<TARGET_FILE:gtst>' '${dat}/m4x4_b1_ev.grf' gtst_ev.txt && \
  '$<TARGET_FILE:gtst>' '${dat}/m4x4_b1_elv.grf' gtst_elv.txt && \
  '${CMAKE_COMMAND}' -E compare_files gtst_ev.txt gtst_elv.txt" )

#check_prog_full
add_test(NAME prg_full_1 COMMAND $<TARGET_FILE:gord> ${dat}/bump.grf bump_f.ord -vt)
add_test(NAME prg_full_2 COMMAND $<TARGET_FILE:gpart> 3 ${dat}/bump.grf bump_f.map -vt)
add_test(NAME prg_full_3 COMMAND $<TARGET_FILE:gpart> -o 3 ${dat}/bump.grf bump_fo.map -vt)

####################################
#  Tests of the parallel programs  #
####################################

if(BUILD_PTSCOTCH)
  # check_prog_dgord
  add_test(NAME dgord_1 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:dgord> ${dat}/bump.grf ${dev_null} -vt)
  add_test(NAME dgord_2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:dgord> ${dat}/bump_b100000.grf ${dev_null} -vt)
  add_test(NAME dgord_3 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:dgord> ${dat}/small_%r.grf small.ord -vt)

  # check_prog_dgpart
  add_test(NAME dgpart_1 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:dgpart> 9 ${dat}/bump.grf ${dev_null} -vmt)
  add_test(NAME dgpart_2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:dgpart> 9 ${dat}/bump_b100000.grf ${dev_null} -vmt)
  add_test(NAME dgpart_3 COMMAND ${BASH} -c "'${MPIEXEC_EXECUTABLE}' ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4 '$<TARGET_FILE:dgpart>' 50 '${dat}/small_%r.grf' ${dev_null} -vmt")

  # check_prog_dgscat-dggath
  add_test(NAME test_dgscat COMMAND
    ${BASH} -c "'${MPIEXEC_EXECUTABLE}' ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 '$<TARGET_FILE:dgscat>' \
      -r1 -c '${dat}/m16x16.grf' out_dgscat_%r.grf ; \
      '${MPIEXEC_EXECUTABLE}' ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 '$<TARGET_FILE:dggath>' \
      -r1 -c out_dgscat_%r.grf out_dggat.grf ; \
      '${CMAKE_COMMAND}' -E compare_files  '${dat}/m16x16.grf' out_dggat.grf")

  add_test(NAME test_dgscat_b1 COMMAND
    ${BASH} -c "'${MPIEXEC_EXECUTABLE}' ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 '$<TARGET_FILE:dgscat>' \
      -c '${dat}/m16x16_b1.grf' out_b1_dgscat_%r.grf ; \
      '${MPIEXEC_EXECUTABLE}' ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3 '$<TARGET_FILE:dggath>' \
      -c out_b1_dgscat_%r.grf out_b1_ddgat.grf ; \
      '${CMAKE_COMMAND}' -E compare_files '${dat}/m16x16_b1.grf' out_b1_ddgat.grf")

  # check_prog_dgtst
  add_test(NAME dgtst COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:dgtst> ${dat}/small_%r.grf)

  # check_prog_dfull
  add_test(NAME dfull_1 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:dgord> ${dat}/bump.grf bump_f.ord -vt)
  add_test(NAME dfull_2 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:dgpart> 3 ${dat}/bump.grf bump_f.map -vmt)
endif(BUILD_PTSCOTCH)
