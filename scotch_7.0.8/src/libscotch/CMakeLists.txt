## Copyright 2014-2016,2021-2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
##                                 to     10 jun 2025     ##
##                                                        ##
############################################################

################
#  dummysizes  #
################

add_executable(dummysizes dummysizes.c)
set_target_properties(dummysizes PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if(BUILD_PTSCOTCH)
  find_package(MPI COMPONENTS C)
  if(NOT MPI_C_FOUND)
    message(FATAL_ERROR "MPI required to compile PT-Scotch")
  endif()
  add_executable(ptdummysizes dummysizes.c)
  set_target_properties(ptdummysizes PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(ptdummysizes PRIVATE MPI::MPI_C)
  add_dependencies(ptdummysizes scotch_h)
  set_target_properties(ptdummysizes PROPERTIES
    COMPILE_FLAGS -DSCOTCH_PTSCOTCH)
endif()

##############
#  scotch.h  #
##############

# Generate scotch.h
add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/scotch.h
  COMMAND $<TARGET_FILE:dummysizes> "-s${SCOTCH_NAME_SUFFIX}"
  ${CMAKE_CURRENT_SOURCE_DIR}/library.h ${GENERATED_INCLUDE_DIR}/scotch.h DEPENDS dummysizes
  DEPENDS library.h
  COMMENT "Generate scotch.h")
add_custom_target(scotch_h
  DEPENDS "${GENERATED_INCLUDE_DIR}/scotch.h")

# Generate scotchf.h
add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/scotchf.h
  COMMAND $<TARGET_FILE:dummysizes>  "-s${SCOTCH_NAME_SUFFIX}"
  ${CMAKE_CURRENT_SOURCE_DIR}/library_f.h ${GENERATED_INCLUDE_DIR}/scotchf.h DEPENDS dummysizes
  DEPENDS library.h
  COMMENT "Generate scotchf.h")
add_custom_target(scotchf_h
  DEPENDS "${GENERATED_INCLUDE_DIR}/scotchf.h")

################
#  ptscotch.h  #
################

# Generate ptscotch.h
if(BUILD_PTSCOTCH)
  add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/ptscotch.h
    COMMAND $<TARGET_FILE:ptdummysizes> "-s${SCOTCH_NAME_SUFFIX}"
    ${CMAKE_CURRENT_SOURCE_DIR}/library_pt.h ${GENERATED_INCLUDE_DIR}/ptscotch.h DEPENDS ptdummysizes
    DEPENDS library_pt.h scotch_h
    COMMENT "Generate ptscotch.h")
  add_custom_target(ptscotch_h
    DEPENDS ${GENERATED_INCLUDE_DIR}/ptscotch.h)

# Generate ptscotchf.h
add_custom_command(OUTPUT ${GENERATED_INCLUDE_DIR}/ptscotchf.h
  COMMAND $<TARGET_FILE:ptdummysizes> "-s${SCOTCH_NAME_SUFFIX}"
  ${CMAKE_CURRENT_SOURCE_DIR}/library_pt_f.h ${GENERATED_INCLUDE_DIR}/ptscotchf.h DEPENDS ptdummysizes
  DEPENDS library.h
  COMMENT "Generate ptscotchf.h")
add_custom_target(ptscotchf_h
  DEPENDS "${GENERATED_INCLUDE_DIR}/ptscotchf.h")
endif()

#################
#  parser_yy.c  #
#################

# Generate parser_yy.c using Bison
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/parser_yy.c ${CMAKE_CURRENT_BINARY_DIR}/parser_ly.h
  COMMAND ${BISON_EXECUTABLE} "-pscotchyy${SCOTCH_NAME_SUFFIX}" --defines="parser_ly.h" --output="parser_yy.c" ${CMAKE_CURRENT_SOURCE_DIR}/parser_yy.y
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/parser_yy.y module.h common.h parser.h parser_yy.h
  COMMENT "Generate parser_yy.c and parser_ly.h using Bison")
add_custom_target(parser_yy_c DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/parser_yy.c ${CMAKE_CURRENT_BINARY_DIR}/parser_ly.h)

#################
#  parser_ll.c  #
#################

# Generate parser_ll.c using Flex
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/parser_ll.c
  ${CMAKE_CURRENT_BINARY_DIR}/parser_lh.h
  COMMAND ${FLEX_EXECUTABLE} "-Pscotchyy${SCOTCH_NAME_SUFFIX}"  "-oparser_ll.c"
  "--header-file=parser_lh.h"
  ${CMAKE_CURRENT_SOURCE_DIR}/parser_ll.l
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/parser_ll.l
  ${CMAKE_CURRENT_BINARY_DIR}/parser_lh.h module.h common.h parser.h parser_ll.h parser_ly.h
  COMMENT "Generate parser_ll.c using Flex")
add_custom_target(parser_ll_c DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/parser_ll.c)

###############
#  libScotch  #
###############

set(SCOTCH_C_SOURCES
  arch.c
  arch.h
  arch_build.c
  arch_build.h
  arch_build2.c
  arch_build2.h
  arch_cmplt.c
  arch_cmplt.h
  arch_cmpltw.c
  arch_cmpltw.h
  arch_deco.c
  arch_deco.h
  arch_deco2.c
  arch_deco2.h
  arch_dist.c
  arch_dist.h
  arch_hcub.c
  arch_hcub.h
  arch_mesh.c
  arch_mesh.h
  arch_tleaf.c
  arch_tleaf.h
  arch_sub.c
  arch_sub.h
  arch_torus.c
  arch_torus.h
  arch_vcmplt.c
  arch_vcmplt.h
  arch_vhcub.c
  arch_vhcub.h
  bgraph.c
  bgraph.h
  bgraph_bipart_bd.c
  bgraph_bipart_bd.h
  bgraph_bipart_df.c
  bgraph_bipart_df.h
  bgraph_bipart_ex.c
  bgraph_bipart_ex.h
  bgraph_bipart_fm.c
  bgraph_bipart_fm.h
  bgraph_bipart_gg.c
  bgraph_bipart_gg.h
  bgraph_bipart_gp.c
  bgraph_bipart_gp.h
  bgraph_bipart_ml.c
  bgraph_bipart_ml.h
  bgraph_bipart_st.c
  bgraph_bipart_st.h
  bgraph_bipart_zr.c
  bgraph_bipart_zr.h
  bgraph_cost.c
  bgraph_check.c
  bgraph_cost.c
  bgraph_store.c
  common.c
  common.h
  common_context.c
  common_error.c
  common_file.c
  common_file.h
  common_file_compress.c
  common_file_compress.h
  common_file_decompress.c
  common_integer.c
  common_memory.c
  common_string.c
  common_stub.c
  common_thread.c
  common_thread_system.c
  common_thread_system.h
  common_timer.c
  common_timer.h
  common_values.c
  common_values.h
  context.c
  context.h
  fibo.c
  fibo.h
  gain.c
  gain.h
  geom.c
  geom.h
  graph.c
  graph.h
  graph_band.c
  graph_base.c
  graph_check.c
  graph_clone.c
  graph_coarsen.c
  graph_coarsen.h
  graph_diam.c
  graph_ielo.c
  library_graph_induce.c
  graph_induce.c
  graph_induce.h
  graph_io.c
  graph_io.h
  graph_io_chac.c
  graph_io_habo.c
  graph_io_habo.h
  graph_io_mmkt.c
  graph_io_mmkt.h
  graph_io_scot.c
  graph_io_scot.h
  graph_list.c
  graph_list.h
  graph_match.c
  graph_match.h
  hall_order_hd.c
  hall_order_hd.h
  hall_order_hf.c
  hall_order_hf.h
  hall_order_hx.c
  hall_order_hx.h
  hgraph.c
  hgraph.h
  hgraph_check.c
  hgraph_induce.c
  hgraph_induce.h
  hgraph_order_bl.c
  hgraph_order_bl.h
  hgraph_order_cc.c
  hgraph_order_cc.h
  hgraph_order_cp.c
  hgraph_order_cp.h
  hgraph_order_gp.c
  hgraph_order_gp.h
  hgraph_order_hd.c
  hgraph_order_hd.h
  hgraph_order_hf.c
  hgraph_order_hf.h
  hgraph_order_hx.c
  hgraph_order_hx.h
  hgraph_order_kp.c
  hgraph_order_kp.h
  hgraph_order_nd.c
  hgraph_order_nd.h
  hgraph_order_si.c
  hgraph_order_si.h
  hgraph_order_st.c
  hgraph_order_st.h
  hmesh.c
  hmesh.h
  hmesh_check.c
  hmesh_hgraph.c
  hmesh_hgraph.h
  hmesh_induce.c
  hmesh_mesh.c
  hmesh_order_bl.c
  hmesh_order_bl.h
  hmesh_order_cp.c
  hmesh_order_cp.h
  hmesh_order_gp.c
  hmesh_order_gp.h
  hmesh_order_gr.c
  hmesh_order_gr.h
  hmesh_order_hd.c
  hmesh_order_hd.h
  hmesh_order_hf.c
  hmesh_order_hf.h
  hmesh_order_hx.c
  hmesh_order_hx.h
  hmesh_order_nd.c
  hmesh_order_nd.h
  hmesh_order_si.c
  hmesh_order_si.h
  hmesh_order_st.c
  hmesh_order_st.h
  kgraph.c
  kgraph.h
  kgraph_band.c
  kgraph_check.c
  kgraph_cost.c
  kgraph_map_bd.c
  kgraph_map_bd.h
  kgraph_map_cp.c
  kgraph_map_cp.h
  kgraph_map_df.c
  kgraph_map_df.h
  kgraph_map_ex.c
  kgraph_map_ex.h
  kgraph_map_fm.c
  kgraph_map_fm.h
  kgraph_map_ml.c
  kgraph_map_ml.h
  kgraph_map_rb.c
  kgraph_map_rb.h
  kgraph_map_rb_map.c
  kgraph_map_rb_map.h
  kgraph_map_rb_part.c
  kgraph_map_rb_part.h
  kgraph_map_st.c
  kgraph_map_st.h
  kgraph_store.c
  library.h
  library_arch.c
  library_arch_build.c
  library_arch_build_f.c
  library_arch_dom.c
  library_arch_dom_f.c
  library_arch_f.c
  library_common.c
  library_common_f.c
  library_context.c
  library_context_graph.c
  library_context_graph_f.c
  library_context_f.c
  #library_errcom.c
  #library_error.c
  #library_error_exit.c
  library_f.h
  library_geom.c
  library_geom_f.c
  library_graph.c
  library_graph_base.c
  library_graph_base_f.c
  library_graph_check.c
  library_graph_check_f.c
  library_graph_coarsen.c
  library_graph_coarsen_f.c
  library_graph_color.c
  library_graph_color_f.c
  library_graph_diam.c
  library_graph_diam_f.c
  library_graph_f.c
  library_graph_io_chac.c
  library_graph_io_chac_f.c
  library_graph_io_habo.c
  library_graph_io_habo_f.c
  library_graph_io_mmkt.c
  library_graph_io_mmkt_f.c
  library_graph_io_scot.c
  library_graph_io_scot_f.c
  library_graph_map.c
  library_graph_map_f.c
  library_graph_map_io.c
  library_graph_map_io.h
  library_graph_map_io_f.c
  library_graph_map_view.c
  library_graph_map_view.h
  library_graph_map_view_f.c
  library_graph_order.c
  library_graph_order_f.c
  library_graph_part_ovl.c
  library_graph_part_ovl_f.c
  library_mapping.c
  library_mapping.h
  library_memory.c
  library_mesh.c
  library_mesh_f.c
  library_mesh_graph.c
  library_mesh_graph_f.c
  library_mesh_io_habo.c
  library_mesh_io_habo_f.c
  library_mesh_io_scot.c
  library_mesh_io_scot_f.c
  library_mesh_order.c
  library_mesh_order_f.c
  library_order.c
  library_order.h
  library_parser.c
  library_parser_f.c
  library_pt.h
  library_pt_f.h
  library_random.c
  library_random_f.c
  library_strat.c
  library_version.c
  library_version_f.c
  mapping.c
  mapping.h
  mapping_check.c
  mapping_io.c
  mapping_io.h
  mesh.c
  mesh.h
  mesh_check.c
  mesh_coarsen.c
  mesh_coarsen.h
  mesh_graph.c
  mesh_graph.h
  mesh_induce_sepa.c
  mesh_induce_sepa.h
  mesh_io.c
  mesh_io.h
  mesh_io_habo.c
  mesh_io_scot.c
  module.h
  order.c
  order.h
  order_check.c
  order_io.c
  parser.c
  parser.h
  ${CMAKE_CURRENT_BINARY_DIR}/parser_ll.c
  ${CMAKE_CURRENT_BINARY_DIR}/parser_ly.h
  ${CMAKE_CURRENT_BINARY_DIR}/parser_yy.c
  parser_yy.h
  ${GENERATED_INCLUDE_DIR}/scotch.h
  vgraph.c
  vgraph.h
  vgraph_check.c
  vgraph_separate_bd.c
  vgraph_separate_bd.h
  vgraph_separate_df.c
  vgraph_separate_df.h
  vgraph_separate_es.c
  vgraph_separate_es.h
  vgraph_separate_fm.c
  vgraph_separate_fm.h
  vgraph_separate_gg.c
  vgraph_separate_gg.h
  vgraph_separate_gp.c
  vgraph_separate_gp.h
  vgraph_separate_ml.c
  vgraph_separate_ml.h
  vgraph_separate_st.c
  vgraph_separate_st.h
  vgraph_separate_th.c
  vgraph_separate_th.h
  vgraph_separate_vw.c
  vgraph_separate_vw.h
  vgraph_separate_zr.c
  vgraph_separate_zr.h
  vgraph_store.c
  vmesh.c
  vmesh.h
  vmesh_check.c
  vmesh_separate_fm.c
  vmesh_separate_fm.h
  vmesh_separate_gg.c
  vmesh_separate_gg.h
  vmesh_separate_gr.c
  vmesh_separate_gr.h
  vmesh_separate_ml.c
  vmesh_separate_ml.h
  vmesh_separate_st.c
  vmesh_separate_st.h
  vmesh_separate_zr.c
  vmesh_separate_zr.h
  vmesh_store.c
  wgraph.c
  wgraph.h
  wgraph_check.c
  wgraph_part_fm.c
  wgraph_part_fm.h
  wgraph_part_es.c
  wgraph_part_es.h
  wgraph_part_ml.c
  wgraph_part_ml.h
  wgraph_part_rb.c
  wgraph_part_rb.h
  wgraph_part_st.c
  wgraph_part_st.h
  wgraph_part_zr.c
  wgraph_part_zr.h
  wgraph_store.c)

add_library(scotch
  ${SCOTCH_C_SOURCES})
set_target_properties(scotch PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
  SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})
add_dependencies(scotch parser_yy_c parser_ll_c)

# scotcherr lib
if(NOT "${LIBSCOTCHERR}" STREQUAL "")
  target_link_libraries(scotch PUBLIC ${LIBSCOTCHERR})
endif()

if(NOT WIN32)
  target_link_libraries(scotch PUBLIC m)
endif()

target_include_directories(scotch PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(APPLE)
  set_target_properties(scotch PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

if(OS_IS_FREEBSD)
  target_compile_definitions(scotch PUBLIC cpu_set_t=cpuset_t __BSD_VISIBLE)
endif()

if(SCOTCH_DETERMINISTIC STREQUAL FULL)
  target_compile_definitions(scotch PUBLIC COMMON_RANDOM_FIXED_SEED)
  target_compile_definitions(scotch PUBLIC SCOTCH_DETERMINISTIC)
  message(STATUS "Determinism: full")
elseif(SCOTCH_DETERMINISTIC STREQUAL FIXED_SEED)
  target_compile_definitions(scotch PUBLIC COMMON_RANDOM_FIXED_SEED)
  message(STATUS "Determinism: fixed random seed")
else()
  message(STATUS "Determinism: none")
endif()

# Add (de)compression libraries
if(ZLIB_FOUND)
  target_compile_definitions(scotch PRIVATE COMMON_FILE_COMPRESS_GZ)
  target_include_directories(scotch PRIVATE ${ZLIB_INCLUDE_DIRS})
  target_link_libraries(scotch PRIVATE ${ZLIB_LIBRARIES})
endif()

if(BZIP2_FOUND)
  target_compile_definitions(scotch PRIVATE COMMON_FILE_COMPRESS_BZ2)
  target_include_directories(scotch PRIVATE ${BZIP2_INCLUDE_DIRS})
  target_link_libraries(scotch PRIVATE ${BZIP2_LIBRARIES})
endif()

if(LIBLZMA_FOUND)
  target_compile_definitions(scotch PRIVATE COMMON_FILE_COMPRESS_LZMA)
  target_include_directories(scotch PRIVATE ${LIBLZMA_INCLUDE_DIRS})
  target_link_libraries(scotch PRIVATE ${LIBLZMA_LIBRARIES})
endif()

# Add thread library
if(Threads_FOUND)
  target_compile_definitions(scotch PRIVATE COMMON_PTHREAD SCOTCH_PTHREAD)
  target_link_libraries(scotch PRIVATE Threads::Threads ${THREADS_PTHREADS_WIN32_LIBRARY})
endif()

# Include files
add_dependencies(scotch scotch_h scotchf_h)
if(BUILD_PTSCOTCH)
  add_dependencies(scotch ptscotch_h ptscotchf_h)
endif()

# Error management libraries
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  target_compile_definitions(scotch PRIVATE SCOTCH_DEBUG_ALL)
endif()

add_library(scotcherr library_error.c)
set_target_properties(scotcherr PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
  SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})
target_include_directories(scotcherr PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_dependencies(scotcherr scotch_h)

add_library(scotcherrexit library_error_exit.c)
set_target_properties(scotcherrexit PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
  SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})
target_include_directories(scotcherrexit PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_dependencies(scotcherrexit scotch_h)

#################
#  libPTScotch  #
#################

set(PTSCOTCH_C_SOURCES
  arch_mpi.c
  arch_mpi.h
  arch_cmplt_mpi.c
  arch_cmplt_mpi.h
  arch_cmpltw_mpi.c
  arch_cmpltw_mpi.h
  arch_deco_mpi.c
  arch_deco_mpi.h
  arch_deco2_mpi.c
  arch_deco2_mpi.h
  arch_dist_mpi.c
  arch_dist_mpi.h
  arch_hcub_mpi.c
  arch_hcub_mpi.h
  arch_mesh_mpi.c
  arch_mesh_mpi.h
  arch_sub_mpi.c
  arch_sub_mpi.h
  arch_tleaf_mpi.c
  arch_tleaf_mpi.h
  arch_torus_mpi.h
  arch_vcmplt_mpi.c
  arch_vcmplt_mpi.h
  arch_vhcub_mpi.c
  arch_vhcub_mpi.h
  bdgraph.c
  bdgraph.h
  bdgraph_bipart_bd.c
  bdgraph_bipart_bd.h
  bdgraph_bipart_df.c
  bdgraph_bipart_df.h
  bdgraph_bipart_ex.c
  bdgraph_bipart_ex.h
  bdgraph_bipart_ml.c
  bdgraph_bipart_ml.h
  bdgraph_bipart_sq.c
  bdgraph_bipart_sq.h
  bdgraph_bipart_st.c
  bdgraph_bipart_st.h
  bdgraph_bipart_zr.c
  bdgraph_bipart_zr.h
  bdgraph_check.c
  bdgraph_gather_all.c
  bdgraph_store.c
  comm.c
  comm.h
  dgraph.c
  dgraph.h
  dgraph_allreduce.c
  dgraph_allreduce.h
  dgraph_band.c
  dgraph_band_grow.h
  dgraph_build.c
  dgraph_build.h
  dgraph_build_grid3d.c
  dgraph_build_grid3d.h
  dgraph_build_hcub.c
  dgraph_check.c
  dgraph_coarsen.c
  dgraph_coarsen.h
  dgraph_compact.c
  dgraph_fold.c
  dgraph_fold.h
  dgraph_fold_comm.c
  dgraph_fold_comm.h
  dgraph_fold_dup.c
  dgraph_fold_dup.h
  dgraph_gather.c
  dgraph_gather_all.c
  dgraph_ghst.c
  dgraph_ghst.h
  dgraph_halo.c
  dgraph_halo.h
  dgraph_induce.c
  dgraph_io_load.c
  dgraph_io_load.h
  dgraph_io_save.c
  dgraph_match.c
  dgraph_match.h
  dgraph_match_check.c
  dgraph_match_sync_coll.c
  dgraph_match_sync_ptop.c
  dgraph_redist.c
  dgraph_redist.h
  dgraph_scatter.c
  dgraph_view.c
  dmapping.c
  dmapping.h
  dmapping_io.c
  dorder.c
  dorder.h
  dorder_gather.c
  dorder_gather.h
  dorder_io.c
  dorder_io_block.c
  dorder_io_tree.c
  dorder_perm.c
  dorder_perm.h
  dorder_tree_dist.c
  hdgraph.c
  hdgraph.h
  hdgraph_check.c
  hdgraph_fold.c
  hdgraph_fold.h
  hdgraph_gather.c
  hdgraph_induce.c
  hdgraph_order_nd.c
  hdgraph_order_nd.h
  hdgraph_order_si.c
  hdgraph_order_si.h
  hdgraph_order_sq.c
  hdgraph_order_sq.h
  hdgraph_order_st.c
  hdgraph_order_st.h
  kdgraph.c
  kdgraph.h
  kdgraph_gather.c
  kdgraph_map_rb.c
  kdgraph_map_rb.h
  kdgraph_map_rb_map.c
  kdgraph_map_rb_map.h
  kdgraph_map_rb_part.c
  kdgraph_map_rb_part.h
  kdgraph_map_st.c
  kdgraph_map_st.h
  library_context_dgraph.c
  library_context_dgraph_f.c
  library_dgraph.c
  library_dgraph_band.c
  library_dgraph_band_f.c
  library_dgraph_build.c
  library_dgraph_build_f.c
  library_dgraph_build_grid3d.c
  library_dgraph_build_grid3d_f.c
  library_dgraph_check.c
  library_dgraph_check_f.c
  library_dgraph_coarsen.c
  library_dgraph_coarsen_f.c
  library_dgraph_f.c
  library_dgraph_gather.c
  library_dgraph_gather_f.c
  library_dgraph_grow.c
  library_dgraph_halo.c
  library_dgraph_halo_f.c
  library_dgraph_induce.c
  library_dgraph_induce_f.c
  library_dgraph_io_load.c
  library_dgraph_io_load_f.c
  library_dgraph_io_save.c
  library_dgraph_io_save_f.c
  library_dgraph_map.c
  library_dgraph_map_f.c
  library_dgraph_map_stat.c
  library_dgraph_map_stat_f.c
  library_dgraph_map_view.c
  library_dgraph_map_view_f.c
  library_dgraph_order.c
  library_dgraph_order_f.c
  library_dgraph_order_gather.c
  library_dgraph_order_gather_f.c
  library_dgraph_order_io.c
  library_dgraph_order_io_block.c
  library_dgraph_order_io_block_f.c
  library_dgraph_order_io_f.c
  library_dgraph_order_perm.c
  library_dgraph_order_perm_f.c
  library_dgraph_order_tree_dist.c
  library_dgraph_order_tree_dist_f.c
  library_dgraph_redist.c
  library_dgraph_redist_f.c
  library_dgraph_scatter.c
  library_dgraph_scatter_f.c
  library_dgraph_stat.c
  library_dgraph_stat.h
  library_dgraph_stat_f.c
  library_dmapping.c
  library_dmapping.h
  library_dorder.c
  ${GENERATED_INCLUDE_DIR}/scotch.h
  ${GENERATED_INCLUDE_DIR}/ptscotch.h
  vdgraph.c
  vdgraph.h
  vdgraph_check.c
  vdgraph_gather_all.c
  vdgraph_separate_bd.c
  vdgraph_separate_bd.h
  vdgraph_separate_df.c
  vdgraph_separate_df.h
  vdgraph_separate_ml.c
  vdgraph_separate_ml.h
  vdgraph_separate_sq.c
  vdgraph_separate_sq.h
  vdgraph_separate_st.c
  vdgraph_separate_st.h
  vdgraph_separate_zr.c
  vdgraph_separate_zr.h
  vdgraph_store.c)

if(BUILD_PTSCOTCH)
  add_library(ptscotch
    ${PTSCOTCH_C_SOURCES})

  set_target_properties(ptscotch PROPERTIES
    VERSION ${SCOTCH_VERSION_LONG}
    SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE}
    COMPILE_FLAGS -DSCOTCH_PTSCOTCH)

  target_link_libraries(ptscotch PUBLIC scotch MPI::MPI_C)

  target_include_directories(ptscotch PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  add_dependencies(ptscotch ptscotch_h scotch_h)

  if(APPLE)
    set_target_properties(ptscotch PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif()

  if(Threads_FOUND)
    target_compile_definitions(ptscotch PRIVATE SCOTCH_PTHREAD COMMON_PTHREAD)
    target_link_libraries(ptscotch PRIVATE Threads::Threads ${THREADS_PTHREADS_WIN32_LIBRARY})
  endif()

  add_library(ptscotcherr library_error.c)
  set_target_properties(ptscotcherr PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
    SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})
  target_include_directories(ptscotcherr PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  add_dependencies(ptscotcherr ptscotch_h scotch_h)

  add_library(ptscotcherrexit library_error_exit.c)
  set_target_properties(ptscotcherrexit PROPERTIES VERSION ${SCOTCH_VERSION_LONG}
    SOVERSION ${SCOTCH_VERSION}.${SCOTCH_RELEASE})
  target_include_directories(ptscotcherrexit PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

  add_dependencies(ptscotcherrexit ptscotch_h scotch_h)
  if(NOT "${LIBPTSCOTCHERR}" STREQUAL "")
    target_link_libraries(ptscotch PUBLIC ${LIBPTSCOTCHERR})
  endif()
endif()

set_target_properties(scotch PROPERTIES PUBLIC_HEADER "${GENERATED_INCLUDE_DIR}/scotch.h;${GENERATED_INCLUDE_DIR}/scotchf.h")
if(BUILD_PTSCOTCH)
  set_target_properties(ptscotch PROPERTIES PUBLIC_HEADER "${GENERATED_INCLUDE_DIR}/ptscotch.h;${GENERATED_INCLUDE_DIR}/ptscotchf.h")
endif()

set(TARGETS_LIST "scotch;scotcherr;scotcherrexit")

foreach(_target ${TARGETS_LIST})
  install(EXPORT ${_target}Targets
    COMPONENT libscotch
    FILE ${_target}Targets.cmake
    NAMESPACE SCOTCH::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scotch)

  install(TARGETS ${_target}
    COMPONENT libscotch
    EXPORT ${_target}Targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  if(BUILD_PTSCOTCH)
    install(EXPORT pt${_target}Targets
      COMPONENT libptscotch
      FILE pt${_target}Targets.cmake
      NAMESPACE SCOTCH::
      DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scotch)

    install(TARGETS pt${_target}
      COMPONENT libptscotch
      EXPORT pt${_target}Targets
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
endforeach()
