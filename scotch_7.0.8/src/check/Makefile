## Copyright 2011,2012,2014-2016,2018-2021 IPB, Universite de Bordeaux, INRIA & CNRS
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

TMPDIR			= /tmp

SCOTCHBINDIR		= ../../bin
SCOTCHINCDIR		= ../../include
SCOTCHLIBDIR		= ../../lib
SCOTCHTGTDIR		= ../../tgt

LIBESMUMPS		= -lesmumps -lscotch -lscotcherr
LIBPTSCOTCH		= -lptscotch -lscotch -lptscotcherr
LIBSCOTCH		= -lscotch -lscotcherr

DEPESMUMPS		= $(SCOTCHLIBDIR)/libesmumps$(LIB)
DEPPTSCOTCH		= $(SCOTCHLIBDIR)/libptscotch$(LIB)
DEPSCOTCH		= $(SCOTCHLIBDIR)/libscotch$(LIB)

include ../Makefile.inc

ifeq ($(strip $(SCOTCH_NAME_SUFFIX)),)
  ifneq ($(findstring -DSCOTCH_NAME_SUFFIX=,$(CFLAGS)),)
    scotch_name_suffix_defs = $(filter -DSCOTCH_NAME_SUFFIX=%,$(CFLAGS))
    scotch_name_suffix_def = $(lastword $(scotch_name_suffix_defs))
    export SCOTCH_NAME_SUFFIX = $(subst -DSCOTCH_NAME_SUFFIX=,,$(scotch_name_suffix_def))
  endif
endif
ifneq ($(strip $(SCOTCH_NAME_SUFFIX)),)
  $(info Using suffix "$(SCOTCH_NAME_SUFFIX)")
  ifeq ($(findstring -DSCOTCH_NAME_SUFFIX=,$(CFLAGS)),)
    CFLAGS += -DSCOTCH_NAME_SUFFIX=$(SCOTCH_NAME_SUFFIX)
  endif
  ifeq ($(findstring -DSCOTCH_NAME_SUFFIX=,$(CCDFLAGS)),)
    CCDFLAGS += -DSCOTCH_NAME_SUFFIX=$(SCOTCH_NAME_SUFFIX)
  endif
endif

ifeq ($(strip $(SCOTCH_METIS_PREFIXFU)),)
  ifneq ($(findstring -DSCOTCH_METIS_PREFIXFU=,$(CFLAGS)),)
    scotch_metis_prefixfu_defs = $(filter -DSCOTCH_METIS_PREFIXFU=%,$(CFLAGS))
    scotch_metis_prefixfu_def = $(lastword $(scotch_metis_prefixfu_defs))
    export SCOTCH_METIS_PREFIXFU = $(subst -DSCOTCH_METIS_PREFIXFU=,,$(scotch_metis_prefixfu_def))
  endif
endif
ifeq ($(strip $(SCOTCH_METIS_PREFIXFU)),)
  ifeq ($(findstring -DSCOTCH_METIS_PREFIX,$(CFLAGS)),)
    export SCOTCH_METIS_PREFIXFU =
  else
    export SCOTCH_METIS_PREFIXFU = SCOTCHF
  endif
endif
ifneq ($(strip $(SCOTCH_METIS_PREFIXFU)),)
  $(info Using MeTiS uppercase Fortran prefix "$(SCOTCH_METIS_PREFIXFU)")
endif

EXECPPROG		?= mpirun
EXECPFLAGS		?= $(if $(findstring Open,$(shell $(EXECPPROG) --version)),--oversubscribe)
# EXECPFLAGS		+= valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes
EXECS			?=
# EXECS			+= valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes
EXECP3			?= $(EXECPPROG) -n 3 $(EXECPFLAGS)
EXECP4			?= $(EXECPPROG) -n 4 $(EXECPFLAGS)

PROGDIFF		?= diff -q
PROGSED			?= sed

##
##  General inference rules.
##

%$(OBJ)	:	%.c
				$(CC) $(CFLAGS) $(CLIBFLAGS) -I$(SCOTCHINCDIR) -c $(<) -o $(@)

%_v3$(EXE) :	%.c
				$(CC) $(CFLAGS) -I$(SCOTCHINCDIR) -L$(SCOTCHLIBDIR) $(<) -o $(@) -lscotchmetisv3 $(SCOTCHLIBS) $(LDFLAGS)

%_v5$(EXE) :	%.c
				$(CC) $(CFLAGS) -I$(SCOTCHINCDIR) -L$(SCOTCHLIBDIR) $(<) -o $(@) -lscotchmetisv5 $(SCOTCHLIBS) $(LDFLAGS)

%$(EXE)	:	%.c
				$(CC) $(CFLAGS) -I$(SCOTCHINCDIR) -L$(SCOTCHLIBDIR) $(<) -o $(@) $(SCOTCHLIBS) $(LDFLAGS)

%.f90	:	%.f90.in
				../libscotch/dummysizes$(EXE) -s$(SCOTCH_NAME_SUFFIX) $(<) $(@).tmp
				$(PROGSED) -e 's/SCOTCHMETISNAMEFU\s*(\([A-Za-z0-8_]*\))/'"$(SCOTCH_METIS_PREFIXFU)"'\1/g' $(@).tmp > $(@)

%_v3$(EXE):	%.f90
				$(FC) $(FFLAGS) -I$(SCOTCHINCDIR) -L$(SCOTCHLIBDIR) $(<) -o $(@) -lscotchmetisv3 $(SCOTCHLIBS) $(LDFLAGS)

%_v5$(EXE):	%.f90
				$(FC) $(FFLAGS) -I$(SCOTCHINCDIR) -L$(SCOTCHLIBDIR) $(<) -o $(@) -lscotchmetisv5 $(SCOTCHLIBS) $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	check	ptcheck	clean	realclean

check				:	realclean
					$(MAKE) CC="$(CCS)" CCD="$(CCS)" SCOTCHLIBS="$(LIBSCOTCH)" $(CHECKSCOTCH)

check_full			:	realclean
					$(MAKE) CC="$(CCS) -DSCOTCH_DEBUG_FULL" CCD="$(CCS)" SCOTCHLIBS="$(LIBSCOTCH)" check_prog_full

ptcheck				:	realclean
					$(MAKE) CFLAGS="$(CFLAGS) -DSCOTCH_PTSCOTCH" CC="$(CCP)" SCOTCHLIBS="$(LIBPTSCOTCH)" $(CHECKPTSCOTCH)

ptcheck_full			:	realclean
					$(MAKE) CFLAGS="$(CFLAGS) -DSCOTCH_PTSCOTCH" CC="$(CCP)" SCOTCHLIBS="$(LIBPTSCOTCH)" check_prog_dfull

escheck				:	realclean
					$(MAKE) CFLAGS="$(CFLAGS)" CC="$(CCP)" SCOTCHLIBS="$(LIBESMUMPS)" $(CHECKESMUMPS)

clean				:
					-$(RM) *~ *$(OBJ)

realclean			:	clean
					-$(RM) 				\
					test_common_file_compress	\
					test_common_random		\
					test_common_thread		\
					test_fibo			\
					test_libesmumps			\
					test_libmetis_v3		\
					test_libmetis_v5		\
					test_libmetis_dual_v3		\
					test_libmetis_dual_v5		\
					test_libmetis_dual_f.f90	\
					test_libmetis_dual_f.f90.tmp	\
					test_libmetis_dual_f_v3		\
					test_libmetis_dual_f_v5		\
					test_scotch_arch		\
					test_scotch_arch_deco		\
					test_scotch_context		\
					test_scotch_graph_coarsen	\
					test_scotch_graph_color		\
					test_scotch_graph_diam		\
					test_scotch_graph_dump2		\
					test_scotch_graph_dump2.c	\
					test_scotch_graph_induce	\
					test_scotch_graph_map		\
					test_scotch_graph_map_copy	\
					test_scotch_graph_order		\
					test_scotch_graph_part_ovl	\
					test_scotch_dgraph_band		\
					test_scotch_dgraph_check	\
					test_scotch_dgraph_coarsen	\
					test_scotch_dgraph_grow		\
					test_scotch_dgraph_induce	\
					test_scotch_dgraph_redist	\
					test_scotch_mesh_graph		\
					test_strat_par			\
					test_strat_seq

##
##  Test cases dependencies.
##

CHECKSCOTCH = 		check_common_file_compress		\
			check_common_random			\
			check_common_thread			\
			check_fibo				\
			check_strat_seq				\
			check_scotch_arch			\
			check_scotch_arch_deco			\
			check_scotch_context			\
			check_scotch_graph_coarsen		\
			check_scotch_graph_color		\
			check_scotch_graph_diam			\
			check_scotch_graph_dump			\
			check_scotch_graph_induce		\
			check_scotch_graph_map			\
			check_scotch_graph_map_copy		\
			check_scotch_graph_order		\
			check_scotch_graph_part_ovl		\
			check_scotch_mesh_graph			\
			check_libmetis_v3			\
			check_libmetis_v5			\
			check_libmetis_dual_v3			\
			check_libmetis_dual_v5			\
			check_libmetis_dual_f_v3		\
			check_libmetis_dual_f_v5		\
			check_prog_amk				\
			check_prog_amk_grf			\
			check_prog_acpl				\
			check_prog_gbase			\
			check_prog_gmap				\
			check_prog_gmk				\
			check_prog_gmk_msh			\
			check_prog_gord				\
			check_prog_gord_block			\
			check_prog_gotst			\
			check_prog_gout				\
			check_prog_gpart			\
			check_prog_gpart_clustering		\
			check_prog_gpart_overlap		\
			check_prog_gpart_remap			\
			check_prog_gscat			\
			check_prog_gtst				\
			check_prog_mmk				\
			check_prog_mord

CHECKPTSCOTCH = 	check_strat_par				\
			check_scotch_dgraph_check		\
			check_scotch_dgraph_band		\
			check_scotch_dgraph_coarsen		\
			check_scotch_dgraph_grow		\
			check_scotch_dgraph_induce      	\
			check_scotch_dgraph_redist		\
			check_prog_dgord			\
			check_prog_dgpart			\
			check_prog_dgscat-dggath		\
			check_prog_dgtst

CHECKESMUMPS = 		check_libesmumps

##
##  Sequential library tests.
##

check_common_file_compress	:	test_common_file_compress
					-cp "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf"
					$(EXECS) ./test_common_file_compress "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf.gz"
					$(EXECS) ./test_common_file_compress "$(TMPDIR)/bump_b1.grf.gz" "$(TMPDIR)/bump_b1.grf"
					$(PROGDIFF) "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf"
					-cp "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf"
					$(EXECS) ./test_common_file_compress "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf.bz2"
					$(EXECS) ./test_common_file_compress "$(TMPDIR)/bump_b1.grf.bz2" "$(TMPDIR)/bump_b1.grf"
					$(PROGDIFF) "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf"
					-cp "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf"
					$(EXECS) ./test_common_file_compress "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf.lzma"
					$(EXECS) ./test_common_file_compress "$(TMPDIR)/bump_b1.grf.lzma" "$(TMPDIR)/bump_b1.grf"
					$(PROGDIFF) "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf"

test_common_file_compress	:	test_common_file_compress.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_common_random		:	test_common_random
					$(EXECS) ./test_common_random "$(TMPDIR)/rand.dat" 0
					$(EXECS) ./test_common_random "$(TMPDIR)/rand.dat" 1

test_common_random		:	test_common_random.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_common_thread		:	test_common_thread
					$(EXECS) ./test_common_thread

test_common_thread		:	test_common_thread.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_fibo			:	test_fibo
					$(EXECS) ./test_fibo

test_fibo			:	test_fibo.c			\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_libesmumps		:	test_libesmumps
					$(EXECS) ./test_libesmumps data/bump.grf
					$(EXECS) ./test_libesmumps data/bump_imbal_32.grf

test_libesmumps			:	test_libesmumps.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)	\
					$(SCOTCHLIBDIR)/libesmumps$(LIB)

##

check_libmetis_v3		:	test_libmetis_v3
					$(EXECS) ./test_libmetis_v3 data/bump.grf
					$(EXECS) ./test_libmetis_v3 data/bump_b1.grf
					$(EXECS) ./test_libmetis_v3 data/bump_imbal_32.grf

test_libmetis_v3		:	test_libmetis.c			\
					$(SCOTCHLIBDIR)/libscotch$(LIB)	\
					$(SCOTCHLIBDIR)/libscotchmetisv3$(LIB)

check_libmetis_v5		:	test_libmetis_v5
					$(EXECS) ./test_libmetis_v5 data/bump.grf
					$(EXECS) ./test_libmetis_v5 data/bump_b1.grf
					$(EXECS) ./test_libmetis_v5 data/bump_imbal_32.grf

test_libmetis_v5		:	test_libmetis.c			\
					$(SCOTCHLIBDIR)/libscotch$(LIB)	\
					$(SCOTCHLIBDIR)/libscotchmetisv5$(LIB)

##

check_libmetis_dual_v3		:	test_libmetis_dual_v3
					$(EXECS) ./test_libmetis_dual_v3
					$(EXECS) ./test_libmetis_dual_v3 1

test_libmetis_dual_v3		:	test_libmetis_dual.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)	\
					$(SCOTCHLIBDIR)/libscotchmetisv3$(LIB)

check_libmetis_dual_v5		:	test_libmetis_dual_v5
					$(EXECS) ./test_libmetis_dual_v5
					$(EXECS) ./test_libmetis_dual_v5 1

test_libmetis_dual_v5		:	test_libmetis_dual.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)	\
					$(SCOTCHLIBDIR)/libscotchmetisv5$(LIB)

##

check_libmetis_dual_f_v3	:	test_libmetis_dual_f_v3
					$(EXECS) ./test_libmetis_dual_f_v3

test_libmetis_dual_f_v3		:	test_libmetis_dual_f.f90	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)	\
					$(SCOTCHLIBDIR)/libscotchmetisv3$(LIB)

check_libmetis_dual_f_v5	:	test_libmetis_dual_f_v5
					$(EXECS) ./test_libmetis_dual_f_v5

test_libmetis_dual_f_v5		:	test_libmetis_dual_f.f90	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)	\
					$(SCOTCHLIBDIR)/libscotchmetisv5$(LIB)

##

check_scotch_arch		:	test_scotch_arch
					$(EXECS) ./test_scotch_arch "$(TMPDIR)/arch.tgt"

test_scotch_arch		:	test_scotch_arch.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_arch_deco		:	test_scotch_arch_deco
					$(EXECS) ./test_scotch_arch_deco data/m4x4.grf "$(TMPDIR)/arch2.tgt"
#					$(EXECS) ./test_scotch_arch_deco data/m4x4_b100000.grf "$(TMPDIR)/arch2.tgt"

test_scotch_arch_deco		:	test_scotch_arch_deco.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_context		:	test_scotch_context
					$(EXECS) ./test_scotch_context data/bump_b1.grf

test_scotch_context		:	test_scotch_context.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_coarsen	:	test_scotch_graph_coarsen
					$(EXECS) ./test_scotch_graph_coarsen data/bump.grf
					$(EXECS) ./test_scotch_graph_coarsen data/bump_b100000.grf

test_scotch_graph_coarsen	:	test_scotch_graph_coarsen.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_color	:	test_scotch_graph_color
					$(EXECS) ./test_scotch_graph_color data/bump.grf
					$(EXECS) ./test_scotch_graph_color data/bump_b100000.grf

test_scotch_graph_color		:	test_scotch_graph_color.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_diam		:	test_scotch_graph_diam
					$(EXECS) ./test_scotch_graph_diam data/bump_b100000.grf
					$(EXECS) ./test_scotch_graph_diam data/m16x16_b100000.grf
					$(EXECS) ./test_scotch_graph_diam data/p2-p2.grf

test_scotch_graph_diam		:	test_scotch_graph_diam.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_dump		:	test_scotch_graph_dump2
					$(EXECS) ./test_scotch_graph_dump2 $(TMPDIR)/m16x16_b1.grf
					$(PROGDIFF) "data/bump_b1.grf" "$(TMPDIR)/bump_b1.grf"

test_scotch_graph_dump2		:	test_scotch_graph_dump2.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

test_scotch_graph_dump2.c	:	test_scotch_graph_dump.c
					$(EXECS) $(SCOTCHBINDIR)/gdump -ptest data/m16x16_b1.grf $(TMPDIR)/m16x16_b1.c
					$(CAT) test_scotch_graph_dump.c $(TMPDIR)/m16x16_b1.c > $(@)


##

check_scotch_graph_induce	:	test_scotch_graph_induce
					$(EXECS) ./test_scotch_graph_induce data/bump.grf
					$(EXECS) ./test_scotch_graph_induce data/bump_b100000.grf

test_scotch_graph_induce	:	test_scotch_graph_induce.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_map		:	test_scotch_graph_map
					$(EXECS) ./test_scotch_graph_map data/m4x4.grf
					$(EXECS) ./test_scotch_graph_map data/m4x4_b100000.grf
					$(EXECS) ./test_scotch_graph_map data/m16x16.grf
					$(EXECS) ./test_scotch_graph_map data/m16x16_b100000.grf

test_scotch_graph_map		:	test_scotch_graph_map.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_map_copy	:	test_scotch_graph_map_copy
					$(EXECS) ./test_scotch_graph_map_copy data/bump.grf
					$(EXECS) ./test_scotch_graph_map_copy data/bump_b100000.grf

test_scotch_graph_map_copy	:	test_scotch_graph_map_copy.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_order	:	test_scotch_graph_order
					$(EXECS) ./test_scotch_graph_order data/bump.grf
					$(EXECS) ./test_scotch_graph_order data/bump_b100000.grf

test_scotch_graph_order		:	test_scotch_graph_order.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_part_ovl	:	test_scotch_graph_part_ovl
					$(EXECS) ./test_scotch_graph_part_ovl 4 data/m16x16.grf /dev/null
					$(EXECS) ./test_scotch_graph_part_ovl 4 data/m16x16_b100000.grf /dev/null

test_scotch_graph_part_ovl	:	test_scotch_graph_part_ovl.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_mesh_graph		:	test_scotch_mesh_graph
					$(EXECS) ./test_scotch_mesh_graph data/ship001.msh

test_scotch_mesh_graph		:	test_scotch_mesh_graph.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_strat_seq			:	test_strat_seq
					$(EXECS) ./test_strat_seq

test_strat_seq			:	test_strat_seq.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##
##  Parallel library tests.
##

check_strat_par			:	test_strat_par
					$(EXECS) ./test_strat_par

test_strat_par			:	test_strat_par.c		\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_band	:	test_scotch_dgraph_band
					$(EXECP3) ./test_scotch_dgraph_band data/bump.grf $(TMPDIR)/test_scotch_dgraph_band.map
					$(EXECP3) ./test_scotch_dgraph_band data/bump_b100000.grf $(TMPDIR)/test_scotch_dgraph_band.map

test_scotch_dgraph_band		:	test_scotch_dgraph_band.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_coarsen	:	test_scotch_dgraph_coarsen
					$(EXECP3) ./test_scotch_dgraph_coarsen data/bump.grf
					$(EXECP3) ./test_scotch_dgraph_coarsen data/bump_b100000.grf
					$(EXECP3) ./test_scotch_dgraph_coarsen data/m4x4_b1.grf

test_scotch_dgraph_coarsen	:	test_scotch_dgraph_coarsen.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_check	:	test_scotch_dgraph_check
					$(EXECP3) ./test_scotch_dgraph_check data/bump.grf
					$(EXECP3) ./test_scotch_dgraph_check data/bump_b100000.grf

test_scotch_dgraph_check	:	test_scotch_dgraph_check.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_grow	:	test_scotch_dgraph_grow
					$(EXECP3) ./test_scotch_dgraph_grow data/bump.grf $(TMPDIR)/test_scotch_dgraph_grow.map
					$(EXECP3) ./test_scotch_dgraph_grow data/bump_b100000.grf $(TMPDIR)/test_scotch_dgraph_grow.map

test_scotch_dgraph_grow		:	test_scotch_dgraph_grow.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_induce	:	test_scotch_dgraph_induce
					$(EXECP3) ./test_scotch_dgraph_induce data/bump.grf
					$(EXECP3) ./test_scotch_dgraph_induce data/bump_b100000.grf

test_scotch_dgraph_induce	:	test_scotch_dgraph_induce.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_redist	:	test_scotch_dgraph_redist
					$(EXECP3) ./test_scotch_dgraph_redist data/bump.grf
					$(EXECP3) ./test_scotch_dgraph_redist data/bump_b100000.grf

test_scotch_dgraph_redist	:	test_scotch_dgraph_redist.c

##
##  Sequential program tests.
##

check_prog_acpl			:	check_prog_amk
					$(EXECS) $(SCOTCHBINDIR)/amk_m2 5 12 | $(SCOTCHBINDIR)/acpl | $(SCOTCHBINDIR)/atst

check_prog_amk			:	check_prog_atst
					$(EXECS) $(SCOTCHBINDIR)/amk_ccc 5 | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_fft2 5 | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_hy 5 | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_m2 5 12 | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_m2 -mo 5 12 | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_p2 5 7 | $(SCOTCHBINDIR)/atst

check_prog_amk_grf		:	check_prog_atst
					$(EXECS) $(SCOTCHBINDIR)/amk_grf data/m16x16.grf | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_grf -2 data/m16x16.grf | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_grf -ldata/m4x4_vertlist.txt data/m4x4.grf | $(SCOTCHBINDIR)/atst
					$(EXECS) $(SCOTCHBINDIR)/amk_grf -ldata/m4x4_vertlist.txt -2 data/m4x4.grf | $(SCOTCHBINDIR)/atst

check_prog_atst			:
					$(EXECS) $(SCOTCHBINDIR)/atst $(SCOTCHTGTDIR)/t4x4x4.tgt

check_prog_gbase		:
					$(EXECS) $(SCOTCHBINDIR)/gbase 1 data/m4x4.grf $(TMPDIR)/m4x4_b1.grf
					$(PROGDIFF) "data/m4x4_b1.grf" "$(TMPDIR)/m4x4_b1.grf"
					$(EXECS) $(SCOTCHBINDIR)/gbase 1 data/m4x4_b1.grf $(TMPDIR)/m4x4_b1.grf
					$(PROGDIFF) "data/m4x4_b1.grf" "$(TMPDIR)/m4x4_b1.grf"
					$(EXECS) $(SCOTCHBINDIR)/gbase 0 data/m4x4_b1.grf $(TMPDIR)/m4x4.grf
					$(PROGDIFF) "data/m4x4.grf" "$(TMPDIR)/m4x4.grf"

check_prog_gmk			:
					$(EXECS) $(SCOTCHBINDIR)/gmk_hy 5 | $(SCOTCHBINDIR)/gtst
					$(EXECS) $(SCOTCHBINDIR)/gmk_m2 5 1 | $(SCOTCHBINDIR)/gtst
					$(EXECS) $(SCOTCHBINDIR)/gmk_m2 -b1 -e 5 4 -g/dev/null | $(SCOTCHBINDIR)/gtst
					$(EXECS) $(SCOTCHBINDIR)/gmk_m2 -t 5 4 | $(SCOTCHBINDIR)/gtst
					$(EXECS) $(SCOTCHBINDIR)/gmk_m2 -b1 -e -t 5 2 | $(SCOTCHBINDIR)/gtst
					$(EXECS) $(SCOTCHBINDIR)/gmk_m2 -b1 -e -t 5 1 | $(SCOTCHBINDIR)/gtst
					$(EXECS) $(SCOTCHBINDIR)/gmk_m3 -b1 5 4 3 -g/dev/null | $(SCOTCHBINDIR)/gtst
					$(EXECS) $(SCOTCHBINDIR)/gmk_ub2 5 | $(SCOTCHBINDIR)/gtst

check_prog_gmk_msh		:	check_prog_mmk
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh $(TMPDIR)/m2_5_1.msh $(TMPDIR)/m2_5_1.grf
					$(EXECS) $(SCOTCHBINDIR)/gtst $(TMPDIR)/m2_5_1.grf
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh $(TMPDIR)/m2_5_4.msh $(TMPDIR)/m2_5_4.grf
					$(EXECS) $(SCOTCHBINDIR)/gtst $(TMPDIR)/m2_5_4.grf
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh $(TMPDIR)/m2_5_4_3.msh $(TMPDIR)/m2_5_4_3.grf
					$(EXECS) $(SCOTCHBINDIR)/gtst $(TMPDIR)/m2_5_4_3.grf
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh data/ship001.msh $(TMPDIR)/ship001.grf
					$(EXECS) $(SCOTCHBINDIR)/gtst $(TMPDIR)/ship001.grf
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh "data/small_55.msh" $(TMPDIR)/small_55_d1.grf -d1
					$(PROGDIFF) "data/small_55_d1.grf" $(TMPDIR)/small_55_d1.grf
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh "data/small_55.msh" $(TMPDIR)/small_55_d2.grf -d2
					$(PROGDIFF) "data/small_55_d2.grf" $(TMPDIR)/small_55_d2.grf
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh "data/small_55.msh" $(TMPDIR)/small_55_d3.grf -d3
					$(PROGDIFF) "data/small_55_d2.grf" $(TMPDIR)/small_55_d3.grf
					$(EXECS) $(SCOTCHBINDIR)/gmk_msh "data/cube_8.msh" $(TMPDIR)/cube_8_d3.grf -d3
					$(PROGDIFF) "data/cube_8_d3.grf" $(TMPDIR)/cube_8_d3.grf

check_prog_mmk                  :
					$(EXECS) $(SCOTCHBINDIR)/mmk_m2 5 1 $(TMPDIR)/m2_5_1.msh
					$(EXECS) $(SCOTCHBINDIR)/mmk_m2 5 4 $(TMPDIR)/m2_5_4.msh
					$(EXECS) $(SCOTCHBINDIR)/mmk_m3 5 4 3 $(TMPDIR)/m2_5_4_3.msh
					$(EXECS) $(SCOTCHBINDIR)/mtst $(TMPDIR)/m2_5_1.msh
					$(EXECS) $(SCOTCHBINDIR)/mtst $(TMPDIR)/m2_5_4.msh
					$(EXECS) $(SCOTCHBINDIR)/mtst $(TMPDIR)/m2_5_4_3.msh

check_prog_mcv		        :	check_prog_mmk
					$(EXECS) $(SCOTCHBINDIR)/mcv $(TMPDIR)/m2_5_1.msh $(TMPDIR)/output_m2_5_1.msh -is -os
					$(EXECS) $(SCOTCHBINDIR)/mcv $(TMPDIR)/m2_5_4.msh $(TMPDIR)/output_m2_5_4.msh -is -os
					$(EXECS) $(SCOTCHBINDIR)/mcv $(TMPDIR)/m2_5_4_3.msh $(TMPDIR)/output_m2_5_4_3.msh -is -os
					$(PROGDIFF) $(TMPDIR)/m2_5_1.msh $(TMPDIR)/output_m2_5_1.msh
					$(PROGDIFF) $(TMPDIR)/m2_5_4.msh $(TMPDIR)/output_m2_5_4.msh
					$(PROGDIFF) $(TMPDIR)/m2_5_4_3.msh $(TMPDIR)/output_m2_5_4_3.msh


check_prog_mord                 :	check_prog_gmk_msh
					$(EXECS) $(SCOTCHBINDIR)/mord $(TMPDIR)/m2_5_1.msh | $(SCOTCHBINDIR)/gotst $(TMPDIR)/m2_5_1.grf
					$(EXECS) $(SCOTCHBINDIR)/mord $(TMPDIR)/m2_5_4.msh | $(SCOTCHBINDIR)/gotst $(TMPDIR)/m2_5_4.grf
					$(EXECS) $(SCOTCHBINDIR)/mord $(TMPDIR)/m2_5_4_3.msh | $(SCOTCHBINDIR)/gotst $(TMPDIR)/m2_5_4_3.grf
					$(EXECS) $(SCOTCHBINDIR)/mord data/ship001.msh | $(SCOTCHBINDIR)/gotst $(TMPDIR)/ship001.grf

check_prog_gmap			:
					$(EXECS) $(SCOTCHBINDIR)/gmap data/bump.grf $(SCOTCHTGTDIR)/h3.tgt $(TMPDIR)/bump_h3.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gmap data/small.grf $(SCOTCHTGTDIR)/m11x13.tgt $(TMPDIR)/small_m11x13.map -vmt

check_prog_gord			:	check_prog_gmk_msh
					$(EXECS) $(SCOTCHBINDIR)/gord data/bump.grf $(TMPDIR)/bump.ord -Cu -vt
					$(EXECS) $(SCOTCHBINDIR)/gord data/bump_b1.grf $(TMPDIR)/bump_b1.ord -Cd -vt
					$(EXECS) $(SCOTCHBINDIR)/gord data/bump_b1.grf '-On{sep=(/((vert)>(240))?(e{strat=g}|g|h);),ole=d{cmin=15,cmax=100000,frat=0},ose=g{pass=3}}' $(TMPDIR)/bump_b1_2.ord -Cr -vt
					$(EXECS) $(SCOTCHBINDIR)/gord data/bump_imbal_32.grf $(TMPDIR)/bump_imbal_32.ord -Cf -vt
					$(EXECS) $(SCOTCHBINDIR)/gord data/nocoarsen.grf $(TMPDIR)/nocoarsen.ord -vt
					$(EXECS) $(SCOTCHBINDIR)/gord data/compress-coarsen.grf $(TMPDIR)/compress-coarsen.ord -vt
					$(EXECS) $(SCOTCHBINDIR)/gord $(TMPDIR)/ship001.grf $(TMPDIR)/ship001.ord -vt

check_prog_gotst		:	check_prog_gord
					$(EXECS) $(SCOTCHBINDIR)/gotst data/bump.grf $(TMPDIR)/bump.ord
					$(EXECS) $(SCOTCHBINDIR)/gotst data/bump_b1.grf $(TMPDIR)/bump_b1.ord
					$(EXECS) $(SCOTCHBINDIR)/gotst data/bump_imbal_32.grf $(TMPDIR)/bump_imbal_32.ord
					$(EXECS) $(SCOTCHBINDIR)/gotst data/nocoarsen.grf $(TMPDIR)/nocoarsen.ord
					$(EXECS) $(SCOTCHBINDIR)/gotst data/compress-coarsen.grf $(TMPDIR)/compress-coarsen.ord
					$(EXECS) $(SCOTCHBINDIR)/gotst $(TMPDIR)/ship001.grf $(TMPDIR)/ship001.ord

check_prog_gord_block		:	check_prog_gord
					$(EXECS) $(SCOTCHBINDIR)/gord data/m4x4_b1_ev.grf $(TMPDIR)/m4x4_b1.ord '-ob{cmin=3,strat=s}'
					$(EXECS) $(SCOTCHBINDIR)/gord data/m4x4_b1_ev.grf $(TMPDIR)/m4x4_b1_ev.ord '-ob{cmin=2,strat=s}'
					$(EXECS) $(SCOTCHBINDIR)/gord data/m4x4_b1_ev.grf $(TMPDIR)/m4x4_b1_ev.ord '-ob{cmin=1,strat=s}'

check_prog_gout			:	check_prog_gord		\
					check_prog_gpart
					$(EXECS) $(SCOTCHBINDIR)/gout '-Oi{c,r}' data/bump.grf data/bump.xyz $(TMPDIR)/bump_k9.map $(TMPDIR)/bump_k9.iv
					$(EXECS) $(SCOTCHBINDIR)/gout '-Om{e}' -gn data/bump.grf /dev/null $(TMPDIR)/bump.ord $(TMPDIR)/bump_ord.eps
					$(EXECS) $(SCOTCHBINDIR)/gout '-Op{d,e,g,r,s}' data/bump.grf data/bump.xyz $(TMPDIR)/bump_k9.map $(TMPDIR)/bump_k9.eps
					$(EXECS) $(SCOTCHBINDIR)/gout '-Ot{b,c,r}' -mn data/bump.grf data/bump.xyz /dev/null $(TMPDIR)/bump_k9.tlp
					$(EXECS) $(SCOTCHBINDIR)/gout '-Ov{v}' data/bump.grf data/bump.xyz $(TMPDIR)/bump_k9.map $(TMPDIR)/bump_k9.vtk

check_prog_gpart		:
					$(EXECS) $(SCOTCHBINDIR)/gpart 9 data/bump.grf $(TMPDIR)/bump_k9.map -Cd -vmt
					$(EXECS) $(SCOTCHBINDIR)/gpart 9 data/bump_b100000.grf $(TMPDIR)/bump_k9.map -Cu -vmt
					$(EXECS) $(SCOTCHBINDIR)/gpart 9 data/nocoarsen.grf $(TMPDIR)/nocoarsen_k9.map -Cr -vmt
					$(EXECS) $(SCOTCHBINDIR)/gpart 5 data/bump_b100000.grf $(TMPDIR)/bump_k5.map -Cf -vmt '-Mr{bal=0.05,sep=m{vert=120,low=h{pass=10}f{bal=0.05,move=120},asc=b{bnd=f{bal=0.05,move=120},org=f{bal=0.05,move=120}}}}|m{vert=1000,low=r{bal=0.05,sep=m{vert=120,low=h{pass=10}f{bal=0.05,move=120},asc=b{bnd=f{bal=0.05,move=120},org=f{bal=0.05,move=120}}}},asc=b{bnd=f{bal=0.05,move=120},org=f{bal=0.05,move=120}}}'

check_prog_gpart_clustering	:
					$(EXECS) $(SCOTCHBINDIR)/gpart -q 1 data/bump.grf $(TMPDIR)/bump_part_cls_9.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gpart -q 1 data/bump_b100000.grf $(TMPDIR)/bump_b100000_part_cls_9.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gmtst data/bump_b100000.grf $(SCOTCHTGTDIR)/vcmplt.tgt $(TMPDIR)/bump_b100000_part_cls_9.map
					$(EXECS) $(SCOTCHBINDIR)/gpart -q 1 data/bump_imbal_32.grf $(TMPDIR)/bump_imbal_32_part_cls_32.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gmtst data/bump_imbal_32.grf $(SCOTCHTGTDIR)/vhcub.tgt $(TMPDIR)/bump_imbal_32_part_cls_32.map

check_prog_gpart_overlap	:
					$(EXECS) $(SCOTCHBINDIR)/gpart -o 9 data/bump.grf $(TMPDIR)/bump_part_ovl_9.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gpart -o 9 data/bump_b100000.grf $(TMPDIR)/bump_b100000_part_ovl_9.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gmtst -o data/bump_b100000.grf $(SCOTCHTGTDIR)/k9.tgt $(TMPDIR)/bump_b100000_part_ovl_9.map
					$(EXECS) $(SCOTCHBINDIR)/gpart -o 32 data/bump_imbal_32.grf $(TMPDIR)/bump_imbal_32_part_ovl_32.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gmtst -o data/bump_imbal_32.grf $(SCOTCHTGTDIR)/k32.tgt $(TMPDIR)/bump_imbal_32_part_ovl_32.map
					$(EXECS) $(SCOTCHBINDIR)/gpart -o 6 data/bump_imbal_32.grf $(TMPDIR)/bump_imbal_32_part_ovl_32_6.map '-me{strat=m{asc=b{width=3,bnd=d{pass=40,dif=1,rem=0}f{move=80,pass=-1,bal=0.01},org=f{move=80,pass=-1,bal=0.01}},low=r{job=t,bal=0.01,map=t,poli=S,sep=(m{asc=b{bnd=f{move=120,pass=-1,bal=0.01,type=b},org=f{move=120,pass=-1,bal=0.01,type=b},width=3},low=h{pass=10}f{move=120,pass=-1,bal=0.01,type=b},vert=120,rat=0.8}|m{asc=b{bnd=f{move=120,pass=-1,bal=0.01,type=b},org=f{move=120,pass=-1,bal=0.01,type=b},width=3},low=h{pass=10}f{move=120,pass=-1,bal=0.01,type=b},vert=120,rat=0.8})},vert=10000,rat=0.8,type=h}}f' -vmt
					$(EXECS) $(SCOTCHBINDIR)/gpart -o 24 data/m4x4.grf $(TMPDIR)/bump_m4x4_part_ovl_24.map -vmt
					$(EXECS) $(SCOTCHBINDIR)/gmtst -o data/m4x4.grf $(SCOTCHTGTDIR)/k24.tgt $(TMPDIR)/bump_m4x4_part_ovl_24.map

check_prog_gpart_remap		:
					$(EXECS) $(SCOTCHBINDIR)/gpart 32 data/bump_imbal_32.grf /dev/null -rodata/bump_old.map -vmt

check_prog_gscat		:
					$(EXECS) $(SCOTCHBINDIR)/gscat 4 data/small.grf $(TMPDIR)/small_%r.grf
					$(EXECS) $(SCOTCHBINDIR)/gscat 4 -i1 data/bump.grf $(TMPDIR)/bump4_%r.grf

check_prog_gtst			:
					$(EXECS) $(SCOTCHBINDIR)/gtst data/m4x4_b1_ev.grf $(TMPDIR)/gtst_ev.txt
					$(EXECS) $(SCOTCHBINDIR)/gtst data/m4x4_b1_elv.grf $(TMPDIR)/gtst_elv.txt
					$(PROGDIFF) "$(TMPDIR)/gtst_ev.txt" "$(TMPDIR)/gtst_elv.txt"

##

check_prog_full		:
					$(EXECS) $(SCOTCHBINDIR)/gord data/bump.grf $(TMPDIR)/bump_f.ord -vt
					$(EXECS) $(SCOTCHBINDIR)/gpart 3 data/bump.grf $(TMPDIR)/bump_f.map -vt
					$(EXECS) $(SCOTCHBINDIR)/gpart -o 3 data/bump.grf $(TMPDIR)/bump_fo.map -vt

##
##  Parallel program tests.
##

check_prog_dgord		:
					$(EXECP3) $(SCOTCHBINDIR)/dgord data/bump.grf /dev/null -Cu -vt
					$(EXECP3) $(SCOTCHBINDIR)/dgord data/bump_b100000.grf /dev/null -Cd -vt
					$(EXECP4) $(SCOTCHBINDIR)/dgord data/small_%r.grf $(TMPDIR)/small.ord -vt

check_prog_dgpart		:
					$(EXECP3) $(SCOTCHBINDIR)/dgpart 9 data/bump.grf /dev/null -Cu -vmt
					$(EXECP3) $(SCOTCHBINDIR)/dgpart 9 data/bump_b100000.grf /dev/null -Cd -vmt
					$(EXECP4) $(SCOTCHBINDIR)/dgpart 50 data/small_%r.grf /dev/null -vmt

check_prog_dgscat-dggath	:	TMP := $(shell mktemp -u)
check_prog_dgscat-dggath	:
					$(EXECP3) $(SCOTCHBINDIR)/dgscat -r1 -c data/m16x16.grf "$(TMP)1_%r.grf"
					$(EXECP3) $(SCOTCHBINDIR)/dggath -r1 -c "$(TMP)1_%r.grf" "$(TMP)1.grf"
					$(PROGDIFF) data/m16x16.grf "$(TMP)1.grf"
					$(EXECP3) $(SCOTCHBINDIR)/dgscat -c data/m16x16_b100000.grf "$(TMP)2_%r.grf"
					$(EXECP3) $(SCOTCHBINDIR)/dggath -c "$(TMP)2_%r.grf" "$(TMP)2.grf"
					$(PROGDIFF) data/m16x16_b100000.grf "$(TMP)2.grf"

check_prog_dgtst		:
					$(EXECP4) $(SCOTCHBINDIR)/dgtst data/small_%r.grf
##

check_prog_dfull		:
					$(EXECP4) $(SCOTCHBINDIR)/dgord data/bump.grf $(TMPDIR)/bump_f.ord -vt
					$(EXECP4) $(SCOTCHBINDIR)/dgpart 3 data/bump.grf $(TMPDIR)/bump_f.map -vmt
