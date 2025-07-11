## Copyright 2007-2010,2018,2019,2021 IPB, Universite de Bordeaux, INRIA & CNRS
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

bindir		= ../../bin
includedir	= ../../include
libdir		= ../../lib

V		= v$(SCOTCH_METIS_VERSION)

##
##  General inference rules.
##

include ../Makefile.inc

%_$(V)$(OBJ):	%.c
				$(CC) $(CFLAGS) $(CLIBFLAGS) -I$(includedir) -I../libscotch -c $(<) -o $(@)

%$(EXE)	:	%.c
		 		$(CC) $(CFLAGS) -I$(includedir) -I../libscotch $(<) -o $(@) -L$(libdir) $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	ptscotch	scotch	ptinstall	install	clean	realclean

scotch				:
					$(MAKE)						\
					CC="$(CCS)"					\
					SCOTCH_METIS_VERSION=3				\
					metis.h						\
					libscotchmetisv3$(LIB)
					$(MAKE)						\
					CC="$(CCS)"					\
					SCOTCH_METIS_VERSION=5				\
					metis.h						\
					metisf.h					\
					libscotchmetisv5$(LIB)

ptscotch			:
					$(MAKE)						\
					CC="$(CCP)"					\
					CFLAGS="$(CFLAGS) -DSCOTCH_PTSCOTCH"		\
					SCOTCH_METIS_VERSION=3				\
					SCOTCH_PARMETIS_VERSION=3			\
					parmetis.h					\
					libptscotchparmetisv3$(LIB)

install				:	scotch
					-$(CP) metis.h metisf.h $(includedir)
					-$(CP) libscotchmetisv3$(LIB) libscotchmetisv5$(LIB) $(libdir)

ptinstall			:	ptscotch
					-$(CP) parmetis.h $(includedir)
					-$(CP) libptscotchparmetisv3$(LIB) $(libdir)

clean				:
					-$(RM) *~ *$(OBJ) lib*$(LIB) metis.h metisf.h parmetis.h

realclean			:	clean

##
##  Todo list.
##

metis.h				:	../libscotch/dummysizes$(EXE)		\
					library_metis.h
					../libscotch/dummysizes$(EXE) library_metis.h metis.h

metisf.h			:	../libscotch/dummysizes$(EXE)		\
					library_metis_f.h
					../libscotch/dummysizes$(EXE) library_metis_f.h metisf.h

metis_graph_dual_$(V)$(OBJ)	:	metis_graph_dual.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h					\
					metis_graph_dual.h

metis_graph_dual_f_$(V)$(OBJ)	:	metis_graph_dual_f.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_graph_order_$(V)$(OBJ)	:	metis_graph_order.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_graph_order_f_$(V)$(OBJ)	:	metis_graph_order_f.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_graph_part_$(V)$(OBJ)	:	metis_graph_part.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h					\
					metis_graph_part.h

metis_graph_part_f_$(V)$(OBJ)	:	metis_graph_part_f.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_graph_part_dual_$(V)$(OBJ):	metis_graph_part_dual.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h					\
					metis_graph_dual.h			\
					metis_graph_part.h

metis_graph_part_dual_f_$(V)$(OBJ):	metis_graph_part_dual_f.c		\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_options_$(V)$(OBJ)	:	metis_options.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_options_f_$(V)$(OBJ)	:	metis_options_f.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

parmetis.h			:	../libscotch/dummysizes$(EXE)		\
					library_parmetis.h
					../libscotch/dummysizes$(EXE) library_parmetis.h parmetis.h

parmetis_dgraph_order_$(V)$(OBJ):	parmetis_dgraph_order.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

parmetis_dgraph_order_f_$(V)$(OBJ):	parmetis_dgraph_order_f.c		\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

parmetis_dgraph_part_$(V)$(OBJ)	:	parmetis_dgraph_part.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

parmetis_dgraph_part_f_$(V)$(OBJ):	parmetis_dgraph_part_f.c		\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

libptscotchparmetis$(V)$(LIB)	:	parmetis_dgraph_order_$(V)$(OBJ)	\
					parmetis_dgraph_order_f_$(V)$(OBJ)	\
					parmetis_dgraph_part_$(V)$(OBJ)		\
					parmetis_dgraph_part_f_$(V)$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(^)
					-$(RANLIB) $(@)

libscotchmetis$(V)$(LIB)	:	metis_graph_dual_$(V)$(OBJ)		\
					metis_graph_dual_f_$(V)$(OBJ)		\
					metis_graph_order_$(V)$(OBJ)		\
					metis_graph_order_f_$(V)$(OBJ)		\
					metis_graph_part_$(V)$(OBJ)		\
					metis_graph_part_f_$(V)$(OBJ)		\
					metis_graph_part_dual_$(V)$(OBJ)	\
					metis_graph_part_dual_f_$(V)$(OBJ)	\
					metis_options_$(V)$(OBJ)		\
					metis_options_f_$(V)$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(^)
					-$(RANLIB) $(@)
