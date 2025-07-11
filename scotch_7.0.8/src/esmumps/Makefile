## Copyright 2004,2007-2010,2012,2015,2020,2022 IPB, Universite de Bordeaux, INRIA & CNRS
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
srclibscotchdir	= ../libscotch

##
##  General inference rules.
##

include ../Makefile.inc

%$(OBJ)	:	%.c
				$(CC) $(CFLAGS) $(CLIBFLAGS) -I$(srclibscotchdir) -I$(includedir) -c $(<) -o $(@)

%$(EXE)	:	%.c
				$(CC) $(CFLAGS) -I$(srclibscotchdir) -I$(includedir) $(<) -o $(@) -L. -lesmumps -L$(libdir) -lscotch -lscotcherrexit $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	ptscotch	scotch	ptinstall	install	clean	realclean

scotch				:
					$(MAKE) CC="$(CCS)" CCD="$(CCS)"	\
					libesmumps$(LIB)			\
					test_esmumps$(EXE)			\
					test_fax$(EXE)

ptscotch			:	scotch

install				:	scotch					\
					library.h
					-$(CP) library.h $(includedir)/esmumps.h
					-$(CP) libesmumps$(LIB) $(libdir)

ptinstall			:	scotch					\
					library.h
					-$(CP) library.h $(includedir)/esmumps.h
					-$(CP) libesmumps$(LIB) $(libdir)/libptesmumps$(LIB)

clean				:
					-$(RM) *~ common.h *$(OBJ) lib*$(LIB) test_esmumps$(EXE) test_fax$(EXE)

realclean			:	clean

##
##  Library dependencies.
##

LIBESMUMPSDEPS =			dof$(OBJ)				\
					esmumps$(OBJ)				\
					graph_graph$(OBJ)			\
					library_esmumps$(OBJ)			\
					library_esmumps_f$(OBJ)			\
					library_esmumps_strats$(OBJ)		\
					order$(OBJ)				\
					order_check$(OBJ)			\
					order_io$(OBJ)				\
					order_scotch_graph$(OBJ)		\
					symbol$(OBJ)				\
					symbol_check$(OBJ)			\
					symbol_cost$(OBJ)			\
					symbol_fax_graph$(OBJ)

##
##  Todo list.
##

module.h			:	$(srclibscotchdir)/module.h

dof$(OBJ)			:	dof.c					\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(includedir)/scotch.h			\
					graph.h					\
					dof.h

esmumps$(OBJ)			:	esmumps.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(includedir)/scotch.h			\
					graph.h					\
					dof.h					\
					symbol.h				\
					order.h					\
					fax.h					\
					esmumps.h

graph_graph$(OBJ)		:	graph_graph.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(includedir)/scotch.h			\
					graph.h					\

library_esmumps$(OBJ)		:	library_esmumps.c			\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(srclibscotchdir)/scotch.h		\
					esmumps.h				\
					library.h

library_esmumps_f$(OBJ)		:	library_esmumps_f.c			\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(includedir)/scotch.h			\
					library.h

library_esmumps_strats$(OBJ)	:	library_esmumps_strats.c		\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(includedir)/scotch.h			\
					library.h

order$(OBJ)			:	order.c					\
					module.h				\
					$(srclibscotchdir)/common.h		\
					graph.h					\
					order.h

order_check$(OBJ)		:	order_check.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					graph.h					\
					order.h

order_scotch_graph$(OBJ)	:	order_scotch_graph.c			\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(includedir)/scotch.h			\
					graph.h					\
					order.h

order_io$(OBJ)			:	order_io.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					graph.h					\
					order.h

symbol$(OBJ)			:	symbol.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					symbol.h

symbol_check$(OBJ)		:	symbol_check.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					symbol.h

symbol_cost$(OBJ)		:	symbol_cost.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					symbol.h				\
					symbol_cost.h

symbol_fax$(OBJ)		:	symbol_fax.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					graph.h					\
					symbol.h				\
					order.h					\
					fax.h

symbol_fax_graph$(OBJ)		:	symbol_fax_graph.c			\
					symbol_fax.c				\
					module.h				\
					$(srclibscotchdir)/common.h		\
					$(includedir)/scotch.h			\
					graph.h					\
					symbol.h				\
					order.h					\
					fax.h

libesmumps$(LIB)		:	$(LIBESMUMPSDEPS)
					$(AR) $(ARFLAGS) $(@) $(?)
					-$(RANLIB) $(@)

test_esmumps$(EXE)		:	test_esmumps.c				\
					$(srclibscotchdir)/common.h		\
					module.h				\
					graph.h					\
					order.h					\
					symbol.h				\
					esmumps.h				\
					libesmumps$(LIB)			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)

test_fax$(EXE)			:	test_fax.c				\
					$(srclibscotchdir)/common.h		\
					module.h				\
					graph.h					\
					order.h					\
					dof.h					\
					symbol.h				\
					fax.h					\
					esmumps.h				\
					libesmumps$(LIB)			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)
