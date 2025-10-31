## Copyright 2021 IPB, Universite de Bordeaux, INRIA & CNRS
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

rootdir		= ../..
includedir	= $(rootdir)/include
libdir		= $(rootdir)/lib

include $(rootdir)/src/Makefile.inc

##
##  Project rules.
##

.PHONY				:	clean	default	all

default				:	all

required			:
					(cd ../../src; make scotch)

all				:	required		\
					test_d2color_seq$(EXE)	\
					test_d2color_shm$(EXE)	\
					test_d2match_seq$(EXE)	\
					test_d2match_shm$(EXE)

clean				:
					-$(RM) test_d2color*$(OBJ) test_d2color_seq$(EXE) test_d2color_shm$(EXE) test_d2match_seq$(EXE) test_d2match_shm$(EXE)

test_d2color_seq$(EXE)		:	test_d2color_seq.c
					$(CC) $(CFLAGS) -Wall -g $(<) $(hamfobj) -I$(includedir) -I$(rootdir)/src/libscotch -L$(libdir) -lscotch -lscotcherr -lm -pthread -o $(@)

test_d2color_shm$(EXE)		:	test_d2color_shm.c
					$(CC) $(CFLAGS) -Wall -g $(<) $(hamfobj) -I$(includedir) -I$(rootdir)/src/libscotch -L$(libdir) -lscotch -lscotcherr -lm -pthread -o $(@)

test_d2match_seq$(EXE)		:	test_d2match_seq.c
					$(CC) $(CFLAGS) -Wall -g $(<) $(hamfobj) -I$(includedir) -I$(rootdir)/src/libscotch -L$(libdir) -lscotch -lscotcherr -lm -pthread -o $(@)

test_d2match_shm$(EXE)		:	test_d2match_shm.c
					$(CC) $(CFLAGS) -Wall -g $(<) $(hamfobj) -I$(includedir) -I$(rootdir)/src/libscotch -L$(libdir) -lscotch -lscotcherr -lm -pthread -o $(@)
