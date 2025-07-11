/* Copyright 2004,2007,2012,2018,2024 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : hgraph_order_st.h                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the data           **/
/**                declarations for the main graph         **/
/**                ordering routine.                       **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 31 oct 1996     **/
/**                                 to   : 29 aug 1998     **/
/**                # Version 3.3  : from : 02 oct 1998     **/
/**                                 to   : 07 sep 2001     **/
/**                # Version 4.0  : from : 29 dec 2001     **/
/**                                 to   : 15 jan 2003     **/
/**                # Version 6.0  : from : 17 oct 2012     **/
/**                                 to   : 05 apr 2018     **/
/**                # Version 7.0  : from : 11 sep 2024     **/
/**                                 to   : 07 nov 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The type definitions.
*/

/*+ Method types. +*/

typedef enum HgraphOrderStMethodType_ {
  HGRAPHORDERSTMETHBL = 0,                        /*+ Block splitting post-processing       +*/
  HGRAPHORDERSTMETHCC,                            /*+ Splitting of connected components     +*/
  HGRAPHORDERSTMETHCP,                            /*+ Graph compression                     +*/
  HGRAPHORDERSTMETHGP,                            /*+ Gibbs-Poole-Stockmeyer                +*/
  HGRAPHORDERSTMETHHD,                            /*+ Block Halo Approximate Minimum Degree +*/
  HGRAPHORDERSTMETHHF,                            /*+ Block Halo Approximate Minimum Fill   +*/
  HGRAPHORDERSTMETHKP,                            /*+ K-way block partitioning              +*/
  HGRAPHORDERSTMETHND,                            /*+ Nested Dissection                     +*/
  HGRAPHORDERSTMETHSI,                            /*+ Simple                                +*/
  HGRAPHORDERSTMETHNBR                            /*+ Number of methods                     +*/
} HgraphOrderStMethodType;

/*+ Method function pointer type. +*/

typedef int (* HgraphOrderFunc) (Hgraph * restrict const, Order * restrict const, const Gnum, OrderCblk * restrict const, const void * const);

/*
**  The external declarations.
*/

extern StratTab             hgraphorderststratab;

/*
**  The function prototypes.
*/

int                         hgraphOrderSt       (Hgraph * restrict const, Order * restrict const, const Gnum, OrderCblk * restrict const, const Strat * const);
