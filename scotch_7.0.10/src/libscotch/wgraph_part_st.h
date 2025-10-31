/* Copyright 2007-2010,2018,2020,2021,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph_part_st.h                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Charles-Edmond BICHOT (v5.1b)           **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declara-  **/
/**                tions for vertex overlapped graph par-  **/
/**                titioning strategy and method tables.   **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to   : 31 may 2018     **/
/**                # Version 6.1  : from : 25 aug 2020     **/
/**                                 to   : 02 dec 2021     **/
/**                # Version 7.0  : from : 07 nov 2024     **/
/**                                 to   : 07 nov 2024     **/
/**                                                        **/
/************************************************************/

/*
** The type definitions.
*/

/** Method types. **/

typedef enum WgraphPartStMethodType_ {
  WGRAPHPARTSTMETHES = 0,                         /*+ Edge separation        +*/
  WGRAPHPARTSTMETHFM,                             /*+ Fiduccia-Mattheyses    +*/
  WGRAPHPARTSTMETHML,                             /*+ Multi-level separation +*/
  WGRAPHPARTSTMETHRB,                             /*+ Recursive bisection    +*/
  WGRAPHPARTSTMETHZR,                             /*+ Zero method            +*/
  WGRAPHPARTSTMETHNBR                             /*+ Number of methods      +*/
} WgraphPartStMethodType;

/*+ Method function pointer type. +*/

typedef int (* WgraphPartFunc) (Wgraph * restrict const, const void * const);

/*
**  The external declarations.
*/

extern StratTab             wgraphpartststratab;

/*
**  The function prototypes.
*/

int                         wgraphPartSt        (Wgraph * restrict const, const Strat * restrict const);
