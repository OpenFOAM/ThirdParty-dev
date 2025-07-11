/* Copyright 2007-2010,2018,2019,2021,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph.h                                **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                Charles-Edmond BICHOT (5.1b)            **/
/**                                                        **/
/**   FUNCTION   : This module contains the data declarat- **/
/**                ions for the vertex overlapped graph    **/
/**                partitioning.                           **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 01 dec 2007     **/
/**                                 to   : 01 jul 2008     **/
/**                # Version 6.0  : from : 05 nov 2009     **/
/**                                 to   : 31 may 2018     **/
/**                # Version 6.1  : from : 23 nov 2021     **/
/**                                 to   : 02 dec 2021     **/
/**                # Version 7.0  : from : 23 aug 2019     **/
/**                                 to   : 16 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The graph structure. +*/

typedef struct Wgraph_ {
  Graph                     s;                    /*+ Source graph                            +*/
  Anum                      partnbr;              /*+ Current number of parts                 +*/
  Gnum                      fronnbr;              /*+ Number of frontier vertices             +*/
  Gnum                      fronload;             /*+ load for frontier                       +*/
  Gnum *                    frontab;              /*+ Array of frontier vertex numbers        +*/
  Gnum *                    compload;             /*+ Array of part loads                     +*/
  Gnum *                    compsize;             /*+ Array of part number of vertives        +*/
  Anum *                    parttax;              /*+ Part array; can be allocated separately +*/
  INT                       levlnum;              /*+ Coarsening level                        +*/
  Context *                 contptr;              /*+ Execution context                       +*/
} Wgraph;

/*+ The save graph structure. +*/

typedef struct WgraphStore_ {
  Anum                      partnbr;              /*+ Number of parts             +*/
  Gnum                      fronnbr;              /*+ Number of frontier vertices +*/
  Gnum                      fronload;             /*+ Frontier load               +*/
  byte *                    datatab;              /*+ Variable-sized data array   +*/
} WgraphStore;

/*+ This structure stores part lists. +*/

typedef struct WgraphPartList_ {
  Gnum                      vertnum;              /*+ Number of vertex of which part is neighbor +*/
  Anum                      nextidx;              /*+ Pointer to index of next recorded neighbor +*/
} WgraphPartList;

/*
**  The function prototypes.
*/

void                        wgraphInit          (Wgraph * const, const Graph * restrict const, const Anum);
void                        wgraphExit          (Wgraph * const);
int                         wgraphAlloc         (Wgraph * const);
void                        wgraphZero          (Wgraph * const);
int                         wgraphCheck         (const Wgraph * const);
int                         wgraphCost          (Wgraph * const);

int                         wgraphStoreInit     (const Wgraph * const, WgraphStore * const);
void                        wgraphStoreExit     (WgraphStore * const);
void                        wgraphStoreSave     (const Wgraph * const, WgraphStore * const);
void                        wgraphStoreUpdt     (Wgraph * const, const WgraphStore * const);
