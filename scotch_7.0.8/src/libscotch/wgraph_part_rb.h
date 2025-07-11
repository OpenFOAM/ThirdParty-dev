/* Copyright 2010,2018,2019,2021 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : wgraph_part_rb.h                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Jun-Ho HER (v6.0)                       **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the vertex overlapped graph partit- **/
/**                ioning based on recursive bipartitioni- **/
/**                ng approach.                            **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 16 mar 2010     **/
/**                                 to   : 31 may 2018     **/
/**                # Version 6.1  : from : 21 nov 2021     **/
/**                                 to   : 23 nov 2021     **/
/**                # Version 7.0  : from : 23 aug 2019     **/
/**                                 to   : 26 nov 2021     **/
/**                                                        **/
/**   NOTES      : # This code derives from the code of    **/
/**                  kgraph_map_rb_part.h for the vertex   **/
/**                  overlapped graph partitioning.        **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct WgraphPartRbParam_ {
  Strat *                   straptr;              /*+ Bipartitioning strategy used +*/
} WgraphPartRbParam;

/*+ This structure holds global data. +*/

typedef struct WgraphPartRbData_ {
  const Graph *             grafptr;              /*+ Pointer to top-level graph          +*/
  Anum *                    parttax;              /*+ Pointer to top-level part array     +*/
  Gnum *                    frontab;              /*+ Pointer to top-level frontier array +*/
  Gnum                      fronnbr;              /*+ Current number of frontier vertices +*/
  Strat *                   straptr;              /*+ Bipartitioning strategy used        +*/
#ifdef SCOTCH_PTHREAD
  pthread_mutex_t           mutedat;              /*+ Mutex for frontier updates          +*/
#endif /* SCOTCH_PTHREAD */
} WgraphPartRbData;

/*+ This structure holds the splitting parameters. +*/

typedef struct WgraphPartRbSplit2_ {
  Gnum                      vertnbr;              /*+ Number of vertices in part   +*/
  Anum                      domnnum;              /*+ Initial domain number to map +*/
  Anum                      domnsiz;              /*+ Number of domains to map     +*/
} WgraphPartRbSplit2;

typedef struct WgraphPartRbSplit_ {
  WgraphPartRbSplit2        splttab[2];           /*+ Array of induced subgraph data  +*/
  WgraphPartRbData *        dataptr;              /*+ Pointer to global data          +*/
  Graph *                   grafptr;              /*+ Graph to induce and bipartition +*/
  Gnum *                    frontab;              /*+ Graph frontier array            +*/
  Gnum                      fronnbr;              /*+ Number of frontier vertices     +*/
  GraphPart *               parttax;              /*+ Graph part array                +*/
  int *                     revaptr;              /*+ Pointer to return value         +*/
} WgraphPartRbSplit;

/*
**  The function prototypes.
*/

int                         wgraphPartRb        (Wgraph * restrict const, const WgraphPartRbParam * restrict const);
