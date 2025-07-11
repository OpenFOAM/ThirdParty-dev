/* Copyright 2004,2007,2011-2013,2015,2018-2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : graph_coarsen.h                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the source graph coarsening         **/
/**                functions.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 02 dec 1992     **/
/**                                 to   : 18 may 1993     **/
/**                # Version 1.3  : from : 30 apr 1994     **/
/**                                 to   : 18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to   : 18 aug 1994     **/
/**                # Version 3.0  : from : 07 jul 1995     **/
/**                                 to   : 28 sep 1995     **/
/**                # Version 3.1  : from : 28 nov 1995     **/
/**                                 to   : 28 nov 1995     **/
/**                # Version 3.2  : from : 07 sep 1996     **/
/**                                 to   : 17 sep 1998     **/
/**                # Version 4.0  : from : 13 dec 2001     **/
/**                                 to   : 05 dec 2004     **/
/**                # Version 6.0  : from : 09 mar 2011     **/
/**                                 to   : 30 aug 2020     **/
/**                # Version 7.0  : from : 28 jul 2018     **/
/**                                 to   : 19 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines.
*/

#if (! defined SCOTCH_PTHREAD) && (! defined GRAPHCOARSENNOTHREAD)
#define GRAPHCOARSENNOTHREAD
#endif /* (! defined SCOTCH_PTHREAD) && (! defined GRAPHCOARSENNOTHREAD) */

/*+ Graph option flags. Their values must be equal
    to those defined in library.h and library_f.h  +*/

#define GRAPHCOARSENNONE            0x0000        /* No options set */

#define GRAPHCOARSENDSTMATE         0x0001        /* Matching/fine-to-coarse array destination provided */
#define GRAPHCOARSENDSTMULT         0x0002        /* Multinode array destination provided               */
#define GRAPHCOARSENHASMULT         0x0004        /* Multinode array provided                           */
#define GRAPHCOARSENUSEMATE         0x0008        /* Matching array data provided                       */

#define GRAPHCOARSENNOCOMPACT       0x1000        /* Create a non-compact graph                         */
#define GRAPHCOARSENDETERMINISTIC   0x2000        /* Use deterministic algorithms only                  */
#define GRAPHCOARSENNOMERGE         0x4000        /* Do not merge isolated vertices                     */

/*+ Prime number for hashing vertex numbers. +*/

#define GRAPHCOARSENHASHPRIME       1049          /*+ Prime number +*/

/*
**  The type and structure definitions.
*/

/*+ Here are the edge matching function types for coarsening. +*/

typedef enum GraphCoarsenType_ {
  GRAPHCOARHEM,                                   /*+ Heavy-edge matching       +*/
  GRAPHCOARSCN,                                   /*+ Scanning (first) matching +*/
  GRAPHCOARNBR                                    /*+ Number of matching types  +*/
} GraphCoarsenType;

/*+ The multinode table element, which contains
    pairs of based indices of collapsed vertices.
    Both values are equal for uncollapsed vertices.
    As the base values of the fine and coarse graphs
    may be different, the values of the collapsed
    vertices are set with respect to the base value
    of the fine graph.                               +*/

typedef struct GraphCoarsenMulti_ {
  Gnum                      vertnum[2];           /*+ Numbers of the collapsed vertices of a multinode +*/
} GraphCoarsenMulti;

/*+ A table made of such elements is used during
    coarsening to build the edge array of the new
    graph, after the labeling of the vertices.    +*/

typedef struct GraphCoarsenHash_ {
  Gnum                      vertorgnum;           /*+ Origin vertex (i.e. pass) number +*/
  Gnum                      vertendnum;           /*+ Other end vertex number          +*/
  Gnum                      edgenum;              /*+ Number of corresponding edge     +*/
} GraphCoarsenHash;

/*+ The thread-specific data block. +*/

typedef struct GraphCoarsenThread_ {
  GraphCoarsenHash *        coarhashtab;          /*+ End vertex hash table (may be local)                    +*/
  Gnum                      coarvertnnd;          /*+ After-last coarse vertex number                         +*/
  Gnum                      coarvertbas;          /*+ Minimum coarse vertex number                            +*/
  Gnum                      coarvertnbr;          /*+ Number of coarse vertices to date                       +*/
  Gnum                      coaredgebas;          /*+ Minimum coarse edge number                              +*/
  Gnum                      coaredloadj;          /*+ Local coarse edge load sum adjust                       +*/
  Gnum                      coardegrmax;          /*+ Local maximum degree                                    +*/
  Gnum                      finevertbas;          /*+ Start of fine vertex range                              +*/
  Gnum                      finevertnnd;          /*+ End of fine vertex range                                +*/
  Gnum *                    finequeutab;          /*+ Queue array (may be global)                             +*/
  Gnum                      finequeudlt;          /*+ Multiplicative factor for queue slot size on first pass +*/
  Gnum                      finequeunbr;          /*+ Number of vertices in (local) queue                     +*/
  Gnum                      scantab[2];
} GraphCoarsenThread;

/*+ The matching and coarsening routine
    parameter structure. It contains the
    thread-independent data.             +*/

typedef struct GraphCoarsenData_ {
  int                       flagval;              /*+ Flags for controlling matching and coarsening   +*/
  const Graph *             finegrafptr;          /*+ Fine graph to perform matching on               +*/
  const Anum *              fineparotax;          /*+ Old part array                                  +*/
  const Anum *              finepfixtax;          /*+ Array of fixed vertices                         +*/
  Gnum                      finevfixnbr;          /*+ Number of fine fixed vertices                   +*/
  Gnum *                    finematetax;          /*+ Fine mate array / fine-to-coarse array          +*/
  Graph *                   coargrafptr;          /*+ Coarse graph to build                           +*/
  Gnum                      coarvertmax;          /*+ Maximum number of vertices to get               +*/
  Gnum                      coarvertnbr;          /*+ Global number of coarse vertices after matching +*/
  GraphCoarsenMulti *       coarmulttab;          /*+ Multinode array                                 +*/
  Gnum                      coarmultsiz;          /*+ Size of multinode array allocated in graph      +*/
  Gnum                      coarhashmsk;          /*+ Hash table mask                                 +*/
  int *                     finelocktax;          /*+ Global matching lock array (if any)             +*/
  GraphCoarsenThread *      thrdtab;              /*+ Array of thread-specific data                   +*/
  int                       fumaval;              /*+ Index of mating routine in function array       +*/
  volatile int              retuval;              /*+ Return value                                    +*/
  Context *                 contptr;              /*+ Execution context                               +*/
} GraphCoarsenData;

/*
**  The function prototypes.
*/

#ifdef SCOTCH_GRAPH_COARSEN
#ifndef GRAPHCOARSENNOTHREAD
static void                 graphCoarsenEdgeCt  (const GraphCoarsenData * restrict const, GraphCoarsenThread * restrict const);
#endif /* GRAPHCOARSENNOTHREAD */
static void                 graphCoarsenEdgeLl  (const GraphCoarsenData * restrict const, GraphCoarsenThread * restrict const);
static void                 graphCoarsenEdgeLu  (const GraphCoarsenData * restrict const, GraphCoarsenThread * restrict const);
#endif /* SCOTCH_GRAPH_COARSEN */

int                         graphCoarsen        (const Graph * restrict const, Graph * restrict const, Gnum * restrict * restrict const, GraphCoarsenMulti * restrict * restrict const, const Gnum, const double, const Gnum, const Anum * restrict const, const Anum * restrict const, const Gnum, Context * restrict const);
int                         graphCoarsenMatch   (const Graph * restrict const, Gnum * restrict * restrict const, Gnum * restrict const, const double, const Gnum, const Anum * restrict const, const Anum * restrict const, const Gnum, Context * restrict const);
int                         graphCoarsenBuild   (const Graph * restrict const, Graph * restrict const, Gnum * restrict const, GraphCoarsenMulti * restrict * restrict const, const Gnum, Context * restrict const);
