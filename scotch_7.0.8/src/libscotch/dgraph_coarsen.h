/* Copyright 2007-2009,2012,2014,2018,2020,2021,2023,2024 Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_coarsen.h                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                                                        **/
/**   FUNCTION   : This file implements the distributed    **/
/**                graph coarsening method.                **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 27 Jul 2005     **/
/**                                 to   : 24 feb 2007     **/
/**                # Version 5.1  : from : 11 nov 2008     **/
/**                                 to   : 26 may 2009     **/
/**                # Version 6.0  : from : 18 sep 2012     **/
/**                                 to   : 07 jun 2018     **/
/**                # Version 6.1  : from : 17 jun 2021     **/
/**                                 to   : 17 jun 2021     **/
/**                # Version 7.0  : from : 14 jan 2020     **/
/**                                 to   : 31 jul 2024     **/
/**                                                        **/
/************************************************************/

/*
** The defines.
*/

/*+ Graph option flags. Their values must be equal
    to those defined in library.h and library_f.h  +*/

#define DGRAPHCOARSENNONE           0x0000        /* No options set                 */
#define DGRAPHCOARSENFOLD           0x0100        /* Do folding without duplication */
#define DGRAPHCOARSENFOLDDUP        0x0300        /* Do folding with duplication    */
#define DGRAPHCOARSENNOMERGE        0x4000        /* Do not merge isolated vertices */

/*
** The type and structure definitions.
*/

/*+ The multinode table element, which contains
    pairs of based indices of collapsed vertices.
    Both values are equal for uncollapsed vertices. +*/

typedef struct DgraphCoarsenMulti_  {
  Gnum                      vertglbnum[2];        /*+ Global indices of the collapsed vertices of a multinode +*/
} DgraphCoarsenMulti;

/* This structure defines vertex data exchange
   cell. It can have many uses depending on the
   type of algorithm.                           */

typedef struct DgraphCoarsenVert_ {
  Gnum                      datatab[2];           /*+ Two values +*/
} DgraphCoarsenVert;

/* This structure defines the inter-process vertex and
   edge count structure. At matching time, it records
   the amount of data to be sent to each neighbor process. */

typedef struct DgraphCoarsenCount_ {
  Gnum                      vertsndnbr;           /*+ Number of vertices to send to this neighbor +*/
  Gnum                      edgesndnbr;           /*+ Number of edges to send to this neighbor    +*/
  Gnum                      vertlocnbr;           /*+ Number of local vertices in this neighbor   +*/
} DgraphCoarsenCount;

/*+ A table made of such elements is used during
    coarsening to build the edge array of the new
    graph, after the labeling of the vertices.    +*/

typedef struct DgraphCoarsenHash_ {
  Gnum                      vertorgnum;           /*+ Origin vertex (i.e. pass) number +*/
  Gnum                      vertendnum;           /*+ Other end vertex number          +*/
  Gnum                      edgelocnum;           /*+ Number of corresponding edge     +*/
} DgraphCoarsenHash;


/*+ The thread-specific data block, which
    has two uses. Firstly, it helps define
    the repartition across threads of the
    adjacency lists received from the
    different neighbors. Secondly, it
    allows reduction of graph properties.  +*/

typedef struct DgraphCoarsenThread_  {
  Gnum                      edgelocsum;           /*+ Upper bound on sum of local (fine) edges for this thread and below  +*/
  Gnum                      velolocsum;           /*+ Sum of vertex loads for each thread, for final reduce-sum operation +*/
  Gnum                      edgelocnbr;           /*+ Number of edges for each thread, for final reduce-sum operation     +*/
  Gnum                      degrlocmax;           /*+ Maximum degree for each thread, for final reduce-max operation      +*/
  int                       retuval;              /*+ Return value                                                        +*/
} DgraphCoarsenThread;

/*+ This structure gathers all data necessary
    to the proper execution of the coarsening
    and matching routines.                    +*/

typedef struct DgraphCoarsenData_ {
  int                       flagval;              /*+ Flag value                                                   +*/
  Dgraph *                  finegrafptr;          /*+ Pointer to fine graph                                        +*/
  Dgraph *                  coargrafptr;          /*+ Pointer to coarse graph which is built                       +*/
  int *                     coarprvptr;           /*+ Pointer to coarse private data to free in case of error      +*/
  DgraphCoarsenVert *       vrcvdattab;           /*+ Area reserved for receiving vertex messages                  +*/
  DgraphCoarsenVert *       vsnddattab;           /*+ Area reserved for sending vertex messages                    +*/
  int *                     vrcvcnttab;           /*+ Count data for vertex receive sub-arrays                     +*/
  int *                     vsndcnttab;           /*+ Count data for vertex send sub-arrays                        +*/
  int *                     vrcvdsptab;           /*+ Displacement for vertex receive sub-arrays [+1]              +*/
  int *                     vsnddsptab;           /*+ Displacement data for vertex send sub-arrays [+1]            +*/
  int *                     nrcvidxtab;           /*+ Count array for neighbor receive sub-arrays                  +*/
  int *                     nsndidxtab;           /*+ Count array for neighbor send sub-arrays                     +*/
  MPI_Request *             nrcvreqtab;           /*+ Request array for receive requests                           +*/
  MPI_Request *             nsndreqtab;           /*+ TRICK: nsndreqtab = (nrcvreqtab + procngbnbr)                +*/
  int *                     procgsttax;           /*+ Array giving the neighbor process index of each ghost vertex +*/
  int                       procngbnxt;           /*+ Index of first neighbor of higher rank than current process  +*/
  DgraphCoarsenCount *      dcntloctab;           /*+ Count array for sending vertices and edges                   +*/
  DgraphCoarsenCount *      dcntglbtab;           /*+ Count array for receiving vertices and edges                 +*/
  Gnum *                    coargsttax;           /*+ Fine-to-coarse vertex index array                            +*/
  DgraphCoarsenMulti *      multloctmp;           /*+ Pointer to multloctab structure to free (if any)             +*/
  DgraphCoarsenMulti *      multloctab;           /*+ Structure which contains the result of the matching          +*/
  Gnum                      multlocnbr;           /*+ Index of next multinode to be created                        +*/
#ifdef SCOTCH_DEBUG_DGRAPH2
  Gnum                      multlocsiz;           /*+ Size of allocated multinode array                            +*/
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  Gnum                      vertrcvnbr;           /*+ Number of fine vertices to be received                       +*/
  Gnum                      edgercvnbr;           /*+ Number of fine edges to be received                          +*/
  Gnum                      edgekptnbr;           /*+ Upper bound on number of edges kept from finer graph         +*/
  Gnum                      vertsndnbr;           /*+ Number of fine vertices to be sent                           +*/
  Gnum                      edgesndnbr;           /*+ Number of fine edges to be sent                              +*/
  Gnum *                    ercvdattab;           /*+ Pointer to adjacency exchange receive array                  +*/
  int *                     ercvdsptab;           /*+ Displacement data for adjacency exchange receive array       +*/
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  int *                     ercvcnttab;           /*+ Count data for adjacency exchange receive array              +*/
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
  Gnum                      degrlocmax;           /*+ Local maximum degree of coarse graph                         +*/
  Gnum                      coarhashmsk;          /*+ Mask value of adjacency hash table (TRICK: size - 1)         +*/
  DgraphCoarsenThread *     thrdtab;              /*+ Optional array for adjacency splitting across threads        +*/
  int                       thrdmin;              /*+ Minimum number of threads across all processes               +*/
  Context *                 contptr;              /*+ Execution context                                            +*/
} DgraphCoarsenData;

/*
** The function prototypes.
*/

#ifdef SCOTCH_DGRAPH_COARSEN
static int                  dgraphCoarsenInit   (DgraphCoarsenData * restrict const, Dgraph * restrict const, Dgraph * restrict const);
static void                 dgraphCoarsenExit   (DgraphCoarsenData * restrict const);
static int                  dgraphCoarsenBuild  (DgraphCoarsenData * restrict const);
#endif /* SCOTCH_DGRAPH_COARSEN */

int                         dgraphCoarsen       (Dgraph * restrict const, Dgraph * restrict const, DgraphCoarsenMulti * restrict * const, const Gnum, const Gnum, const double, const int, Context * restrict const);
Gnum                        dgraphCoarsenVertLocMax (const Dgraph * restrict const, const int);
