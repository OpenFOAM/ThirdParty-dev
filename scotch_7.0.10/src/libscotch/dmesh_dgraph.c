/* Copyright 2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dmesh_dgraph.c                          **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Marc FUENTES                            **/
/**                                                        **/
/**   FUNCTION   : This module contains the source         **/
/**                distributed mesh to distributed graph   **/
/**                conversion function.                    **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 20 jan 2023     **/
/**                                 to   : 27 aug 2025     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dmesh.h"
#include "dmesh_dgraph.h"

/* This routine increases the size of the edge
** hash table array.
** It returns:
** - 0   : if resizing succeeded.
** - !0  : if out of memory.
*/

static
int
dmeshDgraphDualHashEdgeResize (
DmeshDgraphDualHashEdge * restrict *  hashtabptr,
Gnum *                                hashsizptr,
Gnum *                                hashmaxptr,
Gnum *                                hashmskptr,
const Gnum                            vertnum)    /* Current vertex number for used slots */
{
  DmeshDgraphDualHashEdge * restrict  hashtab;
  Gnum                                hashold;
  Gnum                                hashbas;
  Gnum                                hashnnd;
  Gnum                                hashsiz;
  Gnum                                hashmsk;

  hashold = *hashsizptr;
  hashsiz = hashold * 2;
  if ((hashtab = memRealloc (*hashtabptr, hashsiz * sizeof (DmeshDgraphDualHashEdge))) == NULL) /* Resize hash table */
    return (1);

  memSet (hashtab + hashold, ~0, hashold * sizeof (DmeshDgraphDualHashEdge)); /* Initialize second half */

  for (hashbas = hashold - 1; hashtab[hashbas].vertnum == vertnum; hashbas --) ; /* Find start index of last block */
  hashnnd = hashold;                              /* First segment to reconsider ends at the end of the old array  */
  hashmsk = hashsiz - 1;
  while (hashnnd != hashbas) {                    /* For each of the two segments to consider */
    Gnum                hashnum;

    for (hashnum = hashbas; hashnum < hashnnd; hashnum ++) { /* Re-compute position in new table */
      if (hashtab[hashnum].vertnum == vertnum) {  /* If hash slot used in this pass              */
        Gnum                vertend;
        Gnum                hashnew;

        vertend = hashtab[hashnum].vertend;       /* Get hash key value */
        for (hashnew = (vertend * DMESHDGRAPHHASHPRIME) & hashmsk; ; hashnew = (hashnew + 1) & hashmsk) {
          if (hashnew == hashnum)                 /* If hash slot is the same */
            break;                                /* There is nothing to do   */
          if (hashtab[hashnew].vertnum != vertnum) { /* If new slot is empty  */
            hashtab[hashnew] = hashtab[hashnum];  /* Copy data to new slot    */
            hashtab[hashnum].vertnum = ~0;        /* Mark old slot as empty   */
            break;
          }
        }                                         /* Go on searching */
      }
    }

    hashnnd = hashbas;                            /* End of second segment to consider is start of first one    */
    hashbas = 0;                                  /* Start of second segment is beginning of array              */
  }                                               /* After second segment, hashbas = hashnnd = 0 and loop stops */

  *hashtabptr  = hashtab;
  *hashsizptr  = hashsiz;
  *hashmaxptr *= 2;                               /* Adjust maximum capacity proportionally, whatever it is */
  *hashmskptr  = hashmsk;

  return (0);
}

static
int
dmeshDgraphDualProcNum (
Gnum                        vertglbnbr,
Gnum                        vertlocnum,           /* Un-based vertex number */
int                         procglbnbr)           /* Number of processes    */
{
  Gnum                bloksizmin;
  Gnum                bloksizmax;
  Gnum                blokrmnval;
  Gnum                proclocnum;

  blokrmnval = 1 + (vertglbnbr - 1) % (Gnum) procglbnbr; /* Index of first process of size bloksizmin (TRICK: or procglbnbr if it would be 0) */
  bloksizmax = ((vertglbnbr - 1) + (Gnum) procglbnbr) / (Gnum) procglbnbr; /* Maximum size of block per process (procnum < blokrmnval)        */
  bloksizmin = bloksizmax - 1;                    /* Minimum size of block per process (proclocnum >= blokrmnval)                             */

  proclocnum = vertlocnum / bloksizmax;           /* Get estimated index of block */
  if ((proclocnum >= blokrmnval) && (bloksizmin > 0))
    proclocnum = blokrmnval + (vertlocnum - blokrmnval * bloksizmax) / bloksizmin;

  return ((int) proclocnum);
}

/*************************/
/*                       */
/* The distributed graph */
/* building routines.    */
/*                       */
/*************************/

/* This routine builds a distributed dual graph
** (that is, an element graph) from the given
** distributed mesh. An edge is built between any two
** element vertices if these two elements e1 and e2
** have at least noconbr nodes in common.
** It returns:
** - 0  : if the distributed graph has been successfully built.
** - 1  : on error.
*/

int
dmeshDgraphDual (
const Dmesh * restrict const  meshptr,            /*+ Original mesh                                           +*/
Dgraph * restrict const       grafptr,            /*+ Graph to build                                          +*/
const Gnum                    noconbr)            /*+ Number of common points to define adjacency of elements +*/
{
  Gnum                                velmlocnum;
  Gnum                                eelmlocnum;
  int * restrict                      vnflloctax; /* Flag array for sending node vertices to processes     */
  Gnum                                vnodglbnbr;
  Gnum                                vnodglbmax; /* Highest node global index across all processes        */
  Gnum                                vnodlocmax; /* Highest node global index on local process            */
#ifdef SCOTCH_DEBUG_DMESH2
  Gnum                                vnodglbmin; /* Lowest node global number across all processes        */
  Gnum                                vnodlocmin; /* Lowest node global number on local process            */
#endif /* SCOTCH_DEBUG_DMESH2 */
  Gnum                                vnodlocbas; /* First global node index on local process              */
  Gnum                                vnodlocnnd; /* After-last global node index on local process         */
  Gnum                                vnodlocnbr; /* Number of node vertices on local process              */
  Gnum * restrict                     vnodloctax; /* Vertex index array for local node vertices            */
  Gnum                                vnodlocnum;
  Gnum                                vnodlocadj; /* Temporary variable to adjust node vertex indices      */
  Gnum * restrict                     enodloctax; /* Edge element adjacency array for local node vertices  */
  Gnum * restrict                     prfrloctab; /* Start of per-process linked list of edge data to send */ 
  Gnum * restrict                     eeneloctax; /* After-start linked list of per-process edge data      */
  int * restrict                      ercvcnttab; /* Edge receive count array                              */
  int * restrict                      ercvdsptab; /* Edge receive displacement array                       */
  Gnum * restrict                     ercvdattab; /* Receive array for element-node edges                  */
  int                                 ercvdatsiz; /* Size of edge receive array                            */
  int                                 ercvdatidx; /* Index to traverse edge receive array entirely         */
  int * restrict                      esndcnttab; /* Edge send count array                                 */
  int * restrict                      esnddsptab; /* Edge send displacement array                          */
  int                                 esnddatsiz; /* Size of edge send array                               */
  Gnum * restrict                     esnddattab; /* Send array for element-node edges                     */
#define nrcvcnttab                  ercvcnttab    /* Node adjacency receive count array; TRICK: reused     */
  int * restrict                      nrcvdsptab; /* Node adjacency receive displacement array             */
  int                                 nrcvdatsiz; /* Size of node adjacency receive array                  */
  int                                 nrcvdatidx; /* Index to traverse node receive array entirely         */
  Gnum * restrict                     nrcvdattab; /* Receive array for node adjacencies                    */
#define nsndcnttab                  esndcnttab    /* Node adjacency send count array; TRICK: reused        */
#define nsnddsptab                  esnddsptab    /* Node adjacency send displacement array; TRICK: reused */
  int                                 nsnddspval; /* Node adjacency send displacement value                */
  int                                 nsnddatsiz; /* Size of node adjacency send array                     */
  Gnum * restrict                     nsnddattab; /* Send array for node adjacencies                       */
  Gnum * restrict                     vertloctax; /* Vertex array of distributed dual graph                */
  Gnum * restrict                     edgeloctax; /* Edge array of distributed dual graph                  */
  Gnum                                edgelocnnd; /* End marker for current size of edge array             */
  Gnum                                edgelocnum;
  Gnum                                degrlocmax; /* Local maximum degree                                  */
  DmeshDgraphDualHashEdge * restrict  helmtab;    /* Hash table for building dual graph edge adjacencies   */
  Gnum                                helmsiz;
  Gnum                                helmmsk;
  Gnum                                helmmax;
  Gnum                                helmnbr;
  Gnum                                hnodnbr;    /* Expected number of nodes in node hash table           */
  DmeshDgraphDualHashNode * restrict  hnodtab;    /* Hash table for locating node adjacencies in array     */
  Gnum                                hnodsiz;
  Gnum                                hnodmsk;
  int                                 procnum;

  const int                         procglbnbr = meshptr->procglbnbr;
  const int                         proclocnum = meshptr->proclocnum;
  const Gnum                        baseval    = meshptr->baseval;
  const Gnum                        velmlocbas = meshptr->prelvrttab[proclocnum];
  const Gnum                        velmlocnnd = meshptr->prelvrttab[proclocnum + 1];
  const Gnum                        velmlocnbr = velmlocnnd - velmlocbas;
  const Gnum * restrict       const velmloctax = meshptr->velmloctab - velmlocbas; /* TRICK: base with respect to global indices */
  const Gnum * restrict       const eelmloctax = meshptr->eelmloctab - baseval;
  const Gnum                        eelmlocnnd = velmloctax[velmlocnnd]; /* Because element vertex array is compact */
  const Gnum                        eelmlocnbr = eelmlocnnd - baseval;

#ifdef SCOTCH_DEBUG_DMESH2
  vnodlocmin = GNUMMAX;
#endif /* SCOTCH_DEBUG_DMESH2 */
  vnodlocmax = 0;
  for (eelmlocnum = baseval; eelmlocnum < eelmlocnnd; eelmlocnum ++) { /* Search for highest node vertex */
    Gnum                vnodlocend;

    vnodlocend = eelmloctax[eelmlocnum];
#ifdef SCOTCH_DEBUG_DMESH2
    if (vnodlocmin > vnodlocend)
      vnodlocmin = vnodlocend;
#endif /* SCOTCH_DEBUG_DMESH2 */
    if (vnodlocmax < vnodlocend)
      vnodlocmax = vnodlocend;
  }

  if (MPI_Allreduce (&vnodlocmax, &vnodglbmax, 1, GNUM_MPI, MPI_MAX, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshDgraphDual: communication error (1)");
    return (1);
  }
#ifdef SCOTCH_DEBUG_DMESH2
  if (MPI_Allreduce (&vnodlocmin, &vnodglbmin, 1, GNUM_MPI, MPI_MIN, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshDgraphDual: communication error (2)");
    return (1);
  }
  if (vnodglbmin != baseval) {
    errorPrint ("dmeshDgraphDual: invalid parameters");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DMESH2 */
  vnodglbnbr = vnodglbmax - baseval + 1;

  if (((ercvcnttab = memAlloc (procglbnbr * sizeof (int) * 5))  == NULL) || /* TRICK: allocate all cnt/dsp arrays as well      */
      ((prfrloctab = memAlloc (procglbnbr * sizeof (Gnum)))     == NULL) || /* Heads of linked lists to send to each process   */
      ((eeneloctax = memAlloc (eelmlocnbr * sizeof (Gnum) * 2)) == NULL)) { /* TRICK: double array with element and edge index */
    errorPrint ("dmeshDgraphDual: out of memory (1)");
    return (1);
  }
  eeneloctax -= baseval;
  ercvdsptab  = ercvcnttab + procglbnbr;
  esndcnttab  = ercvdsptab + procglbnbr;
  esnddsptab  = esndcnttab + procglbnbr;
  nrcvdsptab  = esnddsptab + procglbnbr;
  memSet (esndcnttab,  0, procglbnbr * sizeof (int)); /* Initialize data send count array */
  memSet (prfrloctab, ~0, procglbnbr * sizeof (Gnum)); /* TRICK: set all values to -1     */

  for (velmlocnum = velmlocbas, eelmlocnum = baseval; /* Chain all edges to send per destination process */
       velmlocnum < velmlocnnd; velmlocnum ++) {
    Gnum                eelmlocnnd;

    for (eelmlocnnd = velmloctax[velmlocnum + 1]; eelmlocnum < eelmlocnnd; eelmlocnum ++) {
      Gnum                vnodlocend;             /* End node value      */
      int                 procnum;                /* Destination process */

      vnodlocend = eelmloctax[eelmlocnum];
      procnum = dmeshDgraphDualProcNum (vnodglbnbr, vnodlocend - baseval, procglbnbr); /* Compute destination process for given (un-based) edge end */

      esndcnttab[procnum] += 2;                   /* Count edge data to send to given process         */
      eeneloctax[2 * eelmlocnum]     = velmlocnum; /* Record element vertex number                    */
      eeneloctax[2 * eelmlocnum + 1] = prfrloctab[procnum]; /* Chain edge number of node data to send */
      prfrloctab[procnum] = eelmlocnum;
    }
  }

  if (MPI_Alltoall (esndcnttab, 1, MPI_INT,       /* Exchange counts of data to be sent/received to/from each process */
                    ercvcnttab, 1, MPI_INT, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshDgraphDual: communication error (3)");
    return (1);
  }

  for (procnum = 0, esnddatsiz = ercvdatsiz = 0; procnum < procglbnbr; procnum ++) {
    esnddsptab[procnum] = esnddatsiz;             /* Set start indices for building send arrays */
    esnddatsiz += esndcnttab[procnum];
    ercvdsptab[procnum] = ercvdatsiz;             /* Set start indices for building receive arrays */
    ercvdatsiz += ercvcnttab[procnum];            /* Sum-up size of data receive array             */
  }

  vnodlocbas = DATASCAN (vnodglbnbr, procglbnbr, proclocnum)     + baseval;
  vnodlocnnd = DATASCAN (vnodglbnbr, procglbnbr, proclocnum + 1) + baseval;
  vnodlocnbr = vnodlocnnd - vnodlocbas;

  if (((vnodloctax = memAlloc ((vnodlocnbr + 1) * sizeof (Gnum))) == NULL) ||
      ((vnflloctax = memAlloc (vnodlocnbr * sizeof (int)))        == NULL) ||
      ((enodloctax = memAlloc (ercvdatsiz * sizeof (Gnum) / 2))   == NULL) || /* Node edge array size is half the received data pairs */
      ((ercvdattab = memAlloc (ercvdatsiz * sizeof (Gnum)))       == NULL) ||
      ((esnddattab = memAlloc (esnddatsiz * sizeof (Gnum)))       == NULL)) { /* Will be freed first */
    errorPrint ("dmeshDgraphDual: out of memory (2)");
    return (1);
  }

  for (procnum = 0; procnum < procglbnbr; procnum ++) { /* Build per-process send arrays from edge linked lists */
    Gnum                esnddspidx;

    esnddspidx = esnddsptab[procnum];
    for (eelmlocnum = prfrloctab[procnum]; eelmlocnum != -1; eelmlocnum = eeneloctax[2 * eelmlocnum + 1]) { /* For all edges in linked list */
      esnddattab[esnddspidx]     = eelmloctax[eelmlocnum]; /* Send node number                 */
      esnddattab[esnddspidx + 1] = eeneloctax[2 * eelmlocnum]; /* Send neighbor element number */
      esnddspidx += 2;                            /* Two more data items added to send array   */
    }
#ifdef SCOTCH_DEBUG_DMESH2
    if ((procnum < (procglbnbr - 1)) &&
        (esnddspidx != esnddsptab[procnum + 1])) {
      errorPrint ("dmeshDgraphDual: internal error (1)");
      return (1);
    }
#endif /* SCOTCH_DEBUG_DMESH2 */
  }

  memFree (eeneloctax + baseval);                 /* Free linked list of element edges */
  memFree (prfrloctab);

  memSet (vnodloctax, 0, vnodlocnbr * sizeof (Gnum)); /* Initialize neighbor count array                          */
  vnodloctax -= vnodlocbas;                       /* TRICK: base node vertex array with respect to global indices */
  enodloctax -= baseval;                          /* Node edge array is commonly based with respect to base value */

  if (MPI_Alltoallv (esnddattab, esndcnttab, esnddsptab, GNUM_MPI,
                     ercvdattab, ercvcnttab, ercvdsptab, GNUM_MPI, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshDgraphDual: communication error (4)");
    return (1);
  }

  memFree (esnddattab);                           /* Free send array for first phase */

  for (ercvdatidx = 0; ercvdatidx < ercvdatsiz; ercvdatidx += 2) { /* Count adjacency of each node vertex */
    Gnum                vnodlocnum;

    vnodlocnum = ercvdattab[ercvdatidx];          /* Get node global index in edge */
#ifdef SCOTCH_DEBUG_DMESH2
    if ((vnodlocnum <  vnodlocbas) ||             /* If received node number is out of bounds for process */
        (vnodlocnum >= vnodlocnnd)) {
      errorPrint ("dmeshDgraphDual: internal error (2)");
      return (1);
    }
#endif /* SCOTCH_DEBUG_DMESH2 */

    vnodloctax[vnodlocnum] ++;                    /* One more edge neighbor for node vertex */
  }

  for (vnodlocnum = vnodlocbas, vnodlocadj = baseval; vnodlocnum < vnodlocnnd; vnodlocnum ++) { /* Set indices to edge array */
    Gnum                vnodloctmp;

    vnodloctmp = vnodloctax[vnodlocnum];
    vnodloctax[vnodlocnum] = vnodlocadj;
    vnodlocadj += vnodloctmp;
  }

  for (ercvdatidx = 0; ercvdatidx < ercvdatsiz; ercvdatidx += 2) { /* Build adjacency of each node vertex */
    Gnum                velmlocnum;
    Gnum                vnodlocnum;

    vnodlocnum = ercvdattab[ercvdatidx];          /* Get edge pair of global indices */
    velmlocnum = ercvdattab[ercvdatidx + 1];

    enodloctax[vnodloctax[vnodlocnum] ++] = velmlocnum; /* Add edge to node adjacency */
  }

  memMov (vnodloctax + vnodlocbas + 1, vnodloctax + vnodlocbas, vnodlocnbr * sizeof (Gnum)); /* Rebuild indices in node vertex array */
  vnodloctax[vnodlocbas] = baseval;

  memSet (vnflloctax, ~0, vnodlocnbr * sizeof (int)); /* Initialize node vertex flag array with invalid process number */
  vnflloctax -= vnodlocbas;                       /* TRICK: base node vertex flag array with respect to global indices */

  for (procnum = 0, ercvdatidx = 0; procnum < procglbnbr; procnum ++) { /* For all processes holding elements */
    Gnum                ercvdatnnd;               /* End index of edge data to consider for current process   */
    int                 nsndcntval;               /* Amount of node adjacency data to send to current process */

    nsndcntval = 0;                               /* No nodes to send to this process yet */
    for (ercvdatnnd = (procnum < (procglbnbr - 1)) ? ercvdsptab[procnum + 1] : ercvdatsiz;
         ercvdatidx < ercvdatnnd; ercvdatidx += 2) { /* Build adjacency of each node vertex */
      Gnum                vnodlocnum;

      vnodlocnum = ercvdattab[ercvdatidx];        /* Get number of node to send back to given process */

      if (vnflloctax[vnodlocnum] != procnum) {    /* If flag not set for this destination process       */
        nsndcntval += 2 + (int) (vnodloctax[vnodlocnum + 1] - vnodloctax[vnodlocnum]); /* Add data size */
        vnflloctax[vnodlocnum] = procnum;         /* Flag node vertex as already being sent             */
      }
    }

    nsndcnttab[procnum] = nsndcntval;
  }

  if (MPI_Alltoall (nsndcnttab, 1, MPI_INT,       /* Exchange counts of data to be sent/received to/from each process */
                    nrcvcnttab, 1, MPI_INT, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshDgraphDual: communication error (5)");
    return (1);
  }

  for (procnum = 0, nsnddatsiz = nrcvdatsiz = 0; procnum < procglbnbr; procnum ++) {
    nsnddsptab[procnum] = nsnddatsiz;             /* Set start indices for building send arrays */
    nsnddatsiz += nsndcnttab[procnum];
    nrcvdsptab[procnum] = nrcvdatsiz;             /* Set start indices for building receive arrays */
    nrcvdatsiz += nrcvcnttab[procnum];            /* Sum-up size of data receive array             */
  }

  edgelocnnd = nrcvdatsiz + baseval + 4;          /* Estimate on number of edges to be built; TRICK: "+4" for 25% increase to work */
  if (((vertloctax = memAlloc ((velmlocnbr + 1)       * sizeof (Gnum))) == NULL) || /* Dual graph vertex array                     */
      ((edgeloctax = memAlloc ((edgelocnnd - baseval) * sizeof (Gnum))) == NULL) ||
      ((nrcvdattab = memAlloc (nrcvdatsiz             * sizeof (Gnum))) == NULL) ||
      ((nsnddattab = memAlloc (nsnddatsiz             * sizeof (Gnum))) == NULL)) { /* Will be freed first */
    errorPrint ("dmeshDgraphDual: out of memory (3)");
    return (1);
  }
  vertloctax -= velmlocbas;                       /* TRICK: base vertex array with respect to global indices */
  edgeloctax -= baseval;

  for (procnum = 0, ercvdatidx = 0, nsnddspval = 0; /* For all processes holding elements */
       procnum < procglbnbr; procnum ++) {
    Gnum                ercvdatnnd;               /* End index of edge data to consider for current process   */

#ifdef SCOTCH_DEBUG_DMESH2
    if (nsnddspval != nsnddsptab[procnum]) {
      errorPrint ("dmeshDgraphDual: internal error (3)");
      return (1);
    }
#endif /* SCOTCH_DEBUG_DMESH2 */
    for (ercvdatnnd = (procnum < (procglbnbr - 1)) ? ercvdsptab[procnum + 1] : ercvdatsiz;
         ercvdatidx < ercvdatnnd; ercvdatidx += 2) { /* Build adjacency of each node vertex */
      Gnum                vnodlocnum;

      vnodlocnum = ercvdattab[ercvdatidx];        /* Get number of node to send back to given process */

      if (vnflloctax[vnodlocnum] != (procglbnbr + procnum)) { /* If flag not set for this destination process; TRICK: new range */
        Gnum                degrval;

        degrval = vnodloctax[vnodlocnum + 1] - vnodloctax[vnodlocnum];
#ifdef SCOTCH_DEBUG_DMESH2
        if ((nsnddspval + 2 + degrval) > nsnddatsiz) {
          errorPrint ("dmeshDgraphDual: internal error (4)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_DMESH2 */
        nsnddattab[nsnddspval ++] = vnodlocnum;   /* Send node number  */
        nsnddattab[nsnddspval ++] = degrval;      /* Send degree       */
        memCpy (nsnddattab + nsnddspval, enodloctax + vnodloctax[vnodlocnum], degrval * sizeof (Gnum)); /* Send neighbors */
        nsnddspval += degrval;                    /* Update start of area                          */
        vnflloctax[vnodlocnum] = (procglbnbr + procnum); /* Flag node vertex as already being sent */
      }
    }
  }
#ifdef SCOTCH_DEBUG_DMESH2
  if (nsnddspval != nsnddatsiz) {
    errorPrint ("dmeshDgraphDual: internal error (5)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DMESH2 */

  memFree (ercvdattab);                           /* Free receive array for first phase */
  memFree (enodloctax + baseval);
  memFree (vnflloctax + vnodlocbas);
  memFree (vnodloctax + vnodlocbas);

  if (MPI_Alltoallv (nsnddattab, nsndcnttab, nsnddsptab, GNUM_MPI,
                     nrcvdattab, nrcvcnttab, nrcvdsptab, GNUM_MPI, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshDgraphDual: communication error (6)");
    return (1);
  }

  memFree (nsnddattab);                           /* Free send array for second phase           */
  memFree (ercvcnttab);                           /* Free process count and displacement arrays */

 for (nrcvdatidx = 0, hnodnbr = 0; nrcvdatidx < nrcvdatsiz; ) { /* Count received nodes to set node hash table size */
    Gnum                degrval;

    degrval     = nrcvdattab[nrcvdatidx + 1];     /* Get node degree                  */
    nrcvdatidx += 2 + degrval;                    /* Skip adjacency to next node data */
    hnodnbr ++;                                   /* One more slot in node hash table */
  }

  for (hnodsiz = 32; hnodsiz < hnodnbr; hnodsiz *= 2) ; /* Get upper power of two                           */
  hnodsiz *= 2;                                   /* Load node hash table at 50% capacity since large array */

  if ((hnodtab = memAlloc (hnodsiz * sizeof (DmeshDgraphDualHashNode))) == NULL) {
    errorPrint ("dmeshDgraphDual: out of memory (4)");
    return (1);
  }
  memSet (hnodtab, ~0, hnodsiz * sizeof (DmeshDgraphDualHashNode));
  hnodmsk = hnodsiz - 1;

 for (nrcvdatidx = 0; nrcvdatidx < nrcvdatsiz; ) { /* Record all node locations in node hash table */
    Gnum                vnodlocnum;
    Gnum                degrval;
    Gnum                hnodnum;

    vnodlocnum = nrcvdattab[nrcvdatidx ++];       /* Get node number */
    degrval    = nrcvdattab[nrcvdatidx];          /* Get node degree */

    for (hnodnum = (vnodlocnum * DMESHDGRAPHHASHPRIME) & hnodmsk; ; hnodnum = (hnodnum + 1) & hnodmsk) {
      if (hnodtab[hnodnum].vnodnum == ~0) {       /* If hash slot empty    */
        hnodtab[hnodnum].vnodnum = vnodlocnum;    /* Record node location  */
        hnodtab[hnodnum].dataidx = nrcvdatidx;    /* Point to degree value */
        break;
      }
    }

    nrcvdatidx += 1 + degrval;                    /* Skip adjacency to next node data */
  }

  helmmax = 64;
  helmsiz = helmmax * 4;                          /* Load edge hash table at 1/4 capacity */
  if ((helmtab = memAlloc (helmsiz * sizeof (DmeshDgraphDualHashEdge))) == NULL) {
    errorPrint ("dmeshDgraphDual: out of memory (5)");
    return (1);
  }
  memSet (helmtab, ~0, helmsiz * sizeof (DmeshDgraphDualHashEdge));
  helmmsk = helmsiz - 1;

  degrlocmax = 0;
  for (velmlocnum = velmlocbas, edgelocnum = baseval; velmlocnum < velmlocnnd; velmlocnum ++) {
    Gnum                eelmlocnum;
    Gnum                eelmlocnnd;
    Gnum                degrlocval;
    Gnum                helmnum;

    vertloctax[velmlocnum] = edgelocnum;

    helmnum = (velmlocnum * DMESHDGRAPHHASHPRIME) & helmmsk; /* Prevent adding loop edge */
    helmtab[helmnum].vertnum = velmlocnum;
    helmtab[helmnum].vertend = velmlocnum;
    helmtab[helmnum].nghbnbr = 0;                 /* Loop edge never created as boundary already crossed */
    helmnbr = 1;                                  /* One entry added to edge hash table                  */

    eelmlocnum = velmloctax[velmlocnum];
    eelmlocnnd = velmloctax[velmlocnum + 1];

    for ( ; eelmlocnum < eelmlocnnd; eelmlocnum ++) { /* For all element edges */
      Gnum                vnodlocnum;
      int                 nrcvdatidx;
      int                 enoddatnum;
      int                 enoddatnnd;
      Gnum                hnodnum;

      vnodlocnum = eelmloctax[eelmlocnum];        /* Get number of node neighbor */
      for (hnodnum = (vnodlocnum * DMESHDGRAPHHASHPRIME) & hnodmsk; ; hnodnum = (hnodnum + 1) & hnodmsk) {
#ifdef SCOTCH_DEBUG_DMESH2
        if (hnodtab[hnodnum].vnodnum == ~0) { /* If hash slot empty */
          errorPrint ("dmeshDgraphDual: internal error (6)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_DMESH2 */
        if (hnodtab[hnodnum].vnodnum == vnodlocnum) /* If hash slot found */
          break;
      }
      nrcvdatidx = hnodtab[hnodnum].dataidx;      /* Get index of node degree then adjacency  */
      enoddatnum = nrcvdatidx + 1;                /* Start at first index of adjacency        */
      enoddatnnd = enoddatnum + nrcvdattab[nrcvdatidx]; /* Compute end index from node degree */

      for ( ; enoddatnum < enoddatnnd; enoddatnum ++) {
        Gnum                velmlocend;
        Gnum                helmend;

        velmlocend = nrcvdattab[enoddatnum];      /* Get number of end element vertex */

        for (helmend = (velmlocend * DMESHDGRAPHHASHPRIME) & helmmsk; ; helmend = (helmend + 1) & helmmsk) {
          Gnum                nghbnbr;

          if (helmtab[helmend].vertnum != velmlocnum) { /* If edge not yet created      */
            if (helmnbr >= helmmax)               /* If edge hash table full, resize it */
              dmeshDgraphDualHashEdgeResize (&helmtab, &helmsiz, &helmmax, &helmmsk, velmlocnum);

            helmtab[helmend].vertnum = velmlocnum; /* Record new edge */
            helmtab[helmend].vertend = velmlocend;
            helmtab[helmend].nghbnbr =            /* One instance recorded to date */
            nghbnbr                  = noconbr - 1;
            helmnbr ++;                           /* One more item added to hash table              */
            goto test;                            /* Check if one instance is enough to create edge */
          }
          if (helmtab[helmend].vertend == velmlocend) { /* If hash slot found                     */
            nghbnbr = helmtab[helmend].nghbnbr;   /* Get number of times neighbor element met yet */
            if (nghbnbr > 0) {                    /* If edge not already created                  */
              helmtab[helmend].nghbnbr = -- nghbnbr; /* One more instance of neighbor element met */
test:         if (nghbnbr <= 0) {                 /* If new instance allows us to reach threshold */
                if (edgelocnum >= edgelocnnd) {   /* If edge array already full                   */
                  Gnum                  edgelocmax;
                  Gnum * restrict       edgeloctmp;

                  edgelocmax = edgelocnnd - baseval; /* Increase size by 25 % */
                  edgelocmax = edgelocmax + (edgelocmax >> 2);

                  if ((edgeloctmp = memRealloc (edgeloctax + baseval, edgelocmax * sizeof (Gnum))) == NULL) {
                    errorPrint ("dmeshDgraphDual: out of memory (6)");
                    return (1);
                  }

                  edgeloctax = edgeloctmp - baseval;
                  edgelocnnd = edgelocmax + baseval;
                }

                edgeloctax[edgelocnum ++] = velmlocend; /* Create edge */
              }
            }
            break;
          }
        }
      }
    }

    degrlocval = edgelocnum - vertloctax[velmlocnum]; /* Compute local maximum degree */
    if (degrlocval > degrlocmax)
      degrlocmax = degrlocval;
  }
  vertloctax[velmlocnum] = edgelocnum;            /* Set end of vertex array */

  memFree (helmtab);
  memFree (hnodtab);
  memFree (nrcvdattab);

  edgeloctax  = memRealloc (edgeloctax + baseval, (edgelocnum - baseval) * sizeof (Gnum)); /* Reduce size of edge array */
  edgeloctax -= baseval;

  if (dgraphBuild2 (grafptr, baseval,             /* Build distributed graph */
                    meshptr->velmlocnbr, meshptr->velmlocnbr,
                    vertloctax + velmlocbas - baseval, vertloctax + velmlocbas - baseval + 1,
                    NULL, meshptr->velmlocnbr, NULL, NULL,
                    edgelocnum - baseval, edgelocnum - baseval, edgeloctax, NULL, NULL, degrlocmax) != 0) {
    errorPrint ("dmeshDgraphDual: cannot build dual graph");
    memFree    (edgeloctax + baseval);
    memFree    (vertloctax + velmlocbas);
    return (1);
  }
  grafptr->flagval |= DGRAPHFREETABS;

#ifdef SCOTCH_DEBUG_DMESH2
  if (dgraphCheck (grafptr) != 0) {
    errorPrint ("dmeshDgraphDual: internal error (7)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DMESH2 */

  return (0);
}
