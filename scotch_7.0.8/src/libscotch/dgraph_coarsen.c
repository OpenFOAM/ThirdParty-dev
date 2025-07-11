/* Copyright 2007-2012,2014,2018-2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_coarsen.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Cedric CHEVALIER (v5.0)                 **/
/**                Selmane LEBDAOUI (v7.0, b_dcoarsen)     **/
/**                                                        **/
/**   FUNCTION   : This file implements the coarsening     **/
/**                phase of the multi-level method.        **/
/**                The implementation uses several         **/
/**                processes, which could have several     **/
/**                threads each (3 at this time).          **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 27 jul 2005     **/
/**                                 to   : 15 may 2008     **/
/**                # Version 5.1  : from : 23 jun 2008     **/
/**                                 to   : 20 feb 2011     **/
/**                # Version 6.0  : from : 11 sep 2012     **/
/**                                 to   : 07 jun 2018     **/
/**                # Version 6.1  : from : 17 jun 2021     **/
/**                                 to   : 27 dec 2021     **/
/**                # Version 7.0  : from : 14 jan 2020     **/
/**                                 to   : 31 jul 2024     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define SCOTCH_DGRAPH_COARSEN

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_allreduce.h"
#include "dgraph_coarsen.h"
#include "dgraph_match.h"

/******************************/
/*                            */
/* Graph coarsening routines. */
/*                            */
/******************************/

static
int
dgraphCoarsenInit (
DgraphCoarsenData * restrict const  coarptr,      /*+ Coarsening data structure +*/
Dgraph * restrict const             finegrafptr,  /*+ Graph to coarsen          +*/
Dgraph * restrict const             coargrafptr)  /*+ Coarse graph to build     +*/
{
  int                 procglbnbr;
  int                 procglbnum;
  int                 procngbnbr;
  int                 procngbnum;
  int                 procngbnxt;
  int                 vertrcvnbr;
  int                 vertsndnbr;
  Gnum                vertlocnbr;
  Gnum                vertgstnbr;
  int                 vdsprcvnum;
  int                 vdspsndnum;
  byte *              bufftab;
  size_t              buffsiz;

  const int * restrict const  fineprocngbtab = finegrafptr->procngbtab;
  const int * restrict const  fineprocrcvtab = finegrafptr->procrcvtab;
  const int * restrict const  fineprocsndtab = finegrafptr->procsndtab;

  vertlocnbr = finegrafptr->vertlocnbr;
  vertgstnbr = finegrafptr->vertgstnbr;
  procglbnbr = finegrafptr->procglbnbr;
  procngbnbr = finegrafptr->procngbnbr;
  vertrcvnbr = vertgstnbr - vertlocnbr;
  vertsndnbr = finegrafptr->procsndnbr;

  coarptr->coarprvptr = NULL;                     /* Assume nothing to free on error */
  coarptr->multloctmp = NULL;
  coarptr->nsndidxtab = NULL;
  coarptr->nrcvidxtab = NULL;
  coarptr->thrdtab    = NULL;

  if ((coarptr->coarprvptr = memAllocGroup ((void **) (void *) /* Allocate distributed coarse graph private data */
                                            &coargrafptr->procdsptab, (size_t) ((procglbnbr + 1) * sizeof (Gnum)),
                                            &coargrafptr->proccnttab, (size_t) (procglbnbr       * sizeof (Gnum)),
                                            &coargrafptr->procngbtab, (size_t) (procglbnbr       * sizeof (int)),
                                            &coargrafptr->procrcvtab, (size_t) (procglbnbr       * sizeof (int)),
                                            &coargrafptr->procsndtab, (size_t) (procglbnbr       * sizeof (int)), NULL)) == NULL) {
    errorPrint ("dgraphCoarsenInit: out of memory (1)");
    return (1);
  }
  coargrafptr->procvrttab = coargrafptr->procdsptab; /* Coarse graph has no holes */

  if (coarptr->multloctab == NULL) {              /* If no multinode array provided */
    Gnum                multlocsiz;               /* Size of local multinode array  */

    multlocsiz = vertlocnbr;                      /* Size of multinode array of plain coarsened graph, not taking folding into account */
#ifdef SCOTCH_DEBUG_DGRAPH2
    coarptr->multlocsiz = multlocsiz;
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    if ((coarptr->multloctab = memAlloc (multlocsiz * sizeof (DgraphCoarsenMulti))) == NULL) {
      errorPrint        ("dgraphCoarsenInit: out of memory (2)");
      dgraphCoarsenExit (coarptr);
      return (1);
    }
    coarptr->multloctmp = coarptr->multloctab;    /* Array will have to be freed on error */
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  else                                            /* User-provided multinode array */
    coarptr->multlocsiz = dgraphCoarsenVertLocMax (finegrafptr, DGRAPHCOARSENNONE); /* Impossible to know: use prescribed bound to test it */
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (memAllocGroup ((void **) (void *)           /* Data used up to edge exchange phase at coarse graph build time */
                     &coarptr->nrcvidxtab, (size_t) (procngbnbr * sizeof (int)),
                     &coarptr->vrcvdsptab, (size_t) ((procglbnbr + 1) * sizeof (int)), /* TRICK: "+1" for size count */
                     &coarptr->coargsttax, (size_t) (vertgstnbr * sizeof (Gnum)),
                     &coarptr->procgsttax, (size_t) (vertrcvnbr * sizeof (int)), /* TRICK: Only purely ghost part of array will be used */
                     &coarptr->vrcvdattab, (size_t) (vertrcvnbr * sizeof (DgraphCoarsenVert)), NULL) == NULL) {
    errorPrint        ("dgraphCoarsenInit: out of memory (3)");
    dgraphCoarsenExit (coarptr);
    return (1);
  }

  buffsiz = 2 * MAX ((procngbnbr * sizeof (MPI_Request)), (procglbnbr * sizeof (int)));
  if (memAllocGroup ((void **) (void *)           /* Data released after coarse vertex index exchange phase */
                     &coarptr->nsndidxtab, (size_t) (procngbnbr * sizeof (int)),
                     &coarptr->vsnddsptab, (size_t) ((procglbnbr + 1) * sizeof (int)), /* TRICK: "+1" for size count check */
                     &bufftab,             (size_t) buffsiz,
                     &coarptr->dcntloctab, (size_t) (procglbnbr * sizeof (DgraphCoarsenCount)),
                     &coarptr->dcntglbtab, (size_t) (procglbnbr * sizeof (DgraphCoarsenCount)),
                     &coarptr->vsnddattab, (size_t) (vertsndnbr * sizeof (DgraphCoarsenVert)), NULL) == NULL) {
    errorPrint        ("dgraphCoarsenInit: out of memory (4)");
    dgraphCoarsenExit (coarptr);
    return (1);
  }
  coarptr->nrcvreqtab = (MPI_Request *) (void *) bufftab; /* TRICK: point-to-point requests and collective arrays share same space */
  coarptr->nsndreqtab = coarptr->nrcvreqtab + procngbnbr;
  coarptr->vrcvcnttab = (int *) (void *) bufftab;
  coarptr->vsndcnttab = coarptr->vrcvcnttab + procglbnbr;

  for (procglbnum = 0, vdsprcvnum = vdspsndnum = 0; /* Build communication index arrays */
       procglbnum < procglbnbr; procglbnum ++) {
    coarptr->vrcvdsptab[procglbnum] = vdsprcvnum;
    coarptr->vsnddsptab[procglbnum] = vdspsndnum;
    vdsprcvnum += fineprocrcvtab[procglbnum];
    vdspsndnum += fineprocsndtab[procglbnum];
  }
  coarptr->vrcvdsptab[procglbnum] = vdsprcvnum;   /* Mark end of communication index arrays */
  coarptr->vsnddsptab[procglbnum] = vdspsndnum;

  for (procngbnum = procngbnxt = 0; procngbnum < procngbnbr; procngbnum ++) {
    if ((procngbnxt == 0) && (fineprocngbtab[procngbnum] > finegrafptr->proclocnum)) { /* Find index of first neighbor of higher rank */
      procngbnxt = procngbnum;
      break;
    }
  }
  coarptr->procngbnxt = procngbnxt;

  coarptr->coargsttax -= finegrafptr->baseval;
  coarptr->finegrafptr = finegrafptr;
  coarptr->coargrafptr = coargrafptr;

  memSet (coarptr->dcntloctab, 0, procglbnbr * sizeof (DgraphCoarsenCount));

  memSet (coarptr->procgsttax, ~0, vertrcvnbr * sizeof (int)); /* Values have not yet been computed                       */
  coarptr->procgsttax -= vertlocnbr + finegrafptr->baseval; /* TRICK: base array such that only purely ghost part is used */

  coarptr->edgekptnbr = 0;

  return (0);
}

static
void
dgraphCoarsenExit (
DgraphCoarsenData * restrict const    coarptr)    /*+ Coarsening data structure +*/
{
  if (coarptr->nsndidxtab != NULL)                /* Auxiliary array is released after first phase of coarse graph building */
    memFree (coarptr->nsndidxtab);
  if (coarptr->nrcvidxtab != NULL)
    memFree (coarptr->nrcvidxtab);
  if (coarptr->multloctmp != NULL)                /* If multinode array not provided nor passed back to calling routine */
    memFree (coarptr->multloctmp);
  if (coarptr->coarprvptr != NULL)                /* If ownership of coarse graph private data not yet transferred to it */
    memFree (coarptr->coarprvptr);
  if (coarptr->thrdtab != NULL)
    memFree (coarptr->thrdtab);
}

/*************************************/
/*                                   */
/* Multinode data exchange routines. */
/*                                   */
/*************************************/

static
int
dgraphCoarsenBuildColl (
DgraphCoarsenData * restrict const  coarptr)
{
  int                 procngbnum;

  MPI_Comm                      proccomm   = coarptr->finegrafptr->proccomm;
  Dgraph * restrict const       grafptr    = coarptr->finegrafptr;
  const Gnum                    vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval;
  const int                     procngbnbr = grafptr->procngbnbr;
  const int * restrict const    procngbtab = grafptr->procngbtab;
  Gnum * restrict const         coargsttax = coarptr->coargsttax;
  int * restrict const          vsndcnttab = coarptr->vsndcnttab;
  int * restrict const          vrcvdsptab = coarptr->coargrafptr->procrcvtab; /* TRICK: use coarse graph procrcvtab and procsndtab */
  int * restrict const          vsnddsptab = coarptr->coargrafptr->procsndtab;
  int * restrict const          nrcvidxtab = coarptr->nrcvidxtab;
  int * restrict const          nsndidxtab = coarptr->nsndidxtab;

  memSet (vsndcnttab, 0, grafptr->procglbnbr * sizeof (int));
  memSet (vrcvdsptab, 0, grafptr->procglbnbr * sizeof (int));
  memSet (vsnddsptab, 0, grafptr->procglbnbr * sizeof (int));
  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) {
    int                 procglbnum;

    procglbnum = procngbtab[procngbnum];
    vsndcnttab[procglbnum] = 2 * (nsndidxtab[procngbnum] - coarptr->vsnddsptab[procglbnum]);
    vrcvdsptab[procglbnum] = 2 * coarptr->vrcvdsptab[procglbnum];
    vsnddsptab[procglbnum] = 2 * coarptr->vsnddsptab[procglbnum];
  }

  if (MPI_Alltoall (vsndcnttab, 1, MPI_INT, coarptr->vrcvcnttab, 1, MPI_INT, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuildColl: communication error (1)");
    return (1);
  }
  if (MPI_Alltoallv (coarptr->vsnddattab, vsndcnttab,          vsnddsptab, GNUM_MPI,
                     coarptr->vrcvdattab, coarptr->vrcvcnttab, vrcvdsptab, GNUM_MPI, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuildColl: communication error (2)");
    return (1);
  }

  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) { /* For all received data chunks */
    int                 vrcvidxnnd;
    int                 vrcvidxnum;
    int                 procglbnum;
    int                 statsiz;

    const DgraphCoarsenVert * restrict const  vrcvdattab = coarptr->vrcvdattab; /* After data is received */

    procglbnum = procngbtab[procngbnum];
    statsiz = coarptr->vrcvcnttab[procglbnum];
    for (vrcvidxnum = coarptr->vrcvdsptab[procglbnum], vrcvidxnnd = vrcvidxnum + (statsiz / 2); /* TRICK: each message item costs 2 Gnum's */
         vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
      Gnum                vertglbnum;             /* Our global number (the one seen as mate by sender) */
      Gnum                vertlocnum;             /* Our local number (the one seen as mate by sender)  */
      Gnum                multglbnum;             /* Global number of coarse vertex                     */

      vertglbnum = vrcvdattab[vrcvidxnum].datatab[0];
      multglbnum = vrcvdattab[vrcvidxnum].datatab[1];
      vertlocnum = vertglbnum - vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((vertlocnum <  grafptr->baseval) ||     /* If matching request is not directed towards our process */
          (vertlocnum >= grafptr->vertlocnnd)) {
        errorPrint ("dgraphCoarsenBuildColl: internal error");
        return (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      coargsttax[vertlocnum] = multglbnum;
    }
    nrcvidxtab[procngbnum] = vrcvidxnnd;          /* Keep receive end index for preparing edge arrays */
  }

  return (0);
}

static
int
dgraphCoarsenBuildPtop (
DgraphCoarsenData * restrict const  coarptr)
{
  int                 procngbnum;
  int                 vrcvreqnbr;

  MPI_Comm                      proccomm   = coarptr->finegrafptr->proccomm;
  Dgraph * restrict const       grafptr    = coarptr->finegrafptr;
  const Gnum                    vertlocadj = grafptr->procvrttab[grafptr->proclocnum] - grafptr->baseval;
  const int                     procngbnbr = grafptr->procngbnbr;
  const int * restrict const    procngbtab = grafptr->procngbtab;
  Gnum * restrict const         coargsttax = coarptr->coargsttax;
  const int * restrict const    vrcvdsptab = coarptr->vrcvdsptab;
  const int * restrict const    vsnddsptab = coarptr->vsnddsptab;
  int * restrict const          nrcvidxtab = coarptr->nrcvidxtab;
  int * restrict const          nsndidxtab = coarptr->nsndidxtab;

  if (procngbnbr > 0) {                           /* If we have neighbors to communication with */
    procngbnum = coarptr->procngbnxt;             /* Post receives in descending order          */
    do {
      int                 procglbnum;

      procngbnum = (procngbnum + (procngbnbr - 1)) % procngbnbr; /* Pre-decrement neighbor rank */
      procglbnum = procngbtab[procngbnum];
      if (MPI_Irecv (coarptr->vrcvdattab + vrcvdsptab[procglbnum], 2 * (vrcvdsptab[procglbnum + 1] - vrcvdsptab[procglbnum]), GNUM_MPI,
                     procglbnum, TAGCOARSEN, proccomm, &coarptr->nrcvreqtab[procngbnum]) != MPI_SUCCESS) {
        errorPrint ("dgraphCoarsenBuildPtop: communication error (1)");
        return (1);
      }
    } while (procngbnum != coarptr->procngbnxt);

    procngbnum = coarptr->procngbnxt;             /* Post sends in ascending order */
    do {
      int                 procglbnum;

      procglbnum = procngbtab[procngbnum];
      if (MPI_Isend (coarptr->vsnddattab + vsnddsptab[procglbnum], 2 * (nsndidxtab[procngbnum] - vsnddsptab[procglbnum]), GNUM_MPI,
                     procglbnum, TAGCOARSEN, proccomm, &coarptr->nsndreqtab[procngbnum]) != MPI_SUCCESS) {
        errorPrint ("dgraphCoarsenBuildPtop: communication error (2)");
        return (1);
      }
      procngbnum = (procngbnum + 1) % procngbnbr; /* Post-increment neighbor rank */
    } while (procngbnum != coarptr->procngbnxt);
  }

  for (vrcvreqnbr = procngbnbr; vrcvreqnbr > 0; vrcvreqnbr --) { /* For all pending receive requests */
    int                 vrcvidxnnd;
    int                 vrcvidxnum;
    int                 procngbnum;
    MPI_Status          statdat;
    int                 statsiz;
    int                 o;

#ifdef SCOTCH_DEBUG_DGRAPH2
    procngbnum = vrcvreqnbr - 1;                  /* Receive messages in order */
    o = MPI_Wait (&coarptr->nrcvreqtab[procngbnum], &statdat);
#else /* SCOTCH_DEBUG_DGRAPH2 */
    o = MPI_Waitany (procngbnbr, coarptr->nrcvreqtab, &procngbnum, &statdat);
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    if ((o != MPI_SUCCESS) ||
        (MPI_Get_count (&statdat, GNUM_MPI, &statsiz) != MPI_SUCCESS)) {
      errorPrint ("dgraphCoarsenBuildPtop: communication error (3)");
      return (1);
    }
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (statdat.MPI_SOURCE != procngbtab[procngbnum]) {
      errorPrint ("dgraphCoarsenBuildPtop: internal error (1)");
      return (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    {
      const DgraphCoarsenVert * restrict const  vrcvdattab = coarptr->vrcvdattab; /* After data is received */

      for (vrcvidxnum = vrcvdsptab[procngbtab[procngbnum]], vrcvidxnnd = vrcvidxnum + (statsiz / 2); /* TRICK: each message item costs 2 Gnum's */
           vrcvidxnum < vrcvidxnnd; vrcvidxnum ++) {
        Gnum                vertglbnum;           /* Our global number (the one seen as mate by sender) */
        Gnum                vertlocnum;           /* Our local number (the one seen as mate by sender)  */
        Gnum                multglbnum;           /* Global number of coarse vertex                     */

        vertglbnum = vrcvdattab[vrcvidxnum].datatab[0];
        multglbnum = vrcvdattab[vrcvidxnum].datatab[1];
        vertlocnum = vertglbnum - vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
        if ((vertlocnum <  grafptr->baseval) ||   /* If matching request is not directed towards our process */
            (vertlocnum >= grafptr->vertlocnnd)) {
          errorPrint ("dgraphCoarsenBuildPtop: internal error (2)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        coargsttax[vertlocnum] = multglbnum;
      }
      nrcvidxtab[procngbnum] = vrcvidxnnd;        /* Keep receive end index for preparing edge arrays */
    }
  }

  if (MPI_Waitall (procngbnbr, coarptr->nsndreqtab, MPI_STATUSES_IGNORE) != MPI_SUCCESS) { /* Wait for send requests to complete */
    errorPrint ("dgraphCoarsenBuildPtop: communication error (4)");
    return (1);
  }

  return (0);
}

/******************************/
/*                            */
/* The coarse graph adjacency */
/* building routine.          */
/*                            */
/******************************/

/* This routine merges the adjacency of the two
** fine vertices of a multinode into a coarse
** vertex, using a hash table for merging. The
** first multinode vertex is always local, while
** the second one can be either local or remote.
** It returns:
** - >= 0  : new coaredgenum value.
*/

static
Gnum
dgraphCoarsenBuildAdj (
const Dgraph * restrict const       finegrafptr,
DgraphCoarsenMulti * restrict const multloctax,
const Gnum                          coarvertlocnum,
const Gnum                          coarvertglbnum,
Gnum * restrict const               coarveloloctax,
Gnum * restrict const               coaredgeloctax,
Gnum                                coaredgelocnum,
Gnum * restrict const               coaredloloctax,
const Gnum                          vertlocadj,   /* Fine vertex local adjustment */
const Gnum * restrict const         coargsttax,
int * restrict const                ercvdsptab,
const Gnum * restrict const         ercvdattab,
const int * restrict const          procgsttax,
DgraphCoarsenHash * restrict const  coarhashtab,
const Gnum                          coarhashmsk)
{
  Gnum                coarvelolocval;
  Gnum                vertlocnum;
  int                 i;

  const Gnum * restrict const       vertloctax = finegrafptr->vertloctax;
  const Gnum * restrict const       vendloctax = finegrafptr->vendloctax;
  const Gnum * restrict const       veloloctax = finegrafptr->veloloctax;
  const Gnum * restrict const       edgeloctax = finegrafptr->edgeloctax;
  const Gnum * restrict const       edgegsttax = finegrafptr->edgegsttax;
  const Gnum * restrict const       edloloctax = finegrafptr->edloloctax;

  i = 0;
  coarvelolocval = 0;
  vertlocnum = multloctax[coarvertlocnum].vertglbnum[0] - vertlocadj;
  while (1) {                                     /* Pseudo-infinite loop on both vertices of the multinode */
    Gnum                vertglbnum;
    Gnum                edgelocnum;
    Gnum                edgelocnnd;
    Gnum                degrlocval;
    int                 procngbnum;
    int                 ercvidxnum;

    coarvelolocval += (veloloctax != NULL) ? veloloctax[vertlocnum] : 1;
    for (edgelocnum = vertloctax[vertlocnum], edgelocnnd = vendloctax[vertlocnum]; /* Loop on edges of first (and sometimes second) local mate */
         edgelocnum < edgelocnnd; edgelocnum ++) {
      Gnum                coarvertglbend;
      Gnum                edlolocval;
      Gnum                h;

      coarvertglbend = coargsttax[edgegsttax[edgelocnum]];
      if (coarvertglbend == coarvertglbnum)       /* If end of collapsed edge */
        continue;

      edlolocval = (edloloctax != NULL) ? edloloctax[edgelocnum] : 1;
      for (h = (coarvertglbend * COARHASHPRIME) & coarhashmsk; ; h = (h + 1) & coarhashmsk) {
        if (coarhashtab[h].vertorgnum != coarvertglbnum) { /* If old slot           */
          coarhashtab[h].vertorgnum = coarvertglbnum; /* Mark it in reference array */
          coarhashtab[h].vertendnum = coarvertglbend;
          coarhashtab[h].edgelocnum = coaredgelocnum;
          coaredgeloctax[coaredgelocnum] = coarvertglbend; /* One more edge created */
          coaredloloctax[coaredgelocnum] = edlolocval;
          coaredgelocnum ++;
          break;                                  /* Give up hashing */
        }
        if (coarhashtab[h].vertendnum == coarvertglbend) { /* If coarse edge already exists */
          coaredloloctax[coarhashtab[h].edgelocnum] += edlolocval;
          break;                                  /* Give up hashing */
        }
      }
    }

    if (i ++ > 0)                                 /* If second local vertex has been processed, exit */
      break;

    vertglbnum = multloctax[coarvertlocnum].vertglbnum[1];

    if (vertglbnum >= 0) {                        /* If second multinode vertex is local */
      if ((vertglbnum - vertlocadj) == vertlocnum) /* If single multinode                */
        break;
      vertlocnum = (vertglbnum - vertlocadj);
      continue;
    }

    edgelocnum = -2 - vertglbnum;                 /* Get edge number associated with remote second multinode vertex    */
    multloctax[coarvertlocnum].vertglbnum[1] = edgeloctax[edgelocnum]; /* Set global number of second multinode vertex */
    procngbnum = procgsttax[edgegsttax[edgelocnum]]; /* Get number of neighbor process to which owns end vertex        */
    ercvidxnum = ercvdsptab[procngbnum];          /* Get start index of remote vertex adgacency list                   */
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (ercvidxnum < 0) {                         /* If we attempt to read an adjacency that has not been sent */
      errorPrint ("dgraphCoarsenBuildAdj: internal error (1)");
      return (GNUMMAX);
    }
    if (ercvdattab[ercvidxnum ++] != edgeloctax[edgelocnum]) { /* If we are not reading the proper adjacency data */
      errorPrint ("dgraphCoarsenBuildAdj: internal error (2)");
      return (GNUMMAX);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    degrlocval = ercvdattab[ercvidxnum ++];       /* Read vertex degree                                */
    coarvelolocval += (veloloctax != NULL) ? ercvdattab[ercvidxnum ++] : 1; /* Read vertex load if any */

    while (degrlocval -- > 0) {                   /* For all edges */
      Gnum                coarvertglbend;
      Gnum                edlolocval;
      Gnum                h;

      coarvertglbend = ercvdattab[ercvidxnum ++]; /* Get coarse end vertex number                 */
      edlolocval = (edloloctax != NULL) ? ercvdattab[ercvidxnum ++] : 1; /* Read edge load if any */
      if (coarvertglbend == coarvertglbnum)       /* If loop coarse edge, skip it                 */
        continue;

      for (h = (coarvertglbend * COARHASHPRIME) & coarhashmsk; ; h = (h + 1) & coarhashmsk) {
        if (coarhashtab[h].vertorgnum != coarvertglbnum) { /* If old slot           */
          coarhashtab[h].vertorgnum = coarvertglbnum; /* Mark it in reference array */
          coarhashtab[h].vertendnum = coarvertglbend;
          coarhashtab[h].edgelocnum = coaredgelocnum;
          coaredgeloctax[coaredgelocnum] = coarvertglbend; /* One more edge created */
          coaredloloctax[coaredgelocnum] = edlolocval;
          coaredgelocnum ++;
          break;                                  /* Give up hashing */
        }
        if (coarhashtab[h].vertendnum == coarvertglbend) { /* If coarse edge already exists */
          coaredloloctax[coarhashtab[h].edgelocnum] += edlolocval;
          break;                                  /* Give up hashing */
        }
      }
    }

    ercvdsptab[procngbnum] = ercvidxnum;          /* Write back updated receive index       */
    break;                                        /* Exit loop after processing remote mate */
  }
  coarveloloctax[coarvertlocnum] = coarvelolocval;

  return (coaredgelocnum);
}

/***********************************/
/*                                 */
/* Coarse graph building routines. */
/*                                 */
/***********************************/

/* This routine finalizes the creation of the
** distributed coarse graph in a sequential way.
** The distributed coarse graph is compact.
** It returns:
** - 0   : if coarse graph was created.
** - !0  : on error.
*/

static
int
dgraphCoarsenBuildSeq (
DgraphCoarsenData * restrict const  coarptr)
{
  Gnum                          coarvertlocnum;
  Gnum                          coarvertlocnnd;
  Gnum                          coarvelolocsum;
  Gnum                          coardegrlocmax;
  Gnum                          coaredgelocnum;
  DgraphCoarsenHash * restrict  coarhashtab;      /* Table for merging vertex edges to same multinode */
  size_t                        coarhashsiz;      /* Size of hash table                               */

  Dgraph * restrict const             finegrafptr    = coarptr->finegrafptr;
  Dgraph * restrict const             coargrafptr    = coarptr->coargrafptr;
  const Gnum * restrict const         coargsttax     = coarptr->coargsttax;
  const Gnum                          vertlocadj     = finegrafptr->procvrttab[finegrafptr->proclocnum] - finegrafptr->baseval;
  const int * restrict const          procgsttax     = coarptr->procgsttax;
  const Gnum * restrict const         ercvdattab     = coarptr->ercvdattab;
  int * restrict const                ercvdsptab     = coarptr->ercvdsptab; /* Use global displacement data */
  DgraphCoarsenMulti * restrict const multloctax     = coarptr->multloctab - finegrafptr->baseval;
  const Gnum                          multlocadj     = coargrafptr->procdsptab[finegrafptr->proclocnum] - finegrafptr->baseval;
  Gnum * restrict const               coarvertloctax = coargrafptr->vertloctax;
  Gnum * restrict const               coarveloloctax = coargrafptr->veloloctax;
  Gnum * restrict const               coaredgeloctax = coargrafptr->edgeloctax;
  Gnum * restrict const               coaredloloctax = coargrafptr->edloloctax;
  const Gnum                          coarhashmsk    = coarptr->coarhashmsk;

  coarhashsiz = (coarhashmsk + 1) * sizeof (DgraphCoarsenHash); /* TRICK: (coarhashmsk + 1) is power of two */
  if ((coarhashtab = memAlloc (coarhashsiz)) == NULL) {
    errorPrint ("dgraphCoarsenBuildSeq: out of memory");
    return (1);
  }
  memSet (coarhashtab, ~0, coarhashsiz);

  coarvelolocsum = 0;
  coardegrlocmax = 0;
  for (coarvertlocnum = coaredgelocnum = finegrafptr->baseval, coarvertlocnnd = coarvertlocnum + coargrafptr->vertlocnbr;
       coarvertlocnum < coarvertlocnnd; coarvertlocnum ++) {
    coarvertloctax[coarvertlocnum] = coaredgelocnum;

    coaredgelocnum = dgraphCoarsenBuildAdj (finegrafptr, multloctax, coarvertlocnum, coarvertlocnum + multlocadj,
                                            coarveloloctax, coaredgeloctax, coaredgelocnum, coaredloloctax,
                                            vertlocadj, coargsttax, ercvdsptab, ercvdattab, procgsttax,
                                            coarhashtab, coarhashmsk);
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (coaredgelocnum > (coargrafptr->edgelocsiz + coargrafptr->baseval)) { /* Number of local edges can be reached, not exceeded */
      errorPrint ("dgraphCoarsenBuildSeq: internal error");
      return (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    coarvelolocsum += coarveloloctax[coarvertlocnum];
    if (coardegrlocmax < (coaredgelocnum - coarvertloctax[coarvertlocnum]))
      coardegrlocmax = (coaredgelocnum - coarvertloctax[coarvertlocnum]);
  }
  coarvertloctax[coarvertlocnum] = coaredgelocnum; /* Set end of compact edge array */

  coargrafptr->degrglbmax = coardegrlocmax;       /* Save local maximum degree before subsequent reduction */
  coargrafptr->velolocsum = coarvelolocsum;
  coargrafptr->edgelocnbr =
  coargrafptr->edgelocsiz = coaredgelocnum - coargrafptr->baseval; /* Compact edge array */

  memFree (coarhashtab);

  return (0);
}

/* This routine aggregates a sum and max
** reduction of partial coarse graph
** parameters computed by multiple
** threads.
*/

#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
static
void
dgraphCoarsenBuildThrReduce (
DgraphCoarsenThread * restrict const  tlocptr,    /* Pointer to local thread block  */
DgraphCoarsenThread * restrict const  tremptr,    /* Pointer to remote thread block */
const void * const                    globptr)    /* Unused                         */
{
  tlocptr->velolocsum += tremptr->velolocsum;     /* Sum vertex loads        */
  tlocptr->edgelocnbr += tremptr->edgelocnbr;     /* Sum number of edges     */
  if (tremptr->degrlocmax > tlocptr->degrlocmax)  /* Take maximum of degrees */
    tlocptr->degrlocmax = tremptr->degrlocmax;
  tlocptr->retuval |= tremptr->retuval;           /* Agregate return values */
}

/* This routine finalizes the creation of the
** distributed coarse graph in a multi-threaded way.
** The distributed coarse graph is not compact, so
** coargrafptr->vendloctab must have been allocated.
** To manage the case where processes may have
** different numbers of threads, only thrdmin threads
** are used per process, instead of thrdnbr.
** It returns:
** - 0   : if coarse graph was created.
** - !0  : on error.
*/

static
void
dgraphCoarsenBuildThr (
ThreadDescriptor * restrict const   descptr,
DgraphCoarsenData * restrict const  coarptr)
{
  Gnum                          coarvertlocnum;
  Gnum                          coarvertlocnnd;
  Gnum                          coarvelolocsum;
  Gnum                          coaredgelocbas;   /* Start index of edge sub-array for current thread */
  Gnum                          coaredgelocnum;
  Gnum                          coardegrlocmax;
  DgraphCoarsenHash * restrict  coarhashtab;      /* Table for merging vertex edges to same multinode */
  size_t                        coarhashsiz;      /* Size of hash table                               */
  int * restrict                ercvdsptab;
  int                           procngbnum;
  int                           o;

  const int                           thrdmin        = coarptr->thrdmin; /* All processes must have same number of threads */
  const int                           thrdnum        = threadNum (descptr);
  Dgraph * restrict const             finegrafptr    = coarptr->finegrafptr;
  Dgraph * restrict const             coargrafptr    = coarptr->coargrafptr;
  const Gnum * restrict const         coargsttax     = coarptr->coargsttax;
  const Gnum                          vertlocadj     = finegrafptr->procvrttab[finegrafptr->proclocnum] - finegrafptr->baseval;
  const int * restrict const          procgsttax     = coarptr->procgsttax;
  const int                           procngbnbr     = finegrafptr->procngbnbr;
  const Gnum * restrict const         ercvdattab     = coarptr->ercvdattab;
  DgraphCoarsenMulti * restrict const multloctax     = coarptr->multloctab - finegrafptr->baseval;
  const Gnum                          multlocadj     = coargrafptr->procdsptab[finegrafptr->proclocnum] - finegrafptr->baseval;
  Gnum * restrict const               coarvertloctax = coargrafptr->vertloctax;
  Gnum * restrict const               coarvendloctax = coargrafptr->vendloctax;
  Gnum * restrict const               coarveloloctax = coargrafptr->veloloctax;
  Gnum * restrict const               coaredgeloctax = coargrafptr->edgeloctax;
  Gnum * restrict const               coaredloloctax = coargrafptr->edloloctax;
  const Gnum                          coarhashmsk    = coarptr->coarhashmsk;

  if (thrdnum >= thrdmin) {                       /* If extra, idle thread */
    coarptr->thrdtab[thrdnum].velolocsum = 0;     /* Nothing to do         */
    coarptr->thrdtab[thrdnum].edgelocnbr = 0;
    coarptr->thrdtab[thrdnum].degrlocmax = 0;
    o = 0;                                        /* Everything went well */
    goto abort2;                                  /* Skip all work        */
  }

  o = 1;                                          /* Assume an error */

  coarhashsiz = (coarhashmsk + 1) * sizeof (DgraphCoarsenHash); /* TRICK: (coarhashmsk + 1) is power of two */
  if (memAllocGroup ((void **) (void *)
                     &coarhashtab, (size_t) coarhashsiz,
                     &ercvdsptab,  (size_t) (procngbnbr * sizeof (int)), NULL) == NULL) {
    errorPrint ("dgraphCoarsenBuildThr: out of memory");
    goto abort2;
  }
  memSet (coarhashtab, ~0, coarhashsiz);

  coaredgelocbas = finegrafptr->baseval;
  if (thrdnum == 0) {                             /* First thread does not have a thread adjacency start index */
    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) {
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (coarptr->ercvcnttab[procngbnum] <= 0)   /* If no data received from this neighbor       */
        ercvdsptab[procngbnum] = -1;              /* Set canary value                             */
      else                                        /* In the general case, computation is harmless */
#endif /* SCOTCH_DEBUG_DGRAPH2 */
        ercvdsptab[procngbnum] = coarptr->ercvdsptab[procngbnum] + 2 * (thrdmin - 1); /* Skip thread index sub-array to reach start of array for thread 0 */
    }
  }
  else {                                          /* Other threads have a thread adjacency start index */
    coaredgelocbas += coarptr->thrdtab[thrdnum].edgelocsum; /* Skip all local edges for this thread    */
    for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) {
      if (coarptr->ercvcnttab[procngbnum] > 0) {  /* If we have received data from this neighbor */
        int                 ercvdspidx;           /* Start index of neighbor sub-array           */

        ercvdspidx = coarptr->ercvdsptab[procngbnum];
        ercvdsptab[procngbnum] = ercvdspidx + ercvdattab[ercvdspidx + (2 * thrdnum) - 1]; /* Get start index for thread slot       */
        coaredgelocbas += ercvdattab[ercvdspidx + (2 * thrdnum) - 2]; /* Skip sum of edges received from all neighbors before slot */
      }
#ifdef SCOTCH_DEBUG_DGRAPH2
      else                                        /* Set canary value */
        ercvdsptab[procngbnum] = -1;
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    }
  }

  coarvelolocsum = 0;
  coardegrlocmax = 0;
  coarvertlocnum = finegrafptr->baseval + DATASCAN (coargrafptr->vertlocnbr, thrdmin, thrdnum);
  coarvertlocnnd = finegrafptr->baseval + DATASCAN (coargrafptr->vertlocnbr, thrdmin, thrdnum + 1);
  for (coaredgelocnum = coaredgelocbas; coarvertlocnum < coarvertlocnnd; coarvertlocnum ++) {
    coarvertloctax[coarvertlocnum] = coaredgelocnum; /* Set start of vertex adjacency sub-array */

    coaredgelocnum = dgraphCoarsenBuildAdj (finegrafptr, multloctax, coarvertlocnum, coarvertlocnum + multlocadj,
                                            coarveloloctax, coaredgeloctax, coaredgelocnum, coaredloloctax,
                                            vertlocadj, coargsttax, ercvdsptab, ercvdattab, procgsttax,
                                            coarhashtab, coarhashmsk);
#ifdef SCOTCH_DEBUG_DGRAPH2
    if (coaredgelocnum > (coargrafptr->edgelocsiz + coargrafptr->baseval)) { /* Number of local edges can be reached, not exceeded */
      errorPrint ("dgraphCoarsenBuildThr: internal error");
      goto abort1;
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    coarvendloctax[coarvertlocnum] = coaredgelocnum; /* Set end of vertex adjacency sub-array */
    coarvelolocsum += coarveloloctax[coarvertlocnum];
    if (coardegrlocmax < (coaredgelocnum - coarvertloctax[coarvertlocnum]))
      coardegrlocmax = (coaredgelocnum - coarvertloctax[coarvertlocnum]);
  }

  o = 0;                                          /* Everything went well */
  coarptr->thrdtab[thrdnum].velolocsum = coarvelolocsum;
  coarptr->thrdtab[thrdnum].edgelocnbr = coaredgelocnum - coaredgelocbas;
  coarptr->thrdtab[thrdnum].degrlocmax = coardegrlocmax;
  if (thrdnum == (thrdmin - 1))
    coargrafptr->edgelocsiz = coaredgelocnum - finegrafptr->baseval; /* For non-compact edge array, array size is end of last edge sub-array */

#ifdef SCOTCH_DEBUG_DGRAPH2
abort1:
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  memFree (coarhashtab);                          /* Free group leader */
abort2:
  coarptr->thrdtab[thrdnum].retuval = o;

  threadReduce (descptr, coarptr->thrdtab, sizeof (DgraphCoarsenThread), (ThreadReduceFunc) dgraphCoarsenBuildThrReduce, 0, NULL); /* Sum edges and get maximum of degrmax */

  return;
}
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

/* This routine performs the coarsening of edges
** with respect to the coarmulttax array computed
** by dgraphMatch. All data must be available when
** running (all receptions done). This function is
** inspired by libscotch/src/graph_coarsen_edge.c.
*/

DGRAPHALLREDUCEMAXSUMOP (4, 1)

static
int
dgraphCoarsenBuild (
DgraphCoarsenData * restrict const  coarptr)
{
  Gnum                          vendlocsiz;       /* Size of local vertex end array (if any) */
  Gnum                          edgelocsiz;       /* Size of local coarse edge array         */
  int                           ercvdatsiz;       /* Size of adjacency receive data array    */
  int                           esnddatsiz;       /* Size of adjacency send data array       */
  Gnum * restrict               ercvdattab;       /* Adjacency receive data array            */
  Gnum * restrict               esnddattab;       /* Adjacency send data array               */
  int * restrict                ercvcnttab;       /* Adjacency receive count array           */
  int * restrict                esndcnttab;       /* Adjacency send count array              */
  int * restrict                ercvdsptab;       /* Adjacency receive displacement array    */
  int * restrict                esnddsptab;       /* Adjacency send displacement array       */
  int                           ercvdspidx;       /* Current receive displacement index      */
  int                           esnddspidx;       /* Current send displacement index         */
  Gnum                          multlocnum;
  int                           procngbnum;
  int                           procnum;
  Gnum                          coarhashnbr;      /* Size of adjacency hash table            */
  int                           vertmltval;
  int                           edgemltval;
  Gnum                          reduloctab[5];
  Gnum                          reduglbtab[5];
  int                           cheklocval;
  int                           chekglbval;
#ifdef SCOTCH_DEBUG_DGRAPH2
  Gnum                          vertlocnum;
  int * restrict                ercvdbgtab;
#endif /* SCOTCH_DEBUG_DGRAPH2 */
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  Gnum                          coarlocnnd;       /* Boundary for current thread slot in neighbor    */
  Gnum                          edgelocsum;       /* Cumulative number of edges for this thread slot */
  int                           procrcvnbr;       /* Number of neighbors actually sending data to us */
  int                           thrdnum;
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  const int                                 thrdnbr     = contextThreadNbr (coarptr->contptr);
  const int                                 thrdmin     = coarptr->thrdmin;
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
  MPI_Comm                                  proccomm    = coarptr->finegrafptr->proccomm;
  Dgraph * restrict const                   finegrafptr = coarptr->finegrafptr;
  Dgraph * restrict const                   coargrafptr = coarptr->coargrafptr;
  Gnum * restrict const                     coargsttax  = coarptr->coargsttax;
  const int                                 procngbnbr  = finegrafptr->procngbnbr;
  const int * restrict const                procngbtab  = finegrafptr->procngbtab;
  const int * restrict const                procgsttax  = coarptr->procgsttax;
  const Gnum * restrict const               vertloctax  = finegrafptr->vertloctax;
  const Gnum * restrict const               vendloctax  = finegrafptr->vendloctax;
  const Gnum * restrict const               veloloctax  = finegrafptr->veloloctax;
  const Gnum * restrict const               edgeloctax  = finegrafptr->edgeloctax;
  const Gnum * restrict const               edgegsttax  = finegrafptr->edgegsttax;
  const Gnum * restrict const               edloloctax  = finegrafptr->edloloctax;
  const Gnum                                vertlocadj  = finegrafptr->procvrttab[finegrafptr->proclocnum] - finegrafptr->baseval;
  const Gnum                                multlocadj  = coargrafptr->procdsptab[finegrafptr->proclocnum];
  const DgraphCoarsenMulti * restrict const multloctab  = coarptr->multloctab;
  DgraphCoarsenVert * const                 vrcvdattab  = coarptr->vrcvdattab; /* [norestrict:async] */
  DgraphCoarsenVert * restrict const        vsnddattab  = coarptr->vsnddattab;
  int * restrict const                      nsndidxtab  = coarptr->nsndidxtab;

#ifdef SCOTCH_DEBUG_DGRAPH2
  memSet (coargsttax + finegrafptr->baseval, ~0, finegrafptr->vertgstnbr * sizeof (Gnum)); /* To ckeck no data will be left uninitialized */
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) /* Reset indices for sending messages */
    nsndidxtab[procngbnum] = coarptr->vsnddsptab[procngbtab[procngbnum]];

#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  if ((coarptr->thrdtab = memAlloc (thrdnbr * sizeof (DgraphCoarsenThread))) == NULL) { /* TRICK: thrdnbr and not thrdmin for reduction across all threads */
    errorPrint ("dgraphCoarsenBuild: out of memory (1)");
    return (1);                                   /* TODO: cleaner exit */
  }

  coarptr->thrdtab[0].edgelocsum = 0;             /* No local edges exist before thread 0 */
  thrdnum = 1;
  coarlocnnd = DATASCAN (coarptr->multlocnbr, thrdmin, 1); /* Compute first thread slot boundary */
  edgelocsum = 0;
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
  for (multlocnum = 0; multlocnum < coarptr->multlocnbr; multlocnum ++) { /* Loop on local multinode data */
    Gnum                ver0locnum;
    Gnum                ver1locnum;
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    Gnum                deg0locval;               /* Degree of first fine multinode vertex */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    while (multlocnum >= coarlocnnd) {            /* If crossed boundary of current thread    */
      coarptr->thrdtab[thrdnum].edgelocsum = edgelocsum; /* Record sum of local edges to date */
      thrdnum ++;
      coarlocnnd = DATASCAN (coarptr->multlocnbr, thrdmin, thrdnum); /* Compute next thread slot boundary */
    }
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

    ver0locnum = multloctab[multlocnum].vertglbnum[0] - vertlocadj; /* Compute local value of (always local) first multinode vertex */ 
#ifdef SCOTCH_DEBUG_DGRAPH2
    if ((ver0locnum <  finegrafptr->baseval) ||
        (ver0locnum >= finegrafptr->vertlocnnd)) {
      errorPrint ("dgraphCoarsenBuild: internal error (1)");
      return (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
    coargsttax[ver0locnum] = multlocnum + multlocadj; /* Un-based number with base adjustment */
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    deg0locval  = vendloctax[ver0locnum] - vertloctax[ver0locnum];
    edgelocsum += deg0locval;                     /* Account for local edges of first multinode vertex */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

    ver1locnum = multloctab[multlocnum].vertglbnum[1];
    if (ver1locnum >= 0) {                        /* If second vertex is local */
      ver1locnum -= vertlocadj;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((ver1locnum <  finegrafptr->baseval) ||
          (ver1locnum >= finegrafptr->vertlocnnd)) {
        errorPrint ("dgraphCoarsenBuild: internal error (2)");
        return (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      coargsttax[ver1locnum] = multlocnum + multlocadj; /* Valid if single multinode or not */
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
      if (ver1locnum != ver0locnum) {             /* If multinode does not comprise a single vertex                          */
        edgelocsum += vendloctax[ver1locnum] - vertloctax[ver1locnum]; /* Account for edges of second local multinode vertex */
        if (deg0locval > 0)                       /* If multinode is not created by merging an isolated vertex to some other */
          edgelocsum -= 2;                        /* Remove collapsed arcs                                                   */
      }
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
    }
    else {                                        /* Second vertex is not local             */
      Gnum                edgelocnum;             /* Fine edge number to remote fine vertex */
      Gnum                ver1glbnum;             /* Global number of fine end vertex       */
      Gnum                ver1gstnum;             /* Local ghost number of fine end vertex  */
      int                 coarsndidx;             /* Index in request send array            */
      int                 procngbnum;             /* Number of target process               */

#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
      edgelocsum --;                              /* Remove collapsed ghost edge */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
      edgelocnum = -2 - ver1locnum;
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((edgelocnum <   finegrafptr->baseval) ||
          (edgelocnum >= (finegrafptr->baseval + finegrafptr->edgelocsiz))) {
        errorPrint ("dgraphCoarsenBuild: internal error (3)");
        return (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      ver1glbnum = edgeloctax[edgelocnum];
      ver1gstnum = edgegsttax[edgelocnum];

      procngbnum = procgsttax[ver1gstnum];        /* Find neighbor owner process */
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (procngbnum < 0) {                       /* If neighbor had not been computed */
        errorPrint ("dgraphCoarsenBuild: internal error (4)");
        return (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

      coarsndidx = nsndidxtab[procngbnum] ++;     /* Get position of message in send array */
#ifdef SCOTCH_DEBUG_DGRAPH2
      if (coarsndidx >= coarptr->vsnddsptab[procngbtab[procngbnum] + 1]) {
        errorPrint ("dgraphCoarsenBuild: internal error (5)");
        return (1);
      }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      vsnddattab[coarsndidx].datatab[0] = ver1glbnum; /* Send fine global remote vertex number                                    */
      vsnddattab[coarsndidx].datatab[1] = multlocnum + multlocadj; /* Send coarse global multinode value; TRICK: sorted ascending */
    }
  }
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  for ( ; thrdnum < thrdmin; thrdnum ++)          /* Fill remaining thread slots (excluding idle ones) */
    coarptr->thrdtab[thrdnum].edgelocsum = edgelocsum; /* Record sum of local edges to date            */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

  if ((((finegrafptr->flagval & DGRAPHCOMMPTOP) != 0) ? dgraphCoarsenBuildPtop : dgraphCoarsenBuildColl) (coarptr) != 0)
    return (1);
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Barrier (proccomm) != MPI_SUCCESS) {    /* Make sure all-to-all communication did complete */
    errorPrint ("dgraphCoarsenBuild: communication error (1)");
    return (1);
  }

  for (vertlocnum = finegrafptr->baseval; vertlocnum < finegrafptr->vertlocnnd; vertlocnum ++) { /* Make sure all multinode data has been exchanged */
    if (coargsttax[vertlocnum] < 0) {
      errorPrint ("dgraphCoarsenBuild: invalid matching");
      return (1);
    }
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (dgraphHaloSync (finegrafptr, coargsttax + finegrafptr->baseval, GNUM_MPI) != 0) {
    errorPrint ("dgraphCoarsenBuild: cannot propagate multinode indices");
    return (1);
  }

  ercvcnttab = coargrafptr->procrcvtab;           /* TRICK: re-use some private coarse graph arrays after vertex exchange phase */
  ercvdsptab = coargrafptr->procsndtab;
  vertmltval = ((veloloctax != NULL) ? 2 : 1);
#ifdef SCOTCH_DEBUG_DGRAPH2
  vertmltval ++;                                  /* Send vertex indices as well, to check inconsistencies */
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  edgemltval = ((edloloctax != NULL) ? 2 : 1);    /* Multiplication factor for edge array size */
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  procrcvnbr = 0;                                 /* No sending neighbors yet */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
  for (procnum = 0, ercvdspidx = 0;               /* TRICK: dcntglbtab array no longer needed afterwards; can be freed */
       procnum < finegrafptr->procglbnbr; procnum ++) {
    ercvdsptab[procnum] = ercvdspidx;
    ercvcnttab[procnum] = coarptr->dcntglbtab[procnum].vertsndnbr * vertmltval +
                          coarptr->dcntglbtab[procnum].edgesndnbr * edgemltval;
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    if (coarptr->dcntglbtab[procnum].vertsndnbr > 0) { /* If process is a neighbor that actually sends data */
      ercvcnttab[procnum] += 2 * (thrdmin - 1);   /* Add space for thread index array                       */
      procrcvnbr ++;                              /* One more sending neighbor for us to receive            */
    }
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
    ercvdspidx += ercvcnttab[procnum];
  }

  memFree (coarptr->nsndidxtab);                  /* Free now useless work memory  */
  coarptr->nsndidxtab = NULL;                     /* This block won't be reclaimed */

  edgelocsiz = coarptr->edgekptnbr + coarptr->edgercvnbr - coarptr->vertrcvnbr; /* TRICK: remote edge to local vertex will always collapse */
  ercvdatsiz = coarptr->vertrcvnbr + coarptr->edgercvnbr; /* Basic size: degrees plus edge data                                            */
  esnddatsiz = coarptr->vertsndnbr + coarptr->edgesndnbr;
  if (finegrafptr->veloloctax != NULL) {          /* Add vertex loads if necessary */
    ercvdatsiz += coarptr->vertrcvnbr;
    esnddatsiz += coarptr->vertsndnbr;
  }
  if (finegrafptr->edloloctax != NULL) {          /* Add edge loads if necessary */
    ercvdatsiz += coarptr->edgercvnbr;
    esnddatsiz += coarptr->edgesndnbr;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  ercvdatsiz += coarptr->vertrcvnbr;              /* Send vertex indices, to check for inconsistencies */
  esnddatsiz += coarptr->vertsndnbr;
#endif /* SCOTCH_DEBUG_DGRAPH2 */
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  ercvdatsiz += 2 * (thrdmin - 1) * procrcvnbr;   /* Set thread slot space only for neighbors from which we receive */
  esnddatsiz += 2 * (thrdmin - 1) * finegrafptr->procngbnbr; /* Cannot compute exactly so take upper bound          */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (ercvdspidx != ercvdatsiz) {
    errorPrint ("dgraphCoarsenBuild: internal error (6)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  cheklocval = 0;
  coargrafptr->flagval = DGRAPHFREETABS | DGRAPHFREEPRIV | DGRAPHVERTGROUP; /* Coarse graph is not yet based */
  coarptr->coarprvptr = NULL;                     /* Transfer ownership of private arrays to coarse graph    */
  esndcnttab = NULL;                              /* In case of memory allocation error                      */
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  if (thrdmin > 1)                                /* If more than one thread    */
    vendlocsiz = coarptr->multlocnbr;             /* Create a non-compact graph */
  else
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
    vendlocsiz = 1;                               /* Else create a compact graph */
  if ((memAllocGroup ((void **) (void *)
                      &coargrafptr->vertloctax, (size_t) ((coarptr->multlocnbr + vendlocsiz) * sizeof (Gnum)),
                      &coargrafptr->veloloctax, (size_t) ( coarptr->multlocnbr               * sizeof (Gnum)), NULL) == NULL) ||
      ((coargrafptr->edgeloctax = memAlloc (edgelocsiz * sizeof (Gnum))) == NULL) ||
      ((coargrafptr->edloloctax = memAlloc (edgelocsiz * sizeof (Gnum))) == NULL) ||
      (memAllocGroup ((void **) (void *)
                      &esndcnttab,  (size_t) (finegrafptr->procglbnbr * sizeof (int)),
                      &esnddsptab,  (size_t) (finegrafptr->procglbnbr * sizeof (int)),
                      &esnddattab,  (size_t) (esnddatsiz * sizeof (Gnum)),
                      &ercvdattab,  (size_t) (ercvdatsiz * sizeof (Gnum)),
#ifdef SCOTCH_DEBUG_DGRAPH2
                      &ercvdbgtab,  (size_t) (finegrafptr->procglbnbr * sizeof (int)),
#endif /* SCOTCH_DEBUG_DGRAPH2 */
                      NULL) == NULL)) {
    errorPrint ("dgraphCoarsenBuild: out of memory (2)");
    cheklocval = 1;
  }
  coarptr->nsndidxtab = esndcnttab;               /* TRICK: allow data array to be released on error    */
#ifdef SCOTCH_DEBUG_DGRAPH1                       /* Communication cannot be overlapped by a useful one */
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_SUM, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuild: communication error (2)");
    chekglbval = 1;
  }
#else /* SCOTCH_DEBUG_DGRAPH1 */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_DGRAPH1 */
  if (chekglbval != 0) {
    dgraphFree (coargrafptr);
    return (1);
  }

  coargrafptr->baseval     = finegrafptr->baseval;
  coargrafptr->vertlocnnd  = coargrafptr->baseval + coargrafptr->vertlocnbr;
  coargrafptr->vertloctax -= coargrafptr->baseval;
  coargrafptr->vendloctax  = coargrafptr->vertloctax + vendlocsiz; /* Graph may or may not be compact */
  coargrafptr->veloloctax -= coargrafptr->baseval;
  coargrafptr->edgeloctax -= coargrafptr->baseval;
  coargrafptr->edloloctax -= coargrafptr->baseval;
  coargrafptr->edgelocsiz  = edgelocsiz;          /* Temporary size of edge array */
  coargrafptr->veloglbsum  = finegrafptr->veloglbsum;

  for (procngbnum = procnum = 0, esnddspidx = 0;  /* Fill-in adjacency arrays to be sent to neighbor processes */
       procngbnum < procngbnbr; procngbnum ++, procnum ++) {
    int                 procglbnum;
    int                 vrcvidxnnd;
    int                 vrcvidxnum;
    int                 esnddspbas;
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    Gnum                coarngbnbr;               /* Number of coarse vertices in neighbor process         */
    Gnum                coarngbadj;               /* Adjustment value for global coarse vertex number      */
    Gnum                coarngbnnd;               /* Boundary for current thread slot in neighbor          */
    Gnum                edgengbsum;               /* Cumulative number of edges to be sent to thread slots */
    int                 thrdnum;
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

    procglbnum = procngbtab[procngbnum];
    while (procnum < procglbnum) {                /* Fill empty slots */
      esnddsptab[procnum] = esnddspidx;
      esndcnttab[procnum] = 0;
      procnum ++;
    }
    esnddsptab[procnum] = esnddspidx;             /* Record start of send array for this (neighbor) process */

    vrcvidxnum = coarptr->vrcvdsptab[procglbnum]; /* Get boundaries of neighbor message */
    vrcvidxnnd = coarptr->nrcvidxtab[procngbnum];
    if (vrcvidxnum >= vrcvidxnnd) {               /* If no data requested by this neighbor */
      esndcnttab[procnum] = 0;                    /* Nothing to send back                  */
#ifdef SCOTCH_DEBUG_DGRAPH2
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
      esnddatsiz -= 2 * (thrdmin - 1);            /* Lower boundary by discarding thread slot index space not used by this neighbor */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      continue;                                   /* Skip to next neighbor */
    }

    esnddspbas = esnddspidx;
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    esnddspidx += 2 * (thrdmin - 1);              /* Reserve space for thread slot index */

    thrdnum = 1;
    coarngbnbr = coargrafptr->proccnttab[procglbnum]; /* Get initial value and number of coarse vertices for this neighbor */
    coarngbadj = coargrafptr->procdsptab[procglbnum];
    coarngbnnd = coarngbadj + DATASCAN (coarngbnbr, thrdmin, 1); /* Compute first thread slot boundary */
    edgengbsum = 0;
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
    do {                                          /* For all multinode requests received, in order  */
      Gnum                vertlocnum;             /* Local number of vertex whose adjacency to send */
      Gnum                edgelocnum;             /* Start index of fine adjacency array            */
      Gnum                edgelocnnd;             /* End index of fine adjacency array              */

      vertlocnum = vrcvdattab[vrcvidxnum].datatab[0] - vertlocadj; /* Get number of fine local vertex whose adjacency is wanted */
      edgelocnum = vertloctax[vertlocnum];
      edgelocnnd = vendloctax[vertlocnum];
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
      while (vrcvdattab[vrcvidxnum].datatab[1] >= coarngbnnd) { /* As long as sender not in the proper thread slot         */
        esnddattab[esnddspbas + (2 * thrdnum) - 2] = edgengbsum; /* Record number of indices for new slot                  */
        esnddattab[esnddspbas + (2 * thrdnum) - 1] = (Gnum) (esnddspidx - esnddspbas); /* Record start index for new slot  */
        thrdnum ++;
        coarngbnnd = coarngbadj + DATASCAN (coarngbnbr, thrdmin, thrdnum); /* Compute next thread slot boundary */
      }
      edgengbsum += edgelocnnd - edgelocnum - 1;  /* Account for edges to be sent in this slot; TRICK : "-1" for collapsed ghost edge */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
#ifdef SCOTCH_DEBUG_DGRAPH2
      if ((esnddspidx + vertmltval + edgemltval * (edgelocnnd - edgelocnum)) > esnddatsiz) {
        errorPrint ("dgraphCoarsenBuild: internal error (7)");
        return (1);
      }
      esnddattab[esnddspidx ++] = vertlocnum + vertlocadj; /* Write global fine vertex number; accounted for in computation just above */
#endif /* SCOTCH_DEBUG_DGRAPH2 */
      esnddattab[esnddspidx ++] = (edgelocnnd - edgelocnum); /* Write degree */
      if (veloloctax != NULL)
        esnddattab[esnddspidx ++] = veloloctax[vertlocnum];
      if (edloloctax != NULL) {                   /* If coarse graph has edge loads              */
        for ( ; edgelocnum < edgelocnnd; edgelocnum ++) { /* Send coarse adjacency and edge load */
          esnddattab[esnddspidx ++] = coargsttax[edgegsttax[edgelocnum]];
          esnddattab[esnddspidx ++] = edloloctax[edgelocnum];
        }
      }
      else {                                      /* Else only send coarse adjacency */
        for ( ; edgelocnum < edgelocnnd; edgelocnum ++)
          esnddattab[esnddspidx ++] = coargsttax[edgegsttax[edgelocnum]];
      }
    } while (++ vrcvidxnum < vrcvidxnnd);
    esndcnttab[procnum] = esnddspidx - esnddspbas; /* Record amount of data to be sent to this process */
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    for ( ; thrdnum < thrdmin; thrdnum ++) {      /* Fill remaining thread slots */
      esnddattab[esnddspbas + (2 * thrdnum) - 2] = edgengbsum;
      esnddattab[esnddspbas + (2 * thrdnum) - 1] = (Gnum) (esnddspidx - esnddspbas);
    }
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
  }
  while (procnum < finegrafptr->procglbnbr) {     /* Complete fill-in of empty slots */
    esnddsptab[procnum] = esnddspidx;
    esndcnttab[procnum] = 0;
    procnum ++;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if (esnddspidx != esnddatsiz) {                 /* Check that all expected data has been encoded */
    errorPrint ("dgraphCoarsenBuild: internal error (8)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (MPI_Alltoall (esndcnttab, 1, MPI_INT, ercvdbgtab, 1, MPI_INT, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuild: communication error (3)");
    return (1);
  }
  for (procnum = 0; procnum < finegrafptr->procglbnbr; procnum ++) {
    if (ercvdbgtab[procnum] != ercvcnttab[procnum]) {
      errorPrint ("dgraphCoarsenBuild: internal error (9)");
      return (1);
    }
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */
  if (MPI_Alltoallv (esnddattab, esndcnttab, esnddsptab, GNUM_MPI,
                     ercvdattab, ercvcnttab, ercvdsptab, GNUM_MPI, proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsenBuild: communication error (4)");
    return (1);
  }

  for (procngbnum = 0; procngbnum < procngbnbr; procngbnum ++) { /* Make receive array store only neighbor-related data */
    ercvdsptab[procngbnum] = ercvdsptab[procngbtab[procngbnum]];
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
    ercvcnttab[procngbnum] = ercvcnttab[procngbtab[procngbnum]]; /* Count array also needed for threaded version */
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
  }
  coarptr->ercvdattab = ercvdattab;
  coarptr->ercvdsptab = ercvdsptab;
#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  coarptr->ercvcnttab = ercvcnttab;
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */

  for (coarhashnbr = 32; coarhashnbr < finegrafptr->degrglbmax; coarhashnbr = coarhashnbr * 2) ; /* Compute size of adjacency hash table */
  coarptr->coarhashmsk = coarhashnbr * 4 - 1;     /* TRICK: size is power of two */

#if (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD)
  if (thrdmin > 1) {                              /* If multithreading, coarse graph shall not be compact */
    contextThreadLaunch (coarptr->contptr, (ThreadFunc) dgraphCoarsenBuildThr, (void *) coarptr);
    if (coarptr->thrdtab[0].retuval != 0) {
      errorPrint ("dgraphCoarsenBuild: could not compute adjacency (1)");
      return (1);
    }

    coargrafptr->velolocsum = coarptr->thrdtab[0].velolocsum; /* Get reduced vertex load sum       */
    coargrafptr->edgelocnbr = coarptr->thrdtab[0].edgelocnbr; /* Get reduced local number of edges */
    coargrafptr->degrglbmax = coarptr->thrdtab[0].degrlocmax; /* Get reduced maximum local degree  */
  }
  else
#endif /* (defined SCOTCH_PTHREAD) && (! defined DGRAPHCOARSENNOTHREAD) */
  {
    if (dgraphCoarsenBuildSeq (coarptr)) {
      errorPrint ("dgraphCoarsenBuild: could not compute adjacency (2)");
      return (1);
    }
  }

  coargrafptr->edgeloctax  = memRealloc (coargrafptr->edgeloctax + coargrafptr->baseval, coargrafptr->edgelocsiz * sizeof (Gnum));
  coargrafptr->edgeloctax -= coargrafptr->baseval;
  coargrafptr->edloloctax  = memRealloc (coargrafptr->edloloctax + coargrafptr->baseval, coargrafptr->edgelocsiz * sizeof (Gnum));
  coargrafptr->edloloctax -= coargrafptr->baseval;

  reduloctab[0] = coargrafptr->vertlocnbr;        /* Compute maximum of local vertices across all processes  */
  reduloctab[1] = coargrafptr->edgelocnbr;        /* Compute maximum of local edges across all processes     */
  reduloctab[2] = coargrafptr->edgelocsiz;        /* Compute maximum of edge array size across all processes */
  reduloctab[3] = coargrafptr->degrglbmax;        /* Compute maximum vertex degree across all processes      */
  reduloctab[4] = coargrafptr->edgelocnbr;        /* Compute sum of edges across all processes               */

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 4, 1, proccomm) != 0) {
    errorPrint ("dgraphCoarsenBuild: communication error (5)");
    return (1);
  }

  coargrafptr->vertglbmax = reduglbtab[0];
  coargrafptr->edgeglbmax = reduglbtab[1];
  coargrafptr->edgeglbsmx = reduglbtab[2];
  coargrafptr->degrglbmax = reduglbtab[3];        /* It is now a real global maximum degree */
  coargrafptr->edgeglbnbr = reduglbtab[4];

  memFree (coarptr->thrdtab);
  coarptr->thrdtab = NULL;

  return (0);
}

/***************************/
/*                         */
/* The coarsening routine. */
/*                         */
/***************************/

/* This routine coarsens the given fine distributed
** graph, as long as the coarsening ratio remains
** below some threshold value and the coarsened graph
** is not too small.
** If a multinode array is provided (*multlocptr != NULL),
** it must be of a size sufficient to hold multinode data
** in any configuration, including in the case of folding
** with duplication, where folded data is spread across
** floor(P/2) processes.
** It returns:
** - 0  : if the graph has been coarsened.
** - 1  : if the graph could not be coarsened.
** - 2  : on error.
*/

int
dgraphCoarsen (
Dgraph * restrict const               finegrafptr, /*+ Graph to coarsen                   +*/
Dgraph * restrict const               coargrafptr, /*+ Coarse graph to build              +*/
DgraphCoarsenMulti * restrict * const multlocptr, /*+ Pointer to un-based multinode array +*/
const Gnum                            passnbr,    /*+ Number of coarsening passes to go   +*/
const Gnum                            coarnbr,    /*+ Minimum number of coarse vertices   +*/
const double                          coarrat,    /*+ Maximum contraction ratio           +*/
const int                             flagval,    /*+ Flag value                          +*/
Context * restrict const              contptr)    /*+ Execution context                   +*/
{
  DgraphMatchData           matedat;              /* Matching state data; includes coarsening handling data   */
  Gnum                      vertrcvnbr;           /* Overall number of vertices to be received from neighbors */
  Gnum                      edgercvnbr;           /* Overall number of edges to be received from neighbors    */
  Gnum                      vertsndnbr;           /* Overall number of vertices to be sent to neighbors       */
  Gnum                      edgesndnbr;           /* Overall number of edges to be sent to neighbors          */
  int                       thrdlocmin;
  int                       thrdglbmin;
  int                       cheklocval;
  Gnum                      coarvertmax;
  Gnum                      passnum;
  int                       procnum;
  int                       o;

#ifdef SCOTCH_DEBUG_DGRAPH1
  if (coarrat < 0.5L)                             /* If impossible coarsening ratio wanted */
    return (1);                                   /* We will never succeed                 */
#endif /* SCOTCH_DEBUG_DGRAPH1 */

  coarvertmax = (Gnum) ((double) finegrafptr->vertglbnbr * coarrat); /* Maximum number of coarse vertices */
  if (coarvertmax < coarnbr)                      /* If there are too few vertices in graph               */
    return (1);                                   /* It is useless to go any further                      */

  if (dgraphGhst (finegrafptr) != 0) {            /* Compute ghost edge array of fine graph if not already present */
    errorPrint ("dgraphCoarsen: cannot compute ghost edge array");
    return (2);
  }

  matedat.c.flagval    = flagval;
  matedat.c.multloctab = *multlocptr;             /* Propagate the provided multinode array or NULL if it has to be allocated */
  matedat.c.contptr    = contptr;                 /* Set execution context as it might be needed for initialization           */
  cheklocval  = dgraphCoarsenInit (&matedat.c, finegrafptr, coargrafptr);
  cheklocval |= dgraphMatchInit   (&matedat, 0.5F);

  thrdlocmin = (cheklocval != 0) ? -1 : contextThreadNbr (contptr); /* TRICK: negative value on error */
  if (MPI_Allreduce (&thrdlocmin, &thrdglbmin, 1, MPI_INT, MPI_MIN, finegrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsen: communication error (1)");
    return (2);
  }
  if (thrdglbmin < 0) {                           /* On error */
    dgraphMatchExit   (&matedat);
    dgraphCoarsenExit (&matedat.c);
    return (2);
  }
  matedat.c.thrdmin = thrdglbmin;                 /* Minimum number of threads across all processes */

  for (passnum = 0; passnum < passnbr; passnum ++) {
    ((passnum == 0) ? dgraphMatchHl : dgraphMatchHy) (&matedat); /* If first pass, process lightest vertices first */

    if ((((finegrafptr->flagval & DGRAPHCOMMPTOP) != 0) ? dgraphMatchSyncPtop : dgraphMatchSyncColl) (&matedat) != 0) {
      errorPrint        ("dgraphCoarsen: cannot perform matching");
      dgraphMatchExit   (&matedat);
      dgraphCoarsenExit (&matedat.c);
      return (2);
    }
  }
  dgraphMatchLy (&matedat);                       /* All remaining vertices are matched locally */

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphMatchCheck (&matedat) != 0) {
    errorPrint        ("dgraphCoarsen: invalid matching");
    dgraphMatchExit   (&matedat);
    dgraphCoarsenExit (&matedat.c);
    return (2);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  dgraphMatchExit (&matedat);

  vertsndnbr =
  edgesndnbr = 0;
  for (procnum = 0; procnum < finegrafptr->procglbnbr; procnum ++) {
    vertsndnbr += matedat.c.dcntloctab[procnum].vertsndnbr;
    edgesndnbr += matedat.c.dcntloctab[procnum].edgesndnbr;
    matedat.c.dcntloctab[procnum].vertlocnbr = matedat.c.multlocnbr;
  }
  matedat.c.vertsndnbr = vertsndnbr;
  matedat.c.edgesndnbr = edgesndnbr;

  if (MPI_Alltoall (matedat.c.dcntloctab, 3, GNUM_MPI, matedat.c.dcntglbtab, 3, GNUM_MPI, finegrafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphCoarsen: communication error (2)");
    return (2);
  }

  vertrcvnbr =
  edgercvnbr = 0;
  coargrafptr->procdsptab[0] = finegrafptr->baseval; /* Build vertex-to-process array */
  for (procnum = 0; procnum < finegrafptr->procglbnbr; procnum ++) {
    Gnum                proccntval;

    vertrcvnbr += matedat.c.dcntglbtab[procnum].vertsndnbr;
    edgercvnbr += matedat.c.dcntglbtab[procnum].edgesndnbr;
    proccntval  = matedat.c.dcntglbtab[procnum].vertlocnbr;
    coargrafptr->proccnttab[procnum] = proccntval;
    coargrafptr->procdsptab[procnum + 1] = coargrafptr->procdsptab[procnum] + proccntval;
  }
  coargrafptr->vertlocnbr = matedat.c.multlocnbr;
  coargrafptr->vertglbnbr = coargrafptr->procdsptab[finegrafptr->procglbnbr] - finegrafptr->baseval;
  matedat.c.vertrcvnbr = vertrcvnbr;
  matedat.c.edgercvnbr = edgercvnbr;

  if (coargrafptr->vertglbnbr > coarvertmax) {    /* If coarsening ratio not met */
    dgraphCoarsenExit (&matedat.c);
    return (1);
  }

  if (matedat.c.multloctmp != NULL) {             /* If we allocated the multinode array */
    matedat.c.multloctmp =
    matedat.c.multloctab = memRealloc (matedat.c.multloctab, matedat.c.multlocnbr * sizeof (DgraphCoarsenMulti)); /* Resize multinode array */
  }

  if (dgraphCoarsenBuild (&matedat.c) != 0) {     /* Build coarse graph */
    dgraphCoarsenExit (&matedat.c);
    return (2);
  }

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (dgraphCheck (coargrafptr) != 0) {           /* Check graph consistency */
    errorPrint ("dgraphCoarsen: inconsistent graph data");
    dgraphFree (coargrafptr);
    return (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  matedat.c.multloctmp = NULL;                    /* So that it will not be freed    */
  dgraphCoarsenExit (&matedat.c);                 /* Free all other temporary arrays */

  o = 0;                                          /* Assume everything is now all right   */
  if (((flagval & DGRAPHCOARSENFOLDDUP) != 0) &&  /* If some form of folding is requested */
      (coargrafptr->procglbnbr >= 2)) {           /* And if there is need to it           */
    Dgraph                coargrafdat;            /* Coarse graph data before folding     */
    DgraphCoarsenMulti *  coarmultptr;            /* Pointer to folded multinode array    */
    MPI_Datatype          coarmultype;

    MPI_Type_contiguous (2, GNUM_MPI, &coarmultype); /* Define type for MPI transfer */
    MPI_Type_commit (&coarmultype);               /* Commit new type                 */

    coargrafdat = *coargrafptr;                   /* Copy unfolded coarse graph data to save area */
    coarmultptr = NULL;                           /* Assume we will not get a multinode array     */
    if ((flagval & DGRAPHCOARSENFOLDDUP) == DGRAPHCOARSENFOLD) { /* Do a simple folding           */
      memSet (coargrafptr, 0, sizeof (Dgraph));   /* Also reset procglbnbr for unused processes   */
      o = dgraphFold (&coargrafdat, 0, coargrafptr, (void *) matedat.c.multloctab, (void **) (void *) &coarmultptr, coarmultype);
    }
    else {                                        /* Do a duplicant-folding */
      int               loopval;

      o = dgraphFoldDup (&coargrafdat, coargrafptr, (void *) matedat.c.multloctab, (void **) (void *) &coarmultptr, coarmultype, contptr);
      loopval = contextIntRandVal (contptr, finegrafptr->proclocnum + contextIntRandVal (contptr, finegrafptr->proclocnum * 2 + 1) + 1);
      while (loopval --)                          /* Desynchronize pseudo-random generator across processes */
        contextIntRandVal (contptr, 2);
    }
    dgraphExit    (&coargrafdat);                 /* Free unfolded graph */
    MPI_Type_free (&coarmultype);
    if (*multlocptr == NULL)                      /* If unfolded multinode array was not user-provided, free it */
      memFree (matedat.c.multloctab);
    *multlocptr = coarmultptr;                    /* Return folded multinode array or NULL */
  }
  else                                            /* No folding at all                                                   */
    *multlocptr = matedat.c.multloctab;           /* Return un-based pointer (maybe the same as initially user-provided) */

  return (o);
}
