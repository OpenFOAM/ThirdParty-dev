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
/**   NAME       : dmesh_io_load.c                         **/
/**                                                        **/
/**   AUTHOR     : Marc FUENTES                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module loads an auxiliary          **/
/**                distributed mesh.                       **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 11 jan 2024     **/
/**                                 to   : 09 aug 2025     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define SCOTCH_DMESH_IO_LOAD

#include "module.h"
#include "common.h"
#include "graph.h"
#include "dgraph_allreduce.h"
#include "dmesh.h"
#include "dmesh_io_load.h"

/* This routine loads a distributed source
** graph from the given stream(s). Either
** one processor holds a non-NULL stream
** of a centralized graph, or all of them
** hold valid streams to either a centralized
** or a distributed graph.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

DGRAPHALLREDUCEMAXSUMOP (6, 3)
DGRAPHALLREDUCEMAXSUMOP (12, 2)

int
dmeshLoad (
Dmesh * restrict const      meshptr,              /* Distributed mesh structure to fill         */
FILE * const                fileptr,              /* Local distributed stream to read data from */
const Gnum                  baseval,              /* Base value (-1 means keep file base)       */
const GraphLoadFlag         flagval)              /* Mesh loading flags (unused)                */
{
  Gnum                reduloctab[9];
  Gnum                reduglbtab[9];
  Gnum                versval;

#ifdef SCOTCH_DEBUG_DMESH2
  if (MPI_Barrier (meshptr->proccomm) != MPI_SUCCESS) { /* Synchronize for debugging */
    errorPrint ("dmeshLoad: communication error (1)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DMESH2 */

  reduloctab[0] = baseval;                        /* Exchange baseval to check it is the same for all */
  reduloctab[1] = - baseval;
  reduloctab[2] = (Gnum) flagval;                 /* Exchange flagval to check it is the same for all */
  reduloctab[3] = - (Gnum) flagval;
  reduloctab[4] = 1;                              /* Set uneffective values for versval */
  reduloctab[5] = -3;
  reduloctab[6] =                                 /* Assume everything will be fine */
  reduloctab[7] =                                 /* Assume does not have a stream  */
  reduloctab[8] = 0;
  if (fileptr != NULL) {
    if (intLoad (fileptr, &versval) != 1) {       /* Read version number */
      errorPrint ("dmeshLoad: bad input");
      versval       = 0;
      reduloctab[6] = 1;
    }
    else if ((versval != 1) && (versval != 3)) {  /* If not a mesh format */
      errorPrint ("dmeshLoad: not a mesh format");
      reduloctab[6] = 1;
    }
    reduloctab[4] = versval;
    reduloctab[5] = - versval;
    reduloctab[7] = 1;                            /* One more process involved in loading */
    reduloctab[8] = meshptr->proclocnum;
  }

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 6, 3, meshptr->proccomm) != 0) {
    errorPrint ("dmeshLoad: communication error (2)");
    return (1);
  }

  if (reduglbtab[6] != 0)                         /* Return from previous errors */
    return (1);

  if ((reduglbtab[0] != - reduglbtab[1])) {
    errorPrint ("dmeshLoad: inconsistent base value");
    return (1);
  }
  if ((reduglbtab[2] != - reduglbtab[3])) {
    errorPrint ("dmeshLoad: inconsistent flag value");
    return (1);
  }
  if ((reduglbtab[7] != 0) &&
      (reduglbtab[4] != - reduglbtab[5])) {
    errorPrint ("dmeshLoad: inconsistent mesh file version value");
    return (1);
  }

  if (reduglbtab[4] == 3) {                       /* If distributed mesh format             */
    if (reduglbtab[7] == meshptr->procglbnbr)     /* If as many input streams as processors */
      return (dmeshLoadAdm (meshptr, fileptr, baseval, flagval)); /* Read distributed mesh  */
  }
  else {                                          /* If centralized mesh format */
    if (reduglbtab[7] == 1)                       /* If only one reader stream  */
      return (dmeshLoadCent (meshptr, fileptr, baseval, flagval, reduglbtab[8])); /* Distribute centralized mesh from known root */
    else if (reduglbtab[7] == meshptr->procglbnbr)
      return (dmeshLoadMulti (meshptr, fileptr, baseval, flagval)); /* Read multi-centralized mesh */
  }

  errorPrint ((reduglbtab[7] == 0)
              ? "dmeshLoad: no input stream provided"
              : "dmeshLoad: invalid number of input streams");
  return (1);
}

/* This routine loads a centralized source
** mesh from a single stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
dmeshLoadCent (
Dmesh * restrict const      meshptr,              /* Distributed mesh to load             */
FILE * const                fileptr,              /* One single centralized stream        */
Gnum                        baseval,              /* Base value (-1 means keep file base) */
const DmeshFlag             flagval,              /* Mesh loading flags                   */
const int                   protnum)              /* Root process number                  */
{
  errorPrint ("dmeshLoadCent: not implemented");
  return (1);
}

/* This routine loads a distributed source
** mesh from a distributed source auxiliary
** distributed mesh file spread across all
** of the streams.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
dmeshLoadAdm (
Dmesh * restrict const      meshptr,              /* Distributed mesh structure to fill         */
FILE * const                fileptr,              /* Local distributed stream to read data from */
Gnum                        baseval,              /* Base value (-1 means keep file base)       */
const DmeshFlag             flagval)              /* Mesh loading flags                         */
{
  Gnum                  baseadj;
  Gnum                  proclocnum;
  Gnum                  vnodglbnbr;               /* Global number of node vertices                           */
  Gnum                  vnodglbnnd;
  Gnum                  velmglbnbr;               /* Global number of element vertices                        */
  Gnum * restrict       velmloctab;               /* Un-based pointer to element vertex index array           */
  Gnum                  velmlocnbr;
  Gnum                  velmlocnum;
  Gnum * restrict       eelmloctax;               /* Based pointer to element edge array to node vertices     */
  Gnum                  eelmlocnbr;
  Gnum                  eelmlocnnd;
  Gnum                  eelmlocnum;
  Gnum * restrict       prelvrttab;
  Gnum                  preldspadj;               /* Adjustment value for element-to-process array            */
  Gnum                  reduloctab[14];
  Gnum                  reduglbtab[14];
  int                   procnum;
  int                   o;

  reduloctab[0] = 0;                              /* Assume everything will be fine */
  o  = intLoad (fileptr, &reduloctab[1]);         /* Read procglbnbr                */
  o += intLoad (fileptr, &proclocnum);            /* Read proclocnum                */
  o += intLoad (fileptr, &reduloctab[3]);         /* Read velmglbnbr                */
  o += intLoad (fileptr, &reduloctab[5]);         /* Read vnodglbnbr                */
  o += intLoad (fileptr, &reduloctab[10]);        /* Read velmlocnbr                */
  o += intLoad (fileptr, &reduloctab[11]);        /* Read eelmlocnbr                */
  o += intLoad (fileptr, &reduloctab[7]);         /* Read baseval                   */
  o += intLoad (fileptr, &reduloctab[9]);         /* Read flagval; not used yet     */
  if ((o != 8)            ||
      (reduloctab[9] < 0) ||
      (reduloctab[9] > 111)) {
    errorPrint ("dmeshLoadAdm: bad input (1)");
    reduloctab[0] = 2;                            /* Immediate abort has maximum value so as to be propagated by MAX reduce */
  }
  reduloctab[2]  = - reduloctab[1];
  reduloctab[4]  = - reduloctab[3];
  reduloctab[6]  = - reduloctab[5];
  reduloctab[8]  = - reduloctab[7];
  reduloctab[12] =   reduloctab[10];              /* Compute sums of local numbers */
  reduloctab[13] =   reduloctab[11];

  if ((int) proclocnum != meshptr->proclocnum) {  /* If fragment is not read by proper process */
    errorPrint ("dmeshLoadAdm: wrong process number to read auxiliary distributed mesh fragment");
    reduloctab[0] |= 1;
  }
  if ((int) reduloctab[1] != meshptr->procglbnbr) {
    errorPrint ("dmeshLoadAdm: wrong number of processes to read auxiliary distributed mesh");
    reduloctab[0] = 2;
  }

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 12, 2, meshptr->proccomm) != 0) {
    errorPrint ("dmeshLoadAdm: communication error (1)");
    reduglbtab[0] = 2;
  }
  if (reduglbtab[0] >= 2)                         /* If has to abort immediately */
    return (1);

  if ((reduglbtab[2]  != - reduglbtab[1]) ||
      (reduglbtab[4]  != - reduglbtab[3]) ||
      (reduglbtab[6]  != - reduglbtab[5]) ||
      (reduglbtab[8]  != - reduglbtab[7]) ||
      (reduglbtab[12] !=   reduloctab[3])) {
    errorPrint ("dmeshLoadAdm: inconsistent distributed mesh headers");
    return (1);
  }

  if (baseval == -1) {                            /* If keep file graph base     */
    baseval = reduglbtab[7];                      /* Set graph base as file base */
    baseadj = 0;                                  /* No base adjustment needed   */
  }
  else                                            /* If set graph base  */
    baseadj = baseval - reduglbtab[7];            /* Update base adjust */

  velmglbnbr = reduloctab[3];
  vnodglbnbr = reduloctab[5];
  velmlocnbr = reduloctab[10];
  eelmlocnbr = reduloctab[11];

  if ((prelvrttab = memAlloc ((meshptr->procglbnbr + 1) * sizeof (Gnum))) == NULL) {
    errorPrint ("dmeshLoadAdm: out of memory (1)");
    reduloctab[0] = 1;
  }
  if (MPI_Allreduce (&reduloctab[0], &reduglbtab[0], 1, GNUM_MPI, MPI_MAX, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshLoadAdm: communication error (2)");
    return (1);
  }
  if (reduglbtab[0] == 1) {
    if (prelvrttab != NULL)
      memFree (prelvrttab);
    return (1);
  }

  if ((velmloctab = memAlloc ((velmlocnbr + 1) * sizeof (Gnum))) == NULL) { /* Vertex array will remain unbased */
    errorPrint ("dmeshLoadAdm: out of memory (2)");
    reduloctab[10] = -1;                          /* TRICK: update velmlocnbr to be broadcast */
  }
  if ((eelmloctax = memAlloc (eelmlocnbr * sizeof (Gnum))) == NULL) { /* Edge array will be based */
    errorPrint ("dmeshLoadAdm: out of memory (3)");
    reduloctab[10] = -1;
  }

  if (MPI_Allgather (&reduloctab[10], 1, GNUM_MPI,
                     prelvrttab,      1, GNUM_MPI, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshLoadAdm: communication error (3)");
    return (1);
  }

  preldspadj = baseval;                           /* Use the new base to build displacement array  */
  for (procnum = 0; procnum < meshptr->procglbnbr; procnum ++) { /* Build element-to-process array */
    Gnum                preldspval;

    preldspval = prelvrttab[procnum];
    if (preldspval < 0) {                         /* If error notified by us or by another process */
      if (eelmloctax != NULL)
        memFree (eelmloctax);                     /* Array is not yet based */
      if (velmloctab != NULL)
        memFree (velmloctab);
      memFree (prelvrttab);
      return  (1);
    }

    prelvrttab[procnum] = preldspadj;
    preldspadj += preldspval;
  }
  prelvrttab[procnum] = preldspadj;               /* Set end of displacement array */
  if (preldspadj != (velmglbnbr + baseval)) {
    errorPrint ("dmeshLoadAdm: bad input (2)");
    reduloctab[0] = 1;
    goto abort;
  }

  vnodglbnnd  = vnodglbnbr + baseval;
  eelmloctax -= baseval;                          /* Base edge array with new base */
  eelmlocnnd  = eelmlocnbr + baseval;
  for (velmlocnum = 0, eelmlocnum = baseval; velmlocnum < velmlocnbr; velmlocnum ++) {
    Gnum                degrval;                  /* Number of neighbors for element */

    velmloctab[velmlocnum] = eelmlocnum;          /* Set start of neighbor list for new vertex */

    if ((intLoad (fileptr, &degrval) != 1) ||
        (degrval < 0)) {
      errorPrint ("dmeshLoadAdm: bad input (3)");
      reduloctab[0] = 1;
      goto abort;
    }
    while (degrval -- > 0) {                      /* For all neighbors              */
      if ((eelmlocnum >= eelmlocnnd)                        || /* If out of bounds  */
          (intLoad (fileptr, &eelmloctax[eelmlocnum]) != 1) ||
          (eelmloctax[eelmlocnum] += baseadj,     /* Adjust end element vertex value */
           eelmloctax[eelmlocnum] <  baseval)               ||
          (eelmloctax[eelmlocnum] >= vnodglbnnd)) {
        errorPrint ("dmeshLoadAdm: bad input (4)");
        reduloctab[0] = 1;
        goto abort;
      }

      eelmlocnum ++;                              /* One more edge recorded */
    }
  }
  velmloctab[velmlocnum] = eelmlocnum;            /* Set end of vertex array */
  if (eelmlocnum != eelmlocnnd) {
    errorPrint ("dmeshLoadAdm: bad input (5)");
    reduloctab[0] = 1;
  }

abort:
  if (MPI_Allreduce (&reduloctab[0], &reduglbtab[0], 1, GNUM_MPI, MPI_MAX, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshLoadAdm: communication error (4)");
    return (1);
  }
  if (reduglbtab[0] != 0) {
    memFree (eelmloctax + baseval);
    memFree (velmloctab);
    memFree (prelvrttab);
    return (1);
  }

  meshptr->flagval   |= DMESHFREETABS;
  meshptr->baseval    = baseval;
  meshptr->velmglbnbr = reduglbtab[3];
  meshptr->velmglbmax = reduglbtab[10];
  meshptr->velmlocnbr = reduloctab[10];
  meshptr->velmloctab = velmloctab;
  meshptr->eelmglbnbr = reduglbtab[13];
  meshptr->eelmglbmax = reduglbtab[11];
  meshptr->eelmlocnbr = reduloctab[11];
  meshptr->eelmloctab = eelmloctax + baseval;
  meshptr->vnodglbnbr = reduloctab[5];
  meshptr->prelvrttab = prelvrttab;

  return (0);
}

/* This routine loads a distributed source
** mesh from a centralized source mesh
** file replicated on all of the streams.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
int
dmeshLoadMulti (
Dmesh * restrict const      meshptr,              /* Distributed mesh to load             */
FILE * const                fileptr,              /* Duplicated centralized streams       */
Gnum                        baseval,              /* Base value (-1 means keep file base) */
const DmeshFlag             flagval)              /* Mesh loading flags                   */
{
  errorPrint ("dmeshLoadMulti: not implemented");
  return (1);
}
