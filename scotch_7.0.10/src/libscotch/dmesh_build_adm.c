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
/**   NAME       : dmesh_build_adm.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the distributed    **/
/**                mesh data structure handling routines.  **/
/**                                                        **/
/**    DATES     : # Version 7.0  : from : 15 aug 2025     **/
/**                                 to   : 15 aug 2025     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "dmesh.h"
#include "dgraph_allreduce.h"

/*************************************/
/*                                   */
/* These routines handle distributed */
/* source meshes.                    */
/*                                   */
/*************************************/

/* This routine builds an auxiliary distributed
** mesh from the local arrays that are passed to
** it.
** As for all routines that build meshes, the
** private fields of the Dmesh structure have to
** be initialized if they are not already.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

DGRAPHALLREDUCEMAXSUMOP (7, 2)

int
dmeshBuildAdm (
Dmesh * restrict const      meshptr,              /* Distributed mesh structure            */
const Gnum                  baseval,              /* Base for indexing                     */
const Gnum                  velmlocnbr,           /* Number of local elements              */
Gnum * const                velmloctab,           /* Local element vertex begin array      */
const Gnum                  eelmlocnbr,           /* Number of local element-to-node edges */
Gnum * const                eelmloctab,           /* Local edge array                      */
const Gnum                  vnodglbnbr)           /* Global number of node vertices        */
{
  Gnum                preldspadj;                 /* Adjustment value for element-to-process array */
  Gnum                reduloctab[9];
  Gnum                reduglbtab[9];
  int                  procnum;

  reduloctab[0] = baseval;                        /* Exchange baseval to check it is the same for all */
  reduloctab[1] = - baseval;
  reduloctab[2] = vnodglbnbr;                     /* Exchange baseval to check it is the same for all */
  reduloctab[3] = - vnodglbnbr;
  reduloctab[4] = velmlocnbr;
  reduloctab[5] = eelmlocnbr;
  reduloctab[6] = 0;                              /* Assume everything will be fine */
  reduloctab[7] = velmlocnbr;
  reduloctab[8] = eelmlocnbr;

  if (meshptr->prelvrttab == NULL) {
    if ((meshptr->prelvrttab = memAlloc ((meshptr->procglbnbr + 1) * sizeof (Gnum))) == NULL) {
      errorPrint ("dmeshBuildAdm: out of memory");
      reduloctab[6] = 1;
    }
    else
      meshptr->flagval |= DMESHFREEPRIV;
  }

  if (dgraphAllreduceMaxSum (reduloctab, reduglbtab, 7, 2, meshptr->proccomm) != 0) {
    errorPrint ("dmeshBuildAdm: communication error (1)");
    return (1);
  }
  if (reduglbtab[6] != 0)
    return (1);

  if ((reduglbtab[1]  != - reduglbtab[0]) ||
      (reduglbtab[3]  != - reduglbtab[2])) {
    errorPrint ("dmeshBuildAdm: inconsistent parameters (1)");
    return (1);
  }

  if (MPI_Allgather (&reduloctab[4],      1, GNUM_MPI,
                     meshptr->prelvrttab, 1, GNUM_MPI, meshptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dmeshBuildAdm: communication error (2)");
    return (1);
  }

  preldspadj = baseval;                           /* Use the new base to build displacement array  */
  for (procnum = 0; procnum < meshptr->procglbnbr; procnum ++) { /* Build element-to-process array */
    Gnum                preldspval;

    preldspval = meshptr->prelvrttab[procnum];
    meshptr->prelvrttab[procnum] = preldspadj;
    preldspadj += preldspval;
  }
  meshptr->prelvrttab[procnum] = preldspadj;      /* Set end of displacement array */
  if (preldspadj != (reduglbtab[7] + baseval)) {
    errorPrint ("dmeshBuildAdm: inconsistent parameters (2)");
    return (1);
  }

  meshptr->baseval    = baseval;
  meshptr->velmglbnbr = reduglbtab[7];
  meshptr->velmglbmax = reduglbtab[4];
  meshptr->velmlocnbr = velmlocnbr;
  meshptr->velmloctab = velmloctab;
  meshptr->eelmglbnbr = reduglbtab[8];
  meshptr->eelmglbmax = reduglbtab[5];
  meshptr->eelmlocnbr = eelmlocnbr;
  meshptr->eelmloctab = eelmloctab;
  meshptr->vnodglbnbr = vnodglbnbr;

  return (0);
}
