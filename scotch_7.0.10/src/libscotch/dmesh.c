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
/**   NAME       : dmesh.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the distributed    **/
/**                mesh data structure handling routines.  **/
/**                                                        **/
/**    DATES     : # Version 7.0  : from : 09 aug 2025     **/
/**                                 to   : 09 aug 2025     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "dmesh.h"

/*************************************/
/*                                   */
/* These routines handle distributed */
/* source meshes.                    */
/*                                   */
/*************************************/

/* This routine initializes a distributed mesh
** structure. In order to avoid collective
** communication whenever possible, the allocation
** of send and receive index arrays is not performed
** in the routine itself, but rather delegated to
** subsequent routines such as dmeshBuild.
** However, these arrays will not be freed by
** dmeshFree, but by dmeshExit.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dmeshInit (
Dmesh * restrict const      meshptr,              /* Distributed mesh structure                     */
MPI_Comm                    proccomm)             /* Communicator to be used for all communications */
{
  memSet (meshptr, 0, sizeof (Dmesh));            /* Clear public and private mesh fields */

  meshptr->proccomm = proccomm;                   /* Set private fields    */
  MPI_Comm_size (proccomm, &meshptr->procglbnbr); /* Get communicator data */
  MPI_Comm_rank (proccomm, &meshptr->proclocnum);

  return (0);
}

/* This routine frees the public and private data
** of the given distributed mesh, but not its
** communicator.
** Private data could have been kept and freed only
** in dmeshExit(). Yet, freeing it along with the
** mesh public data is a way to avoid memory
** fragmentation.
** It is not a collective routine, as no communication
** is needed to perform the freeing of memory structures.
** It returns:
** - VOID  : in all cases.
*/

static
void
dmeshFree2 (
Dmesh * restrict const      meshptr)
{
  if ((meshptr->flagval & DMESHFREEPRIV) != 0) {
    if (meshptr->prelvrttab != NULL)
      memFree (meshptr->prelvrttab);
  }
  if ((meshptr->flagval & DMESHFREETABS) != 0) { /* If local arrays must be freed */
    if (meshptr->velmloctab != NULL)
      memFree (meshptr->velmloctab);
    if (meshptr->eelmloctab != NULL)
      memFree (meshptr->eelmloctab);
  }
}

void
dmeshFree (
Dmesh * restrict const      meshptr)
{
  DmeshFlag           flagval;
  MPI_Comm            proccomm;                   /* Data for temporarily saving private data */
  int                 procglbnbr;
  int                 proclocnum;

  dmeshFree2 (meshptr);                           /* Free all user fields */

  flagval    = meshptr->flagval & DMESHFREECOMM;
  proccomm   = meshptr->proccomm;                 /* Save private fields only */
  procglbnbr = meshptr->procglbnbr;
  proclocnum = meshptr->proclocnum;

  memSet (meshptr, 0, sizeof (Dmesh));           /* Reset mesh structure */

  meshptr->flagval    = flagval;                  /* Restore private fields */
  meshptr->proccomm   = proccomm;
  meshptr->procglbnbr = procglbnbr;
  meshptr->proclocnum = proclocnum;

  return;
}

/* This routine destroys a distributed mesh structure.
** It is not a collective routine, as no communication
** is needed to perform the freeing of memory structures.
** Private data are always destroyed. If this is not
** wanted, use dmeshFree() instead.
** It returns:
** - VOID  : in all cases.
*/

void
dmeshExit (
Dmesh * restrict const     meshptr)
{
  DmeshFlag          flagval;

  flagval = meshptr->flagval;
  if ((flagval & DMESHFREECOMM) != 0)            /* If communicator has to be freed */
    MPI_Comm_free (&meshptr->proccomm);           /* Free it                         */

  dmeshFree2 (meshptr);

#ifdef SCOTCH_DEBUG_DMESH1
  memSet (meshptr, ~0, sizeof (Dmesh));
#endif /* SCOTCH_DEBUG_DMESH1 */
  meshptr->flagval = flagval & ~DMESHBITSUSED;   /* A subsequent dmeshExit() will have no effect */
}
