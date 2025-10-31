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
/**   NAME       : library_dmesh.c                         **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted source mesh handling routines of  **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 12 aug 2025     **/
/**                                 to   : 12 aug 2025     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "context.h"
#include "graph.h"                                /* For graphPtscotch() */
#include "dmesh.h"
#include "ptscotch.h"

/****************************************/
/*                                      */
/* These routines are the C API for the */
/* distributed mesh handling routines.  */
/*                                      */
/****************************************/

/*+ This routine reserves a memory area
*** of a size sufficient to store a
*** SCOTCH_Dmesh structure.
*** It returns:
*** - !NULL  : if the allocation succeeded.
*** - NULL   : on error.
+*/

SCOTCH_Dmesh *
SCOTCH_dmeshAlloc ()
{
  return ((SCOTCH_Dmesh *) memAlloc (sizeof (SCOTCH_Dmesh)));
}

/*+ This routine returns the size, in bytes,
*** of a SCOTCH_Dmesh structure.
*** It returns:
*** - > 0  : in all cases.
+*/

int
SCOTCH_dmeshSizeof ()
{
  return (sizeof (SCOTCH_Dmesh));
}

/*+ This routine initializes the opaque
*** distributed mesh structure used to
*** handle distributed meshs in the
*** Scotch library.
*** It returns:
*** - 0   : if the initialization succeeded.
*** - !0  : on error.
+*/

int
SCOTCH_dmeshInit (
SCOTCH_Dmesh * const        meshptr,
MPI_Comm                    proccomm)             /* Communicator to be used for all communications */
{
#ifdef SCOTCH_DEBUG_DMESH2
  if (sizeof (SCOTCH_Num) != sizeof (Gnum)) {
    errorPrint (STRINGIFY (SCOTCH_dmeshInit) ": internal error (1)");
    return (1);
  }
  if (sizeof (SCOTCH_Dmesh) < sizeof (Dmesh)) {
    errorPrint (STRINGIFY (SCOTCH_dmeshInit) ": internal error (2)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DMESH2 */

  return (dmeshInit ((Dmesh *) meshptr, proccomm));
}

/*+ This routine frees the contents of the
*** given opaque mesh structure.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dmeshExit (
SCOTCH_Dmesh * const        meshptr)
{
  if (! contextContainerTrue (meshptr))
    dmeshExit ((Dmesh *) meshptr);
}

/*+ This routine frees the contents of the
*** given opaque mesh structure but does
*** not free its private data.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dmeshFree (
SCOTCH_Dmesh * const        meshptr)
{
  if (! contextContainerTrue (meshptr))
    dmeshFree ((Dmesh *) meshptr);
}

/*+ This routine accesses mesh size data.
*** NULL pointers on input indicate unwanted
*** data.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dmeshSize (
const SCOTCH_Dmesh * const  libmeshptr,
SCOTCH_Num * const          velmglbnbr,
SCOTCH_Num * const          velmlocnbr,
SCOTCH_Num * const          eelmglbnbr,
SCOTCH_Num * const          eelmlocnbr,
SCOTCH_Num * const          vnodglbnbr)
{
  const Dmesh * const      srcmeshptr = (Dmesh *) CONTEXTOBJECT (libmeshptr);

  if (velmglbnbr != NULL)
    *velmglbnbr = (SCOTCH_Num) (srcmeshptr->velmglbnbr);
  if (velmlocnbr != NULL)
    *velmlocnbr = (SCOTCH_Num) (srcmeshptr->velmlocnbr);
  if (eelmglbnbr != NULL)
    *eelmglbnbr = (SCOTCH_Num) (srcmeshptr->eelmglbnbr);
  if (eelmlocnbr != NULL)
    *eelmlocnbr = (SCOTCH_Num) (srcmeshptr->eelmlocnbr);
  if (vnodglbnbr != NULL)
    *vnodglbnbr = (SCOTCH_Num) (srcmeshptr->vnodglbnbr);
}

/*+ This routine accesses all of the graph data.
*** NULL pointers on input indicate unwanted
*** data. NULL pointers on output indicate
*** unexisting arrays.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_dmeshData (
const SCOTCH_Dmesh * const  libmeshptr,           /* Graph structure to read             */
SCOTCH_Num * const          baseptr,              /* Base value                          */
SCOTCH_Num * const          velmglbnbr,           /* Number of global element vertices   */
SCOTCH_Num * const          velmlocnbr,           /* Number of local element vertices    */
SCOTCH_Num ** const         velmloctab,           /* Element vertex array [velmlocnbr+1] */
SCOTCH_Num * const          eelmglbnbr,           /* Number of global element edges      */
SCOTCH_Num * const          eelmlocnbr,           /* Number of local element edges       */
SCOTCH_Num ** const         eelmloctab,           /* Element edge array [eelmlocnbr]     */
SCOTCH_Num * const          vnodglbnbr,           /* Number of global node vertices      */
MPI_Comm * const            commptr)              /* MPI Communicator                    */
{
  const Dmesh * const      srcmeshptr = (Dmesh *) CONTEXTOBJECT (libmeshptr);

  if (baseptr != NULL)
    *baseptr = srcmeshptr->baseval;
  if (velmglbnbr != NULL)
    *velmglbnbr = (SCOTCH_Num) (srcmeshptr->velmglbnbr);
  if (velmlocnbr != NULL)
    *velmlocnbr = (SCOTCH_Num) (srcmeshptr->velmlocnbr);
  if (velmloctab != NULL)
    *velmloctab = (SCOTCH_Num *) (srcmeshptr->velmloctab);
  if (eelmglbnbr != NULL)
    *eelmglbnbr = (SCOTCH_Num) (srcmeshptr->eelmglbnbr);
  if (eelmlocnbr != NULL)
    *eelmlocnbr = (SCOTCH_Num) (srcmeshptr->eelmlocnbr);
  if (eelmloctab != NULL)
    *eelmloctab = (SCOTCH_Num *) (srcmeshptr->eelmloctab);
  if (vnodglbnbr != NULL)
    *vnodglbnbr = (SCOTCH_Num) (srcmeshptr->vnodglbnbr);
  if (commptr != NULL)
    *commptr = srcmeshptr->proccomm;
}
