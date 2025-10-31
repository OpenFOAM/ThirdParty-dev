/* Copyright 2011,2012,2014,2018-2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_coarsen.c                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the          **/
/**                distributed graph coarsening routine of **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 07 aug 2011     **/
/**                                 to   : 07 aug 2011     **/
/**                # Version 6.0  : from : 11 sep 2012     **/
/**                                 to   : 25 apr 2018     **/
/**                # Version 7.0  : from : 27 aug 2019     **/
/**                                 to   : 12 aug 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "context.h"
#include "dgraph.h"
#include "dgraph_fold_comm.h"
#include "dgraph_coarsen.h"
#include "ptscotch.h"

/************************************/
/*                                  */
/* This routine is the C API for    */
/* the distributed graph coarsening */
/* routine.                         */
/*                                  */
/************************************/

/*+ This routine provides an upper bound of the
*** number of local vertices, on the calling
*** process, of the coarse graph that would be
*** created by calling SCOTCH_dgraphCoarsen().
*** It is notably used to allocate the external
*** multinode array.
*** It returns:
*** - >= 0  : in all cases.
+*/

SCOTCH_Num
SCOTCH_dgraphCoarsenVertLocMax (
const SCOTCH_Dgraph * restrict const  flibgrafptr, /* Fine graph structure */
const SCOTCH_Num                      flagval)    /* Coarsening type       */
{
  CONTEXTDECL                  (flibgrafptr);

  Dgraph * restrict const   finegrafptr = (Dgraph *) CONTEXTOBJECT (flibgrafptr);

  return (dgraphCoarsenVertLocMax (finegrafptr, (int) flagval));
}

/*+ This routine creates a distributed coarse graph
*** from the given fine graph, unless the coarse graph
*** is smaller than some threshold size or the
*** coarsening ratio is above some other threshold.
*** If the coarse graph is created, a coarse-to-fine
*** vertex array is created, that contains a pair of
*** fine indices for each coarse index. It is the
*** user's responsibility to provide and free this
*** array whenever it is no longer needed.
*** It returns:
*** - 0  : if the graph has been coarsened (and folded).
*** - 1  : if the graph could not be coarsened.
*** - 2  : on error.
+*/

int
SCOTCH_dgraphCoarsen (
SCOTCH_Dgraph * const       flibgrafptr,          /* Fine graph structure              */
const SCOTCH_Num            coarnbr,              /* Minimum number of coarse vertices */
const double                coarrat,              /* Maximum contraction ratio         */
const SCOTCH_Num            flagval,              /* Flag for coarsening type          */
SCOTCH_Dgraph * const       clibgrafptr,          /* Coarse graph structure to fill    */
SCOTCH_Num * const          multloctab)           /* Pointer to multinode array        */
{
  DgraphCoarsenMulti *      multlocptr;           /* Pointer to multinode array used by dgraphCoarsen() */
  int                       foldval;              /* Flag of folding type                               */
  CONTEXTDECL              (flibgrafptr);
  int                       o;

  Dgraph * const            coargrafptr = (Dgraph *) CONTEXTOBJECT (clibgrafptr);
#ifdef SCOTCH_DEBUG_LIBRARY1
  Dgraph * const            finegrafptr = (Dgraph *) CONTEXTOBJECT (flibgrafptr);

  MPI_Comm_compare (finegrafptr->proccomm, coargrafptr->proccomm, &o);
  if ((o != MPI_IDENT) && (o != MPI_CONGRUENT)) {
    errorPrint (STRINGIFY (SCOTCH_dgraphCoarsen) ": communicators are not congruent");
    return (2);
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  if (CONTEXTINIT (flibgrafptr) != 0) {
    errorPrint (STRINGIFY (SCOTCH_dgraphCoarsen) ": cannot initialize context");
    return (2);
  }

  foldval = flagval & (DGRAPHCOARSENFOLD | DGRAPHCOARSENFOLDDUP); /* Only consider folding flags   */
  multlocptr = (foldval == DGRAPHCOARSENNONE)     /* If plain coarsening                           */
               ? ((DgraphCoarsenMulti * restrict) multloctab) /* Use user-provided array (or NULL) */
               : NULL;                            /* Else dgraphCoarsen() will create the array    */

  o = dgraphCoarsen (CONTEXTGETOBJECT (flibgrafptr), coargrafptr, &multlocptr,
                     5, coarnbr, coarrat, (int) flagval, CONTEXTGETDATA (flibgrafptr));
  if (o > 2)
    o = 2;

  if (multloctab == NULL) {                       /* If user did not provide a multinode array      */
    if (multlocptr != NULL)                       /* If coarsening (and folding) went well          */
      memFree (multlocptr);                       /* Free allocated multinode array (folded or not) */
  }
  else {                                          /* User provided a multinode array                */
    if (foldval != DGRAPHCOARSENNONE) {           /* If folding took place, user array was not used */
#ifdef SCOTCH_DEBUG_LIBRARY2
      SCOTCH_Num          multlocsiz;             /* Size that the user may have reserved for the multinode folded array */

      multlocsiz = SCOTCH_dgraphCoarsenVertLocMax (flibgrafptr, flagval);
      if (multlocsiz < coargrafptr->vertlocnbr) { /* If advised size is smaller than real size */
        errorPrint (STRINGIFY (SCOTCH_dgraphCoarsen) ": invalid estimated multinode array size");
        o = 2;
      }                                           /* Go on copying anyway as the user may have reserved more */
#endif /* SCOTCH_DEBUG_LIBRARY2 */
      memCpy  (multloctab, multlocptr, coargrafptr->vertlocnbr * sizeof (DgraphCoarsenMulti)); /* Update array with folded multinode data */
      memFree (multlocptr);                       /* Free allocated folded multinode array */
    }
  }

  CONTEXTEXIT (flibgrafptr);
  return (o);
}
