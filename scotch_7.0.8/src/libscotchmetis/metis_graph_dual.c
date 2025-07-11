/* Copyright 2020,2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_dual.c                      **/
/**                                                        **/
/**   AUTHOR     : Marc FUENTES                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the MeTiS partitioning      **/
/**                routines containing routines relative   **/
/**                to dual graphs                          **/
/**                                                        **/
/**   DATES      : # Version 6.1  : from : 01 sep 2020     **/
/**                                 to   : 28 may 2021     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 11 aug 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "graph.h"
#include "mesh.h"
#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */
#include "metis_graph_dual.h"                     /* Our "metis.h" file */

/* This routine creates a (SCOTCH_)Mesh structure
** from the partial mesh connectivity data that is
** passed to it.
** It returns:
** - METIS_OK      : if the mesh has been successfully built.
** - METIS_ERROR*  : on error.
*/

int
_SCOTCH_METIS_MeshToDual2 (
SCOTCH_Mesh * const         meshptr,              /*+ Mesh structure to fill                        +*/
const SCOTCH_Num            baseval,              /*+ Base value                                    +*/
const SCOTCH_Num            vnodnbr,              /*+ Number of nodes in mesh                       +*/
const SCOTCH_Num            velmnbr,              /*+ Number of elements in mesh                    +*/
const SCOTCH_Num * const    verttab,              /*+ Array of start indices of elements in edgetab +*/
const SCOTCH_Num * const    edgetab)              /*+ Array of elements                             +*/
{
  Gnum * restrict     srcverttax;
  Gnum * restrict     srcedgetax;
  Gnum                degrmax;
  Gnum                edgenbr;
  Gnum                edgenum;
  Gnum                vertnum;

  Mesh * const                srcmeshptr = (Mesh *) meshptr; /* Use structure as source mesh */
  const Gnum * restrict const verttax = verttab - baseval;
  const Gnum * restrict const edgetax = edgetab - baseval;
  const Gnum                  velmnnd = velmnbr + baseval;

#ifdef SCOTCH_DEBUG_LIBRARY1
  for (vertnum = baseval; vertnum < velmnnd; vertnum ++) { /* For all element vertices */
    Gnum                edgenum;
    Gnum                edgennd;

    edgenum = verttax[vertnum];
    edgennd = verttax[vertnum + 1];
    if (edgennd < edgenum) {
      SCOTCH_errorPrint ("_SCOTCH_METIS_MeshToDual2: invalid input indices (1)");
      return            (METIS_ERROR_INPUT);
    }
    for ( ; edgenum < edgennd; edgenum ++) {
      if ((edgetax[edgenum] <  baseval) ||
          (edgetax[edgenum] > (baseval + vnodnbr))) {
        SCOTCH_errorPrint ("_SCOTCH_METIS_MeshToDual2: invalid input indices (2)");
        return            (METIS_ERROR_INPUT);
      }
    }
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  srcmeshptr->flagval = MESHFREEEDGE | MESHFREEVERT;
  srcmeshptr->baseval = baseval;
  srcmeshptr->velmbas = baseval;                  /* Elements are placed first */
  srcmeshptr->velmnbr = velmnbr;
  srcmeshptr->velmnnd = velmnnd;
  srcmeshptr->vnodbas = velmnnd;                  /* Nodes are placed after */
  srcmeshptr->vnodnbr = vnodnbr;
  srcmeshptr->vnodnnd = vnodnbr + velmnnd;
  srcmeshptr->velotax = NULL;                     /* No element loads */
  srcmeshptr->velosum = velmnbr;
  srcmeshptr->vnlotax = NULL;                     /* No node loads */
  srcmeshptr->vnlosum = vnodnbr;

  if ((srcverttax = memAlloc ((velmnbr + vnodnbr + 1) * sizeof (Gnum))) == NULL) { /* Allocate compact vertex array for mesh vertices */
    SCOTCH_errorPrint ("_SCOTCH_METIS_MeshToDual2: out of memory (1)");
    return            (METIS_ERROR_MEMORY);
  }
  memSet (srcverttax + velmnbr, 0, vnodnbr * sizeof (Gnum)); /* Initialize node part of array as node vertex degrees */
  srcverttax -= baseval;                          /* Array is now based                                              */
  srcmeshptr->verttax = srcverttax;
  srcmeshptr->vendtax = srcverttax + 1;

  degrmax = 0;
  edgenbr = 0;
  for (vertnum = baseval; vertnum < velmnnd; vertnum ++) { /* For all element vertices */
    Gnum                edgenum;
    Gnum                edgennd;
    Gnum                degrval;

    edgenum = verttax[vertnum];
    edgennd = verttax[vertnum + 1];
    degrval = edgennd - edgenum;
    if (degrval > degrmax)                        /* Compute vertex maximum degree for element vertices */
      degrmax = degrval;
    edgenbr += degrval;                           /* Accumulate number of element-node arcs */

    for ( ; edgenum < edgennd; edgenum ++)
      srcverttax[velmnbr + edgetax[edgenum]] ++;  /* Accumulate degrees of end node vertices */
  }
  edgenbr *= 2;                                   /* Overall number of arcs is twice the number of element-node arcs */
  srcmeshptr->edgenbr = edgenbr;

  if (verttax[baseval] == baseval)                /* If vertex indices are aligned with graph base value                              */
    memCpy (srcverttax + baseval, verttab, velmnbr * sizeof (Gnum)); /* Copy element adjacency list into element part of vertex array */
  else {
    Gnum                edgeadj;

    edgeadj = verttax[baseval] - baseval;
    for (vertnum = baseval; vertnum < velmnnd; vertnum ++)
      srcverttax[vertnum] = verttax[vertnum] + edgeadj;
  }

  for (vertnum = velmnnd, edgenum = verttax[velmnnd]; vertnum < velmnnd + vnodnbr; vertnum ++) { /* For all node vertex indices in mesh */
    Gnum                degrval;

    degrval = srcverttax[vertnum];
    if (degrval > degrmax)                        /* Compute vertex maximum degree for node vertices */
      degrmax = degrval;

    srcverttax[vertnum] = edgenum;                /* Fill node part of vertex array */
    edgenum += degrval;
  }
  srcverttax[vertnum] = edgenum;                  /* Mark end of vertex array */
  srcmeshptr->degrmax = degrmax;

  if ((srcedgetax = memAlloc (edgenbr * sizeof (Gnum))) == NULL) {
    SCOTCH_errorPrint ("_SCOTCH_METIS_MeshToDual2: out of memory (2)");
    memFree (srcverttax + baseval);
    return  (METIS_ERROR_MEMORY);
  }
  srcedgetax -= baseval;
  srcmeshptr->edgetax = srcedgetax;

  for (edgenum = baseval; edgenum < verttax[velmnnd]; edgenum ++) /* Copy skewed edge array for element vertices */
    srcedgetax[edgenum] = edgetax[edgenum] + velmnbr;

  for (vertnum = baseval; vertnum < velmnnd; vertnum ++) {
    Gnum                edgenum;

    for (edgenum = verttax[vertnum]; edgenum < verttax[vertnum + 1]; edgenum ++) /* Build edge array for node vertices */
      srcedgetax[srcverttax[edgetax[edgenum] + velmnbr] ++] = vertnum;
  }

  memMov (srcverttax + velmnnd + 1, srcverttax + velmnnd, (vnodnbr - 1) * sizeof (Gnum));  /* Re-build node part of vertex array */
  srcverttax[velmnnd] = verttax[velmnnd];

  return (METIS_OK);
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_MeshToDual) (
const SCOTCH_Num * const    ne,
const SCOTCH_Num * const    nn,
const SCOTCH_Num * const    eptr,
const SCOTCH_Num * const    eind,
const SCOTCH_Num * const    ncommon,
const SCOTCH_Num * const    nuimflag,
SCOTCH_Num ** const         xadj,
SCOTCH_Num ** const         adjncy)
{
  SCOTCH_Mesh         meshdat;
  SCOTCH_Graph        grafdat;
  SCOTCH_Num          baseval;
  SCOTCH_Num          vertnbr;
  SCOTCH_Num *        verttab;
  SCOTCH_Num *        vendtab;
  SCOTCH_Num          edgenbr;
  SCOTCH_Num *        edgetab;
  int                 o;

  *xadj = NULL;                                   /* Assume something will go wrong */

  SCOTCH_meshInit  (&meshdat);
  SCOTCH_graphInit (&grafdat);

  o = _SCOTCH_METIS_MeshToDual2 (&meshdat, *nuimflag, *nn, *ne, eptr, eind);
  if (o != METIS_OK) {
    SCOTCH_errorPrint ("METIS_MeshToDual: cannot create mesh");
    return (o);
  }
  o = SCOTCH_meshGraphDual (&meshdat, &grafdat, *ncommon);
  SCOTCH_meshExit (&meshdat);                     /* Mesh structure is no longer needed */
  if (o != 0) {
    SCOTCH_errorPrint ("METIS_MeshToDual: cannot create graph from mesh");
    return (o);
  }

  SCOTCH_graphData (&grafdat, &baseval, &vertnbr, &verttab, &vendtab, NULL, NULL, &edgenbr, &edgetab, NULL);

  if (((*xadj   = malloc ((vertnbr + 1) * sizeof (SCOTCH_Num))) == NULL) || /* Do not use libScotch memory allocation as freed by user */
      ((*adjncy = malloc (edgenbr       * sizeof(SCOTCH_Num)))  == NULL)) {
    SCOTCH_errorPrint ("METIS_MeshToDual: out of memory");
    if (*xadj != NULL)
      free (*xadj);
    SCOTCH_graphExit (&grafdat);
    return           (METIS_ERROR_MEMORY);
  }

  memCpy (*xadj,   verttab, (vertnbr + 1) * sizeof (SCOTCH_Num));
  memCpy (*adjncy, edgetab,  edgenbr      * sizeof (SCOTCH_Num));
  SCOTCH_graphExit (&grafdat);

  return (METIS_OK);
}

/*
**
*/

int
SCOTCHMETISNAMEC (METIS_MeshToDual) (
const SCOTCH_Num * const    ne,
const SCOTCH_Num * const    nn,
const SCOTCH_Num * const    eptr,
const SCOTCH_Num * const    eind,
const SCOTCH_Num * const    ncommon,
const SCOTCH_Num * const    nuimflag,
SCOTCH_Num ** const         xadj,
SCOTCH_Num ** const         adjncy)
{
  return (SCOTCHMETISNAMES (METIS_MeshToDual) (ne, nn, eptr, eind, ncommon, nuimflag, xadj, adjncy));
}
