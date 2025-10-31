/* Copyright 2021,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_compact.c                        **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file creates compact graphs from   **/
/**                non-compact ones.                       **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 14 sep 2021     **/
/**                                 to   : 17 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "dgraph.h"

/******************************/
/*                            */
/* Graph compacting routines. */
/*                            */
/******************************/

/* This routine creates grouped compact arrays for
** vertex and edge arrays.
** It returns:
** - 0   : if compact arrays have been created.
** - !0  : on error.
*/

int
dgraphCompact2 (
const Dgraph * restrict const orggrafptr,           /*+ Graph to compact                            +*/
Gnum * restrict * const       cmpvertlocptr,        /*+ Pointer to compact vertex array to build    +*/
Gnum * restrict * const       cmpedgelocptr,        /*+ Pointer to compact edge array to build      +*/
Gnum * restrict * const       cmpedlolocptr)        /*+ Pointer to compact edge load array to build +*/
{
  Gnum                datasiz;
  Gnum                cmpvertlocnum;
  Gnum * restrict     cmpvertloctax;
  Gnum * restrict     cmpedgeloctax;
  Gnum * restrict     cmpedloloctax;
  Gnum                cmpedgelocnum;

  const Gnum * restrict const orgvertloctax = orggrafptr->vertloctax;
  const Gnum * restrict const orgvendloctax = orggrafptr->vendloctax;
  const Gnum * restrict const orgedgeloctax = orggrafptr->edgeloctax;
  const Gnum * restrict const orgedloloctax = orggrafptr->edloloctax;

  datasiz = (orggrafptr->vertlocnbr + 1) + orggrafptr->edgelocnbr + ((orgedloloctax != NULL) ? orggrafptr->edgelocnbr : 0);

  if ((cmpvertloctax = memAlloc (datasiz * sizeof (Gnum))) == NULL) {
    errorPrint ("dgraphCompact2: out of memory");
    return (1);
  }
  cmpvertloctax -= orggrafptr->baseval;
  cmpedgeloctax  = cmpvertloctax + (orggrafptr->vertlocnbr + 1);
  cmpedloloctax  = (orgedloloctax != NULL) ? (cmpedgeloctax + orggrafptr->edgelocnbr) : NULL;

  for (cmpvertlocnum = cmpedgelocnum = orggrafptr->baseval; cmpvertlocnum < orggrafptr->vertlocnnd; cmpvertlocnum ++) { /* For all local vertices */
    Gnum                degrval;
    Gnum                orgedgelocnum;

    cmpvertloctax[cmpvertlocnum] = cmpedgelocnum; /* Set beginning of local edge sub-array */
    orgedgelocnum = orgvertloctax[cmpvertlocnum];
    degrval = orgvendloctax[cmpvertlocnum] - orgedgelocnum;

#ifdef SCOTCH_DEBUG_DGRAPH2
    if ((cmpedgelocnum + degrval) > (orggrafptr->edgelocnbr + orggrafptr->baseval)) {
      errorPrint ("dgraphCompact2: internal error");
      return (1);
    }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

    memCpy (cmpedgeloctax + cmpedgelocnum, orgedgeloctax + orgedgelocnum, degrval * sizeof (Gnum)); /* Copy edge sub-array */
    if (orgedloloctax != NULL)
      memCpy (cmpedloloctax + cmpedgelocnum, orgedloloctax + orgedgelocnum, degrval * sizeof (Gnum)); /* Copy edge load sub-array */

    cmpedgelocnum += degrval;                     /* Advance compact edge index by number of end vertices */
  }
  cmpvertloctax[cmpvertlocnum] = cmpedgelocnum;   /* Set end of compact edge array */

  *cmpvertlocptr = cmpvertloctax;
  *cmpedgelocptr = cmpedgeloctax;
  *cmpedlolocptr = cmpedloloctax;

  return (0);
}
