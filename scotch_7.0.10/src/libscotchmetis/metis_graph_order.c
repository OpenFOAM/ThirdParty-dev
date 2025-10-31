/* Copyright 2007,2008,2010,2012,2018,2019,2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_order.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Amaury JACQUES (v6.0)                   **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the MeTiS ordering          **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 08 sep 2006     **/
/**                                 to   : 07 jun 2007     **/
/**                # Version 5.1  : from : 30 jun 2010     **/
/**                                 to   : 30 jun 2010     **/
/**                # Version 6.0  : from : 13 sep 2012     **/
/**                                 to   : 18 may 2019     **/
/**                # Version 6.1  : from : 20 jun 2021     **/
/**                                 to   : 20 jun 2021     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 11 aug 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "metis.h"                                /* Our "metis.h" file */

/************************************/
/*                                  */
/* These routines are the C API for */
/* MeTiS graph ordering routines.   */
/*                                  */
/************************************/

static
int
_SCOTCH_METIS_Node (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  SCOTCH_Graph        grafdat;                    /* Scotch graph object to interface with libScotch    */
  SCOTCH_Ordering     ordedat;                    /* Scotch ordering object to interface with libScotch */
  SCOTCH_Strat        stradat;
  int                 o;

  o = METIS_ERROR;                                /* Assume an error */

  SCOTCH_graphInit (&grafdat);

  if (SCOTCH_graphBuild (&grafdat,
                         *numflag, *n, xadj, xadj + 1, vwgt, NULL,
                         xadj[*n] - *numflag, adjncy, NULL) == 0) {
    SCOTCH_stratInit (&stradat);
#ifdef SCOTCH_DEBUG_ALL
    if (SCOTCH_graphCheck (&grafdat) == 0)        /* TRICK: next instruction called only if graph is consistent */
#endif /* SCOTCH_DEBUG_ALL */
    {
      if (SCOTCH_graphOrderInit (&grafdat, &ordedat, iperm, perm, /* MeTiS and Scotch have opposite definitions for (inverse) permutations */
                                 NULL, NULL, NULL) == 0) {
        if (SCOTCH_graphOrderCompute (&grafdat, &ordedat, &stradat) == 0)
          o = METIS_OK;
        SCOTCH_graphOrderExit (&grafdat, &ordedat);
      }
    }
    SCOTCH_stratExit (&stradat);
  }
  SCOTCH_graphExit (&grafdat);

  return (o);
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V3_EdgeND) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  return (_SCOTCH_METIS_Node (n, xadj, adjncy, NULL, numflag, options, perm, iperm));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V3_NodeND) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  return (_SCOTCH_METIS_Node (n, xadj, adjncy, NULL, numflag, options, perm, iperm));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V3_NodeWND) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  return (_SCOTCH_METIS_Node (n, xadj, adjncy, vwgt, numflag, options, perm, iperm));
}

/*
**
*/

int
SCOTCHMETISNAMES (METIS_V5_NodeND) (
const SCOTCH_Num * const    nvtxs,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  SCOTCH_Num          baseval;

  baseval = ((options != NULL) && (options != xadj)) ? options[METIS_OPTION_NUMBERING] : 0;

  return (_SCOTCH_METIS_Node (nvtxs, xadj, adjncy, vwgt, &baseval, options, perm, iperm));
}

/*******************/
/*                 */
/* MeTiS v3 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 3)

int
SCOTCHMETISNAMEC (METIS_EdgeND) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  return (SCOTCHMETISNAMES (METIS_V3_EdgeND) (n, xadj, adjncy, numflag, options, perm, iperm));
}

/*
**
*/

int
SCOTCHMETISNAMEC (METIS_NodeND) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  return (SCOTCHMETISNAMES (METIS_V3_NodeND) (n, xadj, adjncy, numflag, options, perm, iperm));
}

/*
**
*/

int
SCOTCHMETISNAMEC (METIS_NodeWND) (
const SCOTCH_Num * const    n,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    numflag,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  return (SCOTCHMETISNAMES (METIS_V3_NodeWND) (n, xadj, adjncy, vwgt, numflag, options, perm, iperm));
}

#endif /* (SCOTCH_METIS_VERSION == 3) */

/*******************/
/*                 */
/* MeTiS v5 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 5)

int
SCOTCHMETISNAMEC (METIS_NodeND) (
const SCOTCH_Num * const    nvtxs,
const SCOTCH_Num * const    xadj,
const SCOTCH_Num * const    adjncy,
const SCOTCH_Num * const    vwgt,
const SCOTCH_Num * const    options,
SCOTCH_Num * const          perm,
SCOTCH_Num * const          iperm)
{
  return (SCOTCHMETISNAMES (METIS_V5_NodeND) (nvtxs, xadj, adjncy, vwgt, options, perm, iperm));
}

#endif /* (SCOTCH_METIS_VERSION == 5) */
