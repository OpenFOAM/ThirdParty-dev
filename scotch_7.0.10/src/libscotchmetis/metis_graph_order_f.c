/* Copyright 2007,2010,2012,2015,2019,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_order_f.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API of   **/
/**                the compatibility library for the       **/
/**                MeTiS ordering routines.                **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 10 sep 2006     **/
/**                                 to   : 10 sep 2006     **/
/**                # Version 5.1  : from : 30 jun 2010     **/
/**                                 to   : 30 jun 2010     **/
/**                # Version 6.0  : from : 13 sep 2012     **/
/**                                 to   : 18 may 2019     **/
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

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the graph ordering routines.   */
/*                                    */
/**************************************/

FORTRAN (                              \
SCOTCHMETISNAMESU (METIS_V3_EDGEND),   \
SCOTCHMETISNAMESL (metis_v3_edgend), ( \
const SCOTCH_Num * const    n,         \
const SCOTCH_Num * const    xadj,      \
const SCOTCH_Num * const    adjncy,    \
const SCOTCH_Num * const    numflag,   \
const SCOTCH_Num * const    options,   \
SCOTCH_Num * const          perm,      \
SCOTCH_Num * const          iperm,     \
int * const                 revaptr),  \
(n, xadj, adjncy, numflag, options, perm, iperm, revaptr))
{
  *revaptr = SCOTCHMETISNAMES (METIS_V3_EdgeND) (n, xadj, adjncy, numflag, options, perm, iperm);
}

/*
**
*/

FORTRAN (                              \
SCOTCHMETISNAMESU (METIS_V3_NODEND),   \
SCOTCHMETISNAMESL (metis_v3_nodend), ( \
const SCOTCH_Num * const    n,         \
const SCOTCH_Num * const    xadj,      \
const SCOTCH_Num * const    adjncy,    \
const SCOTCH_Num * const    numflag,   \
const SCOTCH_Num * const    options,   \
SCOTCH_Num * const          perm,      \
SCOTCH_Num * const          iperm,     \
int * const                 revaptr),  \
(n, xadj, adjncy, numflag, options, perm, iperm, revaptr))
{
  *revaptr = SCOTCHMETISNAMES (METIS_V3_NodeND) (n, xadj, adjncy, numflag, options, perm, iperm);
}

/*
**
*/

FORTRAN (                               \
SCOTCHMETISNAMESU (METIS_V3_NODEWND),   \
SCOTCHMETISNAMESL (metis_v3_nodewnd), ( \
const SCOTCH_Num * const    n,          \
const SCOTCH_Num * const    xadj,       \
const SCOTCH_Num * const    adjncy,     \
const SCOTCH_Num * const    vwgt,       \
const SCOTCH_Num * const    numflag,    \
const SCOTCH_Num * const    options,    \
SCOTCH_Num * const          perm,       \
SCOTCH_Num * const          iperm,      \
int * const                 revaptr),   \
(n, xadj, adjncy, vwgt, numflag, options, perm, iperm, revaptr))
{
  *revaptr = SCOTCHMETISNAMES (METIS_V3_NodeWND) (n, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}

/*
**
*/

FORTRAN (                              \
SCOTCHMETISNAMESU (METIS_V5_NODEND),   \
SCOTCHMETISNAMESL (metis_v5_nodend), ( \
const SCOTCH_Num * const    nvtxs,     \
const SCOTCH_Num * const    xadj,      \
const SCOTCH_Num * const    adjncy,    \
const SCOTCH_Num * const    vwgt,      \
const SCOTCH_Num * const    options,   \
SCOTCH_Num * const          perm,      \
SCOTCH_Num * const          iperm,     \
int * const                 revaptr),  \
(nvtxs, xadj, adjncy, vwgt, options, perm, iperm, revaptr))
{
  *revaptr = SCOTCHMETISNAMES (METIS_V5_NodeND) (nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
}

/*******************/
/*                 */
/* MeTiS v3 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 3)

FORTRAN (                            \
SCOTCHMETISNAMEFU (METIS_EDGEND),    \
SCOTCHMETISNAMEFL (metis_edgend), (  \
const SCOTCH_Num * const    n,       \
const SCOTCH_Num * const    xadj,    \
const SCOTCH_Num * const    adjncy,  \
const SCOTCH_Num * const    numflag, \
const SCOTCH_Num * const    options, \
SCOTCH_Num * const          perm,    \
SCOTCH_Num * const          iperm),  \
(n, xadj, adjncy, numflag, options, perm, iperm))
{
  SCOTCHMETISNAMEC (METIS_EdgeND) (n, xadj, adjncy, numflag, options, perm, iperm);
}

/*
**
*/

FORTRAN (                            \
SCOTCHMETISNAMEFU (METIS_NODEND),    \
SCOTCHMETISNAMEFL (metis_nodend), (  \
const SCOTCH_Num * const    n,       \
const SCOTCH_Num * const    xadj,    \
const SCOTCH_Num * const    adjncy,  \
const SCOTCH_Num * const    numflag, \
const SCOTCH_Num * const    options, \
SCOTCH_Num * const          perm,    \
SCOTCH_Num * const          iperm),  \
(n, xadj, adjncy, numflag, options, perm, iperm))
{
  SCOTCHMETISNAMEC (METIS_NodeND) (n, xadj, adjncy, numflag, options, perm, iperm);
}

/*
**
*/

FORTRAN (                            \
SCOTCHMETISNAMEFU (METIS_NODEWND),   \
SCOTCHMETISNAMEFL (metis_nodewnd), ( \
const SCOTCH_Num * const    n,       \
const SCOTCH_Num * const    xadj,    \
const SCOTCH_Num * const    adjncy,  \
const SCOTCH_Num * const    vwgt,    \
const SCOTCH_Num * const    numflag, \
const SCOTCH_Num * const    options, \
SCOTCH_Num * const          perm,    \
SCOTCH_Num * const          iperm),  \
(n, xadj, adjncy, vwgt, numflag, options, perm, iperm))
{
  SCOTCHMETISNAMEC (METIS_NodeWND) (n, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}

#endif /* (SCOTCH_METIS_VERSION == 3) */

/*******************/
/*                 */
/* MeTiS v5 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 5)

FORTRAN (                            \
SCOTCHMETISNAMEFU (METIS_NODEND),    \
SCOTCHMETISNAMEFL (metis_nodend), (  \
const SCOTCH_Num * const    nvtxs,   \
const SCOTCH_Num * const    xadj,    \
const SCOTCH_Num * const    adjncy,  \
const SCOTCH_Num * const    vwgt,    \
const SCOTCH_Num * const    options, \
SCOTCH_Num * const          perm,    \
SCOTCH_Num * const          iperm),  \
(nvtxs, xadj, adjncy, vwgt, options, perm, iperm))
{
  SCOTCHMETISNAMEC (METIS_NodeND) (nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
}

#endif /* (SCOTCH_METIS_VERSION == 5) */
