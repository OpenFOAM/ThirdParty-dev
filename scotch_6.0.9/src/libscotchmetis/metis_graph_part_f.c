/* Copyright 2007,2010,2012,2015,2019 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : metis_graph_part_f.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Amaury JACQUES                          **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API of   **/
/**                the compatibility library for the       **/
/**                MeTiS partitioning routines.            **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 10 sep 2006     **/
/**                                 to     07 jun 2007     **/
/**                # Version 5.1  : from : 30 jun 2010     **/
/**                                 to     30 jun 2010     **/
/**                # Version 6.0  : from : 13 sep 2012     **/
/**                                 to     18 may 2019     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

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

FORTRAN (                                                       \
SCOTCH_METIS_V3_PARTGRAPHKWAY, scotch_metis_v3_partgraphkway, ( \
const SCOTCH_Num * const    n,                                  \
const SCOTCH_Num * const    xadj,                               \
const SCOTCH_Num * const    adjncy,                             \
const SCOTCH_Num * const    vwgt,                               \
const SCOTCH_Num * const    adjwgt,                             \
const SCOTCH_Num * const    wgtflag,                            \
const SCOTCH_Num * const    numflag,                            \
const SCOTCH_Num * const    nparts,                             \
const SCOTCH_Num * const    options,                            \
SCOTCH_Num * const          edgecut,                            \
SCOTCH_Num * const          part,                               \
int * const                 revaptr),                           \
(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part, revaptr))
{
  *revaptr = SCOTCH_METIS_V3_PartGraphKway (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

/*
**
*/

FORTRAN (                                                                 \
SCOTCH_METIS_V3_PARTGRAPHRECURSIVE, scotch_metis_v3_partgraphrecursive, ( \
const SCOTCH_Num * const    n,                                            \
const SCOTCH_Num * const    xadj,                                         \
const SCOTCH_Num * const    adjncy,                                       \
const SCOTCH_Num * const    vwgt,                                         \
const SCOTCH_Num * const    adjwgt,                                       \
const SCOTCH_Num * const    wgtflag,                                      \
const SCOTCH_Num * const    numflag,                                      \
const SCOTCH_Num * const    nparts,                                       \
const SCOTCH_Num * const    options,                                      \
SCOTCH_Num * const          edgecut,                                      \
SCOTCH_Num * const          part,                                         \
int * const                 revaptr),                                     \
(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part, revaptr))
{
  *revaptr = SCOTCH_METIS_V3_PartGraphRecursive (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

/*
**
*/

FORTRAN (                                                         \
SCOTCH_METIS_V3_PARTGRAPHVKWAY, scotch_metis_v3_partgraphvkway, ( \
const SCOTCH_Num * const    n,                                    \
const SCOTCH_Num * const    xadj,                                 \
const SCOTCH_Num * const    adjncy,                               \
const SCOTCH_Num * const    vwgt,                                 \
const SCOTCH_Num * const    vsize,                                \
const SCOTCH_Num * const    wgtflag,                              \
const SCOTCH_Num * const    numflag,                              \
const SCOTCH_Num * const    nparts,                               \
const SCOTCH_Num * const    options,                              \
SCOTCH_Num * const          volume,                               \
SCOTCH_Num * const          part,                                 \
int * const                 revaptr),                             \
(n, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part, revaptr))
{
  *revaptr = SCOTCH_METIS_V3_PartGraphVKway (n, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}

/*
**
*/

FORTRAN (                                                       \
SCOTCH_METIS_V5_PARTGRAPHKWAY, scotch_metis_v5_partgraphkway, ( \
const SCOTCH_Num * const    nvtxs,                              \
const SCOTCH_Num * const    ncon,                               \
const SCOTCH_Num * const    xadj,                               \
const SCOTCH_Num * const    adjncy,                             \
const SCOTCH_Num * const    vwgt,                               \
const SCOTCH_Num * const    vsize,                              \
const SCOTCH_Num * const    adjwgt,                             \
const SCOTCH_Num * const    nparts,                             \
const SCOTCH_Num * const    tpwgts,                             \
const double * const        ubvec,                              \
const SCOTCH_Num * const    options,                            \
SCOTCH_Num * const          objval,                             \
SCOTCH_Num * const          part,                               \
int * const                 revaptr),                           \
(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part, revaptr))
{
  *revaptr = SCOTCH_METIS_V5_PartGraphKway (nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part);
}

/*
**
*/

FORTRAN (                                                                 \
SCOTCH_METIS_V5_PARTGRAPHRECURSIVE, scotch_metis_v5_partgraphrecursive, ( \
const SCOTCH_Num * const    nvtxs,                                        \
const SCOTCH_Num * const    ncon,                                         \
const SCOTCH_Num * const    xadj,                                         \
const SCOTCH_Num * const    adjncy,                                       \
const SCOTCH_Num * const    vwgt,                                         \
const SCOTCH_Num * const    vsize,                                        \
const SCOTCH_Num * const    adjwgt,                                       \
const SCOTCH_Num * const    nparts,                                       \
const SCOTCH_Num * const    tpwgts,                                       \
const double * const        ubvec,                                        \
const SCOTCH_Num * const    options,                                      \
SCOTCH_Num * const          objval,                                       \
SCOTCH_Num * const          part,                                         \
int * const                 revaptr),                                     \
(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part, revaptr))
{
  *revaptr = SCOTCH_METIS_V5_PartGraphRecursive (nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part);
}

/*******************/
/*                 */
/* MeTiS v3 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 3)

FORTRAN (                                                             \
METISNAMEU (METIS_PARTGRAPHKWAY), METISNAMEL (metis_partgraphkway), ( \
const SCOTCH_Num * const    n,                                        \
const SCOTCH_Num * const    xadj,                                     \
const SCOTCH_Num * const    adjncy,                                   \
const SCOTCH_Num * const    vwgt,                                     \
const SCOTCH_Num * const    adjwgt,                                   \
const SCOTCH_Num * const    wgtflag,                                  \
const SCOTCH_Num * const    numflag,                                  \
const SCOTCH_Num * const    nparts,                                   \
const SCOTCH_Num * const    options,                                  \
SCOTCH_Num * const          edgecut,                                  \
SCOTCH_Num * const          part),                                    \
(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part))
{
  METISNAMEU (METIS_PartGraphKway) (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

/*
**
*/

FORTRAN (                                                                       \
METISNAMEU (METIS_PARTGRAPHRECURSIVE), METISNAMEL (metis_partgraphrecursive), ( \
const SCOTCH_Num * const    n,                                                  \
const SCOTCH_Num * const    xadj,                                               \
const SCOTCH_Num * const    adjncy,                                             \
const SCOTCH_Num * const    vwgt,                                               \
const SCOTCH_Num * const    adjwgt,                                             \
const SCOTCH_Num * const    wgtflag,                                            \
const SCOTCH_Num * const    numflag,                                            \
const SCOTCH_Num * const    nparts,                                             \
const SCOTCH_Num * const    options,                                            \
SCOTCH_Num * const          edgecut,                                            \
SCOTCH_Num * const          part),                                              \
(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part))
{
  METISNAMEU (METIS_PartGraphRecursive) (n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

/*
**
*/

FORTRAN (                                                               \
METISNAMEU (METIS_PARTGRAPHVKWAY), METISNAMEL (metis_partgraphvkway), ( \
const SCOTCH_Num * const    n,                                          \
const SCOTCH_Num * const    xadj,                                       \
const SCOTCH_Num * const    adjncy,                                     \
const SCOTCH_Num * const    vwgt,                                       \
const SCOTCH_Num * const    vsize,                                      \
const SCOTCH_Num * const    wgtflag,                                    \
const SCOTCH_Num * const    numflag,                                    \
const SCOTCH_Num * const    nparts,                                     \
const SCOTCH_Num * const    options,                                    \
SCOTCH_Num * const          volume,                                     \
SCOTCH_Num * const          part),                                      \
(n, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part))
{
  METISNAMEU (METIS_PartGraphVKway) (n, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}

#endif /* (SCOTCH_METIS_VERSION == 3) */

/*******************/
/*                 */
/* MeTiS v5 stubs. */
/*                 */
/*******************/

#if (SCOTCH_METIS_VERSION == 5)

FORTRAN (                                                             \
METISNAMEU (METIS_PARTGRAPHKWAY), METISNAMEL (metis_partgraphkway), ( \
const SCOTCH_Num * const    nvtxs,                                    \
const SCOTCH_Num * const    ncon,                                     \
const SCOTCH_Num * const    xadj,                                     \
const SCOTCH_Num * const    adjncy,                                   \
const SCOTCH_Num * const    vwgt,                                     \
const SCOTCH_Num * const    vsize,                                    \
const SCOTCH_Num * const    adjwgt,                                   \
const SCOTCH_Num * const    nparts,                                   \
const SCOTCH_Num * const    tpwgts,                                   \
const double * const        ubvec,                                    \
const SCOTCH_Num * const    options,                                  \
SCOTCH_Num * const          objval,                                   \
SCOTCH_Num * const          part),                                    \
(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part))
{
  METISNAMEU (METIS_PartGraphKway) (nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part);
}

/*
**
*/

FORTRAN (                                                                       \
METISNAMEU (METIS_PARTGRAPHRECURSIVE), METISNAMEL (metis_partgraphrecursive), ( \
const SCOTCH_Num * const    nvtxs,                                              \
const SCOTCH_Num * const    ncon,                                               \
const SCOTCH_Num * const    xadj,                                               \
const SCOTCH_Num * const    adjncy,                                             \
const SCOTCH_Num * const    vwgt,                                               \
const SCOTCH_Num * const    vsize,                                              \
const SCOTCH_Num * const    adjwgt,                                             \
const SCOTCH_Num * const    nparts,                                             \
const SCOTCH_Num * const    tpwgts,                                             \
const double * const        ubvec,                                              \
const SCOTCH_Num * const    options,                                            \
SCOTCH_Num * const          objval,                                             \
SCOTCH_Num * const          part),                                              \
(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part))
{
  METISNAMEU (METIS_PartGraphRecursive) (nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part);
}

#endif /* (SCOTCH_METIS_VERSION == 5) */
