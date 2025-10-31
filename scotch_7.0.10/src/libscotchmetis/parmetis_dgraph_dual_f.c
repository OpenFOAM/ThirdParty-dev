/* Copyright 2024-2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : parmetis_dgraph_dual_f.c                **/
/**                                                        **/
/**   AUTHOR     : Marc FUENTES                            **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API of   **/
/**                the compatibility library for the       **/
/**                ParMeTiS dual mesh building routine.    **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 28 mar 2022     **/
/**                                 to   : 28 jul 2025     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "ptscotch.h"
#include "parmetis.h"                             /* Our "parmetis.h" file */

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the distributed graph ordering */
/* routines.                          */
/*                                    */
/**************************************/

FORTRAN (                                    \
SCOTCHMETISNAMESU (PARMETIS_V3_MESH2DUAL),   \
SCOTCHMETISNAMESL (parmetis_v3_mesh2dual), ( \
const SCOTCH_Num * const    elmdist,         \
SCOTCH_Num * const          eptr,            \
SCOTCH_Num * const          eind,            \
const SCOTCH_Num * const    numflag,         \
const SCOTCH_Num * const    ncommonnodes,    \
SCOTCH_Num **               xadj,            \
SCOTCH_Num **               adjncy,          \
const MPI_Fint * const      commptr,         \
int * const                 revaptr),        \
(elmdist, eptr, eind, numflag, ncommonnodes, xadj, adjncy, commptr, revaptr))
{
  MPI_Comm            commdat;

  commdat = MPI_Comm_f2c (*commptr);
  *revaptr = SCOTCHMETISNAMES (ParMETIS_V3_Mesh2Dual) (elmdist, eptr, eind, numflag, ncommonnodes, xadj, adjncy, &commdat);
}

/*******************/
/*                 */
/* MeTiS v3 stubs. */
/*                 */
/*******************/

#if (SCOTCH_PARMETIS_VERSION == 3)

FORTRAN (                                    \
SCOTCHMETISNAMEFU (PARMETIS_V3_MESH2DUAL),   \
SCOTCHMETISNAMEFL (parmetis_v3_mesh2dual), ( \
const SCOTCH_Num * const    elmdist,         \
SCOTCH_Num * const          eptr,            \
SCOTCH_Num * const          eind,            \
const SCOTCH_Num * const    numflag,         \
const SCOTCH_Num * const    ncommonnodes,    \
SCOTCH_Num **               xadj,            \
SCOTCH_Num **               adjncy,          \
const MPI_Fint * const      commptr),        \
(elmdist, eptr, eind, numflag, ncommonnodes, xadj, adjncy, commptr))
{
  MPI_Comm            commdat;

  commdat = MPI_Comm_f2c (*commptr);
  SCOTCHMETISNAMES (ParMETIS_V3_Mesh2Dual) (elmdist, eptr, eind, numflag, ncommonnodes, xadj, adjncy, &commdat);
}

#endif /* (SCOTCH_PARMETIS_VERSION == 3) */
