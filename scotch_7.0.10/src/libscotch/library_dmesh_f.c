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
/**   NAME       : library_dmesh_f.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This file contains the Fortran API for  **/
/**                the distributed source mesh handling    **/
/**                routines of the libSCOTCH library.      **/
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
#include "ptscotch.h"

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the mesh handling routines.    */
/*                                    */
/**************************************/

/*
**
*/

SCOTCH_FORTRAN (                      \
DMESHINIT, dmeshinit, (               \
SCOTCH_Dmesh * const        meshptr,  \
const MPI_Fint * const      commptr,  \
int * const                 revaptr), \
(meshptr, commptr, revaptr))
{
  MPI_Comm            commdat;

  commdat = MPI_Comm_f2c (*commptr);
  *revaptr = SCOTCH_dmeshInit (meshptr, commdat);
}

/*
**
*/

SCOTCH_FORTRAN (                      \
DMESHEXIT, dmeshexit, (               \
SCOTCH_Dmesh * const        meshptr), \
(meshptr))
{
  SCOTCH_dmeshExit (meshptr);
}

/*
**
*/

SCOTCH_FORTRAN (                      \
DMESHSIZEOF, dmeshsizeof, (           \
int * const                 sizeptr), \
(sizeptr))
{
  *sizeptr = SCOTCH_dmeshSizeof ();
}

/*
**
*/

SCOTCH_FORTRAN (                         \
DMESHSIZE, dmeshsize, (                  \
const SCOTCH_Dmesh * const  meshptr,     \
SCOTCH_Num * const          velmglbptr,  \
SCOTCH_Num * const          velmlocptr,  \
SCOTCH_Num * const          eelmglbptr,  \
SCOTCH_Num * const          eelmlocptr,  \
SCOTCH_Num * const          vnodglbptr), \
(meshptr, velmglbptr, velmlocptr, eelmglbptr, eelmlocptr, vnodglbptr))
{
  SCOTCH_dmeshSize (meshptr, velmglbptr, velmlocptr, eelmglbptr, eelmlocptr, vnodglbptr);
}

/*
**
*/

SCOTCH_FORTRAN (                         \
DMESHDATA, dmeshdata, (                  \
const SCOTCH_Dmesh * const  meshptr,     \
const SCOTCH_Num * const    indxptr,     \
SCOTCH_Num * const          baseptr,     \
SCOTCH_Num * const          velmglbptr,  \
SCOTCH_Num * const          velmlocptr,  \
SCOTCH_Idx * const          velmlocidx,  \
SCOTCH_Num * const          eelmglbptr,  \
SCOTCH_Num * const          eelmlocptr,  \
SCOTCH_Idx * const          eelmlocidx,  \
SCOTCH_Num * const          vnodglbptr,  \
MPI_Fint * const            commptr),    \
(meshptr, indxptr, baseptr,              \
 velmglbptr, velmlocptr, velmlocidx,     \
 eelmglbptr, eelmlocptr, eelmlocidx,     \
 vnodglbptr, commptr))
{
  SCOTCH_Num *        velmloctab;                 /* Pointer to mesh arrays */
  SCOTCH_Num *        eelmloctab;
  MPI_Comm            commdat;

  SCOTCH_dmeshData (meshptr,    baseptr,
                    velmglbptr, velmlocptr, &velmloctab,
                    eelmglbptr, eelmlocptr, &eelmloctab,
                    vnodglbptr, &commdat);

  *velmlocidx = (SCOTCH_Idx) (velmloctab - indxptr) + 1; /* Add 1 since Fortran indices start at 1 */
  *eelmlocidx = (SCOTCH_Idx) (eelmloctab - indxptr) + 1;

  *commptr = MPI_Comm_c2f (commdat);
}
