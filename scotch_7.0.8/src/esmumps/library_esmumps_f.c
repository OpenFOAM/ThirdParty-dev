/* Copyright 2004,2007,2020 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_esmumps_f.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains Fortran MUMPS      **/
/**                stubs for the ordering routines of the  **/
/**                libSCOTCH + Emilio libfax libraries.    **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 16 may 2001     **/
/**                                 to   : 17 may 2001     **/
/**                # Version 6.0  : from : 22 jan 2020     **/
/**                                 to   : 22 jan 2020     **/
/**                # Version 6.1  : from : 05 sep 2020     **/
/**                                 to   : 05 sep 2020     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "library.h"

/**************************************/
/*                                    */
/* These routines are the Fortran API */
/* for the MUMPS ordering routines.   */
/*                                    */
/**************************************/

ESMUMPS_FORTRAN (                    \
ESMUMPSF, esmumpsf, (                \
const SCOTCH_Num * const    n,       \
const SCOTCH_Num * const    iwlen,   \
SCOTCH_Num * const          petab,   \
const SCOTCH_Num * const    pfree,   \
SCOTCH_Num * const          lentab,  \
SCOTCH_Num * const          iwtab,   \
SCOTCH_Num * const          nvtab,   \
SCOTCH_Num * const          elentab, \
SCOTCH_Num * const          lasttab, \
SCOTCH_Num * const          ncmpa),  \
(n, iwlen, petab, pfree, lentab, iwtab, nvtab, elentab, lasttab, ncmpa))
{
  *ncmpa = esmumps (*n, *iwlen, petab, *pfree, lentab, iwtab, nvtab, elentab, lasttab);
}

/*
**
*/

ESMUMPS_FORTRAN (                    \
ESMUMPSVF, esmumpsvf, (              \
const SCOTCH_Num * const    n,       \
const SCOTCH_Num * const    iwlen,   \
SCOTCH_Num * const          petab,   \
const SCOTCH_Num * const    pfree,   \
SCOTCH_Num * const          lentab,  \
SCOTCH_Num * const          iwtab,   \
SCOTCH_Num * const          nvtab,   \
SCOTCH_Num * const          elentab, \
SCOTCH_Num * const          lasttab, \
SCOTCH_Num * const          ncmpa),  \
(n, iwlen, petab, pfree, lentab, iwtab, nvtab, elentab, lasttab, ncmpa))
{
  *ncmpa = esmumpsv (*n, *iwlen, petab, *pfree, lentab, iwtab, nvtab, elentab, lasttab);
}
