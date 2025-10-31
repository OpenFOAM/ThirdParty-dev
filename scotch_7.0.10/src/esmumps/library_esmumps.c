/* Copyright 2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_esmumps.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains a MUMPS interface  **/
/**                for the ordering routines of the        **/
/**                libSCOTCH + Emilio libfax libraries.    **/
/**                                                        **/
/**   DATES      : # Version 6.1  : from : 05 sep 2020     **/
/**                                 to   : 05 sep 2020     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/**   NOTES      : # This code derives from that of the    **/
/**                  original "esmumps.c" file.            **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "esmumps.h"
#include "library.h"

/**************************************/
/*                                    */
/* This routine acts as an interface  */
/* between ordering software such as  */
/* MUMPS and Scotch+Emilio.           */
/*                                    */
/**************************************/

int
esmumps (
const SCOTCH_Num            n,
const SCOTCH_Num            iwlen,                /*+ Not used, just here for consistency      +*/
SCOTCH_Num * restrict const petab,
const SCOTCH_Num            pfree,
SCOTCH_Num * restrict const lentab,
SCOTCH_Num * restrict const iwtab,
SCOTCH_Num * restrict const nvtab,
SCOTCH_Num * restrict const elentab,              /*+ Permutations computed for debugging only +*/
SCOTCH_Num * restrict const lasttab)              /*+ Permutations computed for debugging only +*/
{
  return (esmumps2 (n, pfree, petab, lentab, iwtab, NULL, nvtab, elentab, lasttab));
}

/*
**
*/

int
esmumpsv (
const SCOTCH_Num            n,
const SCOTCH_Num            iwlen,                /*+ Not used, just here for consistency      +*/
SCOTCH_Num * restrict const petab,
const SCOTCH_Num            pfree,
SCOTCH_Num * restrict const lentab,
SCOTCH_Num * restrict const iwtab,
SCOTCH_Num * restrict const nvtab,
SCOTCH_Num * restrict const elentab,              /*+ Permutations computed for debugging only +*/
SCOTCH_Num * restrict const lasttab)              /*+ Permutations computed for debugging only +*/
{
  return (esmumps2 (n, pfree, petab, lentab, iwtab, nvtab, nvtab, elentab, lasttab));
}
