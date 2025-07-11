/* Copyright 2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_dgraph_map_stat_f.c             **/
/**                                                        **/
/**   AUTHOR     : Clement BARTHELEMY                      **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the Fortran API for the  **/
/**                distributed mapping handling routines   **/
/**                of the libSCOTCH library.               **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 10 sep 2024     **/
/**                                 to   : 19 sep 2024     **/
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
/* for the mapping routines.          */
/*                                    */
/**************************************/

/*
**
*/

SCOTCH_FORTRAN (                            \
DGRAPHMAPSTAT, dgraphmapstat, (             \
SCOTCH_Dgraph * const         grafptr,      \
const SCOTCH_Dmapping * const mapptr,       \
SCOTCH_Num * const            tgtnbrptr,    \
SCOTCH_Num * const            mapnbrptr,    \
SCOTCH_Num * const            mapminptr,    \
SCOTCH_Num * const            mapmaxptr,    \
double * const                mapavgptr,    \
double * const                mapdltptr,    \
SCOTCH_Num * const            ngbsumptr,    \
SCOTCH_Num * const            ngbminptr,    \
SCOTCH_Num * const            ngbmaxptr,    \
SCOTCH_Num * const            cdstmaxptr,   \
SCOTCH_Num                    cdsttab[256], \
SCOTCH_Num * const            cmlosumptr,   \
SCOTCH_Num * const            cmdisumptr,   \
SCOTCH_Num * const            cmexsumptr),  \
(grafptr, mapptr, tgtnbrptr,                \
 mapnbrptr, mapminptr, mapmaxptr,           \
 mapavgptr, mapdltptr,		            \
 ngbsumptr, ngbminptr, ngbmaxptr,           \
 cdstmaxptr, cdsttab, cmlosumptr, cmdisumptr, cmexsumptr))
{
  SCOTCH_dgraphMapStat (grafptr, mapptr, tgtnbrptr,
                        mapnbrptr, mapminptr, mapmaxptr, mapavgptr, mapdltptr,
                        ngbsumptr, ngbminptr, ngbmaxptr,
                        cdstmaxptr, cdsttab,
                        cmlosumptr, cmdisumptr, cmexsumptr);
}
