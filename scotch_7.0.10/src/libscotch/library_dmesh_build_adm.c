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
/**   NAME       : library_dmesh_build_adm.c               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the distri-  **/
/**                buted source mesh building routine of   **/
/**                the libSCOTCH library.                  **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 15 aug 2025     **/
/**                                 to   : 15 aug 2025     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "context.h"
#include "graph.h"
#include "dmesh.h"
#include "ptscotch.h"

/****************************************/
/*                                      */
/* These routines are the C API for the */
/* distributed mesh handling routines.  */
/*                                      */
/****************************************/

/*+ This routine fills the contents of the given
*** opaque distributed mesh structure with the
*** data provided by the user. The base value
*** allows the user to set the graph base to 0 or 1.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_dmeshBuildAdm (
SCOTCH_Dmesh * const        meshptr,
const SCOTCH_Num            baseval,
const SCOTCH_Num            velmlocnbr,           /* Number of local elements              */
SCOTCH_Num * const          velmloctab,           /* Local element vertex begin array      */
const SCOTCH_Num            eelmlocnbr,           /* Number of local element-to-node edges */
SCOTCH_Num * const          eelmloctab,           /* Local edge array                      */
const SCOTCH_Num            vnodglbnbr)
{
  if ((baseval < -1) || (baseval > 1)) {
    errorPrint (STRINGIFY (SCOTCH_dmeshBuildAdm) ": invalid base parameter");
    return (1);
  }

  return (dmeshBuildAdm ((Dmesh * const) CONTEXTOBJECT (meshptr), baseval,
                         velmlocnbr, velmloctab, eelmlocnbr, eelmloctab, vnodglbnbr));
}
