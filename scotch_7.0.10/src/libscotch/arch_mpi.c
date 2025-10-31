/* Copyright 2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_mpi.c                              **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the MPI-related     **/
/**                generic target architecture functions.  **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 18 feb 2018     **/
/**                                 to   : 19 feb 2018     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define ARCH_MPI

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "arch_mpi.h"
#include "arch_cmplt_mpi.h"
#include "arch_cmpltw_mpi.h"
#include "arch_deco_mpi.h"
#include "arch_deco2_mpi.h"
#include "arch_dist_mpi.h"
#include "arch_hcub_mpi.h"
#include "arch_mesh_mpi.h"
#include "arch_sub_mpi.h"
#include "arch_tleaf_mpi.h"
#include "arch_torus_mpi.h"
#include "arch_vcmplt_mpi.h"
#include "arch_vhcub_mpi.h"

/*
**  The static definitions.
*/

static const ArchMpiClass   archMpiClassTab[] = { ARCHMPICLASSBLOCK (Cmplt),
                                                  ARCHMPICLASSBLOCK (Cmpltw),
                                                  ARCHMPICLASSBLOCK (Deco),
                                                  ARCHMPICLASSBLOCK (Deco2), /* Hidden, type-2 decomposition-defined architecture */
                                                  ARCHMPICLASSBLOCK (Dist),
                                                  ARCHMPICLASSBLOCK (Hcub),
                                                  ARCHMPICLASSBLOCK (Tleaf),
                                                  ARCHMPICLASSBLOCK (Ltleaf),
                                                  ARCHMPICLASSBLOCK (Mesh2),
#ifdef SCOTCH_DEBUG_ARCH3
                                                  ARCHMPICLASSBLOCK (Mesh2o),
                                                  ARCHMPICLASSBLOCK (Mesh2u),
#endif /* SCOTCH_DEBUG_ARCH3 */
                                                  ARCHMPICLASSBLOCK (Mesh3),
                                                  ARCHMPICLASSBLOCK (MeshX),
                                                  ARCHMPICLASSBLOCK (Sub),
                                                  ARCHMPICLASSBLOCK (Torus2),
                                                  ARCHMPICLASSBLOCK (Torus3),
                                                  ARCHMPICLASSBLOCK (TorusX),
                                                  ARCHMPICLASSBLOCK (Vcmplt),
                                                  ARCHMPICLASSBLOCK (Vhcub),
                                                  ARCHMPICLASSBLOCKNULL };

/****************************************/
/*                                      */
/* These are the entry points for the   */
/* MPI-related generic domain routines. */
/*                                      */
/****************************************/

/* This routine creates a MPI datatype for the
** domain of the given target architecture.
** It returns:
** - 0   : MPI type has been created and commited.
** - !0  : on error.
*/

int
archDomMpiType (
const Arch * restrict const   archptr,
MPI_Datatype * restrict const typeptr)
{
  MPI_Datatype        typedat;
  int                 o;

  o = archMpiClassTab [archClassNum (archptr->clasptr)].domMpiType (&archptr->data, &typedat); /* Get base datatype */
  if (o == 0) {
    if ((MPI_Type_create_resized (typedat, 0, sizeof (ArchDom), typeptr) != MPI_SUCCESS) ||
        (MPI_Type_commit (typeptr) != MPI_SUCCESS))
      o = 1;
  }
  MPI_Type_free (&typedat);                       /* Free intermediate datatype */

  return (o);
}

/* This routine creates a MPI datatype of
** a domain comprising a given number of
** consecutive Anum's.
** It returns:
** - 0   : MPI type has been created and commited.
** - !0  : on error.
*/

int
archDomMpiTypeAnum (
int                         anumnbr,
MPI_Datatype * const        typeptr)
{
  return ((MPI_Type_contiguous (anumnbr, ANUM_MPI, typeptr) == MPI_SUCCESS) ? 0 : 1); /* Created MPI types have to be committed */
}
