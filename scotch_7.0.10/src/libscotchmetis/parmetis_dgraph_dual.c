/* Copyright 2023-2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : parmetis_graph_dual.c                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Marc FUENTES                            **/
/**                                                        **/
/**   FUNCTION   : This module is the compatibility        **/
/**                library for the ParMeTiS partitioning   **/
/**                routines containing routines related    **/
/**                to dual graphs                          **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 01 sep 2023     **/
/**                                 to   : 08 aug 2025     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define LIBRARY

#include <mpi.h>
#include "module.h"
#include "common.h"
#include "ptscotch.h"
#include "parmetis.h"                                /* Our "parmetis.h" file */

/* This routine builds a distributed dual graph
** (that is, an element graph) from the given
** distributed mesh data. An edge is built between
** any two element vertices if these two elements
** e1 and e2 have at least *nocoptr nodes in common.
** It returns:
** - METIS_OK   : if the distributed graph has been successfully built.
** - !METIS_OK  : on error.
*/

int
SCOTCHMETISNAMES (ParMETIS_V3_Mesh2Dual) (
const SCOTCH_Num * const    elemdsptab,
const SCOTCH_Num * const    velmloctab,
const SCOTCH_Num * const    eelmloctab,
const SCOTCH_Num * const    baseptr,
const SCOTCH_Num * const    nocoptr,
SCOTCH_Num ** const         vertlocptr,
SCOTCH_Num ** const         edgelocptr,
MPI_Comm *                  commptr)
{
  SCOTCH_Dgraph       grafdat;                    /* Dual graph data structure */
  SCOTCH_Dmesh        meshdat;                    /* Mesh data structure       */
  SCOTCH_Num          vnodlocmax;
  SCOTCH_Num          vnodglbmax;
  SCOTCH_Num          vnodglbnbr;
  SCOTCH_Num          eelmlocnum;
  SCOTCH_Num *        vertloctab;
  SCOTCH_Num *        vertloctmp;
  SCOTCH_Num          vertlocnbr;
  SCOTCH_Num *        edgeloctab;
  SCOTCH_Num *        edgeloctmp;
  SCOTCH_Num          edgelocnbr;
  int                 proclocnum;
  int                 cheklocval;
  int                 chekglbval;
  int                 o;

  MPI_Comm_rank (*commptr, &proclocnum);

  const SCOTCH_Num                  velmlocbas = elemdsptab[proclocnum];
  const SCOTCH_Num                  velmlocnnd = elemdsptab[proclocnum + 1];
  const SCOTCH_Num                  velmlocnbr = velmlocnnd - velmlocbas;
  const SCOTCH_Num * restrict const velmloctax = velmloctab - velmlocbas; /* TRICK: base with respect to global indices */
  const SCOTCH_Num * restrict const eelmloctax = eelmloctab - *baseptr;
  const SCOTCH_Num                  eelmlocnnd = velmloctax[velmlocnnd]; /* Because element vertex array is compact     */
  const SCOTCH_Num                  eelmlocnbr = eelmlocnnd - *baseptr;

  vnodlocmax = 0;
  for (eelmlocnum = *baseptr; eelmlocnum < eelmlocnnd; eelmlocnum ++) { /* Search for highest node vertex index */
    SCOTCH_Num          vnodlocend;

    vnodlocend = eelmloctax[eelmlocnum];
    if (vnodlocmax < vnodlocend)
      vnodlocmax = vnodlocend;
  }

  if (MPI_Allreduce (&vnodlocmax, &vnodglbmax, 1, SCOTCH_NUM_MPI, MPI_MAX, *commptr) != MPI_SUCCESS) {
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: communication error (1)");
    return            (METIS_ERROR);
  }
  vnodglbnbr = vnodglbmax + 1 - *baseptr;

  o = METIS_OK;                                   /* Assume everything will go well */
  if (SCOTCH_dmeshInit (&meshdat, *commptr) != 0) {
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: cannot initialize mesh");
    return            (METIS_ERROR);
  }
  if (SCOTCH_dmeshBuildAdm (&meshdat, *baseptr, velmlocnbr, (SCOTCH_Num *) velmloctab,
                            eelmlocnbr, (SCOTCH_Num *) eelmloctab, vnodglbnbr) != 0) {
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: cannot build mesh");
    o = METIS_ERROR;
    goto abort2;
  }

  if (SCOTCH_dgraphInit (&grafdat, *commptr) != 0) {
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: cannot initialize dual graph");
    o = METIS_ERROR;
    goto abort2;
  }
  if (SCOTCH_dmeshDgraphDual (&meshdat, &grafdat, *nocoptr) != 0) {
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: cannot build dual graph");
    o = METIS_ERROR;
    goto abort1;
  }
  SCOTCH_dgraphData (&grafdat, NULL, NULL, &vertlocnbr, NULL, NULL, &vertloctmp, NULL, NULL, NULL,
                     NULL, &edgelocnbr, NULL, &edgeloctmp, NULL, NULL, NULL);

  cheklocval = 0;                               /* Assume everything will go well                                                   */
  if ((vertloctab = malloc ((vertlocnbr + 1) * sizeof (SCOTCH_Num))) == NULL) { /* Use plain malloc() because of user-managed array */
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: out of memory (1)");
    cheklocval = 1;
  }
  if ((edgeloctab = malloc (edgelocnbr * sizeof (SCOTCH_Num))) == NULL) { /* Use plain malloc() because of user-managed array */
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: out of memory (2)");
    cheklocval = 1;
  }
#ifdef SCOTCH_DEBUG_ALL
  if (MPI_Allreduce (&cheklocval, &chekglbval, 1, MPI_INT, MPI_MAX, *commptr) != MPI_SUCCESS) {
    SCOTCH_errorPrint ("SCOTCH_ParMETIS_V3_Mesh2Dual: communication error");
    return            (METIS_ERROR);
  }
#else /* SCOTCH_DEBUG_ALL */
  chekglbval = cheklocval;
#endif /* SCOTCH_DEBUG_ALL */
  if (chekglbval != 0) {
    if (edgeloctab != NULL)
      memFree (edgeloctab);
    if (vertloctab != NULL)
      memFree (vertloctab);
    o = METIS_ERROR_MEMORY;
    goto abort1;
  }

  memCpy (vertloctab, vertloctmp, (vertlocnbr + 1) * sizeof (SCOTCH_Num)); /* Copy dual graph data to user-managed arrays */
  memCpy (edgeloctab, edgeloctmp, edgelocnbr * sizeof (SCOTCH_Num));

  *vertlocptr = vertloctab;
  *edgelocptr = edgeloctab;

abort1:
  SCOTCH_dgraphExit (&grafdat);
abort2:
  SCOTCH_dmeshExit  (&meshdat);

  return (o);
}

/**********************/
/*                    */
/* ParMeTiS v3 stubs. */
/*                    */
/**********************/

#if (SCOTCH_PARMETIS_VERSION == 3)

int
SCOTCHMETISNAMEC (ParMETIS_V3_Mesh2Dual) (
const SCOTCH_Num * const    elemdsptab,
const SCOTCH_Num * const    velmloctab,
const SCOTCH_Num * const    eelmloctab,
const SCOTCH_Num * const    baseptr,
const SCOTCH_Num * const    nocoptr,
SCOTCH_Num ** const         vertlocptr,
SCOTCH_Num ** const         edgelocptr,
MPI_Comm *                  commptr)
{
  return (SCOTCHMETISNAMES (ParMETIS_V3_Mesh2Dual) (elemdsptab, velmloctab, eelmloctab, baseptr, nocoptr,
                                                    vertlocptr, edgelocptr, commptr));
}

#endif /* (SCOTCH_PARMETIS_VERSION == 3) */
